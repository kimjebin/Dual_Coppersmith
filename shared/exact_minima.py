"""Certified lambda1/lambda2 via checkpointed, parallel integer enumeration.

fpLLL reduces the basis and its root-orthogonal quotient. A GMP C++ worker then
uses fraction-free Bareiss data for every sphere exclusion and exact nearest
integer lifts. No floating pruning, candidate limit, or root-shortest assumption
is used. An unfinished run has ``exact=False`` and no lambda2 result.
"""

from __future__ import annotations

from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from contextlib import contextmanager
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import threading
import time

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

SOURCE = Path(__file__).with_name("exact_enum.cpp")
GIB = 1 << 30
PREPARATION_VERSION = "exact-integer-minima-preparation-v2"


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for part in iter(lambda: handle.read(1 << 20), b""):
            digest.update(part)
    return digest.hexdigest()


def _digest(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                    separators=(",", ":")).encode()).hexdigest()


def _atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, sort_keys=True, separators=(",", ":"))
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


@contextmanager
def _lock(folder: Path):
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / ".lock").open("a+") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"Another exact search is using {folder}") from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def _build_native(folder: Path) -> Path:
    compiler = shutil.which("g++") or shutil.which("c++")
    if compiler is None:
        raise RuntimeError("A C++17 compiler and GMP development headers are required")
    source_hash = _hash_file(SOURCE)
    native_folder = folder / ".native" / (platform.system() + "-" + platform.machine())
    native_folder.mkdir(parents=True, exist_ok=True)
    binary = native_folder / ("exact_enum-" + source_hash[:16])
    if not binary.exists():
        temporary = binary.with_suffix(".tmp")
        command = [compiler, "-O3", "-std=c++17", str(SOURCE),
                   "-lgmpxx", "-lgmp", "-o", str(temporary)]
        attempt = subprocess.run(command + ["-fopenmp"], capture_output=True, timeout=180)
        if attempt.returncode:
            print("OpenMP build unavailable; compiling serial preparation", flush=True)
            subprocess.run(command, check=True, timeout=180)
        os.replace(temporary, binary)
    return binary.resolve()


def _nearest(numerator: int, denominator: int) -> int:
    return (2 * numerator + denominator) // (2 * denominator)


def _witness(data: dict, coefficients: list[int]) -> tuple[int, list[int]]:
    if len(coefficients) != len(data["tail"]) or not any(coefficients):
        raise ArithmeticError("Invalid nonzero quotient witness")
    dot = sum(int(c) * int(y) for c, y in zip(data["inner"], coefficients))
    k = _nearest(dot, data["root_sq"])
    vector = [sum(int(y) * int(row[j]) for y, row in
                  zip(coefficients, data["tail"])) - k * int(data["root"][j])
              for j in range(len(data["root"]))]
    return sum(value * value for value in vector), vector


def _matrix_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def _proved_lll(basis, delta):
    return basis.LLL(delta=delta, algorithm="fpLLL:proved", transformation=True)


def _validated_basis(basis):
    from sage.all import ZZ, matrix
    raw_rows = basis.rows() if hasattr(basis, "rows") else basis
    rows = []
    for row in raw_rows:
        checked = []
        for value in row:
            if (isinstance(value, float) or
                (hasattr(value, "parent") and not value.parent().is_exact()) or ZZ(value) != value):
                raise ValueError("Exact minima require integral, non-floating matrix entries")
            checked.append(ZZ(value))
        rows.append(checked)
    return matrix(ZZ, rows)


def _save_document(path, doc):
    doc["payload_sha256"] = _digest(doc)
    _atomic_json(path, doc)
    return doc


def _read_document(path):
    doc = json.loads(Path(path).read_text())
    if doc.get("payload_sha256") != _digest({k: v for k, v in doc.items()
                                              if k != "payload_sha256"}):
        raise ValueError(f"Prepared metadata hash mismatch: {path}")
    return doc


def _make_preparation(basis, folder: Path, binary: Path, *, quotient_delta: float = 0.99,
                      preparation_threads: int = 1, reuse: Path | None = None) -> dict:
    """Persist every LLL stage; mathematical identities do not bind host binaries.

    A legacy final basis may warm-start a NEW ambient LLL, but is never relabelled
    as the untouched output of the old ambient LLL. Its equivalence is verified.
    No LLL or Bareiss preparation here has a wall-clock timeout.
    """
    from sage.all import ZZ, identity_matrix, matrix, vector
    from sage.env import SAGE_VERSION

    original = matrix(ZZ, basis)
    input_rows = _matrix_rows(original)
    input_hash = _digest(input_rows)
    metadata = folder / "prepared.json"
    engine_hash = _digest([_hash_file(SOURCE), PREPARATION_VERSION, quotient_delta])
    identity = _digest([input_hash, engine_hash])
    if metadata.exists():
        doc = _read_document(metadata)
        if doc["identity"] != identity:
            raise ValueError("Checkpoint belongs to another basis or engine version; use a new directory for legacy import")
        if doc["preparation_sha256"] != _hash_file(folder / "integer.prepared"):
            raise ValueError("Integer preparation hash mismatch")
        for stage in ("input_basis", "ambient_lll", "quotient_lll"):
            if doc[stage + "_file_sha256"] != _hash_file(folder / (stage + ".json")):
                raise ValueError(f"Preparation stage hash mismatch: {stage}")
        return _decode_preparation(doc)
    started = time.monotonic()
    n = original.nrows()
    if n < 2 or original.ncols() < n:
        raise ValueError("A full row-rank integer basis of rank at least two is required")
    input_path = folder / "input_basis.json"
    if input_path.exists():
        if _read_document(input_path)["rows"] != input_rows:
            raise ValueError("Saved input basis mismatch")
    else:
        _save_document(input_path, {"schema": "integer-basis-v1", "rows": input_rows,
                                   "rank": n, "columns": original.ncols()})
    legacy = None
    if reuse is not None and (reuse / "prepared.json").exists():
        legacy = _read_document(reuse / "prepared.json")
        if legacy["input_basis_sha256"] != input_hash:
            raise ValueError("Legacy input basis mismatch")
        if legacy["preparation_sha256"] != _hash_file(reuse / "integer.prepared"):
            raise ValueError("Legacy integer preparation hash mismatch")
        if legacy.get("native_source_sha256") is not None:
            if legacy["native_source_sha256"] != _hash_file(SOURCE):
                raise ValueError("Legacy native source is incompatible")
        elif legacy["binary_sha256"] != _hash_file(binary):
            raise ValueError("Legacy v1 import requires the matching audited native implementation")
    ambient_path = folder / "ambient_lll.json"
    if ambient_path.exists():
        ambient = _read_document(ambient_path)
        if ambient["input_basis_sha256"] != input_hash or ambient["delta"] != 0.8:
            raise ValueError("Ambient LLL identity mismatch")
        reduced = matrix(ZZ, ambient["reduced_rows"])
        warm = matrix(ZZ, ambient["lll_input_rows"])
        preconditioner = matrix(ZZ, ambient["preconditioner_rows"])
        transform = matrix(ZZ, ambient["transform_rows"])
        if (preconditioner * original != warm or transform * warm != reduced or
                abs(preconditioner.det()) != 1 or abs(transform.det()) != 1):
            raise ArithmeticError("Invalid saved ambient transformations")
    else:
        warm = original
        preconditioner = identity_matrix(ZZ, n)
        provenance = "original_public_basis"
        if legacy is not None:
            warm = matrix(ZZ, [legacy["root"]] + legacy["tail"])
            preconditioner = original.solve_left(warm).change_ring(ZZ)
            if preconditioner * original != warm or abs(preconditioner.det()) != 1:
                raise ArithmeticError("Legacy basis is not integrally unimodularly equivalent")
            provenance = "verified_legacy_final_basis_warm_start"
        print("EXACT_STAGE " + json.dumps({"stage": "ambient_lll", "rows": n,
                                            "columns": original.ncols(), "input": provenance}), flush=True)
        stage_started = time.monotonic()
        reduced, transform = _proved_lll(warm, 0.8)
        seconds = time.monotonic() - stage_started
        if transform * warm != reduced or abs(transform.det()) != 1:
            raise ArithmeticError("Invalid fpLLL ambient transformation")
        ambient = _save_document(ambient_path, {
            "schema": "ambient-fplll-stage-v1", "input_basis_sha256": input_hash,
            "backend": "fpLLL:proved", "delta": 0.8, "sage_version": SAGE_VERSION,
            "input_provenance": provenance, "preconditioner_rows": _matrix_rows(preconditioner),
            "lll_input_rows": _matrix_rows(warm), "transform_rows": _matrix_rows(transform),
            "reduced_rows": _matrix_rows(reduced), "lll_seconds": seconds,
            "legacy_payload_sha256": legacy["payload_sha256"] if legacy else None})
    order = sorted(range(n), key=lambda i: reduced.row(i).dot_product(reduced.row(i)))
    root = vector(ZZ, reduced.row(order[0]))
    root_sq = int(root.dot_product(root))
    if root_sq <= 0:
        raise ValueError("The basis is not of full row rank")
    tail = matrix(ZZ, [reduced.row(i) for i in order[1:]])
    inner = tail * root
    projection = root_sq * tail - matrix(ZZ, [c * root for c in inner])
    coordinate_content = math.gcd(*(int(value) for value in projection.list()))
    if coordinate_content <= 0:
        raise ValueError("The quotient basis is not of full row rank")
    projection = projection.apply_map(lambda value: value // coordinate_content)
    quotient_path = folder / "quotient_lll.json"
    if quotient_path.exists():
        quotient = _read_document(quotient_path)
        if (quotient["ambient_file_sha256"] != _hash_file(ambient_path) or
                quotient["delta"] != quotient_delta or quotient["input_rows"] != _matrix_rows(projection)):
            raise ValueError("Quotient LLL identity mismatch")
        transform = matrix(ZZ, quotient["transform_rows"])
        projected_reduced = matrix(ZZ, quotient["reduced_rows"])
    else:
        print("EXACT_STAGE " + json.dumps({"stage": "quotient_lll", "rows": n - 1,
                                            "columns": projection.ncols()}), flush=True)
        stage_started = time.monotonic()
        projected_reduced, transform = _proved_lll(projection, quotient_delta)
        seconds = time.monotonic() - stage_started
        if transform * projection != projected_reduced or abs(transform.det()) != 1:
            raise ArithmeticError("Invalid fpLLL quotient transformation")
        quotient = _save_document(quotient_path, {
            "schema": "quotient-fplll-stage-v1", "backend": "fpLLL:proved",
            "delta": quotient_delta, "sage_version": SAGE_VERSION,
            "ambient_file_sha256": _hash_file(ambient_path), "ambient_row_order": order,
            "coordinate_content": str(coordinate_content), "input_rows": _matrix_rows(projection),
            "reduced_rows": _matrix_rows(projected_reduced), "transform_rows": _matrix_rows(transform),
            "lll_seconds": seconds})
    if transform * projection != projected_reduced or abs(transform.det()) != 1:
        raise ArithmeticError("Invalid saved quotient transformation")
    tail = transform * tail
    root_shifts = [_nearest(int(row.dot_product(root)), root_sq) for row in tail.rows()]
    rows = [row - k * root for row, k in zip(tail.rows(), root_shifts)]
    tail = matrix(ZZ, rows)
    inner = [int(row.dot_product(root)) for row in rows]
    gram = root_sq * (tail * tail.transpose()) - matrix(
        ZZ, len(inner), len(inner), lambda i, j: inner[i] * inner[j])
    content = math.gcd(*(int(value) for value in gram.list()))
    if content <= 0:
        raise ValueError("The quotient basis is not of full row rank")
    gram = gram.apply_map(lambda value: value // content)
    best_index = min(range(n - 1), key=lambda i: rows[i].dot_product(rows[i]))
    best_y = [int(i == best_index) for i in range(n - 1)]
    best = int(rows[best_index].dot_product(rows[best_index]))
    same_coordinates = legacy is not None and all((
        legacy["root"] == list(map(str, root)), legacy["tail"] == _matrix_rows(tail),
        legacy["inner"] == list(map(str, inner)), legacy["content"] == str(content),
        legacy["root_sq"] == str(root_sq)))
    temporary = folder / "integer.input"
    with temporary.open("w") as handle:
        handle.write(f"EXACT-GMP-INPUT-1 {n - 1} {root_sq} {content}\n")
        handle.write(" ".join(map(str, inner)) + "\n")
        for row in gram.rows():
            handle.write(" ".join(map(str, row)) + "\n")
    if same_coordinates:
        shutil.copyfile(reuse / "integer.prepared", folder / "integer.prepared.tmp")
    else:
        print("EXACT_STAGE " + json.dumps({"stage": "integer_preparation", "dimension": n - 1}), flush=True)
        subprocess.run([str(binary), "prepare", str(temporary), str(folder / "integer.prepared.tmp"),
                        str(preparation_threads)], check=True)
    os.replace(folder / "integer.prepared.tmp", folder / "integer.prepared")
    doc = {
        "schema": PREPARATION_VERSION, "identity": identity,
        "input_basis_sha256": input_hash, "engine_sha256": engine_hash,
        "native_source_sha256": _hash_file(SOURCE), "binary_sha256": _hash_file(binary),
        "root": list(map(str, root)), "root_sq": str(root_sq), "tail": _matrix_rows(tail),
        "inner": list(map(str, inner)), "content": str(content),
        "root_shifts": list(map(str, root_shifts)), "initial_best": str(best), "initial_y": best_y,
        "preparation_sha256": _hash_file(folder / "integer.prepared"),
        "preparation_seconds": time.monotonic() - started,
        "ambient_lll_seconds": ambient["lll_seconds"], "quotient_lll_seconds": quotient["lll_seconds"],
        "quotient_delta": quotient_delta, "preparation_threads_requested": preparation_threads,
        "legacy_native_reused": same_coordinates,
        **{stage + "_file_sha256": _hash_file(folder / (stage + ".json"))
           for stage in ("input_basis", "ambient_lll", "quotient_lll")},
    }
    doc["frontier_identity"] = _digest(["exact-frontier-v2", input_hash,
                                        doc["root"], doc["tail"], doc["inner"], doc["content"],
                                        doc["preparation_sha256"], doc["native_source_sha256"]])
    # SQLite backup includes WAL contents and never alters the legacy campaign.
    if same_coordinates and (reuse / "frontier.sqlite3").exists() and not (folder / "frontier.sqlite3").exists():
        source_conn = sqlite3.connect(f"file:{reuse / 'frontier.sqlite3'}?mode=ro", uri=True)
        descriptor, import_path = tempfile.mkstemp(prefix="frontier-import-", suffix=".sqlite3", dir=folder)
        os.close(descriptor)
        target_conn = sqlite3.connect(import_path)
        try:
            source_conn.backup(target_conn)
            previous = target_conn.execute("SELECT value FROM meta WHERE key='identity'").fetchone()
            if previous is None or previous[0] != legacy.get("frontier_identity", legacy["payload_sha256"]):
                raise ValueError("Legacy frontier identity mismatch")
            with target_conn:
                target_conn.execute("UPDATE meta SET value=? WHERE key='identity'", (doc["frontier_identity"],))
            target_conn.execute("PRAGMA journal_mode=DELETE")
        finally:
            source_conn.close()
            target_conn.close()
        os.replace(import_path, folder / "frontier.sqlite3")
    _save_document(metadata, doc)
    return _decode_preparation(doc)


def _decode_preparation(doc: dict) -> dict:
    return {"raw": doc, "root": list(map(int, doc["root"])),
            "root_sq": int(doc["root_sq"]),
            "tail": [list(map(int, row)) for row in doc["tail"]],
            "inner": list(map(int, doc["inner"])), "content": int(doc["content"])}


def save_minima_input(basis, checkpoint, instance=None):
    """Save actual public input and integer basis, not just a PRNG seed."""
    folder = Path(checkpoint).resolve()
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / "input.json"
    rows = _matrix_rows(basis)
    if path.exists():
        previous = _read_document(path)
        if previous["basis_rows"] != rows or (instance is not None and _digest(previous["instance"]) != _digest(instance)):
            raise ValueError("Saved public input differs from the requested instance")
        return previous
    return _save_document(path, {"schema": "portable-minima-input-v1", "basis_rows": rows,
                                 "input_basis_sha256": _digest(rows), "instance": instance,
                                 "lattice_dimension": basis.nrows(), "ambient_columns": basis.ncols(),
                                 "enumeration_dimension": basis.nrows() - 1})


def load_minima_input(checkpoint):
    from sage.all import ZZ, matrix
    doc = _read_document(Path(checkpoint) / "input.json")
    if doc["input_basis_sha256"] != _digest(doc["basis_rows"]):
        raise ValueError("Portable input basis hash mismatch")
    basis = matrix(ZZ, doc["basis_rows"])
    if (basis.nrows() != doc["lattice_dimension"] or basis.ncols() != doc["ambient_columns"] or
            doc["enumeration_dimension"] != basis.nrows() - 1):
        raise ValueError("Portable input dimensions mismatch")
    return basis, doc["instance"]


def verify_minima_bundle(checkpoint):
    """Validate a completed transferred bundle, using exact algebra, never LLL."""
    from sage.all import ZZ, matrix, vector
    folder = Path(checkpoint)
    bundle = _read_document(folder / "bundle.json")
    if not bundle.get("preparation_complete") or bundle["schema"] != "portable-minima-bundle-v1":
        raise ValueError("Incomplete or unsupported portable bundle")
    required = {"input.json", "pass1/input_basis.json", "pass1/ambient_lll.json",
                "pass1/quotient_lll.json", "pass1/prepared.json", "pass1/integer.prepared"}
    if not required.issubset(bundle["files"]):
        raise ValueError("Portable bundle omits required files")
    for name, expected in bundle["files"].items():
        if Path(name).is_absolute() or ".." in Path(name).parts:
            raise ValueError("Unsafe bundle file path")
        if _hash_file(folder / name) != expected:
            raise ValueError(f"Bundle file hash mismatch: {name}")
    original, _ = load_minima_input(folder)
    first = folder / "pass1"
    prepared = _read_document(first / "prepared.json")
    ambient = _read_document(first / "ambient_lll.json")
    quotient = _read_document(first / "quotient_lll.json")
    input_hash = _digest(_matrix_rows(original))
    if (prepared["schema"] != PREPARATION_VERSION or prepared["input_basis_sha256"] != input_hash or
            prepared["native_source_sha256"] != _hash_file(SOURCE)):
        raise ValueError("Portable engine or original basis identity mismatch")
    engine_hash = _digest([_hash_file(SOURCE), PREPARATION_VERSION, prepared["quotient_delta"]])
    if (prepared["engine_sha256"] != engine_hash or prepared["identity"] != _digest([input_hash, engine_hash]) or
            prepared["preparation_sha256"] != _hash_file(first / "integer.prepared")):
        raise ValueError("Portable mathematical preparation identity mismatch")
    frontier_identity = _digest(["exact-frontier-v2", input_hash, prepared["root"], prepared["tail"],
                                 prepared["inner"], prepared["content"], prepared["preparation_sha256"],
                                 prepared["native_source_sha256"]])
    if "frontier_identity" in prepared and prepared["frontier_identity"] != frontier_identity:
        raise ValueError("Portable frontier identity mismatch")
    for stage in ("input_basis", "ambient_lll", "quotient_lll"):
        if prepared[stage + "_file_sha256"] != _hash_file(first / (stage + ".json")):
            raise ValueError("Portable stage hash mismatch")
    if _read_document(first / "input_basis.json")["rows"] != _matrix_rows(original):
        raise ValueError("Portable original basis mismatch")
    warm = matrix(ZZ, ambient["lll_input_rows"])
    preconditioner = matrix(ZZ, ambient["preconditioner_rows"])
    transform = matrix(ZZ, ambient["transform_rows"])
    reduced = matrix(ZZ, ambient["reduced_rows"])
    if (ambient["input_basis_sha256"] != input_hash or ambient["backend"] != "fpLLL:proved" or
            ambient["delta"] != 0.8 or preconditioner * original != warm or transform * warm != reduced or
            abs(preconditioner.det()) != 1 or abs(transform.det()) != 1):
        raise ArithmeticError("Invalid portable ambient transformations")
    order = sorted(range(original.nrows()), key=lambda i: reduced.row(i).dot_product(reduced.row(i)))
    root = vector(ZZ, reduced.row(order[0]))
    root_sq = int(root.dot_product(root))
    tail = matrix(ZZ, [reduced.row(i) for i in order[1:]])
    projected = root_sq * tail - matrix(ZZ, [c * root for c in tail * root])
    coordinate_content = math.gcd(*(int(value) for value in projected.list()))
    projected = projected.apply_map(lambda value: value // coordinate_content)
    quotient_transform = matrix(ZZ, quotient["transform_rows"])
    if (quotient["ambient_file_sha256"] != _hash_file(first / "ambient_lll.json") or
            quotient["ambient_row_order"] != order or quotient["coordinate_content"] != str(coordinate_content) or
            quotient["backend"] != "fpLLL:proved" or quotient["delta"] != prepared["quotient_delta"] or
            matrix(ZZ, quotient["input_rows"]) != projected or
            quotient_transform * projected != matrix(ZZ, quotient["reduced_rows"]) or
            abs(quotient_transform.det()) != 1):
        raise ArithmeticError("Invalid portable quotient transformation")
    tail = quotient_transform * tail
    shifts = [_nearest(int(row.dot_product(root)), root_sq) for row in tail.rows()]
    tail = matrix(ZZ, [row - k * root for row, k in zip(tail.rows(), shifts)])
    inner = [int(value) for value in tail * root]
    gram = root_sq * (tail * tail.transpose()) - matrix(
        ZZ, len(inner), len(inner), lambda i, j: inner[i] * inner[j])
    content = math.gcd(*(int(value) for value in gram.list()))
    if (prepared["root"] != list(map(str, root)) or prepared["root_sq"] != str(root_sq) or
            prepared["tail"] != _matrix_rows(tail) or prepared["inner"] != list(map(str, inner)) or
            prepared["root_shifts"] != list(map(str, shifts)) or prepared["content"] != str(content)):
        raise ArithmeticError("Portable final lift basis mismatch")
    data = _decode_preparation(prepared)
    if _witness(data, prepared["initial_y"])[0] != int(prepared["initial_best"]):
        raise ArithmeticError("Invalid portable initial incumbent")
    with (first / "integer.prepared").open() as handle:
        header = handle.readline().split()
        native_inner = list(map(int, handle.readline().split()))
    if header != ["EXACT-GMP-PREP-1", str(original.nrows() - 1), str(root_sq), str(content)] or native_inner != inner:
        raise ValueError("Portable native preparation header mismatch")
    n = original.nrows()
    if (bundle["lattice_dimension"] != n or bundle["ambient_columns"] != original.ncols() or
            bundle["enumeration_dimension"] != n - 1 or bundle["quotient_lll_rows"] != n - 1):
        raise ValueError("Portable bundle dimensions mismatch")
    frontier = first / "frontier.sqlite3"
    if frontier.exists():
        conn = sqlite3.connect(frontier.resolve().as_uri() + "?mode=ro", uri=True)
        try:
            meta = dict(conn.execute("SELECT key,value FROM meta"))
            if meta.get("identity") != prepared.get("frontier_identity", prepared["payload_sha256"]):
                raise ValueError("Portable SQLite frontier identity mismatch")
            incumbent = json.loads(meta["best"])
            if _witness(data, incumbent["y"])[0] != int(incumbent["value"]):
                raise ArithmeticError("Portable frontier incumbent mismatch")
        finally:
            conn.close()
    return bundle


def prepare_successive_minima(basis, checkpoint, *, instance=None, workers=0,
                             quotient_delta=0.99, reuse=None):
    """Unbounded LLL/Bareiss preparation only; this is NOT a lambda2 certificate."""
    if workers < 0 or not 0.25 < quotient_delta < 1:
        raise ValueError("Invalid preparation parameters")
    basis = _validated_basis(basis)
    folder = Path(checkpoint).resolve()
    with _lock(folder):
        save_minima_input(basis, folder, instance)
        binary = _build_native(folder)
        first = folder / "pass1"
        first.mkdir(exist_ok=True)
        data = _make_preparation(basis, first, binary, quotient_delta=quotient_delta,
                                 preparation_threads=min(workers or _cpu_limit(), _cpu_limit()),
                                 reuse=Path(reuse).resolve() / "pass1" if reuse else None)
        result = {"schema": "portable-minima-bundle-v1", "preparation_complete": True,
                  "exact": False, "lambda2_sq": None, "lattice_dimension": basis.nrows(),
                  "ambient_columns": basis.ncols(), "quotient_lll_rows": basis.nrows() - 1,
                  "enumeration_dimension": basis.nrows() - 1,
                  "preparation_seconds": data["raw"]["preparation_seconds"],
                  "ambient_lll_seconds": data["raw"]["ambient_lll_seconds"],
                  "quotient_lll_seconds": data["raw"]["quotient_lll_seconds"],
                  "ambient_input_provenance": _read_document(first / "ambient_lll.json")["input_provenance"],
                  "legacy_native_reused": data["raw"]["legacy_native_reused"],
                  "files": {str(path.relative_to(folder)): _hash_file(path)
                            for path in [folder / "input.json", first / "input_basis.json",
                                         first / "ambient_lll.json", first / "quotient_lll.json",
                                         first / "prepared.json", first / "integer.prepared"]}}
        _save_document(folder / "bundle.json", result)
        return result


class _NativeWorker:
    def __init__(self, binary: Path, prepared: Path, register=None):
        self.proc = subprocess.Popen([str(binary), "serve", str(prepared)],
                                     stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, text=True, bufsize=1)
        try:
            if register is not None:
                register(self)
            header = self._line().split()
            if len(header) != 2 or header[0] != "READY":
                raise RuntimeError("Invalid native worker handshake")
            self.dimension = int(header[1])
        except BaseException:
            self.close()
            raise

    def _line(self) -> str:
        line = self.proc.stdout.readline()
        if not line:
            code = self.proc.wait()
            error = self.proc.stderr.read()
            raise RuntimeError(f"Native enumeration stopped ({code}): {error.strip()}")
        return line.rstrip("\n")

    def initial(self, best: int) -> list[str]:
        self.proc.stdin.write(f"INIT {best}\n")
        self.proc.stdin.flush()
        fields = self._line().split()
        if len(fields) != 2 or fields[0] != "FRAMES":
            raise RuntimeError("Invalid native frontier header")
        frames = [self._line() for _ in range(int(fields[1]))]
        if self._line() != "END":
            raise RuntimeError("Incomplete native frontier")
        return frames

    def chunk(self, frames: list[str], best: int, nodes: int, seconds: float) -> dict:
        self.proc.stdin.write(f"JOB {best} {nodes} {seconds} {len(frames)}\n")
        self.proc.stdin.write("\n".join(frames) + "\n")
        self.proc.stdin.flush()
        fields = self._line().split()
        if len(fields) != 6 or fields[0] != "RESULT":
            raise RuntimeError("Invalid native chunk header")
        value, visited, leaves, count, size = map(int, fields[1:])
        if visited < 0 or leaves < 0 or count < 0 or size not in (0, self.dimension):
            raise RuntimeError("Invalid native chunk statistics")
        y = list(map(int, self._line().split())) if size else []
        if len(y) != size:
            raise RuntimeError("Invalid native witness length")
        frontier = [self._line() for _ in range(count)]
        if self._line() != "END":
            raise RuntimeError("Incomplete native chunk")
        return {"best": value, "y": y, "nodes": visited,
                "leaves": leaves, "frames": frontier}

    def close(self):
        if self.proc.poll() is None:
            self.proc.terminate()
            try:
                self.proc.wait(timeout=3)
            except subprocess.TimeoutExpired:
                self.proc.kill()
                self.proc.wait()
        for stream in (self.proc.stdin, self.proc.stdout, self.proc.stderr):
            if not stream.closed:
                stream.close()


def _rss(pid: int) -> int:
    try:
        for line in Path(f"/proc/{pid}/status").read_text().splitlines():
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) * 1024
    except (OSError, ValueError):
        pass
    return 0


def _cpu_limit() -> int:
    try:
        cpus = len(os.sched_getaffinity(0))
    except AttributeError:
        cpus = os.cpu_count() or 1
    try:
        quota, period = Path("/sys/fs/cgroup/cpu.max").read_text().split()
        if quota != "max":
            cpus = min(cpus, max(1, int(quota) // int(period)))
    except (OSError, ValueError):
        pass
    return cpus


def _resource_plan(prepared: Path, workers: int, memory_gb: float) -> dict:
    cpus = _cpu_limit()
    available = 2 * GIB
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                available = int(line.split()[1]) * 1024
        limit = Path("/sys/fs/cgroup/memory.max").read_text().strip()
        if limit != "max":
            current = int(Path("/sys/fs/cgroup/memory.current").read_text())
            available = min(available, max(0, int(limit) - current))
    except (OSError, ValueError):
        pass
    overhead = _rss(os.getpid())
    budget = min(int(memory_gb * GIB) if memory_gb else int(available * 0.8),
                 overhead + int(available * 0.9))
    estimate = 16 * (1 << 20) + 5 * prepared.stat().st_size
    slots = (budget - overhead) // estimate
    if slots < 1:
        raise MemoryError("Memory budget cannot fit the coordinator and one estimated worker")
    count = max(1, min(workers or cpus, cpus, slots))
    return {"workers": count, "memory_budget_bytes": budget,
            "estimated_worker_bytes": estimate, "available_bytes": available,
            "estimate_is_not_a_hard_bound": True}


def _split_frontier(frames: list[str], target: int) -> list[str]:
    """Partition unvisited coefficient intervals without descending the tree.

    A paused SE frame contains two intervals around its already visited center.
    Splitting these intervals creates disjoint subtrees even at level zero,
    where returning ancestor frames alone would leave just one busy worker.
    """
    if len(frames) >= target:
        return frames
    pieces = []
    for frame in frames:
        fields = frame.split()
        low, high, left, right = map(int, fields[3:7])
        for start, end, descending in ((low, min(left, high), True),
                                        (max(right, low), high, False)):
            if start <= end:
                pieces.append((fields, start, end, descending))
    while len(pieces) < target:
        if not pieces:
            break
        index = max(range(len(pieces)), key=lambda i: pieces[i][2] - pieces[i][1])
        fields, start, end, descending = pieces[index]
        if start == end:
            break
        middle = (start + end) // 2
        pieces[index] = (fields, start, middle, descending)
        pieces.append((fields, middle + 1, end, descending))
    output = []
    for fields, start, end, descending in pieces:
        clone = fields.copy()
        # Only the unvisited interval is retained. Visit its endpoint nearest
        # the original center first, without requiring the unserialized pivot.
        left, right = (end, end + 1) if descending else (start - 1, start)
        clone[3:7] = list(map(str, (start, end, left, right)))
        output.append(" ".join(clone))
    return output


def _search(data: dict, folder: Path, binary: Path, *, workers: int, memory_gb: float,
            max_seconds: float, max_nodes: int, chunk_nodes: int,
            chunk_seconds: float, progress_seconds: float) -> dict:
    prepared = (folder / "integer.prepared").resolve()
    if data["raw"].get("native_source_sha256", _hash_file(SOURCE)) != _hash_file(SOURCE):
        raise ValueError("Native worker source version changed since preparation")
    if data["raw"]["schema"].endswith("-v1") and data["raw"]["binary_sha256"] != _hash_file(binary):
        raise ValueError("Compiled worker identity changed since preparation")
    plan = _resource_plan(prepared, workers, memory_gb)
    conn = sqlite3.connect(folder / "frontier.sqlite3")
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT NOT NULL)")
    conn.execute("CREATE TABLE IF NOT EXISTS jobs(id INTEGER PRIMARY KEY,payload TEXT NOT NULL,status TEXT NOT NULL)")
    identity = data["raw"].get("frontier_identity", data["raw"]["payload_sha256"])
    existing = conn.execute("SELECT value FROM meta WHERE key='identity'").fetchone()
    if existing:
        if existing[0] != identity:
            conn.close()
            raise ValueError("Frontier identity mismatch")
        with conn:
            conn.execute("UPDATE jobs SET status='pending' WHERE status='running'")
    else:
        best = int(data["raw"]["initial_best"])
        worker = _NativeWorker(binary, prepared)
        try:
            frames = _split_frontier(worker.initial(best), 2 * plan["workers"])
        finally:
            worker.close()
        with conn:
            for key, value in (("identity", identity),
                               ("best", json.dumps({"value": str(best), "y": data["raw"]["initial_y"]})),
                               ("nodes", "0"), ("leaves", "0"), ("elapsed", "0")):
                conn.execute("INSERT INTO meta VALUES(?,?)", (key, value))
            for frame in frames:
                conn.execute("INSERT INTO jobs(payload,status) VALUES(?,'pending')", (json.dumps([frame]),))
    incumbent = json.loads(conn.execute("SELECT value FROM meta WHERE key='best'").fetchone()[0])
    best, best_y = int(incumbent["value"]), incumbent["y"]
    if _witness(data, best_y)[0] != best:
        conn.close()
        raise ArithmeticError("Checkpoint incumbent does not match its witness")
    nodes = int(conn.execute("SELECT value FROM meta WHERE key='nodes'").fetchone()[0])
    leaves = int(conn.execute("SELECT value FROM meta WHERE key='leaves'").fetchone()[0])
    elapsed_before = float(conn.execute("SELECT value FROM meta WHERE key='elapsed'").fetchone()[0])
    node_start = nodes
    began = last_progress = time.monotonic()
    reason = "complete"
    stopping = False
    local = threading.local()
    native_workers = []
    registry_lock = threading.Lock()
    abort = threading.Event()

    def register(worker):
        with registry_lock:
            native_workers.append(worker)
            if abort.is_set():
                worker.proc.terminate()

    def work(frames, bound):
        if not hasattr(local, "worker"):
            local.worker = _NativeWorker(binary, prepared, register=register)
        return local.worker.chunk(frames, bound, chunk_nodes, chunk_seconds)

    def save_meta():
        values = {"best": json.dumps({"value": str(best), "y": best_y}),
                  "nodes": str(nodes), "leaves": str(leaves),
                  "elapsed": str(elapsed_before + time.monotonic() - began)}
        for key, value in values.items():
            conn.execute("UPDATE meta SET value=? WHERE key=?", (value, key))

    print(f"Exact GMP search: workers={plan['workers']}, memory budget="
          f"{plan['memory_budget_bytes'] / GIB:.2f} GiB", flush=True)
    print("EXACT_STAGE " + json.dumps({"stage": "enumeration", "dimension": len(data["inner"]),
                                        "workers": plan["workers"]}), flush=True)
    pool = ThreadPoolExecutor(max_workers=plan["workers"])
    futures = {}
    try:
        while True:
            elapsed = time.monotonic() - began
            rss = _rss(os.getpid()) + sum(_rss(w.proc.pid) for w in native_workers)
            if max_seconds and elapsed >= max_seconds:
                stopping, reason = True, "time_limit"
            if max_nodes and nodes - node_start >= max_nodes:
                stopping, reason = True, "node_limit"
            if rss > plan["memory_budget_bytes"]:
                stopping, reason = True, "memory_limit"
            if not stopping:
                pending = conn.execute("SELECT id,payload FROM jobs WHERE status='pending' ORDER BY id LIMIT ?",
                                       (plan["workers"] - len(futures),)).fetchall()
                for job_id, payload in pending:
                    with conn:
                        conn.execute("UPDATE jobs SET status='running' WHERE id=?", (job_id,))
                    futures[pool.submit(work, json.loads(payload), best)] = job_id
            if not futures:
                remaining = conn.execute("SELECT count(*) FROM jobs WHERE status!='done'").fetchone()[0]
                if remaining == 0:
                    reason = "complete"
                    break
                if stopping:
                    break
            done, _ = wait(futures, timeout=0.2, return_when=FIRST_COMPLETED)
            for future in done:
                job_id = futures.pop(future)
                result = future.result()
                if result["y"]:
                    value, _ = _witness(data, result["y"])
                    if value != result["best"]:
                        raise ArithmeticError("Worker incumbent failed exact witness verification")
                    if value < best:
                        best, best_y = value, result["y"]
                elif result["best"] < best:
                    raise ArithmeticError("Worker decreased its bound without a witness")
                nodes += result["nodes"]
                leaves += result["leaves"]
                # Replace consumed work atomically by its disjoint unfinished frames.
                frames = _split_frontier(result["frames"], 2 * plan["workers"])
                groups = min(len(frames), max(4, 2 * plan["workers"]))
                with conn:
                    conn.execute("UPDATE jobs SET status='done' WHERE id=?", (job_id,))
                    for group in range(groups):
                        conn.execute("INSERT INTO jobs(payload,status) VALUES(?,'pending')",
                                     (json.dumps(frames[group::groups]),))
                    save_meta()
            if time.monotonic() - last_progress >= progress_seconds:
                remaining = conn.execute("SELECT count(*) FROM jobs WHERE status!='done'").fetchone()[0]
                print(f"nodes={nodes:,} frontier={remaining} "
                      f"candidate_log2={math.log2(best) / 2:.6f} RSS={rss / GIB:.2f} GiB", flush=True)
                last_progress = time.monotonic()
    except KeyboardInterrupt:
        reason = "interrupted"
        abort.set()
        for worker in native_workers:
            worker.proc.terminate()
    except BaseException:
        abort.set()
        for worker in native_workers:
            worker.proc.terminate()
        with conn:
            save_meta()
            conn.execute("UPDATE jobs SET status='pending' WHERE status='running'")
        raise
    finally:
        pool.shutdown(wait=True, cancel_futures=True)
        for worker in native_workers:
            worker.close()
        with conn:
            save_meta()
            conn.execute("UPDATE jobs SET status='pending' WHERE status='running'")
        remaining = conn.execute("SELECT count(*) FROM jobs WHERE status!='done'").fetchone()[0]
        conn.close()
    best_sq, witness = _witness(data, best_y)
    if best_sq != best:
        raise ArithmeticError("Final incumbent failed exact verification")
    result = {"complete": remaining == 0 and reason == "complete", "reason": reason,
              "off_line_sq": str(best), "off_line_witness": list(map(str, witness)),
              "quotient_coefficients": list(map(str, best_y)), "nodes": nodes,
              "leaves": leaves, "frontier_jobs": remaining,
              "search_seconds": elapsed_before + time.monotonic() - began,
              "invocation_search_seconds": time.monotonic() - began,
              "invocation_nodes": nodes - node_start,
              "resource_plan": plan, "prepared_sha256": identity}
    _atomic_json(folder / "search_result.json", result)
    return result


def exact_successive_minima(basis, checkpoint, *, workers=0, memory_gb=0,
                           max_seconds=0, max_nodes=0, chunk_nodes=20000,
                           chunk_seconds=2.0, progress_seconds=10.0,
                           quotient_delta=0.99) -> dict:
    """Return certified minima, or a resumable incomplete result.

    Resource/time limits apply to search; basis construction and preparation
    are separately reported. ``workers=0`` uses the CPU/RAM-aware maximum.
    Never infer a shortest root from the first reduced row. If the first pass
    finds a shorter off-line vector, rerun around that certified shortest vector.
    """
    from sage.all import ZZ, matrix

    if workers < 0 or max_nodes < 0 or chunk_nodes < 1:
        raise ValueError("Invalid worker or node limit")
    if any(not math.isfinite(value) or value < 0
           for value in (memory_gb, max_seconds)) or any(
               not math.isfinite(value) or value <= 0
               for value in (chunk_seconds, progress_seconds)):
        raise ValueError("Invalid resource or timing limit")
    if not 0.25 < quotient_delta < 1:
        raise ValueError("quotient_delta must lie strictly between 0.25 and 1")
    basis = _validated_basis(basis)
    folder = Path(checkpoint).resolve()
    with _lock(folder):
        binary = _build_native(folder)
        first = folder / "pass1"
        first.mkdir(exist_ok=True)
        preparation_threads = min(workers or _cpu_limit(), _cpu_limit())
        data = _make_preparation(basis, first, binary, quotient_delta=quotient_delta,
                                 preparation_threads=preparation_threads)
        options = dict(workers=workers, memory_gb=memory_gb, max_seconds=max_seconds,
                       max_nodes=max_nodes, chunk_nodes=chunk_nodes,
                       chunk_seconds=chunk_seconds, progress_seconds=progress_seconds)
        search = _search(data, first, binary, **options)
        passes = [search]
        root_sq = data["root_sq"]
        if not search["complete"]:
            result = {"exact": False, "status": search["reason"],
                      "lambda1_sq": None, "lambda2_sq": None, "passes": passes}
        elif int(search["off_line_sq"]) >= root_sq:
            result = {"exact": True, "status": "complete", "lambda1_sq": str(root_sq),
                      "lambda2_sq": search["off_line_sq"],
                      "lambda1_witness": list(map(str, data["root"])),
                      "lambda2_witness": search["off_line_witness"], "passes": passes}
        else:
            # The first off-line minimizer is a certified ambient shortest vector.
            y = list(map(int, search["quotient_coefficients"]))
            dot = sum(c * z for c, z in zip(data["inner"], y))
            coefficients = [-_nearest(dot, root_sq)] + y
            if math.gcd(*coefficients) != 1:
                raise ArithmeticError("Certified shortest vector is not primitive in the lattice")
            smith, left, right = matrix(ZZ, [coefficients]).smith_form()
            extension = right.inverse().change_ring(ZZ)
            extension[0] *= ZZ(left[0, 0] * smith[0, 0])
            if list(extension.row(0)) != coefficients or abs(extension.det()) != 1:
                raise ArithmeticError("Could not extend the shortest vector to a lattice basis")
            complete_basis = matrix(ZZ, [data["root"]] + data["tail"])
            second_basis = extension * complete_basis
            second = folder / "pass2"
            second.mkdir(exist_ok=True)
            second_data = _make_preparation(second_basis, second, binary, quotient_delta=quotient_delta,
                                            preparation_threads=preparation_threads)
            if second_data["root_sq"] != int(search["off_line_sq"]):
                raise ArithmeticError("Second pass lost the certified shortest vector")
            if max_seconds:
                options["max_seconds"] = max(0.001, max_seconds - search["invocation_search_seconds"])
            if max_nodes:
                options["max_nodes"] = max(1, max_nodes - search["invocation_nodes"])
            second_search = _search(second_data, second, binary, **options)
            passes.append(second_search)
            exact = second_search["complete"]
            result = {"exact": exact, "status": "complete" if exact else second_search["reason"],
                      "lambda1_sq": search["off_line_sq"],
                      "lambda1_witness": search["off_line_witness"],
                      "lambda2_sq": second_search["off_line_sq"] if exact else None,
                      "lambda2_witness": second_search["off_line_witness"] if exact else None,
                      "passes": passes}
            if exact and int(result["lambda2_sq"]) < int(result["lambda1_sq"]):
                raise ArithmeticError("Second minimum is smaller than the certified first minimum")
        result.update(schema="exact-successive-minima-v1", reduction_backend="fpLLL:proved",
                      enumeration_backend="C++ GMP integer Bareiss DFS",
                      checkpoint=str(folder), native_source_sha256=_hash_file(SOURCE),
                      native_binary_sha256=_hash_file(binary),
                      preparation_seconds=sum(json.loads((folder / f"pass{i+1}" / "prepared.json").read_text())
                                              ["preparation_seconds"] for i in range(len(passes))))
        if result["exact"]:
            result["log2_lambda1"] = math.log2(int(result["lambda1_sq"])) / 2
            result["log2_lambda2"] = math.log2(int(result["lambda2_sq"])) / 2
        _atomic_json(folder / "result.json", result)
        return result
