"""True SVP lengths of the full and first-jet root-orthogonal lattices.

The original experiments and their certificates are read-only. New geometry,
integer embeddings and exact-enumeration checkpoints live in projected_minima/.
For an integer projected basis P and positive rational scale s, embed
Z e_0 orthogonally to P. Its second minimum is exactly lambda_1(P), because
all nonzero integer vectors have squared length >= 1. The existing certified
successive-minimum engine therefore certifies the requested projected SVP.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import time

import mpmath as mp
from experiment import (
    HERE, GRID, build, variant, name, rows, norm2, matrix, vector, ZZ, QQ,
    diagonal_matrix, gcd,
)
from shared.dual_moment import snf_congruence_kernel
from shared import exact_minima as engine
from certify import proved_lll_xd
from coefficient_certificate import restrict

OUT = HERE / "projected_minima"
NATIVE_SEARCH = engine._search
engine._proved_lll = proved_lll_xd
mp.mp.dps = 80


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp-" + str(os.getpid()))
    temporary.write_text(json.dumps(value, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def log2_int(value):
    return mp.log(int(value), 2)


def prepare_geometry(index):
    cfg = GRID[index]
    folder = OUT / name(cfg)
    target = folder / "geometry.json"
    if target.exists():
        return read(target)
    began = time.monotonic()
    data = build(cfg)
    n = len(data["keys"])
    ell = cfg[3]
    _, exponents, scalar, root = variant(data, "B")
    Q = diagonal_matrix(ZZ, [data["N"] ** (cfg[2] - e) for e in exponents]) * data["forward"] * data["jet"].transpose()
    coefficients = snf_congruence_kernel(Q, data["N"] ** cfg[2])
    assert coefficients.nrows() == coefficients.ncols() == ell + 1
    assert all(z % data["N"] ** cfg[2] == 0 for z in (coefficients * Q.transpose()).list())
    balanced = data["jet"] * diagonal_matrix(ZZ, data["scales"])
    jet = coefficients * balanced
    coeff_index = abs(coefficients.det())
    gram_det = (balanced * balanced.transpose()).det()
    jet_gram_det = (jet * jet.transpose()).det()
    assert jet_gram_det == coeff_index ** 2 * gram_det > 0
    root_coeff = (vector(QQ, [scalar] + [0] * ell) * coefficients.inverse()).change_ring(ZZ)
    assert root_coeff * jet == root and gcd(root_coeff.list()) == 1
    smith, left, right = matrix(ZZ, [root_coeff]).smith_form()
    extension = right.inverse().change_ring(ZZ)
    extension[0] *= left[0, 0] * smith[0, 0]
    assert extension.row(0) == root_coeff and abs(extension.det()) == 1
    jet = extension * jet
    assert jet.row(0) == root
    root_sq = norm2(root)
    determinant_N_exponent = sum(exponents)
    determinant_R_exponent = sum(cfg[1] - sum(k) for k in data["keys"])
    log_det = determinant_N_exponent * log2_int(data["N"]) + determinant_R_exponent * log2_int(data["R"])
    log_root = log2_int(root_sq) / 2
    log_jet_covol = log2_int(jet_gram_det) / 2

    def ball_log(d):
        return d * mp.log(mp.pi, 2) / 2 - mp.loggamma(mp.mpf(d) / 2 + 1) / mp.log(2)

    gh = (log_det - log_root - ball_log(n - 1)) / (n - 1)
    gh_asym = mp.log((n - 1) / (2 * mp.pi * mp.e), 2) / 2 + (log_det - log_root) / (n - 1)
    jet_pred = mp.log(ell / (2 * mp.pi * mp.e), 2) / 2 + (log_jet_covol - log_root) / ell
    jet_ball = (log_jet_covol - log_root - ball_log(ell)) / ell
    geometry = dict(
        schema="root-projected-geometry-v1", index=index, case=name(cfg),
        application=cfg[0].upper(), t=cfg[1], m=cfg[2], ell=ell, ambient_dimension=n,
        full_projected_rank=n - 1, jet_projected_rank=ell,
        root=rows(matrix(ZZ, [root]))[0], root_sq=str(root_sq),
        jet_basis_rows=rows(jet), jet_coefficient_rows=rows(coefficients),
        jet_gram_determinant=str(jet_gram_det),
        log2_lambda1_root=str(log_root), log2_gh_prediction=str(gh),
        log2_gh_asymptotic_prediction=str(gh_asym),
        log2_jet_prediction=str(jet_pred), log2_jet_exact_ball_prediction=str(jet_ball),
        log2_full_projected_covolume=str(log_det - log_root),
        log2_jet_projected_covolume=str(log_jet_covol - log_root),
        seconds=time.monotonic() - began,
    )
    save(target, geometry)
    print("GEOMETRY", index, name(cfg), "seconds", round(geometry["seconds"], 2), flush=True)
    return geometry


def embedding_from_basis(basis, root, original_scale=ZZ(1)):
    assert basis.row(0) == root
    common = gcd(basis.list())
    assert common > 0
    basis = basis.apply_map(lambda z: z // common)
    root = basis.row(0)
    a = norm2(root)
    tail = basis[1:]
    inner = tail * root
    projected = a * tail - matrix(ZZ, [dot * root for dot in inner])
    content = gcd(projected.list())
    assert content > 0
    projected = projected.apply_map(lambda z: z // content)
    scale = QQ(original_scale * common * content) / a
    embedded = matrix(ZZ, [[1] + [0] * projected.ncols()] + [[0] + list(row) for row in projected.rows()])
    assert embedded.nrows() == basis.nrows()
    return embedded, scale, dict(common_divisor=str(common), scaled_root_sq=str(a), projection_coordinate_content=str(content))


def prepare_embedding(index, kind):
    geometry = prepare_geometry(index)
    folder = OUT / geometry["case"] / kind
    target = folder / "embedding.json"
    if target.exists():
        doc = read(target)
        return matrix(ZZ, doc["rows"]), QQ(doc["scale_numerator"]) / ZZ(doc["scale_denominator"]), doc
    if kind == "jet":
        basis = matrix(ZZ, geometry["jet_basis_rows"])
        root = vector(ZZ, geometry["root"])
        embedded, scale, details = embedding_from_basis(basis, root)
    else:
        source = HERE / ("flag_results" if GRID[index][0] == "hm" else "results") / geometry["case"] / "B.json"
        doc = read(source)
        original = matrix(ZZ, doc["reduced_rows"])
        assert doc["b1_on_root"] and norm2(original.row(0)) == ZZ(geometry["root_sq"])
        common = gcd(original.list())
        basis = original.apply_map(lambda z: z // common)
        root = basis.row(0)
        a = norm2(root)
        projected_numerators = [a * norm2(row) - row.dot_product(root) ** 2 for row in basis.rows()[1:]]
        upper_projected_num = min(projected_numerators)
        # Any projection this short has a nearest lift of norm squared at most
        # upper_projected_num/a + a/4. The exact inverse certificate excludes
        # coefficients for ALL such lifts, not only the original row vectors.
        lift_bound = (QQ(upper_projected_num) / a + QQ(a) / 4).ceil()
        restricted, proof = restrict(basis, lift_bound)
        assert proof["active_rows"][0] == 0 and restricted.row(0) == root
        save(folder / "coefficient_exclusion.json", proof)
        embedded, scale, details = embedding_from_basis(restricted, root, common)
        details.update(source=str(source.relative_to(HERE)), source_sha256=engine._hash_file(source),
                       full_rank=basis.nrows(), restricted_rank=restricted.nrows(),
                       projection_candidate_numerator=str(upper_projected_num), lift_bound_sq=str(lift_bound))
    result = dict(schema="integer-orthogonal-svp-embedding-v1", case=geometry["case"], kind=kind,
                  index=index, rows=rows(embedded), scale_numerator=str(scale.numerator()),
                  scale_denominator=str(scale.denominator()), projected_rank=embedded.nrows() - 1,
                  ambient_columns=embedded.ncols(), **details)
    save(target, result)
    print("EMBEDDING", index, kind, "search_dimension", result["projected_rank"], flush=True)
    return embedded, scale, result


def summarize():
    records = []
    for index, cfg in enumerate(GRID):
        folder = OUT / name(cfg)
        if not (folder / "geometry.json").exists():
            continue
        geo = read(folder / "geometry.json")
        record = {key: geo[key] for key in ("index", "case", "application", "t", "m", "ell", "ambient_dimension",
                                              "log2_gh_prediction", "log2_jet_prediction")}
        for kind in ("full", "jet"):
            paths = [folder / kind / "minimum.json"]
            if kind == "full":
                paths.append(folder / "full_strong" / "minimum.json")
                paths.extend(path / "minimum.json" for path in sorted(folder.glob("full_sieve_*")))
            available = [(path, read(path)) for path in paths if path.exists()]
            if available:
                path, minimum = min(available, key=lambda item: (not item[1].get("exact"), float(item[1].get("log2_upper_bound", "inf"))))
            else:
                path, minimum = paths[0], {}
            record[f"{kind}_exact"] = minimum.get("exact", False)
            record[f"{kind}_log2_minimum"] = minimum.get("log2_length", "") if minimum.get("exact") else ""
            record[f"{kind}_log2_upper_bound"] = minimum.get("log2_upper_bound", "")
            record[f"{kind}_squared_norm_numerator"] = minimum.get("squared_norm_numerator", "") if minimum.get("exact") else ""
            record[f"{kind}_squared_norm_denominator"] = minimum.get("squared_norm_denominator", "") if minimum.get("exact") else ""
            record[f"{kind}_result_source"] = str(path.relative_to(HERE)) if minimum else ""
        bound_path = folder / "full_bounds.json"
        interval = read(bound_path) if bound_path.exists() else {}
        record["full_value_certified_6dp"] = bool(record["full_exact"] or interval.get("certified_decimal_places", 0) >= 6)
        record["full_log2_display_6dp"] = (format(float(record["full_log2_minimum"]), ".6f") if record["full_exact"]
                                             else interval.get("log2_rounded", ""))
        record["full_log2_lower_bound"] = interval.get("log2_lower_bound", "")
        record["full_log2_upper_bound_from_certificate"] = interval.get("log2_upper_bound", "")
        records.append(record)
    if records:
        OUT.mkdir(exist_ok=True)
        with (OUT / "comparison.csv").open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
        save(OUT / "status.json", dict(inputs=len(records), full_exact=sum(r["full_exact"] for r in records),
                                       full_display_certified_6dp=sum(r["full_value_certified_6dp"] for r in records),
                                       jet_exact=sum(r["jet_exact"] for r in records), records=records))
    return records


def run(index, kind, seconds, workers, backend):
    basis, scale, embedding = prepare_embedding(index, kind)
    folder = OUT / embedding["case"] / kind
    result_path = folder / "minimum.json"
    if result_path.exists() and read(result_path).get("exact"):
        print("MINIMUM_REUSED", index, kind, read(result_path)["log2_length"], flush=True)
        return
    if backend == "partition":
        os.environ.setdefault("PARTITION_CHUNK_SECONDS", "10")
        from partition_search import search
    elif backend == "fast":
        from fast_search import search
    else:
        search = NATIVE_SEARCH
    engine._search = search
    result = engine.exact_successive_minima(basis, folder / "certificate", workers=workers, memory_gb=1,
        max_seconds=seconds, chunk_nodes=20000, chunk_seconds=1, progress_seconds=15)
    assert result.get("lambda1_sq") in (None, "1")
    exact = bool(result["exact"])
    best_sq = result.get("lambda2_sq")
    witness = result.get("lambda2_witness")
    if best_sq is None:
        run_result = result["passes"][0]
        best_sq = run_result["off_line_sq"]
        witness = run_result["off_line_witness"]
    witness = vector(ZZ, witness)
    assert witness[0] == 0 and norm2(witness) == ZZ(best_sq)
    coefficients = basis.solve_left(matrix(ZZ, [witness])).change_ring(ZZ)
    assert coefficients * basis == matrix(ZZ, [witness])
    squared = QQ(best_sq) * scale ** 2
    length_log = (log2_int(squared.numerator()) - log2_int(squared.denominator())) / 2
    minimum = dict(case=embedding["case"], kind=kind, exact=exact,
                   squared_norm_numerator=str(squared.numerator()), squared_norm_denominator=str(squared.denominator()),
                   log2_length=str(length_log) if exact else None, log2_upper_bound=str(length_log),
                   integer_embedding_witness=list(map(str, witness)),
                   projected_witness_numerator=[str(z * scale.numerator()) for z in witness[1:]],
                   projected_witness_denominator=str(scale.denominator()),
                   certificate=str((folder / "certificate").relative_to(HERE)),
                   embedding_sha256=engine._hash_file(folder / "embedding.json"), backend=backend,
                   reason=result["status"])
    save(result_path, minimum)
    summarize()
    print("PROJECTED_MINIMUM", index, kind, "exact", exact, "log2", str(length_log), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "run", "status"))
    parser.add_argument("--kind", choices=("full", "jet"), default="jet")
    parser.add_argument("--indices", default="all")
    parser.add_argument("--seconds", type=float, default=120)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--backend", choices=("native", "fast", "partition"), default="native")
    args = parser.parse_args()
    if args.mode == "status":
        records = summarize()
        print("COUNTS", len(records), "full", sum(r["full_exact"] for r in records), "jet", sum(r["jet_exact"] for r in records))
        return
    indices = range(len(GRID)) if args.indices == "all" else [int(i) for i in args.indices.split(",")]
    for index in indices:
        if args.mode == "prepare":
            prepare_geometry(index)
            summarize()
        else:
            run(index, args.kind, args.seconds, args.workers, args.backend)
    summarize()


if __name__ == "__main__":
    main()
