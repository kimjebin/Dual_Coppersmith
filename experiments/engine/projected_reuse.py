"""Reuse strong quotient bases and their exact Bareiss data for projected SVP.

The normalized integer Gram is unchanged by the orthogonal unit-vector
embedding. Existing Bareiss coefficients can therefore be copied, while the
root/scale header, identities and search frontiers are new. The original
proved quotient-LLL output is linked by exact algebra and source hashes.
No search completion is imported from the original nearest-lift objective.
"""
import argparse
import math
import os
import shutil
import signal
import time
from pathlib import Path
from datetime import datetime, timezone, timedelta
import subprocess
import sys

from sage.all import identity_matrix
from projected_minima import (OUT, HERE, GRID, read, save, prepare_geometry, rows,
    norm2, matrix, vector, ZZ, QQ, gcd, log2_int, summarize)
from certify import proved_lll_xd
from shared import exact_minima as engine

SOURCES = {0: "B-qbkz65s", 1: "B-qbkz77s", 2: "B-qbkz80s", 3: "B-qbkz80s",
           7: "B-bkz40", 8: "B-qbkz80s", 9: "B-flatter-qbkz80s"}
engine._proved_lll = proved_lll_xd


def prepare(index):
    geo = prepare_geometry(index)
    folder = OUT / geo["case"] / "full_strong"
    checkpoint = folder / "certificate"
    stage = checkpoint / "pass1"
    embedding_path = folder / "embedding.json"
    if (stage / "prepared.json").exists():
        embedding = read(embedding_path)
        return matrix(ZZ, embedding["rows"]), QQ(embedding["scale_numerator"]) / ZZ(embedding["scale_denominator"]), embedding
    began = time.monotonic()
    source = HERE / "certificates" / geo["case"] / SOURCES[index]
    source_stage = source / "pass1"
    old = engine._read_document(source_stage / "prepared.json")
    source_input = engine._read_document(source_stage / "input_basis.json")
    source_ambient = engine._read_document(source_stage / "ambient_lll.json")
    source_quotient = engine._read_document(source_stage / "quotient_lll.json")
    assert old["native_source_sha256"] == engine._hash_file(engine.SOURCE)
    assert old["preparation_sha256"] == engine._hash_file(source_stage / "integer.prepared")
    d = len(old["tail"])
    assert d + 1 == geo["ambient_dimension"]
    proof = read(source / "coefficient_exclusion.json")
    assert proof["active_rows"] == list(range(d + 1))
    measurement_folder = "flatter_results" if index == 9 else "results"
    measurement_path = HERE / measurement_folder / geo["case"] / "B.json"
    measured = read(measurement_path)
    B = matrix(ZZ, measured["reduced_rows"])
    physical_common = gcd(B.list())
    expected = B.apply_map(lambda z: z // physical_common)
    if (source / "preconditioning.json").exists():
        W = matrix(ZZ, read(source / "preconditioning.json")["transform_rows"])
        assert abs(W.det()) == 1
        expected = W * expected
    assert source_input["rows"] == rows(expected)
    U = matrix(ZZ, source_ambient["transform_rows"])
    warm = matrix(ZZ, source_ambient["lll_input_rows"])
    pre = matrix(ZZ, source_ambient["preconditioner_rows"])
    red = matrix(ZZ, source_ambient["reduced_rows"])
    assert pre * expected == warm and U * warm == red
    assert abs(pre.det()) == abs(U.det()) == 1
    order = sorted(range(d + 1), key=lambda i: norm2(red.row(i)))
    root = red.row(order[0])
    assert physical_common * root in (vector(ZZ, geo["root"]), -vector(ZZ, geo["root"]))
    a = norm2(root)
    T0 = matrix(ZZ, [red.row(i) for i in order[1:]])
    P0 = a * T0 - matrix(ZZ, [dot * root for dot in T0 * root])
    coordinate = gcd(P0.list())
    V = matrix(ZZ, source_quotient["transform_rows"])
    assert abs(V.det()) == 1 and V * P0.apply_map(lambda z: z // coordinate) == matrix(ZZ, source_quotient["reduced_rows"])
    T = matrix(ZZ, old["tail"])
    inner = vector(ZZ, old["inner"])
    assert vector(ZZ, old["root"]) == root and str(a) == old["root_sq"]
    assert T * root == inner
    assert T == V * T0 - matrix(ZZ, [ZZ(k) * root for k in old["root_shifts"]])
    P = a * T - matrix(ZZ, [dot * root for dot in inner])
    h = gcd(P.list())
    P = P.apply_map(lambda z: z // h)
    assert P == matrix(ZZ, source_quotient["reduced_rows"])
    scale = QQ(physical_common * h) / a
    embedded = matrix(ZZ, [[1] + [0] * P.ncols()] + [[0] + list(row) for row in P.rows()])
    old_content = ZZ(old["content"])
    gamma = QQ(a * old_content) / h ** 2
    assert gamma.denominator() == 1 and gamma > 0
    gamma = ZZ(gamma)
    H = (a * T * T.transpose() - matrix(ZZ, d, d, lambda i, j: inner[i] * inner[j])).apply_map(lambda z: z // old_content)
    assert P * P.transpose() == gamma * H and gcd(H.list()) == 1
    best_y = list(map(int, old["initial_y"]))
    old_result_path = source / "result.json"
    if old_result_path.exists():
        old_result = read(old_result_path)
        if old_result.get("passes"):
            candidate = list(map(int, old_result["passes"][0]["quotient_coefficients"]))
            if norm2(vector(ZZ, candidate) * P) < norm2(vector(ZZ, best_y) * P):
                best_y = candidate
    best = norm2(vector(ZZ, best_y) * P)
    assert best > 0
    stage.mkdir(parents=True, exist_ok=True)
    with (source_stage / "integer.prepared").open() as inp, (stage / "integer.prepared").open("w") as out:
        assert inp.readline().split() == ["EXACT-GMP-PREP-1", str(d), str(a), str(old_content)]
        assert inp.readline().split() == list(map(str, inner))
        out.write(f"EXACT-GMP-PREP-1 {d} 1 {gamma}\n")
        out.write(" ".join(["0"] * d) + "\n")
        shutil.copyfileobj(inp, out)
    input_rows = rows(embedded)
    input_hash = engine._digest(input_rows)
    synthetic_root = ["1"] + ["0"] * P.ncols()
    tail_rows = [["0"] + list(map(str, row)) for row in P.rows()]
    doc = dict(schema="projected-svp-reused-gram-v1", input_basis_sha256=input_hash,
        native_source_sha256=engine._hash_file(engine.SOURCE), root=synthetic_root, root_sq="1", tail=tail_rows,
        inner=["0"] * d, content=str(gamma), root_shifts=["0"] * d, initial_best=str(best), initial_y=best_y,
        preparation_sha256=engine._hash_file(stage / "integer.prepared"), preparation_seconds=time.monotonic() - began,
        source_preparation_sha256=engine._hash_file(source_stage / "prepared.json"),
        source_integer_preparation_sha256=engine._hash_file(source_stage / "integer.prepared"),
        normalized_integer_gram_unchanged=True,
        lll_provenance="Exact scalar multiple of original proved-LLL quotient output; verified again by integer LLL precondition checker")
    embedding = dict(case=geo["case"], index=index, kind="full", rows=input_rows, projected_rank=d,
        scale_numerator=str(scale.numerator()), scale_denominator=str(scale.denominator()),
        original_measurement=str(measurement_path.relative_to(HERE)), original_measurement_sha256=engine._hash_file(measurement_path),
        source_preparation=str(source_stage.relative_to(HERE)), source_preparation_sha256=engine._hash_file(source_stage / "prepared.json"),
        native_gram_identity_verified=True, complete_original_lattice_basis_verified=True,
        construction="Full-rank original lattice quotient; all original coefficient directions retained; exact scaled normalized Gram identity")
    save(embedding_path, embedding)
    engine._save_document(stage / "prepared.json", doc)
    print("STRONG_PROJECTED_PREPARED", index, SOURCES[index], "dimension", d, "seconds", round(time.monotonic() - began, 1), flush=True)
    return embedded, scale, embedding


def run(index, seconds, workers):
    basis, scale, embedding = prepare(index)
    folder = OUT / embedding["case"] / "full_strong"
    if (folder / "minimum.json").exists() and read(folder / "minimum.json").get("exact"):
        return folder
    candidate_path = OUT / "sieving" / embedding["case"] / "best_candidate.json"
    if candidate_path.exists():
        base_stage = folder / "certificate" / "pass1"
        base = engine._read_document(base_stage / "prepared.json")
        candidate = read(candidate_path)
        if candidate["prepared_sha256"] == base["payload_sha256"]:
            y = list(map(int, candidate["quotient_coefficients"]))
            value, _ = engine._witness(engine._decode_preparation(base), y)
            assert value > 0
            d = len(base["tail"])
            gain = (math.log2(int(base["initial_best"])) - math.log2(value)) / 2
            existing_folder=folder
            existing_bound=int(base['initial_best'])
            for previous_folder in (OUT/embedding['case']).glob('full_sieve_*'):
                previous_path=previous_folder/'certificate'/'pass1'/'prepared.json'
                if not previous_path.exists():
                    continue
                previous=engine._read_document(previous_path)
                if (previous.get('previous_prepared_sha256')==base['payload_sha256']
                        and previous.get('input_basis_sha256')==base['input_basis_sha256']
                        and previous.get('preparation_sha256')==base['preparation_sha256']
                        and int(previous['initial_best'])<existing_bound):
                    existing_folder=previous_folder;existing_bound=int(previous['initial_best'])
            incremental_volume_gain=d*(math.log2(existing_bound)-math.log2(value))/2
            # A radius-volume comparison chooses whether to start a new sphere;
            # it is only a cost model, never a shortestness assertion.
            if incremental_volume_gain < 0.5:
                folder=existing_folder
                embedding=read(folder/'embedding.json')
            elif gain * d >= 0.5:
                identity = engine._digest([base["payload_sha256"], str(value)])[:12]
                improved = OUT / embedding["case"] / ("full_sieve_" + identity)
                new_stage = improved / "certificate" / "pass1"
                if not (new_stage / "prepared.json").exists():
                    new_stage.mkdir(parents=True, exist_ok=True)
                    target_native = new_stage / "integer.prepared"
                    if not target_native.exists():
                        try:
                            os.link(base_stage / "integer.prepared", target_native)
                        except OSError:
                            shutil.copyfile(base_stage / "integer.prepared", target_native)
                    save(improved / "initial_candidate.json", candidate)
                    embedding_new = dict(embedding)
                    embedding_new.update(radius_source="exactly_checked_quick_svp_candidate",
                        base_embedding_sha256=engine._hash_file(folder / "embedding.json"),
                        initial_candidate_sha256=engine._hash_file(improved / "initial_candidate.json"),
                        estimated_radius_volume_gain_bits=gain * d)
                    save(improved / "embedding.json", embedding_new)
                    adjusted = {key: val for key, val in base.items() if key != "payload_sha256"}
                    adjusted.update(initial_best=str(value), initial_y=y,
                        previous_prepared_sha256=base["payload_sha256"],
                        previous_initial_best=base["initial_best"],
                        initial_candidate_sha256=engine._hash_file(improved / "initial_candidate.json"))
                    assert engine._hash_file(target_native) == adjusted["preparation_sha256"]
                    engine._save_document(new_stage / "prepared.json", adjusted)
                    print("ADOPTED_SIEVE_RADIUS", index, "length_gain_bits", gain,
                          "estimated_volume_gain_bits", gain * d, flush=True)
                folder = improved
                embedding = read(folder / "embedding.json")
    if (folder / "minimum.json").exists() and read(folder / "minimum.json").get("exact"):
        return folder
    from partition_search import search
    os.environ.setdefault("PARTITION_CHUNK_SECONDS", "10")
    stage = folder / "certificate" / "pass1"
    prepared = engine._decode_preparation(engine._read_document(stage / "prepared.json"))
    with engine._lock(folder / "certificate"):
        stage_result = search(prepared, stage, None, workers=workers, max_seconds=seconds)
    result = dict(exact=stage_result["complete"], status=stage_result["reason"], passes=[stage_result])
    save(folder / "certificate" / "result.json", result)
    squared = QQ(stage_result["off_line_sq"]) * scale ** 2
    length_log = (log2_int(squared.numerator()) - log2_int(squared.denominator())) / 2
    minimum = dict(case=embedding["case"], kind="full", exact=result["exact"],
        squared_norm_numerator=str(squared.numerator()), squared_norm_denominator=str(squared.denominator()),
        log2_length=str(length_log) if result["exact"] else None, log2_upper_bound=str(length_log),
        integer_embedding_witness=stage_result["off_line_witness"],
        certificate=str((folder / "certificate").relative_to(HERE)),
        embedding_sha256=engine._hash_file(folder / "embedding.json"), backend="partition",
        reason=result["status"])
    save(folder / "minimum.json", minimum)
    if os.environ.get("PROJECTED_DEFER_PUBLICATION") != "1":
        summarize()
    print("STRONG_PROJECTED_MINIMUM", index, result["exact"], str(length_log), flush=True)
    return folder


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=("prepare", "check", "run", "whole", "watch", "stop-weak"))
    p.add_argument("--indices", default="7,0,1,2,3,8,9")
    p.add_argument("--seconds", type=float, default=120)
    p.add_argument("--workers", type=int, default=1)
    args = p.parse_args()
    if args.mode == "stop-weak":
        expected = [str(HERE / "projected_minima.py"), "run", "--kind", "full", "--indices",
                    "7,0,1,2,3,8,9", "--seconds", "60", "--workers", "1", "--backend", "partition"]
        processes = {}
        for entry in Path("/proc").iterdir():
            if not entry.name.isdigit():
                continue
            try:
                argv = [s.decode() for s in (entry / "cmdline").read_bytes().split(b"\0") if s]
                status = (entry / "status").read_text()
                parent = int(next(line.split()[1] for line in status.splitlines() if line.startswith("PPid:")))
                processes[int(entry.name)] = (parent, argv)
            except (FileNotFoundError, ProcessLookupError, StopIteration):
                continue
        matches = [pid for pid, (_, argv) in processes.items() if argv[1:] == expected]
        assert len(matches) <= 1, matches
        for pid in matches:
            children = [(child, argv) for child, (parent, argv) in processes.items() if parent == pid]
            for child, argv in children:
                assert argv and "exact_enum-" in argv[0] and any(str(OUT) in arg for arg in argv[1:]), argv
                print("STOP_OWN_WEAK_PREPARATION", child, flush=True)
                os.kill(child, signal.SIGINT)
            print("STOP_OWN_WEAK_BATCH", pid, flush=True)
            os.kill(pid, signal.SIGINT)
        return
    if args.mode == "watch":
        indices = list(map(int, args.indices.split(",")))
        started = datetime.now(timezone(timedelta(hours=9))).isoformat()
        failures = {index: 0 for index in indices}
        cycle = 0
        def finished_and_audited(index):
            case_folder = OUT / prepare_geometry(index)["case"]
            folders = [case_folder / "full_strong"] + list(case_folder.glob("full_sieve_*"))
            for folder in folders:
                minimum_path = folder / "minimum.json"
                audit_path = folder / "audit.json"
                if (minimum_path.exists() and read(minimum_path).get("exact")
                        and audit_path.exists() and read(audit_path).get("passed")
                        and read(audit_path).get("minimum_sha256") == engine._hash_file(minimum_path)):
                    return True
            return False
        while not (OUT / "STOP_PROJECTED_SEARCH").exists():
            cycle += 1
            pending = [index for index in indices if not finished_and_audited(index)]
            if not pending:
                save(OUT / "worker_progress.json", dict(pid=os.getpid(), state="complete", indices=indices,
                     started_kst=started, cycle=cycle, workers=args.workers))
                summarize()
                return
            viable = [index for index in pending if failures[index] < 3]
            if not viable:
                save(OUT / "worker_progress.json", dict(pid=os.getpid(), state="failed", indices=indices,
                     started_kst=started, cycle=cycle, failures=failures, workers=args.workers))
                return
            for index in viable:
                if (OUT / "STOP_PROJECTED_SEARCH").exists():
                    break
                progress = dict(pid=os.getpid(), state="running", index=index, indices=indices, cycle=cycle,
                    started_kst=started, updated_kst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
                    workers=args.workers, seconds_per_case=args.seconds, failures=failures)
                save(OUT / "worker_progress.json", progress)
                try:
                    result_folder = run(index, args.seconds, args.workers)
                    minimum_path = result_folder / "minimum.json"
                    if minimum_path.exists() and read(minimum_path).get("exact"):
                        audit_path = minimum_path.parent / "audit.json"
                        if not audit_path.exists() or read(audit_path).get("minimum_sha256") != engine._hash_file(minimum_path):
                            audit_command=[sys.executable,str(HERE / "validate_projected_reuse.py"),
                                           "--indices",str(index),"--variant",result_folder.name]
                            saved_run=read(result_folder/'certificate'/'result.json')['passes'][0]
                            if 'partition_source_sha256' in saved_run:
                                audit_command.append('--partition')
                            subprocess.run(audit_command,check=True)
                    subprocess.run([sys.executable, str(HERE / "write_projected_comparison.py")], check=True,
                                   stdout=subprocess.DEVNULL)
                    failures[index] = 0
                except Exception as exc:
                    failures[index] += 1
                    progress.update(state="case_failed", error=repr(exc), failures=failures)
                    save(OUT / "worker_progress.json", progress)
                    print("PROJECTED_CASE_FAILED", index, repr(exc), flush=True)
        save(OUT / "worker_progress.json", dict(pid=os.getpid(), state="stopped", indices=indices,
             started_kst=started, cycle=cycle, workers=args.workers))
        return
    for index in map(int, args.indices.split(",")):
        if args.mode == "prepare":
            prepare(index)
        elif args.mode == "check":
            _, _, embedding = prepare(index)
            stage = OUT / embedding["case"] / "full_strong" / "certificate" / "pass1"
            from fast_search import verify_lll_precondition
            data = engine._decode_preparation(engine._read_document(stage / "prepared.json"))
            checked = verify_lll_precondition(data, stage)
            assert engine._witness(data, data["raw"]["initial_y"])[0] == int(data["raw"]["initial_best"])
            print("PROJECTED_PRECONDITIONS_OK", index, checked["exact_integer_coefficient_bound"], flush=True)
        elif args.mode == "whole":
            _, scale, embedding = prepare(index)
            folder=OUT/embedding['case']/'full_strong'
            stage=folder/'certificate'/'pass1'
            from fast_search import search as whole_search
            data=engine._decode_preparation(engine._read_document(stage/'prepared.json'))
            with engine._lock(folder/'certificate'):
                proof=whole_search(data,stage,None,max_seconds=args.seconds)
                if not proof['complete']:
                    print('WHOLE_SVP_UNFINISHED',index,flush=True)
                    continue
                squared=QQ(proof['off_line_sq'])*scale**2
                length_log=(log2_int(squared.numerator())-log2_int(squared.denominator()))/2
                save(folder/'certificate'/'result.json',dict(exact=True,status='complete',passes=[proof]))
                minimum=dict(case=embedding['case'],kind='full',exact=True,
                    squared_norm_numerator=str(squared.numerator()),squared_norm_denominator=str(squared.denominator()),
                    log2_length=str(length_log),log2_upper_bound=str(length_log),
                    integer_embedding_witness=proof['off_line_witness'],
                    certificate=str((folder/'certificate').relative_to(HERE)),
                    embedding_sha256=engine._hash_file(folder/'embedding.json'),backend='whole_error_bounded_fplll',reason='complete')
                save(folder/'minimum.json',minimum)
                subprocess.run([sys.executable,str(HERE/'validate_projected_reuse.py'),'--indices',str(index)],check=True)
                subprocess.run([sys.executable,str(HERE/'write_projected_comparison.py')],check=True,stdout=subprocess.DEVNULL)
                print('WHOLE_SVP_CERTIFIED',index,str(length_log),flush=True)
        else:
            run(index, args.seconds, args.workers)


if __name__ == "__main__":
    main()
