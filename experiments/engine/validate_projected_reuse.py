"""Audit reused-Gram projected SVP with a fresh independent enumeration."""
import argparse
import json
import sqlite3
import subprocess
import re
from projected_minima import OUT, HERE, GRID, read, save, prepare_geometry, matrix, vector, ZZ, QQ, norm2
from shared import exact_minima as engine
from gram_validation import validate as validate_gram
from fast_search import search, verify_lll_precondition


def verify_partition(prepared, stage, run):
    from partition_search import binaries, children
    from suffix_filter import load_filter, proven_empty, candidate_interval
    binary, source_hash, bounder, bounder_hash = binaries()
    assert run["complete"] and run["partition_source_sha256"] == source_hash
    assert run["partition_binary_sha256"] == engine._hash_file(binary)
    assert run["prepared_sha256"] == prepared["raw"]["payload_sha256"]
    verify_lll_precondition(prepared, stage)
    limits = read(stage / "partition_bounds.json")
    assert limits["bounder_source_sha256"] == bounder_hash
    assert run["bound_file_sha256"] == engine._hash_file(stage / "partition_bounds.json")
    checked = subprocess.run([str(bounder), str(stage / "integer.prepared"), str(prepared["raw"]["initial_best"])],
                             capture_output=True, text=True, check=True)
    bounds = list(map(int, checked.stdout.splitlines()[1].split()[1:]))
    assert bounds == list(map(int, limits["bounds"]))
    filter_doc = load_filter(prepared, stage, recheck=True)
    filter_hash = engine._hash_file(stage / "suffix_filter.json")
    assert run["suffix_filter_sha256"] == filter_hash
    assert run["input_sha256"] == engine._hash_file(stage / "partition.input")
    T = matrix(ZZ, prepared["tail"])
    a, c = ZZ(prepared["root_sq"]), ZZ(prepared["content"])
    inner = vector(ZZ, prepared["inner"])
    H = (a * T * T.transpose() - matrix(ZZ, T.nrows(), T.nrows(), lambda i, j: inner[i] * inner[j])).apply_map(lambda z: z // c)
    with (stage / "partition.input").open() as f:
        assert f.readline().split() == ["FAST-LIFT-1", str(T.nrows()), str(a), str(c), str(prepared["raw"]["initial_best"])]
        assert f.readline().split() == list(map(str, inner))
        assert f.readline().split() == list(map(str, prepared["raw"]["initial_y"]))
        for row in H.rows():
            assert f.readline().split() == list(map(str, row))
    db = sqlite3.connect((stage / "partition.sqlite3").resolve().as_uri() + "?mode=ro", uri=True)
    try:
        metadata = dict(db.execute("SELECT key,value FROM metadata"))
        assert all(str(run[key]) == value for key, value in metadata.items())
        jobs = {key: (parent, status, json.loads(value) if value else None)
                for key, parent, status, value in db.execute("SELECT suffix,parent,status,result FROM jobs")}
    finally:
        db.close()
    assert jobs["[]"][0] is None
    children_by_parent = {}
    for key, (parent, status, proof) in jobs.items():
        assert status in ("done", "split")
        if parent is not None:
            assert parent in jobs
            children_by_parent.setdefault(parent, set()).add(key)
    best_values = []
    accelerator_identity = None
    transfer = None
    for key, (parent, status, proof) in jobs.items():
        suffix = json.loads(key)
        if status == "split":
            if proof and proof.get("method") == "dyadic_child_interval":
                lo, hi = candidate_interval(suffix, bounds, filter_doc)
                assert proof["filter_sha256"] == filter_hash
                assert proof["lower"] == str(lo) and proof["upper"] == str(hi) and lo <= hi
                expected = {json.dumps([j] + suffix, separators=(",", ":")) for j in range(lo, hi + 1)}
            else:
                expected = {json.dumps(child, separators=(",", ":")) for child in children(suffix, bounds)}
            assert children_by_parent.get(key, set()) == expected and expected
            continue
        assert proof and not children_by_parent.get(key)
        if proof.get("method") == "original_nearest_lift_subtree_lower_bound":
            if transfer is None:
                from projection_subtree_transfer import Transfer
                transfer = Transfer(prepared, stage)
            transfer.verify(suffix, proof)
        elif proof.get("method") in ("exact_suffix_lower_bound", "exact_empty_child_interval"):
            assert proof["empty_filter_sha256"] == filter_hash
            if proof["method"] == "exact_suffix_lower_bound":
                assert proven_empty(suffix, filter_doc)
            else:
                lo, hi = candidate_interval(suffix, bounds, filter_doc)
                assert lo > hi
        elif proof.get("subtree_enumerator"):
            if accelerator_identity is None:
                from unrolled_search import build
                accelerator, identity = build(len(prepared["tail"]))
                accelerator_identity = (identity, engine._hash_file(accelerator))
            assert proof["subtree_enumerator"] == "unrolled_enumlib"
            assert (proof["unrolled_source_identity"], proof["unrolled_binary_sha256"]) == accelerator_identity
        else:
            assert not proof.get("method"), "Unsupported completed-subtree proof method"
        value, _ = engine._witness(prepared, list(map(int, proof["y"])))
        assert value == int(proof["best"]) <= int(prepared["raw"]["initial_best"])
        if value < int(prepared["raw"]["initial_best"]) and suffix:
            assert list(map(int, proof["y"][-len(suffix):])) == suffix
        best_values.append(value)
    assert best_values and min(best_values) == int(run["off_line_sq"])
    return run


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--indices", default="7")
    parser.add_argument("--partition", action="store_true")
    parser.add_argument("--variant", default="full_strong")
    args = parser.parse_args()
    assert args.variant == "full_strong" or re.fullmatch(r"full_sieve_[0-9a-f]{12}", args.variant)
    for index in map(int, args.indices.split(",")):
        geo = prepare_geometry(index)
        folder = OUT / geo["case"] / args.variant
        minimum_path = folder / "minimum.json"
        minimum = read(minimum_path)
        assert minimum["exact"]
        embedding = read(folder / "embedding.json")
        assert engine._hash_file(folder / "embedding.json") == minimum["embedding_sha256"]
        source = HERE / embedding["source_preparation"]
        assert engine._hash_file(source / "prepared.json") == embedding["source_preparation_sha256"]
        old = engine._read_document(source / "prepared.json")
        measured_path = HERE / embedding["original_measurement"]
        assert engine._hash_file(measured_path) == embedding["original_measurement_sha256"]
        from audit_results import measurement
        from experiment import build, gcd
        _, original_basis = measurement(build(GRID[index]), measured_path, "B")
        B = matrix(ZZ, embedding["rows"])
        scale = QQ(embedding["scale_numerator"]) / ZZ(embedding["scale_denominator"])
        root = vector(ZZ, old["root"])
        T = matrix(ZZ, old["tail"])
        a = norm2(root)
        ratio = ZZ(geo["root_sq"]) // a
        assert ratio.is_square()
        physical_scale = ratio.isqrt()
        assert physical_scale == gcd(original_basis.list())
        assert physical_scale * root in (vector(ZZ, geo["root"]), -vector(ZZ, geo["root"]))
        original_scaled = original_basis.apply_map(lambda z: z // physical_scale)
        full_source_basis = matrix(ZZ, [root] + list(T.rows()))
        transition = original_scaled.solve_left(full_source_basis).change_ring(ZZ)
        assert transition * original_scaled == full_source_basis and abs(transition.det()) == 1
        expected_projection = T.change_ring(QQ) - matrix(QQ, [QQ(row.dot_product(root)) / a * root for row in T.rows()])
        assert B[1:, 1:] * scale == physical_scale * expected_projection
        stage = folder / "certificate" / "pass1"
        prepared = engine._decode_preparation(engine._read_document(stage / "prepared.json"))
        assert prepared["root_sq"] == 1 and not any(prepared["inner"])
        assert prepared["tail"] == [[int(z) for z in row] for row in B[1:].rows()]
        validate_gram(prepared, stage)
        if args.partition:
            result = read(folder / "certificate" / "result.json")
            replay = verify_partition(prepared, stage, result["passes"][0])
        else:
            replay = search(prepared, stage, None, max_seconds=180)
        assert replay["complete"]
        squared = QQ(replay["off_line_sq"]) * scale ** 2
        assert str(squared.numerator()) == minimum["squared_norm_numerator"]
        assert str(squared.denominator()) == minimum["squared_norm_denominator"]
        witness = vector(ZZ, replay["off_line_witness"])
        coefficients = B.solve_left(matrix(ZZ, [witness])).change_ring(ZZ)
        assert coefficients * B == matrix(ZZ, [witness]) and witness[0] == 0
        actual = scale * witness[1:]
        assert actual.dot_product(actual) == squared and actual.dot_product(vector(ZZ, geo["root"])) == 0
        save(folder / "audit.json", dict(passed=True, case=geo["case"], kind="full",
            source_measurement_and_full_unimodular_equivalence_independently_verified=True,
            projection_identity_verified=True, reused_bareiss_verified_against_actual_integer_gram=True,
            complete_independent_fplll_replay=not args.partition and not replay.get('certificate_reused',False),
            complete_error_bounded_fplll_search_verified=not args.partition,
            complete_partition_proof_tree_verified=args.partition,
            minimum_sha256=engine._hash_file(minimum_path)))
        print("REUSED_PROJECTED_AUDIT_OK", index, flush=True)


if __name__ == "__main__":
    main()
