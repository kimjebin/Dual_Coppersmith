"""Independent geometry, embedding and enumeration checks for projected SVP."""
import argparse
import json
import sqlite3

from projected_minima import OUT, read, save, GRID, build, name, norm2, matrix, vector, ZZ, QQ, gcd
from experiment import diagonal_matrix, variant
from shared.dual_moment import snf_congruence_kernel
from shared import exact_minima as engine
from audit_results import validate_stage, measurement
from coefficient_certificate import bounds
from fast_search import search as independent_search


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kind", choices=("jet", "full", "all"), default="all")
    parser.add_argument("--indices", default="all")
    args = parser.parse_args()
    indices = range(len(GRID)) if args.indices == "all" else map(int, args.indices.split(","))
    results = []
    for index in indices:
        cfg = GRID[index]
        case_folder = OUT / name(cfg)
        geo_path = case_folder / "geometry.json"
        if not geo_path.exists():
            continue
        geo = read(geo_path)
        data = build(cfg)
        _, exponents, scalar, original_root = variant(data, "B")
        assert list(map(str, original_root)) == geo["root"]
        assert str(norm2(original_root)) == geo["root_sq"]
        jet = matrix(ZZ, geo["jet_basis_rows"])
        balanced = data["jet"] * diagonal_matrix(ZZ, data["scales"])
        coefficient_matrix = balanced.solve_left(jet).change_ring(ZZ)
        assert coefficient_matrix * balanced == jet
        congruences = diagonal_matrix(ZZ, [data["N"] ** (cfg[2] - e) for e in exponents]) * data["forward"] * data["jet"].transpose()
        canonical_kernel = snf_congruence_kernel(congruences, data["N"] ** cfg[2])
        assert abs(coefficient_matrix.det()) == abs(canonical_kernel.det())
        assert all(z % data["N"] ** cfg[2] == 0 for z in (coefficient_matrix * congruences.transpose()).list())
        assert jet.row(0) == original_root
        assert str((jet * jet.transpose()).det()) == geo["jet_gram_determinant"]
        for kind in (("jet", "full") if args.kind == "all" else (args.kind,)):
            folder = case_folder / kind
            minimum_path = folder / "minimum.json"
            if not minimum_path.exists() or not read(minimum_path).get("exact"):
                continue
            minimum = read(minimum_path)
            embedding = read(folder / "embedding.json")
            assert minimum["embedding_sha256"] == engine._hash_file(folder / "embedding.json")
            embedded = matrix(ZZ, embedding["rows"])
            scale = QQ(embedding["scale_numerator"]) / ZZ(embedding["scale_denominator"])
            if kind == "jet":
                source_basis = jet
            else:
                from experiment import HERE
                source_path = HERE / embedding["source"]
                assert engine._hash_file(source_path) == embedding["source_sha256"]
                _, full_basis = measurement(data, source_path, "B")
                common = gcd(full_basis.list())
                normalized = full_basis.apply_map(lambda z: z // common)
                proof = read(folder / "coefficient_exclusion.json")
                active, checked = bounds(normalized, ZZ(proof["upper_sq"]), proof)
                assert active == proof["active_rows"] and checked["excluded"] == proof["excluded"]
                root = normalized.row(0)
                a = norm2(root)
                candidate = min(a * norm2(row) - row.dot_product(root) ** 2 for row in normalized.rows()[1:])
                assert QQ(proof["upper_sq"]) >= QQ(candidate) / a + QQ(a) / 4
                source_basis = full_basis.matrix_from_rows(active)
            root = source_basis.row(0)
            a = norm2(root)
            tail = source_basis[1:].change_ring(QQ)
            projected = tail - matrix(QQ, [QQ(row.dot_product(root)) / a * root for row in tail.rows()])
            assert embedded[0, 0] == 1 and all(z == 0 for z in embedded.row(0)[1:])
            assert all(z == 0 for z in embedded.column(0)[1:])
            assert embedded[1:, 1:] * scale == projected
            assert embedded.nrows() == projected.nrows() + 1
            checkpoint = folder / "certificate"
            result = read(checkpoint / "result.json")
            assert result["exact"] and result["lambda1_sq"] == "1" and len(result["passes"]) == 1
            prepared = validate_stage(checkpoint / "pass1", embedded)
            run = result["passes"][0]
            assert run["complete"] and run["reason"] == "complete"
            best, witness = engine._witness(prepared, list(map(int, run["quotient_coefficients"])))
            assert str(best) == result["lambda2_sq"] == run["off_line_sq"]
            assert list(map(str, witness)) == minimum["integer_embedding_witness"]
            assert witness[0] == 0
            integer_witness = vector(ZZ, witness)
            coordinates = embedded.solve_left(matrix(ZZ, [integer_witness])).change_ring(ZZ)
            assert coordinates * embedded == matrix(ZZ, [integer_witness])
            length_sq = QQ(best) * scale ** 2
            assert str(length_sq.numerator()) == minimum["squared_norm_numerator"]
            assert str(length_sq.denominator()) == minimum["squared_norm_denominator"]
            explicit = vector(QQ, minimum["projected_witness_numerator"]) / ZZ(minimum["projected_witness_denominator"])
            assert explicit == scale * integer_witness[1:] and explicit.dot_product(root) == 0
            assert explicit.dot_product(explicit) == length_sq
            # An independent, error-bounded fplll sphere enumeration is feasible
            # for every small-dimensional result and checks the GMP DFS answer.
            replay = independent_search(prepared, checkpoint / "pass1", None, max_seconds=60)
            assert replay["complete"] and replay["off_line_sq"] == str(best)
            audit = dict(case=name(cfg), kind=kind, passed=True, exact_squared_norm_numerator=str(length_sq.numerator()),
                         exact_squared_norm_denominator=str(length_sq.denominator()),
                         original_geometry_verified=True, projection_embedding_verified=True,
                         global_coefficient_exclusion_verified=(kind == "full"),
                         native_gmp_and_independent_fplll_agree=True,
                         minimum_sha256=engine._hash_file(minimum_path))
            save(folder / "audit.json", audit)
            results.append(audit)
            print("PROJECTED_AUDIT_OK", index, kind, flush=True)
    save(OUT / ("audit-" + args.kind + ".json"), dict(passed=len(results), results=results))


if __name__ == "__main__":
    main()
