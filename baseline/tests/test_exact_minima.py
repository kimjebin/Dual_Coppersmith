#!/usr/bin/env sage-python
"""Exact-minimum certificates, resumable frontiers, and public-basis regressions."""

from itertools import product
import json
from math import isqrt
from pathlib import Path
import random
import runpy
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from sage.all import QQ, ZZ, PolynomialRing, diagonal_matrix, matrix, prod

from acd.mp_dual import _build_shifts
from second_minimum.instances import _affine_kernel, acd_basis, hm_basis
from shared.dual_moment import _exp_vec, _monomials_up_to, snf_congruence_kernel
from shared.exact_minima import (_build_native, _proved_lll, _read_document,
                                _split_frontier, exact_successive_minima,
                                load_minima_input, prepare_successive_minima, verify_minima_bundle)
from shared.reduction import reduce_integer_basis


def exhaustive_minima(basis):
    """Exhaust a provably enclosing coefficient box for rank <= 4 tests."""
    gram = basis * basis.transpose()
    radius = int(max(gram[i, i] for i in range(gram.nrows())))
    inverse = gram.inverse()
    bounds = [isqrt(int(radius * inverse[i, i])) for i in range(gram.nrows())]
    candidates = []
    for coefficients in product(*(range(-b, b + 1) for b in bounds)):
        if not any(coefficients):
            continue
        row = matrix(ZZ, [coefficients]) * basis
        norm = int(sum(value ** 2 for value in row.list()))
        if norm <= radius:
            candidates.append((norm, row.row(0)))
    candidates.sort(key=lambda item: item[0])
    first, vector = candidates[0]
    second = next(norm for norm, row in candidates
                  if matrix(QQ, [vector, row]).rank() == 2)
    return first, second


def original_kernel(ring, shifts, modulus, degree, bound):
    keys = [_exp_vec(monomial) for monomial in _monomials_up_to(ring, degree)]
    rows = []
    for shift in shifts:
        coefficients = {
            tuple(key) if hasattr(key, "__iter__") else (key,): value
            for key, value in shift.dict().items()
        }
        rows.append([coefficients.get(key, 0) for key in keys])
    kernel = snf_congruence_kernel(matrix(ZZ, rows), modulus)
    return kernel * diagonal_matrix(ZZ, [bound ** (degree - sum(key)) for key in keys])


class ExactMinimaTests(unittest.TestCase):
    def run_minima(self, basis, **options):
        with tempfile.TemporaryDirectory(prefix="exact-minima-test-") as checkpoint:
            result = exact_successive_minima(basis, checkpoint, workers=2,
                                             chunk_nodes=5, **options)
        self.assertTrue(result["exact"])
        return int(result["lambda1_sq"]), int(result["lambda2_sq"])

    def test_nearest_integer_lift_changes_projected_order(self):
        self.assertEqual(self.run_minima([[10, 0, 0], [5, 11, 0], [0, 0, 12]]),
                         (100, 144))

    def test_ties_and_rectangular_basis(self):
        self.assertEqual(self.run_minima([[2, 0], [0, 2]]), (4, 4))
        self.assertEqual(self.run_minima([[1, 1, 0], [0, 1, 1]]), (2, 2))

    def test_shortest_reduced_row_is_not_assumed_shortest(self):
        basis = [[128, 0, 0, 0], [64, 96, 0, 0],
                 [64, 48, 72, 0], [64, 48, 36, 54]]
        self.assertEqual(self.run_minima(basis), (4212, 6516))

    def test_random_lattices_against_complete_spheres(self):
        rng = random.Random(43)
        for rank in (2, 3, 4):
            for _ in range(2):
                basis = matrix(ZZ, [[rng.randint(-2, 2) if j < i else
                                     (rng.randint(2, 5) if i == j else 0)
                                     for j in range(rank)] for i in range(rank)])
                with self.subTest(basis=basis):
                    self.assertEqual(self.run_minima(basis), exhaustive_minima(basis))

    def test_large_integer_entries(self):
        scale = 1 << 4096
        self.assertEqual(self.run_minima([[scale, 0], [0, 2 * scale]]),
                         (scale * scale, 4 * scale * scale))

    def test_limit_resume_and_checkpoint_identity(self):
        basis = [[128, 0, 0, 0], [64, 96, 0, 0],
                 [64, 48, 72, 0], [64, 48, 36, 54]]
        with tempfile.TemporaryDirectory(prefix="exact-resume-test-") as checkpoint:
            result = exact_successive_minima(basis, checkpoint, workers=1,
                                             max_seconds=1e-9, chunk_nodes=1)
            self.assertFalse(result["exact"])
            self.assertIsNone(result["lambda2_sq"])
            result = exact_successive_minima(basis, checkpoint, workers=3, chunk_nodes=1)
            self.assertEqual((int(result["lambda1_sq"]), int(result["lambda2_sq"])),
                             (4212, 6516))
            self.assertEqual(json.loads((Path(checkpoint) / "result.json").read_text())
                             ["exact"], True)
            with self.assertRaisesRegex(ValueError, "another basis"):
                exact_successive_minima([[2, 0], [0, 2]], checkpoint)
            metadata = Path(checkpoint) / "pass1" / "prepared.json"
            document = json.loads(metadata.read_text())
            document["root_sq"] = "1"
            metadata.write_text(json.dumps(document))
            with self.assertRaisesRegex(ValueError, "metadata hash"):
                exact_successive_minima(basis, checkpoint)

    def test_rejects_inexact_inputs_and_nonfinite_limits(self):
        with tempfile.TemporaryDirectory(prefix="exact-input-test-") as checkpoint:
            for value in (float("nan"), float("inf")):
                for option in ("chunk_seconds", "progress_seconds", "memory_gb", "max_seconds"):
                    with self.assertRaises(ValueError):
                        exact_successive_minima([[2, 0], [0, 2]], checkpoint,
                                                **{option: value})
            with self.assertRaises(ValueError):
                exact_successive_minima([[2.0, 0], [0, 2]], checkpoint)

    def test_horizontal_split_is_disjoint_complete(self):
        frame = "0 0 0 -10 20 -3 7"
        pieces = _split_frontier([frame], 12)
        visited = []
        for piece in pieces:
            _, _, _, low, high, left, right = map(int, piece.split())
            visited.extend(range(low, min(left, high) + 1))
            visited.extend(range(max(right, low), high + 1))
        self.assertEqual(sorted(visited), list(range(-10, -2)) + list(range(7, 21)))
        self.assertEqual(len(visited), len(set(visited)))


class PortablePreparationTests(unittest.TestCase):
    def test_copied_bundle_skips_rng_and_all_lll_with_different_binary(self):
        basis = matrix(ZZ, [[2, 0, 0], [0, 3, 0], [0, 0, 4]])
        with tempfile.TemporaryDirectory(prefix="portable-test-") as temporary:
            source, copied = Path(temporary) / "source", Path(temporary) / "copied"
            prepared = prepare_successive_minima(basis, source, instance={"seed": 7}, workers=1)
            self.assertEqual(prepared["enumeration_dimension"], 2)
            self.assertFalse(prepared["exact"])
            shutil.copytree(source, copied, ignore=shutil.ignore_patterns(".native", "exact_enum-*"))
            # Rebuild independently, then deliberately change only the compiled
            # ELF's nonloaded trailing bytes. Same source, different binary hash.
            binary = _build_native(copied)
            with binary.open("ab") as handle:
                handle.write(b"portable provenance test\n")
            with patch("shared.exact_minima._proved_lll", side_effect=AssertionError("LLL repeated")), \
                    patch("second_minimum.instances.set_random_seed", side_effect=AssertionError("RNG repeated")):
                saved, instance = load_minima_input(copied)
                self.assertEqual(instance, {"seed": 7})
                result = exact_successive_minima(saved, copied, workers=1)
            self.assertTrue(result["exact"])
            self.assertEqual((result["lambda1_sq"], result["lambda2_sq"]), ("4", "9"))
            self.assertNotEqual(result["native_binary_sha256"],
                                _read_document(copied / "pass1" / "prepared.json")["binary_sha256"])
            with patch("shared.exact_minima._proved_lll", side_effect=AssertionError("LLL repeated")):
                self.assertTrue(verify_minima_bundle(copied)["preparation_complete"])

    def test_completed_ambient_stage_survives_quotient_interruption(self):
        basis = matrix(ZZ, [[2, 0, 0], [0, 3, 0], [0, 0, 4]])
        calls = []

        def interrupted(matrix_self, *args, **options):
            calls.append(matrix_self.nrows())
            if len(calls) == 2:
                raise KeyboardInterrupt()
            return _proved_lll(matrix_self, *args, **options)

        with tempfile.TemporaryDirectory(prefix="portable-stage-test-") as checkpoint:
            with patch("shared.exact_minima._proved_lll", interrupted):
                with self.assertRaises(KeyboardInterrupt):
                    prepare_successive_minima(basis, checkpoint, workers=1)
            self.assertTrue((Path(checkpoint) / "pass1" / "ambient_lll.json").exists())
            self.assertFalse((Path(checkpoint) / "pass1" / "quotient_lll.json").exists())
            calls.clear()

            def resumed(matrix_self, *args, **options):
                calls.append(matrix_self.nrows())
                self.assertEqual(matrix_self.nrows(), 2)
                return _proved_lll(matrix_self, *args, **options)

            with patch("shared.exact_minima._proved_lll", resumed):
                result = prepare_successive_minima(basis, checkpoint, workers=1)
            self.assertTrue(result["preparation_complete"])
            self.assertEqual(calls, [2])

    def test_stage_tamper_and_changed_public_inputs_are_rejected(self):
        basis = matrix(ZZ, [[2, 0], [0, 3]])
        with tempfile.TemporaryDirectory(prefix="portable-tamper-test-") as checkpoint:
            prepare_successive_minima(basis, checkpoint, instance={"seed": 0}, workers=1)
            with self.assertRaisesRegex(ValueError, "public input"):
                prepare_successive_minima(basis, checkpoint, instance={"seed": 1}, workers=1)
            stage = Path(checkpoint) / "pass1" / "ambient_lll.json"
            stage.write_text(stage.read_text().replace('"delta":0.8', '"delta":0.7'))
            with self.assertRaisesRegex(ValueError, "stage hash"):
                prepare_successive_minima(basis, checkpoint, instance={"seed": 0}, workers=1)

    def test_verified_warm_import_and_atomic_frontier_resume(self):
        import shared.exact_minima as engine
        basis = matrix(ZZ, [[2, 0, 0], [0, 3, 0], [0, 0, 4]])
        with tempfile.TemporaryDirectory(prefix="portable-import-test-") as temporary:
            source, target = Path(temporary) / "source", Path(temporary) / "target"
            prepare_successive_minima(basis, source, instance={"monomial_keys": [(0, 1)]}, workers=1)
            # Tuple/list normalization must not reject identical regenerated metadata.
            prepare_successive_minima(basis, source, instance={"monomial_keys": [(0, 1)]}, workers=1)
            exact_successive_minima(basis, source, workers=1)
            original_save = engine._save_document

            def interrupted_save(path, document):
                if path == target / "pass1" / "prepared.json":
                    raise KeyboardInterrupt()
                return original_save(path, document)

            with patch("shared.exact_minima._save_document", interrupted_save):
                with self.assertRaises(KeyboardInterrupt):
                    prepare_successive_minima(basis, target, workers=1, reuse=source)
            self.assertTrue((target / "pass1" / "frontier.sqlite3").exists())
            self.assertFalse((target / "pass1" / "prepared.json").exists())
            result = prepare_successive_minima(basis, target, workers=1, reuse=source)
            self.assertTrue(result["legacy_native_reused"])
            with patch("shared.exact_minima._proved_lll", side_effect=AssertionError("LLL repeated")):
                self.assertTrue(verify_minima_bundle(target)["preparation_complete"])
                exact = exact_successive_minima(basis, target, workers=1)
            self.assertTrue(exact["exact"])
            self.assertEqual(exact["passes"][0]["invocation_nodes"], 0)


class PublicBasisTests(unittest.TestCase):
    def test_acd_affine_kernel_matches_snf(self):
        for ell, degree, lifting in ((1, 2, 1), (2, 3, 2), (3, 2, 2)):
            ring = PolynomialRing(ZZ, [f"x{i}" for i in range(ell)])
            xs = ring.gens()
            samples = [29 + 7 * i for i in range(ell)]
            modulus, bound = ZZ(35), ZZ(4)
            fast, _ = _affine_kernel(ring, [x + a for x, a in zip(xs, samples)],
                                     lambda key: modulus ** min(sum(key), lifting), bound, degree)
            shifts = _build_shifts(ring, xs, modulus, samples, degree, lifting)
            old = original_kernel(ring, shifts, modulus ** lifting, degree, bound)
            self.assertEqual(fast.hermite_form(), old.hermite_form())

    def test_hm_affine_kernel_matches_snf(self):
        for ell in (1, 2, 3):
            ring = PolynomialRing(ZZ, [f"x{i}" for i in range(ell)])
            xs = ring.gens()
            degree, modulus, bound = 3, ZZ(35), ZZ(4)
            normalized = xs[0] + 9 + sum((i + 2) * xs[i] for i in range(1, ell))
            inverse = [xs[0] - 9 - sum((i + 2) * xs[i] for i in range(1, ell))] + list(xs[1:])
            fast, _ = _affine_kernel(ring, inverse, lambda key: modulus ** key[0], bound, degree)
            shifts = [modulus ** (degree - k) * normalized ** k
                      * prod(xs[i + 1] ** exponents[i] for i in range(ell - 1))
                      for k in range(degree + 1)
                      for exponents in product(range(degree - k + 1), repeat=ell - 1)
                      if sum(exponents) <= degree - k]
            old = original_kernel(ring, shifts, modulus ** degree, degree, bound)
            self.assertEqual(fast.hermite_form(), old.hermite_form())

    def test_generators_are_deterministic(self):
        for builder, arguments in ((acd_basis, (2, 1, 2)), (hm_basis, (2,))):
            first, metadata = builder(*arguments)
            second, repeated = builder(*arguments)
            self.assertEqual(first, second)
            self.assertEqual(metadata, repeated)

    def test_paper_instance_generation_matches_existing_scripts(self):
        acd_original = runpy.run_path(str(ROOT / "second_minimum" / "acd_dense_upper_bound.sage"))
        fast, metadata = acd_basis(2, 1, 2)
        old, monomials, roots, _, gamma = acd_original["build_basis"](2, 1, 2, 0)
        columns = {tuple(monomial.degrees()): i for i, monomial in enumerate(monomials)}
        old = old.matrix_from_columns([columns[tuple(key)] for key in metadata["monomial_keys"]])
        self.assertEqual(fast.hermite_form(), old.hermite_form())
        self.assertEqual(metadata["planted_root"], list(map(int, roots)))
        self.assertEqual(metadata["gamma"], gamma)

        hm_original = runpy.run_path(str(ROOT / "second_minimum" / "hm_geometry.sage"))
        fast, metadata = hm_basis(2)
        _, old, keys, _, roots, _ = hm_original["build_lattice"](2)
        self.assertEqual(metadata["monomial_keys"], keys)
        self.assertEqual(metadata["planted_root"], list(map(int, roots)))
        self.assertEqual(fast.hermite_form(), old.hermite_form())

    def test_legacy_flatter_flag_never_launches_flatter(self):
        with patch("subprocess.run", side_effect=AssertionError("External reducer called")):
            reduced, backend = reduce_integer_basis(matrix(ZZ, [[19, 2], [7, 11]]),
                                                     use_flatter=True, return_backend=True)
        self.assertEqual(backend, "fpLLL:proved")
        self.assertEqual(abs(reduced.det()), 195)


if __name__ == "__main__":
    unittest.main()
