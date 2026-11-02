#!/usr/bin/env sage-python
"""Small deterministic regression tests for the public artifact."""

import sys
import unittest

ROOT = __file__.rsplit("/tests/", 1)[0]
sys.path.insert(0, ROOT)

from sage.all import PolynomialRing, ZZ, matrix

from shared.dual_moment import (
    howell_congruence_kernel,
    small_roots as dual_small_roots,
    snf_congruence_kernel,
)
from shared.reduction import reduce_integer_basis


class CoreTests(unittest.TestCase):
    def test_kernel_constructions_define_same_lattice(self):
        a = matrix(ZZ, [[2, 4, 1], [6, 3, 9]])
        snf = snf_congruence_kernel(a, 12).hermite_form()
        howell = howell_congruence_kernel(a, 12).hermite_form()
        self.assertEqual(snf, howell)

    def test_fplll_preserves_determinant(self):
        basis = matrix(ZZ, [[19, 2], [7, 11]])
        reduced, backend = reduce_integer_basis(
            basis, use_flatter=False, return_backend=True
        )
        self.assertEqual(backend, "fpLLL:proved")
        self.assertEqual(abs(reduced.det()), abs(basis.det()))

    def test_dual_small_roots_accepts_single_polynomial(self):
        ring = PolynomialRing(ZZ, names=("x",))
        x = ring.gen()
        self.assertEqual(
            dual_small_roots(x - 1, 101, [4], m=1, t=1,
                              use_flatter=False, kernel_method="snf"),
            [(1,)],
        )

    def test_dual_small_roots_accepts_two_polynomial_system(self):
        ring = PolynomialRing(ZZ, names=("x", "y"))
        x, y = ring.gens()
        self.assertEqual(
            dual_small_roots([x - 1, y + 2], 101, [4, 4], m=1, t=1,
                              use_flatter=False, kernel_method="snf"),
            [(1, -2)],
        )


if __name__ == "__main__":
    unittest.main()
