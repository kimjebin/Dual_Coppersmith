import itertools
import json
import math
from pathlib import Path
import random
import unittest
from sage.all import ZZ, matrix
from coefficient_certificate import restrict, bounds


def exhaustive(basis, upper):
    inverse = (basis*basis.transpose()).inverse()
    limits = [math.isqrt(int(upper*inverse[i,i])) for i in range(basis.nrows())]
    candidates = []
    for c in itertools.product(*(range(-b, b+1) for b in limits)):
        if not any(c):
            continue
        v = matrix(ZZ, [c])*basis
        sq = sum(z*z for z in v.list())
        if sq <= upper:
            candidates.append((sq, v.row(0)))
    candidates.sort(key=lambda x: x[0])
    one, v = candidates[0]
    two = next(sq for sq, w in candidates if matrix([v,w]).rank() == 2)
    return one, two


class CoefficientCertificateTests(unittest.TestCase):
    def check_basis(self, B):
        order = sorted(sum(z*z for z in row) for row in B.rows())
        upper = order[1]
        restricted, proof = restrict(B, upper)
        active, verified = bounds(B, upper, proof)
        self.assertEqual(active, proof['active_rows'])
        self.assertEqual(verified['excluded'], proof['excluded'])
        self.assertEqual(exhaustive(B, upper), exhaustive(restricted, upper))
        return active

    def test_excludes_long_direction_and_preserves_ties(self):
        self.assertEqual(self.check_basis(matrix(ZZ, [[2,0,0],[1,3,0],[0,0,100]])), [0,1])
        self.assertEqual(self.check_basis(matrix(ZZ, [[2,0,0],[0,2,0],[0,0,2]])), [0,1,2])

    def test_complete_small_spheres(self):
        rng = random.Random(20261028)
        for n in (2,3,4):
            for _ in range(4):
                B = matrix(ZZ, [[rng.randint(2,5) if i == j else
                    (rng.randint(-2,2) if j<i else 0) for j in range(n)] for i in range(n)])
                self.check_basis(B)

    def test_very_large_coordinate_scales(self):
        B = matrix(ZZ, [[2**3000,0,0],[2**2999,3*2**3000,0],[0,0,2**3300]])
        self.assertEqual(self.check_basis(B), [0,1])

    def test_one_remaining_direction_certifies_first_minimum(self):
        B = matrix(ZZ, [[2,0,0],[1,100,0],[0,0,200]])
        restricted, proof = restrict(B, ZZ(4))
        self.assertEqual(proof['active_rows'], [0])
        self.assertEqual(restricted.row(0).dot_product(restricted.row(0)), 4)


if __name__ == '__main__':
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(CoefficientCertificateTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    (Path(__file__).parent/'logs'/'coefficient-tests.json').write_text(json.dumps(
        dict(tests=result.testsRun, failures=len(result.failures), errors=len(result.errors)))+'\n')
    raise SystemExit(not result.wasSuccessful())
