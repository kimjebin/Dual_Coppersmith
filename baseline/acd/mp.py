import logging
import os
import sys
from itertools import product
from math import gcd

from sage.all import ZZ

path = os.path.dirname(os.path.dirname(os.path.realpath(os.path.abspath(__file__))))
if path not in sys.path:
    sys.path.insert(0, path)

from shared import small_roots
from shared.reduction import reduce_integer_basis


def attack(N, a, rho, t=1, k=1, roots_method="groebner", use_flatter=False):
    """
    Solves the ACD problem using the multivariate polynomial approach.
    More information: Galbraith D. S. et al., "Algorithms for the Approximate Common Divisor Problem" (Section 5)
    :param N: N = p * q0
    :param a: the a samples, with ai = p * qi + ri
    :param rho: the bit length of the r values
    :param t: the parameter t (default: 1)
    :param k: the parameter k (default: 1)
    :param roots_method: the method to use to find roots (default: "groebner")
    :return: the secret integer p and a list containing the r values, or None if p could not be found
    """
    assert len(a) > 0, "At least one a value is required."
    assert t >= k, "t must be greater than or equal to k."

    R = 2 ** rho

    pr = ZZ[tuple(f"x{i}" for i in range(len(a)))]
    x = pr.gens()
    X = [R] * len(x)

    logging.debug("Generating shifts...")

    shifts = []
    for i in product(*[range(t + 1) for _ in x]):
        if sum(i) <= t:
            l = max(k - sum(i), 0)
            fi = N ** l
            for m in range(len(i)):
                fi *= (x[m] - a[m]) ** i[m]

            shifts.append(fi)

    B, monomials = small_roots.create_lattice(pr, shifts, X)
    # Use exactly the same fpLLL backend and delta as the dual pipeline.
    B = reduce_integer_basis(B, use_flatter=use_flatter)
    polynomials = small_roots.reconstruct_polynomials(B, None, N ** k, monomials, X)
    for roots in small_roots.find_roots(pr, polynomials, method=roots_method):
        r = [roots[xi] for xi in x]
        if all(-R < ri < R for ri in r):
            # Each a[i] - r[i] = p * q_i, so g = gcd(N, a[0]-r[0], ..., a[m-1]-r[m-1])
            # = p * gcd(q0, q_1, ...) = p * c, where c is a small spurious common factor of
            # the random q_i (usually 1).  p is the large prime, hence the LARGEST prime
            # factor of g; taking it recovers p exactly and unconditionally (never p*c).
            g = gcd(N, *(a[i] - r[i] for i in range(len(a))))
            if not (1 < g < N):
                continue
            p_recovered = int(max(ZZ(g).prime_divisors()))
            return p_recovered, r
