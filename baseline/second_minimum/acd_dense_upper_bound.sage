#!/usr/bin/env sage
"""Compute the reduced-basis upper bound for dense ACD second minima."""

import argparse
from math import comb, floor, log2
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("SAGE_NUM_THREADS", "1")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from sage.all import (  # noqa: E402
    QQ,
    ZZ,
    RealField,
    diagonal_matrix,
    gcd as sage_gcd,
    randint,
    random_prime,
    set_random_seed,
    vector,
)

from acd.mp_dual import _build_shifts  # noqa: E402
from shared import small_roots  # noqa: E402
from shared.reduction import reduce_integer_basis  # noqa: E402


ETA = 100
RHO = 10
RF = RealField(256)


def gamma_max(ell: int, t: int, m: int) -> float:
    dim = comb(t + ell, ell)
    loss = (
        dim * log2(dim)
        + dim * dim * log2(1.02)
        + RHO * dim * (ell * t / (ell + 1))
    )
    coefficient = comb(m + ell, ell) * m / (ell + 1)
    return (m * ETA * dim - loss) / coefficient


def generate_instance(t: int, m: int, ell: int, seed: int):
    set_random_seed(seed)
    gamma = int(floor(0.9 * gamma_max(ell, t, m)))
    p = random_prime(ZZ(2) ** ETA - 1, lbound=ZZ(2) ** (ETA - 1))
    q0 = randint(ZZ(2) ** (gamma - ETA - 1), ZZ(2) ** (gamma - ETA))
    modulus = p * q0
    samples = []
    roots = []
    for _ in range(ell):
        qi = randint(ZZ(2) ** (gamma - ETA - 1), ZZ(2) ** (gamma - ETA))
        ri = randint(-(2 ** RHO) + 1, 2 ** RHO - 1)
        samples.append(p * qi + ri)
        roots.append(ri)
    return modulus, samples, roots, gamma


def build_basis(t: int, m: int, ell: int, seed: int):
    modulus, samples, roots, gamma = generate_instance(t, m, ell, seed)
    ring = ZZ[tuple(f"x{i}" for i in range(ell))]
    variables = ring.gens()
    shifts = _build_shifts(ring, variables, modulus, samples, t, m)
    matrix_a, monomials = small_roots.create_lattice(ring, shifts, [1] * ell)
    modulus_power = ZZ(modulus) ** m
    smith, _left, right = matrix_a.smith_form()
    divisors = [ZZ(smith[i, i]) for i in range(min(matrix_a.nrows(), matrix_a.ncols()))]
    multipliers = []
    for index in range(matrix_a.ncols()):
        divisor = divisors[index] if index < len(divisors) else ZZ(0)
        multipliers.append(
            ZZ(1) if divisor == 0 else modulus_power // sage_gcd(modulus_power, divisor)
        )
    scales = [ZZ(2) ** (RHO * (t - monomial.degree())) for monomial in monomials]
    basis = diagonal_matrix(ZZ, multipliers) * right.transpose()
    return basis * diagonal_matrix(ZZ, scales), monomials, roots, scales, gamma


def norm_squared(row):
    return ZZ(sum(ZZ(value) ** 2 for value in row))


def log2_norm(row) -> float:
    return float(RF(norm_squared(row)).log(2) / 2)


def proportional(left, right) -> bool:
    pivot = next((i for i, value in enumerate(right) if value), None)
    if pivot is None:
        return False
    return all(
        QQ(left[i]) * QQ(right[pivot]) == QQ(right[i]) * QQ(left[pivot])
        for i in range(len(right))
    )


def ratio_decode(row, monomials, scales, roots) -> bool:
    unbalanced = [QQ(row[i]) / QQ(scales[i]) for i in range(len(scales))]
    zero = next(i for i, monomial in enumerate(monomials) if monomial.degree() == 0)
    if unbalanced[zero] == 0:
        return False
    for variable, root in enumerate(roots):
        unit = next(
            i
            for i, monomial in enumerate(monomials)
            if tuple(monomial.degrees())
            == tuple(1 if j == variable else 0 for j in range(len(roots)))
        )
        if unbalanced[unit] / unbalanced[zero] != root:
            return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("t", type=int)
    parser.add_argument("m", type=int)
    parser.add_argument("ell", type=int)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--backend", choices=("fplll", "lll"), default="fplll",
                        help="fpLLL:proved; lll is a compatibility alias")
    args = parser.parse_args()

    basis, monomials, roots, scales, gamma = build_basis(
        args.t, args.m, args.ell, args.seed
    )
    reduced = reduce_integer_basis(basis)
    order = sorted(range(reduced.nrows()), key=lambda i: norm_squared(reduced.row(i)))
    root_row = vector(ZZ, reduced.row(order[0]))
    if not ratio_decode(root_row, monomials, scales, roots):
        raise RuntimeError("the shortest reduced row does not decode to the planted root")
    independent = next(
        vector(ZZ, reduced.row(index))
        for index in order[1:]
        if not proportional(reduced.row(index), root_row)
    )
    print(
        "UPPER_BOUND t=%d m=%d ell=%d dim=%d seed=%d backend=%s "
        "gamma=%d log2_lambda1=%.6f log2_lambda2_ub=%.6f"
        % (
            args.t,
            args.m,
            args.ell,
            reduced.nrows(),
            args.seed,
            "fpLLL:proved",
            gamma,
            log2_norm(root_row),
            log2_norm(independent),
        )
    )


if __name__ == "__main__":
    main()
