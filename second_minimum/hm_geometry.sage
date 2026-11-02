#!/usr/bin/env sage
"""Recompute the HM lambda1, lambda2, ambient-GH, and first-jet columns."""

import argparse
from itertools import product
from math import lgamma, log
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
    RR,
    ZZ,
    PolynomialRing,
    diagonal_matrix,
    gcd,
    inverse_mod,
    lcm,
    matrix,
    prod,
    randint,
    random_prime,
    set_random_seed,
    vector,
)

from shared import dual_moment as dm  # noqa: E402


PBITS = 380
QBITS = 20
X_BITS = 6
DEGREE = 3


def generate_instance(ell: int, seed: int = 0):
    set_random_seed(seed)
    p = random_prime(ZZ(2) ** PBITS - 1, lbound=ZZ(2) ** (PBITS - 1))
    q = random_prime(ZZ(2) ** QBITS - 1, lbound=ZZ(2) ** (QBITS - 1))
    modulus = p * q
    bound = ZZ(2) ** X_BITS
    ring = PolynomialRing(ZZ, [f"x{i}" for i in range(ell)])
    variables = ring.gens()
    root = [randint(-int(bound) + 1, int(bound) - 1) for _ in range(ell)]
    coefficients = [randint(1, int(modulus) - 1) for _ in range(ell)]
    while gcd(coefficients[0], modulus) != 1:
        coefficients[0] = randint(1, int(modulus) - 1)
    constant = -sum(coefficients[i] * root[i] for i in range(ell))
    constant += int(p) * randint(1, int(q))
    polynomial = constant + sum(
        coefficients[i] * variables[i] for i in range(ell)
    )
    return p, q, modulus, polynomial, root, bound, ring, variables, coefficients


def hm_shifts(polynomial, modulus, ring, variables, coefficients):
    normalized = ring(
        {
            exponent: (ZZ(value) * inverse_mod(ZZ(coefficients[0]), modulus))
            % modulus
            for exponent, value in zip(
                polynomial.exponents(), polynomial.coefficients()
            )
        }
    )
    powers = [ring.one()]
    for _ in range(DEGREE):
        powers.append(powers[-1] * normalized)
    shifts = []
    for k in range(DEGREE + 1):
        for exponents in product(
            *[range(DEGREE - k + 1) for _ in range(len(variables) - 1)]
        ):
            if sum(exponents) <= DEGREE - k:
                shifts.append(
                    ZZ(modulus) ** (DEGREE - k)
                    * powers[k]
                    * prod(
                        variables[1 + j] ** exponents[j]
                        for j in range(len(exponents))
                    )
                )
    return shifts


def build_lattice(ell: int, seed: int = 0):
    p, q, modulus, polynomial, root, bound, ring, variables, coefficients = (
        generate_instance(ell, seed)
    )
    monomials = dm._monomials_up_to(ring, DEGREE)
    keys = [dm._exp_vec(monomial) for monomial in monomials]
    index = {key: i for i, key in enumerate(keys)}
    rows = []
    for shift in hm_shifts(polynomial, modulus, ring, variables, coefficients):
        row = [ZZ(0)] * len(monomials)
        for exponent, value in shift.dict().items():
            key = tuple(exponent)
            if key in index:
                row[index[key]] = ZZ(value)
        rows.append(row)
    coefficient_matrix = matrix(ZZ, rows)
    kernel = dm.snf_congruence_kernel(coefficient_matrix, ZZ(modulus) ** DEGREE)

    def bound_power(key):
        return prod(bound ** key[i] for i in range(ell))

    rho = lcm([bound_power(key) for key in keys])
    scales = [rho // bound_power(key) for key in keys]
    balanced = kernel * diagonal_matrix(ZZ, scales)
    log2det = sum(RR(abs(kernel.det())).log(2) for _ in [0])
    log2det += sum(RR(scale).log(2) for scale in scales)
    reduced = dm._reduce(balanced, use_flatter=False)
    return reduced, balanced, keys, scales, root, float(log2det)


def log2_norm(row) -> float:
    norm_squared = sum(ZZ(value) ** 2 for value in row)
    return float(RR(norm_squared).log(2) / 2)


def ambient_gh(log2det: float, log2lambda1: float, dimension: int) -> float:
    projected_dimension = dimension - 1
    log2_unit_ball = (
        (projected_dimension / 2.0) * RR.pi().log(2)
        - RR(lgamma(projected_dimension / 2.0 + 1.0) / log(2.0))
    )
    return float(
        (RR(log2det) - RR(log2lambda1) - log2_unit_ball)
        / projected_dimension
    )


def jet_gh(balanced, keys, scales, root, log2lambda1: float) -> float:
    ell = len(root)
    jet_rows = [[prod(ZZ(root[i]) ** key[i] for i in range(ell)) for key in keys]]
    for derivative in range(ell):
        row = []
        for key in keys:
            if key[derivative] == 0:
                row.append(ZZ(0))
                continue
            exponent = list(key)
            exponent[derivative] -= 1
            row.append(
                ZZ(key[derivative])
                * prod(ZZ(root[i]) ** exponent[i] for i in range(ell))
            )
        jet_rows.append(row)

    jet_space = matrix(QQ, jet_rows)
    perpendicular = jet_space.right_kernel().basis_matrix()
    constraints = []
    for i in range(perpendicular.nrows()):
        row = [QQ(perpendicular[i, j]) / QQ(scales[j]) for j in range(len(keys))]
        denominator = lcm([value.denominator() for value in row])
        constraints.append([ZZ(value * denominator) for value in row])
    coefficient_constraints = matrix(ZZ, constraints) * balanced.transpose()
    kernel = coefficient_constraints.right_kernel().basis_matrix()
    jet_lattice = kernel * balanced
    gram = jet_lattice * jet_lattice.transpose()
    log2_covolume = float(RR(abs(gram.det())).log(2) / 2)
    return 0.5 * float(RR(ell / (2 * RR.pi() * RR(1).exp())).log(2)) + (
        log2_covolume - log2lambda1
    ) / ell


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("ell", type=int, choices=range(2, 10))
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    reduced, balanced, keys, scales, root, log2det = build_lattice(
        args.ell, args.seed
    )
    order = sorted(range(reduced.nrows()), key=lambda i: log2_norm(reduced.row(i)))
    lambda1 = log2_norm(reduced.row(order[0]))
    print(
        "HM_GEOMETRY ell=%d dim=%d seed=%d log2_lambda1=%.6f "
        "log2_lambda_gh=%.6f log2_lambda_jet=%.6f"
        % (
            args.ell,
            reduced.nrows(),
            args.seed,
            lambda1,
            ambient_gh(log2det, lambda1, reduced.nrows()),
            jet_gh(balanced, keys, scales, root, lambda1),
        )
    )


if __name__ == "__main__":
    main()
