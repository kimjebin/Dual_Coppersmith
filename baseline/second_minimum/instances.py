"""Public ACD/HM bases for exact-minimum runs, without expensive integer SNF.

The affine monomial changes below are unimodular over ZZ. Flipping their
graded diagonal constraints constructs exactly the same N^m congruence kernel
as the attack scripts. Secret factors are used only to generate test instances.
"""

from math import comb, floor, log2

from sage.all import (ZZ, PolynomialRing, diagonal_matrix, gcd, inverse_mod,
                      matrix, prod, randint, random_prime, set_random_seed)

from shared.dual_moment import _exp_vec, _monomials_up_to


def _affine_kernel(ring, inverse_coordinates, multipliers, bound, degree):
    monomials = _monomials_up_to(ring, degree)
    keys = [_exp_vec(monomial) for monomial in monomials]
    inverse_rows = []
    for key in keys:
        polynomial = prod(inverse_coordinates[i] ** key[i] for i in range(len(key)))
        coefficients = {
            tuple(exponent) if hasattr(exponent, "__iter__") else (exponent,): value
            for exponent, value in polynomial.dict().items()
        }
        inverse_rows.append([coefficients.get(exponent, 0) for exponent in keys])
    inverse = matrix(ZZ, inverse_rows)
    scales = [ZZ(bound) ** (degree - sum(key)) for key in keys]
    kernel = diagonal_matrix(ZZ, [multipliers(key) for key in keys]) * inverse.transpose()
    return kernel * diagonal_matrix(ZZ, scales), keys


def acd_basis(t, m, ell, *, seed=0, eta=100, rho=10):
    if ell < 1 or m < 1 or t < m or eta < 3 or rho < 1:
        raise ValueError("Require ell>=1, t>=m>=1, eta>=3 and rho>=1")
    dimension = comb(t + ell, ell)
    loss = (dimension * log2(dimension) + dimension * dimension * log2(1.02)
            + rho * dimension * ell * t / (ell + 1))
    coefficient = comb(m + ell, ell) * m / (ell + 1)
    gamma = int(floor(0.9 * (m * eta * dimension - loss) / coefficient))
    if gamma < eta + 2:
        raise ValueError("The requested ACD parameters produce an invalid modulus range")
    set_random_seed(seed)
    p = random_prime(ZZ(2) ** eta - 1, lbound=ZZ(2) ** (eta - 1))
    q0 = randint(ZZ(2) ** (gamma - eta - 1), ZZ(2) ** (gamma - eta))
    modulus = p * q0
    bound = ZZ(2) ** rho
    samples, roots = [], []
    for _ in range(ell):
        qi = randint(ZZ(2) ** (gamma - eta - 1), ZZ(2) ** (gamma - eta))
        ri = randint(-int(bound) + 1, int(bound) - 1)
        samples.append(p * qi + ri)
        roots.append(int(ri))
    ring = PolynomialRing(ZZ, [f"x{i}" for i in range(ell)])
    inverse = [x + a for x, a in zip(ring.gens(), samples)]
    basis, keys = _affine_kernel(
        ring, inverse, lambda key: modulus ** min(sum(key), m), bound, t)
    return basis, {"application": "ACD", "t": t, "m": m, "ell": ell,
                   "seed": seed, "eta": eta, "rho": rho, "gamma": gamma,
                   "dim": basis.nrows(), "modulus": str(modulus),
                   "public_samples": list(map(str, samples)), "bound": str(bound),
                   "planted_root": roots, "monomial_keys": keys}


def hm_basis(ell, *, seed=0, pbits=380, qbits=20, xbits=6, degree=3):
    if ell < 1 or min(pbits, qbits) < 3 or xbits < 1 or degree < 1:
        raise ValueError("Invalid HM instance parameters")
    set_random_seed(seed)
    p = random_prime(ZZ(2) ** pbits - 1, lbound=ZZ(2) ** (pbits - 1))
    q = random_prime(ZZ(2) ** qbits - 1, lbound=ZZ(2) ** (qbits - 1))
    modulus = p * q
    bound = ZZ(2) ** xbits
    roots = [randint(-int(bound) + 1, int(bound) - 1) for _ in range(ell)]
    coefficients = [randint(1, int(modulus) - 1) for _ in range(ell)]
    while gcd(coefficients[0], modulus) != 1:
        coefficients[0] = randint(1, int(modulus) - 1)
    constant = -sum(coefficients[i] * roots[i] for i in range(ell))
    constant += p * randint(1, int(q))
    ring = PolynomialRing(ZZ, [f"x{i}" for i in range(ell)])
    xs = ring.gens()
    unit = inverse_mod(coefficients[0], modulus)
    normalized_constant = (constant * unit) % modulus
    normalized_coefficients = [(a * unit) % modulus for a in coefficients]
    inverse = [xs[0] - normalized_constant
               - sum(normalized_coefficients[i] * xs[i] for i in range(1, ell))]
    inverse.extend(xs[1:])
    basis, keys = _affine_kernel(
        ring, inverse, lambda key: modulus ** key[0], bound, degree)
    return basis, {"application": "HM", "t": degree, "m": degree, "ell": ell,
                   "seed": seed, "pbits": pbits, "qbits": qbits, "xbits": xbits,
                   "dim": basis.nrows(), "modulus": str(modulus),
                   "public_coefficients": list(map(str, coefficients)),
                   "public_constant": str(constant), "bound": str(bound),
                   "normalized_coefficients": list(map(str, normalized_coefficients)),
                   "normalized_constant": str(normalized_constant),
                   "planted_root": list(map(int, roots)), "monomial_keys": keys}
