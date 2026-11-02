"""
Bounded small-root recovery for  f(x_1, ..., x_l) == 0  (mod N),  |x_i| < b_i,
via the SNF modular-orthogonal ("dual") lattice + moment/Veronese decoding.

This is the Gröbner-free route: instead of reducing a bound-scaled *coefficient*
lattice, LLL-reducing it to short auxiliary polynomials, and solving those with a
Gröbner basis, we build the modular-orthogonal lattice of the shift coefficient
matrix and read the root's *monomial-evaluation* (moment) vector directly off the
reduced basis.

The core structure below follows Algorithm 1 line by line and is preserved exactly:

    Line 1  monomial set   M   = { x^a : |a| <= t }
    Line 2  index set      I   = { (j, b) : 1 <= j <= m, |b| + j*deg f <= t }
    Line 3  shift matrix   A    with rows  g_{j,b} = N^{m-j} x^b f^j  over M
    Line 4  dual kernel    L    = { z : A z == 0 (mod N^m) }          <-- SNF (flip divisors)
    Line 5  balance factor rho  = lcm{ b^a : x^a in M }
    Line 6  balanced basis  L_b = diag(rho / b^a) * L
    Line 7  reduce          reduce L_b                                <-- fpLLL:proved
    Lines 8-12  ratio-decode each reduced row  w  ->  r_i = b_i * w_{e_i} / w_0

The Line-7 reducer is Sage's proved fpLLL backend. Everything that defines
the method -- the SNF dual kernel, the moment balancing, and the ratio decoding --
is unchanged by that choice.
"""
import logging
from collections.abc import Sequence
from itertools import combinations_with_replacement
from itertools import product

from sage.all import ZZ, Zmod, matrix, diagonal_matrix, identity_matrix, lcm, prod

from .reduction import reduce_integer_basis


# --------------------------------------------------------------------------------------
# Line 4 core: modular-orthogonal (dual) kernel  L = { z in Z^n : A z == 0 (mod q) }.
#
# Two constructions that yield the SAME lattice (verified by HNF equality):
#   * snf_congruence_kernel  -- the SNF "flip the divisors" form (the professor's picture);
#     mathematically transparent but Sage's Smith form is computed over the *integers* on a
#     matrix whose entries are ~q = N^m, so it is slow at larger dimensions.
#   * howell_congruence_kernel -- the same kernel via the kernel over Z/qZ (Howell form,
#     entries bounded by q) lifted and completed with q*I.  Much faster: measured ~19x at
#     dim 31 and ~158x at dim 43 (SNF 58s -> 0.37s), identical lattice.
# congruence_kernel() dispatches; the fast Howell form is the default.
# --------------------------------------------------------------------------------------
def snf_congruence_kernel(A, modulus):
    """SNF form of the dual kernel: rows of diag(mult_i) * V^T, mult_i = q / gcd(q, d_i).

    From the Smith Normal Form  U A V = S = diag(d_1, .., d_r, 0, ..):  writing z = V y,
    A z == 0 (mod q)  <=>  d_i y_i == 0 (mod q), so y_i is a multiple of q / gcd(q, d_i)
    (a free coordinate d_i = 0 has multiplier 1) -- the "flipped" divisors.
    """
    q = ZZ(modulus)
    n = A.ncols()
    S, _U, V = A.smith_form()
    divs = [ZZ(S[i, i]) for i in range(min(A.nrows(), n))]
    mult = []
    for i in range(n):
        d = divs[i] if i < len(divs) else ZZ(0)
        mult.append(ZZ(1) if d == 0 else q // q.gcd(d))
    return diagonal_matrix(ZZ, mult) * V.transpose()


def howell_congruence_kernel(A, modulus):
    """Same dual kernel L = { z : A z == 0 (mod q) }, computed the fast way.

    Take a basis of ker(A) over Z/qZ (Sage uses the Howell form; entries stay < q), lift it
    to Z, and complete it with q*I (the lattice always contains q*Z^n); the Hermite form of
    that stack is a Z-basis of L.  Gives the identical lattice as snf_congruence_kernel but
    avoids an integer Smith form over ~q-sized entries, so it is dramatically faster.
    """
    q = ZZ(modulus)
    n = A.ncols()
    Kq = A.change_ring(Zmod(q)).right_kernel_matrix()
    B = Kq.change_ring(ZZ).stack(q * identity_matrix(ZZ, n))
    return B.hermite_form(include_zero_rows=False)


def congruence_kernel(A, modulus, method="howell"):
    """Dual kernel basis; method in {"howell" (fast, default), "snf" (SNF flip form)}."""
    if method == "snf":
        return snf_congruence_kernel(A, modulus)
    if method == "howell":
        return howell_congruence_kernel(A, modulus)
    raise ValueError("unknown kernel method %r" % method)


# --------------------------------------------------------------------------------------
# Line 7 reducer: use the same fpLLL backend as the primal implementation.
# --------------------------------------------------------------------------------------
def _reduce(basis, use_flatter=False):
    return reduce_integer_basis(basis, use_flatter=use_flatter)


def _exp_vec(mono):
    """Exponent tuple of a monomial, uniform for univariate and multivariate rings."""
    ex = mono.exponents()[0]
    return tuple(ex) if hasattr(ex, "__iter__") else (ex,)


def _monomials_up_to(pr, t):
    """Line 1 / Line 2 helper: every monomial of total degree <= t, constant first."""
    gens = pr.gens()
    monos = []
    for deg in range(t + 1):
        for combo in combinations_with_replacement(range(len(gens)), deg):
            mono = pr.one()
            for idx in combo:
                mono *= gens[idx]
            monos.append(mono)
    return monos


def small_roots(f, N, bounds, m=1, t=None, use_flatter=False, kernel_method="howell",
                divisor_mode=False, scan_rows=False):
    """Recover bounded common roots of one polynomial or a polynomial sequence.

    :param f: a polynomial, or a nonempty sequence of polynomials, in
        ZZ[x_1, ..., x_l].
    :param N: the (known) modulus, or a composite with a hidden divisor (see divisor_mode).
    :param bounds: sequence (b_1, ..., b_l) with |r_i| < b_i.
    :param m: modulus lifting exponent (shifts vanish modulo N^m at the root).
    :param t: monomial cutoff (defaults to m times the maximum relation degree).
    :param use_flatter: legacy compatibility argument; fpLLL is always used.
    :param divisor_mode: if True, the modulus is an unknown divisor p of N (the
        Herrmann-May / factoring-with-known-bits setting): the shifts vanish modulo
        p^m at the root and the public N^m is used as a proxy, so a candidate root is
        accepted when all relation values share a nontrivial factor with N rather
        than when every relation vanishes modulo N.
    :param scan_rows: decode every reduced row when True. The paper algorithm uses
        only the first reduced row, which is the default.
    :return: list of recovered roots, each a tuple (r_1, ..., r_l).
    """
    relations = list(f) if isinstance(f, Sequence) else [f]
    if not relations:
        raise ValueError("at least one polynomial relation is required")
    pr = relations[0].parent()
    if any(relation.parent() != pr for relation in relations):
        raise ValueError("all polynomial relations must have the same parent")
    if any(relation.is_zero() for relation in relations):
        raise ValueError("polynomial relations must be nonzero")
    x = pr.gens()
    ell = len(x)
    assert len(bounds) == ell, "one bound per variable"
    b = [ZZ(bi) for bi in bounds]
    N = ZZ(N)
    max_degree = max(relation.degree() for relation in relations)
    if t is None:
        t = m * max_degree
    if t < max_degree:
        raise ValueError("t must be at least the maximum relation degree")

    # -- Line 1: monomial set M, indexed by exponent tuple. ---------------------------
    # NOTE: columns are keyed by the exponent tuple, NOT the monomial object, because for
    # univariate Sage polynomials g.coefficients() and g.monomials() use OPPOSITE degree
    # orderings; g.dict() (exponent -> coeff) avoids that pitfall for both ring types.
    monos = _monomials_up_to(pr, t)
    keys = [_exp_vec(mono) for mono in monos]
    idx = {k: i for i, k in enumerate(keys)}
    n = len(monos)

    def bpow_key(k):
        return prod(b[i] ** ZZ(k[i]) for i in range(ell))

    # -- Lines 2-3: total-degree ideal-power shifts. -----------------------------------
    # For every relation multi-index i with 1 <= |i| <= m, include
    # N^(m-|i|) x^beta prod_a f_a^i_a whenever its total degree is at most t.
    # The |i|=0 rows are zero modulo N^m and impose no kernel constraint.
    powers = []
    for relation in relations:
        relation_powers = [pr.one()]
        for _ in range(m):
            relation_powers.append(relation_powers[-1] * relation)
        powers.append(relation_powers)

    rows = []
    relation_count = len(relations)
    for total_power in range(1, m + 1):
        for relation_exponents in product(range(total_power + 1), repeat=relation_count):
            if sum(relation_exponents) != total_power:
                continue
            ideal_power = prod(
                powers[a][relation_exponents[a]] for a in range(relation_count)
            )
            remaining_degree = t - ideal_power.degree()
            if remaining_degree < 0:
                continue
            for beta in _monomials_up_to(pr, remaining_degree):
                g = (N ** (m - total_power)) * beta * ideal_power
                row = [0] * n
                for e, coeff in g.dict().items():
                    k = tuple(e) if hasattr(e, "__iter__") else (e,)
                    row[idx[k]] = coeff
                rows.append(row)
    A = matrix(ZZ, rows)

    # -- Line 4: dual kernel L = { z : A z == 0 (mod N^m) }. ---------------------------
    L = congruence_kernel(A, N ** m, method=kernel_method)

    # -- Line 5: balancing factor rho = lcm{ b^a }. -----------------------------------
    rho = lcm([bpow_key(k) for k in keys])

    # -- Line 6: balanced dual basis L_b = diag(rho / b^a) * L. ------------------------
    scales = [rho // bpow_key(k) for k in keys]
    LB = L * diagonal_matrix(ZZ, scales)

    # -- Line 7: reduce with fpLLL. -------------------------------------------------
    LB = _reduce(LB, use_flatter=use_flatter)

    # -- Lines 8-12: ratio-decode each reduced row. -----------------------------------
    const_key = (0,) * ell
    unit_key = [tuple(1 if a == i else 0 for a in range(ell)) for i in range(ell)]
    seen = set()
    results = []
    reduced_rows = LB.rows() if scan_rows else [LB.row(0)]
    for w in reduced_rows:
        w0 = ZZ(w[idx[const_key]])      # Line 9: constant coordinate, must be nonzero
        if w0 == 0:
            continue
        r = []                          # Line 10: r_i = b_i * w_{e_i} / w_0
        ok = True
        for i in range(ell):
            num = b[i] * ZZ(w[idx[unit_key[i]]])
            if num % w0 != 0:           # integrality
                ok = False
                break
            r.append(num // w0)
        if not ok:
            continue
        # -- Line 11: bounds, Veronese consistency, and the original congruence. -------
        if any(abs(r[i]) >= b[i] for i in range(ell)):
            continue
        # w must really be a scalar multiple of the moment vector: w_a * b^a == w_0 * r^a.
        if not all(
            ZZ(w[i]) * bpow_key(keys[i]) == w0 * prod(ZZ(r[a]) ** keys[i][a] for a in range(ell))
            for i in range(n)
        ):
            continue
        if divisor_mode:
            # Hidden-divisor setting: p | f(r), so a nontrivial gcd with N reveals a factor.
            g = N
            for relation in relations:
                g = g.gcd(ZZ(relation(*r)))
            if not (1 < g < N):
                continue
        elif any(relation(*r) % N != 0 for relation in relations):
            continue
        root = tuple(int(v) for v in r)     # Line 12
        if root not in seen:
            seen.add(root)
            results.append(root)
    return results
