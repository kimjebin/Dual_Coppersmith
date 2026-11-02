# shift_style_audit.sage
# -----------------------------------------------------------------------------
# Audit for Section 4.1 (ACD): the pure-power ACD shift presentation used by
# mp_dual.py and the usual total-degree Coppersmith presentation (with x^beta
# monomial multipliers) define the SAME modular orthogonal (dual kernel) lattice.
#
#   (1) pure_power  (mp_dual._build_shifts, what the paper/experiments use):
#         N^{max(m-|i|,0)} prod_j (x_j-a_j)^{i_j},          |i| <= t
#   (2) usual_total (standard Coppersmith, with x^beta multipliers):
#         N^{m-|u|} x^beta prod_j (x_j-a_j)^{u_j},           1<=|u|<=m, |beta|+|u|<=t
#
# Here m is the multiplicity (mp_dual calls it k).  Both presentations are shift
# families for the same instance; the dual kernel is L = { v : A v == 0 (mod N^m) }
# on the common monomial column set {x^alpha : |alpha| <= t}.
#
# Finding (all tested (ell,t,m)): the two kernels are EQUAL as lattices (identical
# Hermite normal form and determinant), even though usual_total has many more
# (redundant) rows.  Hence the recovered root and lambda_1/lambda_2 are identical;
# only the raw basis, and therefore the reduction runtime, differ.  Runtime (LLL,
# `timing` mode): comparable at t=2 (~1x) but several-fold slower at t=3
# (2-4x, and usual_total times out at (3,2,7) where pure_power finishes), because
# the extra redundant rows yield a worse-conditioned SNF basis of the same-
# dimension lattice.  Mechanism for the kernel equality: x^beta y^u is an integer
# combination of pure-power rows modulo N^m, so its extra congruences are redundant.
#
# Reuses the production functions of acd/mp_dual.py; nothing in the package is
# modified.  Run:  sage acd/shift_style_audit.sage
# -----------------------------------------------------------------------------
import os, sys, time
from itertools import product
from math import comb, log2

os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                  MKL_NUM_THREADS="1", SAGE_NUM_THREADS="1")
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)                      # dual_coppersmith/
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)                          # so `import mp_dual` resolves

from sage.all import (ZZ, PolynomialRing, random_prime, randint, set_random_seed,
                      diagonal_matrix, matrix, prod,
                      alarm, cancel_alarm, AlarmInterrupt)
from shared import small_roots
from shared.reduction import reduce_integer_basis
from mp_dual import snf_congruence_kernel, _moment_roots


# ---- ACD instance in the solvable regime (same sizing as second_minimum) -------
def gen_acd(ell, t, k, eta=100, rho=10, seed=0):
    set_random_seed(seed)
    d = comb(t + ell, ell)
    lhs = d * log2(d) + d * d * log2(1.02) + rho * d * (ell * t / (ell + 1))
    cg = comb(k + ell, ell) * k / (ell + 1)
    gamma = int(0.9 * (k * eta * d - lhs) / cg)
    p = random_prime(ZZ(2) ** eta - 1, lbound=ZZ(2) ** (eta - 1))
    q0 = randint(ZZ(2) ** (gamma - eta - 1), ZZ(2) ** (gamma - eta))
    N = ZZ(p) * q0
    a, r = [], []
    for _ in range(ell):
        ri = randint(-(2 ** rho) + 1, 2 ** rho - 1)
        a.append(int(p * randint(ZZ(2) ** (gamma - eta - 1), ZZ(2) ** (gamma - eta)) + ri))
        r.append(int(ri))
    return N, a, r, int(p)


def monos_le_t(x, t):
    ell = len(x)
    return [prod(x[j] ** e[j] for j in range(ell))
            for e in product(*[range(t + 1) for _ in x]) if sum(e) <= t]


def build_pure(x, N, a, t, k):
    """N^{max(k-|i|,0)} prod_j (x_j-a_j)^{i_j}, |i| <= t  (the used presentation)."""
    ell = len(x)
    shifts = []
    for i in product(*[range(t + 1) for _ in x]):
        if sum(i) <= t:
            g = ZZ(N) ** max(k - sum(i), 0)
            for j in range(ell):
                g *= (x[j] - a[j]) ** i[j]
            shifts.append(g)
    return shifts


def build_usual(x, N, a, t, k):
    """N^{k-|u|} x^beta prod_j (x_j-a_j)^{u_j}, 1<=|u|<=k, |beta|+|u|<=t."""
    ell = len(x)
    y = [x[j] - a[j] for j in range(ell)]
    shifts = []
    for u in product(*[range(k + 1) for _ in x]):
        du = sum(u)
        if du < 1 or du > k:
            continue
        yu = prod(y[j] ** u[j] for j in range(ell))
        Npow = ZZ(N) ** (k - du)
        for beta in product(*[range(t + 1) for _ in x]):
            if sum(beta) + du <= t:
                xb = prod(x[j] ** beta[j] for j in range(ell))
                shifts.append(Npow * xb * yu)
    return shifts


def coeff_matrix(shifts, monos, x):
    ell = len(x)
    idx = {m: i for i, m in enumerate(monos)}
    rows = []
    for g in shifts:
        row = [0] * len(monos)
        for e, c in g.dict().items():
            row[idx[prod(x[j] ** int(e[j]) for j in range(ell))]] += ZZ(c)
        rows.append(row)
    return matrix(ZZ, rows)


def lattice_det(H):
    return abs(prod(H[i, i] for i in range(H.nrows())))


def compare_kernels(ell, t, k, seed=0):
    N, a, r, p = gen_acd(ell, t, k, seed=seed)
    pr = PolynomialRing(ZZ, ['x%d' % i for i in range(ell)])
    x = pr.gens()
    monos = monos_le_t(x, t)
    pure, usual = build_pure(x, N, a, t, k), build_usual(x, N, a, t, k)
    Ap = coeff_matrix(pure, monos, x)
    Au = coeff_matrix(usual, monos, x)
    Hp = snf_congruence_kernel(Ap, ZZ(N) ** k).hermite_form(include_zero_rows=False)
    Hu = snf_congruence_kernel(Au, ZZ(N) ** k).hermite_form(include_zero_rows=False)
    equal = (Hp == Hu)
    print("[kernel] ell=%d t=%d k=%d  n=%d | rows pure/usual=%d/%d  rank=%d/%d | "
          "det match=%s  KERNEL EQUAL=%s"
          % (ell, t, k, len(monos), len(pure), len(usual), Ap.rank(), Au.rank(),
             lattice_det(Hp) == lattice_det(Hu), equal))
    return equal


def run_recovery(shifts, monos, x, N, r, rho, t, k):
    A = coeff_matrix(shifts, monos, x)
    L = snf_congruence_kernel(A, ZZ(N) ** k)
    R = 2 ** rho
    LB = L * diagonal_matrix(ZZ, [ZZ(R) ** (t - mo.degree()) for mo in monos])
    t0 = time.perf_counter()
    LBr, used = reduce_integer_basis(LB, use_flatter=False, return_backend=True)
    dt = time.perf_counter() - t0
    ok = any([roots[xi] for xi in x] == list(r)
             and all(-R < roots[xi] < R for xi in x)
             for roots in _moment_roots(LBr, monos, x, R, t))
    return dt, ok, used


def compare_recovery(ell, t, k, seed=0):
    N, a, r, p = gen_acd(ell, t, k, seed=seed)
    pr = PolynomialRing(ZZ, ['x%d' % i for i in range(ell)])
    x = pr.gens()
    monos = monos_le_t(x, t)
    pure, usual = build_pure(x, N, a, t, k), build_usual(x, N, a, t, k)
    dp, okp, up = run_recovery(pure, monos, x, N, r, 10, t, k)
    du, oku, uu = run_recovery(usual, monos, x, N, r, 10, t, k)
    print("[recover] ell=%d t=%d k=%d  rows=%d/%d | pure %.3fs ok=%s | usual %.3fs ok=%s"
          " (reducer=%s)" % (ell, t, k, len(pure), len(usual), dp, okp, du, oku, up))


def _reduce(LB):
    return reduce_integer_basis(LB, use_flatter=False, return_backend=True)


def time_pipeline(shifts, monos, x, N, r, rho, t, k, cap=600):
    """Full dual pipeline (kernel + balance + reduce + decode) under a wall cap.
    Wall time includes the proved fpLLL reduction."""
    import time
    R = 2 ** rho
    w0 = time.perf_counter()
    try:
        alarm(cap)
        A = coeff_matrix(shifts, monos, x)
        L = snf_congruence_kernel(A, ZZ(N) ** k)
        t_ker = time.perf_counter() - w0
        LB = L * diagonal_matrix(ZZ, [ZZ(R) ** (t - mo.degree()) for mo in monos])
        w1 = time.perf_counter()
        LBr, red = _reduce(LB)
        t_red = time.perf_counter() - w1
        ok = any([roots[xi] for xi in x] == list(r)
                 and all(-R < roots[xi] < R for xi in x)
                 for roots in _moment_roots(LBr, monos, x, R, t))
        cancel_alarm()
        return dict(to=False, ker=t_ker, red=t_red,
                    tot=time.perf_counter() - w0, ok=ok, reducer=red)
    except AlarmInterrupt:
        cancel_alarm()
        return dict(to=True, ker=None, red=None, tot=float(cap), ok=False, reducer="?")


def timing_table(cap=600):
    """pure_power vs usual_total full-pipeline wall time on the tab:acd params."""
    TABLE = [(2, 1, 10), (2, 1, 11), (2, 1, 12), (2, 1, 13),
             (2, 2, 10), (2, 2, 11), (2, 2, 12),
             (3, 2, 5), (3, 2, 6), (3, 2, 7),
             (3, 3, 5), (3, 3, 6), (3, 3, 7)]
    hdr = ("%-11s %4s %11s | %-22s | %-22s | %6s | %s"
           % ("(t,m) ell", "dim", "rows P/U",
              "pure ker+red=tot (s)", "usual ker+red=tot (s)", "ratio", "ok P/U"))
    print(hdr)
    print("-" * len(hdr))
    for (t, m, ell) in TABLE:
        N, a, r, p = gen_acd(ell, t, m)
        pr = PolynomialRing(ZZ, ['x%d' % i for i in range(ell)])
        x = pr.gens()
        monos = monos_le_t(x, t)
        pure, usual = build_pure(x, N, a, t, m), build_usual(x, N, a, t, m)
        rp = time_pipeline(pure, monos, x, N, r, 10, t, m, cap=cap)
        ru = time_pipeline(usual, monos, x, N, r, 10, t, m, cap=cap)

        def fmt(d):
            return ("t/o" if d["to"]
                    else "%.2f+%.2f=%.2f" % (d["ker"], d["red"], d["tot"]))
        ratio = ("%.2fx" % (ru["tot"] / rp["tot"])
                 if not rp["to"] and not ru["to"] and rp["tot"] > 0 else "-")
        print("(%d,%d) l=%-3d %4d %11s | %-22s | %-22s | %6s | %s/%s"
              % (t, m, ell, len(monos), "%d/%d" % (len(pure), len(usual)),
                 fmt(rp), fmt(ru), ratio, rp["ok"], ru["ok"]))


if __name__ == "__main__":
    sys.stdout.reconfigure(line_buffering=True)
    mode = sys.argv[1] if len(sys.argv) > 1 else "check"
    if mode == "timing":
        print("=== TIMING: full dual pipeline, tab:acd params, 600s wall cap ===")
        timing_table(cap=600)
    else:
        print("=== KERNEL EQUALITY (pure_power vs usual_total) ===")
        for prm in [(2, 2, 1), (2, 2, 2), (2, 3, 2), (2, 3, 3),
                    (3, 2, 1), (3, 2, 2), (3, 3, 2)]:
            compare_kernels(*prm)
        print("\n=== RECOVERY + REDUCTION RUNTIME ===")
        for prm in [(2, 2, 2), (3, 2, 2), (3, 3, 2), (5, 3, 2)]:
            compare_recovery(*prm)
