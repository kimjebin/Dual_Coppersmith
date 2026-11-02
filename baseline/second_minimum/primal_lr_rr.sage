# primal_lr_rr.sage
# -----------------------------------------------------------------------------
# Primal lattice-reduction (LR) timing for the "Primal (LR)" / "Primal (RR)"
# columns of the Section-4 tables.  The primal method runs in two stages:
#   LR = lattice reduction (LLL), then  RR = Groebner root recovery.
# This script measures the LR stage (Sage .LLL(), single-threaded wall time)
# for every table row, capped at the 600 s budget.  RR is then the primal total
# minus LR: RR = t_primal - LR on completing rows, and t/o where the primal
# baseline times out; if LR itself reaches 600 s (HM ell=9) the reduction alone
# exhausts the budget and recovery never runs.
#
# The primal totals t_primal come from the per-application baselines
# (acd/mp.py, shared/small_roots for HM); this
# script reproduces only the LR half, which the tables report separately.
#
# Run:  FLATTER unused here (LLL only).  sage primal_lr_rr.sage
# -----------------------------------------------------------------------------
import os, sys, time
from itertools import product
from math import comb, gcd, log2

os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                  MKL_NUM_THREADS="1", SAGE_NUM_THREADS="1")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "acd"))
from sage.all import (ZZ, PolynomialRing, random_prime, randint, set_random_seed,
                      inverse_mod, prod, alarm, cancel_alarm, AlarmInterrupt)
from shared import small_roots

CAP = 600  # per-row budget (seconds); matches the paper timeout


def time_lll(B):
    """LLL reduction wall time. Returns (seconds, timed_out)."""
    t0 = time.perf_counter()
    try:
        alarm(CAP); B.LLL(delta=0.8, algorithm="fpLLL:proved"); cancel_alarm()
        return time.perf_counter() - t0, False
    except AlarmInterrupt:
        cancel_alarm()
        return time.perf_counter() - t0, True


def acd_lattice(t, k, ell, eta=100, rho=10, seed=0):
    """Bound-scaled ACD coefficient lattice (same construction as acd/mp.py)."""
    set_random_seed(seed)
    d = comb(t + ell, ell)
    lhs = d * log2(d) + d * d * log2(1.02) + rho * d * (ell * t / (ell + 1))
    cg = comb(k + ell, ell) * k / (ell + 1)
    gamma = int(0.9 * (k * eta * d - lhs) / cg)
    p = random_prime(ZZ(2) ** eta - 1, lbound=ZZ(2) ** (eta - 1))
    q0 = randint(ZZ(2) ** (gamma - eta - 1), ZZ(2) ** (gamma - eta)); N = int(p * q0)
    a = [int(p * randint(ZZ(2) ** (gamma - eta - 1), ZZ(2) ** (gamma - eta))
             + randint(-(2 ** rho) + 1, 2 ** rho - 1)) for _ in range(ell)]
    pr = ZZ[tuple("x%d" % i for i in range(ell))]; x = pr.gens()
    shifts = []
    for idx in product(*[range(t + 1) for _ in x]):
        if sum(idx) <= t:
            fi = N ** max(k - sum(idx), 0)
            for m in range(len(idx)):
                fi *= (x[m] - a[m]) ** idx[m]
            shifts.append(fi)
    return small_roots.create_lattice(pr, shifts, [2 ** rho] * ell)[0]


def hm_lattice(ell, pbits=380, qbits=20, xbits=6, t=3, seed=0):
    """Bound-scaled Herrmann-May coefficient lattice (t=m=3)."""
    set_random_seed(seed)
    p = random_prime(ZZ(2) ** pbits - 1, lbound=ZZ(2) ** (pbits - 1))
    q = random_prime(ZZ(2) ** qbits - 1, lbound=ZZ(2) ** (qbits - 1)); N = int(p * q)
    X = ZZ(2) ** xbits
    P = PolynomialRing(ZZ, ["x%d" % i for i in range(ell)]); xs = P.gens()
    r = [randint(-int(X) + 1, int(X) - 1) for _ in range(ell)]
    a = [randint(1, N - 1) for _ in range(ell)]
    while gcd(a[0], N) != 1:
        a[0] = randint(1, N - 1)
    a0 = -sum(a[i] * r[i] for i in range(ell)) + int(p) * randint(1, int(q))
    f = a0 + sum(a[i] * xs[i] for i in range(ell))
    fm = P({e: (ZZ(c) * inverse_mod(ZZ(a[0]), N)) % N
            for e, c in zip(f.exponents(), f.coefficients())})
    fpw = [P.one()]
    for _ in range(t):
        fpw.append(fpw[-1] * fm)
    shifts = []
    for kk in range(t + 1):
        for iv in product(*[range(t - kk + 1) for _ in range(ell - 1)]):
            if sum(iv) <= t - kk:
                shifts.append((ZZ(N) ** (t - kk)) * fpw[kk]
                              * prod(xs[1 + j] ** iv[j] for j in range(ell - 1)))
    return small_roots.create_lattice(P, shifts, [X] * ell)[0]


def main():
    sys.stdout.reconfigure(line_buffering=True)
    print("app | row | dim | LR_s | timed_out")
    for t, k, ell in [(2, 1, 10), (2, 1, 11), (2, 1, 12), (2, 1, 13),
                      (2, 2, 10), (2, 2, 11), (2, 2, 12),
                      (3, 2, 5), (3, 2, 6), (3, 2, 7),
                      (3, 3, 5), (3, 3, 6), (3, 3, 7)]:
        B = acd_lattice(t, k, ell); dt, to = time_lll(B)
        print("ACD (%d,%d,%d) | %d | %.2f | %s" % (t, k, ell, B.ncols(), dt, to))
    for ell in range(2, 10):
        B = hm_lattice(ell); dt, to = time_lll(B)
        print("HM l=%d | %d | %.2f | %s" % (ell, B.ncols(), dt, to))


main()
