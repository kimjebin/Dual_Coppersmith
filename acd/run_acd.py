#!/usr/bin/env python
"""Minimal ACD example: primal (Groebner) vs dual (ours) on one instance.

Reproduces one row of Table 1 (Section 4.1): the primal Coppersmith method
(LLL + Groebner root recovery) versus the dual/moment method of this paper,
which reads the root off a short vector of the modular-orthogonal lattice with
no Groebner basis.

Run from the repository root:

    sage -python acd/run_acd.py            # default (t,m,ell) = (2,1,10)
    sage -python acd/run_acd.py 3 2 6      # (t,m,ell) = (3,2,6)

Both methods use Sage's proved fpLLL backend with the same LLL delta.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root (shared/)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))                   # acd/ (mp, mp_dual)

from math import comb

from sage.all import ZZ, random_prime, randint, set_random_seed

import mp
import mp_dual

ETA = 100
RHO = 10


def gamma_max(ell, t, k):
    """Galbraith--Gebregiyorgis--Murphy solvability bound (used at 90%)."""
    from math import log2

    d = comb(t + ell, ell)
    lhs = d * log2(d) + d * d * log2(1.02) + RHO * d * (ell * t / (ell + 1))
    cg = comb(k + ell, ell) * k / (ell + 1)
    return (k * ETA * d - lhs) / cg


def gen_instance(t, k, ell, seed=0):
    set_random_seed(seed)
    gamma = int(0.9 * gamma_max(ell, t, k))
    p = random_prime(ZZ(2) ** ETA - 1, lbound=ZZ(2) ** (ETA - 1))
    q0 = randint(ZZ(2) ** (gamma - ETA - 1), ZZ(2) ** (gamma - ETA))
    N = p * q0
    a = []
    for _ in range(ell):
        qi = randint(ZZ(2) ** (gamma - ETA - 1), ZZ(2) ** (gamma - ETA))
        ri = randint(-(2 ** RHO) + 1, 2 ** RHO - 1)
        a.append(p * qi + ri)
    return int(p), int(N), [int(x) for x in a], gamma


def _recovered(result, p):
    return result is not None and int(result[0]) == p


def main():
    args = sys.argv[1:]
    t, k, ell = (int(args[0]), int(args[1]), int(args[2])) if len(args) == 3 else (2, 1, 10)
    p, N, a, gamma = gen_instance(t, k, ell)
    dim = comb(ell + t, ell)
    print(f"ACD  eta={ETA} rho={RHO}  (t,m)=({t},{k})  ell={ell}  dim={dim}  gamma={gamma}")

    start = time.time()
    dual = mp_dual.attack(N, a, RHO, t=t, k=k, use_flatter=False)
    t_dual = time.time() - start
    print(f"  dual  (ours):  {'OK ' if _recovered(dual, p) else 'FAIL'}  {t_dual:8.2f} s")

    start = time.time()
    primal = mp.attack(N, a, RHO, t=t, k=k)
    t_primal = time.time() - start
    print(f"  primal      :  {'OK ' if _recovered(primal, p) else 'FAIL'}  {t_primal:8.2f} s")


if __name__ == "__main__":
    main()
