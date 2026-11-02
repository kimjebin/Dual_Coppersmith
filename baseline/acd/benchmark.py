#!/usr/bin/env sage-python
"""Run one ACD benchmark cell and emit one JSON object."""

import argparse
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import mp
import mp_dual
from run_acd import ETA, RHO, gen_instance
from shared.reduction import selected_backend


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("t", type=int)
    parser.add_argument("m", type=int, help="multiplicity (called k in the code)")
    parser.add_argument("ell", type=int)
    parser.add_argument("method", choices=("dual", "primal"))
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--no-flatter", action="store_true",
                        help="legacy option; all reductions already use fpLLL")
    args = parser.parse_args()

    p, modulus, samples, gamma = gen_instance(
        args.t, args.m, args.ell, seed=args.seed
    )
    started = time.perf_counter()
    if args.method == "dual":
        result = mp_dual.attack(
            modulus, samples, RHO, t=args.t, k=args.m,
            use_flatter=False,
        )
        backend = selected_backend()
    else:
        result = mp.attack(modulus, samples, RHO, t=args.t, k=args.m)
        backend = selected_backend() + " + Groebner"

    recovered = result is not None and int(result[0]) == p
    record = {
        "eta": ETA,
        "gamma": gamma,
        "rho": RHO,
        "ell": args.ell,
        "t": args.t,
        "m": args.m,
        "method": args.method,
        "seed": args.seed,
        "backend": backend,
        "time_s": time.perf_counter() - started,
        "recovered": bool(recovered),
        "status": "ok" if recovered else "failed",
    }
    print(json.dumps(record, sort_keys=True), flush=True)
    if not recovered:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
