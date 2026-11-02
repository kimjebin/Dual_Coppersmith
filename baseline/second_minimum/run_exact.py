#!/usr/bin/env sage-python
"""Run exact ACD/HM successive minima with parallel, resumable GMP workers."""

import argparse
import json
import os
from pathlib import Path
import sys
import time

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(name, "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from shared.exact_minima import (exact_successive_minima, load_minima_input,
                                prepare_successive_minima, save_minima_input)  # noqa: E402
from second_minimum.instances import acd_basis, hm_basis  # noqa: E402


ACD_GRID = [(2, 1, ell) for ell in range(10, 14)] + [
    (2, 2, ell) for ell in range(10, 13)] + [
    (3, 2, ell) for ell in range(5, 8)] + [
    (3, 3, ell) for ell in range(5, 8)]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("application", choices=("acd", "hm", "all"))
    parser.add_argument("--t", type=int, default=3)
    parser.add_argument("--m", type=int, default=3)
    parser.add_argument("--ell", type=int)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--workers", type=int, default=0,
                        help="native worker count; 0 selects a CPU/RAM-aware maximum")
    parser.add_argument("--memory-gb", type=float, default=0,
                        help="search RSS budget; 0 uses 80%% of available memory")
    parser.add_argument("--max-seconds", type=float, default=0,
                        help="search time per invocation and per case, excluding preparation")
    parser.add_argument("--max-nodes", type=int, default=0,
                        help="search node budget per invocation and per case")
    parser.add_argument("--chunk-nodes", type=int, default=20000)
    parser.add_argument("--chunk-seconds", type=float, default=2.0)
    parser.add_argument("--progress-seconds", type=float, default=10.0)
    parser.add_argument("--quotient-delta", type=float, default=0.99)
    parser.add_argument("--prepare-only", action="store_true",
                        help="save both LLL stages and portable preparation; never enumerate")
    parser.add_argument("--reuse-from", type=Path,
                        help="legacy campaign directory for verified warm preparation")
    parser.add_argument("--out", type=Path, required=True,
                        help="checkpoint directory; repeat the same command to resume")
    parser.add_argument("--eta", type=int, default=100)
    parser.add_argument("--rho", type=int, default=10)
    parser.add_argument("--pbits", type=int, default=380)
    parser.add_argument("--qbits", type=int, default=20)
    parser.add_argument("--xbits", type=int, default=6)
    args = parser.parse_args(argv)
    if args.application == "all" and args.ell is not None:
        parser.error("--ell cannot select one row of the all-paper-rows run")
    if args.application == "acd" and args.ell is None:
        parser.error("acd requires --ell; use all for the complete 21-row paper grid")
    cases = []
    if args.application in ("acd", "all"):
        grid = ACD_GRID if args.application == "all" else [(args.t, args.m, args.ell)]
        cases.extend(("acd", t, m, ell) for t, m, ell in grid)
    if args.application in ("hm", "all"):
        grid = range(2, 10) if args.ell is None else [args.ell]
        degree = 3 if args.application == "all" else args.t
        cases.extend(("hm", degree, degree, ell) for ell in grid)
    incomplete = False
    for app, t, m, ell in cases:
        case = f"{app}-t{t}-m{m}-ell{ell}-seed{args.seed}"
        print(f"CASE {case}: building public lattice", flush=True)
        began = time.monotonic()
        case_folder = args.out / case
        if (case_folder / "input.json").exists():
            basis, instance = load_minima_input(case_folder)
            requested = dict(application=app.upper(), t=t, m=m, ell=ell, seed=args.seed)
            requested.update(dict(eta=args.eta, rho=args.rho) if app == "acd" else
                             dict(pbits=args.pbits, qbits=args.qbits, xbits=args.xbits))
            if any(instance.get(key) != value for key, value in requested.items()):
                raise ValueError("Saved public instance parameters differ from CLI request")
        elif app == "acd":
            basis, instance = acd_basis(t, m, ell, seed=args.seed, eta=args.eta, rho=args.rho)
        else:
            basis, instance = hm_basis(ell, seed=args.seed, pbits=args.pbits,
                                      qbits=args.qbits, xbits=args.xbits, degree=t)
        build_seconds = time.monotonic() - began
        case_folder.mkdir(parents=True, exist_ok=True)
        instance.update(lattice_dimension=basis.nrows(), ambient_columns=basis.ncols(),
                        quotient_dimension=basis.nrows() - 1)
        temporary = case_folder / "instance.json.tmp"
        with temporary.open("w") as handle:
            json.dump(instance, handle, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, case_folder / "instance.json")
        if args.prepare_only:
            result = prepare_successive_minima(
                basis, case_folder, instance=instance, workers=args.workers,
                quotient_delta=args.quotient_delta,
                reuse=args.reuse_from / case if args.reuse_from else None)
            print("PREPARED_RESULT " + json.dumps({key: result[key] for key in
                  ("preparation_complete", "lattice_dimension", "enumeration_dimension",
                   "ambient_lll_seconds", "quotient_lll_seconds", "legacy_native_reused")},
                  sort_keys=True), flush=True)
            continue
        save_minima_input(basis, case_folder, instance)
        result = exact_successive_minima(
            basis, case_folder, workers=args.workers, memory_gb=args.memory_gb,
            max_seconds=args.max_seconds, max_nodes=args.max_nodes,
            chunk_nodes=args.chunk_nodes, chunk_seconds=args.chunk_seconds,
            progress_seconds=args.progress_seconds, quotient_delta=args.quotient_delta)
        result.update(instance=instance, build_seconds=build_seconds,
                      total_invocation_seconds=time.monotonic() - began)
        incomplete |= not result["exact"]
        # Keep the runner record separate from the engine's certified result.
        path = args.out / case / "row.json"
        temporary = path.with_suffix(".tmp")
        with temporary.open("w") as handle:
            json.dump(result, handle, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
        summary = {key: result.get(key) for key in
                   ("exact", "status", "lambda1_sq", "lambda2_sq", "log2_lambda1", "log2_lambda2")}
        summary.update(case=case, checkpoint=str(path.parent),
                       total_seconds=result["total_invocation_seconds"])
        print("EXACT_RESULT " + json.dumps(summary, sort_keys=True), flush=True)
    return 2 if incomplete else 0


if __name__ == "__main__":
    raise SystemExit(main())
