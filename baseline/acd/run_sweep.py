#!/usr/bin/env python3
"""Run the 13 ACD paper configurations with isolated per-cell timeouts."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import time


ROOT = Path(__file__).resolve().parents[1]
WORKER = Path(__file__).with_name("benchmark.py")
CONFIGS = [
    *((2, 1, ell) for ell in range(10, 14)),
    *((2, 2, ell) for ell in range(10, 13)),
    *((3, 2, ell) for ell in range(5, 8)),
    *((3, 3, ell) for ell in range(5, 8)),
]


def run_cell(t, m, ell, method, seed, cap):
    command = [
        "sage", "-python", str(WORKER), str(t), str(m), str(ell), method,
        "--seed", str(seed),
    ]
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", SAGE_NUM_THREADS="1")
    started = time.perf_counter()
    try:
        proc = subprocess.run(
            command, cwd=ROOT, env=env, text=True, capture_output=True,
            timeout=cap, check=False,
        )
    except subprocess.TimeoutExpired:
        return {"t": t, "m": m, "ell": ell, "method": method, "seed": seed,
                "time_s": cap, "recovered": False, "status": "timeout"}

    lines = [line for line in proc.stdout.splitlines() if line.startswith("{")]
    if lines:
        record = json.loads(lines[-1])
        if proc.returncode and record.get("status") == "ok":
            record["status"] = "error"
        return record
    return {"t": t, "m": m, "ell": ell, "method": method, "seed": seed,
            "time_s": time.perf_counter() - started, "recovered": False,
            "status": "error", "stderr_tail": proc.stderr[-500:]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", choices=("dual", "primal", "both"), default="both")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--cap", type=float, default=600.0)
    parser.add_argument(
        "--out", type=Path,
        default=ROOT / "acd" / "results" / "raw" / "acd_rerun.jsonl",
    )
    parser.add_argument("--force", action="store_true", help="overwrite --out")
    args = parser.parse_args()
    if args.out.exists() and not args.force:
        parser.error(f"{args.out} exists; choose another path or pass --force")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    methods = ("dual", "primal") if args.method == "both" else (args.method,)

    with args.out.open("w", encoding="utf-8") as handle:
        for t, m, ell in CONFIGS:
            for method in methods:
                record = run_cell(t, m, ell, method, args.seed, args.cap)
                handle.write(json.dumps(record, sort_keys=True) + "\n")
                handle.flush()
                print(json.dumps(record, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
