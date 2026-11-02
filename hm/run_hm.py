#!/usr/bin/env python3
"""Run the Herrmann--May cells used by the paper and write JSONL records."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import time


PBITS = 380
QBITS = 20
X_BITS = 6
DEGREE = 3
DEFAULT_CAP = 600.0
ROOT = Path(__file__).resolve().parents[1]
WORKER = Path(__file__).with_name("hm_dual.sage")


def run_cell(ell: int, method: str, cap: float) -> dict[str, object]:
    """Run one isolated HM cell and return a normalized record."""
    command = [
        "sage", str(WORKER), str(PBITS), str(QBITS), str(ell),
        str(X_BITS), str(DEGREE), method,
    ]
    env = dict(os.environ)
    env.update(
        OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1", SAGE_NUM_THREADS="1",
    )
    started = time.perf_counter()
    try:
        proc = subprocess.run(
            command, cwd=ROOT, env=env, text=True, capture_output=True,
            timeout=cap, check=False, stdin=subprocess.DEVNULL,
        )
    except subprocess.TimeoutExpired:
        return {
            "schema_version": 1, "ell": ell, "method": method,
            "status": "timeout", "time_s": cap, "recovered": False,
        }

    result_lines = [
        line for line in proc.stdout.splitlines() if line.startswith("RESULT ")
    ]
    if not result_lines:
        return {
            "schema_version": 1, "ell": ell, "method": method,
            "status": "error", "time_s": time.perf_counter() - started,
            "recovered": False, "stderr_tail": proc.stderr[-500:],
        }

    payload = json.loads(result_lines[-1][7:])
    recovered = bool(payload.get("ok")) and proc.returncode == 0
    return {
        "schema_version": 1, "pbits": PBITS, "qbits": QBITS,
        "x_bits": X_BITS, "t": DEGREE, "m": DEGREE, "ell": ell,
        "dim": payload.get("dim"), "method": method,
        "backend": "fpLLL:proved" if method == "dual" else "fpLLL:proved + Groebner",
        "seed": 0,
        "time_s": payload.get("time_s", time.perf_counter() - started),
        "recovered": recovered, "status": "ok" if recovered else "failed",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--method", choices=("dual", "primal", "both"), default="both"
    )
    parser.add_argument("--cap", type=float, default=DEFAULT_CAP)
    parser.add_argument(
        "--out", type=Path,
        default=ROOT / "hm" / "results" / "raw" / "hm_rerun.jsonl",
    )
    parser.add_argument("--force", action="store_true", help="overwrite --out")
    args = parser.parse_args()

    if args.out.exists() and not args.force:
        parser.error(f"{args.out} exists; choose another path or pass --force")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    methods = ("dual", "primal") if args.method == "both" else (args.method,)

    with args.out.open("w", encoding="utf-8") as handle:
        for ell in range(2, 10):
            for method in methods:
                record = run_cell(ell, method, args.cap)
                line = json.dumps(record, sort_keys=True)
                handle.write(line + "\n")
                handle.flush()
                print(line, flush=True)


if __name__ == "__main__":
    main()
