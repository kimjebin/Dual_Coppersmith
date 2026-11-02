#!/usr/bin/env python3
"""Measure the 21 paper rows with isolated whole-case wall-clock deadlines."""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import json
from math import comb, isfinite
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import threading
import time

ROOT = Path(__file__).resolve().parents[1]
WORKER = Path(__file__).with_name("run_exact.py")
ACD_GRID = [(2, 1, ell) for ell in range(10, 14)] + [
    (2, 2, ell) for ell in range(10, 13)] + [
    (3, 2, ell) for ell in range(5, 8)] + [
    (3, 3, ell) for ell in range(5, 8)]
PAPER_CASES = [("acd", *case) for case in ACD_GRID] + [
    ("hm", 3, 3, ell) for ell in range(2, 10)]


def atomic_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w") as handle:
        json.dump(value, handle, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def stop_group(proc, grace=3):
    """Stop this isolated case and its native workers, even if the leader exits."""
    try:
        os.killpg(proc.pid, signal.SIGINT)
    except ProcessLookupError:
        pass
    try:
        proc.wait(timeout=grace)
    except subprocess.TimeoutExpired:
        pass
    finally:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()


class Campaign:
    """Own all running process groups and prevent starts after cancellation."""

    def __init__(self):
        self.lock = threading.Lock()
        self.cancelled = threading.Event()
        self.active = set()

    def start(self, command, **options):
        with self.lock:
            if self.cancelled.is_set():
                raise InterruptedError("Campaign cancelled before case launch")
            proc = subprocess.Popen(command, **options)
            self.active.add(proc)
            return proc

    def release(self, proc):
        with self.lock:
            self.active.discard(proc)

    def abort(self):
        with self.lock:
            self.cancelled.set()
            processes = list(self.active)
        for proc in processes:
            stop_group(proc)


def run_case(case, *, output, timeout, workers, memory_gb, seed=0, sage="sage", campaign=None):
    app, t, m, ell = case
    name = f"{app}-t{t}-m{m}-ell{ell}-seed{seed}"
    folder = output / name
    folder.mkdir(parents=True, exist_ok=True)
    command = [sage, "-python", str(WORKER), app, "--t", str(t), "--m", str(m),
               "--ell", str(ell), "--seed", str(seed), "--workers", str(workers),
               "--memory-gb", str(memory_gb), "--max-seconds", "0",
               "--progress-seconds", "30", "--out", str(output)]
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               SAGE_NUM_THREADS="1", PYTHONUNBUFFERED="1")
    began = time.monotonic()
    deadline_reached = False
    campaign = campaign or Campaign()
    with (folder / "run.log").open("w") as log:
        proc = campaign.start(command, cwd=ROOT, env=env, stdout=log,
                              stderr=subprocess.STDOUT, start_new_session=True)
        try:
            atomic_json(folder / "supervisor.json", {"pid": proc.pid, "command": command,
                        "started_utc": datetime.now(timezone.utc).isoformat(), "timeout_seconds": timeout})
            proc.wait(timeout=max(0.001, timeout - (time.monotonic() - began)))
            # wait() may observe an exit just after its nominal deadline. Do
            # not credit a late certificate merely because no exception arose.
            if time.monotonic() - began >= timeout:
                deadline_reached = True
                stop_group(proc)
        except subprocess.TimeoutExpired:
            deadline_reached = True
            stop_group(proc)
        except BaseException:
            stop_group(proc)
            raise
        finally:
            campaign.release(proc)
    elapsed = time.monotonic() - began
    expected = comb(t + ell, ell)
    row = {"application": app.upper(), "t": t, "m": m, "ell": ell, "seed": seed,
           "expected_dimension": expected, "lattice_dimension": None,
           "quotient_dimension": None, "ambient_lll_entered": False,
           "quotient_lll_entered": False, "enumeration_entered": False,
           "last_stage": "startup_or_build", "exact": False, "log2_lambda2": None,
           "elapsed_seconds": elapsed, "timeout_seconds": timeout,
           "status": "timeout" if deadline_reached else "error",
           "returncode": proc.returncode, "requested_workers": workers,
           "requested_memory_gb": memory_gb, "checkpoint": str(folder)}
    instance_file = folder / "instance.json"
    if instance_file.exists():
        instance = json.loads(instance_file.read_text())
        row.update(lattice_dimension=instance["lattice_dimension"],
                   quotient_dimension=instance["quotient_dimension"],
                   build_seconds=instance["build_seconds"])
        if row["lattice_dimension"] != expected or row["quotient_dimension"] != expected - 1:
            raise ArithmeticError("Built dimensions differ from the paper parameter grid")
    for line in (folder / "run.log").read_text().splitlines():
        if line.startswith("EXACT_STAGE "):
            stage = json.loads(line[len("EXACT_STAGE "):])
            row["last_stage"] = stage["stage"]
            if stage["stage"] in ("ambient_lll", "quotient_lll"):
                row[stage["stage"] + "_entered"] = True
                actual = stage["rows"]
                target = expected if stage["stage"] == "ambient_lll" else expected - 1
                if actual != target:
                    raise ArithmeticError("Actual fpLLL dimensions disagree with the built basis")
            if stage["stage"] == "enumeration":
                if stage["dimension"] != expected - 1:
                    raise ArithmeticError("Actual native enumeration dimension mismatch")
                row.update(enumeration_entered=True, actual_workers=stage["workers"])
    engine_result = folder / "row.json"
    if not deadline_reached and proc.returncode == 0 and engine_result.exists():
        result = json.loads(engine_result.read_text())
        if result["exact"] and result["lambda2_sq"] is not None:
            row.update(status="complete", exact=True, lambda1_sq=result["lambda1_sq"],
                       lambda2_sq=result["lambda2_sq"], log2_lambda2=result["log2_lambda2"],
                       preparation_seconds=result["preparation_seconds"],
                       search_seconds=sum(p["invocation_search_seconds"] for p in result["passes"]),
                       nodes=sum(p["invocation_nodes"] for p in result["passes"]))
    elif not deadline_reached:
        row["error_tail"] = (folder / "run.log").read_text()[-2000:]
    atomic_json(folder / "timing.json", row)
    return row


def write_table(output, rows):
    ordered = sorted(rows, key=lambda row: PAPER_CASES.index(
        (row["application"].lower(), row["t"], row["m"], row["ell"])))
    atomic_json(output / "timings.json", ordered)
    fields = ["application", "t", "m", "ell", "lattice_dimension", "quotient_dimension",
              "ambient_lll_entered", "quotient_lll_entered", "enumeration_entered",
              "status", "elapsed_seconds", "timeout_seconds", "last_stage", "log2_lambda2"]
    temporary = output / "timings.csv.tmp"
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(ordered)
    os.replace(temporary, output / "timings.csv")
    return ordered


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=float, default=600,
                        help="whole-case wall-clock cap, including startup/build/LLL/preparation")
    parser.add_argument("--jobs", type=int, default=1, help="simultaneous isolated cases")
    parser.add_argument("--workers", type=int, default=8, help="native workers per case")
    parser.add_argument("--memory-gb", type=float, default=6, help="search RSS budget per case")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--out", type=Path, required=True, help="a new output directory")
    args = parser.parse_args(argv)
    if not isfinite(args.timeout) or args.timeout <= 0 or args.jobs < 1 or args.workers < 1:
        parser.error("timeout, jobs and workers must be positive")
    if not isfinite(args.memory_gb) or args.memory_gb <= 0:
        parser.error("memory-gb must be finite and positive")
    if shutil.which("sage") is None:
        parser.error("SageMath must be installed")
    if args.out.exists():
        parser.error("--out already exists; use a new directory for cold timing measurements")
    args.out.mkdir(parents=True)
    args.out = args.out.resolve()
    manifest = {"schema": "exact-minima-wall-timing-v1", "seed": args.seed,
                "timeout_seconds": args.timeout, "simultaneous_cases": args.jobs,
                "workers_per_case": args.workers, "search_memory_gib_per_case": args.memory_gb,
                "timing_scope": "startup + build + both LLL reductions + preparation + exact enumeration",
                "started_utc": datetime.now(timezone.utc).isoformat(),
                "cpu_affinity": sorted(os.sched_getaffinity(0)),
                "case_count": len(PAPER_CASES)}
    atomic_json(args.out / "run_manifest.json", manifest)
    rows = []
    campaign = Campaign()
    pool = ThreadPoolExecutor(max_workers=args.jobs)
    futures = {}
    try:
        futures = {pool.submit(run_case, case, output=args.out, timeout=args.timeout,
                               workers=args.workers, memory_gb=args.memory_gb,
                               seed=args.seed, campaign=campaign): case
                   for case in PAPER_CASES}
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            write_table(args.out, rows)
            print("TIMING_RESULT " + json.dumps({key: row[key] for key in
                  ("application", "t", "m", "ell", "status", "elapsed_seconds", "last_stage")}),
                  flush=True)
    except BaseException:
        campaign.abort()
        for future in futures:
            future.cancel()
        raise
    finally:
        pool.shutdown(wait=True, cancel_futures=True)
    write_table(args.out, rows)
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["completed_rows"] = sum(row["exact"] for row in rows)
    manifest["timeout_rows"] = sum(row["status"] == "timeout" for row in rows)
    manifest["error_rows"] = sum(row["status"] == "error" for row in rows)
    atomic_json(args.out / "run_manifest.json", manifest)
    print(f"TABLE {args.out / 'timings.csv'}", flush=True)
    return 1 if manifest["error_rows"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
