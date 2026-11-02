#!/usr/bin/env python3
"""Prepare all 21 paper rows, with NO case/LLL/preparation timeout.

Run with ordinary Python; each isolated case runs in Sage. Existing portable
inputs and completed stages are resumed, never overwritten with new randomness.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import threading
import time

ROOT = Path(__file__).resolve().parents[1]
GRID = ([('acd', 2, 1, ell) for ell in range(10, 14)] +
        [('acd', 2, 2, ell) for ell in range(10, 13)] +
        [('acd', 3, 2, ell) for ell in range(5, 8)] +
        [('acd', 3, 3, ell) for ell in range(5, 8)] +
        [('hm', 3, 3, ell) for ell in range(2, 10)])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--reuse-from', type=Path)
    parser.add_argument('--jobs', type=int, default=3)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args(argv)
    if args.jobs < 1 or args.workers < 1:
        parser.error('jobs and workers must be positive')
    sage = shutil.which('sage')
    if sage is None:
        parser.error('sage executable not found')
    args.out = args.out.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    began = time.monotonic()
    registry = {}
    registry_lock = threading.Lock()
    abort = threading.Event()

    def run(row):
        app, t, m, ell = row
        case = f'{app}-t{t}-m{m}-ell{ell}-seed{args.seed}'
        folder = args.out / case
        folder.mkdir(exist_ok=True)
        command = [sage, '-python', str(ROOT / 'second_minimum' / 'run_exact.py'),
                   app, '--t', str(t), '--m', str(m), '--ell', str(ell),
                   '--seed', str(args.seed), '--workers', str(args.workers),
                   '--prepare-only', '--out', str(args.out)]
        if args.reuse_from:
            command += ['--reuse-from', str(args.reuse_from.resolve())]
        print(f'START {case} (no timeout)', flush=True)
        started = time.monotonic()
        with (folder / 'preparation.log').open('a') as log:
            with registry_lock:
                if abort.is_set():
                    return dict(case=case, returncode=130, preparation_complete=False)
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                           start_new_session=True)
                registry[case] = process
            try:
                returncode = process.wait()
            finally:
                with registry_lock:
                    registry.pop(case, None)
        result = dict(case=case, application=app.upper(), t=t, m=m, ell=ell,
                      lattice_dimension=math.comb(t + ell, ell),
                      enumeration_dimension=math.comb(t + ell, ell) - 1,
                      returncode=returncode, invocation_seconds=time.monotonic() - started)
        if returncode == 0:
            bundle = json.loads((folder / 'bundle.json').read_text())
            result.update({key: bundle[key] for key in
                           ('preparation_complete', 'ambient_lll_seconds', 'quotient_lll_seconds',
                            'preparation_seconds', 'ambient_input_provenance', 'legacy_native_reused')})
        print('FINISH ' + json.dumps(result, sort_keys=True), flush=True)
        return result

    # Start the three known incomplete legacy preparations first, independently.
    priority = {('acd', 2, 1, 13), ('acd', 3, 2, 7), ('hm', 3, 3, 9)}
    rows = []
    pool = ThreadPoolExecutor(max_workers=args.jobs)
    try:
        pending = {pool.submit(run, row): row for row in sorted(GRID, key=lambda row: row not in priority)}
        for future in as_completed(pending):
            rows.append(future.result())
            rows.sort(key=lambda row: GRID.index((row['application'].lower(), row['t'], row['m'], row['ell'])))
            temporary = args.out / 'preparations.json.tmp'
            temporary.write_text(json.dumps(rows, indent=2) + '\n')
            temporary.replace(args.out / 'preparations.json')
    except BaseException:
        abort.set()
        with registry_lock:
            active = list(registry.values())
            for process in active:
                try:
                    os.killpg(process.pid, signal.SIGINT)
                except ProcessLookupError:
                    pass
        deadline = time.monotonic() + 5
        for process in active:
            try:
                process.wait(timeout=max(0.01, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                pass
        # Kill descendant groups too, even if the Sage coordinator exited first.
        for process in active:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        raise
    finally:
        pool.shutdown(wait=True, cancel_futures=True)
    with (args.out / 'preparations.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(key for row in rows for key in row)))
        writer.writeheader()
        writer.writerows(rows)
    manifest = dict(schema='portable-minima-grid-v1', seed=args.seed, jobs=args.jobs,
                    workers_per_case=args.workers, case_timeout_seconds=None,
                    preparation_timeout_seconds=None, elapsed_seconds=time.monotonic() - began,
                    case_count=len(rows), completed=sum(row.get('preparation_complete', False) for row in rows),
                    rows_file='preparations.json')
    (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return int(manifest['completed'] != len(GRID))


if __name__ == '__main__':
    raise SystemExit(main())
