#!/usr/bin/env python3
"""Visit every transferred case sequentially, preserving resumable frontiers."""

import argparse
import json
import math
from pathlib import Path
import shutil
import signal
import subprocess
import time


def stop_group(process):
    """Allow checkpoint cleanup, then stop the entire isolated case group."""
    try:
        import os
        os.killpg(process.pid, signal.SIGINT)
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()
    except ProcessLookupError:
        process.wait()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=Path(__file__).resolve().parent /
                        'results' / 'portable-lll-20261006')
    parser.add_argument('--case', help='optional exact case directory name')
    parser.add_argument('--workers', type=int, default=0, help='0: CPU/RAM-aware selection')
    parser.add_argument('--memory-gb', type=float, default=0, help='0: automatic search budget')
    parser.add_argument('--max-seconds', type=float, default=600,
                        help='search budget per case/visit; 0: unlimited (not a whole-case deadline)')
    args = parser.parse_args(argv)
    if (args.workers < 0 or any(not math.isfinite(value) or value < 0
                              for value in (args.memory_gb, args.max_seconds))):
        parser.error('worker, memory and time limits must be nonnegative and finite')
    sage = shutil.which('sage')
    if not sage:
        parser.error('SageMath executable "sage" is required')
    data = args.data.resolve()
    rows = json.loads((data / 'preparations.json').read_text())
    if args.case:
        rows = [row for row in rows if row['case'] == args.case]
        if not rows:
            parser.error('unknown case: ' + args.case)

    def priority(row):
        path = data / row['case'] / 'result.json'
        certified = path.exists() and json.loads(path.read_text()).get('exact') is True
        return (not certified, row['enumeration_dimension'], row['case'])

    # Validate/show already-completed cases before potentially expensive searches.
    rows.sort(key=priority)
    outcomes = []
    active = None
    try:
        for row in rows:
            case = row['case']
            print('START ' + case, flush=True)
            command = [sage, '-python', str(Path(__file__).with_name('enumerate_prepared.py')),
                       str(data / case), '--workers', str(args.workers),
                       '--memory-gb', str(args.memory_gb), '--max-seconds', str(args.max_seconds)]
            began = time.monotonic()
            active = subprocess.Popen(command, start_new_session=True)
            returncode = active.wait()
            active = None
            outcome = dict(case=case, returncode=returncode,
                           invocation_wall_seconds=time.monotonic() - began)
            outcomes.append(outcome)
            temporary = data / 'enumeration_sweep.json.tmp'
            temporary.write_text(json.dumps(outcomes, indent=2) + '\n')
            temporary.replace(data / 'enumeration_sweep.json')
            print('FINISH ' + json.dumps(outcome, sort_keys=True), flush=True)
            if returncode not in (0, 2):
                return returncode if returncode > 0 else 1
    except BaseException:
        if active is not None:
            stop_group(active)
        raise
    return 2 if any(row['returncode'] == 2 for row in outcomes) else 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print('Interrupted. Repeat the same command to resume committed work.', flush=True)
        raise SystemExit(130)
