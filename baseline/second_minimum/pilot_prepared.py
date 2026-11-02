#!/usr/bin/env sage-python
"""Short fresh-frontier timing pilots, without changing prepared source cases."""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import sys
import threading
import time

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(name, '1')
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.exact_minima import (_atomic_json, _build_native, _decode_preparation,
                                _lock, _read_document, _search, exact_successive_minima,
                                load_minima_input)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prepared_directory', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=30)
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--memory-gb', type=float, default=6)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args(argv)
    if args.seconds <= 0 or args.jobs < 1 or args.workers < 1 or args.memory_gb <= 0:
        parser.error('positive bounded pilot settings are required')
    source, output = args.prepared_directory.resolve(), args.out.resolve()
    output.mkdir(parents=True, exist_ok=args.resume)
    # Sage initializes signal handling and must first be imported on this thread.
    import sage.all  # noqa: F401
    binary = _build_native(output)
    rows = json.loads((source / 'preparations.json').read_text())
    aborted = threading.Event()
    started = time.monotonic()

    def run(row):
        if aborted.is_set():
            return None
        case = row['case']
        folder = output / case
        first = folder / 'pass1'
        if args.resume and (folder / 'pilot.json').exists():
            saved_summary = json.loads((folder / 'pilot.json').read_text())
            if saved_summary['pilot_total_seconds'] is not None and (
                    saved_summary['pilot_total_seconds'] < saved_summary['pilot_search_seconds']):
                saved_summary['pilot_total_seconds'] = None
                saved_summary['wall_timing_note'] = 'Not recorded by interrupted collector; search timing is intact.'
                _atomic_json(folder / 'pilot.json', saved_summary)
            return saved_summary
        first.mkdir(parents=True, exist_ok=args.resume)
        # These source files are immutable. Only pilot-owned search/result files
        # are written; ZIP packaging later follows links into independent bytes.
        for name in ('input.json', 'instance.json', 'bundle.json'):
            target = folder / name
            expected = source / case / name
            if target.exists():
                if target.resolve() != expected.resolve():
                    raise ValueError('Unexpected pilot source: ' + str(target))
            else:
                target.symlink_to(expected)
        for name in ('input_basis.json', 'ambient_lll.json', 'quotient_lll.json',
                     'prepared.json', 'integer.prepared', 'integer.input'):
            target = first / name
            expected = source / case / 'pass1' / name
            if target.exists():
                if target.resolve() != expected.resolve():
                    raise ValueError('Unexpected pilot source: ' + str(target))
            else:
                target.symlink_to(expected)
        doc = _read_document(first / 'prepared.json')
        data = _decode_preparation(doc)
        began = time.monotonic()
        saved = first / 'search_result.json'
        search = json.loads(saved.read_text()) if args.resume and saved.exists() else None
        recovered = search is not None
        if search is None:
            with _lock(folder):
                search = _search(data, first, binary, workers=args.workers, memory_gb=args.memory_gb,
                                 max_seconds=args.seconds, max_nodes=0, chunk_nodes=20000,
                                 chunk_seconds=1, progress_seconds=10)
        if search['prepared_sha256'] != doc.get('frontier_identity', doc['payload_sha256']):
            raise ValueError('Pilot preparation identity changed')
        certified = search['complete'] and int(search['off_line_sq']) >= data['root_sq']
        summary = {key: row[key] for key in ('case', 'application', 't', 'm', 'ell',
                                            'lattice_dimension', 'enumeration_dimension')}
        summary.update(pilot_search_seconds=search['invocation_search_seconds'],
                       pilot_total_seconds=None if recovered else time.monotonic() - began,
                       complete=search['complete'], lambda2_certified=certified,
                       needs_second_pass=search['complete'] and not certified,
                       nodes=search['nodes'], leaves=search['leaves'],
                       frontier_jobs=search['frontier_jobs'], reason=search['reason'],
                       worker_slot_cap=search['resource_plan']['workers'],
                       prepared_identity=search['prepared_sha256'])
        if recovered:
            summary['wall_timing_note'] = 'Not recorded by interrupted collector; search timing is intact.'
        _atomic_json(folder / 'pilot.json', summary)
        print('PILOT_RESULT ' + json.dumps(summary, sort_keys=True), flush=True)
        return summary

    summaries = []
    pool = ThreadPoolExecutor(max_workers=args.jobs)
    futures = []
    try:
        futures = [pool.submit(run, row) for row in rows]
        for future in as_completed(futures):
            summary = future.result()
            if summary is not None:
                # Keep Sage matrix operations on the main thread. The native
                # search itself is safe to run concurrently in worker threads.
                folder = output / summary['case']
                if summary['lambda2_certified'] and not (folder / 'result.json').exists():
                    basis, _ = load_minima_input(folder)
                    doc = _read_document(folder / 'pass1' / 'prepared.json')
                    result = exact_successive_minima(
                        basis, folder, workers=args.workers, memory_gb=args.memory_gb,
                        quotient_delta=doc['quotient_delta'])
                    if not result['exact'] or result['passes'][0]['invocation_nodes'] != 0:
                        raise AssertionError('Completed pilot was unexpectedly repeated')
                summaries.append(summary)
                ordered = sorted(summaries, key=lambda item: next(
                    i for i, row in enumerate(rows) if row['case'] == item['case']))
                _atomic_json(output / 'pilots.json', ordered)
    except BaseException:
        aborted.set()
        for future in futures:
            future.cancel()
        raise
    finally:
        pool.shutdown(wait=True, cancel_futures=True)
    _atomic_json(output / 'manifest.json', dict(
        schema='enumeration-timing-pilots-v1', case_count=len(summaries),
        seconds_per_search=args.seconds, jobs=args.jobs, requested_workers=args.workers,
        memory_gb_per_search=args.memory_gb, last_collector_run_seconds=time.monotonic() - started,
        summed_search_seconds=sum(row['pilot_search_seconds'] for row in summaries),
        certified_cases=sum(row['lambda2_certified'] for row in summaries),
        source_changed=False, limitations='time budget is cooperative and starts after native INIT; '
        'frontier job counts are not completion percentages'))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
