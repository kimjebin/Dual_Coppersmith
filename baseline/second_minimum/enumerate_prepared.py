#!/usr/bin/env sage-python
"""Resume exact lambda2 enumeration from a copied portable case directory."""

import argparse
import json
import os
from pathlib import Path
import sys

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(name, '1')
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shared.exact_minima import exact_successive_minima, load_minima_input, verify_minima_bundle


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case_directory', type=Path)
    parser.add_argument('--workers', type=int, default=0)
    parser.add_argument('--memory-gb', type=float, default=0)
    parser.add_argument('--max-seconds', type=float, default=0,
                        help='enumeration only; default 0 means no time limit')
    parser.add_argument('--max-nodes', type=int, default=0)
    args = parser.parse_args(argv)
    if not (args.case_directory / 'bundle.json').exists():
        parser.error('a completed portable bundle.json is required')
    verify_minima_bundle(args.case_directory)
    basis, _ = load_minima_input(args.case_directory)
    metadata = json.loads((args.case_directory / 'pass1' / 'prepared.json').read_text())
    result = exact_successive_minima(basis, args.case_directory, workers=args.workers,
                                    memory_gb=args.memory_gb, max_seconds=args.max_seconds,
                                    max_nodes=args.max_nodes, quotient_delta=metadata['quotient_delta'])
    print(json.dumps({key: result.get(key) for key in
                     ('exact', 'status', 'lambda1_sq', 'lambda2_sq', 'log2_lambda1', 'log2_lambda2')}))
    return 0 if result['exact'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
