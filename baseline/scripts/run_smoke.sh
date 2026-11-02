#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export SAGE_NUM_THREADS=1

command -v sage >/dev/null
sage -python tests/test_core.py
sage -python tests/test_exact_minima.py
python3 tests/test_exact_grid.py
python3 tests/test_prepare_grid.py
python3 scripts/audit_data.py
python3 second_minimum/b1_bounds.py check
sage hm/hm_dual.sage 80 20 2 4 2 dual

printf 'Smoke checks completed.\n'
