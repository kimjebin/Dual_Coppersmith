# Dual Coppersmith Artifact

Anonymous source artifact for the ACD and Herrmann--May experiments in the
current paper.  It contains the generic total-degree modular-orthogonal
construction, the two application implementations, the primal baselines, and
the scripts and recorded values used by the active experiment tables.

This directory preserves the original recovery artifact. The current projected
lattice experiments, canonical results, and reproduction instructions are in
the repository root. Historical timing tables here retain their original
precision and status; consult `../results/` for the later exact/interval results.

## Quick start

Original requirement: SageMath 10.8; the combined artifact is checked with SageMath 10.9. Both primal and dual now use Sage's
`fpLLL:proved` backend with delta 0.8. Exact minima additionally require a
C++17 compiler and GMP development headers/libraries. Historical dual/flatter
timing records remain unchanged; fresh fpLLL timings are different measurements.

```bash
cd baseline
./scripts/run_smoke.sh
```

The smoke test checks the generic multi-relation implementation, modular-kernel
construction, data manifest, and one small HM recovery instance.  It does not
rerun the long paper grids.

## Layout

| Path | Contents |
|---|---|
| `shared/dual_moment.py` | Generic total-degree dual Coppersmith construction and ratio decoder |
| `shared/reduction.py` | Shared proved fpLLL backend for primal and dual |
| `shared/exact_minima.py`, `shared/exact_enum.cpp` | Checkpointed parallel integer enumeration for actual successive minima |
| `acd/` | ACD primal and dual implementations, table sweep, and shift-style audit |
| `hm/` | Herrmann--May primal/dual worker and isolated sweep runner |
| `second_minimum/` | ACD/HM minima, GH/jet geometry, dense upper bound, LR timing, and Appendix B1--B4 builder |
| `acd/results/raw/acd_paper.jsonl` | The 13 ACD dual timing records used by the paper |
| `second_minimum/data/paper_table_rows.csv` | The 21 active ACD/HM table rows with field-level provenance |
| `scripts/` | Smoke and data-integrity checks |

## Representative runs

Run commands from this directory.

```bash
# One ACD instance: dual recovery and primal Groebner baseline
sage -python acd/run_acd.py 2 1 10

# All 13 ACD configurations, isolated 600-second cells
python3 acd/run_sweep.py --method both --out acd/results/raw/acd_rerun.jsonl

# HM table sweep; write to a new JSONL file
python3 hm/run_hm.py --method both --out hm/results/raw/hm_rerun.jsonl

# Dense ACD rows: rigorous reduced-vector upper bound on lambda_2
sage second_minimum/acd_dense_upper_bound.sage 3 2 5

# Exact lambda_2 for any ACD/HM row, including formerly upper-bound rows
sage -python second_minimum/run_exact.py acd --t 3 --m 2 --ell 5 \
  --workers 32 --memory-gb 64 --out second_minimum/results/exact-rerun
sage -python second_minimum/run_exact.py hm --ell 2 \
  --workers 32 --memory-gb 64 --out second_minimum/results/exact-rerun

# All 21 paper rows; repeat the command to resume the same checkpoints
sage -python second_minimum/run_exact.py all --workers 32 --memory-gb 64 \
  --out second_minimum/results/exact-rerun

# Existing GH/jet columns (formulas unchanged)
sage second_minimum/acd_jet_covolume.sage 3 3 5
sage second_minimum/hm_geometry.sage 2

# Primal LR timing and Appendix B1--B4 derived table
sage second_minimum/primal_lr_rr.sage
python3 second_minimum/b1_bounds.py check
```

Exact enumeration remains exponential and some cells may take hours or longer.
Only an exhausted frontier is labelled exact; limited runs save their frontier
and report no `lambda2_sq`. See [the exact-run guide](second_minimum/README.md)
for worker/RAM controls and the certificate boundary. Always use a new output
directory for new measurements rather than overwriting committed records.

## Data status

Read `DATA.md` before using the numbers.  The ACD dual timing records are raw and
complete for the active 13-row table.  Several original primal, HM, and minimum
logs were not recovered; their surviving paper values are explicitly labelled
as recorded/derived rather than raw observations.

For dense ACD families `(2,1)` and `(3,2)`, the reported `lambda_2` value is a
reduced-basis **upper bound**, not a certified second minimum.  It must not be
used as a lower-bound certificate for the uSVP gap.

## Third-party code

The retained `shared/small_roots` utilities and ACD shift code are adapted from
[`jvdsn/crypto-attacks`](https://github.com/jvdsn/crypto-attacks).  Its MIT
license is preserved in `shared/LICENSE.crypto-attacks`. The historical
[`flatter`](https://github.com/keeganryan/flatter) reducer is no longer invoked.
