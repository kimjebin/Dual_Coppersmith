# Data Inventory and Provenance

This artifact is scoped to the active ACD and Herrmann--May experiment tables.
Files from the commented-out generic root-line appendix and unrelated exploratory
runs are deliberately excluded.

## Raw records

`acd/results/raw/acd_paper.jsonl` contains exactly the 13 seed-0 dual/flatter
timing records used by the ACD table.  The field names were normalized to the
current runner schema; numerical, seed, status, and recovery fields are
unchanged from the recovered records.

The original per-stage primal ACD logs and the original HM sweep log were not
found.  They must not be represented as raw data.  Fresh reruns can be produced
with `acd/run_sweep.py` and `hm/run_hm.py`, but hardware-dependent fresh timing
must not silently replace the recorded paper values.

## Recorded table rows

`second_minimum/data/paper_table_rows.csv` is the machine-readable inventory of
the 21 active paper rows.  It records LR/RR/dual time, `lambda_1`, the reduction
threshold, `lambda_2`, ambient GH, first-jet length, and provenance.  Values for
which the original raw log is unavailable are labelled as recorded/derived.

The second-minimum and jet scripts regenerate the mathematical quantities, but
their original long-run stdout logs were not recovered.  In particular:

- ACD `(2,2)` and `(3,3)` and all HM rows were obtained by projected PARI
  enumeration and are labelled `enumerated`.
- Dense ACD `(2,1)` and `(3,2)` rows are labelled `upper_bound`.  Their value is
  the norm of a concrete independent reduced vector, so it is a rigorous upper
  bound on `lambda_2`; closeness to the true minimum is empirical.
- PARI `qfminim(..., flag=2)` uses high-precision floating-point pruning.  The
  returned candidate norms are checked in exact integer arithmetic, but the
  archived computation is not a formal no-missed-branch certificate.

The current `second_minimum/run_exact.py` replaces that computation for fresh
runs with resumable C++/GMP integer enumeration and exact nearest-integer
lifts. It labels a result exact only after the entire frontier is exhausted.
The historical CSV, its `enumerated`/`upper_bound` labels, and all raw numerical
values remain unchanged. In particular, no exact status is retroactively
claimed for the archived projected-enumeration values. Fresh fpLLL timings
must also be distinguished from the archived dual/flatter measurements.

## Integrity checks

```bash
python3 scripts/audit_data.py
python3 second_minimum/b1_bounds.py check
```

The manifest verifies hashes and record counts.  The audit also checks the
active parameter grids, dimensions, ACD timing roundings, `lambda_2` kinds, and
the `delta_0=1.02` threshold formula.
