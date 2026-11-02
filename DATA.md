# Data and interpretation

`acd/results/raw/acd_paper.jsonl` contains the 13 original ACD dual/flatter
timing records. `second_minimum/data/paper_table_rows.csv` preserves the 21
original paper rows and their original precision and labels. Their hashes and
record counts are in `DATA_MANIFEST.json`; run `python scripts/audit_data.py`.
Original per-stage primal ACD, HM sweep, and long minimum logs were not
recovered. The surviving paper values are recorded/derived data, not raw logs.

`results/minima.json` consolidates the later 21 case records without changing
their exact rational squared norms, prediction conventions, statuses, or
source references. The CSV comparisons are generated from this file.
`results/rhf_base_21.csv` and `results/rhf_table1_acd_13.csv` report actual saved-basis RHFs.
The conventional RHF is `(norm(b1)/det(L)^(1/n))^(1/n)`; the alternative
`n-1` convention and the historical `1.02^(2n)*lambda1` column are kept.

The five unfinished full projected lengths are candidate upper bounds on
their first minimum. They cannot be used as certified lower bounds on the
original second minimum. A projected candidate lifted to an independent
original-lattice vector gives a valid upper bound, but is not independent
evidence of prediction accuracy. The tables keep this bound separate from
original-lattice BKZ witnesses and the revised paper's upper reports.

`results/revised_lambda2_upper_bounds.csv` transcribes printed Table 5,
page 8 (PDF page 9), with its document hash and three-decimal precision.
The PDF is not included. Rounded printed reports and locally verified exact
squared norms have separate sources. Archived exact labels are retained;
the compact reorganization does not rerun full-dimensional searches.

## Companion Release data

The 107 unchanged ZIP assets contain the bases, exact transforms, witnesses,
lower forms, and archived proof inputs. They include read-only backups of
34 completed search-frontier databases. About 3 GiB is needed after unpacking.
Download all ZIPs into a directory and set:

```bash
export DUAL_COPPERSMITH_DATA=/path/to/downloaded-assets
python scripts/check.py --data-dir "$DUAL_COPPERSMITH_DATA"
```

`second_minimum/data/release_assets.json` binds each asset by name, byte count,
SHA-256, member count, and case/group. Whole-archive hashes bind every member,
including names and contents; extraction also checks CRCs and safe paths.
The full per-file repackaging inventory is omitted from this compact source
tree. The scientific data ZIPs and their internal hash references are unchanged.
Source references in result records resolve relative to the unpacked work
directory. Public proof inputs use relative paths; historical native binary
hashes identify original executions, not relocated rebuilds.

Large deterministic Bareiss files can be rebuilt from the saved exact Gram
with `scripts/rebuild_preparation.py`; the archived preparation hash must match.
Search can be exponential. A timeout provides no shortestness proof, and the
best tested basis is not claimed globally optimal.
