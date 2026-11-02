# Dual Coppersmith experiments

Anonymous research artifact for the original ACD and Herrmann--May (HM)
experiments. The 21 seeded cases comprise 13 ACD and 8 HM cases. This repository
contains the original recovery code, reduction measurements, projected-lattice
experiments, exact witnesses, and certified bounds for the second successive
minimum of the original lattice.

Start with [the 21-case comparison](results/tables/all_21_comparison.csv) and
[the five unfinished cases](results/tables/unfinished_5_comparison.md).
All reported lengths in those tables are base-2 logarithms.

| Quantity | Archived status |
| --- | --- |
| Original-lattice first minimum | 21 exact |
| Original-lattice second minimum | 16 exact; 5 certified intervals |
| Full projected-lattice first minimum | 16 exact; 5 checked candidate upper bounds |
| Jet projected-lattice first minimum | 21 exact |
| GH and jet formulas | Heuristic predictions |

The five unfinished ACD cases are `(2,1,11)`, `(2,1,12)`, `(2,1,13)`,
`(3,2,6)`, and `(3,2,7)`. Exact labels describe archived audited searches.
Packaging checks do not rerun every full-dimensional exhaustive search.
No claim of global optimality is made for the optimized bases.

## Contents

| Directory | Purpose |
| --- | --- |
| `baseline/` | Original recovery implementations, minima engine, and paper tables |
| `core/` | Original ACD/HM lattice construction and portable reporting |
| `configs/` | The 21 public case identifiers and parameter settings |
| `experiments/engine/` | Search, sieving, reduction, basis optimization, and proof checkers |
| `results/cases/` | Exact squared norms, status, prediction conventions, and source references |
| `results/tables/` | Human-readable comparisons and actual RHF measurements |
| `results/figures/` | Exportable comparison figure |
| `data/` | Release asset checksums, member manifests, and repackaging provenance |
| `docs/` | Mathematical methods, reproduction instructions, and data interpretation |
| `scripts/`, `tests/` | Portable entry points and checks |
| `third_party/` | Pinned upstream identities and preserved licenses |

## Quick checks

The source and summary checks use only Python's standard library:

```bash
python scripts/verify_submission.py
python scripts/make_tables.py --output /tmp/dual-tables
```

For scientific checks, use Linux or WSL, SageMath 10.9, GMP, fplll, and a C++17
compiler. A version-pinned conda environment is described in `environment.yml`;
`environment.lock.json` records the versions actually used for the packaging
checks. It is a selected version record, not a complete dependency lock.

```bash
micromamba create -f environment.yml
micromamba activate dual-coppersmith
python scripts/smoke.py --work /tmp/dual-smoke
```

The data ZIPs are companion GitHub Release assets. Download all assets into
one directory and set its location; avoid committing the ZIPs to Git history.
In the prepared package, they are in the sibling directory
`dual-coppersmith-release-assets/`, which the extraction script also detects.

```bash
export DUAL_COPPERSMITH_DATA=/path/to/downloaded-assets
python scripts/verify_submission.py --data-dir "$DUAL_COPPERSMITH_DATA"
python scripts/verify_science.py --work /tmp/dual-science
```

The scientific check reconstructs exact membership, independence, changes of
basis, and lower-form inequalities. `--replay-seconds 240` additionally reruns
decisive small-block searches and fails on timeout. It does not recertify every
archived full-dimensional shortest-vector search.

See [reproduction instructions](docs/REPRODUCING.md),
[mathematical methods](docs/METHODS.md), and [data provenance](docs/DATA.md).
New runs write into the chosen work directory. Preserve the repository and
downloaded assets as the reference snapshot.

## Interpretation

`lambda_GH` and `lambda_jet` are predictions for projected geometry, not
measured shortest-vector lengths. A projected candidate lifted into the
original lattice gives a valid upper bound, but that bound and the candidate
are dependent evidence. The tables therefore keep revised-paper upper bounds,
original-lattice BKZ witnesses, and lifted projected witnesses separate.
Revised-paper numbers are rounded printed values; local exact squared norms
are stored separately.

The actual root-Hermite factor is a descriptive reduction measurement. It is
distinct from the Lovasz parameter of LLL and does not replace an SVP proof.
Upstream copyright notices and licenses are retained.
