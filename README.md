# Dual Coppersmith Artifact

Anonymous research code and recorded results for the ACD and Herrmann--May
(HM) experiments. The artifact contains the original primal/dual recovery
implementations and the later reduction, projection, and second-minimum
experiments for the same 21 seeded cases.

## Quick start

Use Linux or WSL with SageMath 10.9, GMP, fplll 5.5.0, and a C++17 compiler.
The supplied `environment.yml` describes the Sage environment. Standard-library
result checks need Python 3.11 or newer.

```bash
python scripts/check.py
sage -python scripts/smoke.py --work /tmp/dual-smoke
```

The smoke suite checks modular kernels, exact minima on small lattices,
coefficient exclusions, corrupted Gram rejection, and one small HM recovery.
It does not repeat the long paper grids.

## Layout

| Path | Contents |
| --- | --- |
| `acd/`, `hm/` | Original primal and dual recovery code and timing runners |
| `shared/` | Modular kernels, reduction, exact-minimum engine, and upstream licenses |
| `second_minimum/` | Paper geometry and timing scripts |
| `second_minimum/experiments/` | RHF, projected minima, sieving, basis optimization, and bound verification |
| `results/` | 21-case comparisons, exact rational values, statuses, and revised-paper upper reports |
| `scripts/`, `tests/` | Experiment entry points, table generation, and scientific checks |

## Representative runs

Run from this directory and use new output directories for new measurements.

```bash
# One original ACD recovery instance and the HM sweep
sage -python acd/run_acd.py 2 1 10
python hm/run_hm.py --method both --out /tmp/hm-rerun.jsonl

# Exact original-lattice minima; repeat to resume
sage -python second_minimum/run_exact.py hm --ell 2 --workers 1 \
  --max-seconds 120 --out /tmp/hm-minima

# GH/jet formulas and original paper-table checks
sage second_minimum/acd_jet_covolume.sage 3 3 5
python scripts/audit_data.py
python second_minimum/b1_bounds.py check

# Regenerate compact comparisons without Sage or large data
python scripts/make_tables.py --output /tmp/dual-tables
```

See [second-minimum experiments](second_minimum/README.md) for the projection,
RHF, sieving, and bound commands. [DATA.md](DATA.md) explains provenance and
the companion Release data.

## Result status

There are 13 ACD and 8 HM cases. The original first minimum is archived exact
on all 21 cases. The original second minimum has 16 archived exact values and
5 certified intervals. The full projected first minimum has 16 archived exact
values and 5 checked candidate upper bounds; all 21 jet projected values are
archived exact. The unfinished cases are ACD `(2,1,11)`, `(2,1,12)`, `(2,1,13)`,
`(3,2,6)`, and `(3,2,7)`.

Start with [all 21 cases](results/all_21_comparison.csv) and
[the five unfinished cases](results/unfinished_5_comparison.md).
Table lengths are base-2 logarithms. GH and jet columns are heuristic
predictions. Revised-paper, original-lattice BKZ, and lifting-derived upper
bounds remain separate; lifting and its projected candidate share evidence.
Historical exact labels describe archived audited searches, not a fresh
full-dimensional recertification. No global basis optimality is claimed.

Historical dual/flatter timing records retain their original provenance.
Fresh runs here use Sage's proved fpLLL backend and have different timings.
The measured RHF is distinct from the LLL Lovasz parameter and is not an SVP
approximation guarantee.

## License

See [LICENSE](LICENSE), `shared/LICENSE.crypto-attacks`, and
[upstream licenses](shared/licenses/README.md). Required upstream author
attribution is preserved.
