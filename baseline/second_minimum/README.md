# Paper geometry and timing scripts

This directory contains the scripts behind the active ACD/HM norm columns and
the Appendix reduction-quality table.

| Script | Purpose |
|---|---|
| `acd_dense_upper_bound.sage` | Reduced-vector upper bound for dense ACD `(2,1)` and `(3,2)` rows |
| `run_exact.py` | Actual exact ACD/HM lambda1/lambda2; all 21 rows supported |
| `run_exact_grid.py` | Cold 21-row timing campaign with a 600-second whole-case deadline |
| `prepare_grid.py` | Unbounded 21-row LLL/preparation campaign with portable stage checkpoints |
| `enumerate_prepared.py` | Enumerate a transferred case without regenerating inputs or repeating saved LLL |
| `portable_minima.py` | Exact-algebra bundle verification and binary-free archive export |
| `instances.py` | Public congruence kernels via exact affine monomial transformations, avoiding SNF |
| `acd_second_minimum.sage` | Compatibility wrapper for the exact ACD engine |
| `acd_jet_covolume.sage` | ACD ambient GH and first-jet length |
| `hm_second_minimum.sage` | Compatibility wrapper for the exact HM engine |
| `hm_geometry.sage` | HM `lambda_1`, ambient GH, and first-jet length |
| `primal_lr_rr.sage` | Fresh primal LLL timing; RR requires the corresponding total-runtime record |
| `b1_bounds.py` | Deterministic Appendix B1--B4 derivation from the 21-row CSV |

Run from the `baseline/` directory:

```bash
sage second_minimum/acd_dense_upper_bound.sage 3 2 5
sage -python second_minimum/run_exact.py acd --t 3 --m 3 --ell 5 \
  --workers 32 --memory-gb 64 --out second_minimum/results/exact-rerun
sage second_minimum/acd_jet_covolume.sage 3 3 5
sage -python second_minimum/run_exact.py hm --ell 2 \
  --workers 32 --memory-gb 64 --out second_minimum/results/exact-rerun
sage second_minimum/hm_geometry.sage 2
sage second_minimum/primal_lr_rr.sage
python3 second_minimum/b1_bounds.py check
```

## Exact enumeration

Use `all` instead of `acd`/`hm` for all 13 ACD and 8 HM rows. Cases are processed
sequentially, with all selected workers assigned to the current case. Omitting
`--ell` for `hm` runs its eight paper rows. Legacy `.sage` commands above also
accept the worker, memory, checkpoint, and limit flags.

For the requested ten-minute timing table, use the isolated supervisor:

```bash
python3 second_minimum/run_exact_grid.py --timeout 600 \
  --jobs 3 --workers 8 --memory-gb 6 \
  --out second_minimum/results/exact-10min-new
```

Its default timeout is **600 seconds per complete case**, including Sage
startup, instance construction, native compilation, both LLL reductions,
Bareiss preparation, and exact enumeration. It terminates the whole process
group at the deadline, allowing up to three seconds for checkpoint cleanup.
This is different from the search-only `--max-seconds` control below.

Use a new output directory: cold timing runs deliberately reject existing
directories instead of silently timing a resumed certificate. The supervisor
saves each case's log, PID/command/deadline, built dimensions, terminal timing
record, and aggregate `timings.csv`/`timings.json`. `--jobs` controls simultaneous
cases; workers and search RAM budgets are per case. Disclose these settings
when comparing times because cases compete for machine resources.

For each square paper basis, ambient fpLLL has rank `n=binomial(t+ell,ell)`;
quotient fpLLL uses `(n-1)` rows in the same `n` ambient coordinates. The
integer enumeration dimension is `n-1`, not `ell` or the jet dimension.
Stage records distinguish dimensions actually entered from intended dimensions
when the deadline occurs in construction or preparation. Timeout rows never
receive an exact lambda2. Old v1 checkpoints should be imported into a new
portable directory with the preparation command below; the archived campaign
is never rewritten.

## Portable, unbounded preparation

Finish all LLL and integer preparation without a case or preparation timeout:

```bash
python3 second_minimum/prepare_grid.py --jobs 3 --workers 8 \
  --reuse-from second_minimum/results/exact-10min-20261006 \
  --out second_minimum/results/portable-new
sage -python second_minimum/portable_minima.py verify second_minimum/results/portable-new
python3 second_minimum/portable_minima.py export second_minimum/results/portable-new \
  --archive second_minimum/results/portable-new.tar.gz
```

Omit `--reuse-from` for fresh cold preparation. Each completed legacy preparation
can warm-start a **new** ambient fpLLL reduction after exact integral unimodular
equivalence verification. This is recorded as a warm run, not the lost untouched
output or timing of the old cold reduction. Legacy native preparation and search
frontiers are imported only when final root/tail coordinates match exactly and
the native implementation is compatible. Otherwise they are rebuilt. The old
v1 format requires the matching audited binary for import; portable v2 bundles
instead bind the native source and semantic format, not a host binary hash.

Each case saves decimal-string integer matrices and the actual public instance:

- `input.json`: seed, parameters, public ACD samples or HM coefficients/constant,
  normalized public relation, monomial order, and the original integer basis.
- `pass1/ambient_lll.json`: actual fpLLL input/output, delta, warm preconditioner,
  unimodular transformation, provenance, Sage version, and stage time.
- `pass1/quotient_lll.json`: integer projected input/output, coordinate scaling,
  row order, unimodular transformation, delta, and stage time.
- `pass1/prepared.json` and `integer.prepared`: root/tail lift basis, integer
  root shifts, exact Gram content, incumbent, hashes, and portable ASCII Bareiss data.
- `bundle.json`: completion, dimensions, stage times, and file hashes.
  `preparations.json`/`.csv` summarize the complete grid.

The grid never starts enumeration: `preparation_complete=true` does **not** mean
`exact=true` or certify a shortest reference root. Completed LLL stages are
saved immediately; an interruption repeats only an unfinished stage. Ctrl-C
stops queued work and active case process groups, retaining completed stages.
Reusing the same output reads the saved basis without invoking its PRNG again.

Extract the archive on another computer, alongside this repository and Sage,
a C++17 compiler, and GMP development headers. The export excludes `.native`
and all machine binaries; the receiving machine compiles its own native worker.
Do not copy host binaries when transferring a case manually. Run, for example:

```bash
sage -python second_minimum/enumerate_prepared.py \
  second_minimum/results/portable-new/hm-t3-m3-ell9-seed0 \
  --workers 32 --memory-gb 64
```

This validates hashes and exact transformation/lift identities, loads the
saved public basis, and skips both saved pass1 LLL stages. Default enumeration
time is unlimited; `--max-seconds`/`--max-nodes` optionally bound an invocation.
If enumeration disproves the reference root's shortestness, a new pass2
preparation is needed (not a repeat of pass1); its enumeration rank is still
`n-1`. Existing matching-coordinate frontiers continue rather than restarting.
Portability tests simulate transfer, a different rebuilt binary hash, no PRNG
or LLL calls, and interruption between preparation/import commits.

The completed seed-0 21-row campaign, dimension/timing table, verified archive
and checksum are linked in
[the dated preparation record](history/experiments/exact-minima/2026-10-06-portable-preparation.md).
For the full streaming archive-content and actual-case transfer regression:

```bash
sage -python tests/check_portable_archive.py \
  second_minimum/results/portable-lll-20261006.tar.gz \
  second_minimum/results/portable-lll-20261006
```

This manual check is separate from smoke because it reads the multi-GiB archive.

The newer [standalone ZIP](results/dual-minima-transfer-20261006.zip) includes
the necessary code as well as all 21 prepared cases, 14 exact certificates and
seven partial ACD frontiers; the receiving machine does not need this checkout.
See [the transfer guide](TRANSFER_README.md) for Linux/WSL2 prerequisites and
[the measured forecast](history/experiments/exact-minima/2026-10-06-enumeration-forecast-transfer.md)
for search-only timing and the seven unknown completion ETAs. After extraction:

```bash
python3 second_minimum/check_transfer.py
python3 second_minimum/run_prepared_all.py --workers 8 --memory-gb 6 --max-seconds 600
```

Unlike `run_exact.py`, this all-case wrapper defaults to 600 search seconds
per case/visit; repeat to resume, or use `--max-seconds 0` for unlimited search.
The saved LLL stages are not repeated. Full package validation is available as
`python3 tests/check_transfer_zip.py PATH.zip` and uses isolated copied code.

- `--workers 0` (default) chooses a CPU-affinity/cgroup-aware worker count,
  limited by estimated RAM. Explicit counts are also capped to these resources.
- `--memory-gb 0` (default) uses 80% of available memory. The explicit value is
  a GiB RSS budget for the coordinator and native search workers, not a promise
  that basis generation or fpLLL preparation fits that budget.
- `--max-seconds` and `--max-nodes` limit search per case and invocation;
  preparation is excluded. Defaults of zero impose no search limit.
- `--chunk-nodes 20000` and `--chunk-seconds 2` bound checkpoint intervals.
  Smaller chunks improve load sharing and interruption latency but add IPC/SQL
  overhead. `--progress-seconds 10` controls progress messages.
- `--quotient-delta 0.99` controls the additional fpLLL quotient reduction.
  Its cost and the initial basis reduction are preparation, not search time.

Limits are cooperative: in-flight chunks may overshoot time/node thresholds.
Measured RSS can also overshoot its estimate before scheduling stops. Memory
pressure stops new chunks and drains active work; it is not an OS-enforced cap.
Set a container/job limit as well if a strict machine-level cap is required.

Repeat the same command with the same `--out` to resume. Worker counts and
resource limits may change; basis, engine, seed, parameters, and quotient delta
must match the saved identity. Changes to the native source or mathematical
format require a new directory; unrelated Python edits do not. A
checkpoint lock prevents two coordinators from owning the same case. Each
case stores `prepared.json`, the hashed Bareiss preparation, SQLite frontier,
`result.json`, and a runner `row.json`; keep the complete directory for resuming.
Ctrl-C during search saves an incomplete checkpoint. During preparation, retry
the command; only the unfinished stage is rebuilt. Exit code 0 means every
requested case completed, 2 means at least one remains incomplete.

The engine enumerates the root-orthogonal quotient but evaluates every
candidate's **nearest integer lift** in the original lattice. Integer Bareiss
data, GMP square roots, and signed integer division certify every excluded
branch; floating-point values are used only for progress and time limits.
Canonical sign classes, cached endpoints, horizontal subtree splitting,
OpenMP preparation, and persistent native workers avoid redundant work.
There is no fixed solution-list capacity and no floating pruning.

The initial reduced row is only a primitive reference vector, not an assumed
shortest vector. If the first pass finds a shorter vector, a second pass around
that certified shortest vector computes actual lambda2. This also handles ties
and non-square full-row-rank bases. Only frontier exhaustion produces
`exact=true` and `lambda2_sq`; otherwise it is null, not a projected minimum or
an upper bound substituted into the exact column. Completed records include
integer squared lengths, witnesses, source/binary hashes, and resource stats.
Worst-case runtime is still exponential; large-row feasibility is not assured.

Read [the current design contract](design/exact-minima.md) for proof details,
validation evidence, and the public-N literature distinction.

## Historical data

`data/paper_table_rows.csv` is a provenance inventory, not a claim that all
original stdout logs survived. Its dense ACD values remain upper bounds; its
historical `enumerated` values came from projected PARI enumeration and are not
relabelled exact by this implementation. See [DATA.md](../DATA.md). Fresh exact
results and fpLLL timings never overwrite the original table or raw records.
