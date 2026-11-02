# Exact minima: current contract

Scope: this standalone ACD/HM artifact, not sibling experimental engines.

## Current status

Primal and dual reduction use `fpLLL:proved` at delta 0.8. Fresh minima use
the checkpointed C++/GMP engine in `shared/exact_minima.py` and
`shared/exact_enum.cpp`; the quotient reduction defaults to delta 0.99.
Exact affine monomial transformations replace expensive SNF only in the new
minimum-run input builder. Attack constructions, HG/RR filters, GH and jet
formulas, and archived datasets are unchanged.

The engine computes actual successive minima, not merely projected norms.
It evaluates integer lifts and does not assume the shortest reduced row is
globally shortest. A shorter reference requires a second exhaustive pass.
An unfinished frontier must return `exact=false` and `lambda2_sq=null`.
Squared lengths and witnesses are exact integers; log2 fields are summaries.

Parallel work must partition, not prune, the remaining frontier. Every committed
chunk atomically replaces consumed work with disjoint unfinished work and an
exactly verified incumbent. Interrupted or failed chunks are replayed. Engine
and basis identities are checked before resume. Resource limits are cooperative,
not a worst-case runtime or hard RSS guarantee.

Unbounded preparation now saves actual public inputs, both fpLLL stage matrices
and transforms, lift data, and portable ASCII Bareiss data. Mathematical resume
identity excludes host binary hashes. Transfers rebuild workers locally and
verify exact transformations without repeating saved LLL or the PRNG. A new
root-refuted second pass can still require new preparation.

## Evidence and next step

Small exhaustive regressions, integer-lift and root-refuted counterexamples,
checkpoint resume/identity checks, SNF-equivalent input bases, and independent
OpenMP/frontier review passed.

The cold campaign had 13 exact completions and eight 600-second whole-case
timeouts. Read [the timing table](../history/experiments/exact-minima/2026-10-06-ten-minute-grid.md)
for raw timings and timeout stages. Unbounded preparation then completed all
21 rows (37.394 minutes; 18 verified warm starts, three cold). Read
[the preparation table](../history/experiments/exact-minima/2026-10-06-portable-preparation.md)
for dimensions, portability and audit boundaries.

Fresh current-coordinate pilots now certify 14 rows, including HM ell9; the
13 previous exact lengths agree. Their native search times are 0.012–2.857s,
excluding validation, compilation and initial loading. Seven ACD rows remain
incomplete after about 30–32 search seconds. No finite completion ETA follows
from censoring, historical timings, worker-slot caps or pending job counts.
The standalone Linux/WSL2 ZIP carries all 21 preparations, 14 certificates,
seven partial frontiers, source and usage instructions. Read
[the forecast/transfer record](../history/experiments/exact-minima/2026-10-06-enumeration-forecast-transfer.md)
for measurements and validation. Next: resume the seven uncertified rows on
the receiving machine, retaining the exact completion gate.

## Change record

- 2026-10-06: Adopted fpLLL and actual integer-lift minima. Read
  [the proof/implementation record](../history/design/exact-minima/2026-10-06-integer-frontier.md)
  for validation and the public-N distinction.
- 2026-10-06: Added the ten-minute timing supervisor. Read
  [the whole-case deadline contract](../history/design/exact-minima/2026-10-06-whole-case-deadline.md)
  for time accounting, dimensional evidence, and timeout cleanup semantics.
- 2026-10-06: Added unbounded warm/cold preparation and transferable checkpoints.
  Read [the portable preparation contract](../history/design/exact-minima/2026-10-06-portable-preparation.md).
- 2026-10-06: Separated enumeration evidence from transfer/runtime overhead;
  certified HM ell9 and packaged source/data with seven resumable frontiers.
