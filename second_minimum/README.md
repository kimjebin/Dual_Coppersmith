# Second-minimum experiments

The original geometry, dense upper-bound, and LR timing scripts remain here.
`run_exact.py` computes original-lattice successive minima with resumable GMP
workers. `experiments/` contains only the dependencies needed by the later
projection, RHF, sieving, and certified-bound experiments.

Case identifiers and indices are in `data/cases.json`. Run from the repository
root with SageMath. Work directories receive generated code, binaries, data,
and reports; use a separate directory for each new measurement.

## Fresh geometry and a small exact projected search

```bash
sage -python scripts/run_experiment.py geometry --case hm-t3-m3-ell2-seed0 \
  --work /tmp/dual-geometry --fresh-geometry
sage -python scripts/run_experiment.py projected --case hm-t3-m3-ell2-seed0 \
  --work /tmp/dual-projection --new-search -- --kind jet --backend native --seconds 120
```

Without `--new-search`, an archived exact result may be reused and the program
prints `MINIMUM_REUSED`. A fresh search omits archived minima and search states
for that case. Fresh geometry uses no saved geometry.

## Archived data, RHF, and original-lattice bounds

Download the companion ZIPs and set `DUAL_COPPERSMITH_DATA` as in `DATA.md`.
The sibling directory `dual-coppersmith-release-assets/` is also detected.

```bash
sage -python scripts/run_experiment.py rhf --work /tmp/dual-rhf
sage -python scripts/verify_science.py --work /tmp/dual-bounds
sage -python scripts/run_experiment.py optimize --case acd-t2-m1-ell11-seed0 \
  --work /tmp/dual-optimize -- --quick
sage -python scripts/run_experiment.py certify-bounds --case acd-t2-m1-ell11-seed0 \
  --work /tmp/dual-optimize --no-unpack -- --help
```

Use `-- --help` for an engine's options. `--no-unpack` continues an initialized
work directory; normal extraction refuses to overwrite changed reference data.
`verify_science.py --indices 1,2,3,8,9` restricts verification to unfinished
cases. `--replay-seconds 240` additionally reruns decisive small-block SVPs;
timeout is a failed replay. Other block/full-dimensional completion claims
still use the archived execution records. This is not a fresh full SVP proof.

For a missing large preparation file:

```bash
sage -python scripts/rebuild_preparation.py /path/to/work/certificates/CASE/STAGE/pass1
```

## Lower-bound method

For the planted shortest primitive vector `r`, orthogonal projection `pi`
satisfies `lambda1(pi(Lambda)) <= lambda2(Lambda)`. A candidate projected upper
bound cannot be substituted for a certified projected lower bound.

Each proposed projected basis is `C*P` with exact integer unimodular `C`.
Numerical primal/dual BKZ models propose `C`; determinant and transformation
checks validate it. For its exact Gram matrix, a dyadic triangular lower form
`H=c*L*diag(D)*L^T` is certified by exact positive diagonal dominance after
an invertible triangular congruence, proving `G >= H`. Contiguous blocks give
safe SVP lower bounds using the last nonzero coefficient block. Completed
unpruned searches or threshold-emptiness certificates strengthen the block
GS fallback, and the saved rational physical scale is restored exactly.

Reversed dual BKZ-60 wins the tested certificates for ACD `(2,1,11)`,
`(2,1,12)`, and `(3,2,6)`; primal BKZ-60 wins `(2,1,13)`. The `(3,2,7)` baseline
remains best among tested alternatives. A higher GS minimum alone does not
guarantee a better certified block bound. No global optimum is claimed.

Frozen native source hashes and upstream notices are retained. Rebuilt
binaries can differ by compiler/path; replay binds to the saved source
identities and retains historical binary hashes as provenance.

## Optional sieving

G6K candidates supply exactly checked vectors, not shortestness certificates.
The pinned upstream version is in `shared/licenses/README.md`. Install GNU
compilers, make, git, autoconf, automake, libtool, and pkg-config, then:

```bash
bash scripts/build_g6k.sh /tmp/dual-sieve
sage -python scripts/run_experiment.py sieve --work /tmp/dual-sieve -- \
  /path/to/projected/certificate/pass1 --algorithm bgj1 --context 80 --seconds 180
```

Hydrate the desired case and rebuild its preparation first. Export the
comparison figure with `python scripts/plot_results.py --output /tmp/dual-figures`.
