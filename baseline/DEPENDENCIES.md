# Dependencies

The repository has no dependency on sibling directories.

- SageMath 10.8 (Python 3.11 or newer), including `fpylll` and PARI/GP.
- For exact minima: a C++17 compiler (GCC recommended), GMP development
  headers, and `libgmp`/`libgmpxx`. The worker compiles into its checkpoint
  directory on first use; no machine-specific binary is bundled.
- OpenMP is used for parallel Bareiss preparation when the compiler supports
  `-fopenmp`. Otherwise preparation is serial; enumeration still uses multiple
  persistent native worker processes.

All lattice reduction uses Sage's `fpLLL:proved`, normally at delta 0.8; the
exact engine's quotient reduction defaults to 0.99. Neither `flatter` nor
`cuso` is required or invoked. Legacy `use_flatter`/`--no-flatter` arguments
remain accepted for compatibility but cannot change the backend.

Historical timing records used dual/flatter, tested at commit
`d2b8026f29b4a69e987b15d4b240f8a5053275d3`. They are archival data, not claims
about the runtime of this fpLLL implementation.
