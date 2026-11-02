# Mathematical methods

All bases consist of row vectors with Euclidean norm. Write `Lambda` for the
original full-rank modular orthogonal lattice, `r` for its planted shortest
root vector, and `pi` for orthogonal projection onto `r` perpendicular.
The full projected lattice is `pi(Lambda)`; the jet lattice is a separate
structured sublattice described by the first-jet construction.

## Modular orthogonal lattices

The ACD shifts are `N^(max(m-|k|,0)) * g_k`; HM shifts are
`N^(m-k_0) * g_k`, indexed by total-degree monomials up to `t`.
There are `ell` variables and total degree at most `t` (HM uses `t=m=3`).
An integer coefficient vector is in the modular orthogonal lattice exactly
when its inner product with every shift coefficient vector vanishes modulo
`N^m`. Coordinate scales balance the monomials. The implementation constructs
the same congruence kernel using exact integer transformations and checks
membership, index, and volume. See `core/lattices.py`,
`baseline/second_minimum/instances.py`, and `baseline/shared/dual_moment.py`.
Public experiments use the original lattice `B` only.

## Projection and successive minima

Because `r` is a shortest primitive lattice vector, every pair of independent
lattice vectors contains a vector outside its line. Projection does not
increase norm, hence `lambda_1(pi(Lambda)) <= lambda_2(Lambda)`.
In particular, a certified lower bound on the projected first minimum is a
lower bound on the original second minimum. A *candidate upper bound* on the
projected first minimum cannot serve as this lower bound.

A lattice vector `v` independent of `r` gives the upper bound
`max(||r||, ||v||)` for `lambda_2(Lambda)`. Original-lattice membership and
positive Gram determinant `||r||^2 ||v||^2 - <r,v>^2` are checked exactly.
Lifting a projected witness and reducing its root coefficient gives another
valid original-lattice upper bound. These two measurements share their witness.

## Certified lower bounds and changed bases

A plain Gram--Schmidt lower bound depends on the basis. The published current
intervals strengthen it with exact block certificates:

1. Form a projected integer basis `P` with a saved rational physical scale `s`.
2. Propose primal or dual BKZ changes of basis using numerically scaled models.
   Retain an exact integer matrix `C` only after verifying `abs(det(C)) = 1`.
   The actual basis is `C*P`, not the floating-point model.
3. For its exact Gram matrix `G`, produce a dyadic triangular lower form
   `H = c*L*diag(D)*L^T`, with positive diagonal `D`.
4. Verify `G-H` is positive definite by exact diagonal dominance after an
   exact invertible triangular congruence. Floating-point Cholesky proposes
   the form; integer arithmetic validates the inequality.
5. Partition the triangular form into contiguous blocks. For a nonzero
   coefficient vector, its last nonzero block contributes at least that
   block's first minimum. Other contributions are nonnegative. Thus the
   minimum of the block bounds is a lower bound for the whole form.
6. Solve manageable block SVPs without pruning, or certify that no vector
   exists below a chosen threshold. Incomplete blocks retain their safe
   Gram--Schmidt fallback. Multiply squared bounds by `s^2`.

For cases `(2,1,11)`, `(2,1,12)`, and `(3,2,6)`, the selected basis uses
reversed dual BKZ-60; `(2,1,13)` uses primal BKZ-60. For `(3,2,7)`, none of the
tested alternatives improves the baseline's certified block lower bound.
The objective is the largest *certified tested* lower bound, not a proof of
global basis optimality. Improving a GS minimum alone need not improve this
objective.

## Exact minima and heuristic candidates

The archived exact searches use coefficient exclusions with exact residual
checks, integer Gram/Bareiss preparation, complete enumeration, and checked
frontier coverage. The fplll-based native subtree interfaces use unpruned
enumeration under validated source contracts. Sieving and BKZ propose vectors;
their norms and membership are checked exactly, but shortestness requires
complete search coverage or a matching certified lower bound.

Frozen native source hashes are preserved. A rebuilt binary may differ due
to compiler versions and installation paths. Public replay checks bind to
the frozen source and keep historical binary hashes as historical metadata;
they do not claim that a relocated binary has the same bytes.

## Predictions and reduction measurements

For rank `d`, the GH column uses
`(det(projected lattice) / volume(unit ball in dimension d))^(1/d)`.
The jet column retains the paper's asymptotic first-jet prediction. The saved
records state these conventions. Small changes from asymptotic ball-volume
formulas should not be mistaken for experimental variation.

For dimension `n`, the measured RHF is
`(||b1|| / det(Lambda)^(1/n))^(1/n)`. It can be below one for planted
lattices. The historical column `1.02^(2n) * lambda_1` remains a comparison
threshold. Substituting a measured RHF is descriptive; it is not by itself a
new rigorous approximation guarantee. The LLL Lovasz parameter (0.8 in the
original runs) and the root-Hermite factor are different quantities.
