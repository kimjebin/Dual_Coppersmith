"""Use the same explicit fpLLL backend for every primal and dual reduction."""

BACKEND = "fpLLL:proved"


def selected_backend(use_flatter=False):
    """Return the actual backend; the legacy argument is compatibility-only."""
    return BACKEND


def reduce_integer_basis(basis, use_flatter=False, lll_delta=0.8,
                         return_backend=False):
    """Reduce an integer row basis with Sage's proved fpLLL implementation.

    ``use_flatter`` remains accepted for old callers, but never selects an
    external executable. Both attacks use the same algorithm and LLL delta.
    """
    reduced = basis.LLL(delta=lll_delta, algorithm=BACKEND)
    return (reduced, BACKEND) if return_backend else reduced
