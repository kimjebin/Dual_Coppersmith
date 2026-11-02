"""Public lattice preconditioning only; BKZ never supplies the certificate."""
from sage.all import ZZ, matrix


def precondition(basis, block_size=40, loops=4, use_strategies=False):
    from fpylll import IntegerMatrix, BKZ
    A = IntegerMatrix.from_matrix([[int(z) for z in row] for row in basis.rows()])
    U = IntegerMatrix.identity(basis.nrows())
    options=dict(block_size=min(block_size, basis.nrows()),max_loops=loops, flags=BKZ.AUTO_ABORT)
    if use_strategies:
        options['strategies']=BKZ.DEFAULT_STRATEGY
    parameters = BKZ.Param(**options)
    BKZ.reduction(A, parameters, U=U, float_type='mpfr', precision=512)
    reduced = matrix(ZZ, [[A[i,j] for j in range(A.ncols)] for i in range(A.nrows)])
    transform = matrix(ZZ, [[U[i,j] for j in range(U.ncols)] for i in range(U.nrows)])
    assert transform*basis == reduced and abs(transform.det()) == 1
    return reduced, transform
