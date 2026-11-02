"""Exact coefficient exclusion from a dyadic approximate inverse.

Floating point proposes Q only. Every accepted exclusion is an exact integer
inequality obtained by recomputing I-AQ, with ||E||_2 <= max(||E||_1,||E||_inf).
"""
import math
from sage.all import ZZ, QQ, matrix


def ceil_sqrt(n):
    r = math.isqrt(int(n))
    return ZZ(r + (r*r < n))


def bounds(basis, upper_sq, proof):
    n = basis.nrows()
    assert basis.ncols() == n
    powers = proof['row_powers']
    P = int(proof['fraction_bits'])
    Q = matrix(ZZ, proof['inverse_numerators'])
    assert Q.nrows() == Q.ncols() == n
    product = basis * Q
    top = max(powers)
    denominator = ZZ(2) ** (top+P)
    E = matrix(ZZ, n, n, lambda i, j:
        ((ZZ(2)**(powers[i]+P) if i == j else 0)-product[i, j]) * ZZ(2)**(top-powers[i]))
    epsilon_num = max(max(sum(abs(c) for c in E.row(i)) for i in range(n)),
                      max(sum(abs(c) for c in E.column(j)) for j in range(n)))
    assert epsilon_num < denominator, 'Approximate inverse residual is too large'
    frobenius = ceil_sqrt(sum(c*c for c in Q.list()))
    active, checks = [], []
    for j in range(n):
        col = ceil_sqrt(sum(c*c for c in Q.column(j)))
        # ||column_j(A^-1)|| <= col/2^P +
        #        (frobenius/2^P)*epsilon/(1-epsilon).
        numerator = col*(denominator-epsilon_num) + frobenius*epsilon_num
        bound_denominator = ZZ(2)**(P+powers[j])*(denominator-epsilon_num)
        excluded = ZZ(upper_sq)*numerator**2 < bound_denominator**2
        if not excluded:
            active.append(j)
        checks.append(bool(excluded))
    assert len(active) >= 1
    return active, dict(residual_numerator=str(epsilon_num), residual_denominator=str(denominator),
                        excluded=checks)


def restrict(basis, upper_sq):
    import numpy as np
    powers = [(ZZ(row.dot_product(row)).nbits()+1)//2 for row in basis.rows()]
    normalized = np.array([[float(QQ(z)/ZZ(2)**s) for z in row]
                           for row, s in zip(basis.rows(), powers)])
    inverse = np.linalg.inv(normalized)
    assert np.isfinite(inverse).all()
    P = 80
    Q = [[str(int(float(z)*2**P)) for z in row] for row in inverse]
    proof = dict(schema='exact-dyadic-coefficient-exclusion-v1', row_powers=powers,
                 fraction_bits=P, inverse_numerators=Q, upper_sq=str(upper_sq))
    active, checked = bounds(basis, upper_sq, proof)
    proof.update(checked, active_rows=active)
    return basis.matrix_from_rows(active), proof
