"""HM symmetric-cube bases adapted to the public constant-coordinate flag."""
from sage.all import ZZ, QQ, matrix, vector, PolynomialRing, prod


def flag_basis(data):
    ell = data['cfg'][3]
    N, R = data['N'], data['R']
    a = list(map(ZZ, data['meta']['normalized_coefficients']))
    c = ZZ(data['meta']['normalized_constant'])
    # Homogeneous solutions to a*x = 0 mod N. No planted root is used.
    H = matrix(ZZ, ell, ell)
    H[0,0] = N
    for i in range(1,ell):
        H[i,0], H[i,i] = -a[i], 1
    H = H.LLL(delta=0.99, algorithm='fpLLL:proved', fp='xd')
    assert abs(H.det()) == N and all(sum(a[j]*row[j] for j in range(ell)) % N == 0 for row in H.rows())
    orthogonal, _ = H.gram_schmidt()
    v = vector(ZZ, [-c] + [0]*(ell-1))
    for j in reversed(range(ell)):
        q = QQ(v.dot_product(orthogonal.row(j))) / orthogonal.row(j).dot_product(orthogonal.row(j))
        nearest = (2*q.numerator()+q.denominator()) // (2*q.denominator())
        v -= nearest*H.row(j)
    assert (c+sum(a[j]*v[j] for j in range(ell))) % N == 0
    L = matrix(ZZ, [[R]+list(v)] + [[0]+list(row) for row in H.rows()])
    assert abs(L.det()) == N*R
    ring = PolynomialRing(ZZ, [f'c{i}' for i in range(ell+1)])
    cs = ring.gens()
    linear = [sum(cs[i]*L[i,j] for i in range(ell+1)) for j in range(ell+1)]
    homogeneous_keys = [(3-sum(k),)+tuple(k) for k in data['keys']]
    index = {k:i for i,k in enumerate(homogeneous_keys)}
    P = matrix(ZZ, len(index), len(index))
    for col, alpha in enumerate(homogeneous_keys):
        polynomial = prod(linear[j]**alpha[j] for j in range(ell+1))
        for beta, coefficient in polynomial.dict().items():
            P[index[tuple(beta)],col] = coefficient
    return P, [k[0] for k in homogeneous_keys]


def intersect(P, h_degree, label, data, exponents):
    if label != 'B': raise ValueError('Only the original lattice B is supported')
    powers = [0]*len(h_degree)
    candidate = matrix(ZZ, [[z*data['N']**power for z in row]
                            for row,power in zip(P.rows(),powers)])
    raw = matrix(ZZ, [[candidate[i,j]//data['scales'][j] for j in range(candidate.ncols())]
                      for i in range(candidate.nrows())])
    assert all(candidate[i,j] % data['scales'][j] == 0
               for i in range(candidate.nrows()) for j in range(candidate.ncols()))
    changed = raw*data['forward'].transpose()
    transform = matrix(ZZ, changed.nrows(), changed.ncols())
    for j,e in enumerate(exponents):
        divisor = data['N']**e
        for i in range(changed.nrows()):
            assert changed[i,j] % divisor == 0
            transform[i,j] = changed[i,j]//divisor
    assert abs(transform.det()) == 1
    return candidate, transform
