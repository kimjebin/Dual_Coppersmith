"""Original-lattice basis and preparation verification helpers."""
from experiment import *
from shared import exact_minima as engine
from coefficient_certificate import bounds
from fast_search import verify_lll_precondition

def validate_stage(folder,expected):
    doc=engine._read_document(folder/'prepared.json')
    ambient=engine._read_document(folder/'ambient_lll.json')
    quotient=engine._read_document(folder/'quotient_lll.json')
    assert engine._read_document(folder/'input_basis.json')['rows']==rows(expected)
    for key in ('input_basis','ambient_lll','quotient_lll'):
        assert doc[key+'_file_sha256']==engine._hash_file(folder/(key+'.json'))
    assert doc['native_source_sha256']==engine._hash_file(engine.SOURCE)
    assert doc['preparation_sha256']==engine._hash_file(folder/'integer.prepared')
    pre=matrix(ZZ,ambient['preconditioner_rows'])
    U=matrix(ZZ,ambient['transform_rows'])
    red=matrix(ZZ,ambient['reduced_rows'])
    warm=matrix(ZZ,ambient['lll_input_rows'])
    assert pre*expected==warm and U*warm==red
    assert abs(pre.det())==abs(U.det())==1
    order=sorted(range(red.nrows()),key=lambda i:norm2(red.row(i)))
    root=red.row(order[0]); a=norm2(root)
    tail=matrix(ZZ,[red.row(i) for i in order[1:]])
    projection=a*tail-matrix(ZZ,[z*root for z in tail*root])
    coordinate=gcd(projection.list())
    projection=projection.apply_map(lambda z:z//coordinate)
    V=matrix(ZZ,quotient['transform_rows'])
    assert abs(V.det())==1
    assert quotient['ambient_row_order']==order
    assert str(coordinate)==quotient['coordinate_content']
    assert quotient['input_rows']==rows(projection)
    assert V*projection==matrix(ZZ,quotient['reduced_rows'])
    tail=V*tail
    shifts=[engine._nearest(int(row.dot_product(root)),int(a)) for row in tail.rows()]
    tail=matrix(ZZ,[row-k*root for row,k in zip(tail.rows(),shifts)])
    inner=list(tail*root)
    gram=a*(tail*tail.transpose())-matrix(ZZ,len(inner),len(inner),lambda i,j:inner[i]*inner[j])
    assert list(map(str,root))==doc['root'] and str(a)==doc['root_sq']
    assert rows(tail)==doc['tail'] and list(map(str,inner))==doc['inner']
    assert list(map(str,shifts))==doc['root_shifts']
    assert str(gcd(gram.list()))==doc['content']
    data=engine._decode_preparation(doc)
    assert engine._witness(data,doc['initial_y'])[0]==int(doc['initial_best'])
    from gram_validation import validate
    validate(data,folder)
    return data


def measurement(data,path,label):
    doc=json.loads(path.read_text())
    assert doc['case']==name(data['cfg']) and doc['variant']==label
    R=matrix(ZZ,doc['reduced_rows']); U=matrix(ZZ,doc['transform_rows'])
    exponents=variant(data,label)[1]
    n=len(data['keys'])
    assert R.nrows()==R.ncols()==doc['dim']==n
    assert all(R[i,j]%data['scales'][j]==0 for i in range(n) for j in range(n))
    raw=matrix(ZZ,[[R[i,j]//data['scales'][j] for j in range(n)] for i in range(n)])
    transformed=raw*data['forward'].transpose()
    for j,e in enumerate(exponents):
        divisor=data['N']**e
        for i in range(n):
            assert transformed[i,j]%divisor==0
            assert transformed[i,j]//divisor==U[i,j]
    # Inclusion from exact affine congruences, equality from exact volume.
    # This independently verifies the lattice without a second determinant of
    # the much larger unreduced-basis transformation matrix.
    common=gcd(R.list())
    volume=data['N']**sum(exponents)*data['R']**sum(data['cfg'][1]-sum(k) for k in data['keys'])
    assert volume%common**n==0
    scaled=R.apply_map(lambda z:z//common)
    assert abs(scaled.det())==volume//common**n
    for i in (0,1):
        assert str(norm2(R.row(i)))==doc[f'b{i+1}_sq']
        assert abs(ln2norm(R.row(i))-doc[f'log2_b{i+1}'])<1e-10
    return doc,R
