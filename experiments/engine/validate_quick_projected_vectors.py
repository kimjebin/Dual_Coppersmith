"""Check quick-SVP vectors against the original public modular lattice.

This validates membership and exact length, not absence of shorter vectors.
"""
import hashlib
import json

from projected_minima import HERE, OUT, GRID, build, variant, name, read, save, matrix, vector, ZZ, QQ, norm2
from shared import exact_minima as engine


def main():
    audits=[]
    for index in (0,1,2,3,8,9):
        case=name(GRID[index]);folder=OUT/'sieving'/case
        candidate_path=folder/'best_candidate.json'
        raw=candidate_path.read_bytes();candidate=json.loads(raw)
        candidate_hash=hashlib.sha256(raw).hexdigest()
        embedding=read(OUT/case/'full_strong'/'embedding.json')
        source=HERE/embedding['source_preparation']
        prepared=engine._read_document(source/'prepared.json')
        data=build(GRID[index]);_,exponents,_,root=variant(data,'B')
        source_root=vector(ZZ,prepared['root']);source_tail=matrix(ZZ,prepared['tail'])
        ratio=QQ(norm2(root))/norm2(source_root)
        assert ratio.denominator()==1 and ratio.numerator().is_square()
        scale=ratio.numerator().isqrt()
        assert scale*source_root in (root,-root)
        y=vector(ZZ,candidate['quotient_coefficients'])
        preimage=scale*(y*source_tail)
        assert all(preimage[j]%data['scales'][j]==0 for j in range(len(preimage)))
        unbalanced=vector(ZZ,[preimage[j]//data['scales'][j] for j in range(len(preimage))])
        affine=unbalanced*data['forward'].transpose()
        assert all(affine[j]%data['N']**exponents[j]==0 for j in range(len(exponents)))
        projection=preimage.change_ring(QQ)-QQ(preimage.dot_product(root))/norm2(root)*root
        integer_witness=vector(ZZ,candidate['integer_embedding_witness'])
        physical_scale=QQ(candidate['scale_numerator'])/ZZ(candidate['scale_denominator'])
        assert integer_witness[0]==0 and projection==physical_scale*integer_witness[1:]
        squared=QQ(candidate['physical_squared_norm_numerator'])/ZZ(candidate['physical_squared_norm_denominator'])
        assert projection.dot_product(root)==0 and projection.dot_product(projection)==squared>0
        audit=dict(case=case,index=index,passed=True,kind='exact_lattice_membership_and_length_only',
            shortestness_certified=candidate['certified_minimum'],
            candidate_sha256=candidate_hash,public_affine_congruences_verified=True,
            exact_original_lattice_preimage=list(map(str,preimage)),
            exact_projection_verified=True,exact_squared_length_verified=True,
            orthogonality_verified=True,source_preparation_sha256=engine._hash_file(source/'prepared.json'),
            physical_squared_norm_numerator=str(squared.numerator()),physical_squared_norm_denominator=str(squared.denominator()))
        save(folder/('vector_audit_'+candidate_hash[:12]+'.json'),audit)
        save(folder/'vector_audit.json',audit)
        audits.append(dict(case=case,passed=True,candidate_sha256=candidate_hash,
                           shortestness_certified=candidate['certified_minimum']))
        print('ORIGINAL_PROJECTED_VECTOR_VALID',index,case,'shortestness_certified',candidate['certified_minimum'],flush=True)
    save(OUT/'sieving'/'vector_validation.json',dict(passed=len(audits),audits=audits,
         role='Each vector is exactly the projection of an integer vector satisfying the original public modular constraints'))


if __name__=='__main__':
    main()
