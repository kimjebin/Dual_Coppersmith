"""Exact audit of lambda2 interval certificates; replay a decisive small SVP."""
import argparse
from pathlib import Path

from sage.all import diagonal_matrix
from lambda2_interval_experiment import (DEST, HERE, OUT, GRID, name, read, save, fraction, check_member,
    verified_quotient, build, matrix, vector, ZZ, QQ, norm2, ceil_sqrt)
from shared import exact_minima as engine
from fast_search import build_binary
from partition_search import binaries, run_job


def audit(index, replay_seconds):
    cfg=GRID[index];folder=DEST/name(cfg);path=folder/"bounds.json"
    doc=read(path); digest=engine._hash_file(path)
    lower,upper=fraction(doc,"lower_squared"),fraction(doc,"upper_squared")
    assert 0<lower<=upper
    data=build(cfg);root=vector(ZZ,read(OUT/name(cfg)/"geometry.json")["root"])
    root_sq=norm2(root);check_member(root,data)
    replay=None
    if doc["exact_lambda2"]:
        source=HERE/doc["source_certificate"]
        old=read(source)
        accepted=next(row for row in read(HERE/"certified_minima.json") if row["case"]==name(cfg) and row["variant"]=="B")
        assert accepted["exact"] and accepted["certificate_sha256"]==engine._hash_file(source)==doc["source_sha256"]
        assert old["exact"] and lower==upper==ZZ(old["lambda2_sq"])
        witness=vector(ZZ,old["lambda2_witness"])
        check_member(witness,data)
        assert norm2(witness)==upper and root_sq*upper-root.dot_product(witness)**2>0
    else:
        P,scale,_,_,_=verified_quotient(index,data)
        assert scale==fraction(doc,"projected_scale")
        upper_doc=read(folder/"upper_witness.json")
        assert engine._hash_file(folder/"upper_witness.json")==doc["upper_witness_sha256"]
        witness=vector(ZZ,upper_doc["vector"])
        check_member(witness,data)
        assert norm2(witness)==ZZ(upper_doc["squared_norm"])
        assert upper==max(root_sq,norm2(witness)) and root_sq*norm2(witness)-root.dot_product(witness)**2>0
        form_path=folder/"lower_form.json";form=read(form_path)
        assert engine._hash_file(form_path)==doc["lower_form_sha256"]
        L,M=matrix(ZZ,form["L_numerators"]),matrix(ZZ,form["M_numerators"])
        D=list(map(ZZ,form["D_numerators"]));den=ZZ(2)**form["fraction_bits"];d=P.nrows()
        assert L.nrows()==L.ncols()==M.nrows()==M.ncols()==len(D)==d
        assert min(D)>0
        for A in (L,M):
            assert all(A[i,i]==den for i in range(d))
            assert all(A[i,j]==0 for i in range(d) for j in range(i+1,d))
        assert all(100*abs(L[i,j])<=51*den for i in range(d) for j in range(i))
        assert all(100*(den**2*D[i]+L[i,i-1]**2*D[i-1])>=99*den**2*D[i-1] for i in range(1,d))
        G=P*P.transpose();H=L*diagonal_matrix(ZZ,D)*L.transpose()
        difference=den**3*G-ZZ(2)**form["gram_scale_power"]*H
        checked=M*difference*M.transpose()
        margins=[checked[i,i]-sum(abs(checked[i,j]) for j in range(d) if j!=i) for i in range(d)]
        assert min(margins)>0 and str(min(margins))==form["positive_diagonal_dominance_margin"]
        form_scale=QQ(ZZ(2)**form["gram_scale_power"])/den**3
        best=QQ(den**2*min(D))*form_scale*scale**2
        native,source_hash=build_binary();native_hash=engine._hash_file(native)
        winner=None
        for partition in doc["partitions"]:
            endpoints=partition["endpoints"]
            assert endpoints[0]==0 and endpoints[-1]==d and sorted(set(endpoints))==endpoints
            values=[];blocks=[]
            for start,end in zip(endpoints,endpoints[1:]):
                block=L[start:end,start:end]
                gram=block*diagonal_matrix(ZZ,D[start:end])*block.transpose()
                n=end-start;initial=min(gram[i,i] for i in range(n))
                initial_i=next(i for i in range(n) if gram[i,i]==initial)
                unit=[int(i==initial_i) for i in range(n)]
                filename=folder/f"block_{start}_{end}.input";proof=read(filename.with_suffix(".json"))
                assert proof["input_sha256"]==engine._hash_file(filename)
                assert proof["native_source_sha256"]==source_hash
                assert len(proof["native_binary_sha256"])==64 and int(proof["native_binary_sha256"],16)>=0
                assert proof["start"]==start and proof["end"]==end and proof["dimension"]==n
                with filename.open() as stream:
                    assert stream.readline().split()==["FAST-LIFT-1",str(n),"1","1",str(initial)]
                    assert stream.readline().split()==["0"]*n and stream.readline().split()==list(map(str,unit))
                    for row in gram.rows(): assert stream.readline().split()==list(map(str,row))
                    assert not stream.read().strip()
                bounds=list(map(ZZ,proof["coefficient_bounds"]))
                assert len(bounds)==n and max(bounds)<2**29
                if proof["coefficient_bound_method"]=="exact_inverse_gram_cauchy_schwarz":
                    inverse=gram.inverse()
                    assert bounds==[ceil_sqrt((initial-1)*inverse[i,i]) for i in range(n)]
                else:
                    assert proof["coefficient_bound_method"]=="exact_triangular_integer_bound"
                    expected=[ZZ(0)]*n
                    for i in reversed(range(n)):
                        weighted=sum(abs(block[j,i])*expected[j] for j in range(i+1,n))
                        expected[i]=ceil_sqrt(QQ(initial-1)/(den**2*D[start+i]))+(weighted+den-1)//den
                    assert bounds==expected
                if proof["complete_svp"]:
                    result=proof["result"];z=vector(ZZ,result["witness"])
                    value=ZZ(result["squared_minimum"])
                    assert 0<value<=initial and z*gram*z==value and result["precision"]>=128
                else:
                    assert proof["result"] is None
                    value=den**2*min(D[start:end])
                assert str(value)==proof["lower_squared_norm"]
                values.append(value);blocks.append((filename,proof))
            value=QQ(min(values))*form_scale*scale**2
            assert value==fraction(partition,"lower_squared")
            if value>=best:
                best=value
                winner=blocks[values.index(min(values))]
        assert lower==best
        # Re-run the decisive full small-block search via the separate subtree
        # interface. A timeout grants no extra proof credit.
        if replay_seconds>0 and winner and winner[1]["complete_svp"]:
            binary,_,_,_=binaries()
            replay=run_job(binary,winner[0],[],replay_seconds)
            assert replay is not None and replay["best"]==winner[1]["lower_squared_norm"]
    assert engine._hash_file(path)==digest
    record=dict(passed=True,index=index,case=name(cfg),bounds_sha256=digest,
        exact_lambda2=doc["exact_lambda2"],public_upper_witness_membership_verified=True,
        public_upper_witness_independence_verified=True,
        lower_form_psd_and_partition_argument_verified=not doc["exact_lambda2"],
        decisive_block_replayed=replay is not None,
        replay_nodes=replay["nodes"] if replay else None)
    save(folder/"audit.json",record)
    print("LAMBDA2_INTERVAL_AUDIT_OK",index,"replay_nodes",record["replay_nodes"],flush=True)
    return record


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--indices",default=",".join(map(str,range(len(GRID)))))
    parser.add_argument("--replay-seconds",type=float,default=180)
    args=parser.parse_args()
    for index in map(int,args.indices.split(",")):
        audit(index,args.replay_seconds)
    results=[]
    for cfg in GRID:
        folder=DEST/name(cfg)
        if (folder/"audit.json").exists() and (folder/"bounds.json").exists():
            checked=read(folder/"audit.json")
            if checked.get("passed") and checked.get("bounds_sha256")==engine._hash_file(folder/"bounds.json"):
                results.append(checked)
    save(DEST/"validation.json",dict(passed=len(results),results=results))
