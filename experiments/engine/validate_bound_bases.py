"""Independent exact reconstruction and replay of optimized basis certificates."""
import os
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[key]="1"
import argparse
import concurrent.futures
import csv
from pathlib import Path
import subprocess
import sys
import mpmath as mp
from sage.all import ZZ,QQ,matrix,diagonal_matrix
from optimize_bound_basis import DEST,INDICES,load
from certify_bound_bases import verify_binary
from lambda2_interval_experiment import fraction,lognorm,ceil_sqrt
from projected_minima import HERE,OUT,GRID,name,read,save
from shared import exact_minima as engine
from fast_search import build_binary
from partition_search import binaries,run_job

mp.mp.dps=90


def parse(path,H,initial,y):
    n=H.nrows()
    with path.open() as stream:
        assert stream.readline().split()==["FAST-LIFT-1",str(n),"1","1",str(initial)]
        assert stream.readline().split()==["0"]*n
        assert stream.readline().split()==list(map(str,y))
        for row in H.rows(): assert stream.readline().split()==list(map(str,row))
        assert not stream.read().strip()


def audit_trial(folder,P,scale,replay_seconds=0):
    path=folder/"bounds.json";doc=read(path);digest=engine._hash_file(path)
    transform=read(folder/"basis_transform.json");C=matrix(ZZ,transform["transform_rows"])
    assert C.nrows()==C.ncols()==P.nrows() and abs(C.det())==1
    assert engine._digest(transform["transform_rows"])==transform["transform_sha256"]==doc["transform_sha256"]
    assert engine._hash_file(HERE/transform["source_embedding"])==transform["source_embedding_sha256"]
    assert engine._hash_file(folder/"basis_transform.json")==doc["transform_file_sha256"]
    form=read(folder/"lower_form.json")
    assert engine._hash_file(folder/"lower_form.json")==doc["lower_form_sha256"]
    L=matrix(ZZ,form["L_numerators"]);M=matrix(ZZ,form["M_numerators"])
    D=list(map(ZZ,form["D_numerators"]));den=ZZ(2)**form["fraction_bits"];d=P.nrows()
    assert L.nrows()==L.ncols()==M.nrows()==M.ncols()==len(D)==d and min(D)>0
    for A in (L,M):
        assert all(A[i,i]==den for i in range(d))
        assert all(A[i,j]==0 for i in range(d) for j in range(i+1,d))
    R=C*P;G=R*R.transpose();power=ZZ(2)**form["gram_scale_power"]
    residual=M*(den**3*G-power*L*diagonal_matrix(ZZ,D)*L.transpose())*M.transpose()
    margin=min(residual[i,i]-sum(abs(residual[i,j]) for j in range(d) if i!=j) for i in range(d))
    assert margin>0 and str(margin)==form["positive_diagonal_dominance_margin"]
    factor=QQ(power)/den**3*scale**2
    best=den**2*min(D)*factor
    assert best==fraction(doc,"gs_lower_squared")
    checker,checker_hash=verify_binary();native,native_hash=build_binary();winner=None
    for partition in doc["partitions"]:
        endpoints=partition["endpoints"]
        assert endpoints[0]==0 and endpoints[-1]==d and sorted(set(endpoints))==endpoints
        values=[];inputs=[]
        for start,end in zip(endpoints,endpoints[1:]):
            bpath=folder/f"block_{start}_{end}.input";proof=read(bpath.with_suffix(".json"));n=end-start
            A=L[start:end,start:end];H0=A*diagonal_matrix(ZZ,D[start:end])*A.transpose()
            W=matrix(ZZ,proof["block_transform_rows"])
            assert W.nrows()==W.ncols()==n and abs(W.det())==1
            H=W*H0*W.transpose();initial=min(H[i,i] for i in range(n))
            i0=next(i for i in range(n) if H[i,i]==initial);unit=[int(i==i0) for i in range(n)]
            parse(bpath,H,initial,unit)
            assert engine._hash_file(bpath)==proof["input_sha256"]
            # Relocation and compiler versions can change binary bytes. Bind
            # the rebuilt programs to frozen source, retaining archived binary
            # identities as historical metadata rather than claiming equality.
            assert proof["checker_source_sha256"]==checker_hash
            assert proof["native_source_sha256"]==native_hash
            assert all(len(proof[key])==64 and int(proof[key],16)>=0
                       for key in ("checker_binary_sha256","native_binary_sha256"))
            checked=subprocess.run([str(checker),str(bpath)],text=True,capture_output=True,check=True).stdout.splitlines()
            assert checked[0].split()[:2]==["GRAM_LLL_OK",str(n)]
            bounds=list(map(ZZ,proof["coefficient_bounds"]))
            if proof["coefficient_bound_method"]=="exact_bareiss_coefficient_bound":
                assert bounds==list(map(ZZ,checked[1].split()))
            else:
                assert proof["coefficient_bound_method"]=="exact_inverse_gram_cauchy_schwarz"
                inverse=H.inverse();assert bounds==[ceil_sqrt((initial-1)*inverse[i,i]) for i in range(n)]
            assert len(bounds)==n and max(bounds)<2**29
            gsnum,gsden=map(ZZ,checked[2].split())
            fallback=max(QQ(gsnum)/gsden,QQ(den**2*min(D[start:end])))
            assert fallback==fraction(proof,"fallback_gs_squared")
            value=fallback;replay_path=bpath;expected_best=None
            threshold=proof.get("threshold_proof")
            if proof["complete_svp"]:
                result=proof["result"];z=matrix(ZZ,[result["y"]]);value=QQ(result["best"])
                assert 0<value<=initial and not z.is_zero() and (z*H*z.transpose())[0,0]==value
                assert result["precision"]>=128;expected_best=str(value)
            elif threshold:
                replay_path=folder/threshold["input"];cutoff=ZZ(threshold["cutoff"])
                assert 0<cutoff<=initial and engine._hash_file(replay_path)==threshold["input_sha256"]
                parse(replay_path,H,cutoff,[0]*n)
                result=threshold["result"]
                assert ZZ(result["best"])==cutoff and not any(map(ZZ,result["y"])) and result["precision"]>=128
                value=max(fallback,QQ(cutoff));expected_best=str(cutoff)
            else: assert proof["result"] is None
            assert value==fraction(proof,"lower_squared")
            values.append(value);inputs.append((replay_path,expected_best,n))
        bound=min(values)*factor
        assert bound==fraction(partition,"lower_squared")
        if bound>best:
            best=bound;winner=inputs[values.index(min(values))]
    assert best==fraction(doc,"lower_squared")
    replay=None
    if replay_seconds and winner and winner[1] is not None:
        binary,_,_,_=binaries();replay=run_job(binary,winner[0],[],replay_seconds)
        assert replay is not None and replay["best"]==winner[1]
    assert engine._hash_file(path)==digest
    checked=dict(passed=True,bounds_sha256=digest,unimodular_and_lower_form_exactly_verified=True,
                 block_grams_and_lll_contracts_exactly_verified=True,decisive_block_replayed=replay is not None,
                 replay_nodes=replay["nodes"] if replay else None)
    save(folder/"audit.json",checked)
    return checked


def case(index,replay_seconds):
    loaded=load(index);folder,P=loaded[:2];baseline=read(folder/"baseline.json")
    scale=QQ(baseline["scale_numerator"])/ZZ(baseline["scale_denominator"])
    result=read(folder/"result.json");digest=engine._hash_file(folder/"result.json")
    best=fraction(result,"best_lower_squared");old=fraction(result,"baseline_lower_squared")
    assert best>=old
    trial_docs=[read(HERE/cert/"bounds.json") for cert in result["tested_certificates"]]
    assert best==max([old]+[fraction(item,"lower_squared") for item in trial_docs])
    winner=result["selected_certificate"]
    audits=[]
    for certificate in result["tested_certificates"]:
        audits.append(audit_trial(HERE/certificate,P,scale,replay_seconds if certificate==winner else 0))
    basis=read(folder/"best_basis.json");C=matrix(ZZ,basis["transform_rows"])
    assert abs(C.det())==1 and C*P==matrix(ZZ,basis["rows"])
    assert QQ(basis["scale_numerator"])/ZZ(basis["scale_denominator"])==scale
    if winner: assert best==fraction(read(HERE/winner/"bounds.json"),"lower_squared")
    else: assert best==old and C.is_one()
    assert result["best_basis_sha256"]==engine._hash_file(folder/"best_basis.json")
    old_bounds=HERE/result["baseline_bounds"];old_audit=read(old_bounds.parent/"audit.json")
    assert old_audit["passed"] and engine._hash_file(old_bounds)==old_audit["bounds_sha256"]==result["baseline_sha256"]
    assert old==fraction(read(old_bounds),"lower_squared")
    assert engine._hash_file(folder/"result.json")==digest
    save(folder/"audit.json",dict(passed=True,result_sha256=digest,independently_reconstructed_trials=len(audits),
        best_basis_sha256=engine._hash_file(folder/"best_basis.json"),
        decisive_new_block_replayed=bool(winner and replay_seconds>0),
        binary_identity_policy="rebuilt from matching frozen source; archived binary digests retained"))
    print("OPTIMIZED_BASIS_AUDIT_OK",index,len(audits),result["best_lower_log2"],flush=True)


def comparison():
    with (OUT/"empirical_comparison"/"paper_lambda2_comparison.csv").open(newline="",encoding="utf-8-sig") as stream:
        paper={row["case"]:row for row in csv.DictReader(stream)}
    rows=[]
    for index in INDICES:
        folder=DEST/name(GRID[index]);result=read(folder/"result.json");audit=read(folder/"audit.json")
        assert audit["passed"] and audit["result_sha256"]==engine._hash_file(folder/"result.json")
        ref=paper[result["case"]];a=mp.mpf(result["baseline_lower_log2"]);b=mp.mpf(result["best_lower_log2"])
        original_bounds=read(HERE/result["baseline_bounds"])
        p=mp.mpf(original_bounds["prediction_log2"]);u=mp.mpf(ref["paper_lambda2_log2"])
        old_form=read((HERE/result["baseline_bounds"]).parent/"lower_form.json")
        baseline=read(folder/"baseline.json");scale=QQ(baseline["scale_numerator"])/ZZ(baseline["scale_denominator"])
        old_gs=QQ(ZZ(2)**old_form["gram_scale_power"])*min(map(ZZ,old_form["D_numerators"]))/ZZ(2)**old_form["fraction_bits"]*scale**2
        selected_gs=old_gs;best_gs=old_gs
        for certificate in result["tested_certificates"]:
            trial=read(HERE/certificate/"bounds.json");value=fraction(trial,"gs_lower_squared")
            best_gs=max(best_gs,value)
            if certificate==result["selected_certificate"]: selected_gs=value
        rows.append(dict(index=index,case=result["case"],projection_dimension=read(folder/"baseline.json")["dimension"],
            baseline_lower_log2=str(a),optimized_lower_log2=str(b),gain_percent=str(100*mp.expm1((b-a)*mp.log(2))),
            baseline_gs_lower_log2=str(lognorm(old_gs)),selected_basis_gs_lower_log2=str(lognorm(selected_gs)),
            best_tested_gs_lower_log2=str(lognorm(best_gs)),
            prediction_log2=str(p),optimized_lower_over_prediction=str(mp.power(2,b-p)),
            paper_lambda2_upper_log2=str(u),paper_upper_over_optimized_lower=str(mp.power(2,u-b)),
            proposed_bases=result["proposed_bases"],certified_trials=result["certified_trials"],
            selected_basis=result["selected_label"],audit_passed=True,global_basis_optimality_proven=False,
            basis_file=str((folder/"best_basis.json").relative_to(HERE))))
    save(DEST/"comparison.json",rows)
    with (DEST/"comparison.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print("OPTIMIZED_COMPARISON_SAVED",str(DEST/"comparison.csv"),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--index",type=int)
    parser.add_argument("--replay-seconds",type=float,default=240);parser.add_argument("--comparison-only",action="store_true")
    args=parser.parse_args()
    if args.comparison_only: comparison()
    elif args.index is not None: case(args.index,args.replay_seconds)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
            futures=[pool.submit(subprocess.run,[sys.executable,str(Path(__file__).resolve()),"--index",str(index),
                     "--replay-seconds",str(args.replay_seconds)],check=True) for index in INDICES]
            for future in concurrent.futures.as_completed(futures): future.result()
        comparison()
