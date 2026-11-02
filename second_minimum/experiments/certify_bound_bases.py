"""Certify improved lambda2 lower bounds for equivalent proposed bases.

The global basis need not be LLL reduced. Each small projected-block Gram is
separately reduced by an exactly checked unimodular congruence, and its LLL
contract is rechecked using integral Bareiss arithmetic before proved SVP.
"""
import os
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[key]="1"
import argparse
import concurrent.futures
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import time
import mpmath as mp
import numpy as np
from sage.all import ZZ,QQ,matrix,diagonal_matrix,identity_matrix
from fpylll import IntegerMatrix,GSO,LLL,BKZ,FPLLL
from optimize_bound_basis import DEST,INDICES,load
from projected_minima import HERE,OUT,GRID,name,read,save,rows
from lambda2_interval_experiment import fraction, fraction_fields,lognorm,ceil_sqrt
from shared import exact_minima as engine
from fast_search import compile_once,build_binary

mp.mp.dps=80
FPLLL.set_precision(256)


def gram_lll(H):
    n=H.nrows();A=IntegerMatrix.from_matrix([[int(z) for z in row] for row in H.rows()])
    U=IntegerMatrix.identity(n)
    M=GSO.Mat(A,U=U,gram=True,float_type="mpfr",flags=GSO.INT_GRAM)
    assert M.update_gso()
    reducer=LLL.Reduction(M,delta=.999,eta=.5001);reducer()
    if n>=20:
        parameters=BKZ.Param(block_size=min(60,n),strategies=BKZ.DEFAULT_STRATEGY,
                             max_loops=4,max_time=30,flags=BKZ.MAX_LOOPS|BKZ.MAX_TIME|BKZ.AUTO_ABORT)
        BKZ.Reduction(M,reducer,parameters)()
        reducer()
    W=matrix(ZZ,[[U[i,j] for j in range(n)] for i in range(n)])
    assert abs(W.det())==1
    return W*H*W.transpose(),W


def lower_form(G,path):
    d=G.nrows();power=max(abs(z) for z in G.list()).nbits();bits=80;den=ZZ(2)**bits
    numeric=np.array([[float(QQ(z)/ZZ(2)**power) for z in row] for row in G.rows()])
    chol=np.linalg.cholesky((numeric+numeric.T)/2)
    unit=chol/np.diag(chol)[None,:];inverse=np.linalg.inv(unit)
    L=matrix(ZZ,d,d,lambda i,j:den if i==j else int(round(float(unit[i,j])*int(den))) if i>j else 0)
    M=matrix(ZZ,d,d,lambda i,j:den if i==j else int(round(float(inverse[i,j])*int(den))) if i>j else 0)
    for safety in (20,16,12,8,4):
        D=[ZZ(math.floor(float(chol[i,i]**2)*(1-2.**(-safety))*int(den))) for i in range(d)]
        assert min(D)>0
        residual=den**3*G-ZZ(2)**power*(L*diagonal_matrix(ZZ,D)*L.transpose())
        Q=M*residual*M.transpose()
        margin=min(Q[i,i]-sum(abs(Q[i,j]) for j in range(d) if i!=j) for i in range(d))
        if margin>0: break
    else: raise RuntimeError("No positive exact lower-form residual")
    doc=dict(schema="exact-arbitrary-basis-lower-form-v1",fraction_bits=bits,gram_scale_power=int(power),
        L_numerators=rows(L),M_numerators=rows(M),D_numerators=list(map(str,D)),
        safety_bits=safety,positive_diagonal_dominance_margin=str(margin),global_basis_need_not_be_lll=True)
    save(path,doc)
    return L,D,QQ(ZZ(2)**power)/den**3


def verify_binary():
    source=HERE/"verify_bound_gram.cpp";checksum=engine._hash_file(source)
    binary=HERE/"bin"/("verify_bound_gram-"+checksum[:16])
    compile_once(binary,["g++","-O3","-std=c++17",str(source),"-lgmpxx","-lgmp"])
    return binary,checksum


def block(L,D,start,end,folder,seconds,threshold_ratio=None):
    began=time.monotonic();A=L[start:end,start:end]
    original=A*diagonal_matrix(ZZ,D[start:end])*A.transpose();n=end-start
    cached_path=folder/f"block_{start}_{end}.json"
    if cached_path.exists():
        old=read(cached_path)
        if old.get("complete_svp"):
            W=matrix(ZZ,old["block_transform_rows"])
            assert abs(W.det())==1 and old["input_sha256"]==engine._hash_file(cached_path.with_suffix(".input"))
            return fraction(old,"lower_squared")
    H,W=gram_lll(original)
    initial=min(H[i,i] for i in range(n));i0=next(i for i in range(n) if H[i,i]==initial)
    y=[int(i==i0) for i in range(n)]
    path=folder/f"block_{start}_{end}.input"
    with path.open("w") as stream:
        stream.write(f"FAST-LIFT-1 {n} 1 1 {initial}\n"+" ".join(["0"]*n)+"\n"+" ".join(map(str,y))+"\n")
        for row in H.rows(): stream.write(" ".join(map(str,row))+"\n")
    checker,checker_hash=verify_binary()
    verified=subprocess.run([str(checker),str(path)],capture_output=True,text=True,check=True)
    lines=verified.stdout.splitlines();tag,dimension,largest=lines[0].split()
    assert tag=="GRAM_LLL_OK" and int(dimension)==n
    bounds=list(map(ZZ,lines[1].split()));gsnum,gsden=map(ZZ,lines[2].split())
    gs_bound=QQ(gsnum)/gsden
    bound_method="exact_bareiss_coefficient_bound"
    if max(bounds)>=2**29:
        inverse=H.inverse();bounds=[ceil_sqrt((initial-1)*inverse[i,i]) for i in range(n)]
        bound_method="exact_inverse_gram_cauchy_schwarz"
    assert max(bounds)<2**29
    native,source_hash=build_binary();result=None;complete=False;threshold_proof=None
    fallback=max(gs_bound,QQ(L[0,0]**2*min(D[start:end])))
    search_path=path;cutoff=initial
    if threshold_ratio is not None and n>45:
        log_gh_squared=(2*mp.loggamma(mp.mpf(n)/2+1)-n*mp.log(mp.pi)+mp.log(int(H.det())))/n
        cutoff=min(initial,ZZ(int(mp.floor(mp.exp(log_gh_squared)*mp.mpf(str(threshold_ratio))**2))))
        assert cutoff>0
        search_path=folder/f"block_{start}_{end}.threshold.input"
        lines=path.read_text().splitlines();lines[0]=f"FAST-LIFT-1 {n} 1 1 {cutoff}"
        lines[2]=" ".join(["0"]*n)
        search_path.write_text("\n".join(lines)+"\n")
    try:
        run=subprocess.run([str(native),str(search_path)],capture_output=True,text=True,check=True,timeout=seconds)
        tag,best,nodes,leaves,native_seconds,precision=run.stdout.splitlines()[0].split()
        assert tag=="COMPLETE"
        z=matrix(ZZ,[list(map(ZZ,run.stdout.splitlines()[1].split()))])
        assert 0<ZZ(best)<=cutoff
        value=max(fallback,QQ(best))
        search_result=dict(best=best,y=rows(z)[0],nodes=int(nodes),leaves=int(leaves),precision=int(precision),native_seconds=float(native_seconds))
        if search_path==path or ZZ(best)<cutoff:
            assert z*H*z.transpose()==matrix(ZZ,[[ZZ(best)]])
            value=QQ(best);complete=True;result=search_result
        else:
            assert z.is_zero() and ZZ(best)==cutoff
            threshold_proof=dict(input=str(search_path.name),input_sha256=engine._hash_file(search_path),
                cutoff=str(cutoff),ratio_is_heuristic_only=threshold_ratio,result=search_result,
                meaning="Complete error-bounded enumeration found no nonzero vector of squared norm < cutoff")
    except subprocess.TimeoutExpired: value=fallback
    doc=dict(start=start,end=end,dimension=n,complete_svp=complete,block_transform_rows=rows(W),
        exact_lll_preconditions_verified=True,checker_source_sha256=checker_hash,checker_binary_sha256=engine._hash_file(checker),
        coefficient_bound_method=bound_method,coefficient_bounds=list(map(str,bounds)),
        input_sha256=engine._hash_file(path),native_source_sha256=source_hash,native_binary_sha256=engine._hash_file(native),
        search_input=search_path.name,search_input_sha256=engine._hash_file(search_path),search_cutoff=str(cutoff),
        result=result,threshold_proof=threshold_proof,**fraction_fields("lower_squared",value),**fraction_fields("fallback_gs_squared",fallback),
        seconds=time.monotonic()-began)
    save(path.with_suffix(".json"),doc)
    print("OPTIMIZED_BLOCK_CERTIFICATE",folder.parent.parent.name,folder.name,start,end,complete,round(doc["seconds"],2),flush=True)
    return value


def trial(index,proposal,blocks,seconds,loaded,threshold_ratio=None):
    case_folder,P,G,_,_,_=loaded
    C=matrix(ZZ,proposal["transform_rows"]);assert abs(C.det())==1
    suffix="" if threshold_ratio is None else "_threshold"+str(threshold_ratio).replace(".","p")
    folder=case_folder/"certificates"/(proposal["transform_sha256"][:16]+suffix)
    folder.mkdir(parents=True,exist_ok=True)
    result_path=folder/"bounds.json"
    old=None
    if result_path.exists():
        old=read(result_path)
        if old.get("complete_campaign") and set(blocks).issubset(old["block_sizes"]):
            print("OPTIMIZED_BOUND_REUSED",index,proposal["label"],flush=True)
            return old
    save(folder/"basis_transform.json",dict(label=proposal["label"],transform_rows=rows(C),
        transform_sha256=proposal["transform_sha256"],source_embedding=read(case_folder/"baseline.json")["source_embedding"],
        source_embedding_sha256=read(case_folder/"baseline.json")["embedding_sha256"],unimodular_verified=True))
    Gnew=C*G*C.transpose()
    L,D,form_scale=lower_form(Gnew,folder/"lower_form.json")
    baseline=read(case_folder/"baseline.json");scale=QQ(baseline["scale_numerator"])/ZZ(baseline["scale_denominator"])
    best=QQ(L[0,0]**2*min(D))*form_scale*scale**2;gs=best;partitions=[]
    if old:
        assert old["lower_form_sha256"]==engine._hash_file(folder/"lower_form.json")
        partitions=old["partitions"]
        best=max(best,fraction(old,"lower_squared"))
    all_blocks=sorted(set(blocks)|{item["block_size"] for item in partitions})
    upper=QQ(baseline["candidate_squared_numerator"])/ZZ(baseline["candidate_squared_denominator"])
    for size in blocks:
        if any(item["block_size"]==size for item in partitions): continue
        rem=P.nrows()%size;endpoints=[0]+([rem] if rem else [])+list(range(rem+size,P.nrows()+1,size))
        values=[]
        for start,end in zip(endpoints,endpoints[1:]): values.append(block(L,D,start,end,folder,seconds,threshold_ratio))
        lower=min(values)*form_scale*scale**2;best=max(best,lower)
        partitions.append(dict(block_size=size,endpoints=endpoints,**fraction_fields("lower_squared",lower)))
        assert best<=upper, "Certified lower bound contradicts an independently checked vector"
        doc=dict(index=index,case=name(GRID[index]),label=proposal["label"],complete_campaign=size==blocks[-1],
            block_sizes=all_blocks,partitions=partitions,transform_sha256=proposal["transform_sha256"],
            transform_file_sha256=engine._hash_file(folder/"basis_transform.json"),
            lower_form_sha256=engine._hash_file(folder/"lower_form.json"),
            **fraction_fields("lower_squared",best),**fraction_fields("gs_lower_squared",gs),
            log2_lower=str(lognorm(best)),log2_gs_lower=str(lognorm(gs)),
            certificate=str(folder.relative_to(HERE)))
        save(result_path,doc)
    print("BASIS_LOWER_CERTIFIED",index,proposal["label"],doc["log2_lower"],flush=True)
    return doc


def choose(records,count):
    available=[record for record in records if record["label"]!="baseline"]
    chosen=[];used=set()
    def add(item):
        if item["transform_sha256"] not in used:
            chosen.append(item);used.add(item["transform_sha256"])
    for item in available[:2]: add(item)
    if available: add(max(available,key=lambda item:item["gs_min_log2_relative"]))
    for prefix in ("sieve_","primal_bkz","dual_bkz","legacy_"):
        matches=[record for record in available if record["label"].startswith(prefix)]
        if matches: add(matches[0])
        if len(chosen)>=count: break
    for item in available:
        if len(chosen)>=count: break
        add(item)
    return chosen[:count]


def case(index,count,blocks,seconds,labels=None,threshold_ratio=None):
    loaded=load(index);folder=loaded[0];ranking=read(folder/"proposal_ranking.json")
    selected=choose(ranking["records"],count)
    if labels:
        selected=[]
        for prefix in labels.split(","):
            matches=[record for record in ranking["records"] if record["label"].startswith(prefix)]
            if matches: selected.append(matches[0])
    results=[]
    for proposal in selected:
        try: results.append(trial(index,proposal,blocks,seconds,loaded,threshold_ratio))
        except Exception as error:
            save(folder/("failure_"+proposal["transform_sha256"][:16]+".json"),dict(label=proposal["label"],error=repr(error)))
            print("BASIS_CERTIFICATE_FAILED",index,proposal["label"],repr(error),flush=True)
    # Keep every completed earlier trial when expanding the search campaign.
    by_certificate={result["certificate"]:result for result in results}
    for path in (folder/"certificates").glob("*/bounds.json"):
        result=read(path)
        if result.get("complete_campaign"): by_certificate[result["certificate"]]=result
    results=list(by_certificate.values())
    baseline_path=HERE/read(folder/"baseline.json")["baseline_bounds"]
    old=read(baseline_path);best=fraction(old,"lower_squared");winner=None
    for result in results:
        if fraction(result,"lower_squared")>best: best=fraction(result,"lower_squared");winner=result
    if winner:
        C=matrix(ZZ,read(HERE/winner["certificate"]/"basis_transform.json")["transform_rows"])
        R=C*loaded[1];assert abs(C.det())==1
        scale=QQ(read(folder/"baseline.json")["scale_numerator"])/ZZ(read(folder/"baseline.json")["scale_denominator"])
        save(folder/"best_basis.json",dict(rows=rows(R),transform_rows=rows(C),
            scale_numerator=str(scale.numerator()),scale_denominator=str(scale.denominator()),
            source_embedding=read(folder/"baseline.json")["source_embedding"],label=winner["label"]))
    else:
        save(folder/"best_basis.json",dict(rows=rows(loaded[1]),transform_rows=rows(identity_matrix(ZZ,loaded[1].nrows())),
            scale_numerator=read(folder/"baseline.json")["scale_numerator"],scale_denominator=read(folder/"baseline.json")["scale_denominator"],
            source_embedding=read(folder/"baseline.json")["source_embedding"],label="baseline"))
    summary=dict(index=index,case=name(GRID[index]),proposed_bases=ranking["count"],certified_trials=len(results),
        selected_label=winner["label"] if winner else "baseline",selected_certificate=winner["certificate"] if winner else None,
        baseline_bounds=str(baseline_path.relative_to(HERE)),baseline_sha256=engine._hash_file(baseline_path),
        best_basis_sha256=engine._hash_file(folder/"best_basis.json"),
        **fraction_fields("baseline_lower_squared",fraction(old,"lower_squared")),**fraction_fields("best_lower_squared",best),
        baseline_lower_log2=str(lognorm(fraction(old,"lower_squared"))),best_lower_log2=str(lognorm(best)),
        improvement_bits=str(lognorm(best)-lognorm(fraction(old,"lower_squared"))),
        tested_certificates=[result["certificate"] for result in results],global_basis_optimality_proven=False)
    save(folder/"result.json",summary)
    print("BEST_TESTED_LOWER_BOUND",index,summary["selected_label"],summary["improvement_bits"],flush=True)


def selftest():
    folder=DEST/"selftest";folder.mkdir(parents=True,exist_ok=True)
    H=matrix(ZZ,[[5,3],[3,2]])
    R,W=gram_lll(H);assert W*H*W.transpose()==R
    L=matrix(ZZ,[[ZZ(2)**80,0],[ZZ(2)**79,ZZ(2)**80]])
    D=[ZZ(7),ZZ(3)]
    lower=block(L,D,0,2,folder,30)
    gram=L*diagonal_matrix(ZZ,D)*L.transpose()
    # The form is 7*(x+y/2)^2+3*y^2, multiplied by 2^160.
    # y=0 costs >=7, |y|=1 costs >=19/4, and |y|>=2 costs >=12.
    assert lower==QQ(19)*ZZ(2)**158
    # The frozen evaluator also proves absence below an integer cutoff. A zero
    # output vector is a sentinel, never a claimed witness of that norm.
    original=folder/"block_0_2.input";lines=original.read_text().splitlines()
    binary,_=build_binary()
    for cutoff,should_find in ((ZZ(lower)//2,False),(ZZ(lower)+1,True)):
        lines[0]=f"FAST-LIFT-1 2 1 1 {cutoff}";lines[2]="0 0"
        test_path=folder/("threshold_find.input" if should_find else "threshold_empty.input")
        test_path.write_text("\n".join(lines)+"\n")
        answer=subprocess.run([str(binary),str(test_path)],text=True,capture_output=True,check=True).stdout.splitlines()
        best=ZZ(answer[0].split()[1]);z=list(map(ZZ,answer[1].split()))
        assert best==(lower if should_find else cutoff)
        assert any(z)==should_find
    print("GRAM_ONLY_BOUND_SELFTEST_OK",flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--index",type=int);parser.add_argument("--count",type=int,default=4)
    parser.add_argument("--blocks",default="55");parser.add_argument("--seconds",type=float,default=180)
    parser.add_argument("--labels")
    parser.add_argument("--threshold-ratio",type=float)
    parser.add_argument("--selftest",action="store_true");args=parser.parse_args();DEST.mkdir(exist_ok=True)
    if args.selftest: selftest();return
    blocks=list(map(int,args.blocks.split(",")))
    if args.index is not None: case(args.index,args.count,blocks,args.seconds,args.labels,args.threshold_ratio);return
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        command=[sys.executable,str(Path(__file__).resolve()),"--count",str(args.count),"--blocks",args.blocks,"--seconds",str(args.seconds)]
        if args.labels: command += ["--labels",args.labels]
        if args.threshold_ratio is not None: command += ["--threshold-ratio",str(args.threshold_ratio)]
        futures=[pool.submit(subprocess.run,command+["--index",str(index)],check=True) for index in INDICES]
        for future in concurrent.futures.as_completed(futures): future.result()


if __name__=="__main__": main()
