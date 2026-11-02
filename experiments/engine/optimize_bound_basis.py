"""Propose equivalent projected bases specifically for strong lower bounds.

Numerical GS/block-volume scores rank proposals only. Every stored proposal is
an exact integer unimodular transform of the unchanged full projection basis.
"""
import os
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key]="1"
import argparse
import concurrent.futures
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import time

import numpy as np
from sage.all import ZZ,QQ,matrix,identity_matrix
from fpylll import IntegerMatrix,GSO,LLL,BKZ,FPLLL
from projected_minima import HERE,OUT,GRID,name,read,save,rows
from shared import exact_minima as engine

DEST=HERE/"basis_bound_optimization"
INDICES=(1,2,3,8,9)
FPLLL.set_precision(256)


def load(index):
    case=name(GRID[index]);folder=DEST/case;folder.mkdir(parents=True,exist_ok=True)
    source=OUT/case/"full_strong"/"embedding.json";emb=read(source)
    old_audit=read(HERE/"lambda2_intervals"/case/"basis_audit.json")
    assert old_audit["passed"] and old_audit["embedding_sha256"]==engine._hash_file(source)
    P=matrix(ZZ,emb["rows"])[1:,1:]
    scale=QQ(emb["scale_numerator"])/ZZ(emb["scale_denominator"])
    baseline_path=HERE/"lambda2_intervals"/case/"bounds.json"
    accepted=read(HERE/"lambda2_intervals"/case/"audit.json")
    assert accepted["passed"] and accepted["bounds_sha256"]==engine._hash_file(baseline_path)
    baseline=read(baseline_path)
    candidate=read(OUT/"sieving"/case/"best_candidate.json")
    upper=QQ(candidate["physical_squared_norm_numerator"])/ZZ(candidate["physical_squared_norm_denominator"])
    G=P*P.transpose();power=max(abs(z) for z in G.list()).nbits()
    numeric=np.array([[float(QQ(z)/ZZ(2)**power) for z in row] for row in G.rows()])
    bits=max(abs(z).nbits() for z in P.list());shift=bits-256
    model=P.apply_map(lambda z:z>>shift if shift>=0 else z<<(-shift))
    target=float(upper/scale**2/ZZ(2)**power)
    if not (folder/"primal_model.json").exists():
        save(folder/"primal_model.json",dict(rows=rows(model)))
        inverse=np.linalg.inv(numeric);chol=np.linalg.cholesky((inverse+inverse.T)/2)
        normalizer=2.**128/np.max(np.abs(chol))
        dual=matrix(ZZ,[[int(round(float(chol[i,j])*normalizer)) if j<=i else 0 for j in range(P.nrows())]
                        for i in range(P.nrows())])
        assert all(dual[i,i]>0 for i in range(P.nrows()))
        save(folder/"dual_model.json",dict(rows=rows(dual)))
    save(folder/"baseline.json",dict(index=index,case=case,dimension=P.nrows(),
        embedding_sha256=engine._hash_file(source),source_embedding=str(source.relative_to(HERE)),
        baseline_bounds_sha256=engine._hash_file(baseline_path),baseline_bounds=str(baseline_path.relative_to(HERE)),
        scale_numerator=str(scale.numerator()),scale_denominator=str(scale.denominator()),
        candidate_squared_numerator=str(upper.numerator()),candidate_squared_denominator=str(upper.denominator()),
        normalized_gram_power=int(power),normalized_candidate_squared=target))
    return folder,P,G,numeric,model,target


def score(G,target,block=55):
    chol=np.linalg.cholesky((G+G.T)/2)
    diagonal=np.diag(chol)**2;d=len(diagonal)
    if not np.isfinite(diagonal).all() or min(diagonal)<=0: raise ValueError("Bad numerical profile")
    remainder=d%block;endpoints=[0]+([remainder] if remainder else [])+list(range(remainder+block,d+1,block))
    estimates=[]
    for start,end in zip(endpoints,endpoints[1:]):
        k=end-start
        loggh=(math.lgamma(k/2+1)-k/2*math.log(math.pi)+.5*sum(np.log(diagonal[start:end])))/k
        logrow=.5*math.log(min(np.sum(chol[start:end,start:end]**2,axis=1)))
        estimates.append(min(loggh,logrow))
    return dict(gs_min_log2_relative=.5*math.log2(min(diagonal)/target),
        block55_model_log2_relative=min(estimates)/math.log(2)-.5*math.log2(target),
        gs_squared_log2_profile=list(map(float,np.log2(diagonal/target))))


def size_reduce(C,model):
    A=IntegerMatrix.from_matrix([[int(z) for z in row] for row in (C*model).rows()])
    U=IntegerMatrix.identity(C.nrows())
    M=GSO.Mat(A,U=U,float_type="mpfr",flags=GSO.INT_GRAM)
    assert M.update_gso()
    LLL.Reduction(M,delta=.999,eta=.5001).size_reduction(0,C.nrows())
    V=matrix(ZZ,[[U[i,j] for j in range(U.ncols)] for i in range(U.nrows)])
    assert abs(V.det())==1
    return V*C


def pivot_order(G,seed=0,spread=1):
    rng=random.Random(seed);remaining=list(range(len(G)));inverse=np.linalg.inv(G);chosen=[]
    while remaining:
        order=np.argsort(np.diag(inverse));limit=min(spread,len(order))
        selected=int(order[rng.randrange(limit)])
        chosen.append(remaining.pop(selected))
        keep=[j for j in range(len(inverse)) if j!=selected]
        if keep:
            col=inverse[keep,selected].copy()
            inverse=inverse[np.ix_(keep,keep)]-np.outer(col,col)/inverse[selected,selected]
            inverse=(inverse+inverse.T)/2
    return list(reversed(chosen))


def seysen(G,moves=1200):
    H=np.linalg.inv(G);G=G.copy();n=len(G);C=identity_matrix(ZZ,n);history=[]
    for step in range(moves):
        gd=np.diag(G);hd=np.diag(H)
        proposal=.5*(H/hd[:,None]-G/gd[None,:])
        coefficient=np.rint(proposal);np.fill_diagonal(coefficient,0)
        coefficient=np.clip(coefficient,-10000,10000)
        product=hd[:,None]*gd[None,:]
        linear=G*hd[:,None]-H*gd[None,:]
        gain=-(2*coefficient*linear+2*coefficient**2*product)
        np.fill_diagonal(gain,0)
        i,j=np.unravel_index(np.argmax(gain),gain.shape);mu=int(coefficient[i,j])
        if mu==0 or gain[i,j]<1e-8: break
        old_i,old_j=G[i,:].copy(),G[j,:].copy();old_gii=G[i,i];old_gij=G[i,j];old_gjj=G[j,j]
        G[i,:]=old_i+mu*old_j;G[:,i]=G[i,:];G[i,i]=old_gii+2*mu*old_gij+mu*mu*old_gjj
        old_i,old_j=H[i,:].copy(),H[j,:].copy();old_hii=H[i,i];old_hij=H[i,j];old_hjj=H[j,j]
        H[j,:]=old_j-mu*old_i;H[:,j]=H[j,:];H[j,j]=old_hjj-2*mu*old_hij+mu*mu*old_hii
        C[i]=C.row(i)+mu*C.row(j)
        history.append((int(i),int(j),mu))
    return C,len(history)


def bkz_family(index,dual,target_block):
    folder=DEST/name(GRID[index]);model=matrix(ZZ,read(folder/("dual_model.json" if dual else "primal_model.json"))["rows"])
    A=IntegerMatrix.from_matrix([[int(z) for z in row] for row in model.rows()]);d=model.nrows()
    total=identity_matrix(ZZ,d);family="dual" if dual else "primal";checkpoints=[]
    for block in (20,30,40,50,60,70,80):
        if block>target_block: break
        U=IntegerMatrix.identity(d)
        parameters=BKZ.Param(block_size=min(block,d),strategies=BKZ.DEFAULT_STRATEGY,max_loops=2,
            max_time=60,flags=BKZ.MAX_LOOPS|BKZ.MAX_TIME|BKZ.AUTO_ABORT)
        BKZ.reduction(A,parameters,U=U,float_type="mpfr",precision=256)
        V=matrix(ZZ,[[U[i,j] for j in range(d)] for i in range(d)])
        total=V*total
        current=matrix(ZZ,[[A[i,j] for j in range(A.ncols)] for i in range(d)])
        assert total*model==current and abs(total.det())==1
        checkpoints.append(dict(block=block,transform_rows=rows(total)))
        save(folder/f"{family}_bkz_checkpoints.json",dict(family=family,checkpoints=checkpoints))
        print("BKZ_PROPOSAL_CHECKPOINT",index,family,block,flush=True)


def propose(index,quick=False):
    folder,P,G,numeric,model,target=load(index);d=P.nrows();records=[];seen=set()
    existing=read(folder/"proposals.json")["records"] if (folder/"proposals.json").exists() else []
    def add(label,C,reduce=True):
        try:
            C=C.change_ring(ZZ)
            if reduce: C=size_reduce(C,model)
            assert abs(C.det())==1
            key=engine._digest(rows(C))
            if key in seen: return
            W=np.array([[float(z) for z in row] for row in C.rows()])
            candidate=W@numeric@W.T
            quality=score(candidate,target)
            entry=dict(label=label,transform_rows=rows(C),transform_sha256=key,
                exact_unimodular=True,score_kind="Heuristic ranking only; lower bounds are certified separately",**quality)
            records.append(entry);seen.add(key)
            save(folder/"proposals.json",dict(index=index,case=name(GRID[index]),records=records))
        except (ValueError,ArithmeticError,AssertionError,np.linalg.LinAlgError) as error:
            print("PROPOSAL_SKIPPED",index,label,type(error).__name__,flush=True)
    add("baseline",identity_matrix(ZZ,d),False)
    add("reverse",identity_matrix(ZZ,d).matrix_from_rows(list(reversed(range(d)))))
    for seed in range(12):
        permutation=pivot_order(numeric,seed,1 if seed==0 else 2 if seed<6 else 4)
        add(f"dual_pivot_{seed}",identity_matrix(ZZ,d).matrix_from_rows(permutation))
    C,moves=seysen(numeric)
    add(f"seysen_{moves}",C)
    # Insert the actual saved G6K lifts while optimizing the lower-bound
    # profile, rather than the unrelated full-enumeration node estimate.
    stage=OUT/name(GRID[index])/"full_strong"/"certificate"/"pass1"
    prepared=engine._read_document(stage/"prepared.json")
    lift_rows={}
    for path in sorted(stage.glob("sieve*.json")):
        doc=read(path)
        if doc.get("prepared_sha256")!=prepared["payload_sha256"]: continue
        for item in doc.get("public_preconditioning_rows",[]):
            coefficients=tuple(map(ZZ,item["quotient_coefficients"]))
            lift_rows[(int(item["index"]),coefficients)]=None
    point=read(OUT/"sieving"/name(GRID[index])/"best_candidate.json")
    lift_rows[(0,tuple(map(ZZ,point["quotient_coefficients"])))]=None
    C=identity_matrix(ZZ,d);quality=score(numeric,target)["block55_model_log2_relative"]
    for round_index in range(8):
        inverse=C.inverse().change_ring(ZZ);W=np.array([[float(z) for z in row] for row in C.rows()])
        current=W@numeric@W.T;best=None
        for (position,coefficients) in lift_rows:
            y=(matrix(ZZ,[coefficients])*inverse).row(0)
            units=[j for j in range(position,d) if abs(y[j])==1]
            for removed in units[-2:]:
                E=identity_matrix(ZZ,d);E[removed]=y
                permutation=list(range(d));permutation.insert(position,permutation.pop(removed))
                E=E.matrix_from_rows(permutation)
                V=np.array([[float(z) for z in row] for row in E.rows()])
                try: estimate=score(V@current@V.T,target)["block55_model_log2_relative"]
                except np.linalg.LinAlgError: continue
                if estimate>quality+1e-5 and (best is None or estimate>best[0]): best=(estimate,E,position,removed)
        if best is None: break
        quality,E,position,removed=best;C=size_reduce(E*C,model)
        add(f"sieve_bound_insert_{round_index+1}",C,False)
        print("BOUND_INSERTION_IMPROVED",index,round_index+1,quality,position,removed,flush=True)
    for old in existing:
        add(old["label"],matrix(ZZ,old["transform_rows"]),False)
    # Existing high-precision public basis experiments are useful candidates
    # even when their enumeration-cost model was worse than the original.
    embedding=read(OUT/name(GRID[index])/"full_strong"/"embedding.json")
    source=HERE/embedding["source_preparation"]
    quotient=engine._read_document(source/"quotient_lll.json")
    for sibling in source.parent.parent.glob(source.parent.name+"-*"):
        path=sibling/"pass1"/"quotient_lll.json"
        if not path.exists(): continue
        other=engine._read_document(path)
        C=matrix(ZZ,other["transform_rows"])*matrix(ZZ,quotient["transform_rows"]).inverse()
        if C.denominator()!=1: continue
        C=C.change_ring(ZZ)
        assert C*P==matrix(ZZ,other["reduced_rows"])
        add("legacy_"+sibling.name,C)
    if not quick:
        for dual in (True,False):
            family="dual" if dual else "primal"
            try:
                subprocess.run([sys.executable,str(Path(__file__).resolve()),"--index",str(index),
                                "--family",family,"--target-block","60"],check=True,timeout=180)
            except (subprocess.TimeoutExpired,subprocess.CalledProcessError):
                print("BKZ_FAMILY_PARTIAL",index,family,flush=True)
            path=folder/f"{family}_bkz_checkpoints.json"
            if not path.exists(): continue
            for checkpoint in read(path)["checkpoints"]:
                V=matrix(ZZ,checkpoint["transform_rows"])
                C=V.inverse().transpose().change_ring(ZZ) if dual else V
                add(f"{family}_bkz{checkpoint['block']}",C)
                if dual:
                    add(f"{family}_bkz{checkpoint['block']}_reverse",C.matrix_from_rows(list(reversed(range(d)))))
    ranked=sorted(records,key=lambda item:(item["block55_model_log2_relative"],item["gs_min_log2_relative"]),reverse=True)
    save(folder/"proposal_ranking.json",dict(index=index,case=name(GRID[index]),count=len(ranked),records=ranked))
    for record in ranked[:5]:
        print("RANKED_BOUND_BASIS",index,record["label"],record["block55_model_log2_relative"],record["gs_min_log2_relative"],flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--index",type=int)
    parser.add_argument("--quick",action="store_true");parser.add_argument("--family",choices=("dual","primal"))
    parser.add_argument("--target-block",type=int,default=60);args=parser.parse_args();DEST.mkdir(exist_ok=True)
    if args.family: bkz_family(args.index,args.family=="dual",args.target_block);return
    if args.index is not None: propose(args.index,args.quick);return
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        futures=[]
        for index in INDICES:
            command=[sys.executable,str(Path(__file__).resolve()),"--index",str(index)]
            if args.quick: command.append("--quick")
            futures.append(pool.submit(subprocess.run,command,check=True))
        for future in concurrent.futures.as_completed(futures): future.result()


if __name__=="__main__": main()
