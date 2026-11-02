"""Change block-prefix subspaces using public short vectors and local swaps."""
import os
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"): os.environ[key]="1"
import argparse
import concurrent.futures
from copy import copy
from pathlib import Path
import subprocess
import sys
import numpy as np
from sage.all import ZZ,matrix,identity_matrix
from optimize_bound_basis import load,score,size_reduce,INDICES
from projected_minima import HERE,OUT,GRID,name,read,save,rows
from shared import exact_minima as engine


def case(index):
    folder,P,G,numeric,model,target=load(index);d=P.nrows();cut=d-55
    old=read(folder/"proposals.json")["records"];seen={item["transform_sha256"] for item in old}
    candidates=[]
    def consider(label,C):
        W=np.array([[float(z) for z in row] for row in C.rows()])
        try: quality=score(W@numeric@W.T,target)
        except np.linalg.LinAlgError: return
        candidates.append((quality["block55_model_log2_relative"],label,C,quality))
    identity=identity_matrix(ZZ,d)
    # Moving a later vector into the prefix changes the projected tail lattice.
    for start in range(max(0,cut-5),cut):
        for end in range(cut,min(d,cut+25)):
            C=copy(identity);C.swap_rows(start,end)
            consider(f"boundary_swap_{start}_{end}",C)
    chol=np.linalg.cholesky(numeric);mu=chol/np.diag(chol)[None,:]
    vectors=[];stage=OUT/name(GRID[index])/"full_strong"/"certificate"/"pass1"
    prepared=engine._read_document(stage/"prepared.json")
    for path in sorted(stage.glob("sieve*.json")):
        doc=read(path)
        if doc.get("prepared_sha256")!=prepared["payload_sha256"]: continue
        for item in doc.get("public_preconditioning_rows",[]):
            z=list(map(ZZ,item["quotient_coefficients"]));position=int(item["index"])
            for i in reversed(range(position)):
                z[i]=ZZ(-round(sum(float(z[j])*mu[j,i] for j in range(i+1,d))))
            vectors.append(z)
    point=read(OUT/"sieving"/name(GRID[index])/"best_candidate.json")
    vectors.append(list(map(ZZ,point["quotient_coefficients"])))
    for number,z in enumerate(vectors):
        for position in sorted({0,max(0,cut-1),max(0,cut-3)}):
            units=[j for j in range(position,d) if abs(z[j])==1]
            for removed in units[-4:]:
                C=copy(identity);C[removed]=z
                permutation=list(range(d));permutation.insert(position,permutation.pop(removed))
                C=C.matrix_from_rows(permutation)
                consider(f"prefix_lift_{number}_{position}_{removed}",C)
    candidates.sort(key=lambda item:item[0],reverse=True)
    selected=[]
    for prefix,limit in (("boundary_swap_",4),("prefix_lift_",4)):
        selected.extend([item for item in candidates if item[1].startswith(prefix)][:limit])
    for _,label,C,quality in selected:
        C=size_reduce(C,model);assert abs(C.det())==1
        key=engine._digest(rows(C))
        if key in seen: continue
        seen.add(key)
        W=np.array([[float(z) for z in row] for row in C.rows()]);quality=score(W@numeric@W.T,target)
        old.append(dict(label=label,transform_rows=rows(C),transform_sha256=key,exact_unimodular=True,
                        score_kind="Heuristic only; requires exact block certification",**quality))
        print("PREFIX_SUBSPACE_PROPOSAL",index,label,quality["block55_model_log2_relative"],flush=True)
    save(folder/"proposals.json",dict(index=index,case=name(GRID[index]),records=old))
    ranked=sorted(old,key=lambda item:(item["block55_model_log2_relative"],item["gs_min_log2_relative"]),reverse=True)
    save(folder/"proposal_ranking.json",dict(index=index,case=name(GRID[index]),count=len(ranked),records=ranked))
    save(folder/"local_search.json",dict(numerically_screened=len(candidates),saved_candidates=len(selected),
        count_is_heuristic_only=True))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--index",type=int);args=parser.parse_args()
    if args.index is not None: case(args.index)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
            futures=[pool.submit(subprocess.run,[sys.executable,str(Path(__file__).resolve()),"--index",str(index)],check=True)
                     for index in INDICES]
            for future in concurrent.futures.as_completed(futures): future.result()
