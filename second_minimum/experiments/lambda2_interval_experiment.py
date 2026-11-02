"""Certified lambda2 intervals and an independent test of min(GH, jet).

For the unfinished cases, a small exact matrix certificate proves G >= LDL^T.
Contiguous blocks of this triangular lower form give rigorous SVP lower bounds
without solving the full projected SVP. An original-lattice witness gives the
upper bound. Exact previously audited lambda2 values remain point intervals.
"""
import argparse
import concurrent.futures
import csv
from datetime import datetime, timedelta, timezone
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

import mpmath as mp
import numpy as np
from sage.all import diagonal_matrix, identity_matrix, RealIntervalField
from projected_minima import HERE, OUT, GRID, build, variant, name, read, save, rows, norm2, matrix, vector, ZZ, QQ, gcd
from shared import exact_minima as engine

DEST = HERE / "lambda2_intervals"
mp.mp.dps = 80


def fraction_fields(prefix, value):
    value = QQ(value)
    return {prefix+"_numerator": str(value.numerator()), prefix+"_denominator": str(value.denominator())}


def fraction(doc, prefix):
    return QQ(doc[prefix+"_numerator"]) / ZZ(doc[prefix+"_denominator"])


def mpq(value):
    value = QQ(value)
    return mp.mpf(int(value.numerator())) / int(value.denominator())


def lognorm(value):
    return mp.log(mpq(value), 2)/2


def check_member(v, data):
    assert all(v[j] % data["scales"][j] == 0 for j in range(len(v)))
    raw = vector(ZZ, [v[j]//data["scales"][j] for j in range(len(v))])
    affine = raw*data["forward"].transpose()
    exponents = variant(data,"B")[1]
    assert all(affine[j] % data["N"]**exponents[j] == 0 for j in range(len(exponents)))


def verified_quotient(index, data):
    """Check the full public basis and every unimodular change of coordinates."""
    from audit_results import measurement
    case = name(GRID[index])
    folder = OUT / case / "full_strong"
    emb = read(folder / "embedding.json")
    source = HERE / emb["source_preparation"]
    old = engine._read_document(source / "prepared.json")
    assert engine._hash_file(source / "prepared.json") == emb["source_preparation_sha256"]
    measured_path = HERE / emb["original_measurement"]
    assert engine._hash_file(measured_path) == emb["original_measurement_sha256"]
    _, B = measurement(data, measured_path, "B")
    g = gcd(B.list())
    expected = B.apply_map(lambda z:z//g)
    preconditioning = source.parent / "preconditioning.json"
    if preconditioning.exists():
        W = matrix(ZZ,read(preconditioning)["transform_rows"])
        assert abs(W.det()) == 1
        expected = W*expected
    assert engine._read_document(source/"input_basis.json")["rows"] == rows(expected)
    ambient = engine._read_document(source/"ambient_lll.json")
    quotient = engine._read_document(source/"quotient_lll.json")
    pre, U = matrix(ZZ,ambient["preconditioner_rows"]), matrix(ZZ,ambient["transform_rows"])
    warm, red = matrix(ZZ,ambient["lll_input_rows"]), matrix(ZZ,ambient["reduced_rows"])
    assert pre*expected == warm and U*warm == red
    assert abs(pre.det()) == abs(U.det()) == 1
    order = sorted(range(red.nrows()),key=lambda i:norm2(red.row(i)))
    r = red.row(order[0]); a = norm2(r)
    physical_root=variant(data,"B")[3]
    assert g*r in (physical_root,-physical_root)
    T0 = matrix(ZZ,[red.row(i) for i in order[1:]])
    projected = a*T0-matrix(ZZ,[z*r for z in T0*r])
    h = gcd(projected.list())
    projected = projected.apply_map(lambda z:z//h)
    V = matrix(ZZ,quotient["transform_rows"])
    assert abs(V.det()) == 1 and quotient["input_rows"] == rows(projected)
    P = V*projected
    assert quotient["reduced_rows"] == rows(P)
    T = matrix(ZZ,old["tail"])
    assert T == V*T0-matrix(ZZ,[ZZ(k)*r for k in old["root_shifts"]])
    assert list(map(str,r)) == old["root"] and list(map(str,T*r)) == old["inner"]
    embedded = matrix(ZZ,emb["rows"])
    assert embedded[0,0] == 1 and all(z == 0 for z in embedded[0,1:].list())
    assert all(z == 0 for z in embedded[1:,0].list()) and embedded[1:,1:] == P
    scale = QQ(emb["scale_numerator"])/ZZ(emb["scale_denominator"])
    assert scale == QQ(g*h)/a
    assert P.nrows()+1 == B.nrows()
    save(DEST/case/"basis_audit.json",dict(passed=True,index=index,case=case,
        full_public_lattice_verified=True,unimodular_projection_equivalence_verified=True,
        embedding_sha256=engine._hash_file(folder/"embedding.json"),
        source_prepared_sha256=engine._hash_file(source/"prepared.json"),
        measured_sha256=engine._hash_file(measured_path),projected_rank=P.nrows()))
    return P, scale, old, g, T


def lower_form(P, path):
    """Exact preconditioned diagonal dominance proves G >= c*L*D*L^T."""
    G = P*P.transpose()
    d = G.nrows(); power = max(abs(z) for z in G.list()).nbits()
    scaled = np.array([[float(QQ(z)/ZZ(2)**power) for z in row] for row in G.rows()])
    chol = np.linalg.cholesky(scaled)
    unit = chol/np.diag(chol)[None,:]
    inverse = np.linalg.inv(unit)
    bits = 80; den = ZZ(2)**bits
    L = matrix(ZZ,d,d,lambda i,j: den if i==j else int(round(float(unit[i,j])*int(den))) if i>j else 0)
    M = matrix(ZZ,d,d,lambda i,j: den if i==j else int(round(float(inverse[i,j])*int(den))) if i>j else 0)
    for safety_bits in (20,16,12,8,4):
        factor = 1-2.**(-safety_bits)
        D = [ZZ(math.floor(float(chol[i,i]**2)*factor*int(den))) for i in range(d)]
        assert min(D)>0
        H = L*diagonal_matrix(ZZ,D)*L.transpose()
        residual = den**3*G-ZZ(2)**power*H
        checked = M*residual*M.transpose()
        margin = min(checked[i,i]-sum(abs(checked[i,j]) for j in range(d) if j!=i) for i in range(d))
        if margin>0:
            break
    else:
        raise RuntimeError("No positive exact lower-form certificate")
    assert all(100*abs(L[i,j]) <= 51*den for i in range(d) for j in range(i))
    assert all(100*(den**2*D[i]+L[i,i-1]**2*D[i-1]) >= 99*den**2*D[i-1] for i in range(1,d))
    doc = dict(schema="exact-triangular-lower-form-v1",fraction_bits=bits,gram_scale_power=int(power),
        L_numerators=rows(L),M_numerators=rows(M),D_numerators=list(map(str,D)),
        safety_fraction=f"1-2^-{safety_bits}",positive_diagonal_dominance_margin=str(margin),
        argument="M*(2^(3b)*G-2^s*L*diag(D)*L^T)*M^T is strictly diagonally dominant positive; L and M have diagonal 2^b")
    save(path,doc)
    return L,D,QQ(ZZ(2)**power)/den**3,doc


def ceil_sqrt(value):
    value=QQ(value)
    integer=(value.numerator()+value.denominator()-1)//value.denominator()
    r=integer.isqrt()
    return r+(r*r<integer)


def block_svp(L,D,start,end,folder,seconds):
    from fast_search import build_binary
    block = L[start:end,start:end]
    H = block*diagonal_matrix(ZZ,D[start:end])*block.transpose()
    d=H.nrows(); den=L[0,0]
    initial=min(H[i,i] for i in range(d)); initial_i=next(i for i in range(d) if H[i,i]==initial)
    y=[int(i==initial_i) for i in range(d)]
    bounds=[ZZ(0)]*d
    for i in reversed(range(d)):
        weighted=sum(abs(block[j,i])*bounds[j] for j in range(i+1,d))
        bounds[i]=ceil_sqrt(QQ(initial-1)/(den**2*D[start+i]))+(weighted+den-1)//den
    bound_method="exact_triangular_integer_bound"
    if max(bounds)>=2**29:
        inverse=H.inverse()
        bounds=[ceil_sqrt((initial-1)*inverse[i,i]) for i in range(d)]
        bound_method="exact_inverse_gram_cauchy_schwarz"
    assert max(bounds)<2**29
    binary, source_hash=build_binary()
    filename=folder/f"block_{start}_{end}.input"
    with filename.open("w") as stream:
        stream.write(f"FAST-LIFT-1 {d} 1 1 {initial}\n"+" ".join(["0"]*d)+"\n"+" ".join(map(str,y))+"\n")
        for row in H.rows(): stream.write(" ".join(map(str,row))+"\n")
    began=time.monotonic(); complete=False; result=None
    try:
        run=subprocess.run([str(binary),str(filename)],capture_output=True,text=True,check=True,timeout=seconds)
        tag,best,nodes,leaves,native_seconds,precision=run.stdout.splitlines()[0].split()
        assert tag=="COMPLETE"
        witness=vector(ZZ,list(map(int,run.stdout.splitlines()[1].split())))
        assert witness*H*witness==ZZ(best)>0
        complete=True
        result=dict(squared_minimum=best,witness=list(map(str,witness)),nodes=int(nodes),
                    leaves=int(leaves),precision=int(precision),native_seconds=float(native_seconds))
        lower=ZZ(best)
    except subprocess.TimeoutExpired:
        lower=den**2*min(D[start:end])
    doc=dict(start=start,end=end,dimension=d,complete_svp=complete,lower_squared_norm=str(lower),
        input_sha256=engine._hash_file(filename),native_source_sha256=source_hash,
        native_binary_sha256=engine._hash_file(binary),coefficient_bound_method=bound_method,
        coefficient_bounds=list(map(str,bounds)),elapsed_seconds=time.monotonic()-began,result=result)
    save(filename.with_suffix(".json"),doc)
    print("BLOCK_LOWER_BOUND",folder.name,start,end,complete,round(doc["elapsed_seconds"],2),flush=True)
    return lower,doc


def case(index, block_sizes, seconds):
    cfg=GRID[index]; case_name=name(cfg); folder=DEST/case_name; folder.mkdir(parents=True,exist_ok=True)
    geo=read(OUT/case_name/"geometry.json")
    certified=next(row for row in read(HERE/"certified_minima.json") if row["case"]==case_name and row["variant"]=="B")
    root=vector(ZZ,geo["root"]); root_sq=norm2(root)
    result=dict(index=index,case=case_name,application=cfg[0].upper(),t=cfg[1],m=cfg[2],ell=cfg[3],
                dimension=geo["ambient_dimension"],prediction_log2=str(min(mp.mpf(geo["log2_gh_prediction"]),mp.mpf(geo["log2_jet_prediction"]))),
                prediction_exact_ball_log2=str(min(mp.mpf(geo["log2_gh_prediction"]),mp.mpf(geo["log2_jet_exact_ball_prediction"]))),
                root_sq=str(root_sq))
    if certified["exact"]:
        source=HERE/"results"/case_name/"B.certificate.json"
        assert certified["certificate_sha256"]==engine._hash_file(source)
        certificate=read(source)
        assert certificate["exact"] and ZZ(certificate["lambda1_sq"])==root_sq
        squared=ZZ(certificate["lambda2_sq"])
        data=build(cfg)
        witness=vector(ZZ,certificate["lambda2_witness"])
        check_member(root,data);check_member(witness,data)
        assert norm2(witness)==squared and norm2(root)*squared-root.dot_product(witness)**2>0
        lower=upper=QQ(squared)
        result.update(exact_lambda2=True,method="previous_independently_audited_exact_lambda2",
            source_certificate=str(source.relative_to(HERE)),source_sha256=engine._hash_file(source))
    else:
        data=build(cfg)
        P,scale,old,g,T=verified_quotient(index,data)
        candidate_path=OUT/"sieving"/case_name/"best_candidate.json"
        candidate=read(candidate_path)
        base=engine._read_document(OUT/case_name/"full_strong"/"certificate"/"pass1"/"prepared.json")
        assert candidate["prepared_sha256"]==base["payload_sha256"]
        y=vector(ZZ,candidate["quotient_coefficients"])
        preimage=g*(y*T)
        check_member(root,data);check_member(preimage,data)
        dot=preimage.dot_product(root)
        k=(2*dot+root_sq)//(2*root_sq)
        witness=preimage-k*root
        assert norm2(root)*norm2(witness)-root.dot_product(witness)**2>0
        projection=preimage.change_ring(QQ)-QQ(dot)/root_sq*root
        assert norm2(y*P)*scale**2==projection.dot_product(projection)==fraction(candidate,"physical_squared_norm")
        upper=QQ(max(root_sq,norm2(witness)))
        save(folder/"upper_witness.json",dict(vector=list(map(str,witness)),root=list(map(str,root)),
            squared_norm=str(norm2(witness)),candidate_sha256=engine._hash_file(candidate_path),
            membership_and_independence_verified=True))
        L,D,form_scale,form=lower_form(P,folder/"lower_form.json")
        best=QQ(L[0,0]**2*min(D))*scale**2*form_scale
        partitions=[]
        for size in block_sizes:
            remainder=P.nrows()%size
            endpoints=[0]+([remainder] if remainder else [])+list(range(remainder+size,P.nrows()+1,size))
            assert endpoints[-1]==P.nrows()
            bounds=[]
            for start,end in zip(endpoints,endpoints[1:]):
                value,doc=block_svp(L,D,start,end,folder,seconds)
                bounds.append(value)
            lower_partition=QQ(min(bounds))*scale**2*form_scale
            best=max(best,lower_partition)
            partitions.append(dict(block_size=size,endpoints=endpoints,
                lower_squared_numerator=str(lower_partition.numerator()),lower_squared_denominator=str(lower_partition.denominator())))
            intermediate=dict(result,exact_lambda2=False,method="exact_lower_form_and_block_svp; verified_nearest_lift_upper_witness",
                partitions=partitions,**fraction_fields("lower_squared",best),**fraction_fields("upper_squared",upper))
            save(folder/"bounds.json",intermediate)
        lower=best
        result.update(exact_lambda2=False,method="exact_lower_form_and_block_svp; verified_nearest_lift_upper_witness",
            partitions=partitions,lower_form_sha256=engine._hash_file(folder/"lower_form.json"),
            upper_witness_sha256=engine._hash_file(folder/"upper_witness.json"),
            **fraction_fields("projected_scale",scale))
    assert 0<lower<=upper
    result.update(**fraction_fields("lower_squared",lower),**fraction_fields("upper_squared",upper))
    save(folder/"bounds.json",result)
    print("LAMBDA2_INTERVAL",index,str(lognorm(lower)),str(lognorm(upper)),flush=True)


def collect():
    records=[]
    for index,cfg in enumerate(GRID):
        path=DEST/name(cfg)/"bounds.json"
        if not path.exists(): continue
        doc=read(path); lower,upper=fraction(doc,"lower_squared"),fraction(doc,"upper_squared")
        lo,hi,pred=lognorm(lower),lognorm(upper),mp.mpf(doc["prediction_log2"])
        ratio_lo,ratio_hi=mp.power(2,lo-pred),mp.power(2,hi-pred)
        row={key:doc[key] for key in ("index","case","application","t","m","ell","dimension","exact_lambda2")}
        audit_path=path.parent/"audit.json"
        audit=read(audit_path) if audit_path.exists() else {}
        row["independent_audit_passed"]=bool(audit.get("passed") and audit.get("bounds_sha256")==engine._hash_file(path))
        row.update(lower_log2=str(lo),upper_log2=str(hi),width_bits=str(hi-lo),
            bound_width_percent=str(100*(mp.sqrt(mpq(upper/lower))-1)),prediction_log2=str(pred),
            true_over_prediction_lower=str(ratio_lo),true_over_prediction_upper=str(ratio_hi),
            prediction_in_interval=lo<=pred<=hi,
            prediction_distance_bits=str(max(lo-pred,pred-hi,mp.mpf(0))))
        for epsilon in (1,5,10,25,50):
            a,b=mp.mpf(1)-mp.mpf(epsilon)/100,mp.mpf(1)+mp.mpf(epsilon)/100
            row[f"accuracy_{epsilon}pct"]=("certified_within" if ratio_lo>=a and ratio_hi<=b
                else "certified_outside" if ratio_lo>b or ratio_hi<a else "inconclusive")
        corrected=mp.mpf(doc["prediction_exact_ball_log2"])
        row.update(exact_ball_prediction_log2=str(corrected),
            exact_ball_true_over_prediction_lower=str(mp.power(2,lo-corrected)),
            exact_ball_true_over_prediction_upper=str(mp.power(2,hi-corrected)))
        records.append(row)
    if records:
        for filename,readable in (("comparison.csv",False),("comparison_readable.csv",True)):
            target=DEST/filename;temporary=target.with_name(target.name+".tmp-"+str(os.getpid()))
            with temporary.open("w",encoding="utf-8-sig",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader()
                for row in records:
                    out=dict(row)
                    if readable:
                        for key in ("lower_log2","upper_log2","width_bits","bound_width_percent","prediction_log2",
                                    "true_over_prediction_lower","true_over_prediction_upper","prediction_distance_bits",
                                    "exact_ball_prediction_log2","exact_ball_true_over_prediction_lower","exact_ball_true_over_prediction_upper"):
                            out[key]=mp.nstr(mp.mpf(out[key]),13)
                    writer.writerow(out)
            temporary.replace(target)
    save(DEST/"status.json",dict(created_kst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
        cases=len(records),exact=sum(row["exact_lambda2"] for row in records),
        independently_audited=sum(row["independent_audit_passed"] for row in records),
        role="Exact matrix and small-block SVP bounds; independent audit recorded separately. Narrow intervals, rather than prediction membership alone, support accuracy claims",
        records=records))
    print("LAMBDA2_INTERVAL_ROWS",len(records),flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--index",type=int)
    parser.add_argument("--indices")
    parser.add_argument("--blocks",default="35,45")
    parser.add_argument("--seconds",type=float,default=120)
    parser.add_argument("--collect",action="store_true")
    args=parser.parse_args();DEST.mkdir(exist_ok=True)
    if args.collect: collect();return
    if args.index is not None:
        case(args.index,list(map(int,args.blocks.split(","))),args.seconds);return
    old=read(HERE/"certified_minima.json")
    exact_cases={row["case"] for row in old if row["variant"]=="B" and row["exact"]}
    indices=sorted(range(len(GRID)),key=lambda index:(name(GRID[index]) not in exact_cases,index))
    if args.indices:
        selected=set(map(int,args.indices.split(",")))
        indices=[index for index in indices if index in selected]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures={pool.submit(subprocess.run,[sys.executable,str(Path(__file__).resolve()),"--index",str(index),
            "--blocks",args.blocks,"--seconds",str(args.seconds)],check=True):index for index in indices}
        for future in concurrent.futures.as_completed(futures):
            future.result();collect()


if __name__=="__main__": main()
