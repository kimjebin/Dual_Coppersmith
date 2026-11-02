"""Directed-interval audit of the prediction's multiplicative approximation."""
import csv
from datetime import datetime, timedelta, timezone
from math import factorial
from sage.all import RealIntervalField
from lambda2_interval_experiment import DEST,HERE,OUT,GRID,name,read,save,build,variant,fraction,ZZ,QQ
from shared import exact_minima as engine

RIF=RealIntervalField(320)


def unit_ball_volume(d):
    k=d//2
    if d%2==0:
        return RIF.pi()**k/factorial(k)
    return RIF.pi()**k*4**(k+1)*factorial(k+1)/factorial(2*k+2)


def model(data,geo):
    n=geo["ambient_dimension"];ell=data["cfg"][3];t=data["cfg"][1]
    exponents=variant(data,"B")[1]
    detlog=sum(exponents)*RIF(data["N"]).log2()+sum(t-sum(key) for key in data["keys"])*RIF(data["R"]).log2()
    rootlog=RIF(ZZ(geo["root_sq"])).log2()/2
    gh=(detlog-rootlog-unit_ball_volume(n-1).log2())/(n-1)
    jet=(RIF(ell)/(2*RIF.pi()*RIF(1).exp())).log2()/2+(RIF(ZZ(geo["jet_gram_determinant"])).log2()/2-rootlog)/ell
    return RIF(min(gh.lower(),jet.lower()),min(gh.upper(),jet.upper()))


records=[]
for index,cfg in enumerate(GRID):
    folder=DEST/name(cfg);doc=read(folder/"bounds.json");audit=read(folder/"audit.json")
    assert audit["passed"] and audit["bounds_sha256"]==engine._hash_file(folder/"bounds.json")
    data=build(cfg);geo=read(OUT/name(cfg)/"geometry.json")
    prediction=model(data,geo)
    assert abs(float(prediction.center())-float(doc["prediction_log2"]))<1e-9
    lower,upper=fraction(doc,"lower_squared"),fraction(doc,"upper_squared")
    bias=RIF((RIF(lower).log2()/2-prediction).lower(),(RIF(upper).log2()/2-prediction).upper())
    ratio=(bias*RIF(2).log()).exp()
    row=dict(index=index,case=name(cfg),prediction_log2_interval=str(prediction),
        true_over_prediction_interval=str(ratio),error_bits_interval=str(bias),
        ratio_lower_display=float(ratio.lower()),ratio_upper_display=float(ratio.upper()),
        error_bits_lower_display=float(bias.lower()),error_bits_upper_display=float(bias.upper()))
    for factor,label in ((QQ(2),"factor_2"),(QQ(14)/5,"factor_2_8"),(QQ(3),"factor_3")):
        row[label+"_certified"]=(ratio.lower()>=RIF(1/factor).upper() and ratio.upper()<=RIF(factor).lower())
    for epsilon in (1,5,10,25,50):
        lo=RIF(QQ(100-epsilon)/100);hi=RIF(QQ(100+epsilon)/100)
        row[f"accuracy_{epsilon}pct"]=("certified_within" if ratio.lower()>=lo.upper() and ratio.upper()<=hi.lower()
            else "certified_outside" if ratio.lower()>hi.upper() or ratio.upper()<lo.lower() else "inconclusive")
    records.append(row)
    print("RIGOROUS_MODEL_RATIO",index,row["ratio_lower_display"],row["ratio_upper_display"],
          "factor_2_8",row["factor_2_8_certified"],flush=True)
assert len(records)==21
summary=dict(created_kst=datetime.now(timezone(timedelta(hours=9))).isoformat(),passed=21,
    directed_interval_precision_bits=320,prediction="min(full GH with exact unit ball, paper asymptotic jet)",
    factor_2_8_certified_all=all(row["factor_2_8_certified"] for row in records),
    factor_3_certified_all=all(row["factor_3_certified"] for row in records),
    finite_sample_scope="The 21 original B lattices, one seed-0 instance per parameter setting; no claim for other instances",
    interpretation="These are lambda2/prediction enclosures. Prediction membership in a wide lambda2 interval alone does not establish accuracy.",
    records=records)
for epsilon in (1,5,10,25,50):
    key=f"accuracy_{epsilon}pct"
    summary[key]={status:sum(row[key]==status for row in records)
                  for status in ("certified_within","certified_outside","inconclusive")}
save(DEST/"model_validation.json",summary)
with (DEST/"rigorous_model_comparison.csv").open("w",encoding="utf-8-sig",newline="") as stream:
    writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
print("MODEL_ALL_CASES_FACTOR_2_8",summary["factor_2_8_certified_all"],flush=True)
