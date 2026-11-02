"""Empirical prediction comparison: exact projected minima or checked sieve candidates.

This compares point values, not the certified lambda2 lower/upper intervals.
An uncertified sieve candidate remains an upper bound on the projected minimum.
"""
import csv
from datetime import datetime, timedelta, timezone
from fractions import Fraction
import os
import mpmath as mp

from projected_minima import HERE, OUT, GRID, name, read, save
from shared import exact_minima as engine

DEST=OUT/"empirical_comparison"
mp.mp.dps=80


def squared_norm(doc):
    return Fraction(int(doc["squared_norm_numerator"]),int(doc["squared_norm_denominator"]))


def log_length(value):
    return (mp.log(value.numerator,2)-mp.log(value.denominator,2))/2


def main():
    DEST.mkdir(exist_ok=True)
    records=[]
    for index,cfg in enumerate(GRID):
        case=name(cfg);folder=OUT/case;geo_path=folder/"geometry.json";geo=read(geo_path)
        gh,jet=mp.mpf(geo["log2_gh_prediction"]),mp.mpf(geo["log2_jet_prediction"])
        prediction=min(gh,jet)
        audited=[]
        candidates=[folder/"full",folder/"full_strong"]+sorted(folder.glob("full_sieve_*"))
        for candidate_folder in candidates:
            result_path=candidate_folder/"minimum.json";audit_path=candidate_folder/"audit.json"
            if not result_path.exists() or not audit_path.exists(): continue
            result,audit=read(result_path),read(audit_path)
            if result.get("exact") and audit.get("passed") and audit.get("minimum_sha256")==engine._hash_file(result_path):
                audited.append((squared_norm(result),result_path,audit_path))
        if audited:
            assert len({value for value,_,_ in audited})==1
            value,source,audit_path=audited[0]
            certified=True;kind="certified_projected_minimum"
        else:
            source=OUT/"sieving"/case/"best_candidate.json"
            point=read(source);audit_path=source.parent/"vector_audit.json";audit=read(audit_path)
            assert audit["passed"] and audit["candidate_sha256"]==engine._hash_file(source)
            assert audit["public_affine_congruences_verified"] and audit["exact_projection_verified"]
            assert audit["exact_squared_length_verified"]
            value=Fraction(int(point["physical_squared_norm_numerator"]),int(point["physical_squared_norm_denominator"]))
            assert value==Fraction(int(audit["physical_squared_norm_numerator"]),int(audit["physical_squared_norm_denominator"]))
            certified=False;kind="checked_sieving_candidate_upper_bound"
        assert value>0
        observed=log_length(value);delta=observed-prediction
        ratio=mp.exp(delta*mp.log(2));gap=100*mp.expm1(delta*mp.log(2))
        record=dict(index=index,case=case,application=cfg[0].upper(),t=cfg[1],m=cfg[2],ell=cfg[3],
            projected_dimension=geo["full_projected_rank"],comparison_kind=kind,shortestness_certified=certified,
            selected_prediction="GH" if gh<=jet else "jet",prediction_log2=str(prediction),
            compared_projected_length_log2=str(observed),log2_difference_bits=str(delta),
            compared_length_over_prediction=str(ratio),excess_over_prediction_percent=str(gap),
            within_5pct_of_prediction=bool(abs(gap)<=5),
            source=str(source.relative_to(HERE)),source_sha256=engine._hash_file(source),
            vector_or_minimum_audit=str(audit_path.relative_to(HERE)),audit_sha256=engine._hash_file(audit_path),
            geometry_sha256=engine._hash_file(geo_path))
        save(DEST/(case+".json"),record);records.append(record)
        print("EMPIRICAL_COMPARISON",case,"prediction",mp.nstr(prediction,16),
              "point",mp.nstr(observed,16),"gap_percent",mp.nstr(gap,10),"certified",certified,flush=True)
    for filename,readable in (("comparison.csv",False),("comparison_readable.csv",True)):
        target=DEST/filename;temporary=target.with_name(target.name+".tmp-"+str(os.getpid()))
        with temporary.open("w",encoding="utf-8-sig",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader()
            for record in records:
                display=dict(record)
                if readable:
                    for key in ("prediction_log2","compared_projected_length_log2","log2_difference_bits"):
                        display[key]=format(float(record[key]),".6f")
                    for key in ("compared_length_over_prediction","excess_over_prediction_percent"):
                        display[key]=format(float(record[key]),".6f")
                writer.writerow(display)
        temporary.replace(target)
    quick=[record for record in records if not record["shortestness_certified"]]
    save(DEST/"validation.json",dict(created_kst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
        cases=len(records),certified_minima=len(records)-len(quick),uncertified_candidates=len(quick),
        all_sources_bound_to_passing_audits=True,
        compared_value="lambda1 of the root projection when certified, otherwise a checked sieve-vector length",
        excess_definition="100*(compared_length/prediction-1), with prediction=min(GH,jet)",
        sieve_candidates_all_within_5pct=all(record["within_5pct_of_prediction"] for record in quick),
        total_within_5pct=sum(record["within_5pct_of_prediction"] for record in records),
        scope="21 seed-0 base B instances; agreement of uncertified candidates is empirical, not a certified error bound"))


if __name__=="__main__": main()
