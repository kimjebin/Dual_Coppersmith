"""Regenerate public comparison tables from the canonical original-case records."""
from __future__ import annotations
import csv,json
from decimal import Decimal,localcontext
from pathlib import Path
from paths import ROOT


def log_length(squared):
    with localcontext() as context:
        context.prec=80
        value=Decimal(squared["numerator"])/Decimal(squared.get("denominator","1"))
        return str(value.ln()/Decimal(2).ln()/2)


def case_records():
    return json.loads((ROOT/"results/minima.json").read_text())


def regenerate(destination=None):
    destination=Path(destination) if destination else ROOT/"results"
    destination.mkdir(parents=True,exist_ok=True);records=case_records();rows=[]
    for case in records:
        original=case["original"];upper=original["lambda2"]["upper_candidates"]
        by_kind={item["origin"]:item for item in upper}
        local=[item for item in upper if item.get("squared_norm")]
        best=min(local,key=lambda item:Decimal(item["log2_length"]))
        projected=case["projection"]["full"]
        rows.append(dict(case=case["case_id"],application=case["application"],t=case["t"],m=case["m"],ell=case["ell"],dimension=case["dimension"],
            original_lambda1_log2=original["lambda1"]["log2_length"],lambda2_certified_lower_log2=original["lambda2"]["lower"]["log2_length"],
            prediction_GH_log2=case["predictions"]["GH_log2"],prediction_jet_log2=case["predictions"]["jet_log2"],
            projected_value_log2=projected["log2_length"],projected_value_kind=projected["value_kind"],
            original_lambda2_exact_log2=original["lambda2"]["exact"]["log2_length"] if original["lambda2"]["exact"] else "",
            original_paper_lambda2_log2=by_kind["original_paper_reported"]["log2_length"],
            revised_pdf_lambda2_upper_log2=by_kind.get("revised_report_candidate",{}).get("log2_length",""),
            local_BKZ_lambda2_upper_log2=by_kind.get("original_lattice_BKZ_candidate",{}).get("log2_length",""),
            lifted_projection_lambda2_upper_log2=by_kind.get("lifted_projection_candidate",{}).get("log2_length",""),
            best_local_lambda2_upper_log2=best["log2_length"],best_local_upper_origin=best["origin"],
            jet_projected_exact_log2=case["projection"]["jet"]["log2_length"],
            original_lambda2_status="archived_audited_exact" if original["lambda2"]["exact"] else "certified_interval; exact_search_unfinished"))
    assert len(rows)==21
    for filename,subset in (("all_21_comparison.csv",rows),("unfinished_5_comparison.csv",[row for row in rows if not row["original_lambda2_exact_log2"]])):
        with (destination/filename).open("w",newline="",encoding="utf-8") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(subset)
    columns=[("case","Case"),("original_lambda1_log2","Original lambda1"),("lambda2_certified_lower_log2","Certified lower"),
             ("prediction_GH_log2","GH prediction"),("prediction_jet_log2","Jet prediction"),
             ("projected_value_log2","Projected value"),("projected_value_kind","Projected status"),("revised_pdf_lambda2_upper_log2","Revised PDF upper")]
    lines=["# Unfinished original-case comparison","","All numerical lengths below are base-2 logarithms. PDF values retain their published three-decimal precision.","",
           "| "+" | ".join(label for _,label in columns)+" |","|"+"---|"*len(columns)]
    for row in rows:
        if row["original_lambda2_exact_log2"]:continue
        values=[]
        for key,_ in columns:
            value=str(row[key])
            if key.endswith("log2") and value:
                value=format(Decimal(value),".3f" if key=="revised_pdf_lambda2_upper_log2" else ".6f")
            values.append(value)
        lines.append("| "+" | ".join(values)+" |")
    lines.extend(["","Projected candidate lengths are upper bounds on the projected first minimum; they are not certified shortest lengths.",
        "The CSV also separates the original paper, revised PDF, original-lattice BKZ, and lifting-derived upper candidates.",
        "A lifting-derived upper candidate is not an independent validation of the prediction against the projected vector used to construct it.",""])
    (destination/"unfinished_5_comparison.md").write_text("\n".join(lines),encoding="utf-8")
    status=dict(original_cases=21,original_lambda1_certified=21,original_lambda2_archived_exact=sum(bool(row["original_lambda2_exact_log2"]) for row in rows),
        original_lambda2_unfinished=sum(not bool(row["original_lambda2_exact_log2"]) for row in rows),
        projected_full_archived_exact=sum(case["projection"]["full"]["value_kind"]=="archived_audited_exact" for case in records),
        projected_jet_archived_exact=sum(case["projection"]["jet"]["value_kind"]=="archived_audited_exact" for case in records),
        global_basis_optimality_proven=False,full_searches_rerun_for_packaging=False)
    return status
