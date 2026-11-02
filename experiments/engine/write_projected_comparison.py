"""Publish explicit exact/interval/unresolved projected-minimum comparisons."""
import csv
from datetime import datetime, timedelta, timezone
from decimal import Decimal

from projected_minima import OUT, GRID, name, read, save, summarize
from shared import exact_minima as engine


def rounded(value, places=6):
    return format(Decimal(str(value)), f".{places}f")


def main():
    records = summarize()
    readable = []
    jet_audited = full_audited = interval_only = 0
    for record in records:
        folder = OUT / record["case"]
        jet = read(folder / "jet" / "minimum.json")
        jet_audit = read(folder / "jet" / "audit.json")
        assert jet["exact"] and jet_audit["passed"]
        assert jet_audit["minimum_sha256"] == engine._hash_file(folder / "jet" / "minimum.json")
        jet_audited += 1
        if record["full_exact"]:
            from pathlib import Path
            from experiment import HERE
            minimum_path = HERE / record["full_result_source"]
            audit = read(minimum_path.parent / "audit.json")
            assert audit["passed"] and audit["minimum_sha256"] == engine._hash_file(minimum_path)
            full_audited += 1
            full_value = rounded(record["full_log2_minimum"])
            full_status = "exact_squared_norm; independently audited"
            assert (int(record["full_squared_norm_numerator"]) * int(jet["squared_norm_denominator"])
                    <= int(jet["squared_norm_numerator"]) * int(record["full_squared_norm_denominator"]))
        elif record["full_value_certified_6dp"]:
            from experiment import HERE
            interval = read(folder / "full_bounds.json")
            assert engine._hash_file(HERE / interval["source_certificate"]) == interval["source_certificate_sha256"]
            assert engine._hash_file(HERE / "audits" / (record["case"] + ".json")) == interval["source_independent_audit_sha256"]
            full_value = record["full_log2_display_6dp"]
            full_status = "rigorous interval certifies six decimals; exact squared norm unfinished"
            interval_only += 1
        else:
            full_value = ""
            full_status = "minimum not certified; enumeration ongoing"
        candidate_path=OUT/'sieving'/record['case']/'best_candidate.json'
        quick_candidate=read(candidate_path) if candidate_path.exists() else {}
        upper_values=[value for value in (record['full_log2_upper_bound'],quick_candidate.get('best_log2_length','')) if value]
        best_upper=min(upper_values,key=Decimal) if upper_values else ''
        readable.append(dict(
            application=record["application"], t=record["t"], m=record["m"], ell=record["ell"],
            dimension=record["ambient_dimension"],
            predicted_GH_log2=rounded(record["log2_gh_prediction"]), actual_full_projected_log2=full_value,
            full_status=full_status, full_candidate_upper_bound_log2=rounded(best_upper) if best_upper else "",
            quick_svp_candidate_log2=rounded(quick_candidate['best_log2_length']) if quick_candidate else '',
            quick_svp_shortestness_certified=quick_candidate.get('certified_minimum',''),
            predicted_jet_log2=rounded(record["log2_jet_prediction"]), actual_jet_projected_log2=rounded(jet["log2_length"]),
            jet_status="exact_squared_norm; independently audited"))
    with (OUT / "comparison_readable.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(readable[0]))
        writer.writeheader()
        writer.writerows(readable)
    assert len(readable) == 21 and jet_audited == 21
    save(OUT / "comparison_validation.json", dict(
        created_kst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
        inputs=21, jet_exact_independently_audited=jet_audited,
        full_exact_independently_audited=full_audited, full_interval_certified_6dp=interval_only,
        full_values_certified_6dp=full_audited + interval_only,
        unresolved_full=21 - full_audited - interval_only,
        units="base-2 logarithm of Euclidean length",
        definitions=dict(actual_full="lambda_1(pi(L))", actual_jet="lambda_1(pi(L_jet))",
            pi="orthogonal projection perpendicular to the primitive planted-root vector",
            predicted_GH="Paper table convention: exact unit-ball volume, rank n-1",
            predicted_jet="Paper convention: sqrt(ell/(2*pi*e)) times projected jet covolume^(1/ell)"),
        exact_data="comparison.csv and case/kind/minimum.json include exact rational squared norms",
        completeness="Incomplete until all 21 full projected minima are certified; no candidate or heuristic is labelled as an actual minimum."))
    print("VALIDATED_COUNTS", jet_audited, full_audited, interval_only, "unresolved", 21 - full_audited - interval_only)
    for row in readable:
        print(f"{row['application']} ({row['t']},{row['m']},{row['ell']}) {row['dimension']} "
              f"GH={row['predicted_GH_log2']} FULL={row['actual_full_projected_log2'] or 'UNRESOLVED'} "
              f"JET_PRED={row['predicted_jet_log2']} JET_ACTUAL={row['actual_jet_projected_log2']}")


if __name__ == "__main__":
    main()
