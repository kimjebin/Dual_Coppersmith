"""Collect exactly checked quick-SVP witnesses without promoting them to minima."""
import csv
import os
from datetime import datetime, timedelta, timezone
from fractions import Fraction
import mpmath as mp

from projected_minima import HERE, OUT, GRID, name, read, save
from shared import exact_minima as engine

INDICES = (0, 1, 2, 3, 8, 9)
mp.mp.dps = 80


def log_norm(squared):
    return (mp.log(squared.numerator, 2) - mp.log(squared.denominator, 2)) / 2


def collect():
    destination = OUT / "sieving"
    records = []
    for index in INDICES:
        case = name(GRID[index])
        folder = OUT / case / "full_strong"
        stage = folder / "certificate" / "pass1"
        raw = engine._read_document(stage / "prepared.json")
        data = engine._decode_preparation(raw)
        embedding = read(folder / "embedding.json")
        geometry = read(OUT / case / "geometry.json")
        scale = Fraction(int(embedding["scale_numerator"]), int(embedding["scale_denominator"]))
        initial = int(raw["initial_best"])
        best, best_y = initial, list(map(int, raw["initial_y"]))
        best_origin = "initial_BKZ_basis_candidate"
        best_trial = None
        previous_path=destination/case/'best_candidate.json'
        if previous_path.exists():
            previous=read(previous_path)
            if previous.get('prepared_sha256')==raw['payload_sha256']:
                previous_y=list(map(int,previous['quotient_coefficients']))
                previous_value,_=engine._witness(data,previous_y)
                if previous_value<best:
                    best,best_y,best_origin=previous_value,previous_y,previous['best_origin']
        result_paths = [folder / "certificate" / "result.json"]
        result_paths.extend((OUT / case).glob("full_sieve_*/certificate/result.json"))
        for result_path in result_paths:
            if not result_path.exists():
                continue
            result = read(result_path)
            y = list(map(int, result["passes"][0]["quotient_coefficients"]))
            value, _ = engine._witness(data, y)
            assert value == int(result["passes"][0]["off_line_sq"])
            if value < best:
                best, best_y, best_origin = value, y, str(result_path.relative_to(HERE))
        trials = []
        seen = set()
        for path in sorted(stage.glob("sieve*.json")):
            trial = read(path)
            if "quotient_coefficients" not in trial or trial.get("prepared_sha256") != raw["payload_sha256"]:
                continue
            y = list(map(int, trial["quotient_coefficients"]))
            value, witness = engine._witness(data, y)
            assert value == int(trial["squared_length"]) and any(y) and witness[0] == 0
            key = (trial.get("algorithm", "bgj1"), trial["seed"], trial["context"], str(value),trial.get('source_identity'))
            if key in seen:
                continue
            seen.add(key)
            completed = trial.get("completed_context", max([0] + [int(line.split()[1])
                for line in trial.get("completed_sieve_rounds", []) if line.startswith("HEURISTIC_ROUND_COMPLETE ")]))
            trial_info = dict(path=str(path.relative_to(HERE)), algorithm=trial.get("algorithm", "bgj1"),
                seed=trial["seed"], requested_context=trial["context"], completed_context=completed,
                threads=trial.get('threads',1),log2_physical_length=str(log_norm(Fraction(value)*scale**2)),
                kernel_seconds=trial["seconds"], timed_out=trial["time_limit_reached"],
                source_identity=trial.get('source_identity'),
                kernel_returncode=trial.get("kernel_returncode", 0), integer_squared_norm=str(value))
            trials.append(trial_info)
            if value < best:
                best, best_y = value, y
                best_origin = trial_info["path"]
                best_trial = trial_info
        norm, witness = engine._witness(data, best_y)
        assert norm == best and norm > 0 and witness[0] == 0
        numerator = [z * scale.numerator for z in witness[1:]]
        assert sum(z * int(r) for z, r in zip(numerator, geometry["root"])) == 0
        physical = Fraction(best) * scale ** 2
        assert sum(z*z for z in numerator) * physical.denominator == physical.numerator * scale.denominator ** 2
        initial_physical = Fraction(initial) * scale ** 2
        initial_log, best_log = log_norm(initial_physical), log_norm(physical)
        agreeing = sorted({trial["seed"] for trial in trials if int(trial["integer_squared_norm"]) == best})
        certified=False
        for variant_folder in [folder]+list((OUT/case).glob("full_sieve_*")):
            minimum_path=variant_folder/'minimum.json';audit_path=variant_folder/'audit.json'
            if minimum_path.exists() and audit_path.exists():
                minimum=read(minimum_path);audit=read(audit_path)
                if minimum.get('exact') and audit.get('passed') and audit.get('minimum_sha256')==engine._hash_file(minimum_path):
                    exact_sq=Fraction(int(minimum['squared_norm_numerator']),int(minimum['squared_norm_denominator']))
                    assert physical>=exact_sq
                    certified=certified or physical==exact_sq
        record = dict(index=index, case=case, dimension=len(data["tail"]), certified_minimum=certified,
            initial_log2_length=str(initial_log), best_log2_length=str(best_log),
            improvement_bits=str(initial_log-best_log), best_origin=best_origin,
            seed_count_at_best=len(agreeing), seeds_at_best=agreeing, trials=trials,
            best_trial=best_trial, quotient_coefficients=list(map(str,best_y)),
            integer_embedding_witness=list(map(str,witness)),
            scale_numerator=str(scale.numerator), scale_denominator=str(scale.denominator),
            physical_squared_norm_numerator=str(physical.numerator), physical_squared_norm_denominator=str(physical.denominator),
            source_preparation=str(stage.relative_to(HERE)), prepared_sha256=raw["payload_sha256"],
            embedding_sha256=engine._hash_file(folder / "embedding.json"),
            status="Certified projected minimum" if certified else "Exactly checked nonzero projected-lattice vector; shortestness not certified")
        save(destination / case / "best_candidate.json", record)
        records.append(record)
    destination.mkdir(parents=True, exist_ok=True)
    csv_temporary=destination / ("candidate_comparison.csv.tmp-"+str(os.getpid()))
    with csv_temporary.open("w", encoding="utf-8-sig", newline="") as handle:
        fieldnames = ["index", "case", "dimension", "initial_log2_length", "best_log2_length", "improvement_bits",
                      "seed_count_at_best", "best_origin", "certified_minimum"]
        writer=csv.DictWriter(handle,fieldnames=fieldnames,extrasaction="ignore")
        writer.writeheader();writer.writerows(records)
    csv_temporary.replace(destination / "candidate_comparison.csv")
    readable_temporary=destination / ("candidate_comparison_readable.csv.tmp-"+str(os.getpid()))
    with readable_temporary.open("w",encoding="utf-8-sig",newline="") as handle:
        fieldnames=["case","dimension","initial_log2_length","best_log2_length","improvement_bits","seed_count_at_best","certified_minimum"]
        writer=csv.DictWriter(handle,fieldnames=fieldnames,extrasaction="ignore");writer.writeheader()
        for record in records:
            display=dict(record)
            for key in ("initial_log2_length","best_log2_length","improvement_bits"):
                display[key]=format(float(record[key]),".9f")
            writer.writerow(display)
    readable_temporary.replace(destination / "candidate_comparison_readable.csv")
    methods_temporary=destination/('method_trials.csv.tmp-'+str(os.getpid()))
    with methods_temporary.open('w',encoding='utf-8-sig',newline='') as handle:
        fields=['case','dimension','algorithm','seed','threads','requested_context','completed_context','kernel_seconds',
                'timed_out','kernel_returncode','log2_physical_length','certified_minimum','source_identity','path']
        writer=csv.DictWriter(handle,fieldnames=fields,extrasaction='ignore');writer.writeheader()
        for record in records:
            for trial in record['trials']:
                writer.writerow(dict(trial,case=record['case'],dimension=record['dimension'],certified_minimum=False))
    methods_temporary.replace(destination/'method_trials.csv')
    save(destination / "status.json", dict(created_kst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
         cases=len(records), improved=sum(mp.mpf(record["improvement_bits"])>0 for record in records),
         certified_minima=sum(record['certified_minimum'] for record in records), independent_seed_agreement_is_not_a_shortestness_proof=True,
         records=[{key:record[key] for key in ("index","case","dimension","initial_log2_length","best_log2_length","improvement_bits","seed_count_at_best")}
                  for record in records]))
    for record in records:
        print("QUICK_SVP_CANDIDATE", record["index"], record["dimension"],
              mp.nstr(mp.mpf(record["initial_log2_length"]),16), "->", mp.nstr(mp.mpf(record["best_log2_length"]),16),
              "gain_bits",mp.nstr(mp.mpf(record["improvement_bits"]),8),"seeds",record["seeds_at_best"], flush=True)
    return records


if __name__ == "__main__":
    collect()
