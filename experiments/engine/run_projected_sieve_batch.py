"""Finite, two-process alternative-SVP comparison on unresolved projections."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
import subprocess
import sys
import time
import argparse

from projected_minima import HERE, OUT, GRID, name, read, save
from collect_projected_candidates import collect

PROFILES = [(1,"bdgl",11,77,120), (8,"bdgl",11,83,180),
            (2,"bdgl",11,90,180), (3,"bdgl",11,104,180),
            (9,"bgj1",11,119,180), (0,"hk3",12,65,60)]


def execute(profile):
    index, algorithm, seed, context, seconds, *extra = profile
    threads=extra[0] if extra else 1
    stage = OUT / name(GRID[index]) / "full_strong/certificate/pass1"
    output = stage / f"sieve_{algorithm}_seed{seed}_ctx{context}.json"
    if output.exists():
        return dict(index=index, algorithm=algorithm, seed=seed, reused=True, output=str(output.relative_to(HERE)))
    began = time.monotonic()
    command = [sys.executable, str(HERE / "probe_projected_sieve.py"), str(stage), "--context", str(context),
               "--algorithm", algorithm, "--seed", str(seed), "--seconds", str(seconds), "--threads", str(threads)]
    run = subprocess.run(command, text=True, capture_output=True, timeout=seconds+90)
    if run.returncode or not output.exists():
        return dict(index=index, algorithm=algorithm, seed=seed, failed=True, returncode=run.returncode,
                    diagnostics=(run.stderr+run.stdout)[-4000:], elapsed_seconds=time.monotonic()-began)
    result = read(output)
    return dict(index=index, algorithm=algorithm, seed=seed, output=str(output.relative_to(HERE)),
                elapsed_seconds=time.monotonic()-began, log2_length=result["log2_physical_length"],
                completed_context=result["completed_context"], requested_context=result["context"],
                timed_out=result["time_limit_reached"], kernel_returncode=result["kernel_returncode"])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--extended',action='store_true');args=parser.parse_args()
    profiles=[(3,'bgj1',14,104,600,4),(9,'bdgl',14,119,600,4)] if args.extended else PROFILES
    progress_path=OUT/('sieving/followup_progress.json' if args.extended else 'sieving/batch_progress.json')
    started = datetime.now(timezone(timedelta(hours=9))).isoformat()
    completed=[]
    save(progress_path, dict(state="running", started_kst=started, profiles=profiles,
                                                  concurrency=2, completed=completed))
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures={pool.submit(execute,profile):profile for profile in profiles}
        for future in as_completed(futures):
            try:
                result=future.result()
            except Exception as exc:
                result=dict(profile=futures[future], failed=True, error=repr(exc))
            completed.append(result)
            print("SIEVE_PROFILE_RESULT",result,flush=True)
            collect()
            save(progress_path, dict(state="running", started_kst=started, profiles=profiles,
                 concurrency=2, completed=completed))
    save(progress_path, dict(state="complete", started_kst=started,
         finished_kst=datetime.now(timezone(timedelta(hours=9))).isoformat(), concurrency=2,
         profiles=profiles, completed=completed, failed=sum(bool(item.get("failed")) for item in completed)))


if __name__=="__main__":
    main()
