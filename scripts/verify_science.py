"""Reconstruct exact interval algebra and selected optimized-basis proofs.

This validates archived claims and witnesses. Optional block replays rerun
bounded small SVPs; this is not a fresh full-dimensional SVP certification.
"""
from pathlib import Path
import argparse,json,os,subprocess,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from run_experiment import prepare

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--replay-seconds',type=float,default=0)
    p.add_argument('--indices',default='all')
    args=p.parse_args();work=prepare(args.work)
    env=dict(os.environ,DUAL_COPPERSMITH_ROOT=str(ROOT),PYTHONDONTWRITEBYTECODE='1',
        PYTHONPATH=os.pathsep.join((str(ROOT),str(ROOT/'second_minimum/experiments'),str(work))),
        OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    indices=list(range(21)) if args.indices=='all' else list(map(int,args.indices.split(',')))
    subprocess.run([sys.executable,str(work/'validate_lambda2_intervals.py'),'--indices',','.join(map(str,indices)),
        '--replay-seconds',str(args.replay_seconds)],cwd=work,env=env,check=True)
    for index in (1,2,3,8,9):
        if index not in indices:continue
        subprocess.run([sys.executable,str(work/'validate_bound_bases.py'),'--index',str(index),
            '--replay-seconds',str(args.replay_seconds)],cwd=work,env=env,check=True)
    record=dict(indices=indices,exact_interval_algebra_verified=True,
        all_tested_basis_transforms_verified=True,
        decisive_block_replay_requested=args.replay_seconds>0,
        full_dimensional_exhaustive_search_rerun=False)
    (work/'scientific_validation.json').write_text(json.dumps(record,indent=2)+'\n')
    print('SCIENTIFIC_VALIDATION_OK',json.dumps(record),flush=True)

if __name__=='__main__':main()
