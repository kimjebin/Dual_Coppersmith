"""Run meaningful small scientific regressions in an isolated directory."""
from pathlib import Path
import argparse,os,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_experiment import prepare

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work',type=Path,required=True)
    args=parser.parse_args();work=prepare(args.work,hydrate=False)
    env=dict(os.environ,DUAL_COPPERSMITH_ROOT=str(ROOT),PYTHONDONTWRITEBYTECODE='1',
        PYTHONPATH=os.pathsep.join((str(ROOT),str(ROOT/'second_minimum/experiments'),str(work))),
        OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    for filename in ('test_core.py','test_exact_minima.py','test_results.py'):
        subprocess.run([sys.executable,str(ROOT/'tests'/filename)],cwd=work,env=env,check=True)
    for filename in ('test_coefficient_certificate.py','validate_gram_checker.py'):
        subprocess.run([sys.executable,str(work/filename)],cwd=work,env=env,check=True)
    sage=Path(sys.executable).parent/'sage'
    hm=work/'hm/hm_dual.sage';hm.parent.mkdir(exist_ok=True)
    shutil.copyfile(ROOT/'hm/hm_dual.sage',hm)
    subprocess.run([str(sage),str(hm),'80','20','2','4','2','dual'],cwd=work,env=env,check=True)
    print('SCIENTIFIC_SMOKE_OK',flush=True)

if __name__=='__main__':main()
