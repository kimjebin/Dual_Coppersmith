"""Portable entry points for original-case experiments in isolated workspaces."""
from pathlib import Path
import argparse,json,os,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"scripts"))
from unpack_data import unpack

TASKS={
 "geometry":("projected_minima.py",["prepare"]),
 "projected":("projected_minima.py",["run"]),
 "original-minima":("certify_partition.py",[]),
 "rhf":("measure_root_hermite.py",[]),
 "optimize":("optimize_bound_basis.py",[]),
 "certify-bounds":("certify_bound_bases.py",[]),
 "verify-bounds":("validate_bound_bases.py",[]),
 "verify-vectors":("validate_quick_projected_vectors.py",[]),
 "sieve":("probe_projected_sieve.py",[]),
}

def prepare(work,case=None,hydrate=True,new_search=False):
    work=Path(work).resolve();work.mkdir(parents=True,exist_ok=True)
    for source in (ROOT/"experiments/engine").rglob("*"):
        if not source.is_file() or "__pycache__" in source.parts:continue
        relative=source.relative_to(ROOT/"experiments/engine");target=work/relative
        if target.exists() and target.read_bytes()!=source.read_bytes():raise FileExistsError("Use a new work directory; source files differ")
        target.parent.mkdir(parents=True,exist_ok=True)
        if not target.exists():shutil.copyfile(source,target)
    (work/"logs").mkdir(exist_ok=True)
    def skip(name):
        if not new_search:return False
        if name=="results/"+case+"/B.certificate.json":return True
        if name.startswith("certificates/"+case+"/"):return True
        if name.startswith("projected_minima/"+case+"/"):
            return "/certificate/" in name or name.endswith(("/minimum.json","/audit.json"))
        return False
    if hydrate:unpack(work,case,skip=skip)
    return work


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task",choices=TASKS);parser.add_argument("--case");parser.add_argument("--work",type=Path,default=ROOT/".work")
    parser.add_argument("--fresh-geometry",action="store_true",help="Recompute geometry from the seeded original input without hydrating cached geometry")
    parser.add_argument("--new-search",action="store_true",help="Do not hydrate archived minima or search states for the selected case")
    parser.add_argument("--no-unpack",action="store_true",help="Continue in an already initialized work directory without reference extraction")
    args,extra=parser.parse_known_args();extra=extra[1:] if extra[:1]==["--"] else extra
    configs=json.loads((ROOT/"configs/cases_21.json").read_text());index=None
    if args.case:
        matches=[case for case in configs if case["case_id"]==args.case]
        if not matches:parser.error("Unknown original case identifier")
        index=matches[0]["index"]
    if args.task in {"original-minima","optimize","certify-bounds","verify-bounds"} and index is None:parser.error("This task requires --case")
    if args.task in {"rhf","verify-vectors"} and args.case:parser.error("This task checks all recorded cases; omit --case")
    if args.new_search and args.no_unpack:parser.error("A new search requires reference extraction; omit --no-unpack")
    if args.new_search and (not args.case or args.task not in {"original-minima","projected"}):parser.error("--new-search requires a case and a projected/original-minima task")
    work=prepare(args.work,args.case,not (args.fresh_geometry or args.no_unpack),args.new_search)
    filename,prefix=TASKS[args.task];command=[sys.executable,str(work/filename)]+prefix
    if index is not None:
        if args.task in {"geometry","projected"}:command += ["--indices",str(index)]
        elif args.task not in {"rhf","sieve","verify-vectors"}:command += ["--index",str(index)]
    env=dict(os.environ,DUAL_COPPERSMITH_ROOT=str(ROOT),PYTHONPATH=str(ROOT)+os.pathsep+str(ROOT/"baseline")+os.pathsep+os.environ.get("PYTHONPATH",""),
             PYTHONDONTWRITEBYTECODE="1",OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    subprocess.run(command+extra,cwd=work,env=env,check=True)


if __name__=="__main__":main()
