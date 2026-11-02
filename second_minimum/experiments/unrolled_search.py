"""Research compilation only; not wired into accepted certification."""
from shared_bootstrap import *
import hashlib,json,subprocess,sys,time
from pathlib import Path
from sage.all import ZZ,matrix
from experiment import HERE,dump
from fast_search import compile_once,verify_lll_precondition
from shared.exact_minima import _hash_file,_witness
def build(d):
    if not (HERE/'generated/fixed_enumlib.h').exists():
        raise FileNotFoundError('Restore the shipped frozen generated headers before compiling')
    paths=[HERE/'unrolled_subtree.cpp',HERE/'fast_lift.cpp',HERE/'generated'/'fixed_enumlib.h',
           HERE/'generated'/'fixed_fplll_types.h']
    frozen=HERE/'unrolled_validated_sources.json'
    if frozen.exists():
        for name,checksum in json.loads(frozen.read_text())['files'].items():
            assert _hash_file(HERE/name)==checksum,'Validated proof source changed: '+name
    identity=hashlib.sha256((str(d)+''.join(_hash_file(path) for path in paths)).encode()).hexdigest()
    binary=HERE/'bin'/('unrolled_subtree-'+identity[:16]);prefix=sys.prefix
    compile_once(binary,['/usr/bin/g++','-O3','-march=native','-std=c++17','-DPROOF_DIM='+str(d),
        str(HERE/'unrolled_subtree.cpp'),'-I'+prefix+'/include','-L'+prefix+'/lib',
        '-Wl,-rpath,'+prefix+'/lib','-lfplll','-lmpfr','-lgmpxx','-lgmp','-pthread'])
    return binary,identity

def search(data,folder,native_binary,**options):
    began=time.monotonic();folder=Path(folder)
    preconditions=verify_lll_precondition(data,folder)
    T=matrix(ZZ,data['tail']);inner=matrix(ZZ,[data['inner']]);a=ZZ(data['root_sq']);c=ZZ(data['content'])
    H=(a*T*T.transpose()-inner.transpose()*inner).apply_map(lambda z:z//c)
    best=int(data['raw']['initial_best']);y=list(map(int,data['raw']['initial_y']))
    path=folder/'unrolled.input'
    path.write_text(f'FAST-LIFT-1 {T.nrows()} {a} {c} {best}\n'+' '.join(map(str,data['inner']))+'\n'+
        ' '.join(map(str,y))+'\n'+''.join(' '.join(map(str,row))+'\n' for row in H.rows()))
    binary,identity=build(T.nrows())
    from partition_search import run_job
    completed=run_job(binary,path,[],options.get('max_seconds') or None)
    if completed: best=int(completed['best']);y=completed['y']
    sq,witness=_witness(data,y);assert sq==best
    result=dict(complete=completed is not None,reason='complete' if completed else 'time_limit',
        off_line_sq=str(best),off_line_witness=list(map(str,witness)),quotient_coefficients=list(map(str,y)),
        frontier_jobs=0 if completed else None,nodes=completed['nodes'] if completed else None,
        leaves=completed['leaves'] if completed else None,search_seconds=time.monotonic()-began,
        resource_plan=dict(workers=1),prepared_sha256=data['raw']['payload_sha256'],
        unrolled_source_identity=identity,unrolled_binary_sha256=_hash_file(binary),preconditions=preconditions,
        backend='RESEARCH: recursive fixed-suffix enumlib and exact lifts')
    dump(folder/'unrolled_search_result.json',result)
    return result
