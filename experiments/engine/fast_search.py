"""fplll's proved error-bound enumeration, with exact original-lattice lifts."""
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
from sage.all import ZZ, matrix
from experiment import HERE, dump
from shared.exact_minima import _witness,_hash_file


def compile_once(binary,command):
    import fcntl
    binary.parent.mkdir(exist_ok=True)
    with binary.with_suffix('.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if binary.exists(): return
        temporary=binary.with_name(binary.name+'.tmp-'+str(os.getpid()))
        subprocess.run(command+['-o',str(temporary)],check=True)
        temporary.replace(binary)


def verify_lll_precondition_reference(data,folder):
    path=Path(folder)/'integer.prepared'
    digest=_hash_file(path)
    assert digest==data['raw']['preparation_sha256']
    with path.open() as f:
        magic,d,a,c=f.readline().split()
        d=int(d)
        assert magic=='EXACT-GMP-PREP-1' and d==len(data['tail'])
        assert int(a)==data['root_sq'] and int(c)==data['content']
        assert list(map(int,f.readline().split()))==data['inner']
        delta=list(map(int,f.readline().split()))
        assert len(delta)==d+1 and delta[0]==1 and all(z>0 for z in delta)
        mu_upper=[[] for _ in range(d)]
        denominator=1<<64
        for i in range(1,d):
            lower=list(map(int,f.readline().split()))
            assert len(lower)==i
            assert all(100*abs(lower[j])<=51*delta[j+1] for j in range(i)), ('size',i)
            assert 100*(delta[i+1]*delta[i-1]+lower[-1]**2)>=99*delta[i]**2,('Lovasz',i)
            mu_upper[i]=[(abs(lower[j])*denominator+delta[j+1]-1)//delta[j+1] for j in range(i)]
        radius=data['root_sq']*(int(data['raw']['initial_best'])-1)//data['content']
        bounds=[0]*d
        for i in reversed(range(d)):
            radicand=(radius*delta[i]+delta[i+1]-1)//delta[i+1]
            term=math.isqrt(radicand)
            term+=term*term<radicand
            weighted=sum(mu_upper[j][i]*bounds[j] for j in range(i+1,d))
            bounds[i]=term+(weighted+denominator-1)//denominator
        assert max(bounds)<2**50, 'Quotient coefficients need more than the supported exact float integer range'
    return dict(exact_lll_preconditions=True,delta='99/100',eta='51/100',
                integer_preparation_sha256=digest,
                exact_integer_coefficient_bound=str(max(bounds)),coefficient_bounds_verified=True)


def verify_lll_precondition(data,folder):
    path=Path(folder)/'integer.prepared'
    digest=_hash_file(path)
    assert digest==data['raw']['preparation_sha256']
    with path.open() as f:
        header=f.readline().split()
        assert header==['EXACT-GMP-PREP-1',str(len(data['tail'])),str(data['root_sq']),str(data['content'])]
        assert f.readline().split()==list(map(str,data['inner']))
    source=HERE/'verify_lll.cpp'
    source_hash=_hash_file(source)
    binary=HERE/'bin'/('verify_lll-'+source_hash[:16])
    if not binary.exists():
        compile_once(binary,['g++','-O3','-std=c++17',str(source),'-lgmpxx','-lgmp'])
    run=subprocess.run([str(binary),str(path),str(data['raw']['initial_best'])],
                       text=True,capture_output=True,check=True)
    tag,dimension,bound=run.stdout.split()
    assert tag=='LLL_OK' and int(dimension)==len(data['tail']) and int(bound)<2**50
    return dict(exact_lll_preconditions=True,delta='99/100',eta='51/100',
        integer_preparation_sha256=digest,exact_integer_coefficient_bound=bound,
        coefficient_bounds_verified=True,verifier_source_sha256=source_hash,
        verifier_binary_sha256=_hash_file(binary))


def build_binary():
    source = HERE/'fast_lift.cpp'
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    binary = HERE/'bin'/('fast_lift-'+digest[:16])
    if not binary.exists():
        prefix = Path(sys.prefix)
        command = ['g++','-O3','-std=c++17',str(source),'-I'+str(prefix/'include'),
                   '-L'+str(prefix/'lib'),'-Wl,-rpath,'+str(prefix/'lib'),
                   '-lfplll','-lmpfr','-lgmpxx','-lgmp']
        compile_once(binary,command)
    return binary,digest


def search(data,folder,native_binary,**options):
    began=time.monotonic()
    preconditions=verify_lll_precondition(data,folder)
    if int(preconditions['exact_integer_coefficient_bound'])>=2**29:
        raise ValueError('Whole-sphere external enumlib integer range is too large; use certify_partition.py')
    binary,source_hash=build_binary()
    saved_path=Path(folder)/'fast_search_result.json'
    if saved_path.exists():
        saved=json.loads(saved_path.read_text())
        if (saved.get('complete') and saved.get('fast_source_sha256')==source_hash and
            saved.get('fast_binary_sha256')==_hash_file(binary) and
            saved.get('prepared_sha256')==data['raw']['payload_sha256']):
            y=list(map(int,saved['quotient_coefficients']))
            assert _witness(data,y)[0]==int(saved['off_line_sq'])
            saved.update(certificate_reused=True,invocation_search_seconds=time.monotonic()-began,
                         invocation_nodes=0,preconditions=preconditions)
            print('FAST_CERTIFICATE_REUSED',str(folder),flush=True)
            return saved
    a,c=ZZ(data['root_sq']),ZZ(data['content'])
    T=matrix(ZZ,data['tail'])
    inner=matrix(ZZ,[data['inner']])
    gram=a*(T*T.transpose())-inner.transpose()*inner
    assert all(z%c==0 for z in gram.list())
    gram=gram.apply_map(lambda z:z//c)
    best=int(data['raw']['initial_best'])
    y=list(map(int,data['raw']['initial_y']))
    assert _witness(data,y)[0]==best
    path=Path(folder)/'fast.input'
    with path.open('w') as f:
        f.write(f'FAST-LIFT-1 {T.nrows()} {a} {c} {best}\n')
        f.write(' '.join(map(str,data['inner']))+'\n')
        f.write(' '.join(map(str,y))+'\n')
        for row in gram.rows(): f.write(' '.join(map(str,row))+'\n')
    limit=options.get('max_seconds',0)
    completed=False
    try:
        run=subprocess.run([str(binary),str(path)],text=True,capture_output=True,
                           timeout=limit or None,check=True)
        lines=run.stdout.splitlines()
        tag,value,nodes,leaves,seconds,precision=lines[0].split()
        assert tag=='COMPLETE'
        best=int(value); y=list(map(int,lines[1].split()))
        completed=True
        details=dict(nodes=int(nodes),leaves=int(leaves),native_seconds=float(seconds),precision=int(precision))
    except subprocess.TimeoutExpired:
        details=dict(nodes=None,leaves=None,native_seconds=None,precision=None)
    sq,witness=_witness(data,y)
    assert sq==best
    result=dict(complete=completed,reason='complete' if completed else 'time_limit',
        off_line_sq=str(best),off_line_witness=list(map(str,witness)),
        quotient_coefficients=list(map(str,y)),frontier_jobs=0 if completed else None,
        search_seconds=time.monotonic()-began,invocation_search_seconds=time.monotonic()-began,
        invocation_nodes=details['nodes'],resource_plan=dict(workers=1),
        prepared_sha256=data['raw']['payload_sha256'],fast_source_sha256=source_hash,
        fast_binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
        backend='fplll 5.5.0 error-bounded complete enumeration and exact integer lifts',**details)
    result['preconditions']=preconditions
    dump(Path(folder)/'fast_search_result.json',result)
    print('FAST_SEARCH',Path(folder),completed,details,flush=True)
    return result


if __name__=='__main__':
    import argparse
    from shared.exact_minima import _read_document,_decode_preparation
    p=argparse.ArgumentParser()
    p.add_argument('preparation',type=Path)
    p.add_argument('--seconds',type=float,default=600)
    args=p.parse_args()
    data=_decode_preparation(_read_document(args.preparation/'prepared.json'))
    search(data,args.preparation,None,max_seconds=args.seconds)
