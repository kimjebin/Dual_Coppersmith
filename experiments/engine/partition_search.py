"""Resumable exact partition of an unpruned SVP sphere into fixed suffixes.

The integer box comes from the independently checked exact Bareiss data.
An internal tree node is replaced atomically by ALL its canonical-sign children.
Every leaf must complete fplll's error-bounded, unpruned enumeration. A timeout
only causes a split; it is never evidence for a lower bound.
"""
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import sys
import time
from sage.all import ZZ, matrix
from experiment import HERE,dump
from shared.exact_minima import _witness,_hash_file
from fast_search import compile_once,verify_lll_precondition
from suffix_filter import load_filter,proven_empty,candidate_interval


def identity_hash(paths):
    return hashlib.sha256(''.join(_hash_file(p) for p in paths).encode()).hexdigest()


def binaries():
    prefix=Path(sys.prefix)
    cpp=identity_hash([HERE/'fast_subtree.cpp',HERE/'fast_lift.cpp'])
    binary=HERE/'bin'/('fast_subtree-'+cpp[:16])
    compile_once(binary,['g++','-O3','-std=c++17',str(HERE/'fast_subtree.cpp'),
        '-I'+str(prefix/'include'),'-L'+str(prefix/'lib'),'-Wl,-rpath,'+str(prefix/'lib'),
        '-lfplll','-lmpfr','-lgmpxx','-lgmp'])
    bounds_hash=identity_hash([HERE/'partition_bounds.cpp',HERE/'verify_lll.cpp'])
    bounder=HERE/'bin'/('partition_bounds-'+bounds_hash[:16])
    compile_once(bounder,['g++','-O3','-std=c++17',str(HERE/'partition_bounds.cpp'),'-lgmpxx','-lgmp'])
    return binary,cpp,bounder,bounds_hash


def children(suffix,bounds):
    if len(suffix)>=len(bounds): return []
    limit=bounds[len(bounds)-len(suffix)-1]
    # Highest nonzero coordinate is positive. Subsequent coordinates have both
    # signs, because fixed nonzero suffixes break the SVP +/- symmetry.
    return [[j]+suffix for j in range(-limit if any(suffix) else 0,limit+1)]


def run_job(binary,path,suffix,seconds):
    try:
        result=subprocess.run([str(binary),str(path)]+list(map(str,suffix)),
            capture_output=True,text=True,check=True,timeout=seconds)
    except subprocess.TimeoutExpired:
        return None
    except subprocess.CalledProcessError as error:
        raise RuntimeError(f'Proved subtree enumerator failed: suffix={suffix}; stderr={error.stderr.strip()}') from error
    lines=result.stdout.splitlines()
    tag,best,nodes,leaves,elapsed,precision=lines[0].split()
    assert tag=='COMPLETE'
    return dict(best=best,y=list(map(int,lines[1].split())),nodes=int(nodes),leaves=int(leaves),
        seconds=float(elapsed),precision=int(precision))


def search(data,folder,native_binary,**options):
    began=time.monotonic(); folder=Path(folder)
    preconditions=verify_lll_precondition(data,folder)
    binary,source_hash,bounder,bounder_hash=binaries()
    bound_file=folder/'partition_bounds.json'
    bound_identity=dict(prepared_sha256=data['raw']['payload_sha256'],bounder_source_sha256=bounder_hash,
                        integer_preparation_sha256=data['raw']['preparation_sha256'])
    if bound_file.exists():
        bound_doc=json.loads(bound_file.read_text())
        assert all(bound_doc[k]==v for k,v in bound_identity.items())
        bounds=list(map(int,bound_doc['bounds']))
    else:
        result=subprocess.run([str(bounder),str(folder/'integer.prepared'),str(data['raw']['initial_best'])],
            capture_output=True,text=True,check=True)
        lines=result.stdout.splitlines()
        assert lines[0].startswith('LLL_OK ') and lines[1].startswith('BOUNDS ')
        bounds=list(map(int,lines[1].split()[1:]))
        dump(bound_file,dict(**bound_identity,bounds=list(map(str,bounds)),
                            bounder_binary_sha256=_hash_file(bounder)))
    assert len(bounds)==len(data['tail']) and max(bounds)<2**50
    a,c=ZZ(data['root_sq']),ZZ(data['content']); T=matrix(ZZ,data['tail']); inner=matrix(ZZ,[data['inner']])
    H=a*(T*T.transpose())-inner.transpose()*inner
    assert all(z%c==0 for z in H.list()); H=H.apply_map(lambda z:z//c)
    initial=int(data['raw']['initial_best']); initial_y=list(map(int,data['raw']['initial_y']))
    assert _witness(data,initial_y)[0]==initial
    input_path=folder/'partition.input'
    content=f'FAST-LIFT-1 {T.nrows()} {a} {c} {initial}\n'
    content+=' '.join(map(str,data['inner']))+'\n'+' '.join(map(str,initial_y))+'\n'
    content+=''.join(' '.join(map(str,row))+'\n' for row in H.rows())
    if input_path.exists(): assert input_path.read_text()==content
    else: input_path.write_text(content)
    identity=dict(prepared_sha256=data['raw']['payload_sha256'],partition_source_sha256=source_hash,
        partition_binary_sha256=_hash_file(binary),input_sha256=_hash_file(input_path),
        bound_file_sha256=_hash_file(bound_file))
    filter_doc=load_filter(data,folder)
    identity['suffix_filter_sha256']=_hash_file(folder/'suffix_filter.json')
    accelerated=None
    if os.environ.get('PARTITION_USE_UNROLLED')=='1' or (HERE/'ENABLE_UNROLLED_ENUMERATION').exists():
        from unrolled_search import build
        accelerated,accelerated_identity=build(T.nrows())
        accelerated_hash=_hash_file(accelerated)
        print('VALIDATED_UNROLLED_SUBTREES',T.nrows(),accelerated_identity,flush=True)
    def execute(suffix):
        result=run_job(accelerated or binary,input_path,suffix,chunk_seconds)
        if result is not None and accelerated is not None:
            result.update(subtree_enumerator='unrolled_enumlib',
                unrolled_source_identity=accelerated_identity,unrolled_binary_sha256=accelerated_hash)
        return result
    db=sqlite3.connect(folder/'partition.sqlite3',timeout=60)
    db.execute('PRAGMA busy_timeout=60000')
    db.execute('PRAGMA synchronous=FULL')
    db.execute('CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY,value TEXT)')
    db.execute('CREATE TABLE IF NOT EXISTS jobs (suffix TEXT PRIMARY KEY,parent TEXT,status TEXT,result TEXT)')
    saved=dict(db.execute('SELECT key,value FROM metadata'))
    if saved:
        assert all(key in identity and identity[key]==value for key,value in saved.items())
        db.executemany('INSERT OR IGNORE INTO metadata VALUES (?,?)',identity.items());db.commit()
    else:
        db.executemany('INSERT INTO metadata VALUES (?,?)',identity.items())
        db.execute('INSERT INTO jobs VALUES (?,?,?,?)',('[]',None,'pending',None));db.commit()
    # Running jobs may have died with a previous invocation/VM. No proof credit.
    db.execute("UPDATE jobs SET status='pending' WHERE status='running'");db.commit()
    # An empty suffix invokes this installation's external enumerator, whose
    # initial top-level task list grows very large at dimension119. Fix the
    # final coordinate first, so every actual search uses the streaming
    # internal enumerator. The same canonical integer box partition is audited.
    root=db.execute("SELECT status FROM jobs WHERE suffix='[]'").fetchone()
    if root and root[0]=='pending':
        assert db.execute('SELECT count(*) FROM jobs').fetchone()[0]==1
        with db:
            db.executemany('INSERT INTO jobs VALUES (?,?,?,?)',
                [(json.dumps(child,separators=(',',':')),'[]','pending',None)
                 for child in children([],bounds)])
            db.execute("UPDATE jobs SET status='split' WHERE suffix='[]'")
    workers=max(1,int(options.get('workers',4)))
    chunk_seconds=float(options.get('partition_chunk_seconds',os.environ.get('PARTITION_CHUNK_SECONDS','120')))
    max_seconds=float(options.get('max_seconds',0))
    futures={}; last_progress=0;last_bulk_prune=0
    def empty_proof(suffix):
        method=None
        if proven_empty(suffix,filter_doc): method='exact_suffix_lower_bound'
        elif len(suffix)<len(bounds):
            lo,hi=candidate_interval(suffix,bounds,filter_doc)
            if lo>hi: method='exact_empty_child_interval'
        if method is None: return None
        return dict(best=str(initial),y=initial_y,nodes=0,leaves=0,seconds=0,precision=None,
                    method=method,empty_filter_sha256=identity['suffix_filter_sha256'])
    def progress():
        counts=dict(db.execute('SELECT status,count(*) FROM jobs GROUP BY status'))
        dump(folder/'partition_progress.json',dict(counts=counts,seconds=time.monotonic()-began,
            prepared_sha256=identity['prepared_sha256'],workers=workers))
        print('PARTITION_PROGRESS',folder,counts,round(time.monotonic()-began,1),flush=True)
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            while True:
                expired=((max_seconds and time.monotonic()-began>=max_seconds)
                    or (folder/'STOP_AFTER_CURRENT_JOBS').exists() or (HERE/'STOP_SUPERVISOR').exists())
                if not expired and time.monotonic()-last_bulk_prune>30:
                    candidates=list(db.execute("SELECT suffix FROM jobs WHERE status='pending' LIMIT 200000"))
                    if len(candidates)>4*workers:
                        proved=[]
                        for (key,) in candidates:
                            proof=empty_proof(json.loads(key))
                            if proof is not None: proved.append((json.dumps(proof),key))
                        if proved:
                            for offset in range(0,len(proved),256):
                                with db:
                                    db.executemany("UPDATE jobs SET status='done',result=? WHERE suffix=? AND status='pending'",
                                                   proved[offset:offset+256])
                            print('EXACT_BULK_EMPTY',len(proved),str(folder),flush=True)
                    last_bulk_prune=time.monotonic()
                if not expired:
                    pending=list(db.execute("SELECT suffix FROM jobs WHERE status='pending' ORDER BY length(suffix),suffix LIMIT ?",
                                            (workers-len(futures),)))
                    for (key,) in pending:
                        suffix=json.loads(key)
                        result=empty_proof(suffix)
                        if result is not None:
                            db.execute("UPDATE jobs SET status='done',result=? WHERE suffix=?",(json.dumps(result),key));db.commit()
                            continue
                        db.execute("UPDATE jobs SET status='running' WHERE suffix=?",(key,));db.commit()
                        futures[pool.submit(execute,suffix)]=key
                if not futures:
                    if not expired and db.execute("SELECT count(*) FROM jobs WHERE status='pending'").fetchone()[0]:
                        continue
                    break
                done,_=concurrent.futures.wait(futures,timeout=1,
                            return_when=concurrent.futures.FIRST_COMPLETED)
                for future in done:
                    key=futures.pop(future); suffix=json.loads(key); result=future.result()
                    with db:
                        if result is not None:
                            assert _witness(data,result['y'])[0]==int(result['best'])<=initial
                            db.execute("UPDATE jobs SET status='done',result=? WHERE suffix=?",(json.dumps(result),key))
                        else:
                            if len(suffix)==len(bounds): raise RuntimeError('A fully fixed integer vector timed out')
                            lo,hi=candidate_interval(suffix,bounds,filter_doc)
                            if lo>hi:
                                proof=empty_proof(suffix);assert proof is not None
                                db.execute("UPDATE jobs SET status='done',result=? WHERE suffix=?",(json.dumps(proof),key))
                            else:
                                split_proof=dict(method='dyadic_child_interval',lower=str(lo),upper=str(hi),
                                                 filter_sha256=identity['suffix_filter_sha256'])
                                db.executemany('INSERT INTO jobs VALUES (?,?,?,?)',
                                    [(json.dumps([j]+suffix,separators=(',',':')),key,'pending',None) for j in range(lo,hi+1)])
                                db.execute("UPDATE jobs SET status='split',result=? WHERE suffix=?",(json.dumps(split_proof),key))
                if time.monotonic()-last_progress>15:
                    progress(); last_progress=time.monotonic()
        progress()
        unfinished=db.execute("SELECT count(*) FROM jobs WHERE status NOT IN ('done','split')").fetchone()[0]
        results=[json.loads(row[0]) for row in db.execute("SELECT result FROM jobs WHERE status='done'")]
        optimum=min(results,key=lambda r:int(r['best'])) if results else dict(best=str(initial),y=initial_y)
        best=int(optimum['best']); y=optimum['y']; sq,witness=_witness(data,y); assert sq==best
        answer=dict(complete=not unfinished,reason='complete' if not unfinished else 'time_limit',
            off_line_sq=str(best),off_line_witness=list(map(str,witness)),quotient_coefficients=list(map(str,y)),
            frontier_jobs=unfinished,nodes=sum(r['nodes'] for r in results),leaves=sum(r['leaves'] for r in results),
            search_seconds=time.monotonic()-began,resource_plan=dict(workers=workers),
            backend='partitioned fplll 5.5.0 error-bounded complete enumeration and exact integer lifts',
            preconditions=preconditions,**identity)
        dump(folder/'partition_search_result.json',answer)
        return answer
    finally:
        db.close()
