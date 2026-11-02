"""Generate and exactly check upper-bound witnesses; never certify minima."""
from shared_bootstrap import *
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from scipy.special import betaincinv
from experiment import HERE,dump
from fast_search import compile_once
from shared import exact_minima as engine

p=argparse.ArgumentParser();p.add_argument('stage',type=Path);p.add_argument('--context',type=int,default=90)
p.add_argument('--threads',type=int,default=1);p.add_argument('--seed',type=int,default=0)
p.add_argument('--seconds',type=float,default=600);args=p.parse_args();stage=args.stage.resolve()
source=HERE/'third_party'/'g6k-c71e084dd003e1ffbfe8c91a7211b0ce9da22245'
archive=source/'kernel'/'.libs'/'libg6k.a'
paths=[HERE/'sieve_candidates.cpp',HERE/'fast_lift.cpp',archive,source/'kernel'/'g6k_config.h']
identity=hashlib.sha256(''.join(engine._hash_file(path) for path in paths).encode()).hexdigest()
binary=HERE/'bin'/('sieve_candidates-'+identity[:16])
prefix=Path(sys.prefix)
compile_once(binary,['/usr/bin/g++','-O3','-std=c++17','-march=native','-pthread',
    str(HERE/'sieve_candidates.cpp'),str(archive),'-I'+str(source),'-I'+str(prefix/'include'),
    '-L'+str(prefix/'lib'),'-Wl,-rpath,'+str(prefix/'lib'),'-lfplll','-lmpfr','-lgmpxx','-lgmp','-lqd'])
doc=engine._read_document(stage/'prepared.json');data=engine._decode_preparation(doc)
assert engine._hash_file(stage/'integer.prepared')==doc['preparation_sha256']
input_path=stage/'partition.input'
if not input_path.exists(): input_path=stage/'fast.input'
assert input_path.exists()
dimension=len(data['tail']);context=min(dimension,args.context)
target=min(2000000,math.ceil(3.2*(4/3)**(context/2)))
bucket=3.2*math.sqrt(target);alpha=math.sqrt(1-betaincinv((context+1)/2,.5,min(.99,bucket/target)))
command=[str(binary),str(input_path),str(source/'spherical_coding'),str(target),str(args.threads),
    str(alpha),str(context),str(args.seed)]
timed_out=False;kernel_returncode=0
try:
    run=subprocess.run(command,text=True,capture_output=True,check=True,timeout=args.seconds)
    output,diagnostics=run.stdout,run.stderr
except subprocess.TimeoutExpired as error:
    timed_out=True
    def decoded(value): return value.decode() if isinstance(value,bytes) else (value or '')
    output,diagnostics=decoded(error.stdout),decoded(error.stderr)
except subprocess.CalledProcessError as error:
    kernel_returncode=error.returncode
    output,diagnostics=error.stdout or '',error.stderr or ''
lines=output.splitlines()
positions=[i for i,line in enumerate(lines) if line.startswith('HEURISTIC_CANDIDATE ')]
if not positions: raise RuntimeError('No candidate checkpoint was emitted: '+diagnostics[-2000:])
position=positions[-1]
tag,squared,size,seconds,checked,improvements=lines[position].split()
y=list(map(int,lines[position+1].split()));norm,witness=engine._witness(data,y)
assert any(y) and norm==int(squared) and norm<=int(doc['initial_best'])
record=dict(certified_minimum=False,role='heuristic upper-bound witness, exactly checked',
    squared_length=str(norm),quotient_coefficients=list(map(str,y)),witness=list(map(str,witness)),
    dimension=dimension,context=context,database_size=int(size),seconds=float(seconds),
    exact_vectors_checked=int(checked),improvements=int(improvements),
    improvement_bits=(math.log2(int(doc['initial_best']))-math.log2(norm))/2,
    prepared_sha256=doc['payload_sha256'],kernel_commit='c71e084dd003e1ffbfe8c91a7211b0ce9da22245',
    source_identity=identity,binary_sha256=engine._hash_file(binary),seed=args.seed,threads=args.threads,
    time_limit_reached=timed_out,kernel_returncode=kernel_returncode,
    completed_sieve_rounds=diagnostics.splitlines())
rows={}
for line in lines:
    if not line.startswith('HEURISTIC_BASIS_ROW '): continue
    _,index,length,*coordinates=line.split();index=int(index);length=float(length)
    coefficients=list(map(int,coordinates))
    assert len(coefficients)==dimension and any(coefficients)
    if index not in rows or length<rows[index]['heuristic_projected_length']:
        rows[index]=dict(index=index,heuristic_projected_length=length,
                         quotient_coefficients=list(map(str,coefficients)))
record['public_preconditioning_rows']=list(rows.values())
completed=stage/'partition_search_result.json'
if completed.exists():
    exact=json.loads(completed.read_text())
    if exact.get('complete'):
        assert norm>=int(exact['off_line_sq']), 'Candidate contradicts a completed certificate'
        record['matches_existing_certified_value']=norm==int(exact['off_line_sq'])
dump(stage/f'sieve_candidate_seed{args.seed}.json',record)
print('SIEVE_UPPER_BOUND',dimension,'improvement_bits',record['improvement_bits'],
      'seconds',record['seconds'],'NOT A MINIMUM CERTIFICATE',flush=True)
