"""Exact, conservative rejection of fixed suffixes outside the search sphere."""
import json
import math
import subprocess
from experiment import HERE,dump
from fast_search import compile_once
from shared.exact_minima import _hash_file


def load_filter(data,folder,recheck=False):
    path=folder/'suffix_filter.json'
    source=HERE/'suffix_filter_wide.cpp'
    if path.exists():
        saved=json.loads(path.read_text())
        matches=[candidate for candidate in (HERE/'suffix_filter.cpp',HERE/'suffix_filter_wide.cpp')
                 if _hash_file(candidate)==saved['source_sha256']]
        assert len(matches)==1
        source=matches[0]
    digest=_hash_file(source)
    binary=HERE/'bin'/('suffix_filter-'+digest[:16])
    compile_once(binary,['g++','-O3','-std=c++17',str(source),'-lgmpxx','-lgmp'])
    identity=dict(source_sha256=digest,binary_sha256=_hash_file(binary),
        preparation_sha256=data['raw']['preparation_sha256'],prepared_sha256=data['raw']['payload_sha256'])
    if path.exists() and not recheck:
        doc=json.loads(path.read_text());assert all(doc[k]==v for k,v in identity.items());return doc
    result=subprocess.run([str(binary),str(folder/'integer.prepared'),str(data['raw']['initial_best'])],
        capture_output=True,text=True,check=True)
    lines=result.stdout.splitlines();tag,d,start,mu_bits,r_bits=lines[0].split()
    assert tag=='SUFFIX-BOUNDS-1' and int(d)==len(data['tail'])
    doc=dict(**identity,dimension=int(d),start=int(start),mu_bits=int(mu_bits),r_bits=int(r_bits),
             diagonal_lower=lines[1].split(),mu_lower=[line.split() for line in lines[2:]])
    assert len(doc['diagonal_lower'])==len(doc['mu_lower'])==int(d)-int(start)
    assert all(len(row)==i for i,row in enumerate(doc['mu_lower']))
    if path.exists(): assert doc==json.loads(path.read_text())
    else: dump(path,doc)
    return doc


def partial_lower(suffix,doc):
    if not suffix: return 0
    # Use only the last 24 fixed coordinates. Omitting earlier nonnegative
    # complete-square terms preserves a lower bound, even for a longer suffix.
    offset=max(len(suffix)-len(doc['diagonal_lower']),0)
    suffix=suffix[offset:]
    first=doc['dimension']-len(suffix)-doc['start']
    mu_den=1<<doc['mu_bits'];normalizer=1<<(2*doc['mu_bits']+doc['r_bits'])
    total=0
    for i,coefficient in enumerate(suffix):
        low=high=coefficient*mu_den
        for j in range(i+1,len(suffix)):
            m=int(doc['mu_lower'][first+j][first+i]);y=suffix[j]
            low+=y*(m if y>=0 else m+1)
            high+=y*(m+1 if y>=0 else m)
        distance=low if low>0 else -high if high<0 else 0
        total+=distance**2*int(doc['diagonal_lower'][first+i])
    return total


def proven_empty(suffix,doc):
    return partial_lower(suffix,doc)>(1<<(2*doc['mu_bits']+doc['r_bits']))


def candidate_interval(suffix,bounds,doc,canonical=True):
    """Conservative integer interval for the next coefficient, proved from
    exact dyadic partial-norm and center bounds. No probability/pruning model.
    """
    dimension=len(bounds)
    assert len(suffix)<dimension
    i=dimension-len(suffix)-1
    lower=-int(bounds[i]);upper=int(bounds[i])
    if canonical and not any(suffix): lower=0
    if i<doc['start']: return lower,upper
    normalizer=1<<(2*doc['mu_bits']+doc['r_bits'])
    remainder=normalizer-partial_lower(suffix,doc)
    if remainder<0: return 1,0
    diagonal=int(doc['diagonal_lower'][i-doc['start']])
    if diagonal==0: return lower,upper
    denominator=1<<doc['mu_bits'];low=high=0
    for j,y in enumerate(suffix,i+1):
        mu=int(doc['mu_lower'][j-doc['start']][i-doc['start']])
        low+=y*(mu if y>=0 else mu+1)
        high+=y*(mu+1 if y>=0 else mu)
    distance=math.isqrt(remainder//diagonal)
    lower=max(lower,-((distance+high)//denominator))
    upper=min(upper,(distance-low)//denominator)
    return lower,upper
