"""Bind every saved exact GSO coefficient to the actual integer lift Gram."""
import json
import subprocess
from sage.all import ZZ,matrix
from experiment import HERE,dump
from fast_search import compile_once
from shared.exact_minima import _hash_file

def validate(data,folder):
    source=HERE/'verify_gram.cpp';source_hash=_hash_file(source)
    binary=HERE/'bin'/('verify_gram-'+source_hash[:16])
    compile_once(binary,['g++','-O3','-std=c++17',str(source),'-lgmpxx','-lgmp'])
    identity=dict(preparation_sha256=data['raw']['preparation_sha256'],
        prepared_sha256=data['raw']['payload_sha256'],source_sha256=source_hash,binary_sha256=_hash_file(binary))
    saved=folder/'gram_validation.json'
    if saved.exists():
        previous=json.loads(saved.read_text())
        if all(previous.get(k)==v for k,v in identity.items()) and previous.get('exact'):
            return previous
    a,c=ZZ(data['root_sq']),ZZ(data['content']);T=matrix(ZZ,data['tail']);inner=matrix(ZZ,[data['inner']])
    H=a*(T*T.transpose())-inner.transpose()*inner
    assert all(z%c==0 for z in H.list());H=H.apply_map(lambda z:z//c)
    gram_file=folder/'audit_gram.input'
    text=f'FAST-LIFT-1 {T.nrows()} {a} {c} {data["raw"]["initial_best"]}\n'
    text+=' '.join(map(str,data['inner']))+'\n'+' '.join(map(str,data['raw']['initial_y']))+'\n'
    text+=''.join(' '.join(map(str,row))+'\n' for row in H.rows())
    gram_file.write_text(text)
    result=subprocess.run([str(binary),str(folder/'integer.prepared'),str(gram_file)],
        text=True,capture_output=True,check=True)
    assert result.stdout.split()==['GRAM_OK',str(T.nrows())]
    record=dict(**identity,exact=True,dimension=T.nrows(),gram_input_sha256=_hash_file(gram_file))
    dump(saved,record)
    return record
