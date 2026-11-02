"""Reject inconsistent saved exact decompositions, not only accept valid ones."""
from shared_bootstrap import *
from pathlib import Path
import json
import subprocess
import tempfile
from sage.all import ZZ,matrix
from shared import exact_minima as engine
from gram_validation import validate

HERE=Path(__file__).resolve().parent
records=[]
with tempfile.TemporaryDirectory(prefix='gram-checker-tests-',dir=HERE) as temporary:
    for rank,scale in [(2,1),(4,1),(4,2**4096)]:
        B=scale*matrix(ZZ,rank,rank,lambda i,j: 11+i if i==j else (2*i-3*j if j<i else 0))
        folder=Path(temporary)/(str(rank)+'-'+str(scale.bit_length()))
        engine.prepare_successive_minima(B,folder)
        stage=folder/'pass1';data=engine._decode_preparation(engine._read_document(stage/'prepared.json'))
        proof=validate(data,stage)
        binary=HERE/'bin'/('verify_gram-'+proof['source_sha256'][:16])
        lines=(stage/'integer.prepared').read_text().splitlines()
        delta=lines[2].split();delta[-1]=str(int(delta[-1])+1)
        lines[2]=' '.join(delta)
        bad=stage/'bad.prepared';bad.write_text('\n'.join(lines)+'\n')
        result=subprocess.run([str(binary),str(bad),str(stage/'audit_gram.input')],capture_output=True,text=True)
        assert result.returncode and 'pivot mismatch' in result.stderr
        if rank>2:
            lines=(stage/'integer.prepared').read_text().splitlines()
            lower=lines[3].split();lower[0]=str(int(lower[0])+1);lines[3]=' '.join(lower)
            bad.write_text('\n'.join(lines)+'\n')
            result=subprocess.run([str(binary),str(bad),str(stage/'audit_gram.input')],capture_output=True,text=True)
            assert result.returncode and 'lower coefficient mismatch' in result.stderr
        records.append(dict(rank=rank,scale_bits=scale.bit_length(),valid_accepted=True,corrupt_rejected=True))
(HERE/'logs'/'gram-checker-tests.json').write_text(json.dumps(records,indent=2)+'\n')
print('GRAM_CHECKER_TESTS',records,flush=True)
