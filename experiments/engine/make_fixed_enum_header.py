"""Derive a separately named, 64-bit single-thread enumlib research header.

The original MIT notice is retained. Original and accepted sources are intact.
"""
import hashlib,json,os
from pathlib import Path
HERE=Path(__file__).resolve().parent
source=HERE/'third_party'/'fplll-5.5.0'/'fplll'/'enum-parallel'
text=(source/'enumeration.h').read_text()
text=text.replace('#include "fplll_types.h"','#include "fixed_fplll_types.h"')
text=text.replace('namespace enumlib','namespace fixed_enumlib')
text=text.replace('//#define SINGLE_THREADED','#define SINGLE_THREADED')
text=text.replace('array<int, N>','array<int64_t, N>')
text=text.replace('inline int myround(', 'inline int64_t myround(')
text=text.replace('return (int)(round', 'return (int64_t)(round')
text=text.replace('int xi        = (int)(yi);','int64_t xi    = (int64_t)(yi);')
text=text.replace('int xi             = (int)(yi);','int64_t xi         = (int64_t)(yi);')
text=text.replace('int xi2        = myround(ci2);','int64_t xi2    = myround(ci2);')
types=(source/'fplll_types.h').read_text().replace('namespace enumlib','namespace fixed_enumlib')
directory=HERE/'generated';directory.mkdir(exist_ok=True)
def write_if_changed(path,content):
    if path.exists() and path.read_text()==content: return
    temporary=path.with_name(path.name+'.tmp-'+str(os.getpid()))
    temporary.write_text(content);temporary.replace(path)
write_if_changed(directory/'fixed_enumlib.h',text)
write_if_changed(directory/'fixed_fplll_types.h',types)
record=dict(upstream_commit='a8dedce384689047daba154bd50d6215e35bf03b',license='MIT',
    modifications=['separate namespace','64-bit coordinate and step storage','single-thread recursive path'],
    files={name:hashlib.sha256((directory/name).read_bytes()).hexdigest()
           for name in ['fixed_enumlib.h','fixed_fplll_types.h']})
write_if_changed(directory/'fixed_enumlib_provenance.json',json.dumps(record,indent=2)+'\n')
print('RESEARCH_ENUM_HEADER_READY',flush=True)
