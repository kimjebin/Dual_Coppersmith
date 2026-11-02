"""Check source integrity, result semantics, and optional Release data."""
from pathlib import Path,PurePosixPath
import argparse,hashlib,json,sys,tempfile,zipfile
from decimal import Decimal
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(200000)
from core.reporting import case_records,log_length,regenerate

def sha(data):return hashlib.sha256(data).hexdigest()
def rational(record):
    value=record['squared_norm']
    result=Fraction(int(value['numerator']),int(value['denominator']))
    assert result>0
    assert abs(Decimal(record['log2_length'])-Decimal(log_length(value)))<Decimal('1e-70')
    return result

def check_results():
    configs=json.loads((ROOT/'configs/cases_21.json').read_text());records=case_records()
    assert len(configs)==len(records)==21
    assert {r['case_id'] for r in configs}=={r['case_id'] for r in records}
    exact=0;projected_exact=0
    for record in records:
        original=record['original'];first=rational(original['lambda1']);second=original['lambda2']
        lo=rational(second['lower']);assert first<=lo
        if second['exact']:
            assert rational(second['exact'])==lo;exact+=1
        for upper in second['upper_candidates']:
            if upper.get('squared_norm'):assert rational(upper)>=lo
        full=record['projection']['full'];jet=record['projection']['jet']
        assert rational(full)>0 and rational(jet)>0
        assert jet['value_kind']=='archived_audited_exact'
        projected_exact+=full['value_kind']=='archived_audited_exact'
        if second['exact'] is None:assert full['value_kind']=='checked_candidate_upper_bound'
        assert record['predictions']['value_kind']=='heuristic_prediction'
        if record.get('optimization'):assert record['optimization']['global_optimality_proven'] is False
    assert exact==projected_exact==16
    with tempfile.TemporaryDirectory(prefix='dual-tables-') as temporary:
        path=Path(temporary);status=regenerate(path)
        for filename in ('all_21_comparison.csv','unfinished_5_comparison.csv','unfinished_5_comparison.md'):
            assert (path/filename).read_bytes()==(ROOT/'results/tables'/filename).read_bytes(),filename
        assert status==json.loads((ROOT/'results/status.json').read_text())
    return dict(cases=21,original_exact=exact,full_projection_exact=projected_exact,unfinished=5)

def check_source():
    count=0
    for path in ROOT.rglob('*.py'):
        if any(part in ('.git','.work','.data','__pycache__') for part in path.relative_to(ROOT).parts):continue
        compile(path.read_text(encoding='utf-8'),str(path.relative_to(ROOT)),'exec');count+=1
    frozen=json.loads((ROOT/'experiments/engine/unrolled_validated_sources.json').read_text())
    for name,digest in frozen['files'].items():assert sha((ROOT/'experiments/engine'/name).read_bytes())==digest,name
    manifest_path=ROOT/'MANIFEST.json'
    if manifest_path.exists():
        source=json.loads(manifest_path.read_text())
        for item in source['files']:
            path=ROOT/item['path'];assert path.is_file() and sha(path.read_bytes())==item['sha256'],item['path']
    return dict(python_files=count,frozen_native_files=len(frozen['files']),source_manifest_checked=manifest_path.exists())

def check_data(directory):
    manifest=json.loads((ROOT/'data/archive_manifest.json').read_text());seen=set();members=0;total=0
    for asset in manifest['archives']:
        name=asset['asset_name'];assert PurePosixPath(name).name==name and name.endswith('.zip')
        path=Path(directory)/name;data=path.read_bytes()
        assert len(data)==asset['bytes'] and sha(data)==asset['sha256'],name
        total+=len(data);expected={m['path']:m for m in asset['members']}
        with zipfile.ZipFile(path) as archive:
            assert len(archive.namelist())==len(expected) and set(archive.namelist())==set(expected)
            for name,item in expected.items():
                relative=PurePosixPath(name)
                assert not relative.is_absolute() and '..' not in relative.parts and '\\' not in name
                assert name not in seen;seen.add(name)
                value=archive.read(name)
                assert len(value)==item['bytes'] and sha(value)==item['sha256'],name
                members+=1
    provenance=json.loads((ROOT/'data/provenance.json').read_text())
    assert seen=={p['path'] for p in provenance['files']}
    return dict(assets=len(manifest['archives']),members=members,compressed_bytes=total)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data-dir',type=Path)
    args=p.parse_args();result=dict(source=check_source(),results=check_results())
    if args.data_dir:result['data']=check_data(args.data_dir)
    print('SUBMISSION_CHECK_OK',json.dumps(result,sort_keys=True),flush=True)

if __name__=='__main__':main()
