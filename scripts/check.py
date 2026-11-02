"""Check source integrity, result semantics, and optional Release data."""
from pathlib import Path,PurePosixPath
import argparse,hashlib,json,sys,tempfile,zipfile
from decimal import Decimal
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'second_minimum/experiments'))
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(200000)
from reporting import case_records,log_length,regenerate

def sha(data):return hashlib.sha256(data).hexdigest()
def rational(record):
    value=record['squared_norm']
    result=Fraction(int(value['numerator']),int(value['denominator']))
    assert result>0
    assert abs(Decimal(record['log2_length'])-Decimal(log_length(value)))<Decimal('1e-70')
    return result

def check_results():
    configs=json.loads((ROOT/'second_minimum/data/cases.json').read_text());records=case_records()
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
            assert (path/filename).read_bytes()==(ROOT/'results'/filename).read_bytes(),filename
        assert status==json.loads((ROOT/'results/status.json').read_text())
    return dict(cases=21,original_exact=exact,full_projection_exact=projected_exact,unfinished=5)

def check_source():
    count=0
    for path in ROOT.rglob('*.py'):
        if any(part in ('.git','.work','.data','__pycache__') for part in path.relative_to(ROOT).parts):continue
        compile(path.read_text(encoding='utf-8'),str(path.relative_to(ROOT)),'exec');count+=1
    frozen=json.loads((ROOT/'second_minimum/experiments/unrolled_validated_sources.json').read_text())
    for name,digest in frozen['files'].items():assert sha((ROOT/'second_minimum/experiments'/name).read_bytes())==digest,name
    manifest=json.loads((ROOT/'DATA_MANIFEST.json').read_text())
    for item in manifest['datasets']:
        assert sha((ROOT/item['path']).read_bytes())==item['sha256'],item['path']
    return dict(python_files=count,frozen_native_files=len(frozen['files']),original_datasets_verified=len(manifest['datasets']))

def check_data(directory):
    from unpack_data import check_archive,asset_index
    seen=set();members=0;total=0
    for asset in asset_index():
        path=Path(directory)/asset['asset_name']
        with check_archive(path,asset) as archive:
            for member in archive.infolist():
                assert member.filename not in seen;seen.add(member.filename)
                archive.read(member)  # Verify CRC in addition to the ZIP SHA-256.
                members+=1
        total+=path.stat().st_size
    return dict(assets=len(asset_index()),members=members,compressed_bytes=total)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data-dir',type=Path)
    args=p.parse_args();result=dict(source=check_source(),results=check_results())
    if args.data_dir:result['data']=check_data(args.data_dir)
    print('ARTIFACT_CHECK_OK',json.dumps(result,sort_keys=True),flush=True)

if __name__=='__main__':main()
