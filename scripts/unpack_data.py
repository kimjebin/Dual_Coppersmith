"""Check companion ZIP hashes and extract selected data into a work directory."""
from pathlib import Path,PurePosixPath
import argparse,hashlib,json,os,stat,zipfile
ROOT=Path(__file__).resolve().parents[1]

def asset_index():
    return json.loads((ROOT/'second_minimum/data/release_assets.json').read_text())['archives']

def check_archive(path,asset):
    path=Path(path);data=path.read_bytes()
    if len(data)!=asset['bytes'] or hashlib.sha256(data).hexdigest()!=asset['sha256']:
        raise ValueError('Release asset digest mismatch: '+path.name)
    archive=zipfile.ZipFile(path)
    try:
        members=archive.infolist();seen=set()
        if len(members)!=asset['member_count']:raise ValueError('Unexpected archive member count')
        for member in members:
            name=member.filename;relative=PurePosixPath(name)
            if (not name or relative.is_absolute() or '..' in relative.parts or '.git' in relative.parts
                or '\\' in name or ':' in name or name!=relative.as_posix() or name in seen
                or member.is_dir() or stat.S_ISLNK(member.external_attr>>16)):
                raise ValueError('Unsafe or duplicate archive member: '+name)
            seen.add(name)
    except BaseException:
        archive.close();raise
    return archive

def unpack(work,case=None,groups=None,skip=None,data_dir=None):
    work=Path(work).resolve()
    if work==ROOT:raise ValueError('Choose an isolated work directory')
    work.mkdir(parents=True,exist_ok=True);count=0
    candidates=([Path(data_dir)] if data_dir else [])
    if os.environ.get('DUAL_COPPERSMITH_DATA'):candidates.append(Path(os.environ['DUAL_COPPERSMITH_DATA']))
    candidates += [ROOT/'.data',ROOT.parent/'dual-coppersmith-release-assets']
    for asset in asset_index():
        group=asset['group']
        if case and group!='index' and not group.startswith(case+'/'):continue
        if groups and group!='index' and group.rsplit('/',1)[-1] not in groups:continue
        source=next((folder/asset['asset_name'] for folder in candidates if (folder/asset['asset_name']).is_file()),None)
        if source is None:raise FileNotFoundError('Set DUAL_COPPERSMITH_DATA or --data-dir; missing '+asset['asset_name'])
        with check_archive(source,asset) as archive:
            for member in archive.infolist():
                name=member.filename
                if skip and skip(name):continue
                target=(work/name).resolve()
                if work not in target.parents:raise ValueError('Extraction path escapes the work directory')
                data=archive.read(member)  # ZIP SHA-256 binds the member; read also checks CRC.
                if target.exists():
                    if hashlib.sha256(target.read_bytes()).digest()!=hashlib.sha256(data).digest():
                        raise FileExistsError('Use a new work directory or --no-unpack; generated data differ: '+name)
                else:
                    target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
                count+=1
    return count

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True)
    p.add_argument('--case');p.add_argument('--data-dir',type=Path);p.add_argument('--groups')
    args=p.parse_args()
    print('DATA_MEMBERS_VERIFIED_AND_UNPACKED',unpack(args.work,args.case,
        set(args.groups.split(',')) if args.groups else None,data_dir=args.data_dir))
