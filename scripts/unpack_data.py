"""Verify and unpack selected scientific data into an explicit work directory."""
from pathlib import Path
import argparse,hashlib,json,os,sys,zipfile
ROOT=Path(__file__).resolve().parents[1]

def sha(data):return hashlib.sha256(data).hexdigest()

def unpack(work,case=None,groups=None,skip=None,data_dir=None):
    work=Path(work).resolve();work.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((ROOT/"data/archive_manifest.json").read_text());count=0
    for archive in manifest["archives"]:
        group=archive["group"]
        if case and group!="index" and not group.startswith(case+"/"):continue
        if groups and group!="index" and group.rsplit("/",1)[-1] not in groups:continue
        candidates=([Path(data_dir)] if data_dir else [])+([Path(os.environ["DUAL_COPPERSMITH_DATA"])] if os.environ.get("DUAL_COPPERSMITH_DATA") else [])+[ROOT/".data",ROOT.parent/"dual-coppersmith-release-assets"]
        source=next((folder/archive["asset_name"] for folder in candidates if (folder/archive["asset_name"]).is_file()),None)
        if source is None:raise FileNotFoundError("Download the Release data assets and use --data-dir or DUAL_COPPERSMITH_DATA; missing "+archive["asset_name"])
        if sha(source.read_bytes())!=archive["sha256"]:raise ValueError("Data archive digest mismatch: "+archive["asset_name"])
        expected={item["path"]:item for item in archive["members"]}
        with zipfile.ZipFile(source) as stream:
            if set(stream.namelist())!=set(expected):raise ValueError("Unexpected data archive member")
            for name,item in expected.items():
                if skip and skip(name):continue
                target=(work/name).resolve()
                if work not in target.parents:raise ValueError("Unsafe data archive path")
                data=stream.read(name)
                if len(data)!=item["bytes"] or sha(data)!=item["sha256"]:raise ValueError("Data member digest mismatch: "+name)
                if target.exists():
                    if sha(target.read_bytes())!=item["sha256"]:raise FileExistsError("Use a new work directory; existing generated data differ: "+name)
                else:
                    target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
                count+=1
    return count

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--work",type=Path,required=True);parser.add_argument("--case")
    parser.add_argument("--data-dir",type=Path)
    parser.add_argument("--groups",help="Comma-separated inputs,projection,bounds,optimization")
    args=parser.parse_args();count=unpack(args.work,args.case,set(args.groups.split(",")) if args.groups else None,data_dir=args.data_dir)
    print("DATA_MEMBERS_VERIFIED_AND_UNPACKED",count)
