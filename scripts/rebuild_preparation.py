"""Rebuild an omitted deterministic Bareiss file from an exact saved Gram."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/"baseline"))
from shared import exact_minima as engine


def rebuild(stage,threads=1):
    from sage.all import ZZ,matrix,vector
    stage=Path(stage).resolve();raw=engine._read_document(stage/"prepared.json")
    target=stage/"integer.prepared"
    if target.exists():
        if engine._hash_file(target)!=raw["preparation_sha256"]:raise ValueError("Existing preparation differs from its archived digest")
        return False
    source=stage/"integer.input"
    if not source.exists():
        T=matrix(ZZ,raw["tail"]);inner=vector(ZZ,raw["inner"]);a=ZZ(raw["root_sq"]);c=ZZ(raw["content"])
        H=a*T*T.transpose()-matrix(ZZ,T.nrows(),T.nrows(),lambda i,j:inner[i]*inner[j])
        assert all(value%c==0 for value in H.list());H=H.apply_map(lambda value:value//c)
        with source.open("w") as stream:
            stream.write(f"EXACT-GMP-INPUT-1 {T.nrows()} {a} {c}\n"+" ".join(map(str,inner))+"\n")
            for row in H.rows():stream.write(" ".join(map(str,row))+"\n")
    native=engine._build_native(stage.parent)
    temporary=stage/("integer.prepared.rebuilt-"+str(os.getpid()))
    subprocess.run([str(native),"prepare",str(source),str(temporary),str(threads)],check=True)
    if engine._hash_file(temporary)!=raw["preparation_sha256"]:
        temporary.unlink();raise ValueError("Regenerated preparation does not match the archived digest")
    temporary.replace(target);return True


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("stage",type=Path);parser.add_argument("--threads",type=int,default=1)
    args=parser.parse_args();print("PREPARATION_REBUILT" if rebuild(args.stage,args.threads) else "PREPARATION_ALREADY_VERIFIED")
