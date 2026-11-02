"""Export optimized projection bases and exact original-lattice lifts."""
import json
import argparse
import hashlib
from pathlib import Path
import zipfile
from sage.all import ZZ,QQ,matrix,vector,identity_matrix
from optimize_bound_basis import DEST,INDICES,load
from projected_minima import HERE,OUT,GRID,name,read,save,rows,norm2
from shared import exact_minima as engine


def export(index):
    folder,P=load(index)[:2];basis=read(folder/"best_basis.json");C=matrix(ZZ,basis["transform_rows"])
    baseline=read(folder/"baseline.json");embedding=read(HERE/baseline["source_embedding"])
    source=HERE/embedding["source_preparation"];old=engine._read_document(source/"prepared.json")
    assert engine._hash_file(source/"prepared.json")==embedding["source_preparation_sha256"]
    r=vector(ZZ,old["root"]);T=matrix(ZZ,old["tail"]);a=norm2(r);d=T.nrows();n=d+1
    physical_root=vector(ZZ,read(OUT/name(GRID[index])/"geometry.json")["root"])
    j=next(j for j,z in enumerate(r) if z);g=abs(physical_root[j]//r[j])
    assert g>0 and g*r in (physical_root,-physical_root)
    ambient=engine._read_document(source/"ambient_lll.json");quotient=engine._read_document(source/"quotient_lll.json")
    U=matrix(ZZ,ambient["transform_rows"]);pre=matrix(ZZ,ambient["preconditioner_rows"])
    red=matrix(ZZ,ambient["reduced_rows"]);V=matrix(ZZ,quotient["transform_rows"])
    expected=matrix(ZZ,engine._read_document(source/"input_basis.json")["rows"])
    assert U*pre*expected==red
    order=sorted(range(n),key=lambda i:norm2(red.row(i)))
    change=identity_matrix(ZZ,n);change[1:,1:]=V
    for i,k in enumerate(old["root_shifts"]): change[i+1,0]=-ZZ(k)
    full=identity_matrix(ZZ,n);full[1:,1:]=C
    transform=full*change*identity_matrix(ZZ,n).matrix_from_rows(order)*U*pre
    assert abs(transform.det())==1
    lifted=matrix(ZZ,[g*r]+list((g*C*T).rows()))
    assert lifted==g*transform*expected
    scale=QQ(basis["scale_numerator"])/ZZ(basis["scale_denominator"])
    CT=C*T
    assert a*scale*matrix(ZZ,basis["rows"])==g*(a*CT-matrix(ZZ,[z*r for z in CT*r]))
    target=folder/"best_original_basis.json"
    save(target,dict(schema="optimized-original-lattice-basis-v1",case=name(GRID[index]),
        rows=rows(lifted),transform_from_source_input_rows=rows(transform),
        source_input=str((source/"input_basis.json").relative_to(HERE)),
        source_input_sha256=engine._hash_file(source/"input_basis.json"),source_input_physical_integer_scale=str(g),
        projection_basis=str((folder/"best_basis.json").relative_to(HERE)),
        projection_basis_sha256=engine._hash_file(folder/"best_basis.json"),
        first_row_is_planted_root_up_to_sign=True,unimodular_and_projection_identities_verified=True))
    print("OPTIMIZED_ORIGINAL_BASIS_EXPORTED",index,n,flush=True)


def bundle(reuse_exports=False):
    files=[]
    for index in INDICES:
        folder=DEST/name(GRID[index]);result=read(folder/"result.json");audit=read(folder/"audit.json")
        assert audit["passed"] and audit["result_sha256"]==engine._hash_file(folder/"result.json")
        if reuse_exports:
            original=read(folder/"best_original_basis.json")
            assert original["unimodular_and_projection_identities_verified"]
            assert original["projection_basis_sha256"]==engine._hash_file(folder/"best_basis.json")
            assert original["source_input_sha256"]==engine._hash_file(HERE/original["source_input"])
        else: export(index)
        files.extend(folder/file for file in ("best_basis.json","best_original_basis.json","result.json","audit.json","baseline.json"))
        if result["selected_certificate"]:
            proof=HERE/result["selected_certificate"]
            files.extend(path for path in proof.iterdir() if path.suffix in (".json",".input"))
    files.extend(DEST/file for file in ("comparison.csv","comparison.json"))
    manifest=dict(schema="optimized-bound-basis-bundle-v1",scope="Five unfinished original B cases, seed 0",
        definitions={"best_basis.json":"Projected lattice basis, multiply rows by the stored positive rational scale",
                     "best_original_basis.json":"Integer basis of the original lattice, first row is the planted root up to sign",
                     "lower_bound":"Certified L <= lambda_1(projection) <= lambda_2(original); not an exact lambda_2 claim"},
        global_optimality_proven=False,files={str(path.relative_to(HERE)):engine._hash_file(path) for path in files})
    save(DEST/"bundle_manifest.json",manifest);files.append(DEST/"bundle_manifest.json")
    target=DEST/"optimized_bases.zip"
    with zipfile.ZipFile(target,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for path in files: archive.write(path,str(path.relative_to(HERE)))
    print("OPTIMIZED_BASES_BUNDLE_SAVED",target,target.stat().st_size,flush=True)


def verify_bundle():
    with zipfile.ZipFile(DEST/"optimized_bases.zip") as archive:
        manifest=json.loads(archive.read("basis_bound_optimization/bundle_manifest.json"))
        for filename,checksum in manifest["files"].items():
            assert hashlib.sha256(archive.read(filename)).hexdigest()==checksum
    print("OPTIMIZED_BASES_ARCHIVE_HASHES_OK",len(manifest["files"]),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--reuse-exports",action="store_true")
    parser.add_argument("--verify-only",action="store_true");args=parser.parse_args()
    if args.verify_only: verify_bundle()
    else: bundle(args.reuse_exports)
