"""Transfer completed nearest-lift subtrees to a smaller projected SVP sphere.

If a=||r||^2 and P=(a*T-<T,r>*r)/h, every quotient vector has a nearest
original-lattice lift of squared norm at most h^2*||yP||^2/a^2+floor(a/2)^2/a.
Thus a proved original subtree lower bound b excludes projected improvements
below integer incumbent U when a^2*b > h^2*(U-1)+a*floor(a/2)^2.
Original artifacts are read-only. The source SQLite database is snapshotted.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
from contextlib import closing

from projected_minima import HERE, OUT, GRID, name, read, save, matrix, vector, ZZ
from shared import exact_minima as engine
from partition_search import binaries
from suffix_filter import load_filter, proven_empty, candidate_interval


METHOD = "original_nearest_lift_subtree_lower_bound"


def body_hash(path):
    checksum = hashlib.sha256()
    with path.open("rb") as stream:
        stream.readline()
        stream.readline()
        while block := stream.read(1024 * 1024):
            checksum.update(block)
    return checksum.hexdigest()


def excludes(a, h_squared, bound, incumbent):
    return a*a*bound > h_squared*(incumbent-1) + a*(a//2)**2


class Transfer:
    def __init__(self, prepared, stage):
        self.stage = Path(stage)
        path = self.stage / "nearest_lift_transfer.json"
        self.doc = read(path)
        self.doc_hash = engine._hash_file(path)
        source = HERE / self.doc["source_stage"]
        assert engine._hash_file(source / "prepared.json") == self.doc["source_prepared_file_sha256"]
        self.old = engine._decode_preparation(engine._read_document(source / "prepared.json"))
        old, raw = self.old, prepared["raw"]
        assert old["raw"]["payload_sha256"] == self.doc["source_payload_sha256"]
        assert raw["previous_prepared_sha256"] == self.doc["base_payload_sha256"]
        assert raw["payload_sha256"] == self.doc["target_payload_sha256"]
        assert prepared["root_sq"] == 1 and not any(prepared["inner"])
        self.a = int(old["root_sq"])
        h_squared, remainder = divmod(self.a*int(old["content"]), int(prepared["content"]))
        assert remainder == 0 and math.isqrt(h_squared)**2 == h_squared > 0
        self.h_squared = h_squared
        self.incumbent = int(raw["initial_best"])
        assert self.doc["h_squared"] == str(h_squared)
        # The exact normalized Gram used by the source must be the same one.
        T, r = matrix(ZZ, old["tail"]), vector(ZZ, old["root"])
        inner = vector(ZZ, old["inner"])
        assert r.dot_product(r) == self.a and T*r == inner
        H = self.a*T*T.transpose() - matrix(ZZ, len(inner), len(inner), lambda i,j: inner[i]*inner[j])
        assert all(x % old["content"] == 0 for x in H.list())
        H = H.apply_map(lambda x: x//old["content"])
        P = matrix(ZZ, prepared["tail"])
        assert P*P.transpose() == prepared["content"]*H
        assert body_hash(source / "integer.prepared") == body_hash(self.stage / "integer.prepared") == self.doc["native_body_sha256"]
        assert engine._hash_file(source / "partition.input") == self.doc["source_input_sha256"]
        with (source / "partition.input").open() as stream:
            assert stream.readline().split() == ["FAST-LIFT-1", str(T.nrows()), str(self.a),
                                                 str(old["content"]), str(old["raw"]["initial_best"])]
            assert stream.readline().split() == list(map(str, old["inner"]))
            assert stream.readline().split() == list(map(str, old["raw"]["initial_y"]))
            for row in H.rows():
                assert stream.readline().split() == list(map(str, row))
            assert not stream.read().strip()
        snapshot = self.stage / "nearest_lift_source.sqlite3"
        assert engine._hash_file(snapshot) == self.doc["snapshot_sha256"]
        with closing(sqlite3.connect(snapshot.resolve().as_uri()+"?mode=ro", uri=True)) as db:
            metadata = dict(db.execute("SELECT key,value FROM metadata"))
            self.jobs = {key: (status, json.loads(proof) if proof else None)
                         for key,status,proof in db.execute("SELECT suffix,status,result FROM jobs")}
        binary, source_hash, _, _ = binaries()
        assert metadata["partition_source_sha256"] == source_hash
        assert metadata["partition_binary_sha256"] == engine._hash_file(binary)
        assert metadata["prepared_sha256"] == old["raw"]["payload_sha256"]
        assert metadata["input_sha256"] == self.doc["source_input_sha256"]
        self.filter = load_filter(old, source)
        self.filter_hash = engine._hash_file(source / "suffix_filter.json")
        assert metadata["suffix_filter_sha256"] == self.filter_hash
        self.accelerator = None
        self.checked = {}
        self.lower_bound_checks = {}

    def source_lower_bound(self, key):
        if key in self.checked:
            return self.checked[key]
        entry = self.jobs.get(key)
        if not entry or entry[0] != "done":
            return None
        proof = entry[1]
        suffix = json.loads(key)
        assert proof and len(suffix) <= len(self.old["tail"])
        if proof.get("method") == "exact_suffix_lower_bound":
            assert proof["empty_filter_sha256"] == self.filter_hash and proven_empty(suffix, self.filter)
        elif proof.get("method") == "exact_empty_child_interval":
            limits = read((HERE / self.doc["source_stage"]) / "partition_bounds.json")
            lo,hi = candidate_interval(suffix, list(map(int, limits["bounds"])), self.filter)
            assert proof["empty_filter_sha256"] == self.filter_hash and lo > hi
        else:
            assert not proof.get("method"), "Unsupported original proof method"
            if proof.get("subtree_enumerator"):
                if self.accelerator is None:
                    from unrolled_search import build
                    binary, identity = build(len(self.old["tail"]))
                    self.accelerator = identity, engine._hash_file(binary)
                assert proof["subtree_enumerator"] == "unrolled_enumlib"
                assert (proof["unrolled_source_identity"], proof["unrolled_binary_sha256"]) == self.accelerator
        bound = int(proof["best"])
        y = list(map(int, proof["y"]))
        initial = int(self.old["raw"]["initial_best"])
        assert 0 < bound <= initial
        # Nearly all completed source subtrees are empty below the same incumbent.
        # Validate that common witness once, without repeating large matrix products.
        witness_key = bound, tuple(y)
        if witness_key not in self.lower_bound_checks:
            assert engine._witness(self.old, y)[0] == bound
            self.lower_bound_checks[witness_key] = True
        if bound < initial and suffix:
            assert y[-len(suffix):] == suffix
        self.checked[key] = bound
        return bound

    def proof(self, suffix, target_y):
        key = json.dumps(suffix, separators=(",", ":"))
        bound = self.source_lower_bound(key)
        if bound is None or not excludes(self.a, self.h_squared, bound, self.incumbent):
            return None
        return dict(method=METHOD, transfer_sha256=self.doc_hash, source_suffix=key,
            source_original_lower_bound=str(bound), best=str(self.incumbent), y=target_y,
            nodes=0, leaves=0, seconds=0, precision=None)

    def verify(self, suffix, proof):
        assert proof["transfer_sha256"] == self.doc_hash
        assert proof["source_suffix"] == json.dumps(suffix, separators=(",", ":"))
        bound = self.source_lower_bound(proof["source_suffix"])
        assert bound is not None and str(bound) == proof["source_original_lower_bound"]
        assert excludes(self.a, self.h_squared, bound, self.incumbent)
        assert int(proof["best"]) == self.incumbent


def prepare(index):
    case = name(GRID[index])
    base_folder = OUT / case / "full_strong"
    base_stage = base_folder / "certificate" / "pass1"
    base = engine._read_document(base_stage / "prepared.json")
    candidate = read(OUT / "sieving" / case / "best_candidate.json")
    assert candidate["prepared_sha256"] == base["payload_sha256"]
    y = list(map(int, candidate["quotient_coefficients"]))
    value = engine._witness(engine._decode_preparation(base), y)[0]
    identity = engine._digest([base["payload_sha256"], str(value)])[:12]
    target = OUT / case / ("full_sieve_" + identity)
    if target.exists():
        raise FileExistsError(target)
    staging = OUT / case / (".transfer-staging-" + str(os.getpid()))
    stage = staging / "certificate" / "pass1"
    stage.mkdir(parents=True)
    os.link(base_stage / "integer.prepared", stage / "integer.prepared")
    save(staging / "initial_candidate.json", candidate)
    embedding = read(base_folder / "embedding.json")
    embedding.update(radius_source="exactly_checked_quick_svp_candidate",
        base_embedding_sha256=engine._hash_file(base_folder / "embedding.json"),
        initial_candidate_sha256=engine._hash_file(staging / "initial_candidate.json"),
        proof_acceleration="Completed original nearest-lift subtrees, with rigorous projection conversion")
    save(staging / "embedding.json", embedding)
    adjusted = {key: val for key,val in base.items() if key != "payload_sha256"}
    adjusted.update(initial_best=str(value), initial_y=y, previous_prepared_sha256=base["payload_sha256"],
        previous_initial_best=base["initial_best"],
        initial_candidate_sha256=engine._hash_file(staging / "initial_candidate.json"))
    engine._save_document(stage / "prepared.json", adjusted)
    raw = engine._read_document(stage / "prepared.json")
    prepared = engine._decode_preparation(raw)
    source = HERE / embedding["source_preparation"]
    snapshot = stage / "nearest_lift_source.sqlite3"
    with closing(sqlite3.connect((source / "partition.sqlite3").resolve().as_uri()+"?mode=ro", uri=True)) as original:
        with closing(sqlite3.connect(snapshot)) as destination:
            original.backup(destination)
    old = engine._read_document(source / "prepared.json")
    h_squared,remainder = divmod(int(old["root_sq"])*int(old["content"]), int(raw["content"]))
    assert remainder == 0
    save(stage / "nearest_lift_transfer.json", dict(schema="exact-nearest-lift-subtree-transfer-v1",
        source_stage=str(source.relative_to(HERE)), source_prepared_file_sha256=engine._hash_file(source / "prepared.json"),
        source_payload_sha256=old["payload_sha256"], base_payload_sha256=base["payload_sha256"],
        target_payload_sha256=raw["payload_sha256"], h_squared=str(h_squared),
        source_input_sha256=engine._hash_file(source / "partition.input"),
        snapshot_sha256=engine._hash_file(snapshot), native_body_sha256=body_hash(stage / "integer.prepared")))
    context = Transfer(prepared, stage)
    # Initialize ordinary bounds, input and provenance without doing a timed search.
    from partition_search import search
    initialized = search(prepared, stage, None, workers=1, max_seconds=0.000001)
    assert not initialized["complete"]
    limits = list(map(int, read(stage / "partition_bounds.json")["bounds"]))
    filter_doc = load_filter(prepared, stage)
    filter_hash = engine._hash_file(stage / "suffix_filter.json")
    nodes, stack, imported = [], [([], None)], 0
    while stack:
        suffix,parent = stack.pop()
        key = json.dumps(suffix, separators=(",", ":"))
        source_entry = context.jobs.get(key)
        proof = context.proof(suffix, y)
        if proof:
            status = "done"
            imported += 1
        elif proven_empty(suffix, filter_doc):
            status = "done"
            proof = dict(method="exact_suffix_lower_bound", empty_filter_sha256=filter_hash,
                         best=str(value), y=y, nodes=0, leaves=0, seconds=0, precision=None)
        elif source_entry and source_entry[0] == "split":
            lo,hi = candidate_interval(suffix, limits, filter_doc)
            if lo > hi:
                status = "done"
                proof = dict(method="exact_empty_child_interval", empty_filter_sha256=filter_hash,
                             best=str(value), y=y, nodes=0, leaves=0, seconds=0, precision=None)
            else:
                status = "split"
                proof = dict(method="dyadic_child_interval", lower=str(lo), upper=str(hi), filter_sha256=filter_hash)
                stack.extend(([j]+suffix,key) for j in range(lo,hi+1))
        else:
            status = "pending"
            proof = None
        nodes.append((key,parent,status,json.dumps(proof) if proof else None))
    with closing(sqlite3.connect(stage / "partition.sqlite3")) as db:
        with db:
            assert not db.execute("SELECT count(*) FROM jobs WHERE status='running'").fetchone()[0]
            db.execute("DELETE FROM jobs")
            db.executemany("INSERT INTO jobs VALUES (?,?,?,?)", nodes)
    save(stage / "transfer_preparation.json", dict(imported_original_subtrees=imported,
        total_nodes=len(nodes), target_incumbent=str(value), transfer_sha256=context.doc_hash))
    staging.rename(target)
    print("ORIGINAL_SUBTREES_TRANSFERRED", index, imported, "nodes", len(nodes), "folder", target.name, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int, required=True)
    args = parser.parse_args()
    prepare(args.index)
