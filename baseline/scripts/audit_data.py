#!/usr/bin/env python3
"""Validate the data used by the active ACD and HM paper tables."""

import argparse
import csv
import hashlib
import json
from math import comb, log2
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ACD_RAW = ROOT / "acd" / "results" / "raw" / "acd_paper.jsonl"
TABLE_ROWS = ROOT / "second_minimum" / "data" / "paper_table_rows.csv"
ACD_GRID = {
    *((2, 1, ell) for ell in range(10, 14)),
    *((2, 2, ell) for ell in range(10, 13)),
    *((3, 2, ell) for ell in range(5, 8)),
    *((3, 3, ell) for ell in range(5, 8)),
}
HM_GRID = {(3, 3, ell) for ell in range(2, 10)}
DENSE_UPPER_BOUND = {
    *((2, 1, ell) for ell in range(10, 14)),
    *((3, 2, ell) for ell in range(5, 8)),
}


def load_jsonl(path: Path) -> list[dict[str, object]]:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise ValueError(f"{path}:{line_number}: {exc}") from exc
    return rows


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def verify_manifest() -> bool:
    failed = False
    manifest = json.loads((ROOT / "DATA_MANIFEST.json").read_text(encoding="utf-8"))
    for dataset in manifest["datasets"]:
        path = ROOT / dataset["path"]
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        records = len(load_jsonl(path)) if path.suffix == ".jsonl" else len(load_csv(path))
        if digest != dataset["sha256"] or records != dataset["records"]:
            failed = True
            print(f"MANIFEST MISMATCH: {dataset['path']}")
    print(f"Data manifest: {'FAILED' if failed else 'all hashes and counts match'}")
    return failed


def audit_acd_raw() -> tuple[bool, dict[tuple[int, int, int], dict[str, object]]]:
    rows = load_jsonl(ACD_RAW)
    by_key = {(int(r["t"]), int(r["m"]), int(r["ell"])): r for r in rows}
    failed = len(rows) != len(by_key) or set(by_key) != ACD_GRID
    for key, row in by_key.items():
        expected_dim = comb(key[2] + key[0], key[2])
        failed |= int(row["dim"]) != expected_dim
        failed |= row.get("method") != "dual" or row.get("backend") != "flatter"
        failed |= row.get("status") != "ok" or row.get("recovered") is not True
    print(f"ACD raw: {len(set(by_key) & ACD_GRID)}/{len(ACD_GRID)} active rows")
    return failed, by_key


def audit_table_rows(acd_raw: dict[tuple[int, int, int], dict[str, object]]) -> bool:
    rows = load_csv(TABLE_ROWS)
    expected = {("ACD", *key) for key in ACD_GRID} | {("HM", *key) for key in HM_GRID}
    actual = {
        (row["application"], int(row["t"]), int(row["m"]), int(row["ell"]))
        for row in rows
    }
    failed = len(rows) != len(actual) or actual != expected
    for row in rows:
        key = (int(row["t"]), int(row["m"]), int(row["ell"]))
        failed |= int(row["dim"]) != comb(key[2] + key[0], key[2])
        expected_kind = (
            "upper_bound"
            if row["application"] == "ACD" and key in DENSE_UPPER_BOUND
            else "enumerated"
        )
        failed |= row["lambda2_kind"] != expected_kind
        lambda1 = float(row["lambda1_log2"])
        threshold = lambda1 + 2 * int(row["dim"]) * log2(1.02)
        failed |= round(threshold, 1) != float(row["threshold_log2"])
        if row["application"] == "ACD":
            failed |= round(float(acd_raw[key]["time_s"]), 1) != float(row["dual_s"])
    print(f"Paper table inventory: {len(actual & expected)}/{len(expected)} active rows")
    return failed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--require-original-raw",
        action="store_true",
        help="fail because the original primal/HM/minimum logs are unavailable",
    )
    args = parser.parse_args()
    failed = verify_manifest()
    acd_failed, acd_raw = audit_acd_raw()
    failed |= acd_failed
    failed |= audit_table_rows(acd_raw)
    print("Known provenance gap: original primal ACD, HM, and long minimum logs unavailable")
    if args.require_original_raw:
        failed = True
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
