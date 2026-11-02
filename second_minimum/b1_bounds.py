#!/usr/bin/env python3
"""Build and validate the Appendix B1--B4 table from paper_table_rows.csv."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import math
from pathlib import Path
import sys
from typing import Iterable, TextIO


DELTA0 = 1.02
ROUNDING_TOLERANCE = 0.051
DEFAULT_INPUT = Path(__file__).with_name("data") / "paper_table_rows.csv"
EXPECTED_ORDER = (
    *(("ACD", 2, 1, ell) for ell in range(10, 14)),
    *(("ACD", 2, 2, ell) for ell in range(10, 13)),
    *(("ACD", 3, 2, ell) for ell in range(5, 8)),
    *(("ACD", 3, 3, ell) for ell in range(5, 8)),
    *(("HM", 3, 3, ell) for ell in range(2, 10)),
)
UPPER_BOUND_KEYS = {
    *(("ACD", 2, 1, ell) for ell in range(10, 14)),
    *(("ACD", 3, 2, ell) for ell in range(5, 8)),
}
REQUIRED_COLUMNS = {
    "application",
    "t",
    "m",
    "ell",
    "dim",
    "primal_lr_s",
    "primal_rr_s",
    "dual_s",
    "lambda1_log2",
    "threshold_log2",
    "lambda2_log2",
    "lambda2_kind",
    "lambda2_gh_log2",
    "lambda2_jet_log2",
    "provenance",
}


@dataclass(frozen=True)
class PaperRow:
    application: str
    t: int
    m: int
    ell: int
    dim: int
    primal_lr_s: str
    primal_rr_s: str
    dual_s: str
    lambda1_log2: float
    threshold_log2: float
    lambda2_log2: float
    lambda2_kind: str
    lambda2_gh_log2: float
    lambda2_jet_log2: float
    provenance: str

    @property
    def key(self) -> tuple[str, int, int, int]:
        return self.application, self.t, self.m, self.ell


@dataclass(frozen=True)
class BoundRow:
    paper: PaperRow
    b1_log2: float
    b2_log2: float
    b3_log2: float
    b4_log2: float


def _parse_time(value: str, *, allow_not_run: bool) -> None:
    allowed = {"timeout"}
    if allow_not_run:
        allowed.add("not_run")
    if value in allowed:
        return
    try:
        seconds = float(value)
    except ValueError as exc:
        raise ValueError(f"invalid timing value {value!r}") from exc
    if not math.isfinite(seconds) or seconds < 0:
        raise ValueError(f"invalid timing value {value!r}")


def _paper_row(record: dict[str, str], line_number: int) -> PaperRow:
    try:
        row = PaperRow(
            application=record["application"],
            t=int(record["t"]),
            m=int(record["m"]),
            ell=int(record["ell"]),
            dim=int(record["dim"]),
            primal_lr_s=record["primal_lr_s"],
            primal_rr_s=record["primal_rr_s"],
            dual_s=record["dual_s"],
            lambda1_log2=float(record["lambda1_log2"]),
            threshold_log2=float(record["threshold_log2"]),
            lambda2_log2=float(record["lambda2_log2"]),
            lambda2_kind=record["lambda2_kind"],
            lambda2_gh_log2=float(record["lambda2_gh_log2"]),
            lambda2_jet_log2=float(record["lambda2_jet_log2"]),
            provenance=record["provenance"].strip(),
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"invalid input row at CSV line {line_number}: {exc}") from exc
    _parse_time(row.primal_lr_s, allow_not_run=False)
    _parse_time(row.primal_rr_s, allow_not_run=True)
    _parse_time(row.dual_s, allow_not_run=False)
    if not row.provenance:
        raise ValueError(f"missing provenance at CSV line {line_number}")
    return row


def load_rows(path: Path) -> list[PaperRow]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"missing CSV columns: {', '.join(sorted(missing))}")
        parsed = [
            _paper_row(record, line_number)
            for line_number, record in enumerate(reader, 2)
        ]

    by_key: dict[tuple[str, int, int, int], PaperRow] = {}
    for row in parsed:
        if row.key in by_key:
            raise ValueError(f"duplicate paper row: {row.key}")
        by_key[row.key] = row
    expected = set(EXPECTED_ORDER)
    actual = set(by_key)
    if actual != expected:
        missing = sorted(expected - actual)
        unexpected = sorted(actual - expected)
        raise ValueError(
            f"paper row set mismatch: missing={missing}, unexpected={unexpected}"
        )
    return [by_key[key] for key in EXPECTED_ORDER]


def calculate_bounds(row: PaperRow) -> BoundRow:
    n = row.dim
    log2_lll_factor = math.log2(4.0 / 3.0)
    log2_delta0 = math.log2(DELTA0)
    projected_gh_constant = 0.5 * math.log2(
        (n - 1) / (2.0 * math.pi * math.e)
    )
    log2_det = row.lambda1_log2 + (n - 1) * (
        row.lambda2_gh_log2 - projected_gh_constant
    )
    log2_det_root = log2_det / n
    return BoundRow(
        paper=row,
        b1_log2=row.lambda1_log2 + ((n - 1) / 2.0) * log2_lll_factor,
        b2_log2=log2_det_root + ((n - 1) / 4.0) * log2_lll_factor,
        b3_log2=log2_det_root + n * log2_delta0,
        b4_log2=row.lambda1_log2 + 2 * n * log2_delta0,
    )


def verify_rows(rows: Iterable[PaperRow]) -> list[BoundRow]:
    bounds = [calculate_bounds(row) for row in rows]
    for bound in bounds:
        row = bound.paper
        expected_dim = math.comb(row.ell + row.t, row.t)
        if row.dim != expected_dim:
            raise ValueError(
                f"dimension mismatch for {row.key}: {row.dim} != {expected_dim}"
            )
        expected_kind = (
            "upper_bound" if row.key in UPPER_BOUND_KEYS else "enumerated"
        )
        if row.lambda2_kind != expected_kind:
            raise ValueError(
                f"lambda2 kind mismatch for {row.key}: "
                f"{row.lambda2_kind!r} != {expected_kind!r}"
            )
        if abs(row.threshold_log2 - bound.b4_log2) > ROUNDING_TOLERANCE:
            raise ValueError(
                f"threshold/B4 mismatch for {row.key}: "
                f"{row.threshold_log2:.1f} != {bound.b4_log2:.6f}"
            )
        values = (
            row.lambda1_log2,
            row.threshold_log2,
            row.lambda2_log2,
            row.lambda2_gh_log2,
            row.lambda2_jet_log2,
            bound.b1_log2,
            bound.b2_log2,
            bound.b3_log2,
            bound.b4_log2,
        )
        if not all(math.isfinite(value) and value > 0 for value in values):
            raise ValueError(f"non-positive or non-finite norm for {row.key}")
    return bounds


def _write_csv(bounds: Iterable[BoundRow], handle: TextIO) -> None:
    fields = (
        "application",
        "t",
        "m",
        "ell",
        "dim",
        "lambda1_log2",
        "lambda2_log2",
        "lambda2_kind",
        "b1_log2",
        "b2_log2",
        "b3_log2",
        "b4_log2",
        "b1_below_lambda2",
        "b2_below_lambda2",
        "b3_below_lambda2",
        "b4_below_lambda2",
    )
    writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for bound in bounds:
        row = bound.paper
        writer.writerow(
            {
                "application": row.application,
                "t": row.t,
                "m": row.m,
                "ell": row.ell,
                "dim": row.dim,
                "lambda1_log2": f"{row.lambda1_log2:.1f}",
                "lambda2_log2": f"{row.lambda2_log2:.1f}",
                "lambda2_kind": row.lambda2_kind,
                "b1_log2": f"{bound.b1_log2:.1f}",
                "b2_log2": f"{bound.b2_log2:.1f}",
                "b3_log2": f"{bound.b3_log2:.1f}",
                "b4_log2": f"{bound.b4_log2:.1f}",
                "b1_below_lambda2": bound.b1_log2 < row.lambda2_log2,
                "b2_below_lambda2": bound.b2_log2 < row.lambda2_log2,
                "b3_below_lambda2": bound.b3_log2 < row.lambda2_log2,
                "b4_below_lambda2": bound.b4_log2 < row.lambda2_log2,
            }
        )


def _latex_value(value: float, lambda2_log2: float) -> str:
    rendered = f"{value:.1f}"
    return rf"\textbf{{{rendered}}}" if value < lambda2_log2 else rendered


def _latex_instance(row: PaperRow) -> str:
    if row.application == "ACD":
        return rf"ACD \(({row.t},{row.m})\), \(\ell{{=}}{row.ell}\)"
    return rf"HM \(\ell{{=}}{row.ell}\)"


def _write_latex(bounds: Iterable[BoundRow], handle: TextIO) -> None:
    rows = list(bounds)
    print(r"\begin{tabular}{l c r r r r r r}", file=handle)
    print(r"\hline", file=handle)
    print(r"Instance & \(n\) & \(\lambda_1\) & \(\lambda_2\)", file=handle)
    print(r" & \(B_1\) & \(B_2\) & \(B_3\) & \(B_4\) \\", file=handle)
    print(r"\hline", file=handle)
    previous_group: tuple[object, ...] | None = None
    for bound in rows:
        row = bound.paper
        group = (
            (row.application, row.t, row.m)
            if row.application == "ACD"
            else ("HM",)
        )
        if previous_group is not None and group != previous_group:
            print(r"\hline", file=handle)
        previous_group = group
        values = " & ".join(
            (
                _latex_instance(row),
                str(row.dim),
                f"{row.lambda1_log2:.1f}",
                f"{row.lambda2_log2:.1f}",
                _latex_value(bound.b1_log2, row.lambda2_log2),
                _latex_value(bound.b2_log2, row.lambda2_log2),
                _latex_value(bound.b3_log2, row.lambda2_log2),
                _latex_value(bound.b4_log2, row.lambda2_log2),
            )
        )
        print(values + r" \\", file=handle)
    print(r"\hline", file=handle)
    print(r"\end{tabular}", file=handle)


def _below_counts(bounds: Iterable[BoundRow]) -> tuple[int, int, int, int]:
    rows = list(bounds)
    return tuple(
        sum(getattr(bound, field) < bound.paper.lambda2_log2 for bound in rows)
        for field in ("b1_log2", "b2_log2", "b3_log2", "b4_log2")
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("format", choices=("check", "csv", "latex"))
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    try:
        bounds = verify_rows(load_rows(args.input))
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

    if args.format == "check":
        counts = _below_counts(bounds)
        print(
            f"verified {len(bounds)} rows; below lambda2: "
            f"B1={counts[0]}, B2={counts[1]}, B3={counts[2]}, B4={counts[3]}"
        )
        return 0

    handle: TextIO
    close_handle = False
    if args.output:
        handle = args.output.open("w", encoding="utf-8", newline="")
        close_handle = True
    else:
        handle = sys.stdout
    try:
        if args.format == "csv":
            _write_csv(bounds, handle)
        else:
            _write_latex(bounds, handle)
    finally:
        if close_handle:
            handle.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
