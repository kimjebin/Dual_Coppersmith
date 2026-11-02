"""Measure RHF on saved bases without changing experiments or certificates.

The conventional RHF uses ||b1|| = delta_n**n * det(L)**(1/n).
The n-1 convention in Lemma 2.2 of the supplied paper is also recorded.
Neither measured determinant-based RHF is an SVP approximation guarantee.
All logarithms are evaluated with 80 decimal digits. Determinants use the
exact N/R exponents in audit.json, cross-checked against basis metadata.
"""
from __future__ import annotations

import csv
import hashlib
import json
from decimal import Decimal, localcontext
from pathlib import Path
from paths import ROOT

HERE = Path(__file__).resolve().parent
VARIANTS = ("B",)
GRID = (
    [("acd", 2, 1, ell) for ell in range(10, 14)]
    + [("acd", 2, 2, ell) for ell in range(10, 13)]
    + [("acd", 3, m, ell) for m in (2, 3) for ell in range(5, 8)]
    + [("hm", 3, 3, ell) for ell in range(2, 10)]
)


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    out = HERE / "rhf_analysis"
    out.mkdir(exist_ok=True)
    certificates_path = HERE / "certified_minima.json"
    certificates = {
        (doc["case"], doc["variant"]): doc for doc in read(certificates_path)
    }
    paper_path = ROOT / "second_minimum/data/paper_table_rows.csv"
    with paper_path.open(encoding="utf-8", newline="") as handle:
        paper = {
            (row["application"].lower(), int(row["t"]), int(row["m"]), int(row["ell"])): row
            for row in csv.DictReader(handle)
        }
    records = []
    source_hashes = {}
    with localcontext() as ctx:
        ctx.prec = 80
        ln2 = Decimal(2).ln()

        def log2(value):
            return Decimal(value).ln() / ln2

        def pow2(value):
            return (value * ln2).exp()

        def fmt(value):
            return "" if value is None else format(value, ".18f")

        for app, t, m, ell in GRID:
            case = f"{app}-t{t}-m{m}-ell{ell}-seed0"
            audit_path = HERE / "results" / case / "audit.json"
            audit = read(audit_path)
            source_hashes[str(audit_path.relative_to(HERE))] = digest(audit_path)
            instance = audit["instance"]
            modulus, bound = int(instance["modulus"]), int(instance["bound"])
            audits = {row["variant"]: row for row in audit["records"]}
            folder = "flag_results" if app == "hm" else "results"
            for variant in VARIANTS:
                path = HERE / folder / case / f"{variant}.json"
                raw = path.read_bytes()
                doc = json.loads(raw)
                source_hashes[str(path.relative_to(HERE))] = hashlib.sha256(raw).hexdigest()
                basis = doc["reduced_rows"]
                n = len(basis)
                assert n == doc["dim"] == audits[variant]["dim"]
                assert all(len(row) == n for row in basis)
                b1_sq = sum(int(value) ** 2 for value in basis[0])
                assert b1_sq == int(doc["b1_sq"])
                a = audits[variant]["determinant_N_exponent"]
                b = audits[variant]["determinant_R_exponent"]
                assert sum(doc["exponents"]) == a
                assert sum(t - sum(key) for key in instance["monomial_keys"]) == b
                det_root_log = (a * log2(modulus) + b * log2(bound)) / n
                b1_log = log2(b1_sq) / 2
                assert abs(b1_log - Decimal(str(doc["log2_b1"]))) < Decimal("1e-8")
                gap = b1_log - det_root_log
                rhf_n = pow2(gap / n)
                rhf_nm1 = pow2(gap / (n - 1))
                cert = certificates[case, variant]
                first_sq = cert.get("lambda1_sq")
                first_log = log2(int(first_sq)) / 2 if first_sq is not None else None
                root_log = log2(int(cert["root_sq"])) / 2
                threshold_n = first_log + 2 * gap if first_log is not None else None
                threshold_nm1 = first_log + Decimal(2 * n) * gap / (n - 1) if first_log is not None else None
                threshold_fixed = first_log + 2 * n * log2(Decimal("1.02")) if first_log is not None else None
                # This is an a posteriori normalization of the measured SVP
                # approximation ratio, not an independent recovery criterion.
                svp_delta = pow2((b1_log - first_log) / (2 * n)) if first_log is not None else None
                svp_threshold = b1_log if first_log is not None else None
                matches = b1_sq == int(first_sq) if first_sq is not None else ""
                records.append(dict(
                    application=app.upper(), t=t, m=m, ell=ell, dim=n, variant=variant,
                    source=str(path.relative_to(HERE)), backend=doc["backend"],
                    lovasz_delta=doc["delta"], determinant_N_exponent=a, determinant_R_exponent=b,
                    log2_b1=fmt(b1_log), log2_det_root=fmt(det_root_log),
                    delta_rhf_n=fmt(rhf_n), delta_rhf_n_minus_1=fmt(rhf_nm1),
                    log2_lambda1_certified=fmt(first_log), b1_equals_certified_lambda1=matches,
                    b1_on_planted_root=doc["b1_on_root"], log2_planted_root=fmt(root_log),
                    log2_threshold_rhf_n=fmt(threshold_n),
                    log2_threshold_rhf_n_minus_1=fmt(threshold_nm1),
                    log2_threshold_fixed_1_02=fmt(threshold_fixed),
                    delta_svp_a_posteriori=fmt(svp_delta),
                    log2_threshold_svp_a_posteriori=fmt(svp_threshold),
                    log2_lambda2_certified=cert.get("log2_lambda2", ""),
                    paper_lambda1_log2=paper[app, t, m, ell]["lambda1_log2"] if variant == "B" else "",
                    paper_threshold_log2=paper[app, t, m, ell]["threshold_log2"] if variant == "B" else "",
                ))
            print(f"READ {case}: one saved original basis", flush=True)
        assert len(records) == 21
        base = [row for row in records if row["variant"] == "B"]
        acd = [row for row in base if row["application"] == "ACD"]
        assert len(base) == 21 and len(acd) == 13
        assert all(row["b1_equals_certified_lambda1"] for row in base)
        write_csv(out / "base_21.csv", base)
        write_csv(out / "table1_acd_13.csv", acd)
        source_hashes[str(certificates_path.relative_to(HERE))] = digest(certificates_path)
        metadata = dict(
            precision_decimal_digits=80,
            conventional_rhf="delta_n=(norm(b1)/det(L)**(1/n))**(1/n)",
            paper_lemma_convention="delta_n_minus_1=(norm(b1)/det(L)**(1/n))**(1/(n-1))",
            threshold="log2(lambda1)+2*n*log2(delta)",
            svp_a_posteriori="delta_svp=(norm(b1)/lambda1)**(1/(2*n)); uses a certified lambda1 and is not predictive",
            determinant="det(L)=N**a * R**b; a,b from audit.json and checked against basis exponents and monomials",
            provenance="Fresh saved fpLLL bases at Lovasz delta=0.8; these are not the original paper's dual/flatter timing outputs.",
            interpretation="Empirical determinant-based RHF can be below 1 on planted lattices. It cannot be substituted into a worst-case SVP approximation bound as a recovery certificate.",
            counts=dict(inputs=21, bases=21, base_b1_matches_certified_lambda1=21),
            sources_sha256=source_hashes,
        )
        (out / "provenance.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        print("\nACD table: (t,m,ell) n delta_n delta_(n-1) lambda1_log2 threshold_n threshold_(n-1) fixed_1.02", flush=True)
        for row in acd:
            print(f"({row['t']},{row['m']},{row['ell']}) {row['dim']} "
                  f"{float(row['delta_rhf_n']):.8f} {float(row['delta_rhf_n_minus_1']):.8f} "
                  f"{float(row['log2_lambda1_certified']):.6f} "
                  f"{float(row['log2_threshold_rhf_n']):.6f} "
                  f"{float(row['log2_threshold_rhf_n_minus_1']):.6f} "
                  f"{float(row['log2_threshold_fixed_1_02']):.6f}")
        print(f"\nWROTE {out}")


if __name__ == "__main__":
    main()
