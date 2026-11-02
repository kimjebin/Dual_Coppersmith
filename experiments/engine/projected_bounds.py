"""Rigorous projection-length intervals from previously certified lambda2.

For a primitive root r, every nonzero projected vector has a lift with root
component at most ||r||/2. Thus lambda2(L)^2-||r||^2/4 is a lower bound on
lambda1(pi(L))^2. Projection of a certified lambda2 witness is an upper bound.
This certifies displayed decimal values when exact projected SVP is unfinished;
it is explicitly not an exact rational squared-norm certificate.
"""
from sage.all import RealIntervalField
from projected_minima import OUT, HERE, GRID, prepare_geometry, build, variant, name, norm2, vector, ZZ, QQ, read, save
from shared import exact_minima as engine


def main():
    intervals = RealIntervalField(256)
    values = []
    for index, cfg in enumerate(GRID):
        geometry = prepare_geometry(index)
        certificate_path = HERE / "results" / name(cfg) / "B.certificate.json"
        if not certificate_path.exists():
            continue
        certificate = read(certificate_path)
        if not certificate.get("exact"):
            continue
        old_audit_path = HERE / "audits" / (name(cfg) + ".json")
        old_audit = next(item for item in read(old_audit_path) if item["variant"] == "B")
        assert old_audit["exact"] and old_audit["certificate_sha256"] == engine._hash_file(certificate_path)
        data = build(cfg)
        _, exponents, _, root = variant(data, "B")
        a = norm2(root)
        assert str(a) == geometry["root_sq"]
        witness = vector(ZZ, certificate["lambda2_witness"])
        S = ZZ(certificate["lambda2_sq"])
        assert norm2(witness) == S and ZZ(certificate["lambda1_sq"]) == a
        assert all(witness[j] % data["scales"][j] == 0 for j in range(len(witness)))
        raw = vector(ZZ, [witness[j] // data["scales"][j] for j in range(len(witness))])
        affine = raw * data["forward"].transpose()
        assert all(affine[j] % data["N"] ** exponents[j] == 0 for j in range(len(exponents)))
        dot = witness.dot_product(root)
        lower_sq = QQ(S) - QQ(a) / 4
        upper_sq = QQ(S) - QQ(dot * dot) / a
        assert 0 < lower_sq <= upper_sq
        lower_log = intervals(lower_sq).log2().lower() / 2
        upper_log = intervals(upper_sq).log2().upper() / 2
        lower_rounded = (lower_log * 10 ** 6 + QQ(1) / 2).floor()
        upper_rounded = (upper_log * 10 ** 6 + QQ(1) / 2).floor()
        assert lower_rounded == upper_rounded
        rounded = f"{int(lower_rounded) // 10 ** 6}.{int(lower_rounded) % 10 ** 6:06d}"
        result = dict(case=name(cfg), index=index, kind="full", exact_squared_norm=False,
                      certified_decimal_places=6, log2_rounded=rounded,
                      log2_lower_bound=str(lower_log), log2_upper_bound=str(upper_log),
                      lower_squared_norm_numerator=str(lower_sq.numerator()),
                      lower_squared_norm_denominator=str(lower_sq.denominator()),
                      upper_squared_norm_numerator=str(upper_sq.numerator()),
                      upper_squared_norm_denominator=str(upper_sq.denominator()),
                      source_certificate=str(certificate_path.relative_to(HERE)),
                      source_certificate_sha256=engine._hash_file(certificate_path),
                      source_independent_audit_sha256=engine._hash_file(old_audit_path),
                      original_lattice_witness_membership_verified=True,
                      argument="Nearest lift has root component <= norm(root)/2; lower=lambda2^2-root_sq/4; upper=norm(projection(lambda2_witness))^2",
                      directed_rounding_precision_bits=256)
        path = OUT / name(cfg) / "full_bounds.json"
        save(path, result)
        values.append(result)
        print("PROJECTED_INTERVAL_CERTIFIED", index, rounded, "width_log2", float(upper_log - lower_log), flush=True)
    save(OUT / "interval_certificates.json", dict(count=len(values), values=values))


if __name__ == "__main__":
    main()
