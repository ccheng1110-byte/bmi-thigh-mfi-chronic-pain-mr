#!/usr/bin/env python
"""
Why does the data-derived implied SD of BMI come out at ~1.23 when the betas
are demonstrably per-SD (variance-explained test)?

implied_sd = SE * sqrt(2 * N * f * (1-f))

Three candidate explanations are separated here:

  E1  the allele frequency used is wrong  -> recompute with the reference-panel
      frequency (column `af1`) instead of the file's own frequency (`af`)
  E2  the per-SNP N is wrong              -> implied_sd / sqrt(N) should then be
      tightly constant; check the spread of implied_sd across N strata
  E3  the reported SEs are inflated by a constant k -> detect with the LDSC
      intercept: for correctly calibrated SEs the intercept of the regression
      of z^2 on the LD score is ~1; uniform SE inflation by k pushes it to 1/k^2.
"""
import gzip
import math
import statistics

FILES = {
    "BMI": r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain\gwas_harmonised\BMI_harmonised_patched.tsv.gz",
    "TMFI": r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain\gwas_harmonised\TMFI_harmonised_patched.tsv.gz",
    "CPI_full": r"D:\SX\TMFI_CPI_MR\04_stage4_sampling_cov\stage4_gwas\CPI_full_harmonised.tsv.gz",
}


def load(path, need_ldsc=True):
    rows = []
    with gzip.open(path, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {c: h.index(c) for c in h}
        for line in fh:
            x = line.rstrip("\n").split("\t")
            try:
                d = dict(beta=float(x[ix["beta"]]), se=float(x[ix["se"]]),
                         z=float(x[ix["z"]]), n=float(x[ix["n"]]),
                         af=float(x[ix["af"]]))
                d["af1"] = float(x[ix["af1"]]) if "af1" in ix else float("nan")
                d["ldsc"] = float(x[ix["ldscore"]]) if "ldscore" in ix else float("nan")
            except (ValueError, KeyError):
                continue
            if d["se"] <= 0 or d["n"] <= 0:
                continue
            rows.append(d)
    return rows


def med(v):
    v = sorted(v)
    return v[len(v) // 2]


def implied(d, afkey="af"):
    f = d[afkey]
    if not (0.005 < f < 0.995):
        return None
    return d["se"] * math.sqrt(2 * d["n"] * f * (1 - f))


def ldsc_intercept(rows):
    """Weighted regression of z^2 on the LD score; intercept ~1 if calibrated."""
    pts = [(r["ldsc"], r["z"] ** 2) for r in rows
           if r["ldsc"] == r["ldsc"] and r["ldsc"] > 1 and 0.01 <= r["af"] <= 0.99
           and abs(r["z"]) < 40]
    if len(pts) < 1000:
        return None
    # step 1: OLS
    n = len(pts)
    mx = sum(p[0] for p in pts) / n
    my = sum(p[1] for p in pts) / n
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    sxy = sum((p[0] - mx) * (p[1] - my) for p in pts)
    b = sxy / sxx
    a = my - b * mx
    # step 2: WLS with chi-square variance weights, two iterations
    for _ in range(3):
        num = den = sn = snx = snxx = snxy = sny = 0.0
        for l, z2 in pts:
            pred = a + b * l
            w = 1.0 / (2 * max(pred, 0.2) ** 2)
            sn += w; snx += w * l; snxx += w * l * l
            sny += w * z2; snxy += w * l * z2
        det = sn * snxx - snx * snx
        b = (sn * snxy - snx * sny) / det
        a = (snxx * sny - snx * snxy) / det
    return a, b, n


def main():
    for name, path in FILES.items():
        try:
            rows = load(path)
        except FileNotFoundError:
            print(f"{name}: file not found")
            continue
        print("=" * 72)
        print(f"{name}  ({len(rows):,} variants)")
        ia = [v for v in (implied(r, "af") for r in rows) if v]
        ib = [v for v in (implied(r, "af1") for r in rows) if v and v == v]
        print(f"  implied_sd using file AF   : median {med(ia):.4f}")
        if ib:
            print(f"  implied_sd using panel AF  : median {med(ib):.4f}")
            daf = sorted(abs(r["af"] - r["af1"]) for r in rows
                         if r["af1"] == r["af1"] and 0.005 < r["af"] < 0.995)
            daf2 = sorted(abs(r["af"] - (1 - r["af1"])) for r in rows
                          if r["af1"] == r["af1"] and 0.005 < r["af"] < 0.995)
            print(f"  |af - af1| median {med(daf):.4f}   "
                  f"|af - (1-af1)| median {med(daf2):.4f}")
        # implied_sd by N stratum -> is N the culprit?
        qs = sorted(r["n"] for r in rows)
        lo, hi = qs[int(len(qs) * 0.1)], qs[int(len(qs) * 0.9)]
        for lab, sel in (("low N", lambda r: r["n"] <= lo),
                         ("high N", lambda r: r["n"] >= hi)):
            v = [x for x in (implied(r) for r in rows if sel(r)) if x]
            if v:
                print(f"  implied_sd {lab:7s}(N<= {lo:.0f} / >= {hi:.0f}): "
                      f"median {med(v):.4f}")
        res = ldsc_intercept(rows)
        if res:
            a, b, n = res
            print(f"  LDSC intercept = {a:.4f}  (slope {b:.4g}, {n:,} SNPs)")
            print(f"  -> implied SE inflation k = 1/sqrt(intercept) = "
                  f"{1/math.sqrt(a):.4f}" if a > 0 else "")
            print(f"  -> implied_sd / k = {med(ia)/(1/math.sqrt(a)):.4f}"
                  if a > 0 else "")
        print()


if __name__ == "__main__":
    main()
