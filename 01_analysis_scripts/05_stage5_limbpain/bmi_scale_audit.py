#!/usr/bin/env python
"""
Audit the measurement scale of the BMI summary statistics used as an exposure.

Two independent tests:

  T1  Known-locus anchor.
      Per-allele effects at FTO / TMEM18 / MC4R / SEC16B are published in
      kg/m^2 and in SD units.  Observing which one the file matches pins the
      unit of `beta`.

  T2  Variance-explained anchor.
      For independent genome-wide-significant loci,
          V = sum 2 f (1-f) beta^2
      equals the *absolute* phenotypic variance they explain, i.e. the
      variance measured in the same units as `beta`.  GIANT (Locke 2015)
      reports that its 97 loci account for ~2.7% of BMI variance.  With
      Var(BMI) ~= 4.8^2 = 23.04 kg^2/m^4, the two hypotheses give
          per-SD :  V ~= 0.027
          kg/m^2 :  V ~= 0.62
      a 23-fold separation.

  T3  Internal consistency of the implied-SD estimator.
          implied_sd = SE * sqrt(2 N f (1-f))
      is the residual SD of the phenotype in the GWAS's own units.  If the
      GWAS standardised BMI, this must be ~1; if it used kg/m^2, ~4.8.
"""
import gzip
import math

BMI = r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain\gwas_harmonised\BMI_harmonised_patched.tsv.gz"

# published per-allele BMI effects (Locke et al. 2015 / GIANT)
KNOWN = {
    "rs1558902": (16, 53820527, "FTO",    0.39, 0.081),   # kg/m2, SD
    "rs2867125": (2,   622827, "TMEM18",  0.31, 0.065),
    "rs571312":  (18, 57839769, "MC4R",   0.22, 0.046),
    "rs543874":  (1,  177913272, "SEC16B", 0.22, 0.046),
    "rs9939609": (16, 53820527, "FTO2",   0.35, 0.073),
}

VAR_BMI_KGM2 = 4.8 ** 2      # SD of BMI in UKB-like samples, kg/m^2
GIANT_R2 = 0.027             # variance explained by the ~97 GIANT loci


def load():
    out = []
    with gzip.open(BMI, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {c: h.index(c) for c in
              ("rsid", "chrom", "position", "beta", "se", "af", "n", "z")}
        for line in fh:
            x = line.rstrip("\n").split("\t")
            try:
                b = float(x[ix["beta"]]); se = float(x[ix["se"]])
                f = float(x[ix["af"]]);   N = float(x[ix["n"]])
                z = float(x[ix["z"]])
            except ValueError:
                continue
            if not (0.005 < f < 0.995) or N <= 0 or se <= 0:
                continue
            out.append(dict(rsid=x[ix["rsid"]], chrom=x[ix["chrom"]],
                            pos=int(float(x[ix["position"]])),
                            beta=b, se=se, af=f, n=N, z=z))
    return out


def main():
    rows = load()
    print(f"usable variants: {len(rows):,}\n")

    # ---------------- T1 known loci ----------------
    print("T1  known-locus anchor")
    print(f"    {'locus':10s} {'beta_obs':>10s} {'kg/m2':>8s} {'SD':>8s}  verdict")
    by_rs = {r["rsid"]: r for r in rows}
    for rs, (ch, po, nm, kg, sd) in KNOWN.items():
        r = by_rs.get(rs)
        if r is None:
            print(f"    {nm:10s} {rs} not in file")
            continue
        b = abs(r["beta"])
        v = "per-SD" if abs(b - sd) < abs(b - kg) else "kg/m2"
        print(f"    {nm:10s} {b:10.4f} {kg:8.2f} {sd:8.3f}  closer to {v}")
    print()

    # ---------------- T2 variance explained ----------------
    gws = [r for r in rows if abs(r["z"]) > 5.45]
    gws.sort(key=lambda r: -abs(r["z"]))
    kept = []
    for r in gws:
        ok = True
        for k in kept:
            if k["chrom"] == r["chrom"] and abs(k["pos"] - r["pos"]) < 1_000_000:
                ok = False
                break
        if ok:
            kept.append(r)
    V = sum(2 * r["af"] * (1 - r["af"]) * r["beta"] ** 2 for r in kept)
    print("T2  variance-explained anchor")
    print(f"    independent GWS loci (|z|>5.45, 1 Mb clumping): {len(kept)}")
    print(f"    V = sum 2f(1-f)beta^2 = {V:.5f}")
    print(f"    hypothesis per-SD (Var=1)     : expected ~{GIANT_R2:.3f}"
          f"   -> obs/exp = {V/GIANT_R2:.2f}")
    exp_kgm2 = GIANT_R2 * VAR_BMI_KGM2
    print(f"    hypothesis kg/m^2 (Var=23.04) : expected ~{exp_kgm2:.3f}"
          f"   -> obs/exp = {V/exp_kgm2:.2f}")
    print()

    # ---------------- T3 implied SD ----------------
    imp_all = sorted(r["se"] * math.sqrt(2 * r["n"] * r["af"] * (1 - r["af"]))
                     for r in rows)
    imp_gws = sorted(r["se"] * math.sqrt(2 * r["n"] * r["af"] * (1 - r["af"]))
                     for r in kept)
    med = imp_all[len(imp_all) // 2]
    print("T3  implied residual SD = SE * sqrt(2 N f (1-f))")
    print(f"    all variants : median {med:.4f}  "
          f"[p25 {imp_all[len(imp_all)//4]:.4f}, p75 {imp_all[3*len(imp_all)//4]:.4f}]")
    print(f"    GWS loci     : median {imp_gws[len(imp_gws)//2]:.4f}")
    print(f"    per-SD expectation : ~1.00 (slightly <1 after covariates)")
    print(f"    kg/m^2 expectation : ~{4.8:.2f}")
    print()

    # ---------------- consequence ----------------
    print("Consequence of treating SD_BMI as the implied SD")
    for cand, lab in ((med, "implied (all)"), (1.0, "true per-SD")):
        print(f"    SD_BMI = {cand:.4f} ({lab})")
    print(f"    ratio implied/1.0 = {med:.4f};  an effect expressed per SD is"
          f" inflated by that factor if SD_BMI is set to the implied value"
          f" while beta is already per SD.")


if __name__ == "__main__":
    main()
