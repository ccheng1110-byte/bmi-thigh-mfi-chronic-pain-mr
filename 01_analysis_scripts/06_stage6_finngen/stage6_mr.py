# -*- coding: utf-8 -*-
"""Stage 6 - FinnGen R12 external replication (outcome-side independence).

Why this stage exists
---------------------
Every outcome used so far (CPI_full, CPI_no_imaging, knee pain, back pain,
figshare knee pain) is drawn from UK Biobank, which is also the source of the
TMFI and BMI exposures.  Reviewers have repeatedly read the resulting overlap
as an unresolvable limitation.  FinnGen R12 shares no samples with UKB, so a
concordant result here turns the overlap argument from theory into evidence.

What is held constant
---------------------
* exposures: TMFI and BMI, the same frozen 42-SNP and 93-SNP instrument sets
* estimators: imported from `stage5_mr` so there is exactly one implementation
* SD_TMFI / SD_BMI: unchanged

What deliberately differs
-------------------------
* NO lambda.  Stage 5 outcomes were BOLT-LMM risk differences and had to be
  divided by mu(1-mu) to reach log-odds.  FinnGen fits REGENIE logistic models;
  `harmonise_finngen.py` verified SE^2*2Nf(1-f) = 1/(mu(1-mu)) to within
  1.7-4.5%, so the betas are already log odds ratios.  Applying lambda here
  would be the same double-conversion bug in the opposite direction.
* outcomes are binary ICD-based diagnoses, not the continuous CPI score, so
  this is a related-phenotype replication and is reported as such.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

S5 = r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain"
sys.path.insert(0, S5)
# one implementation of every estimator, shared with Stage 5
from stage5_mr import (ivw, weighted_median, egger, steiger, mvmr,  # noqa: E402
                       or_per_sd, load_harm, k)

S5H = os.path.join(S5, "gwas_harmonised")          # exposures: TMFI, BMI (UKB)
HERE = os.path.dirname(os.path.abspath(__file__))
S6H = os.path.join(HERE, "gwas_harmonised")        # outcomes: FinnGen R12
RIN = os.path.join(HERE, "r_input")
OUT = os.path.join(HERE, "results")
os.makedirs(RIN, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

INST1 = r"D:\SX\TMFI_CPI_MR\01_stage1_instruments\tmfi_gws_ldclumped_r2_0p001_kb10000.csv"
INST2 = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\mvmr_union_exposures_usable.csv"

SD_TMFI = 2.030884204160846
SD_BMI = 1.0

OUTCOMES = {
    "FG_PAIN": os.path.join(S6H, "PAIN_harmonised.tsv.gz"),
    "FG_LIMBPAIN": os.path.join(S6H, "M13_LIMBPAIN_harmonised.tsv.gz"),
    "FG_FIBROMYALGIA": os.path.join(S6H, "M13_FIBROMYALGIA_harmonised.tsv.gz"),
}


def main():
    tmfi = load_harm(os.path.join(S5H, "TMFI_harmonised_patched.tsv.gz"), "tmfi")
    bmi = load_harm(os.path.join(S5H, "BMI_harmonised_patched.tsv.gz"), "bmi")

    i1 = pd.read_csv(INST1)
    i1["key"] = k(i1["CHR"], i1["POS"])
    i2 = pd.read_csv(INST2)
    i2["key"] = k(i2["CHR"], i2["POS"])
    print("Stage 1 instruments: %d   Stage 2 instruments: %d" % (len(i1), len(i2)))

    rows = []
    for oc, path in OUTCOMES.items():
        out = load_harm(path, "out")
        print("\n" + "=" * 78)
        print("%s   (harmonised variants %d)   scale: log-OR, no lambda"
              % (oc, len(out)))

        # ---------------- Stage 1 : univariable, 42 SNPs -------------------
        d1 = (i1[["SNP", "key"]]
              .merge(tmfi, on="key", how="inner")
              .merge(out, on="key", how="inner"))
        print("  Stage 1: %d/%d instruments present in the outcome" % (len(d1), len(i1)))
        bx = d1["tmfi_beta"].to_numpy(float)
        sx = d1["tmfi_se"].to_numpy(float)
        nx = d1["tmfi_n"].to_numpy(float)
        by = d1["out_beta"].to_numpy(float)
        sy = d1["out_se"].to_numpy(float)
        ny = d1["out_n"].to_numpy(float)
        # orientation sanity: |beta| must match the frozen instrument file
        chk = i1.set_index("key")["BETA"]
        dev = np.abs(np.abs(bx) - np.abs(d1["key"].map(chk).to_numpy(float))).max()
        print("  |beta| vs frozen instrument file: max abs deviation = %.2e" % dev)

        th, se_fe, se_mre, q, df, phi = ivw(bx, by, sy)
        orr = or_per_sd(th, se_mre, SD_TMFI)
        wm, wm_se = weighted_median(bx, by, sy)
        eb, ea, eb_se, ea_se, qe, phie = egger(bx, by, sy)
        st = steiger(bx, sx, nx, by, sy, ny)
        print("  IVW-FE  beta=%+.5f (SE %.5f)   IVW-MRE beta=%+.5f (SE %.5f)"
              % (th, se_fe, th, se_mre))
        print("  Q=%.1f (df %d) phi=%.3f  I2=%.0f%%"
              % (q, df, phi, max(0, (q - df) / q) * 100))
        print("  OR per SD = %.4f (%.4f-%.4f)  P=%.3g   MDE(80%% power) OR/SD = %.4f"
              % (orr[0], orr[1], orr[2], orr[3], orr[4]))
        print("  weighted median = %+.4f (logOR/SD, boot SE %.4f)"
              % (wm * SD_TMFI, wm_se * SD_TMFI))
        print("  MR-Egger slope=%+.5f (SE %.5f)  intercept=%+.5f (SE %.5f, P=%.3g)"
              % (eb, eb_se, ea, ea_se, 2 * stats.norm.sf(abs(ea / ea_se))))
        print("  Steiger: mean R2 exposure=%.2e vs outcome=%.2e; combined z=%+.2f "
              "(P=%.2g); %.0f%% of SNPs exposure-dominant"
              % (st["mean_r2_exposure"], st["mean_r2_outcome"], st["combined_z"],
                 st["p_correct_direction"], st["frac_snp_exposure_dominant"] * 100))
        rows.append({"stage": "1 UVMR", "outcome": oc, "method": "IVW-MRE",
                     "n_snp": len(d1), "beta": th, "se": se_mre,
                     "or_per_sd": orr[0], "ci_low": orr[1], "ci_high": orr[2],
                     "p": orr[3], "mde_or_per_sd": orr[4],
                     "Q": q, "Q_df": df, "phi": phi})
        rows.append({"stage": "1 UVMR", "outcome": oc, "method": "IVW-FE",
                     "n_snp": len(d1), "beta": th, "se": se_fe,
                     "or_per_sd": float(np.exp(th * SD_TMFI)),
                     "ci_low": float(np.exp((th - 1.96 * se_fe) * SD_TMFI)),
                     "ci_high": float(np.exp((th + 1.96 * se_fe) * SD_TMFI)),
                     "p": float(2 * stats.norm.sf(abs(th / se_fe)))})
        rows.append({"stage": "1 UVMR", "outcome": oc, "method": "weighted median",
                     "n_snp": len(d1), "beta": wm, "se": wm_se,
                     "or_per_sd": float(np.exp(wm * SD_TMFI)),
                     "p": float(2 * stats.norm.sf(abs(wm / wm_se)))})
        rows.append({"stage": "1 UVMR", "outcome": oc, "method": "MR-Egger",
                     "n_snp": len(d1), "beta": eb, "se": eb_se,
                     "or_per_sd": float(np.exp(eb * SD_TMFI)),
                     "p": float(2 * stats.norm.sf(abs(eb / eb_se))),
                     "egger_intercept": ea,
                     "egger_intercept_p": float(2 * stats.norm.sf(abs(ea / ea_se)))})
        rows.append({"stage": "1 UVMR", "outcome": oc, "method": "Steiger",
                     "n_snp": len(d1), "beta": st["combined_z"],
                     "p": st["p_correct_direction"],
                     "mean_r2_exposure": st["mean_r2_exposure"],
                     "mean_r2_outcome": st["mean_r2_outcome"]})
        # inputs for the independent base-R re-implementation: already log-OR,
        # so nothing is rescaled here or there
        d1.to_csv(os.path.join(RIN, "R_stage1_%s.tsv" % oc), sep="\t", index=False)

        # ---------------- Stage 2 : multivariable, 93 SNPs -----------------
        d2 = (i2[["SNP", "key"]]
              .merge(tmfi, on="key", how="inner")
              .merge(bmi, on="key", how="inner")
              .merge(out, on="key", how="inner"))
        print("  Stage 2: %d/%d instruments present in the outcome" % (len(d2), len(i2)))
        X = np.column_stack([d2["tmfi_beta"].to_numpy(float),
                             d2["bmi_beta"].to_numpy(float)])
        seX = np.column_stack([d2["tmfi_se"].to_numpy(float),
                               d2["bmi_se"].to_numpy(float)])
        by2 = d2["out_beta"].to_numpy(float)
        sy2 = d2["out_se"].to_numpy(float)
        m = mvmr(X, by2, sy2, seX)
        for nm, j, sd in (("TMFI", 0, SD_TMFI), ("BMI", 1, SD_BMI)):
            orr = or_per_sd(m["theta"][j], m["se"][j], sd)
            print("  MVMR %s: beta=%+.5f (SE %.5f)  OR/SD=%.4f (%.4f-%.4f) "
                  "P=%.3g  MDE OR/SD=%.4f"
                  % (nm, m["theta"][j], m["se"][j], orr[0], orr[1], orr[2],
                     orr[3], orr[4]))
            rows.append({"stage": "2 MVMR", "outcome": oc, "method": "IVW-MRE:%s" % nm,
                         "n_snp": len(d2), "beta": float(m["theta"][j]),
                         "se": float(m["se"][j]), "or_per_sd": orr[0],
                         "ci_low": orr[1], "ci_high": orr[2], "p": orr[3],
                         "mde_or_per_sd": orr[4], "Q": m["q"], "Q_df": m["df"],
                         "phi": m["phi"], "F_univariate": float(m["F_univariate"][j]),
                         "F_conditional_A": float(m["F_conditional_A"][j]),
                         "F_conditional_B": float(m["F_conditional_B"][j])})
        print("  Q=%.1f (df %d) phi=%.3f" % (m["q"], m["df"], m["phi"]))
        print("  univariate mean F: TMFI %.1f, BMI %.1f"
              % (m["F_univariate"][0], m["F_univariate"][1]))
        print("  conditional F (A, R2-adjusted): TMFI %.2f, BMI %.2f"
              % (m["F_conditional_A"][0], m["F_conditional_A"][1]))
        print("  conditional F (B, weighted partial): TMFI %.2f, BMI %.2f"
              % (m["F_conditional_B"][0], m["F_conditional_B"][1]))

        # ---- univariable BMI -> outcome (93-SNP set) as a scale reference --
        bt, st_fe, st_mre, qt, dft, phit = ivw(X[:, 1], by2, sy2)
        orb = or_per_sd(bt, st_mre, SD_BMI)
        print("  UVMR BMI (93 SNPs): beta=%+.5f (SE %.5f)  OR/SD=%.4f "
              "(%.4f-%.4f) P=%.3g"
              % (bt, st_mre, orb[0], orb[1], orb[2], orb[3]))
        rows.append({"stage": "2 MVMR", "outcome": oc, "method": "UVMR BMI",
                     "n_snp": len(d2), "beta": float(bt), "se": float(st_mre),
                     "or_per_sd": orb[0], "ci_low": orb[1], "ci_high": orb[2],
                     "p": orb[3]})

        # ---- leave-one-out influence on the MVMR BMI coefficient ----------
        loo = []
        for i in range(len(X)):
            keep = np.ones(len(X), bool)
            keep[i] = False
            try:
                mm = mvmr(X[keep], by2[keep], sy2[keep], seX[keep])
                loo.append((d2["out_rsid"].iloc[i] if "out_rsid" in d2.columns
                            else str(i), mm["theta"][1]))
            except np.linalg.LinAlgError:
                pass
        if loo:
            loo.sort(key=lambda t: t[1])
            print("  LOO on BMI coefficient: min %+.5f (drop %s), max %+.5f "
                  "(drop %s), full %+.5f"
                  % (loo[0][1], loo[0][0], loo[-1][1], loo[-1][0], m["theta"][1]))
        d2.to_csv(os.path.join(RIN, "R_stage2_%s.tsv" % oc), sep="\t", index=False)

        # ---- which instruments are absent, and why ------------------------
        have = set(d2["out_rsid"]) if "out_rsid" in d2.columns else set()
        miss = [s for s in i2["SNP"] if s not in have]
        if miss:
            print("  instruments absent from this outcome: %s" % ", ".join(miss))

    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, "stage6_mr_results.csv"), index=False)
    print("\n[out] %s" % os.path.join(OUT, "stage6_mr_results.csv"))


if __name__ == "__main__":
    main()
