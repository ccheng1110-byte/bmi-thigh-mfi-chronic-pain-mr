#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Stage 4 - Sampling covariance audit
Step 4: translate cross-trait LDSC intercepts into the bias they induce on the
MR estimates reported in the manuscript, and remove it.

Why this works
--------------
Cross-trait LD Score regression estimates

    intercept_km = N_shared * rho_phen / sqrt(N_k N_m)

Because z = beta/se, the component of the covariance between two summary
statistics at SNP j that is attributable to shared individuals is

    cov(epsilon_Xj, epsilon_Yj) = intercept_KM * se_Xj * se_Yj          (1)

Every IVW-type MR estimator is a ratio of two weighted sums; the numerator is
X'W * beta_Y.  Taking expectations and using (1),

    E[X'W beta_Y] = (X'WX) theta  +  u,
    u_i           = intercept_i * sum_j (se_Xij / se_Yj)                (2)

so the entire effect of sample overlap is

    bias(theta_hat) = (X'WX)^{-1} u                                     (3)

Setting K = 1 recovers the familiar univariable expression

    bias = intercept * sum_j (se_Xj/se_Yj) / sum_j (beta_Xj^2/se_Yj^2)  (4)

Because the bias enters through the *noise* covariance while the signal enters
through beta_X^2, strong instruments mechanically dilute it: equations (3)-(4)
shrink like 1/mean(z_X^2).  This is why the quantity a reviewer actually cares
about is not the raw overlap but the overlap propagated through the reported
instrument set.

Genomic-control rescaling
-------------------------
If a contributing GWAS was genomic-control corrected by a factor lambda before
release, all of its z statistics were divided by sqrt(lambda), and its own
single-trait LDSC intercept estimates 1/lambda rather than 1.  A product of two
such GWAS therefore returns intercept_obs = intercept_true / sqrt(c1*c2),
where c_k is the single-trait intercept of trait k.  We therefore also report
GC-rescaled intercepts, intercept_true = intercept_obs / sqrt(c1*c2), which is
the conservative correction whenever a contributing intercept sits below 1.

This is the analytic core of MRlap (Mounier & Kutalik 2023, eLife), applied
here transparently to the exact instrument sets reported in the manuscript so
that every step can be audited line by line.
"""
import json
import math
import os

import numpy as np
import pandas as pd
from scipy import stats

BASE = r"D:\SX\TMFI_CPI_MR"
OUT_DIR = os.path.join(BASE, "04_stage4_sampling_cov", "stage4_results")
S1 = os.path.join(BASE, "01_stage1_instruments")
S3 = os.path.join(BASE, "03_stage3_mvmr")

os.makedirs(OUT_DIR, exist_ok=True)

OUTCOMES = {
    "CPI_full": {
        "cross_key": "TMFI__CPI_full",
        "bmi_cross_key": "BMI__CPI_full",
        "stage1_outcome_file": "cpi_full_ldclumped_harmonised.csv",
        "stage2_file": "cpi_full_mvmr_harmonised_usable.csv",
    },
    "CPI_no_imaging": {
        "cross_key": "TMFI__CPI_no_imaging",
        "bmi_cross_key": "BMI__CPI_no_imaging",
        "stage1_outcome_file": "cpi_no_imaging_ldclumped_harmonised.csv",
        "stage2_file": "cpi_no_imaging_mvmr_harmonised_usable.csv",
    },
}

# Standardisation scales frozen in the Stage 2 / Stage 3 analyses
SD_TMFI = 2.030884204160846
SD_BMI = 1.0  # per SD of BMI
SD_OUTCOME = {"CPI_full": 2.886503, "CPI_no_imaging": 2.931411}

# Standard errors as actually reported in Table 3 (multiplicative random
# effects, i.e. inflated by the dispersion statistic phi).  The IVW algebra
# reproduced below returns fixed-effect standard errors; quoting both keeps
# the tipping-point analysis conservative, because the larger Table 3 standard
# error demands a far larger intercept before any conclusion moves.
REPORTED_SE_RAW = {
    "CPI_full": {"TMFI": 0.02441, "BMI": 0.07745},
    "CPI_no_imaging": {"TMFI": 0.02733, "BMI": 0.08666},
}


def ivw(bx, by, sey):
    w = 1.0 / sey ** 2
    beta = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    se = np.sqrt(1.0 / np.sum(w * bx ** 2))
    return beta, se


def mvmr_ivw(Xv, yv, sey):
    w = 1.0 / sey ** 2
    Xw = Xv * w[:, None]
    A = Xv.T @ Xw
    theta = np.linalg.solve(A, Xw.T @ yv)
    cov = np.linalg.inv(A)
    return theta, np.sqrt(np.diag(cov)), cov


def main():
    res = json.load(open(os.path.join(OUT_DIR, "stage4_cross_trait_ldsc.json")))
    P = res["pairs"]
    summary = []

    def gc_rescale(key_outcome, outcome_name, tmfi_key, bmi_key):
        """observed -> GC-corrected intercepts for the TMFI and BMI paths."""
        c_tmfi = P[tmfi_key]["trait1"]["intercept"]       # TMFI single-trait
        c_out = P[tmfi_key]["trait2"]["intercept"]        # outcome single-trait
        c_bmi = P[bmi_key]["trait1"]["intercept"]         # BMI single-trait
        c_out2 = P[bmi_key]["trait2"]["intercept"]        # outcome single-trait
        i_tmfi_obs = P[tmfi_key]["cross"]["intercept"]
        i_bmi_obs = P[bmi_key]["cross"]["intercept"]
        s_tmfi = np.sqrt(c_tmfi * c_out)
        s_bmi = np.sqrt(c_bmi * c_out2)
        return {
            "i_tmfi_obs": i_tmfi_obs, "i_tmfi_se": P[tmfi_key]["cross"]["intercept_se"],
            "i_bmi_obs": i_bmi_obs, "i_bmi_se": P[bmi_key]["cross"]["intercept_se"],
            "i_tmfi_gc": i_tmfi_obs / s_tmfi, "i_bmi_gc": i_bmi_obs / s_bmi,
            "c_tmfi": c_tmfi, "c_out": c_out, "c_bmi": c_bmi, "c_out2": c_out2,
        }

    print("=" * 78)
    print("STAGE 4 - sampling covariance audit: from intercept to MR bias")
    print("=" * 78)

    tmfi = pd.read_csv(os.path.join(S1, "tmfi_gws_ldclumped_r2_0p001_kb10000.csv"))
    tmfi = tmfi.rename(columns={"SNP": "snp"})[["snp", "BETA", "SE"]]

    # ------------------------------------------------------------------
    print("\n### 1. Univariable MR (Stage 1, 42 instruments)")
    print("-" * 78)
    for oc, cfg in OUTCOMES.items():
        cross = res["pairs"][cfg["cross_key"]]
        g = gc_rescale(cfg["cross_key"], oc, cfg["cross_key"], cfg["bmi_cross_key"])
        intercept = g["i_tmfi_gc"]
        ise = cross["cross"]["intercept_se"] / np.sqrt(g["c_tmfi"] * g["c_out"])
        N1 = cross["mean_n"]["TMFI"]
        N2 = cross["mean_n"][oc]
        out = pd.read_csv(os.path.join(S1, cfg["stage1_outcome_file"]))
        out = out[["SNP", "BETA_OUT_HARMONISED", "SE_OUT"]] \
            if "BETA_OUT_HARMONISED" in out.columns else out
        d = tmfi.merge(out, left_on="snp", right_on="SNP", how="inner")
        assert len(d) == 42, f"expected 42 instruments, got {len(d)}"

        bx = d["BETA"].to_numpy(float)
        se_x = d["SE"].to_numpy(float)
        by = d["BETA_OUT_HARMONISED"].to_numpy(float)
        se_y = d["SE_OUT"].to_numpy(float)

        beta_obs, se_ivw = ivw(bx, by, se_y)
        # eq (4)
        bias = intercept * np.sum(se_x / se_y) / np.sum(bx ** 2 / se_y ** 2)
        beta_corr = beta_obs - bias
        scale1 = SD_TMFI / SD_OUTCOME[oc]

        print(f"\n[{oc}]  TMFI-outcome cross intercept = {g['i_tmfi_obs']:.4f} obs "
              f"-> {intercept:.4f} GC-rescaled (SE {ise:.4f})")
        print(f"   single-trait intercepts: TMFI {g['c_tmfi']:.4f}, {oc} {g['c_out']:.4f} "
              f"-> rescale factor {np.sqrt(g['c_tmfi']*g['c_out']):.4f}")
        print(f"   N_TMFI~{N1:,.0f}, N_{oc}~{N2:,.0f};  "
              f"N_shared*rho_phen implied = {intercept*np.sqrt(N1*N2):,.0f}")
        print(f"   IVW as reported : beta = {beta_obs:+.5f}  SE = {se_ivw:.5f}  "
              f"({beta_obs - 1.96*se_ivw:+.4f} to {beta_obs + 1.96*se_ivw:+.4f} raw | "
              f"{beta_obs*scale1 - 1.96*se_ivw*scale1:+.4f} to "
              f"{beta_obs*scale1 + 1.96*se_ivw*scale1:+.4f} SD/SD)")
        print(f"   overlap bias    : {bias:+.5f} raw / {bias*scale1:+.5f} SD/SD "
              f"({bias/se_ivw:+.2f} reported SE)")
        print(f"   overlap-corrected: beta = {beta_corr:+.5f} raw / "
              f"{beta_corr*scale1:+.5f} SD/SD  "
              f"({beta_corr - 1.96*se_ivw:+.4f} to {beta_corr + 1.96*se_ivw:+.4f} raw | "
              f"{beta_corr*scale1 - 1.96*se_ivw*scale1:+.4f} to "
              f"{beta_corr*scale1 + 1.96*se_ivw*scale1:+.4f} SD/SD)")
        summary.append({
            "stage": "Stage 1 univariable IVW", "outcome": oc, "scenario": "measured",
            "intercept_tmfi": intercept, "intercept_tmfi_se": ise,
            "intercept_bmi": g["i_bmi_gc"],
            "n_shared_times_rho": intercept * np.sqrt(N1 * N2),
            "beta_observed_raw": beta_obs, "se_raw": se_ivw,
            "bias_raw": bias, "bias_std": bias * scale1,
            "bias_in_se": bias / se_ivw,
            "beta_corrected_raw": beta_corr,
            "beta_obs_std": beta_obs * scale1, "se_std": se_ivw * scale1,
            "beta_corr_std": beta_corr * scale1,
        })

    # ------------------------------------------------------------------
    print("\n### 2. Multivariable MR (Stage 2, 93 instruments)")
    print("-" * 78)
    print("Bias follows eq (3): bias = (X'WX)^-1 u, "
          "u_i = intercept_i * sum_j (se_Xij/se_Yj).")
    for oc, cfg in OUTCOMES.items():
        g = gc_rescale(cfg["cross_key"], oc, cfg["cross_key"], cfg["bmi_cross_key"])
        i_tmfi = g["i_tmfi_gc"]
        i_bmi_obs = g["i_bmi_obs"]
        i_bmi_se = g["i_bmi_se"]
        i_bmi_gc = g["i_bmi_gc"]
        d = pd.read_csv(os.path.join(S3, cfg["stage2_file"]))
        bx1 = d["TMFI_BETA"].to_numpy(float)
        se1 = d["TMFI_SE"].to_numpy(float)
        bx2 = d["BMI_BETA_HARMONISED"].to_numpy(float)
        se2 = d["BMI_SE"].to_numpy(float)
        by = d["OUTCOME_BETA_H"].to_numpy(float)
        sey = d["OUTCOME_SE"].to_numpy(float)
        Xv = np.column_stack([bx1, bx2])

        theta_r, se_th, Ainv = mvmr_ivw(Xv, by, sey)
        scale = np.array([SD_TMFI / SD_OUTCOME[oc], SD_BMI / SD_OUTCOME[oc]])
        theta_s = theta_r * scale
        se_s = se_th * scale

        # eq (2): each exposure's sensitivity to its own overlap
        S1v = float(np.sum(se1 / sey))
        S2v = float(np.sum(se2 / sey))
        Sv = np.array([S1v, S2v])

        print(f"\n[{oc}] instruments = {len(d)}")
        print(f"   reproduced IVW (must match Table 3):  "
              f"TMFI {theta_r[0]:+.5f} -> {theta_s[0]:+.5f} SD/SD (SE {se_s[0]:.5f});  "
              f"BMI {theta_r[1]:+.5f} -> {theta_s[1]:+.5f} SD/SD (SE {se_s[1]:.5f})")
        print(f"   sum_j(se_TMFI/se_Y) = {S1v:.2f};  "
              f"sum_j(se_BMI/se_Y) = {S2v:.2f}")
        print(f"   median |z| of instruments: TMFI "
              f"{np.median(np.abs(bx1/se1)):.1f}, BMI {np.median(np.abs(bx2/se2)):.1f}")

        def bias_from(i1, i2):
            """eq (3)."""
            u = np.array([i1, i2]) * Sv
            return Ainv @ u

        def z_attaining(i1, i2):
            corr = theta_r - bias_from(i1, i2)
            return corr / se_th

        def tipping(target, other_i, se_vec):
            """Smallest |intercept| on exposure `target`, holding the other
            exposure at other_i, that moves the corrected estimate of `target`
            to the p = 0.05 boundary.  The bias is linear in the intercept,
            bias(i) = b0 + i * Sv[t] * Ainv[:, t], so the root is exact."""
            u0 = np.zeros(2)
            u0[1 - target] = other_i * Sv[1 - target]
            b0 = Ainv @ u0
            slope = Sv[target] * Ainv[target, target]
            if slope == 0:
                return np.nan
            base = theta_r[target] - b0[target]
            cands = [(base - sgn * 1.96 * se_vec[target]) / slope
                     for sgn in (-1, 1)]
            return float(min(cands, key=abs))

        # conservative upper bound for the BMI path from its own 95% CI
        i_bmi_ub = (i_bmi_obs + 1.96 * i_bmi_se) / np.sqrt(g["c_bmi"] * g["c_out2"])
        scenarios = [
            ("A. measured intercepts (TMFI and BMI)", i_tmfi, i_bmi_gc),
            ("B. BMI at the upper 95% CI bound of its intercept", i_tmfi, i_bmi_ub),
            ("C. BMI overlapping the outcome as much as TMFI does (>10x measured)",
             i_tmfi, i_tmfi),
            ("D. BMI intercept 3x the TMFI intercept (extreme stress test)",
             i_tmfi, 3 * i_tmfi),
        ]
        for tag, i1, i2 in scenarios:
            bias_r = bias_from(i1, i2)
            bias_s = bias_r * scale
            corr_r = theta_r - bias_r
            corr_s = corr_r * scale
            zt, zb = z_attaining(i1, i2)
            pt = 2 * stats.norm.sf(abs(zt))
            pb = 2 * stats.norm.sf(abs(zb))
            print(f"\n   -- {tag} --")
            print(f"      intercept_TMFI={i1:+.4f}   intercept_BMI={i2:+.4f}")
            print(f"      TMFI bias {bias_r[0]:+.5f} raw / {bias_s[0]:+.5f} SD/SD "
                  f"({bias_r[0]/se_th[0]:+.2f} SE)  ->  corrected "
                  f"{corr_r[0]:+.5f} raw / {corr_s[0]:+.5f} SD/SD, "
                  f"95%CI {corr_s[0]-1.96*se_s[0]:+.4f} to {corr_s[0]+1.96*se_s[0]:+.4f}"
                  f", p={pt:.3g}")
            print(f"      BMI  bias {bias_r[1]:+.5f} raw / {bias_s[1]:+.5f} SD/SD "
                  f"({bias_r[1]/se_th[1]:+.2f} SE)  ->  corrected "
                  f"{corr_r[1]:+.5f} raw / {corr_s[1]:+.5f} SD/SD, "
                  f"95%CI {corr_s[1]-1.96*se_s[1]:+.4f} to {corr_s[1]+1.96*se_s[1]:+.4f}"
                  f", p={pb:.3g}")
            summary.append({
                "stage": "Stage 2 multivariable IVW", "outcome": oc, "scenario": tag,
                "intercept_tmfi": i1, "intercept_bmi": i2,
                "tmfi_bias_std": bias_s[0], "tmfi_bias_in_se": bias_r[0] / se_th[0],
                "tmfi_corrected_std": corr_s[0], "tmfi_se_std": se_s[0],
                "tmfi_p_after": pt,
                "bmi_bias_std": bias_s[1], "bmi_bias_in_se": bias_r[1] / se_th[1],
                "bmi_corrected_std": corr_s[1], "bmi_se_std": se_s[1],
                "bmi_p_after": pb,
            })

        # ---- tipping point: how big must each intercept be to change the call --
        # Two standard errors are carried through: the fixed-effect SE returned
        # by the IVW algebra, and the larger multiplicative-random-effects SE
        # actually quoted in Table 3 (inflated by the dispersion statistic phi).
        # The larger SE demands a larger intercept, so quoting both keeps the
        # analysis conservative.
        se_fe = se_th
        se_rep = np.array([REPORTED_SE_RAW[oc]["TMFI"],
                           REPORTED_SE_RAW[oc]["BMI"]])
        tip_t = tipping(0, i_bmi_gc, se_fe)
        tip_t_rep = tipping(0, i_bmi_gc, se_rep)
        tip_b0 = tipping(1, 0.0, se_fe)
        tip_b0_rep = tipping(1, 0.0, se_rep)
        print(f"\n   tipping point for TMFI (BMI held at its measured intercept): "
              f"intercept_TMFI must reach {tip_t:+.4f} "
              f"({abs(tip_t/i_tmfi):.1f}x its measured value) before the corrected "
              f"TMFI effect attains p = 0.05")
        print(f"      using the larger Table 3 standard error instead: "
              f"{tip_t_rep:+.4f} ({abs(tip_t_rep/i_tmfi):.1f}x measured)")
        print(f"   tipping point for BMI with no TMFI overlap at all: "
              f"intercept_BMI must reach {tip_b0:+.4f} "
              f"({abs(tip_b0/i_bmi_gc):.1f}x its measured value) before the corrected "
              f"BMI effect loses significance at p = 0.05")
        print(f"      using the larger Table 3 standard error instead: "
              f"{tip_b0_rep:+.4f}")
        pair = res["pairs"][cfg["cross_key"]]
        bmi_pair = res["pairs"][cfg["bmi_cross_key"]]
        n_tmfi = pair["mean_n"]["TMFI"]
        n_out = pair["mean_n"][oc]
        n_bmi = bmi_pair["mean_n"]["BMI"]

        def ceiling(nx, ny):
            """Largest cross-trait intercept attainable by sample overlap alone:
            intercept = N_shared*rho/sqrt(Nx*Ny) with |rho|<=1 and
            N_shared <= min(Nx,Ny)  =>  ceiling = sqrt(min/max)."""
            return float(np.sqrt(min(nx, ny) / max(nx, ny)))

        ceil_t = ceiling(n_tmfi, n_out)
        ceil_b = ceiling(n_bmi, n_out)
        print(f"\n   structural ceilings (|rho_phen|<=1, N_shared<=min(N1,N2)):")
        for nm, cval, tp, tpr, meas in (("TMFI", ceil_t, tip_t, tip_t_rep, i_tmfi),
                                        ("BMI ", ceil_b, tip_b0, tip_b0_rep, i_bmi_gc)):
            verdict = ("IMPOSSIBLE - the intercept would have to exceed the "
                       "maximum any amount of overlap can generate"
                       if abs(tp) > cval else
                       f"requires {100*abs(tp)/cval:.1f}% of the theoretical "
                       f"maximum overlap with |rho_phen|=1")
            verdict_r = ("IMPOSSIBLE" if abs(tpr) > cval else "reachable in principle")
            print(f"      {nm}: ceiling {cval:+.4f} | tipping(FE SE) {tp:+.4f} -> "
                  f"{verdict}")
            print(f"           {'':>{len(nm)}}                    | "
                  f"tipping(Table 3 SE) {tpr:+.4f} -> {verdict_r}")
        summary.append({
            "stage": "Stage 2 structural ceiling", "outcome": oc,
            "scenario": "max attainable intercept vs tipping point",
            "intercept_tmfi": ceil_t, "intercept_bmi": ceil_b,
            "tmfi_corrected_std": tip_t, "bmi_corrected_std": tip_b0,
            "tmfi_se_std": tip_t_rep, "bmi_se_std": tip_b0_rep,
        })

    df = pd.DataFrame(summary)
    df.to_csv(os.path.join(OUT_DIR, "stage4_bias_correction.csv"), index=False)
    print(f"\n[out] {os.path.join(OUT_DIR, 'stage4_bias_correction.csv')}")


if __name__ == "__main__":
    main()
