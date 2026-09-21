#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Stage 5 - Site-specific (lower-limb / axial) pain extension
==========================================================
Outcomes:  BACK_PAIN_UKB (ukb-b-9838, n=461,857, 118,471 cases)
           KNEE_PAIN_UKB (ukb-b-16254, n=461,857,  98,704 cases)
           KNEE_PAIN_FIGSHARE (independent knee-pain GWAS, n=429,737)

Stage 1 : univariable MR, the frozen 42-SNP TMFI instrument set
Stage 2 : multivariable MR, the frozen 93-SNP TMFI + BMI set (BMI = internal
          positive control)

All exposure and outcome files are already harmonised against the same
1000 Genomes EUR reference panel, so merging on the variant key gives
consistent allele orientation on both sides.

Binary outcomes -> estimates are log-odds; reported as OR per 1 SD of exposure.

NOTE ON STEIGER: the Python implementation uses the large-sample
R^2 = z^2 / (z^2 + N) approximation and a Fisher-z comparison; the
authoritative Steiger test is produced by the companion R script
(stage5_validation.R, TwoSampleMR::directionality_test).
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

S4 = r"D:\SX\TMFI_CPI_MR\04_stage4_sampling_cov\stage4_gwas"
S5 = r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain\gwas_harmonised"
RIN = r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain\r_input"
OUT = r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain\results"
os.makedirs(RIN, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

INST1 = r"D:\SX\TMFI_CPI_MR\01_stage1_instruments\tmfi_gws_ldclumped_r2_0p001_kb10000.csv"
INST2 = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\mvmr_union_exposures_usable.csv"

SD_TMFI = 2.030884204160846
SD_BMI = 1.2307496342161395

# patched files additionally contain the palindromic instruments, recovered by
# allele-frequency matching (see patch_palindromic.py)
OUTCOMES = {
    "BACK_PAIN_UKB": os.path.join(S5, "BACK_PAIN_UKB_harmonised_patched.tsv.gz"),
    "KNEE_PAIN_UKB": os.path.join(S5, "KNEE_PAIN_UKB_harmonised_patched.tsv.gz"),
    "KNEE_PAIN_FIGSHARE": os.path.join(S5, "KNEE_PAIN_FIGSHARE_harmonised_patched.tsv.gz"),
}


def k(c, p):
    return c.astype(np.int64) * 1_000_000_000 + p.astype(np.int64)


def load_harm(path, tag):
    d = pd.read_csv(path, sep="\t")
    ren = {"beta": f"{tag}_beta", "se": f"{tag}_se", "n": f"{tag}_n",
           "af": f"{tag}_af"}
    if "rsid" in d.columns:
        ren["rsid"] = f"{tag}_rsid"
    d = d.rename(columns=ren)
    keep = [c for c in ["key", "chrom", "position", f"{tag}_rsid",
                        f"{tag}_beta", f"{tag}_se", f"{tag}_n", f"{tag}_af"]
            if c in d.columns]
    return d[keep]


# --------------------------------------------------------------------------
# estimators
# --------------------------------------------------------------------------
def ivw(bx, by, sey):
    w = 1.0 / sey ** 2
    den = np.sum(w * bx ** 2)
    theta = np.sum(w * bx * by) / den
    se = np.sqrt(1.0 / den)
    q = np.sum(w * (by - theta * bx) ** 2)
    df = len(bx) - 1
    phi = max(1.0, q / df)
    return theta, se, se * np.sqrt(phi), q, df, phi


def weighted_median(bx, by, sey, n_boot=2000, seed=7):
    ratio = by / bx
    w = bx ** 2 / sey ** 2
    o = np.argsort(ratio)
    r, ww = ratio[o], w[o]
    c = (np.cumsum(ww) - 0.5 * ww) / np.sum(ww)
    est = np.interp(0.5, c, r)
    rng = np.random.default_rng(seed)
    boot = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, len(r), len(r))
        rr, www = r[idx], ww[idx]
        oo = np.argsort(rr)
        rr, www = rr[oo], www[oo]
        cc = (np.cumsum(www) - 0.5 * www) / np.sum(www)
        boot[i] = np.interp(0.5, cc, rr)
    return float(est), float(np.std(boot, ddof=1))


def egger(bx, by, sey):
    X = np.column_stack([np.ones_like(bx), bx])
    w = 1.0 / sey ** 2
    Xw = X * w[:, None]
    A = X.T @ Xw
    Ainv = np.linalg.inv(A)
    th = Ainv @ (Xw.T @ by)
    res = by - X @ th
    q = np.sum(w * res ** 2)
    df = len(bx) - 2
    phi = max(1.0, q / df)
    se = np.sqrt(np.diag(Ainv) * phi)
    return float(th[1]), float(th[0]), float(se[1]), float(se[0]), float(q), float(phi)


def steiger(bx, sex, nx, by, sey, ny):
    """R^2 = z^2/(z^2+N); Fisher-z comparison of the two correlations."""
    zx = np.abs(bx / sex)
    zy = np.abs(by / sey)
    r2x = zx ** 2 / (zx ** 2 + nx)
    r2y = zy ** 2 / (zy ** 2 + ny)
    rx = np.sqrt(np.clip(r2x, 0, 0.999999))
    ry = np.sqrt(np.clip(r2y, 0, 0.999999))
    fx = np.arctanh(rx)
    fy = np.arctanh(ry)
    se = np.sqrt(1.0 / np.maximum(nx - 3, 1) + 1.0 / np.maximum(ny - 3, 1))
    z = (fx - fy) / se
    zc = float(z.mean() * np.sqrt(len(z)))
    return {"mean_r2_exposure": float(np.mean(r2x)),
            "mean_r2_outcome": float(np.mean(r2y)),
            "combined_z": zc,
            "p_correct_direction": float(2 * stats.norm.sf(zc)) if zc > 0 else 1.0,
            "frac_snp_exposure_dominant": float((z > 0).mean())}


def mvmr(X, by, sey, seX):
    """MVMR IVW with multiplicative random effects; conditional F statistics."""
    k_, p_ = X.shape
    w = 1.0 / sey ** 2
    Xw = X * w[:, None]
    A = X.T @ Xw
    Ainv = np.linalg.inv(A)
    theta = Ainv @ (Xw.T @ by)
    res = by - X @ theta
    q = float(np.sum(w * res ** 2))
    df = k_ - p_
    phi = max(1.0, q / df)
    cov = Ainv * phi
    se = np.sqrt(np.diag(cov))
    # ---- strength: univariate mean F and two candidate conditional F -----
    prec = 1.0 / seX ** 2
    F_uni = np.mean(X ** 2 * prec, axis=0)
    F_cond_A = np.empty(p_)
    F_cond_B = np.empty(p_)
    for j in range(p_):
        others = [i for i in range(p_) if i != j]
        # (A) R^2-adjusted: meanF * (1 - R^2) with R^2 from unweighted OLS
        r2 = np.corrcoef(X[:, j], X[:, others[0]])[0, 1] ** 2 if p_ == 2 else np.nan
        F_cond_A[j] = F_uni[j] * (1 - r2)
        # (B) weighted partial F: residualise beta_j on the other exposures
        Xo = np.column_stack([np.ones(k_)] + [X[:, i] for i in others])
        wj = prec[:, j]
        Ao = Xo.T @ (Xo * wj[:, None])
        bo = Xo.T @ (X[:, j] * wj)
        try:
            proj = np.linalg.solve(Ao, bo)
            rss = X[:, j] - Xo @ proj
            F_cond_B[j] = float(np.sum(rss ** 2 * wj) / (k_ - p_))
        except np.linalg.LinAlgError:
            F_cond_B[j] = np.nan
    return {"theta": theta, "se": se, "q": q, "df": df, "phi": phi,
            "F_univariate": F_uni, "F_conditional_A": F_cond_A,
            "F_conditional_B": F_cond_B}


def or_per_sd(beta, se, sd):
    b = beta * sd
    s = se * sd
    return (float(np.exp(b)),
            float(np.exp(b - 1.96 * s)), float(np.exp(b + 1.96 * s)),
            float(2 * stats.norm.sf(abs(beta / se))),
            float(np.exp((1.96 + 0.8416) * s)))


def main():
    global SEX1, SEX2
    tmfi = load_harm(os.path.join(S5, "TMFI_harmonised_patched.tsv.gz"), "tmfi")
    bmi = load_harm(os.path.join(S5, "BMI_harmonised_patched.tsv.gz"), "bmi")

    i1 = pd.read_csv(INST1)
    i1["key"] = k(i1["CHR"], i1["POS"])
    i2 = pd.read_csv(INST2)
    i2["key"] = k(i2["CHR"], i2["POS"])
    print(f"Stage 1 instruments: {len(i1)}   Stage 2 instruments: {len(i2)}")

    rows = []
    for oc, path in OUTCOMES.items():
        out = load_harm(path, "out")
        print(f"\n{'=' * 78}\n{oc}   (harmonised variants {len(out):,})")

        # ---------------- Stage 1 : univariable, 42 SNPs -------------------
        d1 = (i1[["SNP", "key"]]
              .merge(tmfi, on="key", how="inner")
              .merge(out, on="key", how="inner"))
        print(f"  Stage 1: {len(d1)}/{len(i1)} instruments present in the outcome")
        bx = d1["tmfi_beta"].to_numpy(float)
        sx = d1["tmfi_se"].to_numpy(float)
        nx = d1["tmfi_n"].to_numpy(float)
        by = d1["out_beta"].to_numpy(float)
        sy = d1["out_se"].to_numpy(float)
        ny = d1["out_n"].to_numpy(float)
        # orientation sanity: |beta| must equal the frozen instrument file
        chk = i1.set_index("key")["BETA"]
        dev = np.abs(np.abs(bx) - np.abs(d1["key"].map(chk).to_numpy(float))).max()
        print(f"  |beta| vs frozen instrument file: max abs deviation = {dev:.2e}")

        th, se_fe, se_mre, q, df, phi = ivw(bx, by, sy)
        orr = or_per_sd(th, se_mre, SD_TMFI)
        wm, wm_se = weighted_median(bx, by, sy)
        eb, ea, eb_se, ea_se, qe, phie = egger(bx, by, sy)
        st = steiger(bx, sx, nx, by, sy, ny)
        print(f"  IVW-FE  beta={th:+.5f} (SE {se_fe:.5f})   "
              f"IVW-MRE beta={th:+.5f} (SE {se_mre:.5f})")
        print(f"  Q={q:.1f} (df {df}) phi={phi:.3f}  I2={max(0,(q-df)/q)*100:.0f}%")
        print(f"  OR per SD = {orr[0]:.4f} ({orr[1]:.4f}-{orr[2]:.4f})  P={orr[3]:.3g}"
              f"   MDE(80% power) OR/SD = {orr[4]:.4f}")
        print(f"  weighted median = {wm * SD_TMFI:+.4f} (logOR/SD, boot SE "
              f"{wm_se * SD_TMFI:.4f})")
        print(f"  MR-Egger slope={eb:+.5f} (SE {eb_se:.5f})  "
              f"intercept={ea:+.5f} (SE {ea_se:.5f}, "
              f"P={2 * stats.norm.sf(abs(ea / ea_se)):.3g})")
        print(f"  Steiger: mean R2 exposure={st['mean_r2_exposure']:.2e} vs "
              f"outcome={st['mean_r2_outcome']:.2e}; combined z={st['combined_z']:+.2f}"
              f" (P={st['p_correct_direction']:.2g}); "
              f"{st['frac_snp_exposure_dominant'] * 100:.0f}% of SNPs exposure-dominant")
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
                     "or_per_sd": float(np.exp(wm * SD_TMFI))})
        rows.append({"stage": "1 UVMR", "outcome": oc, "method": "MR-Egger",
                     "n_snp": len(d1), "beta": eb, "se": eb_se,
                     "or_per_sd": float(np.exp(eb * SD_TMFI)),
                     "p": float(2 * stats.norm.sf(abs(eb / eb_se))),
                     "egger_intercept": ea, "egger_intercept_p":
                         float(2 * stats.norm.sf(abs(ea / ea_se)))})
        rows.append({"stage": "1 UVMR", "outcome": oc, "method": "Steiger",
                     "n_snp": len(d1), "beta": st["combined_z"],
                     "p": st["p_correct_direction"],
                     "mean_r2_exposure": st["mean_r2_exposure"],
                     "mean_r2_outcome": st["mean_r2_outcome"]})
        d1.to_csv(os.path.join(RIN, f"R_stage1_{oc}.tsv"), sep="\t", index=False)

        # ---------------- Stage 2 : multivariable, 93 SNPs -----------------
        d2 = (i2[["SNP", "key"]]
              .merge(tmfi, on="key", how="inner")
              .merge(bmi, on="key", how="inner")
              .merge(out, on="key", how="inner"))
        print(f"  Stage 2: {len(d2)}/{len(i2)} instruments present in the outcome")
        X = np.column_stack([d2["tmfi_beta"].to_numpy(float),
                             d2["bmi_beta"].to_numpy(float)])
        SEX1 = d2["tmfi_se"].to_numpy(float)
        SEX2 = d2["bmi_se"].to_numpy(float)
        by2 = d2["out_beta"].to_numpy(float)
        sy2 = d2["out_se"].to_numpy(float)
        seX = np.column_stack([SEX1, SEX2])
        m = mvmr(X, by2, sy2, seX)
        for nm, j, sd in (("TMFI", 0, SD_TMFI), ("BMI", 1, SD_BMI)):
            orr = or_per_sd(m["theta"][j], m["se"][j], sd)
            print(f"  MVMR {nm}: beta={m['theta'][j]:+.5f} (SE {m['se'][j]:.5f})  "
                  f"OR/SD={orr[0]:.4f} ({orr[1]:.4f}-{orr[2]:.4f}) P={orr[3]:.3g}  "
                  f"MDE OR/SD={orr[4]:.4f}")
            rows.append({"stage": "2 MVMR", "outcome": oc, "method": f"IVW-MRE:{nm}",
                         "n_snp": len(d2), "beta": float(m["theta"][j]),
                         "se": float(m["se"][j]), "or_per_sd": orr[0],
                         "ci_low": orr[1], "ci_high": orr[2], "p": orr[3],
                         "mde_or_per_sd": orr[4], "Q": m["q"], "Q_df": m["df"],
                         "phi": m["phi"], "F_univariate": float(m["F_univariate"][j]),
                         "F_conditional_A": float(m["F_conditional_A"][j]),
                         "F_conditional_B": float(m["F_conditional_B"][j])})
        print(f"  Q={m['q']:.1f} (df {m['df']}) phi={m['phi']:.3f}")
        print(f"  univariate mean F: TMFI {m['F_univariate'][0]:.1f}, "
              f"BMI {m['F_univariate'][1]:.1f}")
        print(f"  conditional F (A, R2-adjusted): TMFI {m['F_conditional_A'][0]:.2f}, "
              f"BMI {m['F_conditional_A'][1]:.2f}")
        print(f"  conditional F (B, weighted partial): TMFI "
              f"{m['F_conditional_B'][0]:.2f}, BMI {m['F_conditional_B'][1]:.2f}")

        # ---- univariable BMI -> outcome (93-SNP set) as a scale reference --
        bt, st_fe, st_mre, qt, dft, phit = ivw(X[:, 1], by2, sy2)
        orb = or_per_sd(bt, st_mre, SD_BMI)
        print(f"  UVMR BMI (93 SNPs): beta={bt:+.5f} (SE {st_mre:.5f})  "
              f"OR/SD={orb[0]:.4f} ({orb[1]:.4f}-{orb[2]:.4f}) P={orb[3]:.3g}")
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
            print(f"  LOO on BMI coefficient: min {loo[0][1]:+.5f} "
                  f"(drop {loo[0][0]}), max {loo[-1][1]:+.5f} (drop {loo[-1][0]}), "
                  f"full {m['theta'][1]:+.5f}")
        d2.to_csv(os.path.join(RIN, f"R_stage2_{oc}.tsv"), sep="\t", index=False)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "stage5_mr_results.csv"), index=False)
    print(f"\n[out] {os.path.join(OUT, 'stage5_mr_results.csv')}")


if __name__ == "__main__":
    main()
