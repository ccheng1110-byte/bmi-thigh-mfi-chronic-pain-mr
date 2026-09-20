#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Stage 5 - propagate the measured cross-trait LDSC intercepts into MR bias.

Identical algebra to the Stage 4 audit:
  cov(eps_Xj, eps_Yj) = I * se_Xj * se_Yj
  Stage 1 (IVW)   bias = I * sum(se_X / se_Y) / sum(beta_X^2 / se_Y^2)
  Stage 2 (MVMR)  bias = (X'WX)^{-1} u ,  u_k = I_k * sum_j (se_Xkj / se_Yj)

Structural ceiling: intercept = N_shared*rho/sqrt(N1*N2) with |rho| <= 1 and
N_shared <= min(N1, N2), so |intercept| <= sqrt(min(N1,N2)/max(N1,N2)).
"""
import json
import os

import numpy as np
import pandas as pd

S4 = r"D:\SX\TMFI_CPI_MR\04_stage4_sampling_cov\stage4_gwas"
S5 = r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain\gwas_harmonised"
RES = r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain\results"
INST1 = r"D:\SX\TMFI_CPI_MR\01_stage1_instruments\tmfi_gws_ldclumped_r2_0p001_kb10000.csv"
INST2 = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\mvmr_union_exposures_usable.csv"

SD_TMFI = 2.030884204160846
SD_BMI = 1.0  # GIANT coefficients are per SD of BMI (intercept 0.62 implies ~26% SE inflation, so the implied 1.23 is not the exposure SD)

OUT = {"BACK_PAIN_UKB": "BACK_PAIN_UKB_harmonised_patched.tsv.gz",
       "KNEE_PAIN_UKB": "KNEE_PAIN_UKB_harmonised_patched.tsv.gz",
       "KNEE_PAIN_FIGSHARE": "KNEE_PAIN_FIGSHARE_harmonised_patched.tsv.gz"}

# risk difference -> log odds ratio; identical constants to stage5_mr.py.
# Scaling beta_Y and se_Y by the same constant leaves bias_in_se, the LDSC
# intercept, the structural ceiling and the tipping point unchanged, because
# both bias and theta scale by lambda while D/S scales by 1/lambda.
LAMBDA = {"BACK_PAIN_UKB": 5.24349,
          "KNEE_PAIN_UKB": 5.95101,
          "KNEE_PAIN_FIGSHARE": 5.76064}

ld = json.load(open(os.path.join(
    r"D:\SX\TMFI_CPI_MR\04_stage4_sampling_cov\stage4_results",
    "stage5_cross_trait_ldsc.json")))


def k(c, p):
    return c.astype(np.int64) * 1_000_000_000 + p.astype(np.int64)


def load(path):
    return pd.read_csv(path, sep="\t")


tmfi = load(os.path.join(S5, "TMFI_harmonised_patched.tsv.gz"))
bmi = load(os.path.join(S5, "BMI_harmonised_patched.tsv.gz"))
i1 = pd.read_csv(INST1)
i1["key"] = k(i1["CHR"], i1["POS"])
i2 = pd.read_csv(INST2)
i2["key"] = k(i2["CHR"], i2["POS"])

rows = []
for oc, f in OUT.items():
    out = load(os.path.join(S5, f))
    lam = LAMBDA[oc]
    out["beta"] = out["beta"] * lam
    out["se"] = out["se"] * lam
    tm = out[["key", "beta", "se", "n"]].rename(
        columns={"beta": "out_beta", "se": "out_se", "n": "out_n"})
    x1 = (i1[["SNP", "key"]].merge(tmfi[["key", "beta", "se", "n"]], on="key",
                                   how="inner").merge(tm, on="key", how="inner"))
    x2 = (i2[["SNP", "key"]]
          .merge(tmfi[["key", "beta", "se", "n"]].rename(
              columns={"beta": "t_beta", "se": "t_se", "n": "t_n"}), on="key",
              how="inner")
          .merge(bmi[["key", "beta", "se", "n"]].rename(
              columns={"beta": "b_beta", "se": "b_se", "n": "b_n"}), on="key",
              how="inner")
          .merge(tm, on="key", how="inner"))

    pair = ld["pairs"][f"TMFI__{oc}"]
    i_tmfi = pair["cross"]["intercept"]
    i_tmfi_se = pair["cross"]["intercept_se"]
    n_tmfi = pair["mean_n"]["TMFI"]
    n_out = pair["mean_n"][oc]
    bp = ld["pairs"].get(f"BMI__{oc}")
    i_bmi = bp["cross"]["intercept"] if bp else 0.0
    n_bmi = bp["mean_n"]["BMI"] if bp else np.nan

    print(f"\n{'=' * 76}\n{oc}")
    print(f"  intercept TMFI x outcome = {i_tmfi:+.4f} (SE {i_tmfi_se:.4f}); "
          f"BMI x outcome = {i_bmi:+.4f}")

    # ---------------- Stage 1 -------------------------------------------
    bx = x1["beta"].to_numpy(float)
    sx = x1["se"].to_numpy(float)
    by = x1["out_beta"].to_numpy(float)
    sy = x1["out_se"].to_numpy(float)
    w = 1 / sy ** 2
    theta = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    se_fe = np.sqrt(1 / np.sum(w * bx ** 2))
    q = np.sum(w * (by - theta * bx) ** 2)
    phi = max(1.0, q / (len(bx) - 1))
    se = se_fe * np.sqrt(phi)
    S = np.sum(sx / sy)
    D = np.sum(w * bx ** 2)
    bias = i_tmfi * S / D
    corr = theta - bias
    ceil = np.sqrt(min(n_tmfi, n_out) / max(n_tmfi, n_out))
    tip = min(abs((theta - s * 1.96 * se) * D / S) for s in (-1, 1))
    print(f"  Stage 1: theta={theta:+.6f} (SE {se:.6f})  overlap bias="
          f"{bias:+.6f} ({bias / se:+.3f} SE, {bias * SD_TMFI:+.2e} logOR/SD)")
    print(f"           corrected theta={corr:+.6f} -> OR/SD="
          f"{np.exp(corr * SD_TMFI):.4f}")
    print(f"           ceiling |intercept| <= {ceil:.3f}; intercept needed to "
          f"reach p=0.05: {tip:.3f} "
          f"({'IMPOSSIBLE' if tip > ceil else f'{tip/ceil*100:.0f}% of ceiling'})")
    rows.append({"outcome": oc, "stage": "1 UVMR", "intercept": i_tmfi,
                 "theta": theta, "se": se, "bias": bias,
                 "bias_in_se": bias / se, "theta_corrected": corr,
                 "or_per_sd_corrected": float(np.exp(corr * SD_TMFI)),
                 "ceiling": ceil, "intercept_for_p05": tip})

    # ---------------- Stage 2 -------------------------------------------
    X = np.column_stack([x2["t_beta"].to_numpy(float),
                         x2["b_beta"].to_numpy(float)])
    seX = np.column_stack([x2["t_se"].to_numpy(float),
                           x2["b_se"].to_numpy(float)])
    by2 = x2["out_beta"].to_numpy(float)
    sy2 = x2["out_se"].to_numpy(float)
    w2 = 1 / sy2 ** 2
    A = X.T @ (X * w2[:, None])
    Ainv = np.linalg.inv(A)
    th2 = Ainv @ (X.T @ (by2 * w2))
    res = by2 - X @ th2
    q2 = np.sum(w2 * res ** 2)
    phi2 = max(1.0, q2 / (len(by2) - 2))
    se2 = np.sqrt(np.diag(Ainv) * phi2)
    u = np.array([i_tmfi * np.sum(seX[:, 0] / sy2),
                  i_bmi * np.sum(seX[:, 1] / sy2)])
    bias2 = Ainv @ u
    corr2 = th2 - bias2
    ceil_t = np.sqrt(min(n_tmfi, n_out) / max(n_tmfi, n_out))
    ceil_b = np.sqrt(min(n_bmi, n_out) / max(n_bmi, n_out)) if n_bmi == n_bmi else np.nan
    print(f"  Stage 2: TMFI theta={th2[0]:+.6f} (SE {se2[0]:.6f}) bias="
          f"{bias2[0]:+.6f} ({bias2[0] / se2[0]:+.3f} SE) -> corrected "
          f"{corr2[0]:+.6f} OR/SD={np.exp(corr2[0] * SD_TMFI):.4f}")
    print(f"           BMI  theta={th2[1]:+.6f} (SE {se2[1]:.6f}) bias="
          f"{bias2[1]:+.6f} ({bias2[1] / se2[1]:+.3f} SE) -> corrected "
          f"{corr2[1]:+.6f} OR/SD={np.exp(corr2[1] * SD_BMI):.4f}")
    print(f"           ceilings: TMFI {ceil_t:.3f}, BMI {ceil_b:.3f}")
    for j, nm, sd, cl in ((0, "TMFI", SD_TMFI, ceil_t), (1, "BMI", SD_BMI, ceil_b)):
        rows.append({"outcome": oc, "stage": "2 MVMR",
                     "exposure": nm, "intercept": (i_tmfi, i_bmi)[j],
                     "theta": float(th2[j]), "se": float(se2[j]),
                     "bias": float(bias2[j]), "bias_in_se": float(bias2[j] / se2[j]),
                     "theta_corrected": float(corr2[j]),
                     "or_per_sd_corrected": float(np.exp(corr2[j] * sd)),
                     "ceiling": float(cl)})

pd.DataFrame(rows).to_csv(os.path.join(RES, "stage5_overlap_bias.csv"), index=False)
print(f"\n[out] {os.path.join(RES, 'stage5_overlap_bias.csv')}")
