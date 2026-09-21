# -*- coding: utf-8 -*-
"""Stage 2B: power / minimum detectable effect for the MVMR *direct* effects.

Extends the Stage 2 univariable power analysis to the Stage 3 multivariable
model, using exactly the same estimators so the two stages are comparable.

Stage 2 convention (kept identical here):
    implied_sd(se, n, eaf) = se * sqrt(2 * n * eaf * (1 - eaf))
    scale                  = sd_exposure / sd_outcome
    se_standardized        = se_raw * scale
    MDE(target) solved exactly from the two-sided normal power function
"""
import json
from math import sqrt
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd

ROOT = Path(r"D:\SX\TMFI_CPI_MR")
STAGE1 = ROOT / "01_stage1_instruments"
STAGE2 = ROOT / "02_stage2_power"
OUT = ROOT / "03_stage3_mvmr"
STAGE2B = OUT / "stage2b_power"
STAGE2B.mkdir(parents=True, exist_ok=True)

ALPHA = 0.05
NORMAL = NormalDist()
Z_CRIT = NORMAL.inv_cdf(1 - ALPHA / 2)
EFFECT_GRID = [0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10]


def two_sided_power(effect, se):
    if se <= 0:
        return float("nan")
    mu = abs(effect) / se
    return NORMAL.cdf(-Z_CRIT - mu) + 1 - NORMAL.cdf(Z_CRIT - mu)


def exact_mde(se, target_power):
    lo, hi = 0.0, 10 * se
    for _ in range(100):
        mid = (lo + hi) / 2
        if two_sided_power(mid, se) < target_power:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def implied_sd(se, n, eaf):
    return se * sqrt(2 * n * eaf * (1 - eaf))


def booly(v):
    return str(v).strip().upper() in ("TRUE", "1", "YES")


def main():
    methods = pd.read_csv(OUT / "stage3_mvmr_methods.csv")
    joint = pd.read_csv(OUT / "cpi_full_mvmr_harmonised_usable.csv")
    stage1 = pd.read_csv(STAGE1 / "tmfi_gws_ldclumped_r2_0p001_kb10000.csv")
    # Stage 2 exposure SD, reused verbatim for continuity
    sd_tmfi_stage2 = 2.030884204160846

    specs = [
        ("CPI_full", "cpi_full_mvmr_harmonised_usable.csv", 2.8864741714064026),
        ("CPI_no_imaging", "cpi_no_imaging_mvmr_harmonised_usable.csv", 2.932677728436147),
    ]

    # Implied exposure SDs from the joint instrument set (TMFI- and BMI-GWS subsets)
    tmfi_gws = joint.loc[joint["TMFI_GWS"].map(booly)].copy()
    bmi_gws = joint.loc[joint["BMI_GWS"].map(booly)].copy()
    sd_tmfi_joint = float(np.median(
        tmfi_gws.apply(lambda r: implied_sd(r["TMFI_SE"], r["TMFI_N"], r["TMFI_EAF"]), axis=1)
    ))
    # GIANT coefficients are already per SD of BMI; the implied-SD
    # estimator returns 1.23 only because that release's standard errors
    # are ~26% inflated (single-trait LD score intercept 0.62).
    sd_bmi_joint = 1.0
    sd_bmi_joint_implied = float(np.median(   # audit trail only
        bmi_gws.apply(lambda r: implied_sd(r["BMI_SE"], r["BMI_N"], r["BMI_EAF_HARMONISED"]), axis=1)
    ))

    rows = []
    summaries = {}

    for outcome, fname, sd_y_default in specs:
        dat = pd.read_csv(OUT / fname)
        sd_y = float(np.median(
            dat.apply(lambda r: implied_sd(r["OUTCOME_SE"], r["OUTCOME_N"], r["OUTCOME_EAF_H"]), axis=1)
        ))

        for exposure, sd_exp in (("TMFI", sd_tmfi_stage2), ("BMI", sd_bmi_joint)):
            sub = methods.loc[
                (methods["outcome"] == outcome)
                & (methods["exposure"] == exposure)
                & (methods["method"] == "MVMR_IVW_multiplicative_random_effects")
            ]
            if sub.empty:
                continue
            beta_raw = float(sub["beta"].iloc[0])
            se_raw = float(sub["se"].iloc[0])
            scale = sd_exp / sd_y
            se_std = se_raw * scale
            beta_std = beta_raw * scale

            for e in EFFECT_GRID:
                rows.append({
                    "outcome": outcome,
                    "exposure": exposure,
                    "effect_sd_y_per_sd_exposure": e,
                    "power_two_sided_alpha_0_05": two_sided_power(e, se_std),
                })

            summaries[f"{outcome}|{exposure}"] = {
                "implied_outcome_sd": sd_y,
                "implied_exposure_sd": sd_exp,
                "mvmr_beta_raw": beta_raw,
                "mvmr_se_raw": se_raw,
                "mvmr_beta_standardized_sd_y_per_sd_exposure": beta_std,
                "mvmr_se_standardized": se_std,
                "ci95_standardized": [beta_std - 1.96 * se_std, beta_std + 1.96 * se_std],
                "mde_80_raw": exact_mde(se_raw, 0.80),
                "mde_90_raw": exact_mde(se_raw, 0.90),
                "mde_80_standardized": exact_mde(se_std, 0.80),
                "mde_90_standardized": exact_mde(se_std, 0.90),
            }

    power_curve = pd.DataFrame(rows)
    power_curve.to_csv(STAGE2B / "stage2b_mvmr_power_curve.csv", index=False)

    summary = {
        "parameters": {
            "alpha_two_sided": ALPHA,
            "estimator": "Stage 3 MVMR IVW multiplicative random effects",
            "estimand": "direct effect conditional on the other exposure",
            "power_method": "exact two-sided normal approximation using observed MVMR-MRE SE",
            "scaling_method": "implied SD from SE, N and EAF; scale = sd_exposure / sd_outcome",
        },
        "exposure_sd": {
            "TMFI_stage2_value": sd_tmfi_stage2,
            "TMFI_joint_set_check": sd_tmfi_joint,
            "BMI_joint_set": sd_bmi_joint,
            "BMI_joint_set_implied_not_used": sd_bmi_joint_implied,
        },
        "results": summaries,
    }
    (STAGE2B / "stage2b_mvmr_power_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print("\nPower curve:\n", power_curve.to_string(index=False))


if __name__ == "__main__":
    main()
