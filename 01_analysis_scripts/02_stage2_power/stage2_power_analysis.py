import json
from math import sqrt
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd


ROOT = Path(r"D:\SX\TMFI_CPI_MR")
STAGE1 = ROOT / "01_stage1_instruments"
OUT = ROOT / "02_stage2_power"
OUT.mkdir(parents=True, exist_ok=True)

N_COVARIATES = 24  # age, sex, education, batch, and 20 PCs
ALPHA = 0.05
NORMAL = NormalDist()
Z_CRIT = NORMAL.inv_cdf(1 - ALPHA / 2)


def two_sided_power(effect, se):
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
    return se * np.sqrt(2 * n * eaf * (1 - eaf))


def main():
    instruments = pd.read_csv(STAGE1 / "tmfi_gws_ldclumped_r2_0p001_kb10000.csv")
    instruments["RESIDUAL_DF"] = instruments["N"] - N_COVARIATES - 2
    instruments["R2_PARTIAL"] = instruments["F"] / (instruments["F"] + instruments["RESIDUAL_DF"])
    instruments["IMPLIED_TMFI_SD"] = implied_sd(instruments["SE"], instruments["N"], instruments["AF1"])
    instruments.to_csv(OUT / "stage2_instrument_r2.csv", index=False)

    total_r2 = float(instruments["R2_PARTIAL"].sum())
    k = len(instruments)
    n_exp = int(instruments["N"].median())
    set_f = ((n_exp - k - N_COVARIATES - 1) / k) * (total_r2 / (1 - total_r2))
    sd_exp = float(instruments["IMPLIED_TMFI_SD"].median())

    methods = pd.read_csv(STAGE1 / "formal_mr_methods.csv")
    rows = []
    outcome_summaries = {}
    specs = [
        ("CPI_full", "cpi_full_ldclumped_harmonised.csv"),
        ("CPI_no_imaging", "cpi_no_imaging_ldclumped_harmonised.csv"),
    ]
    standardized_effects = [0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10]

    for outcome, filename in specs:
        dat = pd.read_csv(STAGE1 / filename)
        sd_y_values = implied_sd(dat["SE_OUT"], dat["N_OUT"], dat["EAF_OUT_HARMONISED"])
        sd_y = float(np.median(sd_y_values))
        n_out = int(dat["N_OUT"].median())
        ivw = methods.loc[
            (methods["outcome"] == outcome)
            & (methods["method"] == "IVW_multiplicative_random_effects")
        ].iloc[0]
        se_raw = float(ivw["se"])
        beta_raw = float(ivw["beta"])
        scale = sd_exp / sd_y
        se_std = se_raw * scale
        beta_std = beta_raw * scale
        mde80_raw = exact_mde(se_raw, 0.80)
        mde90_raw = exact_mde(se_raw, 0.90)
        mde80_std = exact_mde(se_std, 0.80)
        mde90_std = exact_mde(se_std, 0.90)

        for effect_std in standardized_effects:
            rows.append(
                {
                    "outcome": outcome,
                    "effect_sd_y_per_sd_tmfi": effect_std,
                    "power_two_sided_alpha_0_05": two_sided_power(effect_std, se_std),
                    "se_standardized_mre": se_std,
                }
            )

        outcome_summaries[outcome] = {
            "median_n_outcome": n_out,
            "implied_outcome_sd": sd_y,
            "ivw_mre_beta_raw": beta_raw,
            "ivw_mre_se_raw": se_raw,
            "ivw_mre_beta_standardized": beta_std,
            "ivw_mre_se_standardized": se_std,
            "mde_80_raw_outcome_units_per_1_tmfi_unit": mde80_raw,
            "mde_90_raw_outcome_units_per_1_tmfi_unit": mde90_raw,
            "mde_80_sd_outcome_per_sd_tmfi": mde80_std,
            "mde_90_sd_outcome_per_sd_tmfi": mde90_std,
        }

    power = pd.DataFrame(rows)
    power.to_csv(OUT / "stage2_power_curve.csv", index=False)

    summary = {
        "parameters": {
            "alpha_two_sided": ALPHA,
            "exposure_covariates": N_COVARIATES,
            "power_method": "exact two-sided normal approximation using observed IVW-MRE SE",
            "r2_method": "sum of per-SNP partial R2 = F/(F+residual_df); valid as an approximation after LD clumping",
        },
        "exposure": {
            "instrument_count": k,
            "median_exposure_n": n_exp,
            "total_partial_r2": total_r2,
            "set_level_f": float(set_f),
            "median_implied_tmfi_sd": sd_exp,
            "implied_tmfi_sd_iqr": [
                float(instruments["IMPLIED_TMFI_SD"].quantile(0.25)),
                float(instruments["IMPLIED_TMFI_SD"].quantile(0.75)),
            ],
        },
        "outcomes": outcome_summaries,
        "unit_status": {
            "paper": "TMFI is the fat-free-muscle-volume-weighted average of four regional MRI-derived TMFI measures.",
            "most_likely_gwas_beta_unit": "raw TMFI units, likely percentage points; implied SD is about 2.03 units.",
            "conflict": "The article states continuous fastGWA-mlm, while the public GWAS.txt contains --fastGWA-mlm-binary. Exact preprocessing is not documented.",
            "reporting_rule": "Use standardized effects for primary power interpretation until authors clarify the GWAS phenotype coding/unit.",
        },
    }
    (OUT / "stage2_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print("\nPower curve:\n", power.to_string(index=False))


if __name__ == "__main__":
    main()
