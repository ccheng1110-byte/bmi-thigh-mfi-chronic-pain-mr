"""Stage 3B robustness analyses for the two-exposure MVMR.

This script intentionally uses only the already harmonised Stage 3 tables and
NumPy/Pandas. It does not run R. The robust fit is labelled exploratory
Huber-IRLS (it is not a substitute for a package-specific publication method).
"""

import json
from math import erfc, sqrt
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"D:\SX\TMFI_CPI_MR")
OUT = ROOT / "03_stage3_mvmr"
P = 2


def normal_p(z):
    return float(erfc(abs(float(z)) / sqrt(2.0)))


def fit_ivw(dat):
    x = dat[["TMFI_BETA", "BMI_BETA_HARMONISED"]].to_numpy(float)
    y = dat["OUTCOME_BETA_H"].to_numpy(float)
    sy = dat["OUTCOME_SE"].to_numpy(float)
    w = 1.0 / sy**2
    xtw = x.T * w
    bread = np.linalg.inv(xtw @ x)
    beta = bread @ (xtw @ y)
    residual = y - x @ beta
    q = float(np.sum(w * residual**2))
    df = len(dat) - P
    phi = max(1.0, q / df)
    se_fixed = np.sqrt(np.diag(bread))
    se_mre = se_fixed * sqrt(phi)
    return {
        "beta": beta,
        "se_fixed": se_fixed,
        "se_mre": se_mre,
        "residual": residual,
        "w": w,
        "q": q,
        "df": df,
        "phi": phi,
    }


def fit_huber(dat, c=1.345, max_iter=100):
    x = dat[["TMFI_BETA", "BMI_BETA_HARMONISED"]].to_numpy(float)
    y = dat["OUTCOME_BETA_H"].to_numpy(float)
    sy = dat["OUTCOME_SE"].to_numpy(float)
    beta = fit_ivw(dat)["beta"].copy()
    robust_weights = np.ones(len(dat))
    for _ in range(max_iter):
        residual = y - x @ beta
        scale = np.median(np.abs(residual / sy)) / 0.6745
        scale = max(float(scale), 1e-8)
        u = residual / (sy * scale)
        robust_weights = np.minimum(1.0, c / np.maximum(np.abs(u), 1e-12))
        w = robust_weights / sy**2
        new_beta = np.linalg.solve((x.T * w) @ x, (x.T * w) @ y)
        if np.max(np.abs(new_beta - beta)) < 1e-10:
            beta = new_beta
            break
        beta = new_beta
    residual = y - x @ beta
    bread = np.linalg.inv((x.T * w) @ x)
    q = float(np.sum(w * residual**2))
    phi = max(1.0, q / (len(dat) - P))
    se = np.sqrt(np.diag(bread) * phi)
    return {"beta": beta, "se": se, "residual": residual, "weights": robust_weights, "phi": phi}


def result_rows(outcome, method, fit, n):
    if "se_mre" in fit:
        se = fit["se_mre"]
    else:
        se = fit["se"]
    rows = []
    for exposure, idx in [("TMFI", 0), ("BMI", 1)]:
        z = fit["beta"][idx] / se[idx]
        rows.append(
            {
                "outcome": outcome,
                "method": method,
                "exposure": exposure,
                "nsnp": int(n),
                "beta": float(fit["beta"][idx]),
                "se": float(se[idx]),
                "p_normal": normal_p(z),
                "ci_low": float(fit["beta"][idx] - 1.96 * se[idx]),
                "ci_high": float(fit["beta"][idx] + 1.96 * se[idx]),
            }
        )
    return rows


def influence_table(dat, outcome):
    fit = fit_ivw(dat)
    x = dat[["TMFI_BETA", "BMI_BETA_HARMONISED"]].to_numpy(float)
    xw = x * np.sqrt(fit["w"])[:, None]
    h = np.sum((xw @ np.linalg.inv(xw.T @ xw)) * xw, axis=1)
    h = np.clip(h, 0, 0.999999)
    std_resid = np.sqrt(fit["w"]) * fit["residual"] / np.sqrt(fit["phi"] * (1 - h))
    cook = std_resid**2 * h / (P * (1 - h))
    out = dat[["SNP", "TMFI_F", "BMI_F", "TMFI_GWS", "BMI_GWS"]].copy()
    out["outcome"] = outcome
    out["leverage"] = h
    out["std_residual"] = std_resid
    out["cook_distance"] = cook
    out["abs_std_residual"] = np.abs(std_resid)
    out["cook_gt_4_over_n"] = cook > 4.0 / len(dat)
    out["abs_stdresid_gt_3"] = np.abs(std_resid) > 3.0
    return out.sort_values(["cook_distance", "abs_std_residual"], ascending=False).reset_index(drop=True)


def loo_table(dat, outcome):
    rows = []
    for snp in dat["SNP"].astype(str):
        subset = dat.loc[dat["SNP"].astype(str) != snp].copy()
        fit = fit_ivw(subset)
        rows.extend(result_rows(outcome, f"leave_one_out:{snp}", fit, len(subset)))
    return pd.DataFrame(rows)


def stratified(dat, outcome):
    d = dat.copy()
    d["tool_group"] = np.select(
        [d["TMFI_GWS"] & d["BMI_GWS"], d["TMFI_GWS"], d["BMI_GWS"]],
        ["both_GWS", "TMFI_only_GWS", "BMI_only_GWS"],
        default="neither_GWS",
    )
    d["strength_group"] = np.select(
        [(d["TMFI_F"] >= 10) & (d["BMI_F"] >= 10), d["TMFI_F"] >= 10, d["BMI_F"] >= 10],
        ["both_F_ge_10", "TMFI_F_ge_10_only", "BMI_F_ge_10_only"],
        default="neither_F_ge_10",
    )
    rows = []
    for grouping in ["tool_group", "strength_group"]:
        for label, subset in d.groupby(grouping, sort=True):
            if len(subset) < 5:
                continue
            fit = fit_ivw(subset)
            rows.extend(result_rows(outcome, f"stratified:{grouping}={label}", fit, len(subset)))
            rows[-2]["group_n"] = int(len(subset))
            rows[-1]["group_n"] = int(len(subset))
    return pd.DataFrame(rows), d[["SNP", "tool_group", "strength_group"]]


def run_outcome(filename, outcome):
    dat = pd.read_csv(OUT / filename)
    influence = influence_table(dat, outcome)
    loo = loo_table(dat, outcome)
    strat, groups = stratified(dat, outcome)
    primary = fit_ivw(dat)
    robust = fit_huber(dat)

    result = result_rows(outcome, "primary_IVW_MRE", primary, len(dat))
    result.extend(result_rows(outcome, "exploratory_Huber_IRLS", robust, len(dat)))
    for method, mask in [
        ("exclude_Cook_gt_4_over_n", influence["cook_gt_4_over_n"]),
        ("exclude_abs_stdresid_gt_3", influence["abs_stdresid_gt_3"]),
    ]:
        excluded = set(influence.loc[mask, "SNP"].astype(str))
        subset = dat.loc[~dat["SNP"].astype(str).isin(excluded)]
        if len(subset) >= 5 and len(subset) < len(dat):
            result.extend(result_rows(outcome, method, fit_ivw(subset), len(subset)))
    for k in [1, 3]:
        excluded = set(influence.head(k)["SNP"].astype(str))
        subset = dat.loc[~dat["SNP"].astype(str).isin(excluded)]
        result.extend(result_rows(outcome, f"exclude_top_{k}_Cook", fit_ivw(subset), len(subset)))

    return pd.DataFrame(result), loo, influence, strat, groups, primary, robust


def main():
    all_results, all_loo, all_influence, all_strat, all_groups = [], [], [], [], []
    summary = {"outcomes": {}, "notes": []}
    for filename, outcome in [
        ("cpi_full_mvmr_harmonised_usable.csv", "CPI_full"),
        ("cpi_no_imaging_mvmr_harmonised_usable.csv", "CPI_no_imaging"),
    ]:
        results, loo, influence, strat, groups, primary, robust = run_outcome(filename, outcome)
        all_results.append(results)
        all_loo.append(loo)
        all_influence.append(influence)
        all_strat.append(strat)
        all_groups.append(groups.assign(outcome=outcome))
        summary["outcomes"][outcome] = {
            "nsnp": int(len(influence)),
            "cook_gt_4_over_n": int(influence["cook_gt_4_over_n"].sum()),
            "abs_stdresid_gt_3": int(influence["abs_stdresid_gt_3"].sum()),
            "top_cook_snps": influence.head(10)[["SNP", "cook_distance", "std_residual"]].to_dict(orient="records"),
            "tool_group_counts": groups["tool_group"].value_counts().to_dict(),
            "strength_group_counts": groups["strength_group"].value_counts().to_dict(),
            "huber_downweighted_n": int(np.sum(robust["weights"] < 0.999999)),
        }
    pd.concat(all_results, ignore_index=True).to_csv(OUT / "stage3b_robustness_results.csv", index=False)
    pd.concat(all_loo, ignore_index=True).to_csv(OUT / "stage3b_leave_one_out.csv", index=False)
    pd.concat(all_influence, ignore_index=True).to_csv(OUT / "stage3b_influence_diagnostics.csv", index=False)
    pd.concat(all_strat, ignore_index=True).to_csv(OUT / "stage3b_stratified_results.csv", index=False)
    pd.concat(all_groups, ignore_index=True).to_csv(OUT / "stage3b_variant_groups.csv", index=False)
    summary["notes"] = [
        "Primary estimates use the same multiplicative-random-effects convention as Stage 3.",
        "Huber-IRLS is exploratory and is not a substitute for a package-specific robust MVMR estimator.",
        "External phenotype databases were not queried automatically; top influence SNPs are listed for manual PhenoScanner/GWAS Catalog review.",
        "No R code was executed in Stage 3B.",
    ]
    (OUT / "stage3b_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
