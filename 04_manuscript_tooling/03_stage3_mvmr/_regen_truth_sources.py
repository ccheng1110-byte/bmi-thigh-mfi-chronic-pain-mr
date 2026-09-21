#!/usr/bin/env python
"""Regenerate the two `single source of truth` artefacts.

Both were left at their pre-P0-1 values while stage3_mvmr_methods.csv moved on,
so the manuscript could quote either. This rebuilds them from the current
(post-P0-1) exports with

  * sd_BMI = 1.00   (GIANT coefficients are per SD of BMI)
  * P values from the standard normal, matching the +/-1.96 SE intervals

and keeps the t-based values in a separate column for reference.

Pure standard library (no numpy/scipy) so it runs in the sandbox.
"""
import csv
import json
import math
from pathlib import Path
from statistics import NormalDist

S3 = Path(r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr")
N = NormalDist()
ZCRIT = 1.959964
SD_TMFI = 2.030884204160846
SD_BMI = 1.0


# ---------------------------------------------------------------- t distribution
def _betacf(a, b, x):
    tiny, eps = 1e-30, 3e-16
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < eps:
            break
    return h


def _betai(a, b, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
             + a * math.log(x) + b * math.log(1 - x))
    bt = math.exp(lbeta)
    if x < (a + 1) / (a + b + 2):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1 - x) / b


def p_t(t, df):
    return _betai(df / 2.0, 0.5, df / (df + t * t))


def two_sided_power(effect, se):
    mu = abs(effect) / se
    return N.cdf(-ZCRIT - mu) + 1 - N.cdf(ZCRIT - mu)


def exact_mde(se, target):
    lo, hi = 0.0, 10 * se
    for _ in range(200):
        mid = (lo + hi) / 2
        if two_sided_power(mid, se) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def load(p):
    return list(csv.DictReader(open(p, newline="", encoding="utf-8")))


def main():
    # outcome SDs: calibrated GWAS, unchanged by P0-1 (they depend only on SE)
    s2b = json.loads((S3 / "stage2b_power" / "stage2b_mvmr_power_summary.json")
                     .read_text(encoding="utf-8"))
    SD_Y = {oc: s2b["results"][f"{oc}|TMFI"]["implied_outcome_sd"]
            for oc in ("CPI_full", "CPI_no_imaging")}

    meth = load(S3 / "stage3_mvmr_methods.csv")
    diag = {r["outcome"]: r for r in load(S3 / "stage3_mvmr_diagnostics.csv")}

    # ------------------------- stage3_mvmr_methods_unified.csv -----------------
    rows = []
    for oc in ("CPI_full", "CPI_no_imaging"):
        sd_y = SD_Y[oc]
        d = diag[oc]
        specs = []
        for meth_name, short in (("MVMR_IVW_multiplicative_random_effects",
                                  "IVW-multiplicative_random_effects"),
                                 ("MVMR_IVW_fixed", "IVW-fixed")):
            for expo in ("TMFI", "BMI"):
                sub = [r for r in meth if r["outcome"] == oc
                       and r["exposure"] == expo and r["method"] == meth_name]
                if sub:
                    specs.append((short, expo, float(sub[0]["beta"]),
                                  float(sub[0]["se"]), 91))
        specs.append(("Egger", "TMFI", float(d["egger_tmfi_beta"]),
                      float(d["egger_tmfi_se_mre"]), 90))
        specs.append(("Egger", "BMI", float(d["egger_bmi_beta"]),
                      float(d["egger_bmi_se_mre"]), 90))
        specs.append(("Egger-intercept", "NA", float(d["egger_intercept"]),
                      float(d["egger_intercept_se_mre"]), 90))
        for name, expo, est, se, df in specs:
            if expo == "TMFI":
                sd_x = SD_TMFI
            elif expo == "BMI":
                sd_x = SD_BMI
            else:
                sd_x = float("nan")
            scale = sd_x / sd_y if sd_x == sd_x else float("nan")
            rows.append({
                "outcome": oc, "method": name, "exposure": expo,
                "nsnp": 93, "df_p": df,
                "beta": repr(est), "se": repr(se),
                "ci_low": repr(est - ZCRIT * se), "ci_high": repr(est + ZCRIT * se),
                "p_t": repr(p_t(abs(est / se), df)),
                "p_z": repr(2 * N.cdf(-abs(est / se))),
                "sd_exposure": "" if sd_x != sd_x else repr(sd_x),
                "sd_outcome": repr(sd_y),
                "beta_standardized": "" if scale != scale else repr(est * scale),
                "se_standardized": "" if scale != scale else repr(se * scale),
            })
    with open(S3 / "stage3_mvmr_methods_unified.csv", "w", newline="",
              encoding="utf-8") as fh:
        wtr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        wtr.writeheader()
        wtr.writerows(rows)
    print(f"wrote stage3_mvmr_methods_unified.csv  ({len(rows)} rows)")

    # ------------------------- stage2b_mvmr_power_summary.json -----------------
    results = {}
    for oc in ("CPI_full", "CPI_no_imaging"):
        sd_y = SD_Y[oc]
        for expo, sd_x in (("TMFI", SD_TMFI), ("BMI", SD_BMI)):
            sub = [r for r in meth if r["outcome"] == oc and r["exposure"] == expo
                   and r["method"] == "MVMR_IVW_multiplicative_random_effects"][0]
            b, se = float(sub["beta"]), float(sub["se"])
            scale = sd_x / sd_y
            b_std, se_std = b * scale, se * scale
            results[f"{oc}|{expo}"] = {
                "implied_outcome_sd": sd_y,
                "implied_exposure_sd": sd_x,
                "mvmr_beta_raw": b,
                "mvmr_se_raw": se,
                "mvmr_beta_standardized_sd_y_per_sd_exposure": b_std,
                "mvmr_se_standardized": se_std,
                "ci95_standardized": [b_std - ZCRIT * se_std,
                                      b_std + ZCRIT * se_std],
                "mde_80_raw": exact_mde(se, 0.80),
                "mde_90_raw": exact_mde(se, 0.90),
                "mde_80_standardized": exact_mde(se_std, 0.80),
                "mde_90_standardized": exact_mde(se_std, 0.90),
                "power_at_observed_effect": two_sided_power(b_std, se_std),
            }
    summary = {
        "parameters": {
            "alpha_two_sided": 0.05,
            "estimator": "Stage 3 MVMR IVW multiplicative random effects",
            "estimand": "direct effect conditional on the other exposure",
            "power_method": "exact two-sided normal approximation using "
                            "observed MVMR-MRE SE",
            "scaling_method": "sd_exposure / sd_outcome; sd from SE, N and EAF "
                              "for TMFI and CPI; sd_BMI = 1.00 because the "
                              "GIANT coefficients are already per SD of BMI",
        },
        "exposure_sd": {
            "TMFI_stage2_value": SD_TMFI,
            "BMI_joint_set": SD_BMI,
            "BMI_joint_set_implied_not_used": 1.2307496342161395,
        },
        "results": results,
    }
    (S3 / "stage2b_power" / "stage2b_mvmr_power_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print("wrote stage2b_mvmr_power_summary.json")

    for k, v in results.items():
        bstd = v["mvmr_beta_standardized_sd_y_per_sd_exposure"]
        m80 = v["mde_80_standardized"]
        pw = v["power_at_observed_effect"]
        print(f"  {k:24s} std={bstd:+.4f}  mde80={m80:.4f}  power={pw:.3f}")


if __name__ == "__main__":
    main()
