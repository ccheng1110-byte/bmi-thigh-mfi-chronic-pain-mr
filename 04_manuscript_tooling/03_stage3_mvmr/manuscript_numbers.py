# -*- coding: utf-8 -*-
"""Single source of truth for every number quoted in the manuscript.

Purpose
-------
The draft mixed two conventions for MVMR P values (z-based Stage 3B exports and
t-based R MVMR output). This script recomputes everything under one declared rule
and prints a formatted table so the manuscript and the CSV outputs cannot drift.

Declared conventions (must be reproduced in Statistical reporting)
------------------------------------------------------------------
* Every estimator: two-sided P from the standard normal distribution applied
  to z = beta / SE. One rule, matching the +/- 1.96 * SE intervals, so that a
  reported interval excludes zero exactly when the reported P is below 0.05.
  (The earlier draft mixed conventions: t with df 91 for MVMR IVW and df 90
  for MVMR-Egger, while the intervals already used 1.96. The t rule is still
  exported as p_t for reference and is the more conservative of the two.)
* Univariable MR: P values are reported exactly as exported by TwoSampleMR
  and are not recomputed; its IVW and weighted-median tests are already
  normal, and only the two MR-Egger tests use t with 40 df (a difference
  below 0.003).
* Confidence intervals: estimate +/- 1.96 * SE throughout.
* Standardised scale:
      implied_sd(se, n, eaf) = se * sqrt(2 * n * eaf * (1 - eaf))
      scale                  = sd_exposure / sd_outcome

Coverage
--------
Stage 1 (univariable), Stage 2B (power / MDE), Stage 3 (MVMR primary, fixed
effects, Egger, diagnostics), Stage 3B (robustness), Stage 4 (sample-overlap
bias correction) and Stage 5 (site-matched secondary analysis, including the
independent base-R reproduction used as the agreement check).

Stage 3B and Stage 5 numbers are read from their frozen CSV exports rather than
recomputed here, because they are produced by separate pipelines
(_recompute_3b.R and stage5_mr.py respectively) that must be re-run after any
change to the harmonised inputs. After the P0-1 palindromic-orientation fix,
Stage 3B was recomputed and rewritten; the values printed here are therefore the
corrected ones. If the orientation or instrument set changes again, re-run
_recompute_3b.R and stage5_mr.py before trusting this output.
"""
import json
from math import sqrt
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(r"D:\SX\TMFI_CPI_MR")
STAGE1 = ROOT / "01_stage1_instruments"
STAGE4 = ROOT / "04_stage4_sampling_cov"
STAGE5 = ROOT / "05_stage5_limbpain"
OUT = ROOT / "03_stage3_mvmr"
STAGE2B = OUT / "stage2b_power"

SD_TMFI = 2.030884204160846  # from Stage 2 exposure SD, reused verbatim
ZCRIT   = 1.959964           # normal quantile; matches the reported intervals


def implied_sd(se, n, eaf):
    return se * sqrt(2 * n * eaf * (1 - eaf))


def main():
    diag = pd.read_csv(OUT / "stage3_mvmr_diagnostics.csv").set_index("outcome")
    joint = pd.read_csv(OUT / "cpi_full_mvmr_harmonised_usable.csv")
    methods = pd.read_csv(OUT / "stage3_mvmr_methods.csv")
    uv_methods = pd.read_csv(STAGE1 / "formal_mr_methods.csv")
    uv_diag = pd.read_csv(STAGE1 / "formal_mr_diagnostics.csv").set_index("outcome")
    s2b = json.loads((STAGE2B / "stage2b_mvmr_power_summary.json").read_text(encoding="utf-8"))

    def booly(v):
        return str(v).strip().upper() in ("TRUE", "1", "YES")

    bmi_joint = joint.loc[joint["BMI_GWS"].map(booly)]
    # GIANT coefficients are already per SD of BMI; the implied-SD
    # estimator returns 1.23 for them only because that release's
    # standard errors are ~26% inflated (single-trait LD score intercept
    # 0.62). Dividing it out gives 0.97-1.02, so sd_BMI = 1.00.
    sd_bmi = 1.0
    sd_bmi_implied = float(np.median(   # audit trail only
        bmi_joint.apply(lambda r: implied_sd(r["BMI_SE"], r["BMI_N"], r["BMI_EAF_HARMONISED"]), axis=1)
    ))

    print("### Data provenance")
    for nm, col_n, col_se, col_eaf, src_note in [
        ("BMI", "BMI_N", "BMI_SE", "BMI_EAF_HARMONISED", "joint set"),
        ("TMFI", "TMFI_N", "TMFI_SE", "TMFI_EAF", "joint set"),
        ("CPI_full", "OUTCOME_N", "OUTCOME_SE", "OUTCOME_EAF_H", "cpi_full set"),
    ]:
        sub = joint.dropna(subset=[col_n, col_se])
        print(f"  {nm:9s} ({src_note}) per-SNP N: {int(sub[col_n].min())}-{int(sub[col_n].max())}"
              f"  implied SD={np.median(sub.apply(lambda r: implied_sd(r[col_se], r[col_n], r[col_eaf]), axis=1)):.6f}")
    joint_ni = pd.read_csv(OUT / "cpi_no_imaging_mvmr_harmonised_usable.csv")
    print(f"  CPI_no_imaging per-SNP N: {int(joint_ni['OUTCOME_N'].min())}-{int(joint_ni['OUTCOME_N'].max())}")
    print(f"  sd_BMI_joint={sd_bmi:.6f}  (implied, not used: {sd_bmi_implied:.6f})  sd_TMFI={SD_TMFI:.6f}")

    print("\n### Table - univariable MR (as exported by TwoSampleMR, nSNP=42)")
    for _, r in uv_methods.iterrows():
        print(f"  {r['outcome']:16s} {r['method']:36s} beta={r['beta']:+.4f} se={r['se']:.4f} "
              f"ci=({r['ci_low']:+.4f},{r['ci_high']:+.4f}) p={r['p']:.3g}")
    for oc, r in uv_diag.iterrows():
        print(f"  {oc:16s} Q={r['cochran_q']:.2f} df={int(r['cochran_q_df'])} p={r['cochran_q_p']:.3g} "
              f"phi={r['phi']:.3f} int={r['egger_intercept']:+.5f}+-{r['egger_intercept_se']:.5f} "
              f"int_p={r['egger_intercept_p']:.3g} steiger={int(r['steiger_correct_count'])}/{int(r['steiger_total'])} "
              f"loo=({r['loo_min']:+.4f},{r['loo_max']:+.4f})")

    print("\n### Standardised univariable ivw-mre effects")
    for oc, se_key in [("CPI_full", "implied_outcome_sd"), ("CPI_no_imaging", None)]:
        row = uv_methods[(uv_methods["outcome"] == oc) &
                         (uv_methods["method"] == "IVW_multiplicative_random_effects")].iloc[0]
        sd_y = s2b["results"][f"{oc}|TMFI"]["implied_outcome_sd"]
        scale = SD_TMFI / sd_y
        est, se = row["beta"] * scale, row["se"] * scale
        print(f"  {oc:16s} beta_std={est:+.4f} se_std={se:.4f} "
              f"ci=({est-1.96*se:+.4f},{est+1.96*se:+.4f})")

    print("\n### Table - MVMR IVW and Egger (two-sided normal P)")
    for oc in ["CPI_full", "CPI_no_imaging"]:
        sd_y = s2b["results"][f"{oc}|TMFI"]["implied_outcome_sd"]
        for expo, sd_x in [("TMFI", SD_TMFI), ("BMI", sd_bmi)]:
            sub = methods[(methods["outcome"] == oc) & (methods["exposure"] == expo) &
                          (methods["method"] == "MVMR_IVW_multiplicative_random_effects")].iloc[0]
            p = 2 * stats.norm.sf(abs(sub["beta"] / sub["se"]))
            est, se = sub["beta"], sub["se"]
            scale = sd_x / sd_y
            print(f"  {oc:16s} MRE  {expo:5s} beta={est:+.4f} se={se:.4f} "
                  f"ci=({est-ZCRIT*se:+.4f},{est+ZCRIT*se:+.4f}) p={p:.3g} | "
                  f"std beta={est*scale:+.4f} ci=({(est-ZCRIT*se)*scale:+.4f},{(est+ZCRIT*se)*scale:+.4f})")
        d = diag.loc[oc]
        for kind, b, se in [("egger_TMFI", d["egger_tmfi_beta"], d["egger_tmfi_se_mre"]),
                            ("egger_BMI", d["egger_bmi_beta"], d["egger_bmi_se_mre"]),
                            ("egger_int", d["egger_intercept"], d["egger_intercept_se_mre"])]:
            p = 2 * stats.norm.sf(abs(b / se))
            print(f"  {oc:16s} {kind:10s} beta={b:+.5f} se={se:.5f} "
                  f"ci=({b-ZCRIT*se:+.5f},{b+ZCRIT*se:+.5f}) p={p:.3g}")

    print("\n### Fixed-effect MVMR (two-sided normal P)")
    for oc in ["CPI_full", "CPI_no_imaging"]:
        for expo in ["TMFI", "BMI"]:
            sub = methods[(methods["outcome"] == oc) & (methods["exposure"] == expo) &
                          (methods["method"] == "MVMR_IVW_fixed")].iloc[0]
            p = 2 * stats.norm.sf(abs(sub["beta"] / sub["se"]))
            print(f"  {oc:16s} FE   {expo:5s} beta={sub['beta']:+.4f} se={sub['se']:.4f} "
                  f"ci=({sub['beta']-ZCRIT*sub['se']:+.4f},{sub['beta']+ZCRIT*sub['se']:+.4f}) p={p:.3g}")

    print("\n### Diagnostics")
    for oc in ["CPI_full", "CPI_no_imaging"]:
        d = diag.loc[oc]
        print(f"  {oc:16s} Q_ivw={d['q_ivw_outcome_only']:.2f}(df={int(d['q_ivw_df'])}) "
              f"phi={d['phi']:.4f} Qval={d['q_validity_modified']:.2f}(df={int(d['q_validity_df'])},"
              f"p={d['q_validity_p']:.3g}) F_TMFI={d['tmfi_conditional_f']:.2f} "
              f"F_BMI={d['bmi_conditional_f']:.2f} egger_Q={d['egger_q']:.2f}(df={int(d['egger_df'])})")

    print("\n### Stage 2B power / MDE (standardised)")
    for k, v in s2b["results"].items():
        print(f"  {k:26s} beta_std={v['mvmr_beta_standardized_sd_y_per_sd_exposure']:+.4f} "
              f"se_std={v['mvmr_se_standardized']:.4f} "
              f"ci=({v['mvmr_beta_standardized_sd_y_per_sd_exposure']-1.959964*v['mvmr_se_standardized']:+.4f},"
              f"{v['mvmr_beta_standardized_sd_y_per_sd_exposure']+1.959964*v['mvmr_se_standardized']:+.4f}) "
              f"mde80={v['mde_80_standardized']:.4f} mde90={v['mde_90_standardized']:.4f}")

    # ------------------------------------------------------------------
    # Persist the single-convention outputs so the manuscript tables and the
    # frozen CSVs cannot drift apart.
    # ------------------------------------------------------------------
    rows = []
    for oc in ["CPI_full", "CPI_no_imaging"]:
        sd_y = s2b["results"][f"{oc}|TMFI"]["implied_outcome_sd"]
        d = diag.loc[oc]
        specs = []
        for meth in ["MVMR_IVW_multiplicative_random_effects", "MVMR_IVW_fixed"]:
            for expo in ["TMFI", "BMI"]:
                sub = methods[(methods["outcome"] == oc) & (methods["exposure"] == expo) &
                              (methods["method"] == meth)]
                if not sub.empty:
                    specs.append((meth.replace("MVMR_IVW_", "IVW-"), expo,
                                  sub.iloc[0]["beta"], sub.iloc[0]["se"], 91))
        specs.append(("Egger", "TMFI", d["egger_tmfi_beta"], d["egger_tmfi_se_mre"], 90))
        specs.append(("Egger", "BMI", d["egger_bmi_beta"], d["egger_bmi_se_mre"], 90))
        specs.append(("Egger-intercept", "NA", d["egger_intercept"], d["egger_intercept_se_mre"], 90))
        for meth, expo, est, se, df in specs:
            sd_x = SD_TMFI if expo == "TMFI" else (sd_bmi if expo == "BMI" else float("nan"))
            scale = sd_x / sd_y if sd_x == sd_x else float("nan")
            rows.append({
                "outcome": oc,
                "method": meth,
                "exposure": expo,
                "nsnp": 93,
                "df_p": df,
                "beta": est,
                "se": se,
                "ci_low": est - 1.959964 * se,
                "ci_high": est + 1.959964 * se,
                "p_t": 2 * stats.t.sf(abs(est / se), df),
                "p_z": 2 * stats.norm.sf(abs(est / se)),
                "sd_exposure": sd_x,
                "sd_outcome": sd_y,
                "beta_standardized": est * scale,
                "se_standardized": se * scale,
            })
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "stage3_mvmr_methods_unified.csv", index=False)
    print(f"\nWrote {OUT / 'stage3_mvmr_methods_unified.csv'}  rows={len(tab)}")

    # Append the unified columns to the original export so it stays in sync.
    orig = methods.copy()
    key = orig.set_index(["outcome", "exposure", "method"])
    orig["p_t_df91"] = [tab.loc[(tab["outcome"] == r["outcome"]) & (tab["exposure"] == r["exposure"]) &
                                (tab["method"] == ("IVW-multiplicative_random_effects"
                                                  if r["method"] == "MVMR_IVW_multiplicative_random_effects"
                                                  else "IVW-fixed")), "p_t"].iloc[0]
                        for _, r in orig.iterrows()]
    orig.to_csv(OUT / "stage3_mvmr_methods.csv", index=False)
    print(f"Updated {OUT / 'stage3_mvmr_methods.csv'} with column p_t_df91")

    # ------------------------------------------------------------------
    # Stage 3B robustness: influence-point and reweighting analyses.
    # The table is recomputed from the corrected (P0-1) harmonised inputs by
    # _recompute_3b.R whenever the orientation fix changes; it is read here so
    # the manuscript cannot quote a stale number.
    # ------------------------------------------------------------------
    print("\n### Table - Stage 3B robustness (influence and reweighting)")
    s3b = pd.read_csv(OUT / "stage3b_robustness_results.csv")
    for _, r in s3b.iterrows():
        print(f"  {r['outcome']:16s} {r['method']:26s} {r['exposure']:5s} n={int(r['nsnp'])} "
              f"beta={r['beta']:+.4f} se={r['se']:.4f} "
              f"ci=({r['ci_low']:+.4f},{r['ci_high']:+.4f}) p={r['p_normal']:.3g}")

    # ------------------------------------------------------------------
    # Stage 4: sample-overlap audit (cross-trait LDSC intercepts propagated
    # into MR bias, with the structural ceiling and the tipping point).
    # ------------------------------------------------------------------
    print("\n### Table 7 - sample-overlap bias correction (Stage 4)")
    s4 = pd.read_csv(STAGE4 / "stage4_results" / "stage4_bias_correction.csv")
    for _, r in s4.iterrows():
        print(f"  {r['stage']:28s} {r['outcome']:16s} {r['scenario']:16s} "
              f"int_tmfi={r['intercept_tmfi']:+.5f} n_shared_rho={r['n_shared_times_rho']:.1f} "
              f"bias_in_se={r['bias_in_se']:.4f} corr_std={r['beta_corr_std']:+.5f} "
              f"| tmfi_corr_std={r['tmfi_corrected_std']:+.5f} bmi_corr_std={r['bmi_corrected_std']:+.5f} "
              f"| tmfi_bias_in_se={r['tmfi_bias_in_se']:.4f} bmi_bias_in_se={r['bmi_bias_in_se']:.4f}")

    # ------------------------------------------------------------------
    # Stage 5: pre-specified site-matched secondary analysis.
    # Panel A = MR estimates, Panel B = power / overlap audit.
    # ------------------------------------------------------------------
    print("\n### Table 8 Panel A - site-matched MR (Stage 5)")
    s5 = pd.read_csv(STAGE5 / "results" / "stage5_mr_results.csv")
    for _, r in s5.iterrows():
        print(f"  {str(r['stage']):12s} {str(r['outcome']):18s} {str(r['method']):8s} n={int(r['n_snp'])} "
              f"OR={r['or_per_sd']:.4f} ci=({r['ci_low']:.4f},{r['ci_high']:.4f}) p={r['p']:.3g} "
              f"mde={r['mde_or_per_sd']:.4f} F_condA={r['F_conditional_A']:.2f} F_condB={r['F_conditional_B']:.2f}")

    print("\n### Table 8 Panel B - overlap audit for the secondary outcomes (Stage 5)")
    s5o = pd.read_csv(STAGE5 / "results" / "stage5_overlap_bias.csv")
    for _, r in s5o.iterrows():
        print(f"  {str(r['outcome']):18s} {str(r['stage']):12s} intercept={r['intercept']:+.5f} "
              f"ceiling={r['ceiling']:.4f} bias_in_se={r['bias_in_se']:.4f} "
              f"intercept_for_p05={r['intercept_for_p05']:.3f}")

    print("\n### Stage 5 independent base-R reproduction (agreement check)")
    s5b = pd.read_csv(STAGE5 / "results" / "stage5_baser_all.csv")
    for _, r in s5b.iterrows():
        print(f"  {str(r['outcome']):18s} {str(r['method']):8s} OR={r['OR']:.7f} "
              f"ci=({r['lo']:.7f},{r['hi']:.7f}) p={r['p']:.6g}")


if __name__ == "__main__":
    main()
