#!/usr/bin/env python
"""
One-pass recomputation of every manuscript number affected by

  (1) SD_BMI: 1.23075 (data-implied, contaminated by the GIANT SE inflation)
              -> 1.0    (the BMI summary statistics are already per SD of BMI)
  (2) the inference rule: P values move from Student t (df 91 / 90) to the
      standard normal, so that they match the +/- 1.96 SE intervals that the
      manuscript already reports.

Pure standard library so it runs anywhere in the project sandbox.

Inputs (all post-P0-1, verified against _backup_before_P0fix)
-------------------------------------------------------------
  03_stage3_mvmr/stage3_mvmr_methods.csv          MVMR IVW (FE and MRE)
  03_stage3_mvmr/stage3_mvmr_diagnostics.csv      MVMR-Egger slopes/intercept
  03_stage3_mvmr/stage2b_power/...json            outcome SDs and MDE machinery
  04_stage4_sampling_cov/stage4_results/stage4_bias_correction.csv
  05_stage5_limbpain/results/stage5_mr_results_scalecorrected.csv
  05_stage5_limbpain/results/stage5_overlap_bias_scalecorrected.csv
"""
import csv
import json
import math
from pathlib import Path
from statistics import NormalDist

ROOT = Path(r"D:\SX\TMFI_CPI_MR")
S3 = ROOT / "03_stage3_mvmr"
S4 = ROOT / "04_stage4_sampling_cov" / "stage4_results"
S5 = ROOT / "05_stage5_limbpain" / "results"

N = NormalDist()
Z_CRIT = 1.959964  # the value the manuscript already uses for its intervals

SD_TMFI = 2.030884204160846
SD_BMI_OLD = 1.2307496342161395
SD_BMI_NEW = 1.0
K_SE = 1.2634  # SE inflation of the GIANT BMI release (1/sqrt(LDSC intercept 0.6265))


def pz(z):
    return 2 * N.cdf(-abs(z))


def two_sided_power(effect, se):
    mu = abs(effect) / se
    return N.cdf(-Z_CRIT - mu) + 1 - N.cdf(Z_CRIT - mu)


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
    # ------------------------------------------------------------------
    # outcome SDs (calibrated GWAS: LDSC intercept 1.0033 / 0.9956)
    # ------------------------------------------------------------------
    s2b = json.loads((S3 / "stage2b_power" / "stage2b_mvmr_power_summary.json")
                     .read_text(encoding="utf-8"))
    SD_Y = {oc: s2b["results"][f"{oc}|TMFI"]["implied_outcome_sd"]
            for oc in ("CPI_full", "CPI_no_imaging")}

    out = []
    w = out.append

    w("## A. Exposure SD")
    w(f"  SD_TMFI      {SD_TMFI:.6f}   (implied; LDSC intercept 0.9956 -> calibrated)")
    w(f"  SD_BMI old   {SD_BMI_OLD:.6f}   (implied; LDSC intercept 0.6265 -> contaminated)")
    w(f"  SD_BMI new   {SD_BMI_NEW:.6f}   (betas are per SD of BMI)")
    w(f"  rescale factor for every BMI quantity = {SD_BMI_NEW/SD_BMI_OLD:.6f}")
    w("")

    # ------------------------------------------------------------------
    # B. Table 3 - MVMR primary
    # ------------------------------------------------------------------
    meth = load(S3 / "stage3_mvmr_methods.csv")
    diag = {r["outcome"]: r for r in load(S3 / "stage3_mvmr_diagnostics.csv")}
    w("## B. Table 3 / Table 5 - MVMR estimates (P: t -> z)")
    w(f"  {'outcome':15s} {'method':6s} {'expo':5s} "
      f"{'beta':>9s} {'se':>8s} {'CI (1.96SE)':>22s} {'P_old(t)':>10s} {'P_new(z)':>11s} "
      f"{'std_old':>8s} {'std_new':>8s} {'stdCI_new':>20s}")
    for oc in ("CPI_full", "CPI_no_imaging"):
        sd_y = SD_Y[oc]
        rows = [r for r in meth if r["outcome"] == oc]
        for r in rows:
            b, se = float(r["beta"]), float(r["se"])
            expo = r["exposure"]
            sd_x_old = SD_TMFI if expo == "TMFI" else SD_BMI_OLD
            sd_x_new = SD_TMFI if expo == "TMFI" else SD_BMI_NEW
            z = b / se
            pt = float(r["p_t_df91"])
            pn = pz(z)
            std_old = b * sd_x_old / sd_y
            std_new = b * sd_x_new / sd_y
            lo, hi = b - Z_CRIT * se, b + Z_CRIT * se
            slo, shi = std_new - Z_CRIT * se * sd_x_new / sd_y, \
                       std_new + Z_CRIT * se * sd_x_new / sd_y
            w(f"  {oc:15s} {r['method'].replace('MVMR_IVW_',''):6s} {expo:5s} "
              f"{b:9.4f} {se:8.4f} {lo:10.4f} to {hi:9.4f} {pt:10.3g} {pn:11.3g} "
              f"{std_old:8.4f} {std_new:8.4f} {slo:9.4f} to {shi:8.4f}".replace(
                  f"{lo:10.4f} to {hi:9.4f}", f"[{lo:.4f},{hi:.4f}]"))
        # Egger
        d = diag[oc]
        for lab, b, se, dfold in (
                ("Egger", float(d["egger_tmfi_beta"]), float(d["egger_tmfi_se_mre"]), 90),
                ("Egger", float(d["egger_bmi_beta"]), float(d["egger_bmi_se_mre"]), 90)):
            expo = "TMFI" if "tmfi" in lab.lower() else "BMI"
            pass
        for expo, b, se in (("TMFI", float(d["egger_tmfi_beta"]), float(d["egger_tmfi_se_mre"])),
                            ("BMI", float(d["egger_bmi_beta"]), float(d["egger_bmi_se_mre"]))):
            sd_x_old = SD_TMFI if expo == "TMFI" else SD_BMI_OLD
            sd_x_new = SD_TMFI if expo == "TMFI" else SD_BMI_NEW
            z = b / se
            lo, hi = b - Z_CRIT * se, b + Z_CRIT * se
            w(f"  {oc:15s} {'Egger':6s} {expo:5s} {b:9.4f} {se:8.4f} "
              f"[{lo:.4f},{hi:.4f}] {'-':>10s} {pz(z):11.3g} "
              f"{b*sd_x_old/sd_y:8.4f} {b*sd_x_new/sd_y:8.4f}")
        b, se = float(d["egger_intercept"]), float(d["egger_intercept_se_mre"])
        w(f"  {oc:15s} {'Egger':6s} {'int':5s} {b:9.5f} {se:8.5f} "
          f"[{b-Z_CRIT*se:.5f},{b+Z_CRIT*se:.5f}] {'-':>10s} {pz(b/se):11.3g}")
        w("")

    # ------------------------------------------------------------------
    # C. Table 6 - power / MDE
    # ------------------------------------------------------------------
    w("## C. Table 6 - power and minimum detectable effect (standardized scale)")
    w(f"  {'outcome':15s} {'expo':5s} {'std_old':>8s} {'std_new':>8s} "
      f"{'CI_new':>20s} {'MDE80_old':>10s} {'MDE80_new':>10s} "
      f"{'MDE90_new':>10s} {'power':>7s}")
    for oc in ("CPI_full", "CPI_no_imaging"):
        sd_y = SD_Y[oc]
        for expo in ("TMFI", "BMI"):
            sub = [r for r in meth if r["outcome"] == oc and r["exposure"] == expo
                   and r["method"] == "MVMR_IVW_multiplicative_random_effects"][0]
            b, se = float(sub["beta"]), float(sub["se"])
            for tag, sdx in (("old", SD_TMFI if expo == "TMFI" else SD_BMI_OLD),
                             ("new", SD_TMFI if expo == "TMFI" else SD_BMI_NEW)):
                sc = sdx / sd_y
                if tag == "new":
                    std, sse = b * sc, se * sc
                    m80, m90 = exact_mde(sse, .80), exact_mde(sse, .90)
                    pw = two_sided_power(std, sse)
                    w(f"  {oc:15s} {expo:5s} {'':8s} {std:8.4f} "
                      f"[{std-Z_CRIT*sse:.4f},{std+Z_CRIT*sse:.4f}] {'':10s} "
                      f"{m80:10.4f} {m90:10.4f} {pw:7.3f}")
                else:
                    std_o = b * sdx / sd_y
                    w(f"  {oc:15s} {expo:5s} {std_o:8.4f} {'':8s} {'':20s} "
                      f"{exact_mde(se*sdx/sd_y,.80):10.4f} {'':10s} {'':10s} {'':7s}")
    w("")

    # ------------------------------------------------------------------
    # D. Table 7 - Stage 4 overlap correction
    # ------------------------------------------------------------------
    w("## D. Table 7 - Stage 4 overlap bias correction (BMI columns rescale)")
    f = SD_BMI_NEW / SD_BMI_OLD
    for r in load(S4 / "stage4_bias_correction.csv"):
        if r["stage"] != "Stage 2 multivariable IVW":
            continue
        oc, sc = r["outcome"], r["scenario"]
        for pre in ("bmi", "tmfi"):
            col = f"{pre}_corrected_std"
            if r[col] in ("NA", ""):
                continue
            old = float(r[col])
            new = old * (f if pre == "bmi" else 1.0)
            bias_old = float(r[f"{pre}_bias_std"]) if r[f"{pre}_bias_std"] != "NA" else None
            bias_new = bias_old * f if (bias_old is not None and pre == "bmi") else bias_old
            bs = "n/a" if bias_old is None else f"{bias_old:+.5f} -> {bias_new:+.5f}"
            w(f"  {oc:15s} {sc[:46]:46s} {pre.upper():4s} "
              f"corrected {old:+.5f} -> {new:+.5f}   bias {bs}   "
              f"P_after {r[f'{pre}_p_after']}")
    w("")

    # ------------------------------------------------------------------
    # E. Table 8 - Stage 5
    # ------------------------------------------------------------------
    w("## E. Table 8 - site-matched secondary analysis")
    w("  OR_new = OR_old ** (SD_BMI_new / SD_BMI_old)  for BMI rows only")
    s5 = load(S5 / "stage5_mr_results_scalecorrected.csv")
    mde_tmfi = {}
    or_bmi = {}
    for r in s5:
        tag = f"{r['stage']} {r['outcome']} {r['method']}"
        if ":BMI" not in r["method"] and r["method"] != "UVMR BMI":
            if r["mde_or_per_sd"]:
                mde_tmfi[(r["outcome"], r["stage"])] = float(r["mde_or_per_sd"])
        if ":BMI" in r["method"]:
            or_bmi[r["outcome"]] = float(r["or_per_sd"])
    for r in s5:
        if ":BMI" not in r["method"] and r["method"] != "UVMR BMI":
            continue
        o_old = float(r["or_per_sd"])
        lo_old, hi_old = float(r["ci_low"]), float(r["ci_high"])
        m_old = float(r["mde_or_per_sd"]) if r["mde_or_per_sd"] else None
        o_new = o_old ** f
        lo_new, hi_new = lo_old ** f, hi_old ** f
        m_new = m_old ** f if m_old else None
        w(f"  {r['stage']} {r['outcome']:20s} {r['method']:12s} "
          f"OR {o_old:.4f} -> {o_new:.4f}   CI [{lo_old:.4f},{hi_old:.4f}] -> "
          f"[{lo_new:.4f},{hi_new:.4f}]   MDE {('%.4f'%m_old) if m_old else '-':>7s} -> "
          f"{('%.4f'%m_new) if m_new else '-':>7s}   P {float(r['p']):.3g}")
    w("")
    w("  Panel B - MDE as a fraction of the BMI effect (log scale)")
    for oc in ("KNEE_PAIN_UKB", "BACK_PAIN_UKB", "KNEE_PAIN_FIGSHARE"):
        m_old = mde_tmfi.get((oc, "1 UVMR"))
        b_old = or_bmi.get(oc)
        if m_old is None:
            continue
        frac_old = math.log(m_old) / math.log(b_old) if b_old else None
        frac_new = math.log(m_old) / (math.log(b_old) * f) if b_old else None
        w(f"    {oc:20s} MDE {m_old:.4f}  BMI {b_old:.4f} -> {b_old**f:.4f}   "
          f"fraction {frac_old:.3f} -> {frac_new:.3f}")
    w("")

    w("## F. Instrument strength (BMI only; conservative because of the SE inflation)")
    w(f"  reported conditional F (Stage 5 closed form)  42.81  -> SE-rescaled "
      f"{42.812489288510605 * K_SE**2:.1f}")
    w(f"  reported conditional F (R MVMR, Table 5)      35.23  -> SE-rescaled "
      f"{35.2337665188432 * K_SE**2:.1f}")
    w(f"  TMFI conditional F is unaffected (LDSC intercept 0.9956)")

    print("\n".join(out))


if __name__ == "__main__":
    main()
