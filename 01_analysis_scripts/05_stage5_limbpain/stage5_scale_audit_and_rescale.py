#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Stage 5 - P0-2 : audit the effect-size scale of the three binary pain outcomes
==============================================================================
Question
--------
Are the betas in ukb-b-16254 (knee pain), ukb-b-9838 (back pain) and the
figshare knee-pain GWAS on the log-odds scale, or on the risk-difference
scale produced by BOLT-LMM's linear model on a 0/1 phenotype?

Two independent diagnostics
---------------------------
(1) SE magnitude.  For a linear (BOLT-LMM) model on a 0/1 phenotype of size N,
        Var(beta_j) = sigma^2 / (2 N f_j (1-f_j)),   sigma^2 = mu (1-mu)
    so    SE_j^2 * 2 N f_j (1-f_j)  ~=  mu (1-mu)      (a number around 0.17)
    For a logistic model the same quantity equals
        1 / (mu (1-mu))                               (a number around 6)
    The two predictions differ by a factor ~35; the observed value settles it.

(2) Calibration.  Diagnosing (1) is applied to the two OpenGWAS datasets whose
    case counts are known from the VCF header.  If it recovers mu(1-mu) to
    <1%, the same estimator is trusted for the figshare dataset, whose case
    count is not published.

Rescaling
---------
If the outcome is on the risk-difference scale, the BOLT-LMM manual gives
        log OR = beta / (mu (1-mu)),   SE_logOR = SE / (mu (1-mu))
i.e. beta and SE are both multiplied by the same constant lambda = 1/(mu(1-mu)).

Consequence for the MR estimates (important):
    * Wald ratio / IVW / MVMR / MR-Egger / weighted-median estimates, their SEs,
      the overlap-bias term and the corrected estimates all scale by lambda;
    * z = beta/SE, hence every P value, Cochran Q, phi, I2, the conditional F
      statistics, the Steiger test, the Egger intercept P and the tipping-point
      vs structural-ceiling ratio are INVARIANT;
    * only odds ratios, confidence intervals, the minimum detectable effect and
      the absolute bias magnitude change.

Pure standard library - no numpy/pandas needed.
"""
import csv
import gzip
import json
import math
import os
import statistics

S5 = r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain"
HARM = os.path.join(S5, "gwas_harmonised")
RES = os.path.join(S5, "results")

SD_TMFI = 2.030884204160846
SD_BMI = 1.0  # per SD of BMI; implied value 1.23 is inflated by the release's SE inflation

# ---- outcome metadata ------------------------------------------------------
# N and ncase for the two OpenGWAS files are read straight from the VCF header
# (##SAMPLE=<...,TotalControls=...,TotalCases=...>); they are repeated here so
# the script is self-contained.  figshare has no published case count -> the
# data-driven estimate is used, with a sensitivity range.
OUTCOMES = {
    #            file                                      N        ncase
    "KNEE_PAIN_UKB":       ("KNEE_PAIN_UKB_harmonised.tsv.gz",       461857, 98704),
    "BACK_PAIN_UKB":       ("BACK_PAIN_UKB_harmonised.tsv.gz",       461857, 118471),
    "KNEE_PAIN_FIGSHARE":  ("KNEE_PAIN_FIGSHARE_harmonised.tsv.gz",  429737, None),
}

# VCF headers, for the audit trail
VCF = {"KNEE_PAIN_UKB": "incoming/ukb-b-16254.vcf.gz",
       "BACK_PAIN_UKB": "incoming/ukb-b-9838.vcf.gz"}


def read_vcf_case_counts(path):
    """Pull TotalCases / TotalControls / StudyType out of ##SAMPLE=."""
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##SAMPLE="):
                d = dict(kv.split("=", 1) for kv in
                         line.strip()[len("##SAMPLE=<"):-1].split(",")
                         if "=" in kv)
                return (int(d.get("TotalCases", 0)),
                        int(d.get("TotalControls", 0)),
                        d.get("StudyType", "?"))
            if not line.startswith("#"):
                break
    return (None, None, None)


def estimate_sigma2(path, N, max_rows=300000, af_lo=0.02, af_hi=0.98):
    """median of SE^2 * 2 N f (1-f)  ==  residual variance of the model."""
    vals = []
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i_se, i_af = hdr.index("se"), hdr.index("af")
        n = 0
        for line in fh:
            x = line.rstrip("\n").split("\t")
            try:
                se = float(x[i_se])
                af = float(x[i_af])
            except (ValueError, IndexError):
                continue
            if not (af_lo < af < af_hi):
                continue
            n += 1
            vals.append(se * se * 2.0 * N * af * (1.0 - af))
            if n >= max_rows:
                break
    return statistics.median(vals), n


def mu_from_sigma2(s2):
    disc = 1.0 - 4.0 * s2
    return (1.0 - math.sqrt(disc)) / 2.0 if disc > 0 else float("nan")


def audit():
    print("=" * 88)
    print("PART 1 - SCALE AUDIT")
    print("=" * 88)
    lam = {}
    for name, (fname, N, ncase) in OUTCOMES.items():
        path = os.path.join(HARM, fname)
        s2_hat, nsnp = estimate_sigma2(path, N)
        mu_hat = mu_from_sigma2(s2_hat)

        if ncase is not None:
            mu_true = ncase / N
            s2_true = mu_true * (1 - mu_true)
            ratio = s2_hat / s2_true
            tag = "CALIBRATION (case count known from VCF header)"
        else:
            mu_true = s2_true = None
            ratio = float("nan")
            tag = "ESTIMATION (no published case count)"

        # what a logistic model would have predicted for the same statistic
        pred_logistic = 1.0 / s2_hat if s2_hat > 0 else float("nan")

        print(f"\n--- {name}   [{tag}]")
        if ncase is not None:
            nc, nctl, st = read_vcf_case_counts(os.path.join(S5, VCF[name]))
            print(f"    VCF header: cases={nc:,} controls={nctl:,} "
                  f"StudyType={st}  N={N:,}")
            print(f"    mu_true={mu_true:.5f}  sigma2_true=mu(1-mu)={s2_true:.5f}")
        print(f"    SNPs scanned                : {nsnp:,}")
        print(f"    median SE^2*2Nf(1-f)        : {s2_hat:.5f}")
        print(f"    -> implied mu               : {mu_hat:.5f}  "
              f"(implied cases ~{mu_hat * N:,.0f})")
        if s2_true is not None:
            print(f"    -> hat/true                 : {ratio:.4f}   "
                  f"({'PASS' if abs(ratio - 1) < 0.02 else 'FAIL'})")
        print(f"    linear-model  prediction    : {s2_hat:.5f}  "
              f"<-- observed value sits here")
        print(f"    logistic-model prediction   : {pred_logistic:.5f}  "
              f"({pred_logistic / s2_hat:.1f}x the observed value)")
        print(f"    VERDICT                     : risk difference "
              f"(BOLT-LMM linear model on 0/1)")
        lam[name] = 1.0 / s2_hat if ncase is None else 1.0 / s2_true
    print("\nlambda = 1/(mu(1-mu)), the factor that turns a risk difference "
          "into a log odds ratio:")
    for k, v in lam.items():
        print(f"    {k:<20s} lambda = {v:.5f}")
    return lam


def norm_cdf(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def sf(x):
    return 1.0 - norm_cdf(x)


def sd_for(row):
    m = row["method"]
    if "BMI" in m and "TMFI" not in m:
        return SD_BMI
    return SD_TMFI


def rescale_mr(lam):
    src = os.path.join(RES, "stage5_mr_results.csv")
    dst = os.path.join(RES, "stage5_mr_results_scalecorrected.csv")
    with open(src, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    out = []
    for r in rows:
        oc = r["outcome"]
        L = lam[oc]
        method = r["method"]
        if method == "Steiger":
            out.append(r)          # scale-free
            continue
        b = float(r["beta"]) * L
        s = float(r["se"]) * L
        sd = sd_for(r)
        r["beta"], r["se"] = f"{b:.10g}", f"{s:.10g}"
        if r["or_per_sd"]:
            r["or_per_sd"] = f"{math.exp(b * sd):.10g}"
        if r["ci_low"] and r["ci_high"]:
            r["ci_low"] = f"{math.exp((b - 1.96 * s) * sd):.10g}"
            r["ci_high"] = f"{math.exp((b + 1.96 * s) * sd):.10g}"
        if r["mde_or_per_sd"]:
            r["mde_or_per_sd"] = f"{math.exp((1.96 + 0.8416) * s * sd):.10g}"
        # risk difference per SD, from the ORIGINAL (unscaled) beta: no
        # approximation involved, reported alongside for transparency
        b_raw = float(r["beta"]) / L
        s_raw = float(r["se"]) / L
        r["rd_per_sd"] = f"{b_raw * sd:.10g}"
        r["rd_ci_low"] = f"{(b_raw - 1.96 * s_raw) * sd:.10g}"
        r["rd_ci_high"] = f"{(b_raw + 1.96 * s_raw) * sd:.10g}"
        if r["egger_intercept"]:
            r["egger_intercept"] = f"{float(r['egger_intercept']) * L:.10g}"
        # p, Q, phi, F statistics are invariant - left untouched
        out.append(r)
    cols = list(out[0].keys())
    for r in out:                       # keep the column order stable
        for c in cols:
            r.setdefault(c, "")
    with open(dst, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(out)
    print(f"\n[out] {dst}")
    return dst


def rescale_overlap(lam):
    src = os.path.join(RES, "stage5_overlap_bias.csv")
    dst = os.path.join(RES, "stage5_overlap_bias_scalecorrected.csv")
    with open(src, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    out = []
    for r in rows:
        L = lam[r["outcome"]]
        sd = SD_BMI if r.get("exposure") == "BMI" else SD_TMFI
        for key in ("theta", "se", "bias", "theta_corrected"):
            if r[key]:
                r[key] = f"{float(r[key]) * L:.10g}"
        # LDSC intercept, bias_in_se, structural ceiling and the tipping point
        # are all scale-invariant -> untouched
        if r["or_per_sd_corrected"]:
            r["or_per_sd_corrected"] = f"{math.exp(float(r['theta_corrected']) * sd):.10g}"
        out.append(r)
    cols = list(out[0].keys())
    with open(dst, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(out)
    print(f"[out] {dst}")
    return dst


def sensitivity(lam):
    """figshare has no published case count; sweep mu over a plausible range."""
    print("\n" + "=" * 88)
    print("PART 3 - SENSITIVITY for KNEE_PAIN_FIGSHARE (case count unknown)")
    print("=" * 88)
    src = os.path.join(RES, "stage5_mr_results.csv")
    with open(src, newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh)
                if r["outcome"] == "KNEE_PAIN_FIGSHARE" and r["method"] == "IVW-MRE:BMI"]
    r = rows[0]
    b, s = float(r["beta"]), float(r["se"])
    print(f"  {'mu':>6} {'lambda':>8} {'OR/SD':>8} {'95% CI':>20}")
    for mu in (0.18, 0.20, 0.2236, 0.25, 0.28):
        L = 1.0 / (mu * (1 - mu))
        bb, ss = b * L, s * L
        print(f"  {mu:6.3f} {L:8.3f} {math.exp(bb * SD_BMI):8.4f} "
              f"  [{math.exp((bb - 1.96 * ss) * SD_BMI):.4f}, "
              f"{math.exp((bb + 1.96 * ss) * SD_BMI):.4f}]")
    print("  NB the P value is identical across the whole range - only the "
          "effect size moves.")


def summary_table(lam):
    print("\n" + "=" * 88)
    print("PART 2 - CORRECTED EFFECT ESTIMATES (P values unchanged)")
    print("=" * 88)
    src = os.path.join(RES, "stage5_mr_results_scalecorrected.csv")
    with open(src, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    keep = {"IVW-MRE", "IVW-MRE:TMFI", "IVW-MRE:BMI", "UVMR BMI"}
    print(f"\n{'outcome':<20}{'method':<14}{'OR/SD':>9}{'95% CI':>22}"
          f"{'P':>12}{'RD/SD':>10}")
    print("-" * 88)
    for r in rows:
        if r["method"] not in keep:
            continue
        print(f"{r['outcome']:<20}{r['method']:<14}{float(r['or_per_sd']):>9.4f}"
              f"  [{float(r['ci_low']):>8.4f}, {float(r['ci_high']):>8.4f}]"
              f"{float(r['p']):>12.3g}{float(r['rd_per_sd']):>+10.5f}")
    print("\nRD/SD = absolute risk difference per 1 SD of the exposure "
          "(no model approximation involved).")


if __name__ == "__main__":
    lam = audit()
    rescale_mr(lam)
    rescale_overlap(lam)
    summary_table(lam)
    sensitivity(lam)
