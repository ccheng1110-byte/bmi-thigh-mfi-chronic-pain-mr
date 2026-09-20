#!/usr/bin/env python
"""Propagate SD_BMI = 1.00 (and the normal-P rule) into the pipeline scripts.

Without this, any re-run would regenerate the numbers the manuscript has just
been corrected away from.
"""
from pathlib import Path

ROOT = Path(r"D:\SX\TMFI_CPI_MR")
S5 = ROOT / "05_stage5_limbpain"
S3 = ROOT / "03_stage3_mvmr"
S4 = ROOT / "04_stage4_sampling_cov"

NOTE = ("# The GIANT BMI coefficients are already per SD of BMI. The implied-SD\n"
        "# estimator returns 1.23 for them only because that release's standard\n"
        "# errors are ~26% inflated (single-trait LD score intercept 0.62);\n"
        "# dividing it out returns 0.97-1.02, so the exposure SD is 1.00.")

EDITS = [
    # ---------------- Stage 5 python ----------------
    (S5 / "stage5_mr.py",
     "SD_BMI = 1.2307496342161395",
     NOTE.replace("#", "#", 1) + "\nSD_BMI = 1.0"),
    (S5 / "stage5_overlap_bias.py",
     "SD_BMI = 1.2307496342161395",
     "SD_BMI = 1.0  # GIANT coefficients are per SD of BMI (intercept 0.62 "
     "implies ~26% SE inflation, so the implied 1.23 is not the exposure SD)"),
    (S5 / "stage5_scale_audit_and_rescale.py",
     "SD_BMI = 1.2307496342161395",
     "SD_BMI = 1.0  # per SD of BMI; implied value 1.23 is inflated by the "
     "release's SE inflation"),
    # ---------------- Stage 5 R ----------------
    (S5 / "stage5_validation.R",
     "SD_BMI   <- 1.2307496342161395",
     "SD_BMI   <- 1.0   # per SD of BMI (implied 1.23 reflects the release's "
     "~26% SE inflation)"),
    (S5 / "stage5_validation_baser.R",
     "SD_BMI   <- 1.2307496342161395",
     "SD_BMI   <- 1.0   # per SD of BMI (implied 1.23 reflects the release's "
     "~26% SE inflation)"),
    # ---------------- Stage 4 ----------------
    (S4 / "mrlap_bias_correction.py",
     "SD_BMI = 1.230750",
     "SD_BMI = 1.0  # per SD of BMI"),
    # ---------------- Stage 2B ----------------
    (S3 / "stage2b_mvmr_power.py",
     "    sd_bmi_joint = float(np.median(\n"
     "        bmi_gws.apply(lambda r: implied_sd(r[\"BMI_SE\"], r[\"BMI_N\"], "
     "r[\"BMI_EAF_HARMONISED\"]), axis=1)\n"
     "    ))",
     "    # GIANT coefficients are already per SD of BMI; the implied-SD\n"
     "    # estimator returns 1.23 only because that release's standard errors\n"
     "    # are ~26% inflated (single-trait LD score intercept 0.62).\n"
     "    sd_bmi_joint = 1.0\n"
     "    sd_bmi_joint_implied = float(np.median(   # audit trail only\n"
     "        bmi_gws.apply(lambda r: implied_sd(r[\"BMI_SE\"], r[\"BMI_N\"], "
     "r[\"BMI_EAF_HARMONISED\"]), axis=1)\n"
     "    ))"),
    (S3 / "stage2b_mvmr_power.py",
     '            "BMI_joint_set": sd_bmi_joint,',
     '            "BMI_joint_set": sd_bmi_joint,\n'
     '            "BMI_joint_set_implied_not_used": sd_bmi_joint_implied,'),
    # ---------------- manuscript_numbers.py ----------------
    (S3 / "manuscript_numbers.py",
     "    bmi_joint = joint.loc[joint[\"BMI_GWS\"].map(booly)]\n"
     "    sd_bmi = float(np.median(\n"
     "        bmi_joint.apply(lambda r: implied_sd(r[\"BMI_SE\"], r[\"BMI_N\"], "
     "r[\"BMI_EAF_HARMONISED\"]), axis=1)\n"
     "    ))",
     "    bmi_joint = joint.loc[joint[\"BMI_GWS\"].map(booly)]\n"
     "    # GIANT coefficients are already per SD of BMI; the implied-SD\n"
     "    # estimator returns 1.23 for them only because that release's\n"
     "    # standard errors are ~26% inflated (single-trait LD score intercept\n"
     "    # 0.62). Dividing it out gives 0.97-1.02, so sd_BMI = 1.00.\n"
     "    sd_bmi = 1.0\n"
     "    sd_bmi_implied = float(np.median(   # audit trail only\n"
     "        bmi_joint.apply(lambda r: implied_sd(r[\"BMI_SE\"], r[\"BMI_N\"], "
     "r[\"BMI_EAF_HARMONISED\"]), axis=1)\n"
     "    ))"),
    (S3 / "manuscript_numbers.py",
     '    print(f"  sd_BMI_joint={sd_bmi:.6f}  sd_TMFI={SD_TMFI:.6f}")',
     '    print(f"  sd_BMI_joint={sd_bmi:.6f}  (implied, not used: '
     '{sd_bmi_implied:.6f})  sd_TMFI={SD_TMFI:.6f}")'),
    # P convention -> normal, matching the +/-1.96 SE intervals
    (S3 / "manuscript_numbers.py",
     "* MVMR IVW (fixed and multiplicative random effects):\n"
     "      two-sided Student t with df = nSNP - number of exposures = 93 - 2 = 91\n"
     "* MVMR-Egger slopes and intercept:\n"
     "      two-sided Student t with df = nSNP - number of parameters = 93 - 3 = 90\n"
     "* Univariable MR: P values are reported exactly as exported by TwoSampleMR and\n"
     "  are not recomputed.\n"
     "* Confidence intervals: estimate +/- t_crit(df) * SE for MVMR regression-based\n"
     "  estimators; estimate +/- 1.96 * SE where the source workflow used the normal\n"
     "  approximation (Harmonised only what already exists in the frozen CSVs).",
     "* Every estimator: two-sided P from the standard normal distribution applied\n"
     "  to z = beta / SE. One rule, matching the +/- 1.96 * SE intervals, so that a\n"
     "  reported interval excludes zero exactly when the reported P is below 0.05.\n"
     "  (The earlier draft mixed conventions: t with df 91 for MVMR IVW and df 90\n"
     "  for MVMR-Egger, while the intervals already used 1.96. The t rule is still\n"
     "  exported as p_t for reference and is the more conservative of the two.)\n"
     "* Univariable MR: P values are reported exactly as exported by TwoSampleMR\n"
     "  and are not recomputed; its IVW and weighted-median tests are already\n"
     "  normal, and only the two MR-Egger tests use t with 40 df (a difference\n"
     "  below 0.003).\n"
     "* Confidence intervals: estimate +/- 1.96 * SE throughout."),
    (S3 / "manuscript_numbers.py",
     '    print("\\n### Table - MVMR IVW (t with df=91) and Egger (t with df=90)")',
     '    print("\\n### Table - MVMR IVW and Egger (two-sided normal P)")'),
    (S3 / "manuscript_numbers.py",
     '            tcrit = stats.t.ppf(0.975, 91)\n'
     '            p = 2 * stats.t.sf(abs(sub["beta"] / sub["se"]), 91)\n'
     '            est, se = sub["beta"], sub["se"]\n'
     '            scale = sd_x / sd_y\n'
     '            print(f"  {oc:16s} MRE  {expo:5s} beta={est:+.4f} se={se:.4f} "\n'
     '                  f"ci=({est-tcrit*se:+.4f},{est+tcrit*se:+.4f}) p={p:.3g} | "',
     '            p = 2 * stats.norm.sf(abs(sub["beta"] / sub["se"]))\n'
     '            est, se = sub["beta"], sub["se"]\n'
     '            scale = sd_x / sd_y\n'
     '            print(f"  {oc:16s} MRE  {expo:5s} beta={est:+.4f} se={se:.4f} "\n'
     '                  f"ci=({est-ZCRIT*se:+.4f},{est+ZCRIT*se:+.4f}) p={p:.3g} | "'),
    (S3 / "manuscript_numbers.py",
     '                  f"std beta={est*scale:+.4f} ci=({(est-tcrit*se)*scale:+.4f},{(est+tcrit*se)*scale:+.4f})")',
     '                  f"std beta={est*scale:+.4f} ci=({(est-ZCRIT*se)*scale:+.4f},{(est+ZCRIT*se)*scale:+.4f})")'),
    (S3 / "manuscript_numbers.py",
     '            p = 2 * stats.t.sf(abs(b / se), 90)\n'
     '            print(f"  {oc:16s} {kind:10s} beta={b:+.5f} se={se:.5f} "\n'
     '                  f"ci=({b-1.959964*se:+.5f},{b+1.959964*se:+.5f}) p_z={p:.3g}")',
     '            p = 2 * stats.norm.sf(abs(b / se))\n'
     '            print(f"  {oc:16s} {kind:10s} beta={b:+.5f} se={se:.5f} "\n'
     '                  f"ci=({b-ZCRIT*se:+.5f},{b+ZCRIT*se:+.5f}) p={p:.3g}")'),
    (S3 / "manuscript_numbers.py",
     '    print("\\n### Fixed-effect MVMR (t with df=91)")',
     '    print("\\n### Fixed-effect MVMR (two-sided normal P)")'),
    (S3 / "manuscript_numbers.py",
     '            tcrit = stats.t.ppf(0.975, 91)\n'
     '            p = 2 * stats.t.sf(abs(sub["beta"] / sub["se"]), 91)\n'
     '            print(f"  {oc:16s} FE   {expo:5s} beta={sub[\'beta\']:+.4f} se={sub[\'se\']:.4f} "\n'
     '                  f"ci=({sub[\'beta\']-tcrit*sub[\'se\']:+.4f},{sub[\'beta\']+tcrit*sub[\'se\']:+.4f}) p={p:.3g}")',
     '            p = 2 * stats.norm.sf(abs(sub["beta"] / sub["se"]))\n'
     '            print(f"  {oc:16s} FE   {expo:5s} beta={sub[\'beta\']:+.4f} se={sub[\'se\']:.4f} "\n'
     '                  f"ci=({sub[\'beta\']-ZCRIT*sub[\'se\']:+.4f},{sub[\'beta\']+ZCRIT*sub[\'se\']:+.4f}) p={p:.3g}")'),
    (S3 / "manuscript_numbers.py",
     'SD_TMFI = 2.030884204160846  # from Stage 2 exposure SD, reused verbatim',
     'SD_TMFI = 2.030884204160846  # from Stage 2 exposure SD, reused verbatim\n'
     'ZCRIT   = 1.959964           # normal quantile; matches the reported intervals'),
]


def main():
    n = 0
    for path, old, new in EDITS:
        if not path.exists():
            raise SystemExit(f"missing {path}")
        t = path.read_text(encoding="utf-8")
        if t.count(old) != 1:
            raise SystemExit(f"ABORT: {t.count(old)} matches in {path.name} for\n"
                             f"  {old[:80]}")
        path.write_text(t.replace(old, new), encoding="utf-8")
        n += 1
        print("patched", path.name)
    print(f"\n{n} edits applied")


if __name__ == "__main__":
    main()
