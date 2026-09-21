#!/usr/bin/env python
"""Apply the SD_BMI (=1.00) and P0-3 (normal P) corrections to the manuscript.

Every replacement is asserted to occur exactly once, so a stale draft fails
loudly instead of being silently half-updated.
"""
from pathlib import Path

MS = Path(r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr"
          r"\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md")

REPL = [
    # ---------------- abstract ----------------
    ("Stage 2 gave a positive BMI direct effect (0.111 SD per SD, "
     "95% CI 0.047–0.176, P=0.00075)",
     "Stage 2 gave a positive BMI direct effect (0.090 SD per SD, "
     "95% CI 0.038–0.143, P=0.00075)"),
    ("Power reached 80% for a TMFI effect of 0.048 SD, 2.3-fold smaller "
     "than the BMI effect recovered.",
     "Power reached 80% for a TMFI effect of 0.048 SD, 1.9-fold smaller "
     "than the BMI effect recovered."),
    ("moved estimates by ≤0.001 SD per SD;",
     "moved estimates by ≤0.0015 SD per SD;"),

    # ---------------- Results ----------------
    ("corresponding to 0.111 SD of pain intensity per SD of BMI "
     "(95% CI 0.047–0.176)",
     "corresponding to 0.090 SD of pain intensity per SD of BMI "
     "(95% CI 0.038–0.143)"),
    ("or 0.125 SD per SD (95% CI 0.054–0.196)",
     "or 0.101 SD per SD (95% CI 0.044–0.159)"),
    ("the BMI direct effect moved by at most 0.005 SD per SD in either "
     "direction",
     "the BMI direct effect moved by at most 0.004 SD per SD in either "
     "direction"),
    ("BMI gave an odds ratio per SD of 1.3487 (1.2223–1.4881, P=2.5×10⁻⁹) "
     "for knee pain and 1.1679 (1.0653–1.2803, P=9.4×10⁻⁴) for back pain",
     "BMI gave an odds ratio per SD of 1.2751 (1.1772–1.3812, P=2.5×10⁻⁹) "
     "for knee pain and 1.1344 (1.0528–1.2223, P=9.4×10⁻⁴) for back pain"),
    ("approximately one quarter and two fifths of the BMI effect that the "
     "same model recovered",
     "approximately 27% and 48% of the BMI effect that the same model "
     "recovered"),
    ("+5.03 (+3.37 to +6.68) and +2.96 (+1.21 to +4.71) for BMI",
     "+4.08 (+2.74 to +5.43) and +2.40 (+0.98 to +3.83) for BMI"),

    # ---------------- Discussion / Conclusion ----------------
    ("recovered a BMI effect of 0.111 SD per SD with approximately 93% "
     "power",
     "recovered a BMI effect of 0.090 SD per SD with approximately 92% "
     "power"),
    ("recovered the BMI effect at approximately 93% power while retaining "
     "80% power for a TMFI effect of 0.048 SD per SD",
     "recovered the BMI effect at approximately 92% power while retaining "
     "80% power for a TMFI effect of 0.048 SD per SD"),
    ("(OR per SD 0.884, 95% CI 0.807–0.968, P=0.008)",
     "(OR per SD 0.905, 95% CI 0.840–0.974, P=0.008)", 2),
    ("a rigorous correction would require the per-cohort λ.",
     "a rigorous correction would require the per-cohort λ. The same "
     "correction inflates the released standard errors by about 27%, which "
     "is why the exposure SD for BMI is set to 1.00 by definition rather "
     "than taken from the implied-SD estimator (Methods), and why the "
     "conditional F statistics reported for BMI are conservative by a "
     "factor of about 1.6."),

    # ---------------- Table 3 ----------------
    ("| IVW multiplicative random effects (primary) | CPI_full | BMI | "
     "0.2610 | 0.0775 | 0.1092 to 0.4128 | 0.0011 | 0.1113 |",
     "| IVW multiplicative random effects (primary) | CPI_full | BMI | "
     "0.2610 | 0.0775 | 0.1092 to 0.4128 | 0.00075 | 0.0904 |"),
    ("| IVW multiplicative random effects (primary) | CPI_no_imaging | BMI | "
     "0.2972 | 0.0864 | 0.1279 to 0.4665 | 0.00088 | 0.1248 |",
     "| IVW multiplicative random effects (primary) | CPI_no_imaging | BMI | "
     "0.2972 | 0.0864 | 0.1279 to 0.4665 | 0.00058 | 0.1014 |"),
    ("| IVW multiplicative random effects (primary) | CPI_no_imaging | TMFI | "
     "−0.0121 | 0.0272 | −0.0655 to 0.0413 | 0.657 | −0.0084 |",
     "| IVW multiplicative random effects (primary) | CPI_no_imaging | TMFI | "
     "−0.0121 | 0.0272 | −0.0655 to 0.0413 | 0.656 | −0.0084 |"),
    ("| IVW fixed effects | CPI_full | BMI | 0.2610 | 0.0563 | "
     "0.1507 to 0.3712 | 1.2×10⁻⁵ | 0.1113 |",
     "| IVW fixed effects | CPI_full | BMI | 0.2610 | 0.0563 | "
     "0.1507 to 0.3712 | 3.5×10⁻⁶ | 0.0904 |"),
    ("| IVW fixed effects | CPI_no_imaging | TMFI | −0.0121 | 0.0202 | "
     "−0.0518 to 0.0275 | 0.550 | −0.0084 |",
     "| IVW fixed effects | CPI_no_imaging | TMFI | −0.0121 | 0.0202 | "
     "−0.0518 to 0.0275 | 0.549 | −0.0084 |"),
    ("| IVW fixed effects | CPI_no_imaging | BMI | 0.2972 | 0.0642 | "
     "0.1714 to 0.4230 | 1.2×10⁻⁵ | 0.1248 |",
     "| IVW fixed effects | CPI_no_imaging | BMI | 0.2972 | 0.0642 | "
     "0.1714 to 0.4230 | 3.6×10⁻⁶ | 0.1014 |"),
    ("| MVMR-Egger slope | CPI_full | BMI | 0.2721 | 0.0793 | "
     "0.1167 to 0.4275 | 0.00091 | 0.1160 |",
     "| MVMR-Egger slope | CPI_full | BMI | 0.2721 | 0.0793 | "
     "0.1167 to 0.4275 | 0.00060 | 0.0943 |"),
    ("| MVMR-Egger slope | CPI_no_imaging | TMFI | −0.0122 | 0.0274 | "
     "−0.0659 to 0.0415 | 0.658 | −0.0084 |",
     "| MVMR-Egger slope | CPI_no_imaging | TMFI | −0.0122 | 0.0274 | "
     "−0.0659 to 0.0415 | 0.657 | −0.0084 |"),
    ("| MVMR-Egger slope | CPI_no_imaging | BMI | 0.2941 | 0.0886 | "
     "0.1204 to 0.4679 | 0.0013 | 0.1235 |",
     "| MVMR-Egger slope | CPI_no_imaging | BMI | 0.2941 | 0.0886 | "
     "0.1204 to 0.4679 | 0.00090 | 0.1003 |"),
    ("| MVMR-Egger intercept | CPI_full | — | −0.0013 | 0.0019 | "
     "−0.0050 to 0.0023 | 0.483 | — |",
     "| MVMR-Egger intercept | CPI_full | — | −0.0013 | 0.0019 | "
     "−0.0050 to 0.0023 | 0.481 | — |"),
    ("| MVMR-Egger intercept | CPI_no_imaging | — | 0.0004 | 0.0021 | "
     "−0.0037 to 0.0044 | 0.864 | — |",
     "| MVMR-Egger intercept | CPI_no_imaging | — | 0.0004 | 0.0021 | "
     "−0.0037 to 0.0044 | 0.863 | — |"),
    ("**Table note:** IVW P values are two-sided Student t with 91 degrees "
     "of freedom; MVMR-Egger P values use 90 degrees of freedom. "
     "Standardized values use the implied SDs listed in Table 6.",
     "**Table note:** P values are two-sided from the standard normal "
     "distribution applied to z = β/SE, matching the ±1.96 × SE intervals. "
     "Standardized values use sd<sub>TMFI</sub>=2.03 and sd<sub>BMI</sub>=1.00 "
     "with the outcome SDs listed in Table 6."),

    # ---------------- Table 6 ----------------
    ("| Stage 2 (MVMR) | CPI_full | BMI | 0.1113 | 0.0465 to 0.1758 | "
     "0.0925 | 0.1070 |",
     "| Stage 2 (MVMR) | CPI_full | BMI | 0.0904 | 0.0378 to 0.1430 | "
     "0.0752 | 0.0870 |"),
    ("| Stage 2 (MVMR) | CPI_no_imaging | BMI | 0.1248 | 0.0537 to 0.1960 | "
     "0.1016 | 0.1175 |",
     "| Stage 2 (MVMR) | CPI_no_imaging | BMI | 0.1014 | 0.0436 to 0.1591 | "
     "0.0826 | 0.0955 |"),
    ("giving 2.03 for TMFI, 1.23 for BMI, 2.89 for CPI_full and 2.93 for "
     "CPI_no_imaging. The observed BMI direct effects correspond to "
     "approximately 93% power on CPI_full and 92% on CPI_no_imaging.",
     "giving 2.03 for TMFI, 2.89 for CPI_full and 2.93 for CPI_no_imaging; "
     "sd<sub>BMI</sub> is 1.00 by definition because the GIANT coefficients "
     "are already per SD (Methods). The observed BMI direct effects "
     "correspond to approximately 92% power on CPI_full and 93% on "
     "CPI_no_imaging."),

    # ---------------- Table 7 ----------------
    ("| Stage 2 (measured intercepts) | CPI_full | BMI | 0.11127 | "
     "−0.00150 | −0.02 | 0.11277 | 2.6×10⁻⁶ |",
     "| Stage 2 (measured intercepts) | CPI_full | BMI | 0.09041 | "
     "−0.00122 | −0.02 | 0.09163 | 2.6×10⁻⁶ |"),
    ("| Stage 2 (measured intercepts) | CPI_no_imaging | BMI | 0.12477 | "
     "−0.00066 | −0.01 | 0.12543 | 3.2×10⁻⁶ |",
     "| Stage 2 (measured intercepts) | CPI_no_imaging | BMI | 0.10135 | "
     "−0.00054 | −0.01 | 0.10191 | 3.2×10⁻⁶ |"),
    ("the corrected BMI direct effect stayed within 0.106 to 0.113 and "
     "0.123 to 0.125",
     "the corrected BMI direct effect stayed within 0.086 to 0.092 and "
     "0.100 to 0.102"),

    # ---------------- Table 8 ----------------
    ("| Stage 2 | Knee pain (UKB) | BMI | 1.3487 | 1.2223–1.4881 | "
     "2.5×10⁻⁹ | 93 | 42.8 |",
     "| Stage 2 | Knee pain (UKB) | BMI | 1.2751 | 1.1772–1.3812 | "
     "2.5×10⁻⁹ | 93 | 42.8 |"),
    ("| Stage 2 | Back pain (UKB) | BMI | 1.1679 | 1.0653–1.2803 | "
     "9.4×10⁻⁴ | 93 | 42.8 |",
     "| Stage 2 | Back pain (UKB) | BMI | 1.1344 | 1.0528–1.2223 | "
     "9.4×10⁻⁴ | 93 | 42.8 |"),
    ("| Knee pain (UKB) | 1.0681 | 1.3487 | 25% |",
     "| Knee pain (UKB) | 1.0681 | 1.2751 | 27% |"),
    ("| Back pain | 1.0621 | 1.1679 | 39% |",
     "| Back pain | 1.0621 | 1.1344 | 48% |"),
]


def main():
    t = MS.read_text(encoding="utf-8")
    n = 0
    for item in REPL:
        old, new = item[0], item[1]
        want = item[2] if len(item) > 2 else 1
        c = t.count(old)
        if c != want:
            raise SystemExit(f"ABORT: {c} occurrences (want {want}) of\n"
                             f"  {old[:90]}...")
        t = t.replace(old, new)
        n += 1
    MS.write_text(t, encoding="utf-8")
    print(f"applied {n} replacements to {MS.name}")


if __name__ == "__main__":
    main()
