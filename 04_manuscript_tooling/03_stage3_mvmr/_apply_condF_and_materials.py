"""Unify conditional F statistics + fix figure numbering / cross-references / deposited wording.

Every replacement carries a uniqueness assertion so a partially-edited manuscript
can never be silently double-patched.
"""
from pathlib import Path

P = Path(r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md")
t = P.read_text(encoding="utf-8")

REPL = [
    # ------------------------------------------------------------------
    # 1. Conditional F: single reported value (R MVMR), closed-form -> S11
    # ------------------------------------------------------------------
    ("(35.23 from the R `MVMR` package, 42.8 from the closed-form approximation)",
     "(35.23)", 1),

    ("Conditional instrument strength was assessed with multivariable conditional F statistics. "
     "Two implementations were used, and both are reported rather than reconciled, because they "
     "answer slightly different questions. The primary analysis (Table 5) uses the conditional F "
     "returned by the R `MVMR` package, which is computed from the multivariable model's own "
     "residual variance. The site-matched secondary analysis (Table 8) uses the closed-form "
     "approximation F = mean F \u00d7 (1 \u2212 R\u00b2), where mean F is the mean univariable F across the "
     "instrument and R\u00b2 is the proportion of exposure variance explained by the other exposure in "
     "the joint model. The two differ in weighting and in how the residual variance is defined, so "
     "their absolute values are not interchangeable; both are far above the conventional threshold "
     "of 10 and lead to the same conclusion, and the corresponding values are quoted side by side "
     "in the table notes.",
     "Conditional instrument strength was assessed with multivariable conditional F statistics "
     "computed by the R `MVMR` package from the joint model's own residual variance [5]; these are "
     "the values reported throughout (Tables 5 and 8), so that a single instrument-strength "
     "statement applies to every analysis. A closed-form approximation, F = mean F \u00d7 (1 \u2212 R\u00b2), "
     "where mean F is the mean univariable F across the instrument and R\u00b2 is the proportion of "
     "exposure variance explained by the other exposure in the joint model, was computed "
     "independently and is reported as a sensitivity analysis in Supplementary Table S11: it "
     "returns 23.6 for TMFI and 42.8 for BMI against 22.12 and 35.23 from `MVMR`. The two "
     "estimators weight the instruments differently and define the residual variance differently, "
     "so their absolute values are not interchangeable, but both are far above the conventional "
     "threshold of 10 and lead to the same conclusion.", 1),

    ("Conditional F statistics here are the closed-form approximation F = mean F \u00d7 (1 \u2212 R\u00b2), "
     "giving 23.6 for TMFI and 42.8 for BMI; the R `MVMR` package returns 22.12 and 35.23 for the "
     "same instrument set (Table 5). The two are not numerically comparable because they differ in "
     "weighting and in the definition of the residual variance, and this is the only diagnostic in "
     "the manuscript whose value depends on the implementation rather than on the data. Both "
     "exceed the threshold of 10 by a wide margin.",
     "Conditional F statistics are those reported in Table 5 (R `MVMR`: 22.12 for TMFI, 35.23 for "
     "BMI); the closed-form approximation returning 23.6 and 42.8 is reported as a sensitivity "
     "analysis in Supplementary Table S11. Both exceed the threshold of 10 by a wide margin.", 1),

    # Table 8 Panel A rows
    ("| Stage 2 | Knee pain (UKB) | TMFI | 1.0146 | 0.9640\u20131.0678 | 0.58 | 93 | 23.6 |",
     "| Stage 2 | Knee pain (UKB) | TMFI | 1.0146 | 0.9640\u20131.0678 | 0.58 | 93 | 22.12 |", 1),
    ("| Stage 2 | Knee pain (UKB) | BMI | 1.2751 | 1.1772\u20131.3812 | 2.5\u00d710\u207b\u2079 | 93 | 42.8 |",
     "| Stage 2 | Knee pain (UKB) | BMI | 1.2751 | 1.1772\u20131.3812 | 2.5\u00d710\u207b\u2079 | 93 | 35.23 |", 1),
    ("| Stage 2 | Back pain (UKB) | TMFI | 1.0005 | 0.9539\u20131.0495 | 0.98 | 93 | 23.6 |",
     "| Stage 2 | Back pain (UKB) | TMFI | 1.0005 | 0.9539\u20131.0495 | 0.98 | 93 | 22.12 |", 1),
    ("| Stage 2 | Back pain (UKB) | BMI | 1.1344 | 1.0528\u20131.2223 | 9.4\u00d710\u207b\u2074 | 93 | 42.8 |",
     "| Stage 2 | Back pain (UKB) | BMI | 1.1344 | 1.0528\u20131.2223 | 9.4\u00d710\u207b\u2074 | 93 | 35.23 |", 1),

    # ------------------------------------------------------------------
    # 2. Figure numbering: second "Figure 3" -> Figure 4; LOO cross-ref 2 -> 3
    # ------------------------------------------------------------------
    ("### Figure 3. Influence-point diagnostics and robustness envelope.",
     "### Figure 4. Influence-point diagnostics and robustness envelope.", 1),
    ("| Leave-one-out estimate range | 0.0035 to 0.0177 | \u22120.0036 to 0.0168 | see Figure 2 | see Figure 2 |",
     "| Leave-one-out estimate range | 0.0035 to 0.0177 | \u22120.0036 to 0.0168 | see Figure 3 | see Figure 3 |", 1),

    # ------------------------------------------------------------------
    # 3. deposited: unify on "will be deposited at submission"
    # ------------------------------------------------------------------
    ("The following files accompany the manuscript and are deposited with the code archive.",
     "The following files accompany the manuscript and will be deposited with the code archive at "
     "the time of submission (Code availability).", 1),

    # ------------------------------------------------------------------
    # 4. New supplementary item S11
    # ------------------------------------------------------------------
    ("11. **Supplementary Table S10** \u2014 Independent re-implementation of the site-matched estimators: "
     "complete Stage 5 output from the Python and base-R implementations side by side, with the "
     "frozen harmonized inputs they share (`stage5_mr_results.csv`, `stage5_baser_all.csv`, "
     "`r_input/R_stage1_*.tsv`, `r_input/R_stage2_*.tsv`).",
     "11. **Supplementary Table S10** \u2014 Independent re-implementation of the site-matched estimators: "
     "complete Stage 5 output from the Python and base-R implementations side by side, with the "
     "frozen harmonized inputs they share (`stage5_mr_results.csv`, `stage5_baser_all.csv`, "
     "`r_input/R_stage1_*.tsv`, `r_input/R_stage2_*.tsv`).\n"
     "12. **Supplementary Table S11** \u2014 Conditional instrument strength under two estimators: the R "
     "`MVMR` conditional F reported in Tables 5 and 8, and the closed-form approximation "
     "F = mean F \u00d7 (1 \u2212 R\u00b2) reported as a sensitivity analysis "
     "(`stage3_mvmr_diagnostics.csv`, `R_validation_strength.csv`).", 1),
]

for old, new, exp in REPL:
    c = t.count(old)
    if c != exp:
        raise SystemExit("ABORT (%d != %d) for:\n  %s" % (c, exp, old[:120]))
    t = t.replace(old, new)

P.write_text(t, encoding="utf-8")
print("applied %d replacements" % len(REPL))
