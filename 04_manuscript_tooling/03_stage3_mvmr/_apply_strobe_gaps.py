"""Add the four STROBE-MR items the manuscript does not yet address:
6d missing data, 6e multiple testing, 9b preregistration, 17 generalizability."""
from pathlib import Path

P = Path(r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md")
t = P.read_text(encoding="utf-8")

ANCHOR_CI = ("a difference below 0.003. Raw-scale coefficients are reported alongside all "
             "standardized values.")

ADD_CI = """
Missing data were handled by complete-case harmonization. A variant was carried into an analysis only if the exposure file, the outcome file and the LD reference panel all supplied a usable allele pair, frequency, effect and standard error for it; variants failing any of these were dropped from that analysis rather than imputed, and the counts resulting from each step are given under Instrument selection and harmonization and in Supplementary Table S1. No multiple-testing correction was applied: the primary contrast is a single pre-specified comparison between two exposures fitted in the same model, the site-matched secondary outcomes were pre-specified as a replication of that contrast rather than as new hypotheses, and the reported P values are interpreted as measures of compatibility with the null rather than as pass/fail decisions at a threshold. The study was not preregistered; the instrument sets, outcome definitions, estimator hierarchy and sensitivity analyses were fixed before the Stage 2 estimates were inspected, and the frozen intermediate files that record this are deposited with the code archive, but no protocol was posted to a public registry."""

ANCHOR_LIM = "Seven boundaries qualify the estimates."

REPL = [
    (ANCHOR_CI, ANCHOR_CI + ADD_CI, 1),
]

# generalizability: append a sentence to the final limitation item if a suitable anchor exists.
GEN_ANCHOR = ("Seventh, the conversion of the site-matched risk differences to the log-odds scale "
              "assumes a constant case fraction across genotypes")
if t.count(GEN_ANCHOR) != 1:
    GEN_ANCHOR = None

for old, new, exp in REPL:
    c = t.count(old)
    if c != exp:
        raise SystemExit("ABORT (%d != %d) for:\n  %s" % (c, exp, old[:100]))
    t = t.replace(old, new, 1)

P.write_text(t, encoding="utf-8")
print("ok; generalizability anchor found:", GEN_ANCHOR is not None)
