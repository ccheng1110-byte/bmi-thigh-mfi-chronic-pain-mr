# -*- coding: utf-8 -*-
"""Add the STROBE-MR reporting statement, supplementary checklist entry and reference."""
from pathlib import Path

P = Path(r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md")
AP = chr(39)  # straight apostrophe

R15_TAIL = ("15. Bycroft C, Freeman C, Petkova D, Band G, Elliott LT, Sharp K, Motyer A, "
            "Vukcevic D, Delaneau O, O%sConnell J, Cortes A, Welsh S, Young A, Effingham M, "
            "McVean G, Leslie S, Allen N, Donnelly P, Marchini J. The UK Biobank resource with "
            "deep phenotyping and genomic data. *Nature*. 2018;562:203\u2013209. "
            "doi:[10.1038/s41586-018-0579-z](https://doi.org/10.1038/s41586-018-0579-z)." % AP)

R16 = ("16. Skrivankova VW, Richmond RC, Woolf BAR, Davies NM, Swanson SA, VanderWeele TJ, "
       "et al. Strengthening the reporting of observational studies in epidemiology using "
       "Mendelian randomization: the STROBE-MR statement. *JAMA*. 2021;326(16):1614\u20131621. "
       "doi:[10.1001/jama.2021.18236](https://doi.org/10.1001/jama.2021.18236).")

REPL = [
    ("The study was not preregistered: instrument sets, outcome definitions, estimator hierarchy "
     "and sensitivity analyses were fixed before the Stage 2 estimates were inspected and the "
     "frozen intermediate files recording this are deposited with the code archive, but no "
     "protocol was posted to a public registry.",
     "The study was not preregistered: instrument sets, outcome definitions, estimator hierarchy "
     "and sensitivity analyses were fixed before the Stage 2 estimates were inspected and the "
     "frozen intermediate files recording this are deposited with the code archive, but no "
     "protocol was posted to a public registry. The study is reported in accordance with the "
     "STROBE-MR statement [16]; the completed 20-item checklist, with every item mapped to its "
     "location in this manuscript, is supplied as Supplementary Table S12."),
    ("12. **Supplementary Table S11** \u2014 Conditional instrument strength under two estimators: "
     "the R `MVMR` conditional F reported in Tables 5 and 8, and the closed-form approximation "
     "F = mean F \u00d7 (1 \u2212 R\u00b2) reported as a sensitivity analysis "
     "(`stage3_mvmr_diagnostics.csv`, `R_validation_strength.csv`).",
     "12. **Supplementary Table S11** \u2014 Conditional instrument strength under two estimators: "
     "the R `MVMR` conditional F reported in Tables 5 and 8, and the closed-form approximation "
     "F = mean F \u00d7 (1 \u2212 R\u00b2) reported as a sensitivity analysis "
     "(`stage3_mvmr_diagnostics.csv`, `R_validation_strength.csv`).\n"
     "13. **Supplementary Table S12** \u2014 STROBE-MR reporting checklist: all 20 items with the "
     "manuscript section and line reference satisfying each one."),
    (R15_TAIL, R15_TAIL + "\n" + R16),
]


def main():
    t = P.read_text(encoding="utf-8")
    bad = [(i, t.count(o), o[:70]) for i, (o, n) in enumerate(REPL) if t.count(o) != 1]
    if bad:
        print("VERIFY FAILED - nothing written")
        for b in bad:
            print("   #%d count=%d :: %s" % b)
        raise SystemExit(1)
    for o, n in REPL:
        t = t.replace(o, n)
    P.write_text(t, encoding="utf-8")
    print("ok; STROBE mentions =", t.count("STROBE"), "; ref16 =", t.count("jama.2021.18236"))


if __name__ == "__main__":
    main()
