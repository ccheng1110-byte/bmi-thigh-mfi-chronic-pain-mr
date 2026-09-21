"""Apply the P1 wording revisions + two stale-number fixes to the manuscript.

Each replacement must occur exactly once; the script aborts otherwise, so a
silent no-op or an over-broad edit cannot slip through.
"""

from pathlib import Path

P = Path(r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md")

REPL = [
    # --- P1-1 absolute-claim about sample overlap -------------------------
    (
        "conditional F statistics were 22.1 and 35.2, with no directional pleiotropy.",
        "conditional F statistics were 22.1 and 35.2, with no evidence of directional pleiotropy.",
    ),
    (
        "Overlap was real (cross-trait intercept 0.052, P=1.2×10⁻¹¹) but moved estimates by "
        "≤0.0015 SD per SD; reversing either finding would require more overlap than any shared "
        "sample can generate.",
        "Overlap was real (cross-trait intercept 0.052, P=1.2×10⁻¹¹) but moved estimates by "
        "≤0.0015 SD per SD under the evaluated scenarios; the intercept needed to reverse either "
        "finding exceeded the structural ceiling attainable under any sharing of participants.",
    ),
    (
        "### Sample overlap was real but could not account for either finding",
        "### Sample overlap was real but left both findings unchanged under the evaluated scenarios",
    ),
    (
        "The absence of exposure-side overlap additionally excludes correlated measurement error "
        "across the columns of the multivariable design as a separate bias source.",
        "A negative intercept estimate does not establish that the true overlap is zero, and these "
        "estimates carry their own uncertainty; they do, however, argue against correlated "
        "measurement error across the columns of the multivariable design as an additional bias source.",
    ),
    (
        "The conclusion is therefore stronger than insensitivity: no degree of UK Biobank sample "
        "overlap, even complete overlap with perfectly correlated phenotypes, could have generated "
        "enough sampling covariance to reverse either Stage 2 finding.",
        "This ceiling argument is a sensitivity analysis under stated assumptions rather than a "
        "proof: it takes the intercepts as estimated without error and treats the reported standard "
        "errors as the relevant reversal threshold. Under those assumptions, no configuration of "
        "overlapping participants — including complete overlap with perfectly correlated phenotypes "
        "— attained an intercept large enough to reverse either Stage 2 finding.",
    ),
    (
        "In every case the required intercept exceeds the ceiling, so no configuration of "
        "overlapping participants could have reversed either Stage 2 finding.",
        "In every case the required intercept exceeds the ceiling. Under the assumptions of this "
        "sensitivity analysis — intercepts taken as estimated and the reported standard errors used "
        "as the reversal threshold — no configuration of overlapping participants attained an "
        "intercept large enough to reverse either Stage 2 finding.",
    ),
    (
        "and neither result could be produced by sample overlap. Sample overlap, the one concern "
        "that previously could be raised only qualitatively against these analyses, has been "
        "quantified here and shown to be incapable of reversing either finding, so it need not be "
        "the focus of that replication.",
        "and neither result changed materially under the overlap scenarios evaluated. Sample "
        "overlap, the one concern that previously could be raised only qualitatively against these "
        "analyses, has been quantified here; under the evaluated scenarios it did not come close to "
        "reversing either finding, so it need not be the focus of that replication.",
    ),

    # --- P1-2 MDE read as an exclusion boundary ---------------------------
    (
        "These data are therefore inconsistent with a TMFI effect of the same order as the "
        "adiposity effect and are compatible only with effects small enough to be of limited "
        "phenotypic consequence.",
        "These data are therefore inconsistent with a TMFI effect of the same order as the "
        "adiposity effect: the confidence intervals exclude effects of that magnitude, while "
        "smaller effects remain compatible with the data.",
    ),

    # --- P1-3 diagnostics written as validity proofs ----------------------
    (
        "the MVMR-Egger intercepts were indistinguishable from zero, so that heterogeneity is not "
        "unbalanced directional pleiotropy; the fixed-effect fit reproduced both point estimates, "
        "so heterogeneity inflates variance rather than shifting estimates; and conditional F "
        "statistics of 22.12 and 35.23 exclude weak-instrument bias as an explanation for either "
        "result.",
        "the MVMR-Egger intercepts were indistinguishable from zero, providing no evidence of "
        "average directional pleiotropy although pleiotropic bias cannot be excluded; the "
        "fixed-effect fit reproduced both point estimates, so heterogeneity inflates variance "
        "rather than shifting estimates; and conditional F statistics of 22.12 and 35.23 reduce "
        "concern about severe conditional weak-instrument bias without excluding it.",
    ),
    (
        "A Stage 2 model was therefore not needed to establish that the TMFI total effect is null;",
        "A Stage 2 model was therefore not needed to show that the data contain no detectable TMFI "
        "total effect;",
    ),
    (
        "while genetically predicted thigh muscle fat infiltration showed no effect, either as a "
        "total effect or as an effect conditional on BMI. The contribution of this paper is not "
        "that it separates one signal from another — Stage 1 already established that the TMFI "
        "total effect is null — but that it converts that null into a quantitative statement.",
        "while genetically predicted thigh muscle fat infiltration showed no evidence of an effect, "
        "either as a total effect or as an effect conditional on BMI. The contribution of this "
        "paper is not that it separates one signal from another — Stage 1 already showed no "
        "detectable TMFI total effect — but that it converts that null into a quantitative statement.",
    ),
    (
        "whereas genetically predicted thigh muscle fat infiltration showed neither a total nor a "
        "conditional direct effect.",
        "whereas genetically predicted thigh muscle fat infiltration showed no evidence of a total "
        "or conditional direct effect.",
    ),
    (
        "genetically predicted thigh muscle fat infiltration showed no effect on knee or back pain,",
        "genetically predicted thigh muscle fat infiltration showed no evidence of an effect on "
        "knee or back pain,",
    ),

    # --- P1-4 secondary outcome called an independent replication ---------
    (
        "An independent knee-pain GWAS from the same cohort (429,737 participants [13]) was used "
        "to replicate the univariable TMFI estimate only.",
        "A second knee-pain GWAS from the same cohort, with a different phenotyping and analysis "
        "pipeline (429,737 participants [13]), was used as a same-cohort concordance check on the "
        "univariable TMFI estimate only; because it draws on the same participants, it tests "
        "specification rather than providing independent replication.",
    ),
    (
        "An independent knee-pain GWAS reproduced the Stage 1 null (0.9965, 0.9561–1.0386, P=0.87).",
        "The second knee-pain GWAS, which shares participants with the first and therefore tests "
        "specification rather than sampling, returned a concordant Stage 1 null (0.9965, "
        "0.9561–1.0386, P=0.87).",
    ),
    (
        "plus an independent knee-pain replication source",
        "plus a second knee-pain GWAS specification from the same cohort",
    ),
    (
        "| Stage 1 | Knee pain (replication; 429,737) | TMFI | 0.9965 |",
        "| Stage 1 | Knee pain (2nd specification; 429,737) | TMFI | 0.9965 |",
    ),
    (
        "| Knee pain (replication) | 1.0609 | not interpreted |",
        "| Knee pain (2nd specification) | 1.0609 | not interpreted |",
    ),
    (
        "1.0179–1.0201 for the replication source",
        "1.0179–1.0201 for the second knee-pain specification",
    ),

    # --- P1-5 genetic correlation read as BMI-mediated / adiposity-derived -
    (
        "That this correlation is not causal is shown by the multivariable model, in which "
        "conditioning on BMI — itself correlated with knee pain at r_g = 0.433 — removes the TMFI "
        "association entirely. The site-specific genetic signal is therefore real but "
        "adiposity-derived.",
        "The multivariable model shows that this correlation does not translate into a detectable "
        "direct effect: conditioning on BMI — itself correlated with knee pain at r_g = 0.433 — "
        "left no detectable TMFI association. The shared genetic architecture with the "
        "site-specific outcomes is therefore real, but it did not yield a detectable direct TMFI "
        "association in these models, and the design does not identify what generates the correlation.",
    ),

    # --- P1-6 reviewer-facing / defensive phrasing ------------------------
    (
        "One estimator did differ, and we report it rather than conceal it.",
        "One estimator did differ between the two implementations.",
    ),

    # --- garbled sentence in the Figure 3 legend --------------------------
    (
        "while TMFI estimates remain within the range excluded by neither significance nor the "
        "power bound.",
        "while TMFI estimates remain close to zero and within the confidence intervals reported in "
        "Table 3.",
    ),

    # --- stale numbers left behind by the SD_BMI correction ----------------
    (
        "with a larger BMI effect at the weight-bearing joint (1.349 vs 1.168).",
        "with a larger BMI effect at the weight-bearing joint (1.275 vs 1.134).",
    ),
    (
        "against +5.03 (+3.37 to +6.68) and +2.96 (+1.21 to +4.71) points per SD of BMI.",
        "against +4.08 (+2.74 to +5.43) and +2.41 (+0.98 to +3.83) points per SD of BMI.",
    ),
]


def main():
    t = P.read_text(encoding="utf-8")
    n = 0
    for old, new in REPL:
        c = t.count(old)
        if c != 1:
            raise SystemExit(f"ABORT: {c} occurrences of:\n  {old[:110]}")
        t = t.replace(old, new)
        n += 1
    P.write_text(t, encoding="utf-8")
    print(f"applied {n} replacements")


if __name__ == "__main__":
    main()
