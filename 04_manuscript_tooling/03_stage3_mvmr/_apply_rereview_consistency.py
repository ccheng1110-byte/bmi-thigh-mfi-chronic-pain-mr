# -*- coding: utf-8 -*-
"""Apply the four unambiguous consistency corrections identified by the
2026-09-20 round-2 review verification.

Each replacement is asserted to occur exactly once, and the pre-change
occurrence count is printed so the edit can be audited.
"""
import io
import sys

MS = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"

REPLACEMENTS = [
    # 1. Table 5 Egger-intercept P values used the Student-t export although the
    #    manuscript declares the normal rule (Table 3 and Results use 0.481/0.863).
    (
        "| Intercept P value | 0.344 | 0.324 | 0.483 | 0.864 |",
        "| Intercept P value | 0.344 | 0.324 | 0.481 | 0.863 |",
        1,
        "Table 5 Egger intercept P -> declared normal rule",
    ),
    # 2. Back-pain BMI risk difference: log OR 0.126084 / lambda 5.24349 = 2.405 pp,
    #    i.e. 2.40, not 2.41.
    (
        "against +4.08 (+2.74 to +5.43) and +2.41 (+0.98 to +3.83) points per SD of BMI",
        "against +4.08 (+2.74 to +5.43) and +2.40 (+0.98 to +3.83) points per SD of BMI",
        1,
        "Results back-pain BMI risk difference 2.41 -> 2.40",
    ),
    # 3. Table 5 note wrongly said Table 8 reports the closed-form conditional F.
    (
        "Table 8 reports the same quantity under the closed-form approximation F = mean F × (1 − R²), which yields 23.6 and 42.8; the two sets should not be compared numerically because they differ in weighting and in the definition of the residual variance, but both exceed the threshold of 10 by a wide margin and support the same interpretation.",
        "The closed-form approximation F = mean F × (1 − R²), computed independently, returns 23.6 for TMFI and 42.8 for BMI and is reported as a sensitivity analysis in Supplementary Table S11; the two estimators should not be compared numerically because they differ in weighting and in the definition of the residual variance, but both exceed the threshold of 10 by a wide margin and support the same interpretation. Tables 5 and 8 therefore quote the same `MVMR` values.",
        1,
        "Table 5 note: correct the cross-reference to Table 8 / S11",
    ),
    # 4. Discussion "about one quarter" contradicts the 35% and 28% stated in Results.
    (
        "The confidence intervals therefore exclude any positive TMFI effect larger than about one quarter of the BMI effect.",
        "The upper confidence limits therefore exclude any positive TMFI effect larger than 35% of the BMI point estimate for `CPI_full` and 28% for `CPI_no_imaging`.",
        1,
        "Discussion: about one quarter -> 35% / 28%",
    ),
]


def main():
    text = io.open(MS, encoding="utf-8").read()
    for old, new, expect, label in REPLACEMENTS:
        n = text.count(old)
        if n != expect:
            print("ABORT: %s -> expected %d occurrence(s), found %d" % (label, expect, n))
            sys.exit(1)
        text = text.replace(old, new)
        print("applied: %s (occurrences=%d)" % (label, n))
    io.open(MS, "w", encoding="utf-8", newline="").write(text)
    print("written:", MS)


if __name__ == "__main__":
    main()
