# -*- coding: utf-8 -*-
"""Small-fixes batch: title hedge, abstract readability, two-stage clarification, Table 7 note."""
import io
import os
import re
import sys

P = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"
raw = io.open(P, "rb").read().decode("utf-8")
crlf = raw.count("\r\n") > 0
txt = raw.replace("\r\n", "\n")

EDITS = [
    # --- 1. title: remove the unqualified claim -------------------------------------------
    ("T1 title hedge",
     "# Body mass index, but not thigh muscle fat infiltration, shows a direct genetic association "
     "with chronic pain intensity: a multivariable Mendelian randomization study with explicit power bounds",
     "# Body mass index, but not thigh muscle fat infiltration, shows a direct genetic association "
     "with chronic pain intensity within the resolution of this design: a multivariable Mendelian "
     "randomization study with explicit power bounds"),

    # --- 2. abstract background: shorter, keep the comparator ------------------------------
    ("T2 abstract background",
     "Thigh muscle fat infiltration (TMFI) and body mass index (BMI) may both influence chronic pain "
     "intensity (CPI), but their effects are difficult to separate. BMI was pre-specified as an internal "
     "comparator, so that a null TMFI estimate could be read against an explicit power bound.",
     "Thigh muscle fat infiltration (TMFI) and body mass index (BMI) may both influence chronic pain "
     "intensity (CPI) but are hard to separate. BMI was pre-specified as an internal comparator, so a "
     "null TMFI estimate could be read against an explicit power bound."),

    # --- 3. abstract methods: non-specialist entry, drop 'Two-stage MR' and CPI_no_imaging --
    ("T3 abstract methods",
     "Two-stage Mendelian randomization using summary statistics for TMFI (n=46,347), BMI (GIANT 2015) "
     "and two UK Biobank CPI outcomes (n=134,627 and 106,697). Stage 1 estimated the TMFI total effect "
     "on 42 clumped variants; Stage 2 fitted two-exposure multivariable MR on 93 variants with "
     "MVMR-Egger, conditional F statistics, influence diagnostics and targeted exclusions. Power was "
     "expressed as the minimum detectable effect. Shared participants were quantified by cross-trait LD "
     "score regression and propagated to every reported estimate.",
     "We used summary statistics for TMFI (n=46,347), BMI (GIANT 2015) and two UK Biobank chronic pain "
     "intensity outcomes (n=134,627 and 106,697). We first estimated the total effect of TMFI using 42 "
     "independent genetic variants, then estimated TMFI and BMI jointly using 93 so each effect is "
     "adjusted for the other, with standard diagnostics for instrument strength and pleiotropy. Because "
     "exposure and outcomes come from the same cohort, that overlap was quantified genetically and "
     "carried through every estimate; results are reported against the smallest effect detectable at 80% "
     "power rather than significance."),

    # --- 4. abstract results: plain language, drop the code-like label and the ceiling jargon
    ("T4 abstract results",
     "No TMFI total effect emerged (0.008 SD per SD, 95% CI \u22120.018 to 0.033, P=0.558; 0.005, "
     "\u22120.022 to 0.032, P=0.721). Stage 2 gave a positive BMI direct effect (0.090 SD per SD, 95% CI "
     "0.038\u20130.143, P=0.00075) against a null TMFI effect (\u22120.002, \u22120.036 to 0.031, "
     "P=0.898), with the same pattern for CPI_no_imaging; conditional F statistics were 22.1 and 35.2, "
     "with no evidence of directional pleiotropy. Power reached 80% for a TMFI effect of 0.048 SD, "
     "1.9-fold smaller than the BMI effect recovered. Overlap was real (cross-trait intercept 0.052, "
     "P=1.2\u00d710\u207b\u00b9\u00b9) but moved estimates by \u22640.0015 SD per SD under the evaluated "
     "scenarios; the intercept needed to reverse either finding exceeded the structural ceiling attainable "
     "under any sharing of participants. Two site-matched secondary outcomes also returned a TMFI null "
     "(knee pain OR 1.015 per SD, P=0.58; back pain 1.001, P=0.98) with a larger BMI effect at the "
     "weight-bearing joint (1.275 vs 1.134).",
     "No TMFI total effect emerged (0.008 SD per SD, 95% CI \u22120.018 to 0.033, P=0.558; 0.005, "
     "\u22120.022 to 0.032, P=0.721). The joint model gave a positive BMI direct effect (0.090 SD per SD, "
     "95% CI 0.038\u20130.143, P=0.00075) against a null TMFI effect (\u22120.002, \u22120.036 to 0.031, "
     "P=0.898), with the same pattern for the narrower pain definition; instruments remained strong after "
     "conditioning (F 22.1 and 35.2) and showed no directional pleiotropy. The study could detect a TMFI "
     "effect of 0.048 SD at 80% power, 1.9-fold smaller than the BMI effect recovered. Exposure and "
     "outcomes share participants (intercept 0.052, P=1.2\u00d710\u207b\u00b9\u00b9); this moved estimates "
     "by \u22640.0015 SD per SD, and no degree of sharing could reverse either finding. Two site-specific "
     "secondary outcomes also returned a TMFI null (knee pain OR 1.015 per SD, P=0.58; back pain 1.001, "
     "P=0.98), with a descriptively larger BMI effect at the weight-bearing joint (1.275 vs 1.134)."),

    # --- 5. abstract conclusions: tighter --------------------------------------------------
    ("T5 abstract conclusions",
     "Genetically predicted adiposity, but not thigh muscle fat infiltration, is associated with chronic "
     "pain intensity. With an explicit power bound, an internal comparator that behaved as expected and "
     "quantified sample overlap, this is a bounded null: effects above roughly 0.05 SD per SD are "
     "unlikely, smaller ones remain open.",
     "Genetically predicted adiposity, but not thigh muscle fat infiltration, is associated with chronic "
     "pain intensity. Against an explicit power bound, an internal comparator that behaved as expected "
     "and quantified participant overlap, this is a bounded null: effects above about 0.05 SD per SD are "
     "unlikely; smaller ones remain open."),

    # --- 6. two-stage clarification (R3-M6) -------------------------------------------------
    ("T6 two-stage clarification",
     "We addressed both weaknesses in a two-stage design.",
     "We addressed both weaknesses in two analysis stages applied to the same summary data \u2014 a stage "
     "here is an analysis step, not an independent sample, and the design is not a two-stage "
     "least-squares fit."),

    # --- 7. Table 7 note: the 'conservative' claim is directionally wrong --------------------
    ("T7 Table 7 note",
     "Bias/SE for Stage 1 is quoted against the standard errors reported in Table 2 and for Stage 2 "
     "against the fixed-effect standard errors returned by the IVW algebra, which are the smaller of the "
     "two choices and therefore the conservative one.",
     "Bias/SE for Stage 1 is quoted against the standard errors reported in Table 2 and for Stage 2 "
     "against the fixed-effect standard errors returned by the IVW algebra. The tipping-point intercept "
     "is proportional to the standard error used, so the smaller fixed-effect errors give the smaller "
     "\u2014 that is, the less conservative \u2014 tipping point; recomputing Panel C with the larger "
     "multiplicative-random-effects errors of Table 3 gives 1.271, 1.242, 1.028 and 1.276 for the four "
     "paths in the order listed, each still above its structural ceiling (by factors of 1.57, 1.18, 1.16 "
     "and 1.39). Every scenario in this table therefore falls on the same side of its ceiling under "
     "either choice of standard error."),
]

fail = 0
for name, old, new in EDITS:
    n = txt.count(old)
    if n != 1:
        print("FAIL %-28s occurrences=%d (expected 1)" % (name, n))
        fail += 1
        continue
    txt = txt.replace(old, new)
    print("ok   %-28s %+d chars" % (name, len(new) - len(old)))

if fail:
    print("\nABORT: nothing written.")
    sys.exit(1)

# abstract ceiling check
lines_new = txt.split("\n")


def nw(s):
    s = re.sub(r"^#+\s*", "", s)
    return len([w for w in s.split() if any(c.isalnum() for c in w)])


body = nw(lines_new[10]) + nw(lines_new[14]) + nw(lines_new[18]) + nw(lines_new[22])
print("\nabstract body = %d words (ceiling 350); incl. keyword line = %d" % (body, body + nw(lines_new[24])))
if body > 345:
    print("ABORT: abstract too long.")
    sys.exit(1)

out = txt.replace("\n", "\r\n") if crlf else txt
bak = P + ".bak_2026-09-20c"
if not os.path.isfile(bak):
    io.open(bak, "wb").write(raw.encode("utf-8"))
io.open(P, "wb").write(out.encode("utf-8"))
print("written", P, "CRLF=", crlf, "| backup", os.path.basename(bak))
