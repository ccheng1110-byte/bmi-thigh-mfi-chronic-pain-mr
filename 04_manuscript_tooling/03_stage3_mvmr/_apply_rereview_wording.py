# -*- coding: utf-8 -*-
"""Apply the three agreed re-review revisions to the manuscript.

1. positive control -> internal comparator (BMI) / consistency check (nested subset)
2. same-cohort outcome wording: reproduced/replicating -> same pattern (nested, not independent)
3. soften the L215 heterogeneity direction claim and the L265 pathway / intervention-target claim
"""
import io
import os
import sys

P = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"

raw = io.open(P, "rb").read().decode("utf-8")
crlf = raw.count("\r\n") > 0
txt = raw.replace("\r\n", "\n")

EDITS = []

# ---------------- Group 1: positive control -> comparator / consistency check -------------
EDITS.append((
    "G1-01 abstract background",
    "BMI was pre-specified as an internal positive control, so that a null TMFI estimate could be read against an explicit power bound.",
    "BMI was pre-specified as an internal comparator, so that a null TMFI estimate could be read against an explicit power bound.",
))
EDITS.append((
    "G1-02 abstract conclusion",
    "With an explicit power bound, a recovered positive control and quantified sample overlap,",
    "With an explicit power bound, an internal comparator that behaved as expected and quantified sample overlap,",
))
EDITS.append((
    "G1-03 intro design",
    "treating the BMI coefficient as an internal positive control and reporting study power as the minimum detectable effect",
    "treating the BMI coefficient as an internal comparator and reporting study power as the minimum detectable effect",
))
EDITS.append((
    "G1-04 intro rationale",
    "BMI instead served as an internal positive control: recovering a BMI effect in the same model demonstrates",
    "BMI instead served as an internal comparator: recovering a BMI effect in the same model demonstrates",
))
EDITS.append((
    "G1-05 methods pipeline check",
    "The internal positive control \u2014 CPI_no_imaging is a nested subset of CPI_full and therefore maximally overlapping \u2014 returned",
    "The nested-subset consistency check \u2014 CPI_no_imaging is a subset of CPI_full and therefore maximally overlapping by design \u2014 returned",
))
EDITS.append((
    "G1-06 results power",
    "That bound becomes interpretable against the positive control obtained in the same model.",
    "That bound becomes interpretable against the BMI comparator estimated in the same model.",
))
EDITS.append((
    "G1-07 results overlap",
    "The internal positive control recovered near-complete overlap with the accuracy needed to trust the substantive estimates: CPI_no_imaging is a nested subset of CPI_full, and the pipeline returned",
    "The nested-subset consistency check recovered near-complete overlap with the accuracy needed to trust the substantive estimates: CPI_no_imaging is a subset of CPI_full, and the pipeline returned",
))
EDITS.append((
    "G1-08 site outcomes",
    "The internal positive control behaved as expected and, usefully, differentiated the two sites:",
    "The internal comparator again behaved as expected and, usefully, differentiated the two sites:",
))
EDITS.append((
    "G1-09 discussion bounded null",
    "Using BMI as an internal positive control shows that the instrument set, harmonization and outcome data were capable of detecting a body-composition effect in this specific contrast; a null obtained alongside such a control is a different quantity from a null obtained in isolation.",
    "Using BMI as an internal comparator shows that the instrument set, harmonization and outcome data were capable of detecting a body-composition effect in this specific contrast; a null obtained alongside such a comparator is a different quantity from a null obtained in isolation. The comparator is not a positive control with a known answer: what it establishes is the sensitivity of the same model in this specific contrast, not that any effect of the size we excluded would necessarily have been recovered.",
))
EDITS.append((
    "G1-10 discussion implication",
    "null results in multivariable MR should be reported with the minimum detectable effect and, where one exists, an internal positive control;",
    "null results in multivariable MR should be reported with the minimum detectable effect and, where one exists, an internal comparator;",
))
EDITS.append((
    "G1-11 table row label",
    "| CPI_full \u00d7 CPI_no_imaging (positive control) |",
    "| CPI_full \u00d7 CPI_no_imaging (nested-subset consistency check) |",
))
EDITS.append((
    "G1-12 table note",
    "**Table note:** audit positive controls on these data \u2014",
    "**Table note:** consistency checks on these data \u2014",
))

# ---------------- Group 2: same-cohort outcome wording ------------------------------------
EDITS.append((
    "G2-01 abstract results",
    "P=0.898), replicating for CPI_no_imaging;",
    "P=0.898), with the same pattern for CPI_no_imaging;",
))
EDITS.append((
    "G2-02 abstract site outcomes",
    "Two site-matched secondary outcomes reproduced the TMFI null",
    "Two site-matched secondary outcomes also returned a TMFI null",
))
EDITS.append((
    "G2-03 results primary",
    "The pattern was reproduced for `CPI_no_imaging`:",
    "The same pattern held for `CPI_no_imaging`, a subset of `CPI_full` and not an independent replication:",
))

# ---------------- Group 3: soften over-stated claims --------------------------------------
EDITS.append((
    "G3-01 heterogeneity direction (L215)",
    "Three observations delimit what this can and cannot explain: MVMR-Egger intercepts were indistinguishable from zero, so the scatter is not a signature of unbalanced directional pleiotropy; the two exposures move in opposite directions, which balanced noise across instruments cannot manufacture; and the fixed-effect fit reproduces both point estimates almost exactly, so heterogeneity inflates variance rather than shifting estimates. The Q result therefore limits a claim of complete instrument validity and motivates caution about absolute magnitudes, but provides no mechanism that would convert a genuinely null TMFI effect into the present data.",
    "Three observations delimit what this can and cannot explain: MVMR-Egger intercepts were indistinguishable from zero, so the scatter shows no signature of unbalanced directional pleiotropy, although a component below this diagnostic's resolution cannot be excluded; the two exposures move in opposite directions, which is not the pattern expected if the instruments acted through a single shared pathway, though this is a consistency argument rather than proof; and the fixed-effect fit reproduces both point estimates almost exactly, indicating that heterogeneity inflates variance rather than shifting estimates. The Q result therefore limits any claim of complete instrument validity and motivates caution about absolute magnitudes. It does not by itself exclude a true TMFI effect, and we have not identified any route among these diagnostics by which heterogeneity alone would be expected to convert such an effect into the null estimates reported here.",
))
EDITS.append((
    "G3-02 pathway / intervention target (L265)",
    "Substantively, for the specific contrast tested here, adiposity-related pathways carry the effect and thigh muscle composition does not. TMFI remains useful as a marker of muscle quality and as a candidate mediator or effect modifier, but these genetic data do not support advancing it as an independent interventional target for chronic pain intensity.",
    "Substantively, for the contrast tested here, genetically predicted systemic adiposity was associated with chronic pain intensity whereas genetically predicted thigh muscle composition was not; because the BMI instrument indexes adiposity as a whole, this localises the difference no further than systemic adiposity versus thigh muscle composition and does not identify a pathway. TMFI remains useful as a marker of muscle quality and as a candidate mediator or effect modifier, but these estimates describe lifelong genetically predicted differences and do not support advancing thigh muscle fat infiltration as an intervention target for chronic pain intensity while adiposity itself is unaddressed.",
))
EDITS.append((
    "G3-03 generic adiposity pathway (L257)",
    "once the adiposity pathway is accounted for.",
    "once systemic adiposity is accounted for.",
))

# ---------------- execute -----------------------------------------------------------------
fail = 0
for name, old, new in EDITS:
    n = txt.count(old)
    if n != 1:
        print("FAIL %-46s occurrences=%d (expected 1)" % (name, n))
        fail += 1
        continue
    txt = txt.replace(old, new)
    print("ok   %-46s +%d chars" % (name, len(new) - len(old)))

if fail:
    print("\nABORT: %d edit(s) unmatched, nothing written." % fail)
    sys.exit(1)

out = txt.replace("\n", "\r\n") if crlf else txt
bak = P + ".bak_2026-09-20b"
if not os.path.isfile(bak):
    io.open(bak, "wb").write(raw.encode("utf-8"))
io.open(P, "wb").write(out.encode("utf-8"))
print("\nwritten: %s (CRLF=%s), backup: %s" % (P, crlf, os.path.basename(bak)))
