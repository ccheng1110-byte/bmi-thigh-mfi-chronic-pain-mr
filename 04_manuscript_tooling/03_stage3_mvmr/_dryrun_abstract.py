# -*- coding: utf-8 -*-
"""Dry-run the abstract rewrite and check the 350-word ceiling before writing anything."""
import io
import re

P = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"
lines = io.open(P, encoding="utf-8").read().replace("\r\n", "\n").split("\n")

BG_OLD = lines[10]
BG_NEW = ("Thigh muscle fat infiltration (TMFI) and body mass index (BMI) may both influence chronic pain "
          "intensity (CPI) but are hard to separate. BMI was pre-specified as an internal comparator, so a "
          "null TMFI estimate could be read against an explicit power bound.")

ME_NEW = ("We used published summary statistics for TMFI (n=46,347), BMI (GIANT 2015) and two UK Biobank "
          "chronic pain intensity outcomes (n=134,627 and 106,697). We first estimated the total effect of "
          "TMFI using 42 independent genetic variants, then estimated TMFI and BMI jointly using 93 variants "
          "so that each effect is adjusted for the other, with standard diagnostics for instrument strength "
          "and pleiotropy. Because exposure and outcomes come from the same cohort, that overlap was "
          "quantified genetically and carried through every estimate; results are reported against the "
          "smallest effect detectable at 80% power rather than significance alone.")

RE_NEW = ("No TMFI total effect emerged (0.008 SD per SD, 95% CI \u22120.018 to 0.033, P=0.558; 0.005, "
          "\u22120.022 to 0.032, P=0.721). The joint model gave a positive BMI direct effect (0.090 SD per "
          "SD, 95% CI 0.038\u20130.143, P=0.00075) against a null TMFI effect (\u22120.002, \u22120.036 to "
          "0.031, P=0.898), with the same pattern for the narrower pain definition; instruments remained "
          "strong after conditioning (F 22.1 and 35.2), with no evidence of directional pleiotropy. The "
          "study could detect a TMFI effect of 0.048 SD at 80% power, 1.9-fold smaller than the BMI effect "
          "recovered. Exposure and outcomes share participants (intercept 0.052, P=1.2\u00d710\u207b\u00b9\u00b9); "
          "this moved estimates by \u22640.0015 SD per SD, and no degree of sharing could reverse either "
          "finding. Two site-specific secondary outcomes also returned a TMFI null (knee pain OR 1.015 per "
          "SD, P=0.58; back pain 1.001, P=0.98), with a descriptively larger BMI effect at the "
          "weight-bearing joint (1.275 vs 1.134).")

CO_NEW = ("Genetically predicted adiposity, but not thigh muscle fat infiltration, is associated with chronic "
          "pain intensity. Against an explicit power bound, an internal comparator that behaved as expected "
          "and quantified participant overlap, this is a bounded null: effects above about 0.05 SD per SD "
          "are unlikely; smaller ones remain open.")


def nw(s):
    s = re.sub(r"^#+\s*", "", s)
    return len([w for w in s.split() if any(c.isalnum() for c in w)])


blocks = [("Background", lines[10], BG_NEW),
          ("Methods", lines[14], ME_NEW),
          ("Results", lines[18], RE_NEW),
          ("Conclusions", lines[22], CO_NEW)]

old_total = new_total = 0
for name, old, new in blocks:
    old_total += nw(old)
    new_total += nw(new)
    print("%-12s %3d -> %3d  (%+d)" % (name, nw(old), nw(new), nw(new) - nw(old)))

print("\nabstract body  %d -> %d  (%+d)" % (old_total, new_total, new_total - old_total))
kw = nw(lines[24])
print("keywords line   %d" % kw)
print("TOTAL incl. keywords: %d  (ceiling 350) -> %s" % (new_total + kw, "OK" if new_total + kw <= 350 else "OVER"))
