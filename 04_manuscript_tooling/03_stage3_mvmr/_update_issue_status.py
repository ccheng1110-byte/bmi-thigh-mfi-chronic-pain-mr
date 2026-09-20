# -*- coding: utf-8 -*-
"""Tick the checklist items completed in the 2026-09-19 finishing pass."""
from pathlib import Path

P = Path(r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Manuscript_Substantive_Issues_and_Corrections_2026-09-19.md")

# 1-based line numbers of items completed in this pass
DONE = [228, 229, 230, 231, 235, 236, 237, 238, 239, 240] + list(range(263, 273)) + list(range(278, 285))

# line number -> note appended after the item text
NOTES = {
    228: " Main text and Tables 5/8 unified on the R `MVMR` values 22.12 / 35.23.",
    229: " The closed-form approximations 23.6 / 42.8 were moved to Supplementary Table S11; the main text mentions them only as a sensitivity analysis.",
    230: " Stated explicitly that the two sets are not numerically comparable and do not validate each other.",
    235: " Figure numbers renumbered to Figure 1-4, with no duplicates.",
    236: " All cross-references to Figures 1-4 have been checked throughout.",
    237: " Changed to 'Seven boundaries', mapped item by item.",
    238: " Stated in Methods: the instrument set, outcome definitions, estimator hierarchy and sensitivity analyses were fixed before the Stage 2 results were inspected.",
    239: " Completed the 20-item checklist (STROBE-MR_Checklist_TMFI_BMI_CPI_2026-09-19.md) and declared it in the main text as Supplementary Table S12.",
    240: " Changed to 'will be deposited at submission'; no longer claimed as already deposited.",
    263: " **[PARTIAL]** Disclosed and switched to the standardised scale; the physical units and the `fastGWA-mlm-binary` flag remain to be confirmed with the source authors.",
    264: " Confirmed (BOLT-LMM risk difference).",
    265: " Recomputed and back-filled using lambda = beta / [mu(1-mu)].",
    266: " Unified on the normal-theory rule.",
    267: " Uniformly back-filled by script, with a residual scan.",
    268: " Revised 27 overreaching statements.",
    269: " Reframed as a within-cohort consistency check.",
    270: " Main text compressed from 9,005 to 8,650 words; MD5 checksums and the inference history moved to the supplement.",
    271: " Figure numbering, declarations and STROBE-MR completed.",
    272: " Completed (Figures 1-4 and Tables 1-8 references consistent; abstract 339 words).",
    278: " **[PARTIAL]** Scale still undetermined; all scale-dependent claims rewritten on the standardised scale.",
    279: " Confirmed (BOLT-LMM risk difference, lambda transformation).",
    280: " Unified (+/-1.96 x SE interval + normal-theory P).",
    281: " Revised.",
    282: " Reframed as scenario-qualified plus sensitivity analysis.",
    283: " Revised.",
    284: " Compressed the technical-defence passage.",
}


def main():
    lines = P.read_text(encoding="utf-8").split("\n")
    n = 0
    for ln in DONE:
        i = ln - 1
        if i >= len(lines):
            print("SKIP out of range", ln)
            continue
        if "[ ]" not in lines[i]:
            print("SKIP no checkbox at line", ln, "::", lines[i][:70])
            continue
        note = NOTES.get(ln, "")
        lines[i] = lines[i].replace("[ ]", "[x]", 1) + note
        n += 1
    P.write_text("\n".join(lines), encoding="utf-8")
    print("ticked %d items" % n)
    print("remaining unchecked in file:", sum(1 for l in lines if "[ ]" in l))


if __name__ == "__main__":
    main()
