# -*- coding: utf-8 -*-
"""Fix the self-contradictory 'intercept above the ceiling' phrasing (L122, L219).

Observed 0.910 (SE 0.010) vs ceiling 0.890 -> the intercept sits AT the structural ceiling
within ~2 SE, not 'against' it. Terminology is also aligned with the rest of the paper,
which uses 'structural ceiling' rather than 'theoretical ceiling'.
"""
import io
import os
import sys

P = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"
raw = io.open(P, "rb").read().decode("utf-8")
crlf = raw.count("\r\n") > 0
txt = raw.replace("\r\n", "\n")

EDITS = [
    ("C-01 methods pipeline check (L122)",
     "returned a cross-trait intercept of 0.910 (SE 0.010) against a theoretical ceiling of "
     "\u221a(105,977/133,721) = 0.890, a genetic correlation of 1.0025, and N_shared\u00b7\u03c1_phen = 108,324 "
     "against a nominal subset size of 106,697, i.e. agreement within 2%.",
     "returned a cross-trait intercept of 0.910 (SE 0.010), at its structural ceiling of "
     "\u221a(105,977/133,721) = 0.890 within ~2 SE; a genetic correlation of 1.0025; and N_shared\u00b7\u03c1_phen = 108,324 "
     "against a nominal subset size of 106,697, i.e. agreement within 2%."),
    ("C-02 results overlap (L219)",
     "returned a cross-trait intercept of 0.910 (SE 0.010) against a theoretical ceiling of 0.890, "
     "and a genetic correlation of 1.0025.",
     "returned a cross-trait intercept of 0.910 (SE 0.010), at its structural ceiling of 0.890 "
     "within ~2 SE, and a genetic correlation of 1.0025."),
]

fail = 0
for name, old, new in EDITS:
    n = txt.count(old)
    if n != 1:
        print("FAIL %-34s occurrences=%d" % (name, n))
        fail += 1
        continue
    txt = txt.replace(old, new)
    print("ok   %-34s +%d chars" % (name, len(new) - len(old)))

if fail:
    print("ABORT, nothing written.")
    sys.exit(1)

out = txt.replace("\n", "\r\n") if crlf else txt
io.open(P, "wb").write(out.encode("utf-8"))
print("written", P, "CRLF=", crlf)
