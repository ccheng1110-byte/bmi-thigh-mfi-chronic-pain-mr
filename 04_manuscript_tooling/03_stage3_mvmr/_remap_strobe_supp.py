# -*- coding: utf-8 -*-
"""Shift STROBE-MR line references >= 244 by +1 (S17 inserted at L243)."""
import io
import re

P = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\STROBE-MR_Checklist_TMFI_BMI_CPI_2026-09-20.md"
BAK = P + ".bak_remap2026-09-20"

with io.open(P, "r", encoding="utf-8", newline="") as f:
    text = f.read()
with io.open(BAK, "w", encoding="utf-8", newline="") as f:
    f.write(text)

THRESH = 244


def bump(m):
    n = int(m.group(1))
    return "L%d" % (n + 1) if n >= THRESH else m.group(0)


new, k = re.subn(r"L(\d+)", bump, text)
print("L-refs shifted:", sum(1 for m in re.finditer(r"L(\d+)", text) if int(m.group(1)) >= THRESH))

with io.open(P, "w", encoding="utf-8", newline="") as f:
    f.write(new)

for m in re.finditer(r"L(2[4-9]\d|3\d\d)", new):
    s = max(0, m.start() - 60)
    print("   ", new[s:m.end() + 20].replace("\n", " "))
