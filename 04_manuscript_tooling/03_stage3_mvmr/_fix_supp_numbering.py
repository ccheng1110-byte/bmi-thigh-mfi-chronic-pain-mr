# -*- coding: utf-8 -*-
"""Renumber the '## Supplementary material' list to 1..N after inserting S17."""
import io
import re

P = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"
with io.open(P, "r", encoding="utf-8", newline="") as f:
    text = f.read()

start = text.index("## Supplementary material")
end = text.index("## Supplementary Methods")
seg = text[start:end]
lines = seg.split("\r\n")

pat = re.compile(r"^(\d+)\.\s")
idx = [i for i, l in enumerate(lines) if pat.match(l)]
print("numbered items found:", len(idx))
for i in idx:
    print("   old %-3s %s" % (pat.match(lines[i]).group(1), lines[i][:70]))

for n, i in enumerate(idx, 1):
    lines[i] = pat.sub("%d. " % n, lines[i], count=1)

print("\nafter renumbering:")
for i in idx:
    print("   new %-3s %s" % (pat.match(lines[i]).group(1), lines[i][:70]))

text = text[:start] + "\r\n".join(lines) + text[end:]
with io.open(P, "w", encoding="utf-8", newline="") as f:
    f.write(text)

with io.open(P, "rb") as f:
    d = f.read()
print("\nCRLF=%d bareLF=%d" % (d.count(b"\r\n"), len(re.findall(rb"(?<!\r)\n", d))))
