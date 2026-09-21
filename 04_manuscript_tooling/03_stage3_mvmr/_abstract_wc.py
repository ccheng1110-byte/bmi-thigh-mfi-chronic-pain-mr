# -*- coding: utf-8 -*-
"""Per-section word count of the abstract, prose-only convention."""
import io, os, re, sys

D = os.path.dirname(os.path.abspath(__file__))
MAN = os.path.join(D, "Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md")
with io.open(MAN, "r", encoding="utf-8") as f:
    t = f.read()

ab = t.split("## Abstract", 1)[1].split("\n## Introduction", 1)[0]
W = lambda s: len(re.findall(r"[^\s]+", s))

mk = re.search(r"\*\*Keywords:\*\*(.*)$", ab, re.S)
kw = W(mk.group(1)) if mk else 0

secs = re.split(r"\n### ", ab)
total = 0
for s in secs:
    if not s.strip():
        continue
    name = s.split("\n", 1)[0].strip()
    if name == "Background" or name.startswith("Background"):
        body = s
    body = s
    # drop the "### Name" heading line
    body = "\n".join(body.split("\n")[1:]) if "\n" in body else ""
    n = W(body)
    total += n
    print("%-14s %4d" % (name, n))
print("%-14s %4d   <-- prose (headings + 'Keywords:' label excluded)" % ("TOTAL", total))
print("%-14s %4d" % ("keywords", kw))
print("cap 350 -> headroom:", 350 - total)
