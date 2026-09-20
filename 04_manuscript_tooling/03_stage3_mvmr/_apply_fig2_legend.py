# -*- coding: utf-8 -*-
"""Update the Figure 2 legend: the figure is now on a common standardized scale."""
import io
import os
import sys

P = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"
raw = io.open(P, "rb").read().decode("utf-8")
crlf = raw.count("\r\n") > 0
txt = raw.replace("\r\n", "\n")

EDITS = [
    ("F2 legend",
     "Blue marks TMFI and orange marks BMI. Because the two exposures are plotted on their own source "
     "scales, cross-exposure comparison should use the standardized values in Table 3 rather than the "
     "plotted heights. The vertical dashed line denotes zero.",
     "Blue marks TMFI and orange marks BMI. Both exposures are plotted on the same standardized scale "
     "\u2014 SD of pain intensity per 1 SD of the exposure \u2014 so plotted heights are comparable "
     "across exposures and correspond to the standardized column of Table 3, where the raw-unit "
     "estimates are also given. The vertical dashed line denotes zero."),
]

fail = 0
for name, old, new in EDITS:
    n = txt.count(old)
    if n != 1:
        print("FAIL %s occurrences=%d" % (name, n))
        fail += 1
        continue
    txt = txt.replace(old, new)
    print("ok   %-12s %+d chars" % (name, len(new) - len(old)))

if fail:
    sys.exit(1)

out = txt.replace("\n", "\r\n") if crlf else txt
io.open(P, "wb").write(out.encode("utf-8"))
print("written", P)
