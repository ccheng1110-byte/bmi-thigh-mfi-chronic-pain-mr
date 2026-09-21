# -*- coding: utf-8 -*-
"""Residual 'positive control' occurrences missed in the first pass (L269 limitations, L449 methods summary)."""
import io
import os
import sys

P = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"
raw = io.open(P, "rb").read().decode("utf-8")
crlf = raw.count("\r\n") > 0
txt = raw.replace("\r\n", "\n")

EDITS = [
    ("F-01 limitations nested-subset check",
     "because the nested-subset positive control recovered near-complete overlap to within 2%",
     "because the nested-subset consistency check recovered near-complete overlap to within 2%"),
    ("F-02 methods summary BMI role",
     "with BMI serving as an internal positive control in the same model",
     "with BMI serving as an internal comparator in the same model"),
]

fail = 0
for name, old, new in EDITS:
    n = txt.count(old)
    if n != 1:
        print("FAIL %-40s occurrences=%d" % (name, n))
        fail += 1
        continue
    txt = txt.replace(old, new)
    print("ok   %-40s +%d chars" % (name, len(new) - len(old)))

if fail:
    print("ABORT")
    sys.exit(1)

out = txt.replace("\n", "\r\n") if crlf else txt
io.open(P, "wb").write(out.encode("utf-8"))
print("written", P, "CRLF=", crlf)
