# -*- coding: utf-8 -*-
import io, re, os
D = os.path.dirname(os.path.abspath(__file__))
files = {
    "PRE  (.bak_origpack2026-09-20)": os.path.join(D, "Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md.bak_origpack2026-09-20"),
    "POST (current)": os.path.join(D, "Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"),
}
for label, p in files.items():
    with io.open(p, "r", encoding="utf-8") as f:
        t = f.read()
    body = t.split("\n## Introduction", 1)[1].split("\n## Main tables", 1)[0]
    meth = t.split("\n## Methods", 1)[1].split("\n## Results", 1)[0]
    res = t.split("\n## Results", 1)[1].split("\n## Discussion", 1)[0]
    disc = t.split("\n## Discussion", 1)[1].split("\n## Conclusion", 1)[0]
    intro = t.split("\n## Introduction", 1)[1].split("\n## Methods", 1)[0]
    print("%-32s body=%5d (intro %4d, methods %5d, results %5d, discussion %5d) lines=%d" % (
        label,
        len(re.findall(r"[^\s]+", body)),
        len(re.findall(r"[^\s]+", intro)),
        len(re.findall(r"[^\s]+", meth)),
        len(re.findall(r"[^\s]+", res)),
        len(re.findall(r"[^\s]+", disc)),
        t.count("\n") + 1))
