# -*- coding: utf-8 -*-
"""Final audit of the JCSM-adapted manuscript after the revision pass."""
import re, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

P = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"
raw = open(P, 'rb').read()
text = raw.decode('utf-8').replace('\r\n', '\n')
W = lambda s: len(re.findall(r"[^\s]+", s))
lines = text.split('\n')

print("lines:", len(lines))
print("CRLF:", raw.count(b'\r\n'), " bare LF:", len(re.findall(rb'(?<!\r)\n', raw)))

title = lines[0].lstrip('#').strip()
ab = text.split('## Abstract', 1)[1].split('\n## Introduction', 1)[0]
ab_body = '\n'.join(l for l in ab.split('\n') if not l.startswith('#') and '**Keywords:**' not in l)
kws = [k.strip() for k in re.search(r"\*\*Keywords:\*\*\s*(.+)", text).group(1).split(';') if k.strip()]
tracked = W(text[text.index('## Introduction'):text.index('## Main tables')])
seg = text[text.index('## Introduction'):text.index('## Main tables')]
strict = W('\n'.join(l for l in seg.split('\n') if not l.startswith('#')))
refs = len([l for l in text[text.index('## References'):].split('\n') if re.match(r"^\d+\.\s", l.strip())])
tabs = len(re.findall(r"^### Table \d", text, re.M))
figs = len(re.findall(r"^### Figure \d", text, re.M))

print("\n===== JCSM HARD LIMITS =====")
for label, val, lim, ok in [
    ("Title (words)", W(title), 17, W(title) <= 17),
    ("Abstract (words)", W(ab_body), 400, W(ab_body) <= 400),
    ("Main text TRACKED (incl. ### heads)", tracked, 4500, tracked <= 4500),
    ("Main text STRICT (prose only)", strict, 4500, strict <= 4500),
    ("Keywords", len(kws), (4, 6), 4 <= len(kws) <= 6),
    ("References", refs, 40, refs <= 40),
    ("Main tables+figures", tabs + figs, 8, tabs + figs <= 8),
]:
    print("  %-36s %-8s limit %-8s %s" % (label, val, lim, "OK" if ok else "*** FAIL ***"))

print("\n  margin, tracked metric:", 4500 - tracked, "words")
print("  margin, strict metric :", 4500 - strict, "words")

# ---------------------------------------------------------------- content integrity
print("\n===== CONTENT INTEGRITY =====")
MUST = ["0.008", "0.558", "0.005", "0.721", "0.090", "0.038–0.143", "0.00075",
        "−0.002", "0.898", "22.12", "35.23", "0.048", "0.052", "1.2×10⁻¹¹",
        "1.0146", "1.0005", "1.2751", "1.1344", "0.983", "1.217", "1.116–1.327",
        "9.0×10⁻⁶", "46,347", "134,627", "106,697", "499,362", "22.75%", "1.61",
        "0.905", "1.26", "3,623", "42/42", "0.031", "0.037", "0.910", "0.890",
        "0.0014", "1.0681", "1.0621", "1.019", "1.001", "0.176", "0.196",
        "descriptively larger", "internal comparator", "bounded null",
        "structural ceiling of 0.890 within ~2 SE", "Supplementary Table S13",
        "Eight boundaries", "in three groups", "Outcome-side triangulation",
        "directionally concordant but imprecise", "conditional on the exposure scale",
        "fixed in advance", "secondary and exploratory", "no overlap bias was expected",
        "Mendelian randomization estimated", "MRlap", "STROBE-MR"]
bad = [m for m in MUST if m not in text]
print("  must-keep: %d checked, %d missing %s" % (len(MUST), len(bad), bad if bad else ""))

GONE = ["Two features make this more than a failed test", "licenses the interpretation",
        "balanced scatter rather than directional pleiotropy",
        "within expectation for three endpoints under the no-correction rule",
        "reproduced both central findings", "adiposity gradient across weight-bearing sites",
        "cross-trait intercept is zero by construction", "inference does not depend on re-clumping",
        "pre-specified as an internal comparator", "pre-specified replications of it",
        "Both contrasts were recovered in FinnGen", "Table 8", "Table 7", "Table 9",
        "Figure 3", "Figure 4", "prespecified"]
left = [g for g in GONE if g in text]
print("  must-go  : %d checked, %d remaining %s" % (len(GONE), len(left), left if left else ""))

# ---------------------------------------------------------------- remaining author tasks
print("\n===== REMAINING AUTHOR TASKS =====")
for tag, pat in [("author placeholders", r"To be completed by the authors"),
                 ("code DOI placeholders", r"will be deposited")]:
    print("  %-26s %d" % (tag, len(re.findall(pat, text))))

print("\n===== STRUCTURE =====")
for ln in lines:
    if ln.startswith('## '):
        print("  ", ln)
