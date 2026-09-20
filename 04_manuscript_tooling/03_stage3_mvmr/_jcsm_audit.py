# -*- coding: utf-8 -*-
"""Audit the manuscript against JCSM hard limits."""
import re, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

P = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"
raw = open(P, 'rb').read()
text = raw.decode('utf-8')
lines = text.split('\n')
print("lines:", len(lines))

W = lambda s: len(re.findall(r"[^\s]+", s))

# ---- title ----
title = lines[0].lstrip('#').strip()
print("\n[TITLE] words =", W(title))
print(title)

# ---- abstract ----
ab = text.split('## Abstract', 1)[1].split('\n## Introduction', 1)[0]
ab_body = '\n'.join(l for l in ab.split('\n') if not l.startswith('#') and '**Keywords:**' not in l)
print("\n[ABSTRACT] body words =", W(ab_body))

# ---- keywords ----
m = re.search(r"\*\*Keywords:\*\*\s*(.+)", text)
kws = [k.strip() for k in m.group(1).split(';') if k.strip()]
print("[KEYWORDS] n =", len(kws), "->", kws)

# ---- body sections ----
body_start = text.index('## Introduction')
body_end = text.index('## References') if '## References' in text else len(text)
body = text[body_start:body_end]

secs = []
cur = None
buf = []
for ln in body.split('\n'):
    if ln.startswith('## '):
        if cur: secs.append((cur, '\n'.join(buf)))
        cur = ln[3:].strip(); buf = []
    else:
        buf.append(ln)
if cur: secs.append((cur, '\n'.join(buf)))

print("\n[SECTIONS] (title lines excluded from counts)")
tot = 0
for name, cont in secs:
    c = '\n'.join(l for l in cont.split('\n') if not l.startswith('#'))
    w = W(c); tot += w
    print("  %-28s %6d" % (name[:28], w))
print("  %-28s %6d" % ("TOTAL BODY", tot))

# ---- references ----
refs = text[text.index('## References'):]
ref_lines = [l for l in refs.split('\n') if re.match(r"^\d+\.\s", l.strip())]
print("\n[REFERENCES] n =", len(ref_lines))

# ---- tables & figures ----
tables = sorted(set(int(x) for x in re.findall(r"^### Table (\d+)", text, re.M)))
figs = sorted(set(int(x) for x in re.findall(r"^### Figure (\d+)", text, re.M)))
if not figs:
    figs = sorted(set(int(x) for x in re.findall(r"\*\*Figure (\d+)", text)))
print("\n[TABLES in main text] n =", len(tables), tables)
print("[FIGURES in main text] n =", len(figs), figs)
print("[MAIN FLOAT TOTAL] =", len(tables) + len(figs), " (JCSM limit: 8)")

# ---- JCSM-countable body: exclude tables, figure legends, abstract, refs ----
countable = 0
for name, cont in secs:
    if name in ('Main tables', 'Figure legends'):
        continue
    c = '\n'.join(l for l in cont.split('\n') if not l.startswith('#'))
    countable += W(c)
print("\n[JCSM-COUNTABLE BODY] (Introduction..Conclusion, excl. tables & fig legends) =", countable)

# ---- JCSM verdict ----
print("\n===== JCSM HARD LIMITS =====")
print("Title      : %d / 17   %s" % (W(title), "OK" if W(title) <= 17 else "OVER by %d" % (W(title)-17)))
print("Abstract   : %d / 400  %s" % (W(ab_body), "OK" if W(ab_body) <= 400 else "OVER"))
print("Body words : %d / 4500 %s" % (countable, "OK" if countable <= 4500 else "OVER by %d (must cut %.0f%%)" % (countable-4500, 100*(countable-4500)/countable)))
print("Keywords   : %d / 4-6  %s" % (len(kws), "OK" if 4 <= len(kws) <= 6 else "OUT OF RANGE"))
print("References : %d / 40   %s" % (len(ref_lines), "OK" if len(ref_lines) <= 40 else "OVER"))
print("Main floats: %d / 8    %s" % (len(tables)+len(figs), "OK" if len(tables)+len(figs) <= 8 else "OVER by %d" % (len(tables)+len(figs)-8)))
