# -*- coding: utf-8 -*-
"""P2 step 1: zero-download diagnostic of LD reference panel coverage.

Answers:
  1. exact panel size / chromosome spread
  2. candidate -> panel coverage (reproduce 22.75%)
  3. profile of the missing candidates: MAF, P, chromosome
  4. how many missing candidates are "locally covered" i.e. have a panel
     variant within X kb (clumping window is 10,000 kb but LD decays fast)
"""
import io, os, re, csv, sys
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(ROOT)
BIM = os.path.join(REPO, "reference", "1000G_EUR.bim")
CAND = os.path.join(ROOT, "tmfi_gws_exact_clean.csv")
CLUMPED = os.path.join(ROOT, "tmfi_gws_ldclumped_r2_0p001_kb10000.csv")

print("=" * 68)
print("P2 STEP 1 : panel coverage diagnostic (no downloads)")
print("=" * 68)

# ---------- panel ----------
panel = {}          # rsid -> (chr, pos)
panel_pos = defaultdict(list)   # chr -> sorted positions
n_bim = 0
with io.open(BIM, "r", encoding="utf-8") as f:
    for line in f:
        p = line.split()
        n_bim += 1
        panel[p[1]] = (p[0], int(p[3]))
        panel_pos[p[0]].append(int(p[3]))
for c in panel_pos:
    panel_pos[c].sort()
print("panel variants            : %d" % n_bim)
print("panel chromosomes         : %s" % ",".join(sorted(panel_pos, key=lambda x: int(x))))
print("panel median gap (chr1)   : %.0f bp" % (
    sorted(panel_pos["1"][i+1] - panel_pos["1"][i] for i in range(len(panel_pos["1"])-1))[len(panel_pos["1"])//2]))
print("panel variants per chrom  : min=%d max=%d" % (
    min(len(v) for v in panel_pos.values()), max(len(v) for v in panel_pos.values())))

# ---------- candidates ----------
with io.open(CAND, "r", encoding="utf-8") as f:
    rd = list(csv.DictReader(f))
print("\ncandidate file columns    : %s" % rd[0].keys())
print("candidate rows            : %d" % len(rd))

# normalise column names
def pick(keys, *alts):
    for a in alts:
        for k in keys:
            if k.lower() == a.lower():
                return k
    return None

krs = pick(rd[0].keys(), "SNP", "snp", "rsid", "Variant", "variant", "ID")
kchr = pick(rd[0].keys(), "chr", "chrom", "chromosome", "CHR")
kpos = pick(rd[0].keys(), "pos", "position", "BP", "bp")
kp = pick(rd[0].keys(), "p", "pval", "P", "p_value", "pvalue")
kmaf = pick(rd[0].keys(), "maf", "MAF", "freq", "eaf")
print("resolved cols             : rsid=%s chr=%s pos=%s p=%s maf=%s" % (krs, kchr, kpos, kp, kmaf))

cands = []
for r in rd:
    rs = r[krs].strip()
    try:
        c = str(int(float(r[kchr])))
    except Exception:
        c = r[kchr].strip()
    try:
        pos = int(float(r[kpos]))
    except Exception:
        pos = None
    try:
        pv = float(r[kp])
    except Exception:
        pv = None
    try:
        maf = float(r[kmaf])
    except Exception:
        maf = None
    cands.append((rs, c, pos, pv, maf))

inpanel = [x for x in cands if x[0] in panel]
missing = [x for x in cands if x[0] not in panel]
print("\nclean strong candidates   : %d" % len(cands))
print("  in panel                : %d (%.2f%%)" % (len(inpanel), 100.0*len(inpanel)/len(cands)))
print("  missing from panel      : %d (%.2f%%)" % (len(missing), 100.0*len(missing)/len(cands)))

# positional double-check (rsid alias check) for a sample
import random
random.seed(20260920)
sample = random.sample(missing, min(400, len(missing)))
pos_hit = 0
for rs, c, pos, pv, maf in sample:
    if c in panel_pos and pos is not None:
        lst = panel_pos[c]
        lo, hi = 0, len(lst)
        while lo < hi:
            mid = (lo+hi)//2
            if lst[mid] < pos: lo = mid+1
            else: hi = mid
        if lo < len(lst) and lst[lo] == pos:
            pos_hit += 1
        elif lo > 0 and abs(lst[lo-1]-pos) <= 2:
            pos_hit += 1
print("  positional re-check (n=%d sample): %d variants DO exist at same bp -> rsid alias issue"
      % (len(sample), pos_hit))

# ---------- profile of missing ----------
def med(v):
    v = sorted([x for x in v if x is not None])
    return v[len(v)//2] if v else None

print("\n--- profile ---")
print("MAF   median  in-panel=%.3f   missing=%.3f" % (med([x[4] for x in inpanel]), med([x[4] for x in missing])))
print("-log10(P) median  in-panel=%.2f   missing=%.2f" % (
    med([-1*__import__("math").log10(x[3]) for x in inpanel if x[3]]),
    med([-1*__import__("math").log10(x[3]) for x in missing if x[3]])))
cm_in = Counter(x[1] for x in inpanel)
cm_mi = Counter(x[1] for x in missing)
print("\nchrom   in-panel   missing   ratio")
for c in sorted(set(list(cm_in)+list(cm_mi)), key=lambda x: int(x)):
    a, b = cm_in.get(c, 0), cm_mi.get(c, 0)
    print("  %-4s  %5d     %5d     %.2f" % (c, a, b, a/float(a+b)))

# MAF strata coverage
print("\nMAF strata coverage")
strata = [(0.0, 0.01), (0.01, 0.05), (0.05, 0.10), (0.10, 0.20), (0.20, 0.45), (0.45, 1.01)]
for lo, hi in strata:
    sub = [x for x in cands if x[4] is not None and lo <= x[4] < hi]
    if not sub:
        continue
    inn = sum(1 for x in sub if x[0] in panel)
    print("  MAF [%.2f,%.2f)  n=%5d  in-panel=%5d  (%.1f%%)" % (lo, hi, len(sub), inn, 100.0*inn/len(sub)))

# ---------- local coverage of missing variants ----------
import bisect
print("\n--- local coverage: is a panel variant near each missing candidate? ---")
for win in (10, 50, 100, 250, 500, 1000):
    cov = 0
    for rs, c, pos, pv, maf in missing:
        if c not in panel_pos or pos is None:
            continue
        lst = panel_pos[c]
        i = bisect.bisect_left(lst, pos)
        ok = (i < len(lst) and lst[i]-pos <= win*1000) or (i > 0 and pos-lst[i-1] <= win*1000)
        if ok:
            cov += 1
    print("  within %5d kb : %5d / %d  (%.1f%%)" % (win, cov, len(missing), 100.0*cov/len(missing)))

# ---------- clumped set ----------
with io.open(CLUMPED, "r", encoding="utf-8") as f:
    cl = list(csv.DictReader(f))
print("\nclumped instruments       : %d" % len(cl))

# ---------- export missing list for remote probe ----------
out = os.path.join(ROOT, "p2_missing_candidates.txt")
with io.open(out, "w", encoding="utf-8") as f:
    f.write("rsid\tchr\tpos\tmaf\tp\n")
    for rs, c, pos, pv, maf in missing:
        f.write("%s\t%s\t%s\t%s\t%s\n" % (rs, c, pos if pos is not None else "",
                                          "" if maf is None else maf, "" if pv is None else pv))
print("missing list written      : %s" % out)
print("\nDONE step1")
