# -*- coding: utf-8 -*-
"""P2 step 7: quantitative impact -- what would the extra instruments buy?

Also resolves the +7 vs +38 discrepancy by checking whether the strongest
"untagged" candidates sit INSIDE a locus that already has an instrument.
"""
import io, os, csv, math
from collections import defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(ROOT)
N_TMFI = 46347
SD_TMFI = 2.030884204160846

print("=" * 68)
print("P2 STEP 7 : quantitative impact of a denser panel")
print("=" * 68)

rows = []
with io.open(os.path.join(ROOT, "tmfi_gws_exact_clean.csv"), "r", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        rows.append({"rs": r["SNP"].strip(), "chr": r["CHR"].strip(),
                     "pos": int(float(r["POS"])), "p": float(r["P"]),
                     "maf": float(r["MAF"]), "beta": float(r["BETA"]),
                     "se": float(r["SE"]), "f": float(r["F"])})
instr = []
with io.open(os.path.join(ROOT, "tmfi_gws_ldclumped_r2_0p001_kb10000.csv"), "r", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        instr.append(r["SNP"].strip())
iset = set(instr)
byrs = {r["rs"]: r for r in rows}

# ---------- loci ----------
def loci_of(gap):
    L = []
    for c in sorted(set(r["chr"] for r in rows), key=lambda x: int(x)):
        sub = sorted([r for r in rows if r["chr"] == c], key=lambda r: r["pos"])
        cur = [sub[0]]
        for r in sub[1:]:
            if r["pos"] - cur[-1]["pos"] > gap:
                L.append(cur); cur = [r]
            else:
                cur.append(r)
        L.append(cur)
    return L

for gap in (100_000, 500_000, 1_000_000):
    L = loci_of(gap)
    ni = [sum(1 for x in l if x["rs"] in iset) for l in L]
    print("gap %7d bp : %2d loci, instruments per locus: %s" % (
        gap, len(L), dict(sorted(defaultdict(int, {k: ni.count(k) for k in set(ni)}).items()))))

# ---------- where do the strongest untagged-but-strong candidates sit? ----------
print("\n--- are the strongest 'untagged' candidates in loci that already have an instrument? ---")
panel = set()
with io.open(os.path.join(REPO, "reference", "1000G_EUR.bim"), "r", encoding="utf-8") as f:
    for line in f:
        panel.add(line.split()[1])
L = loci_of(500_000)
loc_of = {}
for li, l in enumerate(L):
    for x in l:
        loc_of[x["rs"]] = li
ins_loc = {}
for s in instr:
    ins_loc.setdefault(loc_of[s], []).append(s)

strong = sorted([r for r in rows if r["rs"] not in iset and r["p"] < 1e-20],
                key=lambda r: r["p"])
print("candidates with P<1e-20 that are NOT instruments: %d" % len(strong))
inloc = 0
for r in strong[:20]:
    li = loc_of[r["rs"]]
    has = li in ins_loc
    if has: inloc += 1
    print("  %-14s chr%-3s P=%.2e  in panel=%-5s  locus#%d %s" % (
        r["rs"], r["chr"], r["p"], r["rs"] in panel, li,
        ("HAS instrument " + ins_loc[li][0]) if has else "NO instrument"))
print("  ... of these, %d / %d sit in a locus that already has an instrument" % (
    sum(1 for r in strong if loc_of[r["rs"]] in ins_loc), len(strong)))

# ---------- R2 accounting ----------
def r2_of(r):
    return r["f"] / (N_TMFI - 2 + r["f"])

sum42 = sum(r2_of(byrs[s]) for s in instr)
print("\n--- variance explained ---")
print("  42 retained instruments : sum R2 = %.5f  (mean %.5f)" % (sum42, sum42/42))

blind = [l for l in L if not any(x["rs"] in panel for x in l)]
extra = [min(l, key=lambda x: x["p"]) for l in blind]
sum7 = sum(r2_of(x) for x in extra)
print("  7 blind-locus best SNPs : sum R2 = %.5f  (mean %.5f)" % (sum7, sum7/len(extra)))
print("  combined (49)           : sum R2 = %.5f  (+%.1f%%)" % (
    sum42 + sum7, 100.0*sum7/sum42))

# MDE scales as 1/sqrt(sum R2)
mde_now = 0.048
mde_new = mde_now * math.sqrt(sum42 / (sum42 + sum7))
print("\n--- implied change in the minimum detectable effect ---")
print("  MDE with 42 instruments : %.3f SD per SD (as reported)" % mde_now)
print("  MDE with 49 instruments : %.3f SD per SD" % mde_new)
print("  improvement             : %.1f%%  (absolute %.4f SD)" % (
    100.0*(mde_now-mde_new)/mde_now, mde_now-mde_new))

# also: what if a denser panel doubled the instrument set (unlikely upper bound)?
for k in (7, 13, 20, 42):
    s = sum7 * (k/7.0)
    m = mde_now * math.sqrt(sum42/(sum42+s))
    print("  hypothetical +%-2d instruments -> MDE %.4f SD (-%.1f%%)" % (
        k, m, 100.0*(mde_now-m)/mde_now))

print("\nDONE step7")
