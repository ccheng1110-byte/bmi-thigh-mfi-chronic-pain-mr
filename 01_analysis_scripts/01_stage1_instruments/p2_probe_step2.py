# -*- coding: utf-8 -*-
"""P2 step 2: locus-level analysis.

The decisive question is NOT "how many candidate SNPs are missing" but
"how many independent LOCI are invisible to the clumping step".

Method:
  1. cluster all 3,838 candidates into loci by physical gap thresholds
  2. for each locus ask: does it contain >=1 variant present in the panel?
  3. loci with zero panel variants are loci where clumping was blind
     -> upper bound on how many extra instruments a denser panel could add
  4. same analysis restricted to the strongest locus representatives
"""
import io, os, csv, bisect, sys
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(ROOT)
BIM = os.path.join(REPO, "reference", "1000G_EUR.bim")
CAND = os.path.join(ROOT, "tmfi_gws_exact_clean.csv")
CLUMPED = os.path.join(ROOT, "tmfi_gws_ldclumped_r2_0p001_kb10000.csv")
PLINK_LOG = os.path.join(ROOT, "plink_tmfi_r2_0p001_kb10000.log")

print("=" * 68)
print("P2 STEP 2 : locus-level blind-spot analysis")
print("=" * 68)

panel = set()
panel_pos = defaultdict(list)
with io.open(BIM, "r", encoding="utf-8") as f:
    for line in f:
        p = line.split()
        panel.add(p[1])
        panel_pos[p[0]].append(int(p[3]))
for c in panel_pos:
    panel_pos[c].sort()

rows = []
with io.open(CAND, "r", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        rows.append({
            "rs": r["SNP"].strip(),
            "chr": r["CHR"].strip(),
            "pos": int(float(r["POS"])),
            "p": float(r["P"]),
            "maf": float(r["MAF"]),
            "beta": float(r["BETA"]),
            "se": float(r["SE"]),
        })
print("candidates                : %d" % len(rows))

# markers used as PLINK input?
infile = os.path.join(ROOT, "tmfi_gws_for_plink.tsv")
n_plink_input = 0
if os.path.exists(infile):
    with io.open(infile, "r", encoding="utf-8") as f:
        n_plink_input = sum(1 for _ in f) - 1
    print("PLINK input markers       : %d" % n_plink_input)

# ---------- PLINK log: how many top variants were warned missing ----------
warn = 0
if os.path.exists(PLINK_LOG):
    with io.open(PLINK_LOG, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if "is missing from the main dataset" in line:
                warn += 1
    print("PLINK 'missing top variant' warnings : %d" % warn)

# ---------- locus clustering ----------
def cluster(gap_bp):
    """return list of loci; each locus = list of candidate indices"""
    loci = []
    for c in sorted(set(r["chr"] for r in rows), key=lambda x: int(x)):
        sub = sorted([r for r in rows if r["chr"] == c], key=lambda r: r["pos"])
        cur = [sub[0]]
        for r in sub[1:]:
            if r["pos"] - cur[-1]["pos"] > gap_bp:
                loci.append(cur)
                cur = [r]
            else:
                cur.append(r)
        loci.append(cur)
    return loci

print("\n--- locus clustering by physical gap ---")
print("gap(kb)   loci   loci_with_panel   blind_loci   %blind   max_extra_instruments")
summary = {}
for gapkb in (100, 250, 500, 1000, 2000, 5000):
    loci = cluster(gapkb * 1000)
    withp = [L for L in loci if any(x["rs"] in panel for x in L)]
    blind = [L for L in loci if not any(x["rs"] in panel for x in L)]
    # a blind locus can contribute at most 1 extra instrument if it is a
    # genuinely independent signal; count only blind loci whose best P is
    # comparable to what already got in (i.e. below the genome-wide threshold,
    # which all candidates satisfy by construction)
    print("%7d  %5d   %10d       %8d   %5.1f%%   %d" % (
        gapkb, len(loci), len(withp), len(blind), 100.0*len(blind)/len(loci), len(blind)))
    summary[gapkb] = (len(loci), len(withp), len(blind))

# ---------- how many loci did the 42 instruments come from? ----------
cl = []
with io.open(CLUMPED, "r", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        cl.append(r)
print("\nclumped instruments       : %d" % len(cl))
print("clumped columns           : %s" % list(cl[0].keys()))

# ---------- blind loci detail (500 kb gap) ----------
gapkb = 500
loci = cluster(gapkb * 1000)
blind = [L for L in loci if not any(x["rs"] in panel for x in L)]
blind.sort(key=lambda L: -len(L))
print("\n--- blind loci detail (gap=%dkb), top 25 by size ---" % gapkb)
print("chr   span_kb   n_cand   best_P     best_rsid        MAF")
for L in blind[:25]:
    best = min(L, key=lambda r: r["p"])
    print("%-4s  %7.0f   %5d   %.2e   %-15s  %.3f" % (
        best["chr"], (L[-1]["pos"]-L[0]["pos"])/1000.0, len(L), best["p"], best["rs"], best["maf"]))

print("\ntotal candidates inside blind loci : %d (%.1f%% of all)" % (
    sum(len(L) for L in blind), 100.0*sum(len(L) for L in blind)/len(rows)))
print("blind loci with >=5 candidates     : %d" % sum(1 for L in blind if len(L) >= 5))
print("blind loci with 1 candidate        : %d" % sum(1 for L in blind if len(L) == 1))
print("blind loci, best P < 1e-20         : %d" % sum(1 for L in blind if min(x["p"] for x in L) < 1e-20))
print("blind loci, best P < 1e-30         : %d" % sum(1 for L in blind if min(x["p"] for x in L) < 1e-30))

# ---------- strong-locus view: only loci whose best P < 5e-8 (all) ----------
print("\n--- interpretation helper ---")
for gapkb in (250, 500, 1000):
    n, wp, bl = summary[gapkb]
    print("gap %5dkb: %d loci total; %d visible to clumping; %d blind (upper bound of +%d instruments, %.0f%% more)"
          % (gapkb, n, wp, bl, bl, 100.0*bl/max(wp, 1)))
print("\nDONE step2")
