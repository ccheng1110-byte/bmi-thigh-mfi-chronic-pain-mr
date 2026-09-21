# -*- coding: utf-8 -*-
"""P2 step 6: r2 noise floor of the 503-person panel, computed directly
from the .bed (avoids PLINK --r2 inter-chr quirks).

Also computes the empirical distribution of r2 among the 42 instruments
themselves, split into 'within 10 Mb same chromosome' (the pairs PLINK
actually evaluated during clumping) and 'unlinked' pairs.
"""
import io, os, csv, math
from collections import defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(ROOT)
REF = os.path.join(REPO, "reference", "1000G_EUR")

print("=" * 68)
print("P2 STEP 6 : r2 noise floor from raw genotypes")
print("=" * 68)

# ---------- read bim ----------
bim = []
idx = {}
with io.open(REF + ".bim", "r", encoding="utf-8") as f:
    for i, line in enumerate(f):
        p = line.split()
        bim.append((p[0], p[1], int(p[3])))
        idx[p[1]] = i
n_snp = len(bim)
with io.open(REF + ".fam", "r", encoding="utf-8") as f:
    n_ind = sum(1 for _ in f)
print("panel: %d variants x %d individuals" % (n_snp, n_ind))

instr = []
with io.open(os.path.join(ROOT, "tmfi_gws_ldclumped_r2_0p001_kb10000.csv"), "r", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        instr.append(r["SNP"].strip())
print("instruments: %d" % len(instr))

# ---------- read bed (SNP-major) ----------
bytes_per_snp = (n_ind + 3) // 4
want = {idx[s]: s for s in instr if s in idx}
print("instruments found in bed: %d" % len(want))

geno = {}
with open(REF + ".bed", "rb") as f:
    magic = f.read(3)
    assert magic == b"\x6c\x1b\x01", "not a SNP-major PLINK bed: %r" % magic
    for s_i in sorted(want):
        f.seek(3 + s_i * bytes_per_snp)
        buf = f.read(bytes_per_snp)
        g = []
        for b in buf:
            for shift in (0, 2, 4, 6):
                code = (b >> shift) & 3
                # PLINK bed 2-bit: 00=hom A1(0) 01=missing(1) 10=het(2) 11=hom A2(3)
                g.append(-1 if code == 1 else (0 if code == 0 else (1 if code == 2 else 2)))
        geno[want[s_i]] = g[:n_ind]
print("genotypes loaded: %d" % len(geno))

# ---- sanity check: allele frequency vs the GWAS MAF column ----
gw = {r["SNP"].strip(): (float(r["MAF"]), r["A1"].strip(), r["A2"].strip())
      for r in csv.DictReader(io.open(os.path.join(ROOT, "tmfi_gws_exact_clean.csv"),
                                      "r", encoding="utf-8"))}
print("\n--- decode sanity check: panel allele frequency vs GWAS MAF ---")
diffs = []
for s in sorted(geno):
    v = [x for x in geno[s] if x >= 0]
    miss = 1.0 - len(v)/float(len(geno[s]))
    if not v or s not in gw:
        continue
    # dosage of A1 (coded 0 = hom A1) -> freq of A2 allele
    fA2 = sum(v)/(2.0*len(v))
    fA1 = 1.0 - fA2
    maf_panel = min(fA1, fA2)
    maf_gwas = gw[s][0]
    diffs.append(abs(maf_panel - maf_gwas))
print("  n compared          : %d" % len(diffs))
print("  mean |MAF_panel-MAF_gwas| : %.4f" % (sum(diffs)/len(diffs)))
print("  max                 : %.4f" % max(diffs))
print("  missingness (mean)  : %.4f" % (sum(1.0 - len([x for x in geno[s] if x >= 0])/float(len(geno[s]))
                                           for s in geno)/len(geno)))


def r2(g1, g2):
    xs, ys = [], []
    for a, b in zip(g1, g2):
        if a < 0 or b < 0:
            continue
        xs.append(a); ys.append(b)
    n = len(xs)
    if n < 30:
        return None
    mx = sum(xs) / n; my = sum(ys) / n
    sxy = sxx = syy = 0.0
    for a, b in zip(xs, ys):
        dx = a - mx; dy = b - my
        sxy += dx * dy; sxx += dx * dx; syy += dy * dy
    if sxx <= 0 or syy <= 0:
        return None
    r = sxy / math.sqrt(sxx * syy)
    return r * r


pos = {r["SNP"].strip(): (r["CHR"].strip(), int(float(r["POS"])))
       for r in csv.DictReader(io.open(os.path.join(ROOT, "tmfi_gws_exact_clean.csv"),
                                       "r", encoding="utf-8"))}

unlinked, linked = [], []
names = sorted(geno)
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = names[i], names[j]
        v = r2(geno[a], geno[b])
        if v is None:
            continue
        ca, pa = pos[a]; cb, pb = pos[b]
        same_win = (ca == cb and abs(pa - pb) <= 10_000_000)
        (linked if same_win else unlinked).append((v, a, b))

un = sorted(x[0] for x in unlinked)
ln = sorted(x[0] for x in linked)
print("\n--- UNLINKED instrument pairs (different chromosome or >10 Mb) ---")
print("  n pairs        : %d" % len(un))
print("  mean   r2      : %.5f     (theory 1/(N-1) = %.5f)" % (sum(un)/len(un), 1.0/(n_ind-1)))
print("  median r2      : %.5f" % un[len(un)//2])
print("  p90 / p99 /max : %.5f / %.5f / %.5f" % (un[int(len(un)*.9)], un[int(len(un)*.99)], un[-1]))
print("  fraction with r2 > 0.001 (the clumping threshold) : %.1f%%"
      % (100.0*sum(1 for v in un if v > 0.001)/len(un)))
print("  fraction with r2 > 0.002                          : %.1f%%"
      % (100.0*sum(1 for v in un if v > 0.002)/len(un)))

print("\n--- instrument pairs INSIDE the clumping window (same chr, <=10 Mb) ---")
if ln:
    print("  n pairs        : %d" % len(ln))
    print("  median / max r2: %.5f / %.5f" % (ln[len(ln)//2], ln[-1]))
    print("  all below 0.001? %s" % ("YES (clumping succeeded)" if ln[-1] < 0.001 else "NO"))
else:
    print("  none (all 42 instruments are >10 Mb apart or on different chromosomes)")

print("\nINTERPRETATION")
print("  E[r2] for unlinked SNPs in a %d-person panel is 1/(N-1) = %.5f." % (n_ind, 1.0/(n_ind-1)))
if un:
    print("  Observed median %.5f, and %.0f%% of truly independent pairs already exceed"
          % (un[len(un)//2], 100.0*sum(1 for v in un if v > 0.001)/len(un)))
    print("  r2 = 0.001.  The --clump-r2 0.001 threshold therefore sits BELOW the")
    print("  sampling noise floor: pruning is driven partly by chance and is, if")
    print("  anything, over-conservative.  A denser panel would NOT relax this -")
    print("  the limit is the 503 individuals, not the variant count.")
print("\nDONE step6")
