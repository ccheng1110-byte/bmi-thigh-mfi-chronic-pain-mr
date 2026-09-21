# -*- coding: utf-8 -*-
"""P2 step 5: corrected tagging (after fixing --ld-window) + r2 noise floor."""
import io, os, csv, bisect, subprocess, sys, random
from collections import defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(ROOT)
REF = os.path.join(REPO, "reference", "1000G_EUR")
BIM = REF + ".bim"
PLINK = os.path.join(REPO, "tools", "plink19", "plink.exe")
WORK = os.path.join(ROOT, "p2_ldwork")
os.makedirs(WORK, exist_ok=True)

print("=" * 68)
print("P2 STEP 5 : corrected tagging + noise floor")
print("=" * 68)

panel_rows = defaultdict(list)
panel_rs = {}
with io.open(BIM, "r", encoding="utf-8") as f:
    for line in f:
        p = line.split()
        panel_rs[p[1]] = (p[0], int(p[3]))
        panel_rows[p[0]].append((int(p[3]), p[1]))
for c in panel_rows:
    panel_rows[c].sort()

rows = []
with io.open(os.path.join(ROOT, "tmfi_gws_exact_clean.csv"), "r", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        rows.append({"rs": r["SNP"].strip(), "chr": r["CHR"].strip(),
                     "pos": int(float(r["POS"])), "p": float(r["P"]),
                     "maf": float(r["MAF"]), "beta": float(r["BETA"]),
                     "se": float(r["SE"])})
instr = []
with io.open(os.path.join(ROOT, "tmfi_gws_ldclumped_r2_0p001_kb10000.csv"), "r", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        instr.append(r["SNP"].strip())
iset = set(instr)
missing = [r for r in rows if r["rs"] not in panel_rs]
visible = [r for r in rows if r["rs"] in panel_rs]

# ---------- noise floor: inter-chromosomal pairs among the 42 ----------
nf = os.path.join(WORK, "noise2")
cmd = [PLINK, "--bfile", REF, "--extract", os.path.join(WORK, "instruments.txt"),
       "--r2", "inter-chr", "--out", nf]
r = subprocess.run(cmd, capture_output=True, text=True)
vals = []
if r.returncode == 0 and os.path.exists(nf + ".ld"):
    with io.open(nf + ".ld", "r", encoding="utf-8") as f:
        f.readline()
        for line in f:
            p = line.split()
            if len(p) >= 7:
                vals.append(float(p[6]))
if vals:
    vals.sort()
    print("\n--- r2 noise floor (PLINK inter-chr) ---")
    print("  pairs                 : %d" % len(vals))
    print("  mean   r2             : %.5f   (theory 1/(N-1) = %.5f)" % (sum(vals)/len(vals), 1.0/502))
    print("  median r2             : %.5f" % vals[len(vals)//2])
else:
    print("\n(noise floor not available from PLINK here - see step6, computed from raw .bed)")

# ---------- corrected tagging ----------
def nearest_proxy(r):
    c, pos = r["chr"], r["pos"]
    lst = [x[0] for x in panel_rows[c]]
    i = bisect.bisect_left(lst, pos)
    best, bd = None, None
    if i < len(lst):
        best, bd = panel_rows[c][i][1], lst[i]-pos
    if i > 0:
        d = pos - lst[i-1]
        if bd is None or d < bd:
            best, bd = panel_rows[c][i-1][1], d
    return best, bd

proxy = {r["rs"]: nearest_proxy(r) for r in missing}
proxy_ids = sorted(set(v[0] for v in proxy.values() if v[0]))
f_all = os.path.join(WORK, "extract5.txt")
with io.open(f_all, "w", encoding="utf-8") as f:
    f.write("\n".join(sorted(set(proxy_ids) | iset)) + "\n")
out = os.path.join(WORK, "tag5")
cmd = [PLINK, "--bfile", REF, "--extract", f_all, "--r2",
       "--ld-snp-list", os.path.join(WORK, "instruments.txt"),
       "--ld-window-kb", "10000", "--ld-window", "99999",
       "--ld-window-r2", "0.001", "--out", out]
r = subprocess.run(cmd, capture_output=True, text=True)
print("\nplink exit %d" % r.returncode)
if r.returncode != 0:
    print(r.stdout[-1500:]); sys.exit(1)

pairs = defaultdict(float)
with io.open(out + ".ld", "r", encoding="utf-8") as f:
    f.readline()
    for line in f:
        p = line.split()
        if len(p) >= 7:
            b, r2 = p[5], float(p[6])
            if r2 > pairs.get(b, -1):
                pairs[b] = r2

tagged = [x for x in missing if pairs.get(proxy[x["rs"]][0], 0) >= 0.001]
untagged = [x for x in missing if pairs.get(proxy[x["rs"]][0], 0) < 0.001]
print("\n--- corrected tagging of panel-missing candidates ---")
print("  tagged by a retained instrument (r2>0.001) : %d (%.1f%%)" % (len(tagged), 100.0*len(tagged)/len(missing)))
print("  NOT tagged                                 : %d (%.1f%%)" % (len(untagged), 100.0*len(untagged)/len(missing)))

# cluster untagged into loci (500 kb gap)
um = sorted(untagged, key=lambda r: (int(r["chr"]), r["pos"]))
loci, cur = [], [um[0]] if um else []
for r in um[1:]:
    if r["chr"] == cur[-1]["chr"] and r["pos"] - cur[-1]["pos"] <= 500000:
        cur.append(r)
    else:
        loci.append(cur); cur = [r]
if cur: loci.append(cur)
print("\n  untagged candidates cluster into %d loci -> theoretical max +%d instruments" % (len(loci), len(loci)))
print("\n  chr   n   best_P      best_rsid       MAF    z~F")
for L in sorted(loci, key=lambda L: min(x["p"] for x in L)):
    b = min(L, key=lambda x: x["p"])
    z2 = (b["beta"]/b["se"])**2
    print("  %-4s %3d  %.3e  %-14s  %.3f  F=%.1f" % (b["chr"], len(L), b["p"], b["rs"], b["maf"], z2))

# how many untagged loci are actually NEW loci (no visible candidate at all)?
allc = sorted(rows, key=lambda r: (int(r["chr"]), r["pos"]))
big, cur = [], [allc[0]]
for r in allc[1:]:
    if r["chr"] == cur[-1]["chr"] and r["pos"] - cur[-1]["pos"] <= 500000:
        cur.append(r)
    else:
        big.append(cur); cur = [r]
big.append(cur)
blind = [L for L in big if not any(x["rs"] in panel_rs for x in L)]
print("\n--- cross-check with step-2 locus geography ---")
print("  total candidate loci (500kb gap)          : %d" % len(big))
print("  loci with >=1 panel-visible variant       : %d" % (len(big)-len(blind)))
print("  loci with 0 panel-visible variant (blind) : %d  -> hard upper bound +%d" % (len(blind), len(blind)))
print("  candidates living in blind loci           : %d" % sum(len(L) for L in blind))
print("\n  blind loci detail:")
for L in blind:
    b = min(L, key=lambda x: x["p"])
    print("    chr%-3s n=%d best_P=%.2e %s MAF=%.3f" % (b["chr"], len(L), b["p"], b["rs"], b["maf"]))
print("\nDONE step5")
