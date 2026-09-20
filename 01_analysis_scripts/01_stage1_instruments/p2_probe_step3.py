# -*- coding: utf-8 -*-
"""P2 step 3: tagging test -- are the panel-missing candidates already
covered by the 42 retained instruments?

Logic
-----
A missing candidate represents a LOST instrument only if it is NOT in LD
(r2 < 0.001) with any retained instrument.  We cannot compute LD for a
variant that is absent from the panel, so we use the nearest panel variant
as a positional proxy (panel median inter-variant gap is ~807 bp, and 99.9%
of missing candidates have a panel variant within 10 kb).

If the proxy has r2 > 0.001 with a retained instrument, the missing candidate
sits inside that instrument's LD block and would have been clumped away even
with a perfect panel -> no loss.
If no proxy-instrument pair reaches r2 > 0.001, the missing candidate marks a
signal that the current instrument set does not tag -> possible extra
instrument.
"""
import io, os, csv, bisect, subprocess, sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(ROOT)
REF = os.path.join(REPO, "reference", "1000G_EUR")
BIM = REF + ".bim"
PLINK = os.path.join(REPO, "tools", "plink19", "plink.exe")
WORK = os.path.join(ROOT, "p2_ldwork")
os.makedirs(WORK, exist_ok=True)

print("=" * 68)
print("P2 STEP 3 : tagging test (missing candidates vs 42 instruments)")
print("=" * 68)

# ---------- panel index ----------
panel_rs = {}                       # rsid -> (chr,pos)
panel_rows = defaultdict(list)      # chr -> [(pos, rsid)]
with io.open(BIM, "r", encoding="utf-8") as f:
    for line in f:
        p = line.split()
        panel_rs[p[1]] = (p[0], int(p[3]))
        panel_rows[p[0]].append((int(p[3]), p[1]))
for c in panel_rows:
    panel_rows[c].sort()
panel_pos = {c: [x[0] for x in v] for c, v in panel_rows.items()}

# ---------- candidates ----------
rows = []
with io.open(os.path.join(ROOT, "tmfi_gws_exact_clean.csv"), "r", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        rows.append({"rs": r["SNP"].strip(), "chr": r["CHR"].strip(),
                     "pos": int(float(r["POS"])), "p": float(r["P"]),
                     "maf": float(r["MAF"])})
missing = [r for r in rows if r["rs"] not in panel_rs]
print("missing candidates        : %d" % len(missing))

# ---------- instruments ----------
instr = []
with io.open(os.path.join(ROOT, "tmfi_gws_ldclumped_r2_0p001_kb10000.csv"), "r", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        instr.append(r["SNP"].strip())
print("retained instruments      : %d" % len(instr))
assert all(s in panel_rs for s in instr), "instrument not in panel!"

# ---------- nearest panel proxy for each missing candidate ----------
proxy = {}          # missing rs -> (proxy rs, distance bp)
for r in missing:
    c, pos = r["chr"], r["pos"]
    lst = panel_pos.get(c)
    if not lst:
        continue
    i = bisect.bisect_left(lst, pos)
    best, bd = None, None
    if i < len(lst):
        best, bd = panel_rows[c][i][1], lst[i] - pos
    if i > 0:
        d = pos - lst[i-1]
        if bd is None or d < bd:
            best, bd = panel_rows[c][i-1][1], d
    proxy[r["rs"]] = (best, bd)

dists = sorted(v[1] for v in proxy.values())
print("proxy distance bp         : median=%d  p90=%d  max=%d" % (
    dists[len(dists)//2], dists[int(len(dists)*0.9)], dists[-1]))
print("proxies within 10 kb      : %d / %d" % (
    sum(1 for v in proxy.values() if v[1] <= 10000), len(proxy)))

proxy_ids = sorted(set(v[0] for v in proxy.values() if v[0]))
print("unique proxy variants     : %d" % len(proxy_ids))

# ---------- write plink inputs ----------
allset = sorted(set(instr) | set(proxy_ids))
f_all = os.path.join(WORK, "extract.txt")
f_ins = os.path.join(WORK, "instruments.txt")
with io.open(f_all, "w", encoding="utf-8") as f:
    f.write("\n".join(allset) + "\n")
with io.open(f_ins, "w", encoding="utf-8") as f:
    f.write("\n".join(instr) + "\n")

out = os.path.join(WORK, "tag")
cmd = [PLINK, "--bfile", REF, "--extract", f_all, "--r2",
       "--ld-snp-list", f_ins,
       "--ld-window-kb", "10000", "--ld-window-r2", "0.001",
       "--out", out]
print("\nrunning PLINK ...")
r = subprocess.run(cmd, capture_output=True, text=True)
print("plink exit                : %d" % r.returncode)
if r.returncode != 0:
    print(r.stdout[-3000:])
    print(r.stderr[-3000:])
    sys.exit(1)

# ---------- parse ----------
ldfile = out + ".ld"
pairs = defaultdict(float)      # proxy rs -> max r2 to any instrument
nrow = 0
with io.open(ldfile, "r", encoding="utf-8") as f:
    hdr = f.readline()
    for line in f:
        p = line.split()
        if len(p) < 7:
            continue
        nrow += 1
        a, b, r2 = p[2], p[5], float(p[6])
        # a is the ld-snp (instrument); b is the partner
        if r2 > pairs.get(b, -1):
            pairs[b] = r2
print("LD rows reported (r2>=0.001): %d" % nrow)

tagged = 0
untagged = []
for rs, (px, d) in proxy.items():
    if px in pairs and pairs[px] >= 0.001:
        tagged += 1
    else:
        untagged.append(rs)
print("\n--- result ---")
print("missing candidates TAGGED by a retained instrument (r2>0.001) : %d (%.1f%%)"
      % (tagged, 100.0*tagged/len(proxy)))
print("missing candidates NOT tagged -> candidate lost instruments   : %d (%.1f%%)"
      % (len(untagged), 100.0*len(untagged)/len(proxy)))

# ---------- cluster untagged into loci ----------
um = [r for r in rows if r["rs"] in set(untagged)]
um.sort(key=lambda r: (int(r["chr"]), r["pos"]))
loci, cur = [], [um[0]] if um else []
for r in um[1:]:
    if r["chr"] == cur[-1]["chr"] and r["pos"] - cur[-1]["pos"] <= 500000:
        cur.append(r)
    else:
        loci.append(cur); cur = [r]
if cur:
    loci.append(cur)
print("\nuntagged candidates cluster into %d loci (500 kb gap) -> upper bound of +%d instruments"
      % (len(loci), len(loci)))
print("\nchr   n   best_P      best_rsid        MAF     F(z^2)")
import math
for L in sorted(loci, key=lambda L: min(x["p"] for x in L)):
    b = min(L, key=lambda x: x["p"])
    z = abs(0.0)
    print("%-4s %3d  %.3e  %-15s  %.3f" % (b["chr"], len(L), b["p"], b["rs"], b["maf"]))

# ---------- write report ----------
rep = os.path.join(ROOT, "p2_tagging_result.txt")
with io.open(rep, "w", encoding="utf-8") as f:
    f.write("missing_candidates\t%d\n" % len(missing))
    f.write("tagged_by_instrument\t%d\n" % tagged)
    f.write("untagged\t%d\n" % len(untagged))
    f.write("untagged_loci_500kb\t%d\n" % len(loci))
print("\nwritten: %s" % rep)
print("DONE step3")
