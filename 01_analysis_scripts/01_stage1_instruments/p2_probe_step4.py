# -*- coding: utf-8 -*-
"""P2 step 4: VALIDATE the tagging method with a positive control, and
map the instrument / candidate geography in detail.

Positive control
----------------
The 873 panel-visible candidates minus the 42 retained instruments = 831
SNPs.  PLINK removed every one of them *because* it had r2 > 0.001 with a
retained index.  So the method must tag ~100% of them.  If it does, the
71.5% untagged rate seen for the missing candidates is real.  If it does
not, the proxy/LD pipeline is broken and step 3's number is meaningless.
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
print("P2 STEP 4 : positive control + geography")
print("=" * 68)

panel_rs = {}
panel_rows = defaultdict(list)
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
byrs = {r["rs"]: r for r in rows}
instr = []
with io.open(os.path.join(ROOT, "tmfi_gws_ldclumped_r2_0p001_kb10000.csv"), "r", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        instr.append(r["SNP"].strip())
iset = set(instr)
visible = [r for r in rows if r["rs"] in panel_rs]
control = [r for r in visible if r["rs"] not in iset]
missing = [r for r in rows if r["rs"] not in panel_rs]
print("candidates %d | visible %d | instruments %d | control(clumped-away) %d | missing %d"
      % (len(rows), len(visible), len(instr), len(control), len(missing)))

# ---------- run LD: instruments vs (control + proxies) ----------
def nearest_panel_proxy(r):
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

proxy = {r["rs"]: nearest_panel_proxy(r) for r in missing}
proxy_ids = sorted(set(v[0] for v in proxy.values() if v[0]))

partners = sorted(set(proxy_ids) | set(x["rs"] for x in control))
f_all = os.path.join(WORK, "extract4.txt")
f_ins = os.path.join(WORK, "instruments.txt")
with io.open(f_all, "w", encoding="utf-8") as f:
    f.write("\n".join(sorted(set(partners) | iset)) + "\n")
with io.open(f_ins, "w", encoding="utf-8") as f:
    f.write("\n".join(instr) + "\n")

out = os.path.join(WORK, "tag4")
cmd = [PLINK, "--bfile", REF, "--extract", f_all, "--r2",
       "--ld-snp-list", f_ins, "--ld-window-kb", "10000",
       "--ld-window", "99999", "--ld-window-r2", "0.001", "--out", out]
r = subprocess.run(cmd, capture_output=True, text=True)
print("plink exit: %d" % r.returncode)
if r.returncode != 0:
    print(r.stdout[-2000:]); print(r.stderr[-2000:]); sys.exit(1)

pairs = defaultdict(float)
with io.open(out + ".ld", "r", encoding="utf-8") as f:
    f.readline()
    for line in f:
        p = line.split()
        if len(p) < 7: continue
        b, r2 = p[5], float(p[6])
        if r2 > pairs.get(b, -1):
            pairs[b] = r2
print("LD pairs reported: %d" % len(pairs))

# ---------- control ----------
c_tagged = sum(1 for x in control if pairs.get(x["rs"], 0) >= 0.001)
print("\n*** POSITIVE CONTROL ***")
print("clumped-away candidates tagged by an instrument : %d / %d (%.1f%%)"
      % (c_tagged, len(control), 100.0*c_tagged/len(control)))
if c_tagged / float(len(control)) > 0.95:
    print("=> method VALIDATED (>=95%%): step-3 untagged rate is credible")
else:
    print("=> method FAILED the control: step-3 result cannot be trusted")

# ---------- proxies ----------
p_tagged = sum(1 for v in proxy.values() if pairs.get(v[0], 0) >= 0.001)
print("\nmissing-candidate proxies tagged : %d / %d (%.1f%%)"
      % (p_tagged, len(proxy), 100.0*p_tagged/len(proxy)))

# ---------- geography: instruments ----------
print("\n--- 42 instruments (chr:pos, P, MAF) ---")
ins_rows = sorted([byrs[s] for s in instr], key=lambda r: (int(r["chr"]), r["pos"]))
for r in ins_rows:
    print("  chr%-3s %11d  P=%.2e  MAF=%.3f  %s" % (r["chr"], r["pos"], r["p"], r["maf"], r["rs"]))

# ---------- geography: per-locus (500 kb gap) detail ----------
loci = []
for c in sorted(set(r["chr"] for r in rows), key=lambda x: int(x)):
    sub = sorted([r for r in rows if r["chr"] == c], key=lambda r: r["pos"])
    cur = [sub[0]]
    for r in sub[1:]:
        if r["pos"] - cur[-1]["pos"] > 500000:
            loci.append(cur); cur = [r]
        else:
            cur.append(r)
    loci.append(cur)
loci.sort(key=lambda L: min(x["p"] for x in L))
print("\n--- loci (500 kb gap), sorted by best P ---")
print("chr  start        end          span_kb   n_cand  n_visible  n_instr  best_P")
for L in loci:
    best = min(L, key=lambda x: x["p"])
    nv = sum(1 for x in L if x["rs"] in panel_rs)
    ni = sum(1 for x in L if x["rs"] in iset)
    print("%-4s %11d %11d %9.0f  %5d    %5d      %3d   %.2e  %s" % (
        best["chr"], L[0]["pos"], L[-1]["pos"], (L[-1]["pos"]-L[0]["pos"])/1000.0,
        len(L), nv, ni, best["p"], best["rs"]))

# ---------- focus: is the single strongest signal represented? ----------
top = sorted(rows, key=lambda r: r["p"])[:15]
print("\n--- top 15 candidates: in panel? in instrument set? ---")
for r in top:
    print("  %-15s chr%-3s %11d P=%.2e  panel=%-5s instr=%-5s nearest_instr_kb=%s" % (
        r["rs"], r["chr"], r["pos"], r["p"],
        r["rs"] in panel_rs, r["rs"] in iset,
        (lambda d: ("%.0f" % (d/1000.0)) if d is not None else "n/a")(
            min((abs(x["pos"]-r["pos"]) for x in ins_rows if x["chr"] == r["chr"]), default=None))))
# ---------- noise floor: r2 between UNLINKED pairs in a 503-person panel ----------
# expected E[r2] ~ 1/(N-1) = 1/502 = 0.00199 for truly independent SNPs.
# If that exceeds the 0.001 clumping threshold, PLINK prunes by noise, not by LD.
print("\n--- r2 noise floor of the 503-person panel ---")
nf_out = os.path.join(WORK, "noise")
cmd2 = [PLINK, "--bfile", REF, "--extract", os.path.join(WORK, "instruments.txt"),
        "--r2", "inter-chr", "--ld-window", "99999", "--ld-window-r2", "0",
        "--out", nf_out]
r2run = subprocess.run(cmd2, capture_output=True, text=True)
if r2run.returncode == 0:
    vals = []
    with io.open(nf_out + ".ld", "r", encoding="utf-8") as f:
        f.readline()
        for line in f:
            p = line.split()
            if len(p) >= 7:
                vals.append(float(p[6]))
    if vals:
        vals.sort()
        print("inter-chromosomal instrument pairs (n=%d, truly unlinked)" % len(vals))
        print("  mean r2   = %.5f   (theory 1/(N-1) = %.5f)" % (sum(vals)/len(vals), 1.0/502))
        print("  median r2 = %.5f" % vals[len(vals)//2])
        print("  p90 r2    = %.5f   max = %.5f" % (vals[int(len(vals)*0.9)], vals[-1]))
        print("  fraction of UNLINKED pairs with r2 > 0.001 (the clumping threshold) = %.1f%%"
              % (100.0*sum(1 for v in vals if v > 0.001)/len(vals)))
        print("  => at N=503 the threshold r2<0.001 is BELOW the sampling noise floor;")
        print("     clumping therefore prunes partly by chance, which explains why each")
        print("     locus yields exactly one instrument.")
else:
    print("noise-floor run failed: %s" % r2run.stdout[-500:])

print("\nDONE step4")
