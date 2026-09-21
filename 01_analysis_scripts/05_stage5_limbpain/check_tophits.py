"""Verify the identity of the downloaded limb-pain GWAS files by comparing
their genome-wide significant top hits against each other and against known
locus fingerprints (e.g. GDF5 for knee OA / knee pain).

Files:
  ukb-b-9838.vcf.gz   -- claimed Back pain  (OpenGWAS, UKB, EUR)
  ukb-b-16254.vcf.gz  -- claimed Knee pain  (OpenGWAS, UKB, EUR)
  7vsno7.txt.gz       -- figshare knee pain GWAS (n=429,737)
"""
import gzip
import os

D = r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain\incoming"
THRESH = 7.0          # -log10(p) >= 7  => p <= 1e-7
KEEP_PER_CHR = 40

# known locus fingerprints (GRCh37)
FINGER = {
    "rs143383": ("20", 34025733, "GDF5 - knee OA / knee development"),
    "rs7562140": ("1", 0, "placeholder"),
}


def scan_vcf(path, tag):
    hits = []
    finger = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            fmt = f[8].split(":")
            val = f[9].split(":")
            try:
                lp = float(val[fmt.index("LP")])
                es = float(val[fmt.index("ES")])
            except (ValueError, IndexError):
                continue
            rid = f[2]
            if rid in FINGER:
                finger[rid] = (f[0], int(f[1]), es, lp)
            if lp >= THRESH:
                hits.append((lp, f[0], int(f[1]), rid, es))
    hits.sort(reverse=True)
    # thin: keep at most KEEP_PER_CHR per chromosome so the list is readable
    per = {}
    thin = []
    for h in hits:
        c = per.get(h[1], 0)
        if c < KEEP_PER_CHR:
            thin.append(h)
            per[h[1]] = c + 1
    print(f"\n===== {tag}  ({os.path.basename(path)}) =====")
    print(f"   variants with p <= 1e-7 : {len(hits):,}")
    print("   top independent-ish hits (top 30 by p):")
    for lp, c, pos, rid, es in thin[:30]:
        print(f"      chr{c}:{pos:>9}  {rid:<16} ES={es:+.4f}  -log10P={lp:.2f}")
    print("   locus fingerprints:")
    for rid, (c, pos, desc) in FINGER.items():
        if rid in finger:
            fc, fp, es, lp = finger[rid]
            print(f"      {rid:<12} chr{fc}:{fp}  ES={es:+.4f}  -log10P={lp:.2f}"
                  f"   ({desc})")
        else:
            print(f"      {rid:<12} not found in file")
    return {h[3]: h[0] for h in hits}


def scan_flat(path, tag):
    hits = []
    finger = {}
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ic = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or f[0] == "CHR":
                continue
            try:
                p = float(f[ic["P"]])
                beta = float(f[ic["BETA"]])
                if p <= 0:
                    continue
                lp = -__import__("math").log10(p)
            except (ValueError, IndexError):
                continue
            rid = f[ic["SNP"]]
            if rid in FINGER:
                finger[rid] = (f[ic["CHR"]], int(float(f[ic["POS"]])), beta, lp)
            if lp >= THRESH:
                hits.append((lp, f[ic["CHR"]], int(float(f[ic["POS"]])), rid, beta))
    hits.sort(reverse=True)
    per = {}
    thin = []
    for h in hits:
        c = per.get(h[1], 0)
        if c < KEEP_PER_CHR:
            thin.append(h)
            per[h[1]] = c + 1
    print(f"\n===== {tag}  ({os.path.basename(path)}) =====")
    print(f"   variants with p <= 1e-7 : {len(hits):,}")
    print("   top hits (top 30 by p):")
    for lp, c, pos, rid, es in thin[:30]:
        print(f"      chr{c}:{pos:>9}  {rid:<16} BETA={es:+.4f}  -log10P={lp:.2f}")
    print("   locus fingerprints:")
    for rid, (c, pos, desc) in FINGER.items():
        if rid in finger:
            fc, fp, es, lp = finger[rid]
            print(f"      {rid:<12} chr{fc}:{fp}  BETA={es:+.4f}  -log10P={lp:.2f}"
                  f"   ({desc})")
        else:
            print(f"      {rid:<12} not found in file")
    return {h[3]: h[0] for h in hits}


a = scan_vcf(os.path.join(D, "ukb-b-9838.vcf.gz"), "UKB-b:9838 (claimed BACK pain)")
b = scan_vcf(os.path.join(D, "ukb-b-16254.vcf.gz"), "UKB-b:16254 (claimed KNEE pain)")
c = scan_flat(os.path.join(D, "7vsno7.txt.gz"), "figshare 7vsno7 (claimed knee pain)")

print("\n===== pairwise overlap of genome-wide significant loci =====")
for n1, s1 in (("9838(back)", a), ("16254(knee?)", b), ("figshare(knee)", c)):
    for n2, s2 in (("9838(back)", a), ("16254(knee?)", b), ("figshare(knee)", c)):
        if n1 >= n2:
            continue
        inter = set(s1) & set(s2)
        print(f"   {n1:<15} x {n2:<15} shared rsids = {len(inter):,}")
