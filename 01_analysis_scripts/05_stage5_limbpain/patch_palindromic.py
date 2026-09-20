#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Recover the palindromic instruments that the generic harmoniser dropped.

For an A/T or G/C variant the allele pair is uninformative, so orientation is
resolved from allele frequency: the reported A1 is taken to be the allele whose
frequency in the GWAS sample is closest to the frequency of the reference-panel
allele (the same rule TwoSampleMR uses when `action = 2`).
"""
import gzip
import os
import shutil

import numpy as np
import pandas as pd

LDSC = r"D:\SX\TMFI_CPI_MR\04_stage4_sampling_cov\stage4_ldsc\eur_ldscores.tsv.gz"
S5 = r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain\gwas_harmonised"
STD = r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain\standardised"
INC = r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain\incoming"

INST1 = r"D:\SX\TMFI_CPI_MR\01_stage1_instruments\tmfi_gws_ldclumped_r2_0p001_kb10000.csv"
INST2 = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\mvmr_union_exposures_usable.csv"

RAW = {
    "BACK_PAIN_UKB": (os.path.join(STD, "BACK_PAIN_UKB_standard.tsv.gz"),
                      {"chr": 0, "snp": 2, "pos": 1, "a1": 3, "a2": 4,
                       "n": 5, "af": 6, "beta": 7, "se": 8, "p": 9}),
    "KNEE_PAIN_UKB": (os.path.join(STD, "KNEE_PAIN_UKB_standard.tsv.gz"),
                      {"chr": 0, "snp": 2, "pos": 1, "a1": 3, "a2": 4,
                       "n": 5, "af": 6, "beta": 7, "se": 8, "p": 9}),
    "KNEE_PAIN_FIGSHARE": (os.path.join(INC, "7vsno7.txt.gz"),
                           {"chr": 0, "snp": 1, "pos": 2, "a1": 3, "a2": 4,
                            "n": 5, "af": 6, "beta": 7, "se": 8, "p": 9}),
}

i1 = pd.read_csv(INST1)
i2 = pd.read_csv(INST2)
targets = set(i1.loc[i1["PALINDROMIC"].astype(bool), "RSID"].astype(str)) if "RSID" in i1 else set()
targets |= set(i1.loc[i1["PALINDROMIC"].astype(bool), "SNP"].astype(str))
i2_pal = [s for s, a, b in zip(i2["SNP"], i2["TMFI_A1"], i2["TMFI_A2"])
          if {str(a).upper(), str(b).upper()} in ({"A", "T"}, {"C", "G"})]
targets |= set(map(str, i2_pal))
targets = {t for t in targets if t.startswith("rs")}
print(f"palindromic instruments to recover ({len(targets)}): {sorted(targets)}")

ref = pd.read_csv(LDSC, sep="\t")
ref = ref.rename(columns={"chr": "chrom", "pos": "position"})
ref["key"] = ref["chrom"].astype(np.int64) * 1_000_000_000 + ref["position"].astype(np.int64)

for oc, (path, c) in RAW.items():
    found = []
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            if f[c["snp"]] in targets:
                found.append({k: f[j] for k, j in c.items()})
    d = pd.DataFrame(found)
    print(f"\n[{oc}] raw matches: {len(d)} / {len(targets)}")
    if not len(d):
        continue
    for col in ("chr", "pos", "n", "af", "beta", "se", "p"):
        d[col] = pd.to_numeric(d[col], errors="coerce")
    d["chrom"] = d["chr"].astype(int)
    d["position"] = d["pos"].astype(int)
    d["key"] = d["chrom"] * 1_000_000_000 + d["position"]
    d = d.rename(columns={"a1": "fa1", "a2": "fa2"})
    d = d.merge(ref[["key", "a1", "a2", "af1", "ldscore"]], on="key", how="left")
    d = d[d["a1"].notna()]
    fa1 = d["fa1"].str.upper()
    fa2 = d["fa2"].str.upper()
    a1 = d["a1"].str.upper()
    a2 = d["a2"].str.upper()
    f1 = d["af"].astype(float)
    same = (fa1 == a1) & (fa2 == a2)
    swap = (fa1 == a2) & (fa2 == a1)
    # palindromic -> decide by allele-frequency distance
    d_pal = ~same & ~swap
    dist_keep = (f1 - d["af1"]).abs()
    dist_flip = ((1 - f1) - d["af1"]).abs()
    flip = swap | (d_pal & (dist_flip < dist_keep))
    d["beta_h"] = np.where(flip, -d["beta"], d["beta"])
    d["af_h"] = np.where(flip, 1 - f1, f1)
    # Always reported on the reference panel's a1; see patch_exposures.py.
    d["allele_ea"] = a1
    d["allele_nea"] = a2
    d["z"] = d["beta_h"] / d["se"]
    d["maf"] = np.minimum(d["af_h"], 1 - d["af_h"])
    out = d[["key", "chrom", "position", "snp", "allele_ea", "allele_nea",
             "beta_h", "se", "z", "n", "af_h", "af1", "maf", "ldscore"]].copy()
    out.columns = ["key", "chrom", "position", "rsid", "allele_ea", "allele_nea",
                   "beta", "se", "z", "n", "af", "af1", "maf", "ldscore"]
    print(out[["rsid", "chrom", "position", "allele_ea", "allele_nea", "beta",
               "se", "af", "af1", "z"]].to_string(index=False))

    src = os.path.join(S5, f"{oc}_harmonised.tsv.gz")
    dst = os.path.join(S5, f"{oc}_harmonised_patched.tsv.gz")
    base = pd.read_csv(src, sep="\t")
    base = base[~base["key"].isin(set(out["key"]))]
    patched = pd.concat([base, out], ignore_index=True)
    patched.to_csv(dst, sep="\t", index=False, compression="gzip",
                   float_format="%.6g")
    print(f"  patched file: {dst}  rows {len(base):,} -> {len(patched):,}")
    shutil.copy2(dst, os.path.join(
        r"D:\SX\TMFI_CPI_MR\04_stage4_sampling_cov\stage4_gwas",
        f"{oc}_harmonised.tsv.gz"))
    print("  copied into stage4_gwas (for the LDSC audit)")
