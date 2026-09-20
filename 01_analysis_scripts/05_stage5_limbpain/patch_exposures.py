#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Recover the palindromic instruments on the EXPOSURE side (TMFI and BMI).

Stage 4 harmonisation dropped every A/T and G/C variant, so the frozen 42-SNP
and 93-SNP instrument sets could not be matched on the exposure side either.
Here the exposure betas are taken from the frozen instrument files themselves
and oriented against the reference panel by allele frequency.
"""
import os

import numpy as np
import pandas as pd

LDSC = r"D:\SX\TMFI_CPI_MR\04_stage4_sampling_cov\stage4_ldsc\eur_ldscores.tsv.gz"
S4 = r"D:\SX\TMFI_CPI_MR\04_stage4_sampling_cov\stage4_gwas"
S5 = r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain\gwas_harmonised"

INST1 = r"D:\SX\TMFI_CPI_MR\01_stage1_instruments\tmfi_gws_ldclumped_r2_0p001_kb10000.csv"
INST2 = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\mvmr_union_exposures_usable.csv"

PAL = ({"A", "T"}, {"C", "G"})
ref = pd.read_csv(LDSC, sep="\t").rename(columns={"chr": "chrom", "pos": "position"})
ref["key"] = ref["chrom"].astype(np.int64) * 1_000_000_000 + \
    ref["position"].astype(np.int64)

i1 = pd.read_csv(INST1)
i2 = pd.read_csv(INST2)


def is_pal(a, b):
    return {str(a).upper(), str(b).upper()} in PAL


# ---- TMFI rows ------------------------------------------------------------
rows = {}
for _, r in i1[i1["PALINDROMIC"].astype(bool)].iterrows():
    rows[str(r["SNP"])] = dict(snp=r["SNP"], chrom=int(r["CHR"]),
                               pos=int(r["POS"]), a1=r["A1"], a2=r["A2"],
                               eaf=r["AF1"], beta=r["BETA"], se=r["SE"],
                               n=r["N"])
for _, r in i2.iterrows():
    s = str(r["SNP"])
    if s in rows or not is_pal(r["TMFI_A1"], r["TMFI_A2"]):
        continue
    rows[s] = dict(snp=s, chrom=int(r["CHR"]), pos=int(r["POS"]),
                   a1=r["TMFI_A1"], a2=r["TMFI_A2"], eaf=r["TMFI_EAF"],
                   beta=r["TMFI_BETA"], se=r["TMFI_SE"], n=r["TMFI_N"])
tmfi_rows = pd.DataFrame(rows.values())

bmi_rows = []
for _, r in i2.iterrows():
    if not is_pal(r["BMI_A1"], r["BMI_A2"]):
        continue
    bmi_rows.append(dict(snp=r["SNP"], chrom=int(r["CHR"]), pos=int(r["POS"]),
                         a1=r["BMI_A1"], a2=r["BMI_A2"], eaf=r["BMI_EAF"],
                         beta=r["BMI_BETA"], se=r["BMI_SE"], n=r["BMI_N"]))
bmi_rows = pd.DataFrame(bmi_rows)


def orient(df, tag):
    df["key"] = df["chrom"].astype(np.int64) * 1_000_000_000 + \
        df["pos"].astype(np.int64)
    df = df.rename(columns={"a1": "ea1", "a2": "ea2"})
    d = df.merge(ref[["key", "a1", "a2", "af1", "ldscore"]], on="key", how="left")
    d = d[d["a1"].notna()].copy()
    eaf = d["eaf"].astype(float)
    flip = ((1 - eaf) - d["af1"]).abs() < (eaf - d["af1"]).abs()
    d["beta"] = np.where(flip, -d["beta"], d["beta"])
    d["af"] = np.where(flip, 1 - eaf, eaf)
    # The harmonised record is always reported on the reference panel's a1:
    # `af` above is the frequency of that allele by construction, so writing
    # allele_ea = a2 when the strand was flipped would label the variant with
    # the allele whose frequency is 1 - af and contradict the beta sign.
    d["allele_ea"] = d["a1"]
    d["allele_nea"] = d["a2"]
    d["z"] = d["beta"] / d["se"]
    d["maf"] = np.minimum(d["af"], 1 - d["af"])
    out = d[["key", "chrom", "pos", "snp", "allele_ea", "allele_nea", "beta",
             "se", "z", "n", "af", "af1", "maf", "ldscore"]].copy()
    out.columns = ["key", "chrom", "position", "rsid", "allele_ea", "allele_nea",
                   "beta", "se", "z", "n", "af", "af1", "maf", "ldscore"]
    print(f"\n[{tag}] recovered {len(out)} palindromic exposure variants")
    print(out[["rsid", "chrom", "position", "allele_ea", "allele_nea", "beta",
               "se", "af", "af1"]].to_string(index=False))
    return out


for tag, new, src in (("TMFI", orient(tmfi_rows, "TMFI"), "TMFI_harmonised.tsv.gz"),
                      ("BMI", orient(bmi_rows, "BMI"), "BMI_harmonised.tsv.gz")):
    base = pd.read_csv(os.path.join(S4, src), sep="\t")
    base = base[~base["key"].isin(set(new["key"]))]
    patched = pd.concat([base, new], ignore_index=True)
    dst = os.path.join(S5, f"{tag}_harmonised_patched.tsv.gz")
    patched.to_csv(dst, sep="\t", index=False, compression="gzip",
                   float_format="%.6g")
    print(f"  {src}: {len(base):,} -> {len(patched):,}  -> {dst}")
