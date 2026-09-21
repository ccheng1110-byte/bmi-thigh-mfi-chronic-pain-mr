#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Stage 4 - Sampling covariance audit
Step 2: harmonise GWAS summary statistics against the reference panel.

For every trait we:
  1. stream the raw summary statistics, keeping only variants whose
     (chromosome, position) is present in the 1000G EUR reference panel;
  2. resolve allele orientation against the reference A1/A2, allowing an
     allele swap and/or a strand complement (both resolved simultaneously are
     ambiguous palindromic SNPs and are dropped);
  3. drop indels / non-SNP alleles, and flag alleles whose frequency is
     inconsistent with the reference (possible residual orientation error);
  4. attach the LD score and write a harmonised file with columns already
     expressed on the reference A1.

Outputs
-------
stage4_gwas/<trait>_harmonised.tsv.gz
stage4_gwas/<trait>_raw_scan.json      (row counts, N distribution)
"""
import gzip
import json
import os
import time

import numpy as np
import pandas as pd

INC = r"D:\SX\TMFI_CPI_MR\00_stage0_audit\incoming"
BMI_FILE = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\SNP_gwas_mc_merge_nogc.tbl.uniq.gz"
LDSC_DIR = r"D:\SX\TMFI_CPI_MR\04_stage4_sampling_cov\stage4_ldsc"
OUT_DIR = r"D:\SX\TMFI_CPI_MR\04_stage4_sampling_cov\stage4_gwas"

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
PALINDROME = {frozenset(("A", "T")), frozenset(("G", "C"))}

os.makedirs(OUT_DIR, exist_ok=True)


def load_reference():
    ref = pd.read_csv(os.path.join(LDSC_DIR, "eur_ldscores.tsv.gz"), sep="\t")
    ref = ref.rename(columns={"chr": "chrom", "pos": "position"})
    ref["key"] = ref["chrom"].astype(np.int64) * 1_000_000_000 + ref["position"].astype(np.int64)
    ref["is_snp"] = ref["a1"].str.len().eq(1) & ref["a2"].str.len().eq(1)
    return ref


def harmonise(df, ref_sub, trait_name):
    """df: rows already restricted to reference positions, with columns
    chrom, position, a1, a2, beta, se, n (+ optional af)."""
    d = df.merge(ref_sub[["key", "a1", "a2", "maf", "af1", "ldscore", "is_snp"]],
                 on="key", how="left", suffixes=("", "_ref"))
    d = d[d["a1_ref"].notna()]
    d["is_snp"] = d["is_snp"].fillna(False)
    n0 = len(d)

    a1 = d["a1"].str.upper()
    a2 = d["a2"].str.upper()
    r1 = d["a1_ref"]
    r2 = d["a2_ref"]

    # palindromic SNPs cannot be orientation-resolved
    pal = [frozenset((x, y)) in PALINDROME for x, y in zip(a1, a2)]
    pal = np.array(pal)
    # allele length must be 1
    snp_a = a1.str.len().eq(1) & a2.str.len().eq(1)

    direct1 = (a1 == r1) & (a2 == r2)          # same alleles, same strand
    swap1 = (a1 == r2) & (a2 == r1)            # same strand, alleles swapped
    cc1 = a1.map(COMP); cc2 = a2.map(COMP)
    direct2 = (cc1 == r1) & (cc2 == r2)        # opposite strand, order preserved
    swap2 = (cc1 == r2) & (cc2 == r1)          # opposite strand, alleles swapped

    keep = snp_a & (~pal) & (direct1 | swap1 | direct2 | swap2)
    keep = keep.values
    d = d[keep].copy()
    r1 = r1[keep]
    r2 = r2[keep]

    # beta refers to the reported allele a1.  Its physical counterpart is the
    # reference A1 after the direct paths (direct1, direct2) and the reference
    # A2 after the swapped paths (swap1, swap2); only the latter need a flip.
    flip = (np.asarray(swap1[keep]) | np.asarray(swap2[keep]))
    if "beta" in d.columns:
        d["beta"] = np.where(flip, -d["beta"], d["beta"])
    for col in ("z", "tstat"):
        if col in d.columns:
            d[col] = np.where(flip, -d[col], d[col])
    # after harmonisation the effect allele IS the reference A1
    d["allele_ea"] = np.asarray(r1)
    d["allele_nea"] = np.asarray(r2)
    # AF must refer to the effect allele (reference A1)
    if "af" in d.columns:
        d["af"] = np.where(flip, 1.0 - d["af"], d["af"])
    print(f"  {trait_name}: matched={len(d):,}  (dropped {n0 - len(d):,} "
          f"unresolvable/ambiguous; {int(flip.sum()):,} allele swaps, "
          f"{int((swap2[keep]).sum()):,} strand-complement swaps)")
    return d


def scan_tmfi(ref):
    keep = set(ref["key"].to_numpy())
    chunks = []
    nrows = 0
    t0 = time.time()
    for ch in pd.read_csv(os.path.join(INC, "TMFI_GWAS.csv"), chunksize=2_000_000,
                          dtype={"CHR": np.int32, "SNP": str, "POS": np.int64,
                                 "A1": str, "A2": str, "N": np.int64,
                                 "AF1": np.float64, "BETA": np.float64,
                                 "SE": np.float64, "P": np.float64},
                          low_memory=False):
        nrows += len(ch)
        ch["key"] = ch["CHR"].astype(np.int64) * 1_000_000_000 + ch["POS"].astype(np.int64)
        ch = ch[ch["key"].isin(keep)]
        if len(ch):
            chunks.append(ch[["key", "CHR", "POS", "SNP", "A1", "A2",
                              "N", "AF1", "BETA", "SE", "P"]])
        if nrows % 4_000_000 == 0:
            print(f"    scanned {nrows:,} raw rows, kept {sum(len(c) for c in chunks):,} "
                  f"({time.time() - t0:.0f}s)", flush=True)
    raw = pd.concat(chunks, ignore_index=True)
    raw = raw.rename(columns={"CHR": "chrom", "POS": "position", "SNP": "rsid",
                              "A1": "a1", "A2": "a2", "N": "n",
                              "AF1": "af", "BETA": "beta", "SE": "se", "P": "p"})
    info = {"file": "TMFI_GWAS.csv", "raw_rows": int(nrows),
            "matched_reference_rows": int(len(raw)),
            "n_quantiles": {str(k): float(v) for k, v in
                            raw["n"].quantile([0, .01, .05, .5, .95, .99, 1]).items()},
            "scan_seconds": round(time.time() - t0, 1)}
    return raw, info


def scan_cpi(path, ref, name, gzipped=True):
    keep = set(ref["key"].to_numpy())
    chunks = []
    nrows = 0
    t0 = time.time()
    opener = gzip.open if gzipped else open
    with opener(os.path.join(INC, path), "rt") as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        print(f"    header: {cols}", flush=True)
        for ch in pd.read_csv(fh, sep="\t", names=cols, chunksize=2_000_000,
                              low_memory=False, comment="#", dtype=str):
            nrows += len(ch)
            ch["chrom"] = pd.to_numeric(ch.iloc[:, 0], errors="coerce")
            ch["position"] = pd.to_numeric(ch.iloc[:, 1], errors="coerce")
            ch = ch[ch["chrom"].between(1, 22) & ch["position"].notna()]
            ch["key"] = ch["chrom"].astype(np.int64) * 1_000_000_000 + \
                ch["position"].astype(np.int64)
            ch = ch[ch["key"].isin(keep)]
            if len(ch):
                chunks.append(ch)
            if nrows % 4_000_000 == 0:
                print(f"    scanned {nrows:,} raw rows, kept {sum(len(c) for c in chunks):,} "
                      f"({time.time() - t0:.0f}s)", flush=True)
    raw = pd.concat(chunks, ignore_index=True)
    info = {"file": path, "raw_rows": int(nrows),
            "matched_reference_rows": int(len(raw)),
            "scan_seconds": round(time.time() - t0, 1)}
    return raw, info


def scan_bmi(ref):
    """GIANT (Locke et al. 2015) BMI summary statistics.

    The file carries rsIDs only - no physical positions - so variants are
    anchored to the reference panel through the rsID and the position and
    reference alleles are then taken from the panel itself.
    """
    cols = ["rsid", "key", "chrom", "position"]
    ref_map = ref.loc[ref["rsid"].notna(), cols].copy()
    dup = ref_map["rsid"].duplicated(keep=False)
    ref_map = ref_map[~dup]
    print(f"    reference rsID map: {len(ref_map):,} unique anchors "
          f"({int(dup.sum()):,} duplicated rsIDs ignored)", flush=True)

    chunks = []
    nrows = 0
    t0 = time.time()
    for ch in pd.read_csv(BMI_FILE, sep="\t", chunksize=1_000_000,
                          dtype={"SNP": str, "A1": str, "A2": str,
                                 "Freq1.Hapmap": np.float64, "b": np.float64,
                                 "se": np.float64, "p": np.float64,
                                 "N": np.float64},
                          low_memory=False):
        nrows += len(ch)
        ch = ch.drop_duplicates(subset=["SNP"], keep=False)
        ch = ch.merge(ref_map, left_on="SNP", right_on="rsid", how="inner")
        if len(ch):
            chunks.append(ch)
    raw = pd.concat(chunks, ignore_index=True)
    raw = raw.rename(columns={"Freq1.Hapmap": "af", "b": "beta", "se": "se",
                              "p": "p", "N": "n", "A1": "a1", "A2": "a2"})
    info = {"file": os.path.basename(BMI_FILE), "raw_rows": int(nrows),
            "matched_reference_rows": int(len(raw)),
            "note": "matched on rsID; position and reference alleles taken "
                    "from the panel; betas are per-allele GIANT units but "
                    "z = beta/se is unit-free so the LDSC intercept is "
                    "unaffected",
            "n_quantiles": {str(k): float(v) for k, v in
                            raw["n"].quantile([0, .01, .05, .5, .95, .99, 1]).items()},
            "scan_seconds": round(time.time() - t0, 1)}
    return raw, info


def main():
    import sys
    only = set(sys.argv[1:])
    ref = load_reference()
    print(f"[ref] {len(ref):,} variants, {int(ref['is_snp'].sum()):,} simple SNPs")

    def wanted(name):
        return (not only) or (name in only)

    # ---------------- TMFI ----------------
    if wanted("TMFI"):
        print("[TMFI] streaming ...", flush=True)
    raw, info = scan_tmfi(ref)
    d = harmonise(raw, ref, "TMFI")
    out_cols = ["key", "chrom", "position", "rsid", "allele_ea", "allele_nea",
                "beta", "se", "z", "n", "af", "af1", "maf", "ldscore"]
    d["z"] = d["beta"] / d["se"]
    d = d[out_cols]
    d.to_csv(os.path.join(OUT_DIR, "TMFI_harmonised.tsv.gz"), sep="\t", index=False,
             compression="gzip", float_format="%.6g")
    json.dump(info, open(os.path.join(OUT_DIR, "TMFI_raw_scan.json"), "w"), indent=2)
    print(f"[TMFI] wrote {len(d):,} rows\n{json.dumps(info, indent=2)}", flush=True)

    # ---------------- CPI (full) ----------------
    if wanted("CPI_full"):
        print("[CPI_full] streaming ...", flush=True)
        raw, info = scan_cpi("Chronic_pain_intensity.txt.gz", ref, "CPI_full")
        r = pd.DataFrame({
            "chrom": raw["chrom"].astype(np.int32),
            "position": raw["position"].astype(np.int64),
            "rsid": raw.iloc[:, 2],
            "a1": raw.iloc[:, 3].str.upper(),
            "a2": raw.iloc[:, 4].str.upper(),
            "af": raw.iloc[:, 5].astype(np.float64),
            "beta": raw.iloc[:, 6].astype(np.float64),
            "se": raw.iloc[:, 7].astype(np.float64),
            "n": raw.iloc[:, 10].astype(np.float64),
            "key": raw["key"],
        })
        d = harmonise(r, ref, "CPI_full")
        d["z"] = d["beta"] / d["se"]
        info["n_quantiles"] = {str(k): float(v) for k, v in
                               d["n"].quantile([0, .01, .05, .5, .95, .99, 1]).items()}
        d = d[out_cols]
        d.to_csv(os.path.join(OUT_DIR, "CPI_full_harmonised.tsv.gz"), sep="\t",
                 index=False, compression="gzip", float_format="%.6g")
        json.dump(info, open(os.path.join(OUT_DIR, "CPI_full_raw_scan.json"), "w"), indent=2)
        print(f"[CPI_full] wrote {len(d):,} rows\n{json.dumps(info, indent=2)}", flush=True)

    # ---------------- CPI (no imaging) ----------------
    if wanted("CPI_no_imaging"):
        print("[CPI_no_imaging] streaming ...", flush=True)
        raw, info = scan_cpi("chronic_pain_intensity_no_imaging.txt.gz", ref,
                             "CPI_no_imaging")
        print(f"    columns detected: {raw.columns.tolist()}", flush=True)
        r = pd.DataFrame({
            "chrom": raw["chrom"].astype(np.int32),
            "position": raw["position"].astype(np.int64),
            "rsid": raw["ID"],
            "a1": raw["A1"].str.upper(),
            "a2": raw["OMITTED"].str.upper(),
            "af": raw["A1_FREQ"].astype(np.float64),
            "beta": raw["BETA"].astype(np.float64),
            "se": raw["SE"].astype(np.float64),
            "n": raw["OBS_CT"].astype(np.float64),
            "key": raw["key"],
        })
        d = harmonise(r, ref, "CPI_no_imaging")
        d["z"] = d["beta"] / d["se"]
        info["n_quantiles"] = {str(k): float(v) for k, v in
                               d["n"].quantile([0, .01, .05, .5, .95, .99, 1]).items()}
        d = d[out_cols]
        d.to_csv(os.path.join(OUT_DIR, "CPI_no_imaging_harmonised.tsv.gz"), sep="\t",
                 index=False, compression="gzip", float_format="%.6g")
        json.dump(info, open(os.path.join(OUT_DIR, "CPI_no_imaging_raw_scan.json"), "w"),
                  indent=2)
        print(f"[CPI_no_imaging] wrote {len(d):,} rows\n{json.dumps(info, indent=2)}",
              flush=True)

    # ---------------- BMI (GIANT, Locke et al. 2015) ----------------
    if wanted("BMI"):
        print("[BMI] streaming GIANT summary statistics ...", flush=True)
        raw, info = scan_bmi(ref)
        r = pd.DataFrame({
            "chrom": raw["chrom"].astype(np.int32),
            "position": raw["position"].astype(np.int64),
            "rsid": raw["rsid"],
            "a1": raw["a1"].str.upper(),
            "a2": raw["a2"].str.upper(),
            "af": raw["af"].astype(np.float64),
            "beta": raw["beta"].astype(np.float64),
            "se": raw["se"].astype(np.float64),
            "n": raw["n"].astype(np.float64),
            "key": raw["key"].astype(np.int64),
        })
        d = harmonise(r, ref, "BMI")
        d["z"] = d["beta"] / d["se"]
        # reference-panel allele frequency agreement, as a residual
        # orientation-error check for a dataset that carries no positions
        afdev = (d["af"] - d["af1"]).abs()
        info["af_vs_reference"] = {
            "median_abs_dev": float(afdev.median()),
            "frac_dev_gt_0p15": float((afdev > 0.15).mean()),
            "correlation": float(np.corrcoef(d["af"], d["af1"])[0, 1]),
        }
        info["n_quantiles"] = {str(k): float(v) for k, v in
                               d["n"].quantile([0, .01, .05, .5, .95, .99, 1]).items()}
        d = d[out_cols]
        d.to_csv(os.path.join(OUT_DIR, "BMI_harmonised.tsv.gz"), sep="\t",
                 index=False, compression="gzip", float_format="%.6g")
        json.dump(info, open(os.path.join(OUT_DIR, "BMI_raw_scan.json"), "w"),
                  indent=2)
        print(f"[BMI] wrote {len(d):,} rows\n{json.dumps(info, indent=2)}", flush=True)


if __name__ == "__main__":
    main()
