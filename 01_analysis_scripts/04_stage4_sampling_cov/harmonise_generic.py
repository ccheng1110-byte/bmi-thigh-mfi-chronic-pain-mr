#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Stage 5 - Lower-limb / site-specific pain extension
Generic harmoniser: read ANY standard GWAS summary-statistics file for a new
OUTCOME (knee pain, low back pain, ...), match it to the 1000G EUR reference
panel already used by Stage 4, resolve allele orientation, attach the LD score,
and write a harmonised file ready for the MR + overlap-audit pipeline.

Column auto-detection covers the formats produced by figshare / FinnGen /
IEU OpenGWAS / UK Biobank SAIGE/REGENIE exports:
  rsid : SNP, rsid, ID, rs_number, variant
  chr  : CHR, chr, #CHR, chromosome
  pos  : POS, pos, BP, position, bp
  a1   : A1, a1, EA, effect_allele, allele1   (reported EFFECT allele)
  a2   : A2, a2, NE, other_allele, reference_allele, allele2
  beta : BETA, beta, b, effect, log_odds
  se   : SE, se, standard_error, stderr
  p    : P, p, p_value, pval
  n    : N, n, sample_size, total_sample_size, OBS_CT
  af   : AF, af, A1_FREQ, FRQ, freq, effect_allele_freq, EAF
If only n_cases / n_controls are present, n = cases + controls.
If the file has no physical positions (rsid only), it is anchored on rsid like
the GIANT BMI loader.

Usage:
  python harmonise_generic.py <trait_name> <input_file> [--gz]
"""
import argparse
import gzip
import json
import os
import time

import numpy as np
import pandas as pd

LDSC_DIR = r"D:\SX\TMFI_CPI_MR\04_stage4_sampling_cov\stage4_ldsc"
OUT_DIR = r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain\gwas_harmonised"
os.makedirs(OUT_DIR, exist_ok=True)

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
PALINDROME = {frozenset(("A", "T")), frozenset(("G", "C"))}

# ---- column-name synonyms (lower-cased) ----
SYN = {
    "rsid": ["snp", "rsid", "id", "rs_number", "variant", "snpid"],
    "chr": ["chr", "#chr", "chrom", "chromosome"],
    "pos": ["pos", "bp", "position", "bp_hg19", "bp_hg38"],
    "a1": ["a1", "ea", "effect_allele", "allele1", "allele_1"],
    "a2": ["a2", "ne", "other_allele", "reference_allele", "allele2", "allele_2"],
    "beta": ["beta", "b", "effect", "log_odds", "effect_size"],
    "se": ["se", "standard_error", "stderr", "standarderror"],
    "p": ["p", "p_value", "pval", "p_value_gc"],
    "n": ["n", "sample_size", "total_sample_size", "obs_ct", "n_total"],
    "af": ["af", "a1_freq", "af1", "freq", "frq", "effect_allele_freq", "eaf",
           "a1freq", "alt_freq", "a1_freq_hrc"],
    "cases": ["n_cases", "cases", "ncase"],
    "controls": ["n_controls", "controls", "ncontrol"],
}


def _find(columns):
    cols = [c.lower() for c in columns]
    out = {}
    for role, names in SYN.items():
        for cand in names:
            for i, c in enumerate(cols):
                if c == cand or c.endswith("_" + cand) or c.startswith(cand + "_"):
                    out[role] = columns[i]
                    break
            if role in out:
                break
    return out


def load_reference():
    ref = pd.read_csv(os.path.join(LDSC_DIR, "eur_ldscores.tsv.gz"), sep="\t")
    ref = ref.rename(columns={"chr": "chrom", "pos": "position"})
    ref["key"] = ref["chrom"].astype(np.int64) * 1_000_000_000 + ref["position"].astype(np.int64)
    ref["is_snp"] = ref["a1"].str.len().eq(1) & ref["a2"].str.len().eq(1)
    return ref


def harmonise(df, ref_sub, trait_name):
    d = df.merge(ref_sub[["key", "a1", "a2", "maf", "af1", "ldscore", "is_snp"]],
                 on="key", how="left", suffixes=("", "_ref"))
    d = d[d["a1_ref"].notna()]
    d["is_snp"] = d["is_snp"].fillna(False)
    n0 = len(d)
    a1 = d["a1"].str.upper()
    a2 = d["a2"].str.upper()
    r1 = d["a1_ref"]
    r2 = d["a2_ref"]
    pal = np.array([frozenset((x, y)) in PALINDROME for x, y in zip(a1, a2)])
    snp_a = a1.str.len().eq(1) & a2.str.len().eq(1)
    direct1 = (a1 == r1) & (a2 == r2)
    swap1 = (a1 == r2) & (a2 == r1)
    cc1 = a1.map(COMP); cc2 = a2.map(COMP)
    direct2 = (cc1 == r1) & (cc2 == r2)
    swap2 = (cc1 == r2) & (cc2 == r1)
    keep = snp_a & (~pal) & (direct1 | swap1 | direct2 | swap2)
    keep = keep.values
    d = d[keep].copy()
    r1 = r1[keep]; r2 = r2[keep]
    flip = (np.asarray(swap1[keep]) | np.asarray(swap2[keep]))
    if "beta" in d.columns:
        d["beta"] = np.where(flip, -d["beta"], d["beta"])
    d["allele_ea"] = np.asarray(r1)
    d["allele_nea"] = np.asarray(r2)
    if "af" in d.columns:
        d["af"] = np.where(flip, 1.0 - d["af"], d["af"])
    print(f"  {trait_name}: matched={len(d):,}  (dropped {n0 - len(d):,} "
          f"unresolvable/ambiguous; {int(flip.sum()):,} allele swaps, "
          f"{int(np.asarray(swap2[keep]).sum()):,} strand-complement swaps)")
    return d


def stream_and_harmonise(path, trait_name, gz):
    ref = load_reference()
    keep_keys = set(ref["key"].to_numpy())
    ref_rsid = ref.loc[ref["rsid"].notna(), ["rsid", "key", "chrom", "position"]].copy()
    ref_rsid = ref_rsid.drop_duplicates(subset=["rsid"], keep=False)
    opener = gzip.open if gz else open
    chunks = []
    nrows = 0
    t0 = time.time()
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        print(f"    header ({len(header)} cols): {header[:14]}", flush=True)
        mapping = _find(header)
        print(f"    mapped: {mapping}", flush=True)
        has_pos = ("chr" in mapping) and ("pos" in mapping)
        has_rsid = "rsid" in mapping
        if not has_pos and not has_rsid:
            raise SystemExit("No chromosome/position or rsid column detected.")
        usecols = list(mapping.values())
        reader = pd.read_csv(fh, sep="\t", names=header, chunksize=2_000_000,
                             low_memory=False, comment="#", dtype=str,
                             usecols=lambda c: c in usecols)
        for ch in reader:
            nrows += len(ch)
            ch = ch.rename(columns={v: k for k, v in mapping.items()})
            if has_pos:
                ch["chrom"] = pd.to_numeric(ch["chr"], errors="coerce")
                ch["position"] = pd.to_numeric(ch["pos"], errors="coerce")
                ch = ch[ch["chrom"].between(1, 22) & ch["position"].notna()]
                ch["key"] = ch["chrom"].astype(np.int64) * 1_000_000_000 + \
                    ch["position"].astype(np.int64)
                ch = ch[ch["key"].isin(keep_keys)]
            else:  # rsid-only anchoring
                ch = ch.merge(ref_rsid, left_on="rsid", right_on="rsid", how="inner")
            if "cases" in ch.columns and "controls" in ch.columns:
                ch["n"] = pd.to_numeric(ch["cases"], errors="coerce") + \
                    pd.to_numeric(ch["controls"], errors="coerce")
            for c in ("beta", "se", "p", "n", "af"):
                if c in ch.columns:
                    ch[c] = pd.to_numeric(ch[c], errors="coerce")
            ch = ch.dropna(subset=[c for c in ("beta", "se", "p", "n") if c in ch.columns])
            if len(ch):
                chunks.append(ch)
            if nrows % 4_000_000 == 0:
                print(f"    scanned {nrows:,} raw rows, kept "
                      f"{sum(len(c) for c in chunks):,} ({time.time() - t0:.0f}s)",
                      flush=True)
    raw = pd.concat(chunks, ignore_index=True)
    out_cols = ["key", "chrom", "position", "rsid", "allele_ea", "allele_nea",
                "beta", "se", "z", "n", "af", "af1", "maf", "ldscore"]
    d = harmonise(raw, ref, trait_name)
    d["z"] = d["beta"] / d["se"]
    d = d[out_cols]
    out_path = os.path.join(OUT_DIR, f"{trait_name}_harmonised.tsv.gz")
    d.to_csv(out_path, sep="\t", index=False, compression="gzip", float_format="%.6g")
    info = {"file": path, "trait": trait_name, "raw_rows": int(nrows),
            "harmonised_rows": int(len(d)),
            "n_quantiles": {str(k): float(v) for k, v in
                            d["n"].quantile([0, .01, .05, .5, .95, .99, 1]).items()},
            "af_vs_reference": {
                "median_abs_dev": float((d["af"] - d["af1"]).abs().median()) if "af" in d else None,
                "correlation": float(np.corrcoef(d["af"], d["af1"])[0, 1]) if "af" in d else None,
            },
            "scan_seconds": round(time.time() - t0, 1)}
    json.dump(info, open(os.path.join(OUT_DIR, f"{trait_name}_raw_scan.json"), "w"),
              indent=2)
    print(f"[{trait_name}] wrote {len(d):,} rows -> {out_path}\n"
          f"{json.dumps(info, indent=2)}", flush=True)
    return info


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("trait")
    ap.add_argument("input")
    ap.add_argument("--gz", action="store_true")
    args = ap.parse_args()
    stream_and_harmonise(args.input, args.trait, args.gz)
