#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Convert an IEU OpenGWAS harmonised VCF (ukb-b-*.vcf.gz) into the flat TSV
format consumed by harmonise_generic.py.

OpenGWAS VCF conventions:
  - non-effect allele is aligned to the human reference (B37)  -> A2 = REF
  - effect allele is the non-reference allele                  -> A1 = ALT
  - the single sample column carries FORMAT keys ES:SE:LP:AF:ID
  - per-study n comes from the ##SAMPLE meta line
     (TotalCases + TotalControls)

Output columns: CHR POS SNP A1 A2 N AF1 BETA SE P
"""
import gzip
import os
import re
import sys

D = r"D:\SX\TMFI_CPI_MR\05_stage5_limbpain\incoming"
OUT = r"D:\SX/TMFI_CPI_MR/05_stage5_limbpain/standardised"
os.makedirs(OUT, exist_ok=True)

FILES = {
    "BACK_PAIN_UKB": "ukb-b-9838.vcf.gz",
    "KNEE_PAIN_UKB": "ukb-b-16254.vcf.gz",
}


def meta_n(line):
    cases = int(re.search(r"TotalCases=(\d+)", line).group(1))
    controls = int(re.search(r"TotalControls=(\d+)", line).group(1))
    return cases, controls


def convert(trait, fname):
    src = os.path.join(D, fname)
    dst = os.path.join(OUT, f"{trait}_standard.tsv.gz")
    n_total = None
    cases = controls = None
    fmt_idx = None
    kept = scanned = 0
    with gzip.open(src, "rt") as fh, gzip.open(dst, "wt", compresslevel=6) as out:
        out.write("CHR\tPOS\tSNP\tA1\tA2\tN\tAF1\tBETA\tSE\tP\n")
        for line in fh:
            if line.startswith("##"):
                if line.startswith("##SAMPLE"):
                    cases, controls = meta_n(line)
                    n_total = cases + controls
                continue
            if line.startswith("#CHROM"):
                hdr = line.lstrip("#").rstrip("\n").split("\t")
                fmt_idx = {k: i for i, k in enumerate(hdr)}
                continue
            f = line.rstrip("\n").split("\t")
            chrom = f[fmt_idx["CHROM"]]
            if not chrom.isdigit() or not (1 <= int(chrom) <= 22):
                continue
            scanned += 1
            keys = f[fmt_idx["FORMAT"]].split(":")
            vals = f[-1].split(":")          # single sample column is last
            kv = dict(zip(keys, vals))
            es = kv.get("ES")
            se = kv.get("SE")
            lp = kv.get("LP")
            af = kv.get("AF")
            if es in (".", "NA", "") or se in (".", "NA", "") or lp in (".", "NA", ""):
                continue
            try:
                es = float(es)
                se = float(se)
                lp = float(lp)
            except ValueError:
                continue
            if se <= 0:
                continue
            # p from -log10 p; guard underflow
            try:
                p = 10.0 ** (-lp)
            except OverflowError:
                p = 0.0
            out.write(f"{chrom}\t{f[fmt_idx['POS']]}\t{f[fmt_idx['ID']]}\t"
                      f"{f[fmt_idx['ALT']]}\t{f[fmt_idx['REF']]}\t{n_total}\t"
                      f"{af if af not in ('.', 'NA', '') else 'NA'}\t"
                      f"{es:.6g}\t{se:.6g}\t{p:.6g}\n")
            kept += 1
    print(f"[{trait}] cases={cases:,} controls={controls:,} n={n_total:,}")
    print(f"[{trait}] scanned autosomal rows {scanned:,}, wrote {kept:,} -> {dst}")
    return dst


if __name__ == "__main__":
    only = sys.argv[1:] or list(FILES)
    for t in only:
        convert(t, FILES[t])
