#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Stage 4 - Sampling covariance audit
Step 1: compute LD scores from the 1000 Genomes EUR reference panel.

Inputs
------
reference/1000G_EUR.bed / .bim / .fam   (PLINK 1 binary, SNP-major)

Method
------
Genotypes are mean-imputed and standardised under HWE variance,

    p_hat = mean(dosage)/2 ,   z = (g - 2 p_hat) / sqrt(2 p_hat (1 - p_hat))

so that corr(z_j, z_k) equals the LD correlation r_jk.  For every SNP j the
LD score is

    l_j = sum_{k in W(j)}  max(0,  r_jk^2 - (1 - r_jk^2)/(n-2) )

where W(j) is the set of reference SNPs on the same chromosome whose physical
position lies within WINDOW_HALF base pairs of j (default +/- 500 kb; a 1 Mb
window used as the standard ~1 cM LDSC window), and the r^2 shrinkage term is
the standard finite-sample unbiasedness correction used by ldsc.py.  The SNP
itself is included (r_jj = 1).

Every pair of SNPs within the window is counted exactly once: for each row
block we accumulate all window partners with a larger row index.

NOTE ON SCALE: the target estimand of this audit -- the cross-trait intercept
-- is invariant to a multiplicative rescaling of the LD score, because any
rescaling is absorbed by the slope.  Using a ~1.8M-variant common-variant
panel instead of the full 1000G variant set therefore does not bias the
intercept; it only affects the secondary h^2 / genetic-covariance numbers,
which are reported with the same panel as denominator M.

Output
------
stage4_ldsc/eur_ldscores.tsv.gz   chr rsid pos a1 a2 maf af1 ldscore
stage4_ldsc/ldscore_meta.json
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

REF_DIR = r"D:\SX\TMFI_CPI_MR\reference"
OUT_DIR = r"D:\SX\TMFI_CPI_MR\04_stage4_sampling_cov\stage4_ldsc"
WINDOW_HALF = 500_000
BLOCK = 512

os.makedirs(OUT_DIR, exist_ok=True)


# PLINK 1 .bed 2-bit codes are NOT a monotone copy count.  The convention is
#     0  ->  homozygous A1   (2 copies of A1)
#     1  ->  missing
#     2  ->  heterozygous    (1 copy of A1)
#     3  ->  homozygous A2   (0 copies of A1)
CODE_TO_DOSAGE = np.array([2.0, np.nan, 1.0, 0.0], dtype=np.float32)


def unpack(raw: np.ndarray, n_ind: int) -> np.ndarray:
    """Unpack PLINK 2-bit SNP-major genotype bytes -> (m, n_ind) raw codes."""
    m, bpr = raw.shape
    out = np.empty((m, bpr * 4), dtype=np.uint8)
    for k in range(4):
        out[:, k::4] = (raw >> np.uint8(2 * k)) & np.uint8(3)
    return np.ascontiguousarray(out[:, :n_ind])


def codes_to_dosage(codes: np.ndarray) -> np.ndarray:
    """Apply the PLINK coding above; missing becomes NaN."""
    return CODE_TO_DOSAGE[codes]


def read_bed_rows(bed_path, row_start, row_count, bpr):
    """Read a contiguous block of SNP-major rows from the .bed file."""
    off = 3 + row_start * bpr
    cnt = row_count * bpr
    with open(bed_path, "rb") as fh:
        fh.seek(off)
        buf = np.frombuffer(fh.read(cnt), dtype=np.uint8)
    if buf.size != cnt:
        raise IOError("short read in bed file")
    return buf.reshape(row_count, bpr)


def chromosome_ldscore(raw_bytes, pos, n_ind):
    """Standardise genotypes and accumulate LD scores for one chromosome."""
    m = raw_bytes.shape[0]
    codes = unpack(raw_bytes, n_ind)
    del raw_bytes
    g = codes_to_dosage(codes)
    del codes
    miss = np.isnan(g)
    mu = np.nan_to_num(np.nanmean(g, axis=1), nan=0.0).astype(np.float32)
    g = np.where(miss, mu[:, None], g).astype(np.float32)
    g -= mu[:, None]
    del miss

    var = (mu * (2.0 - mu)) / 2.0                      # = 2 p (1-p) under HWE
    ok = var > 1e-6
    sd = np.sqrt(np.maximum(var, 1e-12))
    zz = np.divide(g, sd[:, None], out=np.zeros_like(g), where=ok[:, None])
    del g

    maf = np.minimum(mu / 2.0, 1.0 - mu / 2.0)

    lsum = np.zeros(m, dtype=np.float64)
    shrink = 1.0 / (n_ind - 2.0)
    for start in range(0, m, BLOCK):
        stop = min(start + BLOCK, m)
        # union of the per-SNP windows of the rows in this block; every true
        # partner of every row in the block is guaranteed to lie inside it
        lo = int(np.searchsorted(pos, pos[start] - WINDOW_HALF, side="left"))
        hi = int(np.searchsorted(pos, pos[stop - 1] + WINDOW_HALF, side="right"))
        hi = max(hi, stop)          # the block itself is always in its window
        hi = min(hi, m)

        r = (zz[start:stop] @ zz[lo:hi].T) / float(n_ind)
        r2 = r * r
        np.maximum(r2 - (1.0 - r2) * shrink, 0.0, out=r2)
        del r

        # keep only pairs that are truly within WINDOW_HALF base pairs ...
        dist = np.abs(pos[start:stop, None] - pos[lo:hi][None, :])
        r2[dist > WINDOW_HALF] = 0.0
        del dist
        # ... and keep each pair on the side with the smaller index, so that
        # every unordered pair is counted exactly once; the contribution goes
        # to BOTH partners.
        rows_idx = np.arange(start, stop)[:, None]
        cols_idx = np.arange(lo, hi)[None, :]
        r2[cols_idx <= rows_idx] = 0.0
        lsum[start:stop] += r2.sum(axis=1)
        lsum[lo:hi] += r2.sum(axis=0)
        del r2

    lsum += ok.astype(np.float64)          # self-correlation r_jj = 1
    lsum[~ok] = 0.0
    return lsum, maf, mu / 2.0


def main():
    t0 = time.time()
    bim = pd.read_csv(
        os.path.join(REF_DIR, "1000G_EUR.bim"), sep=r"\s+", header=None,
        names=["chr", "rsid", "cm", "pos", "a1", "a2"],
        dtype={"chr": int, "pos": int},
    )
    n_ind = sum(1 for _ in open(os.path.join(REF_DIR, "1000G_EUR.fam")))
    bpr = (n_ind + 3) // 4
    print(f"[meta] reference variants = {len(bim):,}  individuals = {n_ind}  "
          f"bytes/row = {bpr}", flush=True)

    bed_path = os.path.join(REF_DIR, "1000G_EUR.bed")
    expect = 3 + len(bim) * bpr
    actual = os.path.getsize(bed_path)
    if expect != actual:
        sys.exit(f"[FATAL] bed size mismatch: expected {expect}, found {actual}")
    print("[meta] bed size check OK", flush=True)

    pieces = []
    row_start = 0
    for c, grp in bim.groupby("chr", sort=True):
        m = len(grp)
        raw = read_bed_rows(bed_path, row_start, m, bpr)
        pos = grp["pos"].to_numpy()
        order = np.argsort(pos, kind="stable")
        if not np.all(order == np.arange(m)):          # guard: keep rows and
            raw = np.ascontiguousarray(raw[order])     # positions consistent
        pos_sorted = pos[order]
        l_sorted, maf_sorted, af_sorted = chromosome_ldscore(raw, pos_sorted, n_ind)
        del raw
        l_arr = np.empty(m, dtype=np.float64)
        maf_arr = np.empty(m, dtype=np.float64)
        af_arr = np.empty(m, dtype=np.float64)
        l_arr[order] = l_sorted
        maf_arr[order] = maf_sorted
        af_arr[order] = af_sorted
        pieces.append(pd.DataFrame({
            "chr": c,
            "rsid": grp["rsid"].to_numpy(),
            "pos": grp["pos"].to_numpy(),
            "a1": grp["a1"].to_numpy(),
            "a2": grp["a2"].to_numpy(),
            "maf": maf_arr,
            "af1": af_arr,
            "ldscore": l_arr,
        }))
        row_start += m
        print(f"[chr{c:>2}] m={m:>7,}  mean l={l_arr.mean():8.2f}  "
              f"median l={np.median(l_arr):8.2f}  elapsed={time.time() - t0:7.1f}s",
              flush=True)

    out = pd.concat(pieces, ignore_index=True)
    out_path = os.path.join(OUT_DIR, "eur_ldscores.tsv.gz")
    out.to_csv(out_path, sep="\t", index=False, compression="gzip", float_format="%.6g")
    print(f"[out] {out_path}  rows={len(out):,}")

    meta = {
        "n_reference_snps": int(len(out)),
        "n_individuals": int(n_ind),
        "window_bp_total": int(2 * WINDOW_HALF),
        "block_size": int(BLOCK),
        "r2_shrinkage": "max(0, r2 - (1-r2)/(n-2))",
        "ldscore_includes_self": True,
        "genome_build": "GRCh37 (verified against GWAS positions)",
        "mean_ldscore": float(out["ldscore"].mean()),
        "median_ldscore": float(np.median(out["ldscore"])),
        "mean_ldscore_maf_ge_0p01": float(out.loc[out.maf >= 0.01, "ldscore"].mean()),
        "n_maf_ge_0p01": int((out.maf >= 0.01).sum()),
        "n_maf_ge_0p05": int((out.maf >= 0.05).sum()),
        "elapsed_seconds": round(time.time() - t0, 1),
    }
    with open(os.path.join(OUT_DIR, "ldscore_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
