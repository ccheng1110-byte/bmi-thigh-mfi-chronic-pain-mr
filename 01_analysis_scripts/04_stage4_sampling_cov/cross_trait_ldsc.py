#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Stage 4 - Sampling covariance audit
Step 3: single-trait and cross-trait LD Score regression.

Estimating equations (Bulik-Sullivan et al. 2015, Nat Genet 47:291-295)
-----------------------------------------------------------------------
Single trait k            E[z_kj^2]   = 1 + (N_k h_k^2 / M) * l_j
Cross trait (k, m)        E[z_kj z_mj] = (sqrt(N_k N_m) rho_g / M) * l_j
                                        + N_shared * rho_phen / sqrt(N_k N_m)

so the *intercept* of the cross-trait regression is

    intercept_km = N_shared * rho_phen / sqrt(N_k N_m)

which is exactly the quantity that generates bias in a two-sample MR design
when the two GWAS share individuals.  Note that a multiplicative rescaling of
the LD score changes the slope but NOT the intercept, so the target estimand
is insensitive to how the LD scores were normalised.

Estimation: OLS -> weights -> weighted least squares (refit weights once),
with a delete-one-block jackknife (200 contiguous blocks) for standard errors.
"""
import itertools
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

GWAS_DIR = r"D:\SX\TMFI_CPI_MR\04_stage4_sampling_cov\stage4_gwas"
OUT_DIR = r"D:\SX\TMFI_CPI_MR\04_stage4_sampling_cov\stage4_results"
N_JACK = 200
MAF_MIN = 0.01
AF_TOL = 0.15

os.makedirs(OUT_DIR, exist_ok=True)


# --------------------------------------------------------------------------
# regression machinery
# --------------------------------------------------------------------------
def _blocks(n, n_blocks):
    edges = np.linspace(0, n, n_blocks + 1).astype(int)
    return [(edges[i], edges[i + 1]) for i in range(n_blocks)]


def wls_fit(X, y, w, block_slices):
    """Return (theta, A_inv, per-block A and b contributions)."""
    Xw = X * w[:, None]
    A = X.T @ Xw
    b = Xw.T @ y
    Ainv = np.linalg.inv(A)
    theta = Ainv @ b
    Ab = []
    for s, e in block_slices:
        Xb = X[s:e]
        Xbw = Xb * w[s:e, None]
        Ab.append((Xb.T @ Xbw, Xbw.T @ y[s:e]))
    return theta, Ainv, A, b, Ab


def jackknife_se(theta_full, A, b, Ab, solvefn):
    """Delete-one-block jackknife standard errors (solved products reused)."""
    est = []
    for Ai, bi in Ab:
        Ai2 = A - Ai
        bi2 = b - bi
        try:
            est.append(np.linalg.solve(Ai2, bi2))
        except np.linalg.LinAlgError:
            est.append(np.full(len(bi2), np.nan))
    est = np.array(est)
    est = est[~np.isnan(est).any(axis=1)]
    n = len(est)
    var = (n - 1) / n * ((est - est.mean(axis=0)) ** 2).sum(axis=0)
    return np.sqrt(var), est


def single_trait_regression(l, z, N, M, weights=None, blocks=None):
    X = np.column_stack([np.ones_like(l), N * l / M])
    y = z * z
    if weights is None:
        w = np.ones_like(l)
    else:
        w = weights
    theta, Ainv, A, b, Ab = wls_fit(X, y, w, blocks)
    se, jk = jackknife_se(theta, A, b, Ab, None)
    return {"intercept": theta[0], "intercept_se": se[0],
            "hsq": theta[1], "hsq_se": se[1],
            "_theta": theta, "_A": A, "_b": b, "_Ab": Ab}


def cross_trait_regression(l, z1, z2, N1, N2, M, weights=None, blocks=None):
    X = np.column_stack([np.ones_like(l), np.sqrt(N1 * N2) * l / M])
    y = z1 * z2
    w = np.ones_like(l) if weights is None else weights
    theta, Ainv, A, b, Ab = wls_fit(X, y, w, blocks)
    se, jk = jackknife_se(theta, A, b, Ab, None)
    return {"intercept": theta[0], "intercept_se": se[0],
            "gencov": theta[1], "gencov_se": se[1],
            "_theta": theta, "_A": A, "_b": b, "_Ab": Ab}


# --------------------------------------------------------------------------
# analysis set
# --------------------------------------------------------------------------
def load_trait(name):
    d = pd.read_csv(os.path.join(GWAS_DIR, f"{name}_harmonised.tsv.gz"), sep="\t")
    d = d[np.isfinite(d["beta"]) & np.isfinite(d["se"]) & (d["se"] > 0)]
    d = d[(d["n"] > 0) & np.isfinite(d["n"])]
    d["af_dev"] = (d["af"] - d["af1"]).abs()
    return d


def build_pair(traits, maf_min=MAF_MIN, af_tol=AF_TOL, chi2_cap=True):
    """Merge harmonised traits on the reference panel and apply LDSC filters."""
    if isinstance(traits, str):
        traits = [traits]
    merged = load_trait(traits[0])
    cols = ["key", "chrom", "position", "ldscore", "maf", "af1"]
    ren = {c: c for c in cols}
    for nm in traits:
        d = load_trait(nm)
        keep = cols + [f"{nm}_beta", f"{nm}_se", f"{nm}_z", f"{nm}_n", f"{nm}_af_dev"]
        d = d.rename(columns={"beta": f"{nm}_beta", "se": f"{nm}_se",
                              "z": f"{nm}_z", "n": f"{nm}_n",
                              "af_dev": f"{nm}_af_dev"})
        merged = merged.merge(d[keep], on=cols, how="inner")
    merged = merged.rename(columns=ren)

    m0 = len(merged)
    merged = merged[(merged["ldscore"] > 0) & (merged["maf"] >= maf_min)]
    m1 = len(merged)
    # allele-frequency concordance with the reference panel
    dev_cols = [f"{nm}_af_dev" for nm in traits]
    ok = np.ones(len(merged), bool)
    for c in dev_cols:
        ok &= (merged[c].fillna(0) <= af_tol).to_numpy()
    merged = merged[ok]
    m2 = len(merged)
    # chi-square cap
    if chi2_cap:
        for nm in traits:
            cap = max(80.0, 0.001 * float(merged[f"{nm}_n"].max()))
            merged = merged[(merged[f"{nm}_z"] ** 2) <= cap]
    m3 = len(merged)
    merged = merged.sort_values(["chrom", "position"]).reset_index(drop=True)
    info = {"n_merged": m0, "after_ldsc_maf": m1, "after_af_check": m2,
            "after_chi2_cap": m3}
    return merged, info


def main():
    ld = pd.read_csv(os.path.join(
        r"D:\SX\TMFI_CPI_MR\04_stage4_sampling_cov\stage4_ldsc",
        "eur_ldscores.tsv.gz"), sep="\t")
    M = int((ld["ldscore"] > 0).sum())
    print(f"[M] number of polymorphic reference SNPs used as denominator: {M:,}")

    pairs = [("BACK_PAIN_UKB", "KNEE_PAIN_UKB"),
             ("KNEE_PAIN_UKB", "KNEE_PAIN_FIGSHARE"),
             ("TMFI", "BACK_PAIN_UKB"),
             ("TMFI", "KNEE_PAIN_UKB"),
             ("TMFI", "KNEE_PAIN_FIGSHARE"),
             ("BMI", "BACK_PAIN_UKB"),
             ("BMI", "KNEE_PAIN_UKB")]
    results = {"M": M, "settings": {"jackknife_blocks": N_JACK, "maf_min": MAF_MIN,
                                    "af_tol": AF_TOL}, "pairs": {}}

    for t1, t2 in pairs:
        print(f"\n=== {t1} x {t2} ===", flush=True)
        data, finfo = build_pair([t1, t2])
        print(f"  SNPs used: {len(data):,}   {finfo}")
        l = data["ldscore"].to_numpy(float)
        N1 = data[f"{t1}_n"].to_numpy(float)
        N2 = data[f"{t2}_n"].to_numpy(float)
        z1 = data[f"{t1}_z"].to_numpy(float)
        z2 = data[f"{t2}_z"].to_numpy(float)
        blocks = _blocks(len(l), N_JACK)

        # ---- step 1: OLS single-trait for initial variance parameters ------
        s1 = single_trait_regression(l, z1, N1, M, None, blocks)
        s2 = single_trait_regression(l, z2, N2, M, None, blocks)
        c0 = cross_trait_regression(l, z1, z2, N1, N2, M, None, blocks)

        # ---- step 2: weights, then WLS (weights refreshed once) -----------
        def make_weights():
            y1 = (s1["intercept"] + (s1["hsq"] * N1 / M) * l)
            y2 = (s2["intercept"] + (s2["hsq"] * N2 / M) * l)
            wc = 1.0 / np.maximum(y1 * y2 + (c0["gencov"] * np.sqrt(N1 * N2) / M * l
                                            + c0["intercept"]) ** 2, 1e-9)
            wh1 = 1.0 / np.maximum(4 * y1 - 2, 1e-9)
            wh2 = 1.0 / np.maximum(4 * y2 - 2, 1e-9)
            return wh1, wh2, wc

        wh1, wh2, wc = make_weights()
        s1 = single_trait_regression(l, z1, N1, M, wh1, blocks)
        s2 = single_trait_regression(l, z2, N2, M, wh2, blocks)
        c0 = cross_trait_regression(l, z1, z2, N1, N2, M, wc, blocks)
        wh1, wh2, wc = make_weights()
        sf1 = single_trait_regression(l, z1, N1, M, wh1, blocks)
        sf2 = single_trait_regression(l, z2, N2, M, wh2, blocks)
        sf = cross_trait_regression(l, z1, z2, N1, N2, M, wc, blocks)

        # ---- genetic correlation -----------------------------------------
        rg = sf["gencov"] / np.sqrt(sf1["hsq"] * sf2["hsq"]) \
            if (sf1["hsq"] > 0 and sf2["hsq"] > 0) else np.nan

        def describe(name, s):
            z = s["intercept"] / s["intercept_se"]
            p = 2 * stats.norm.sf(abs(z))
            print(f"  [{name}] intercept={s['intercept']:.4f} "
                  f"(SE {s['intercept_se']:.4f}, p={p:.3g})   "
                  f"h2={s['hsq']:.4f} (SE {s['hsq_se']:.4f})")
            return {"intercept": s["intercept"], "intercept_se": s["intercept_se"],
                    "intercept_p": float(p), "hsq": s["hsq"], "hsq_se": s["hsq_se"]}

        o1 = describe(t1, sf1)
        o2 = describe(t2, sf2)
        ci = sf["intercept"]
        cse = sf["intercept_se"]
        cz = ci / cse
        cp = 2 * stats.norm.sf(abs(cz))
        Nmean1 = float(N1.mean())
        Nmean2 = float(N2.mean())
        n_shared_lb = ci * np.sqrt(Nmean1 * Nmean2)          # rho_phen <= 1
        print(f"  [cross] intercept={ci:.4f} (SE {cse:.4f}, p={cp:.3g})   "
              f"gencov={sf['gencov']:.5f} (SE {sf['gencov_se']:.5f})   rg={rg:.4f}")
        print(f"  N1={Nmean1:,.0f}  N2={Nmean2:,.0f}   "
              f"N_shared*rho_phen >= {n_shared_lb:,.0f}   "
              f"(rho_phen=1 => N_shared ~ {n_shared_lb:,.0f})")
        print(f"  mean chi2 {t1}={np.mean(z1**2):.3f}   {t2}={np.mean(z2**2):.3f}")

        results["pairs"][f"{t1}__{t2}"] = {
            "n_snps": int(len(data)),
            "filter_trace": {k: int(v) for k, v in finfo.items()},
            "mean_n": {t1: Nmean1, t2: Nmean2},
            "trait1": o1, "trait2": o2,
            "cross": {"intercept": float(ci), "intercept_se": float(cse),
                      "intercept_p": float(cp), "gencov": float(sf["gencov"]),
                      "gencov_se": float(sf["gencov_se"]),
                      "rg": float(rg)},
            "ols_cross_intercept": float(c0["intercept"]),
            "n_shared_times_rho": float(n_shared_lb),
            "n_shared_if_rho1": float(n_shared_lb),
            "mean_chi2": {t1: float(np.mean(z1 ** 2)), t2: float(np.mean(z2 ** 2))},
        }

    out_path = os.path.join(OUT_DIR, "stage5_cross_trait_ldsc.json")
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\n[out] {out_path}")


if __name__ == "__main__":
    main()
