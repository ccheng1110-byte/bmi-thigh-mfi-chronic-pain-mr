# -*- coding: utf-8 -*-
"""Stage 6 - Step 2: allele harmonisation + effect-scale audit for FinnGen R12.

Two things are settled here and nothing downstream may revisit them:

1. ALIGNMENT.  FinnGen reports beta on the ALT allele.  Our instrument file
   already fixes the effect allele (TMFI_A1, and BMI harmonised to the same
   orientation), so every variant is oriented to TMFI_A1.  Palindromic
   variants cannot be oriented from the alleles alone (ref/alt strand of a
   palindromic SNP is ambiguous), so they are resolved by allele frequency and
   dropped when the frequency evidence is not decisive.

2. SCALE.  The whole point of this replication is to escape the BOLT-LMM
   risk-difference trap hit in Stage 5.  FinnGen fits REGENIE logistic models,
   so beta should already be a log odds ratio and NO lambda rescaling may be
   applied.  We test that claim rather than assume it:

       logistic  (log OR) :  SE^2 * 2 N f (1-f)  =  1 / (mu (1-mu))
       linear    (risk d) :  SE^2 * 2 N f (1-f)  =      mu (1-mu)

   The two targets differ by 1/(mu(1-mu))^2 -- a factor of 16 for PAIN and
   10,000 for M13_FIBROMYALGIA -- so the verdict cannot be ambiguous.  Only
   non-significant variants are used, because FinnGen back-derives standard
   errors for some significant hits (Firth / SPA) and those break the identity.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import os
import statistics

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXTRACTED = os.path.join(HERE, "extracted")
HARM = os.path.join(HERE, "gwas_harmonised")
UNION = os.path.join(HERE, "instruments_master.csv")   # 109 rsids, both sets
MANIFEST = os.path.join(HERE, "finngen_R12_manifest.tsv")
PANEL = os.path.join(ROOT, "04_stage4_sampling_cov", "stage4_ldsc", "eur_ldscores.tsv.gz")

PHENOS = ["PAIN", "M13_FIBROMYALGIA", "M13_LIMBPAIN"]
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
AF_DECIDE = 0.15     # palindromic: |af_alt - af_exp| must beat this to be usable
AF_CONCORD = 0.20    # Finnish vs UKB drift tolerance for the QC flag


def manifest_n():
    """case/control counts straight from the R12 manifest (no transcription)."""
    out = {}
    with io.open(MANIFEST, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["phenocode"]] = (int(r["num_cases"]), int(r["num_controls"]))
    return out


def load_panel():
    p = {}
    with gzip.open(PANEL, "rt", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            p[r["rsid"]] = float(r["ldscore"])
    return p


def align(ref, alt, af_alt, a1, a2, af_exp):
    """Return (sign, af_aligned, rule).  sign None => unusable."""
    if ref not in COMP or alt not in COMP:
        return None, None, "non-snp"
    palindromic = COMP.get(a1) == a2
    if not palindromic:
        if (alt, ref) == (a1, a2):
            return 1.0, af_alt, "direct"
        if (alt, ref) == (a2, a1):
            return -1.0, 1.0 - af_alt, "flip"
        if (COMP[alt], COMP[ref]) == (a1, a2):
            return 1.0, af_alt, "complement-direct"
        if (COMP[alt], COMP[ref]) == (a2, a1):
            return -1.0, 1.0 - af_alt, "complement-flip"
        return None, None, "allele-mismatch"
    # palindromic: strand of ref/alt is unknowable, decide on allele frequency
    if alt == a1:
        d_same, d_flip = abs(af_alt - af_exp), abs((1.0 - af_alt) - af_exp)
    elif alt == a2:
        d_flip, d_same = abs(af_alt - af_exp), abs((1.0 - af_alt) - af_exp)
    else:
        return None, None, "allele-mismatch"
    if min(d_same, d_flip) > AF_DECIDE:
        return None, None, "palindromic-af-undecidable"
    if d_same < d_flip:
        return 1.0, af_alt, "palindromic-af-direct"
    return -1.0, 1.0 - af_alt, "palindromic-af-flip"


def main():
    os.makedirs(HARM, exist_ok=True)
    man = manifest_n()
    panel = load_panel()
    union = {}
    with io.open(UNION, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            union[r["rsid"]] = r

    report = {}
    for ph in PHENOS:
        raw_path = os.path.join(EXTRACTED, "%s_raw.tsv" % ph)
        if not os.path.isfile(raw_path):
            print("[skip] %s -- no extraction" % ph)
            continue
        ncase, nctrl = man[ph]
        n = ncase + nctrl
        mu = ncase / n
        n_eff = n * mu * (1 - mu)

        rows = list(csv.DictReader(io.open(raw_path, encoding="utf-8"), delimiter="\t"))
        out_rows, dropped, rules, af_dev = [], [], {}, []
        for r in rows:
            rs = r["rsid"]
            u = union[rs]
            a1, a2 = u["A1"], u["A2"]
            try:
                beta = float(r["beta"]); se = float(r["sebeta"])
                af_alt = float(r["af_alt"]); pval = float(r["pval"])
            except (TypeError, ValueError):
                dropped.append((rs, "missing-value")); continue
            af_exp = float(u["EAF"])
            sign, af, rule = align(r["ref"], r["alt"], af_alt, a1, a2, af_exp)
            rules[rule] = rules.get(rule, 0) + 1
            if sign is None:
                dropped.append((rs, rule)); continue
            af_dev.append(abs(af - af_exp))
            chrom, pos = int(u["CHR"]), int(u["POS"])     # b37, from our panel
            out_rows.append({
                "key": chrom * 1_000_000_000 + pos,
                "chrom": chrom, "position": pos, "rsid": rs,
                "allele_ea": a1, "allele_nea": a2,
                "beta": sign * beta, "se": se, "z": sign * beta / se,
                "n": n, "n_eff": round(n_eff, 1),
                "af": af, "af1": af, "maf": min(af, 1 - af),
                "ldscore": panel.get(rs, ""),
                "pval": pval, "rule": rule,
                "fg_chrom": r["chrom"], "fg_pos": r["pos"],
                "af_dev": abs(af - af_exp),
            })

        # ---------------- scale identity ----------------
        use = [r for r in out_rows if r["pval"] > 0.05 and 0.0 < r["maf"] < 0.5]
        vals = [r["se"] ** 2 * 2.0 * n * r["maf"] * (1 - r["maf"]) for r in use]
        med = statistics.median(vals) if vals else float("nan")
        tgt_logit = 1.0 / (mu * (1 - mu))
        tgt_linear = mu * (1 - mu)
        verdict = ("log-OR (logistic)" if abs(med - tgt_logit) < abs(med - tgt_linear)
                   else "RISK DIFFERENCE (linear)")
        ratio_logit = med / tgt_logit

        cols = ["key", "chrom", "position", "rsid", "allele_ea", "allele_nea",
                "beta", "se", "z", "n", "af", "af1", "maf", "ldscore"]
        out_tsv = os.path.join(HARM, "%s_harmonised.tsv.gz" % ph)
        with gzip.open(out_tsv, "wt", encoding="utf-8", newline="") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in out_rows:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        # full audit copy (extra columns) for inspection
        with io.open(os.path.join(HARM, "%s_harmonised_audit.tsv" % ph), "w",
                     encoding="utf-8", newline="") as fh:
            acols = cols + ["pval", "rule", "fg_chrom", "fg_pos", "af_dev", "n_eff"]
            fh.write("\t".join(acols) + "\n")
            for r in out_rows:
                fh.write("\t".join(str(r[c]) for c in acols) + "\n")

        print("\n=== %s" % ph)
        print("    N = %d  (%d cases / %d controls)  mu = %.5f  N_eff = %.0f"
              % (n, ncase, nctrl, mu, n_eff))
        print("    extracted %d -> harmonised %d, dropped %d"
              % (len(rows), len(out_rows), len(dropped)))
        print("    alignment rules: %s" % rules)
        if dropped:
            print("    dropped: %s" % ", ".join("%s(%s)" % d for d in dropped[:12]))
        print("    AF |FinnGen - UKB| : median %.4f  max %.4f  >%.2f: %d"
              % (statistics.median(af_dev), max(af_dev), AF_CONCORD,
                 sum(1 for d in af_dev if d > AF_CONCORD)))
        print("    scale identity on %d non-significant variants:" % len(vals))
        print("        observed median SE^2*2Nf(1-f) = %.4f" % med)
        print("        logistic target  1/(mu(1-mu)) = %.4f   ratio %.3f"
              % (tgt_logit, ratio_logit))
        print("        linear  target      mu(1-mu)  = %.5f  ratio %.3f"
              % (tgt_linear, med / tgt_linear))
        print("        VERDICT: %s" % verdict)
        print("    wrote %s" % out_tsv)

        report[ph] = {
            "n": n, "ncase": ncase, "ncontrol": nctrl, "mu": mu, "n_eff": n_eff,
            "n_extracted": len(rows), "n_harmonised": len(out_rows),
            "dropped": dropped, "rules": rules,
            "af_dev_median": statistics.median(af_dev), "af_dev_max": max(af_dev),
            "scale_identity_median": med,
            "target_logit": tgt_logit, "target_linear": tgt_linear,
            "ratio_to_logit": ratio_logit, "verdict": verdict,
            "n_variants_tested": len(vals),
        }

    with io.open(os.path.join(HARM, "harmonise_report.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)
    print("\nreport -> %s" % os.path.join(HARM, "harmonise_report.json"))


if __name__ == "__main__":
    main()
