# -*- coding: utf-8 -*-
"""Stage 6 - Step 0: one master instrument table covering BOTH frozen sets.

The two frozen instrument sets are not nested:

    42-SNP set  (01_stage1_instruments, TMFI genome-wide significant only)
    93-SNP set  (03_stage3_mvmr union, TMFI + BMI)   <-- only 26 of the 42

Only 26 of the 42 TMFI-significant SNPs survive into the 93-SNP union, so an
extraction keyed on the union alone would silently reduce the Stage 1
univariable analysis from 42 instruments to 26.  We therefore extract the
union of both sets (109 rsids) and keep the 93-SNP union's orientation
wherever a variant appears in both.
"""
import csv
import io
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
I1 = os.path.join(ROOT, "01_stage1_instruments",
                  "tmfi_gws_ldclumped_r2_0p001_kb10000.csv")
I2 = os.path.join(ROOT, "03_stage3_mvmr", "mvmr_union_exposures_usable.csv")
OUT = os.path.join(HERE, "instruments_master.csv")


def main():
    master = {}
    with io.open(I2, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            master[r["SNP"]] = {
                "rsid": r["SNP"], "CHR": r["CHR"], "POS": r["POS"],
                "A1": r["TMFI_A1"], "A2": r["TMFI_A2"], "EAF": r["TMFI_EAF"],
                "source": "union93",
            }
    n_from_union = len(master)
    added = 0
    with io.open(I1, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["SNP"] in master:
                continue
            master[r["SNP"]] = {
                "rsid": r["SNP"], "CHR": r["CHR"], "POS": r["POS"],
                "A1": r["A1"], "A2": r["A2"], "EAF": r["AF1"],
                "source": "tmfi42",
            }
            added += 1
    cols = ["rsid", "CHR", "POS", "A1", "A2", "EAF", "source"]
    with io.open(OUT, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for rs in sorted(master, key=lambda s: (int(master[s]["CHR"]),
                                                int(master[s]["POS"]))):
            w.writerow(master[rs])
    print("93-SNP union: %d" % n_from_union)
    print("added from the 42-SNP set: %d" % added)
    print("master instrument table: %d rsids -> %s" % (len(master), OUT))
    return master


if __name__ == "__main__":
    main()
