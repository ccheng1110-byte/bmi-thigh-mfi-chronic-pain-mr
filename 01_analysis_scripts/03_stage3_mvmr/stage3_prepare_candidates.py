import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"D:\SX\TMFI_CPI_MR")
RAW = ROOT / "00_stage0_audit" / "incoming"
STAGE1 = ROOT / "01_stage1_instruments"
OUT = ROOT / "03_stage3_mvmr"
BMI_PATH = OUT / "SNP_gwas_mc_merge_nogc.tbl.uniq.gz"
P_THRESHOLD = 5e-8
CHUNK = 500_000


def complement(a):
    return str(a).upper().translate(str.maketrans("ACGT", "TGCA"))


def pair_compatible(a1, a2, b1, b2):
    a = {str(a1).upper(), str(a2).upper()}
    b = {str(b1).upper(), str(b2).upper()}
    return a == b or a == {complement(b1), complement(b2)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    bmi = pd.read_csv(BMI_PATH, sep="\t", compression="gzip", dtype={"SNP": "string"})
    bmi["F_BMI"] = (bmi["b"] / bmi["se"]) ** 2
    bmi["RSID"] = bmi["SNP"].astype(str).str.match(r"^rs\d+$")
    bmi = bmi.sort_values(["p", "N"], ascending=[True, False], kind="mergesort").drop_duplicates("SNP")
    bmi_gws = bmi.loc[(bmi["p"] < P_THRESHOLD) & (bmi["F_BMI"] > 10) & bmi["RSID"]].copy()
    bmi_gws.to_csv(OUT / "bmi_gws_candidates.csv", index=False)

    tmfi_clean = pd.read_csv(STAGE1 / "tmfi_gws_exact_clean.csv", dtype={"SNP": "string"})
    candidate_ids = set(tmfi_clean["SNP"].astype(str)) | set(bmi_gws["SNP"].astype(str))

    tmfi_hits = []
    for chunk in pd.read_csv(RAW / "TMFI_GWAS.csv", chunksize=CHUNK, dtype={"SNP": "string"}):
        hit = chunk.loc[chunk["SNP"].astype(str).isin(candidate_ids)].copy()
        if not hit.empty:
            tmfi_hits.append(hit)
    tmfi_lookup_raw = pd.concat(tmfi_hits, ignore_index=True)

    bmi_lookup = bmi.loc[bmi["SNP"].astype(str).isin(candidate_ids)].copy()
    bmi_by_id = bmi_lookup.set_index("SNP", drop=False)
    rows = []
    for snp, group in tmfi_lookup_raw.groupby("SNP", sort=False):
        if snp not in bmi_by_id.index:
            continue
        brow = bmi_by_id.loc[snp]
        if isinstance(brow, pd.DataFrame):
            brow = brow.iloc[0]
        compatible = group.loc[
            group.apply(lambda r: pair_compatible(r["A1"], r["A2"], brow["A1"], brow["A2"]), axis=1)
        ].copy()
        if compatible.empty:
            continue
        compatible = compatible.sort_values(["P", "N"], ascending=[True, False], kind="mergesort")
        trow = compatible.iloc[0]
        rows.append(
            {
                "SNP": snp,
                "CHR": int(trow["CHR"]),
                "POS": int(trow["POS"]),
                "TMFI_A1": str(trow["A1"]).upper(),
                "TMFI_A2": str(trow["A2"]).upper(),
                "TMFI_EAF": float(trow["AF1"]),
                "TMFI_BETA": float(trow["BETA"]),
                "TMFI_SE": float(trow["SE"]),
                "TMFI_P": float(trow["P"]),
                "TMFI_N": int(trow["N"]),
                "BMI_A1": str(brow["A1"]).upper(),
                "BMI_A2": str(brow["A2"]).upper(),
                "BMI_EAF": float(brow["Freq1.Hapmap"]),
                "BMI_BETA": float(brow["b"]),
                "BMI_SE": float(brow["se"]),
                "BMI_P": float(brow["p"]),
                "BMI_N": int(brow["N"]),
            }
        )
    cross = pd.DataFrame(rows)
    cross["TMFI_F"] = (cross["TMFI_BETA"] / cross["TMFI_SE"]) ** 2
    cross["BMI_F"] = (cross["BMI_BETA"] / cross["BMI_SE"]) ** 2
    cross["TMFI_GWS"] = cross["TMFI_P"] < P_THRESHOLD
    cross["BMI_GWS"] = cross["BMI_P"] < P_THRESHOLD
    cross = cross.loc[cross["TMFI_GWS"] | cross["BMI_GWS"]].copy()
    cross["P_SELECT"] = cross[["TMFI_P", "BMI_P"]].min(axis=1)

    panel_ids = set(
        pd.read_csv(
            ROOT / "reference" / "1000G_EUR.bim",
            sep=r"\s+",
            header=None,
            usecols=[1],
            names=["SNP"],
            dtype=str,
        )["SNP"]
    )
    cross["IN_LD_PANEL"] = cross["SNP"].astype(str).isin(panel_ids)
    cross.to_csv(OUT / "mvmr_candidate_exposure_associations.csv", index=False)
    clump_input = cross.loc[cross["IN_LD_PANEL"], ["SNP", "P_SELECT"]].rename(columns={"P_SELECT": "P"})
    clump_input.to_csv(OUT / "mvmr_clump_input.tsv", sep="\t", index=False)

    summary = {
        "bmi_rows": int(len(bmi)),
        "bmi_gws_strong_rsid": int(len(bmi_gws)),
        "tmfi_gws_strong_rsid": int(len(tmfi_clean)),
        "candidate_union_ids": int(len(candidate_ids)),
        "cross_exposure_allele_compatible": int(len(cross)),
        "candidate_tmfi_gws": int(cross["TMFI_GWS"].sum()),
        "candidate_bmi_gws": int(cross["BMI_GWS"].sum()),
        "candidate_gws_both": int((cross["TMFI_GWS"] & cross["BMI_GWS"]).sum()),
        "candidate_in_ld_panel": int(cross["IN_LD_PANEL"].sum()),
    }
    (OUT / "stage3_prepare_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
