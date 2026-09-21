import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"D:\SX\TMFI_CPI_MR")
RAW = ROOT / "00_stage0_audit" / "incoming"
OUT = ROOT / "01_stage1_instruments"
CLUMPED = OUT / "plink_tmfi_r2_0p001_kb10000.clumped"
CHUNK = 500_000


def read_clumped_snps():
    clumped = pd.read_csv(CLUMPED, sep=r"\s+")
    return clumped["SNP"].astype(str).tolist()


def match_and_harmonise(candidates, name, filename, mapping):
    wanted = set(candidates["SNP"].astype(str))
    pieces = []
    for chunk in pd.read_csv(
        RAW / filename,
        sep="\t",
        compression="gzip",
        chunksize=CHUNK,
        dtype={mapping["snp"]: "string"},
    ):
        hit = chunk.loc[chunk[mapping["snp"]].astype(str).isin(wanted)].copy()
        if not hit.empty:
            pieces.append(hit)
    outcome = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()
    outcome.to_csv(OUT / f"{name}_ldclumped_matches_raw.csv", index=False)

    if outcome.empty:
        return {"candidate_rows": len(candidates), "matched_rows": 0, "harmonised": 0}

    renamed = outcome.rename(
        columns={
            mapping["chr"]: "CHR_OUT",
            mapping["pos"]: "POS_OUT",
            mapping["snp"]: "SNP",
            mapping["ea"]: "EA_OUT",
            mapping["oa"]: "OA_OUT",
            mapping["freq"]: "EAF_OUT",
            mapping["beta"]: "BETA_OUT",
            mapping["se"]: "SE_OUT",
            mapping["p"]: "P_OUT",
            mapping["n"]: "N_OUT",
        }
    )
    renamed = renamed[
        ["CHR_OUT", "POS_OUT", "SNP", "EA_OUT", "OA_OUT", "EAF_OUT", "BETA_OUT", "SE_OUT", "P_OUT", "N_OUT"]
    ]
    merged = candidates.merge(renamed, on="SNP", how="left", validate="one_to_one")

    statuses, beta_h, eaf_h = [], [], []
    for row in merged.itertuples(index=False):
        if pd.isna(row.BETA_OUT):
            statuses.append("missing_outcome")
            beta_h.append(np.nan)
            eaf_h.append(np.nan)
            continue
        ea, oa = str(row.A1).upper(), str(row.A2).upper()
        eao, oao = str(row.EA_OUT).upper(), str(row.OA_OUT).upper()
        pal = {ea, oa} in ({"A", "T"}, {"C", "G"}) and len(ea) == len(oa) == 1
        if pal and float(row.MAF) >= 0.42:
            statuses.append("drop_palindromic_intermediate")
            beta_h.append(np.nan)
            eaf_h.append(np.nan)
        elif ea == eao and oa == oao:
            statuses.append("aligned")
            beta_h.append(float(row.BETA_OUT))
            eaf_h.append(float(row.EAF_OUT))
        elif ea == oao and oa == eao:
            statuses.append("swapped")
            beta_h.append(-float(row.BETA_OUT))
            eaf_h.append(1 - float(row.EAF_OUT))
        else:
            statuses.append("allele_mismatch")
            beta_h.append(np.nan)
            eaf_h.append(np.nan)

    merged["HARMONISE_STATUS"] = statuses
    merged["BETA_OUT_HARMONISED"] = beta_h
    merged["EAF_OUT_HARMONISED"] = eaf_h
    merged.to_csv(OUT / f"{name}_ldclumped_harmonised.csv", index=False)
    usable = merged["HARMONISE_STATUS"].isin(["aligned", "swapped"])
    return {
        "candidate_rows": int(len(candidates)),
        "matched_rows": int(len(outcome)),
        "harmonised": int(usable.sum()),
        "status_counts": {str(k): int(v) for k, v in merged["HARMONISE_STATUS"].value_counts().items()},
    }


def main():
    clean = pd.read_csv(OUT / "tmfi_gws_exact_clean.csv")
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
    clean["IN_LD_PANEL"] = clean["SNP"].astype(str).isin(panel_ids)
    index_snps = read_clumped_snps()
    candidates = clean.set_index("SNP").loc[index_snps].reset_index()
    candidates.to_csv(OUT / "tmfi_gws_ldclumped_r2_0p001_kb10000.csv", index=False)

    full = match_and_harmonise(
        candidates,
        "cpi_full",
        "Chronic_pain_intensity.txt.gz",
        {"chr": "CHR", "pos": "BP", "snp": "SNP", "ea": "A1", "oa": "A2", "freq": "MAF", "beta": "BETA", "se": "SE", "p": "P", "n": "N"},
    )
    noimg = match_and_harmonise(
        candidates,
        "cpi_no_imaging",
        "chronic_pain_intensity_no_imaging.txt.gz",
        {"chr": "#CHROM", "pos": "POS", "snp": "ID", "ea": "A1", "oa": "OMITTED", "freq": "A1_FREQ", "beta": "BETA", "se": "SE", "p": "P", "n": "OBS_CT"},
    )
    summary = {
        "parameters": {"p_threshold": 5e-8, "clump_r2": 0.001, "clump_kb": 10000, "reference": "1000G_EUR_GRCh37_503"},
        "panel_overlap": {
            "eligible_strong_rsid": int(len(clean)),
            "present_in_panel": int(clean["IN_LD_PANEL"].sum()),
            "missing_from_panel": int((~clean["IN_LD_PANEL"]).sum()),
            "overlap_fraction": float(clean["IN_LD_PANEL"].mean()),
        },
        "instruments": {
            "ld_clumped": int(len(candidates)),
            "f_min": float(candidates["F"].min()),
            "f_median": float(candidates["F"].median()),
            "f_max": float(candidates["F"].max()),
            "palindromic": int(candidates["PALINDROMIC"].sum()),
        },
        "cpi_full": full,
        "cpi_no_imaging": noimg,
    }
    (OUT / "stage1_formal_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
