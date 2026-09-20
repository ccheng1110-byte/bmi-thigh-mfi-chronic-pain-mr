import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"D:\SX\TMFI_CPI_MR")
RAW = ROOT / "00_stage0_audit" / "incoming"
OUT = ROOT / "01_stage1_instruments"
OUT.mkdir(parents=True, exist_ok=True)

P_THRESHOLD = 5e-8
F_THRESHOLD = 10.0
DISTANCE_BP = 10_000_000
CHUNK = 500_000


def allele_pair_key(row):
    alleles = sorted([str(row.A1).upper(), str(row.A2).upper()])
    return f"{int(row.CHR)}:{int(row.POS)}:{alleles[0]}:{alleles[1]}"


def distance_prune(df):
    selected = []
    selected_positions = {chrom: [] for chrom in range(1, 23)}
    for row in df.sort_values(["P", "CHR", "POS", "SNP"], kind="mergesort").itertuples(index=False):
        chrom = int(row.CHR)
        pos = int(row.POS)
        if all(abs(pos - old) >= DISTANCE_BP for old in selected_positions[chrom]):
            selected.append(row)
            selected_positions[chrom].append(pos)
    return pd.DataFrame(selected, columns=df.columns)


def extract_tmfi():
    path = RAW / "TMFI_GWAS.csv"
    hits = []
    for chunk in pd.read_csv(path, chunksize=CHUNK):
        hit = chunk.loc[chunk["P"] < P_THRESHOLD].copy()
        if not hit.empty:
            hits.append(hit)
    raw = pd.concat(hits, ignore_index=True)
    raw["F"] = (raw["BETA"] / raw["SE"]) ** 2
    raw["MAF"] = np.minimum(raw["AF1"], 1 - raw["AF1"])
    raw["PALINDROMIC"] = ((raw["A1"] == "A") & (raw["A2"] == "T")) | ((raw["A1"] == "T") & (raw["A2"] == "A")) | ((raw["A1"] == "C") & (raw["A2"] == "G")) | ((raw["A1"] == "G") & (raw["A2"] == "C"))
    raw["VARIANT_KEY"] = raw.apply(allele_pair_key, axis=1)
    raw["RSID"] = raw["SNP"].astype(str).str.match(r"^rs\d+$")
    raw.to_csv(OUT / "tmfi_gws_raw.csv", index=False)

    key_counts = raw["VARIANT_KEY"].value_counts()
    id_keys = raw.groupby("SNP")["VARIANT_KEY"].nunique()
    raw["DUPLICATE_VARIANT_KEY"] = raw["VARIANT_KEY"].map(key_counts).gt(1)
    raw["AMBIGUOUS_SNP_ID"] = raw["SNP"].map(id_keys).gt(1)
    clean = raw.sort_values(["P", "N"], ascending=[True, False], kind="mergesort").drop_duplicates("VARIANT_KEY", keep="first")
    clean = clean.loc[(clean["CHR"].between(1, 22)) & clean["RSID"] & (~clean["AMBIGUOUS_SNP_ID"]) & (clean["F"] > F_THRESHOLD)].copy()
    clean.to_csv(OUT / "tmfi_gws_exact_clean.csv", index=False)
    distance = distance_prune(clean)
    distance.to_csv(OUT / "tmfi_gws_distance10mb_provisional.csv", index=False)
    return raw, clean, distance


def match_outcome(candidates, name, filename, mapping):
    wanted = set(candidates["SNP"].astype(str))
    pieces = []
    for chunk in pd.read_csv(RAW / filename, sep="\t", compression="gzip", chunksize=CHUNK, dtype={mapping["snp"]: "string"}):
        hit = chunk.loc[chunk[mapping["snp"]].astype(str).isin(wanted)].copy()
        if not hit.empty:
            pieces.append(hit)
    outcome = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()
    if outcome.empty:
        outcome.to_csv(OUT / f"{name}_candidate_matches_raw.csv", index=False)
        return outcome, pd.DataFrame(), {"matched_rows": 0, "harmonised": 0}

    outcome.to_csv(OUT / f"{name}_candidate_matches_raw.csv", index=False)
    renamed = outcome.rename(columns={
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
    })
    renamed = renamed[["CHR_OUT", "POS_OUT", "SNP", "EA_OUT", "OA_OUT", "EAF_OUT", "BETA_OUT", "SE_OUT", "P_OUT", "N_OUT"]]
    merged = candidates.merge(renamed, on="SNP", how="left", validate="one_to_one")
    statuses = []
    beta_h = []
    eaf_h = []
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
    merged.to_csv(OUT / f"{name}_distance10mb_harmonised.csv", index=False)
    usable = merged.loc[merged["HARMONISE_STATUS"].isin(["aligned", "swapped"])].copy()
    return outcome, usable, {
        "matched_rows": int(len(outcome)),
        "candidate_rows": int(len(candidates)),
        "harmonised": int(len(usable)),
        "status_counts": {str(k): int(v) for k, v in merged["HARMONISE_STATUS"].value_counts().items()},
    }


def main():
    raw, clean, distance = extract_tmfi()
    full_raw, full_use, full_summary = match_outcome(
        distance,
        "cpi_full",
        "Chronic_pain_intensity.txt.gz",
        {"chr": "CHR", "pos": "BP", "snp": "SNP", "ea": "A1", "oa": "A2", "freq": "MAF", "beta": "BETA", "se": "SE", "p": "P", "n": "N"},
    )
    no_raw, no_use, no_summary = match_outcome(
        distance,
        "cpi_no_imaging",
        "chronic_pain_intensity_no_imaging.txt.gz",
        {"chr": "#CHROM", "pos": "POS", "snp": "ID", "ea": "A1", "oa": "OMITTED", "freq": "A1_FREQ", "beta": "BETA", "se": "SE", "p": "P", "n": "OBS_CT"},
    )
    summary = {
        "parameters": {"p_threshold": P_THRESHOLD, "f_threshold": F_THRESHOLD, "provisional_distance_bp": DISTANCE_BP},
        "tmfi": {
            "gws_rows_raw": int(len(raw)),
            "unique_variant_keys": int(raw["VARIANT_KEY"].nunique()),
            "duplicate_variant_key_excess_rows": int(len(raw) - raw["VARIANT_KEY"].nunique()),
            "ambiguous_snp_id_rows": int(raw["AMBIGUOUS_SNP_ID"].sum()),
            "rsid_rows": int(raw["RSID"].sum()),
            "f_gt_10_rows_raw": int((raw["F"] > F_THRESHOLD).sum()),
            "exact_clean_strong_rsid_rows": int(len(clean)),
            "distance10mb_provisional_rows": int(len(distance)),
            "f_min_distance_set": float(distance["F"].min()),
            "f_median_distance_set": float(distance["F"].median()),
            "f_max_distance_set": float(distance["F"].max()),
            "palindromic_rows_distance_set": int(distance["PALINDROMIC"].sum()),
            "palindromic_intermediate_rows_distance_set": int((distance["PALINDROMIC"] & (distance["MAF"] >= 0.42)).sum()),
        },
        "cpi_full": full_summary,
        "cpi_no_imaging": no_summary,
        "important_note": "10 Mb physical distance pruning is a provisional feasibility screen, not LD clumping. Publication-grade instruments require EUR LD clumping.",
    }
    (OUT / "stage1_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
