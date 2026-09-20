import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(r"D:\SX\TMFI_CPI_MR\00_stage0_audit")
INCOMING = BASE / "incoming"
OUT = BASE / "stage0_scan_results.json"
CHUNK_SIZE = 500_000


DATASETS = [
    {
        "name": "TMFI",
        "file": "TMFI_GWAS.csv",
        "sep": ",",
        "chrom": "CHR",
        "pos": "POS",
        "id": "SNP",
        "ea": "A1",
        "oa": "A2",
        "freq": "AF1",
        "beta": "BETA",
        "se": "SE",
        "p": "P",
        "n": "N",
        "stat": None,
        "expected_columns": ["CHR", "SNP", "POS", "A1", "A2", "N", "AF1", "BETA", "SE", "P"],
    },
    {
        "name": "CPI_full",
        "file": "Chronic_pain_intensity.txt.gz",
        "sep": "\t",
        "chrom": "CHR",
        "pos": "BP",
        "id": "SNP",
        "ea": "A1",
        "oa": "A2",
        "freq": "MAF",
        "beta": "BETA",
        "se": "SE",
        "p": "P",
        "n": "N",
        "stat": "Z",
        "expected_columns": ["CHR", "BP", "SNP", "A1", "A2", "MAF", "BETA", "SE", "Z", "P", "N"],
    },
    {
        "name": "CPI_no_imaging",
        "file": "chronic_pain_intensity_no_imaging.txt.gz",
        "sep": "\t",
        "chrom": "#CHROM",
        "pos": "POS",
        "id": "ID",
        "ea": "A1",
        "oa": "OMITTED",
        "freq": "A1_FREQ",
        "beta": "BETA",
        "se": "SE",
        "p": "P",
        "n": "OBS_CT",
        "stat": "T_STAT",
        "expected_columns": ["#CHROM", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF?", "A1", "OMITTED", "A1_FREQ", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE"],
    },
]


def finite_min(current, values):
    values = values[np.isfinite(values)]
    if values.size == 0:
        return current
    value = float(values.min())
    return value if current is None else min(current, value)


def finite_max(current, values):
    values = values[np.isfinite(values)]
    if values.size == 0:
        return current
    value = float(values.max())
    return value if current is None else max(current, value)


def audit_dataset(spec):
    path = INCOMING / spec["file"]
    result = {
        "dataset": spec["name"],
        "file": spec["file"],
        "bytes": path.stat().st_size,
        "rows": 0,
        "columns": None,
        "columns_match_expected": None,
        "chunks": 0,
        "chromosomes": {},
        "position_min": None,
        "position_max": None,
        "missing_or_non_numeric": {},
        "p_min": None,
        "p_max": None,
        "p_out_of_range": 0,
        "se_min": None,
        "se_max": None,
        "se_nonpositive": 0,
        "freq_min": None,
        "freq_max": None,
        "freq_out_of_range": 0,
        "beta_min": None,
        "beta_max": None,
        "beta_nonfinite": 0,
        "n_min": None,
        "n_max": None,
        "invalid_allele_rows": 0,
        "identical_effect_other_allele_rows": 0,
        "palindromic_snp_rows": 0,
        "palindromic_snp_maf_0_42_0_5_rows": 0,
        "rsid_rows": 0,
        "non_rsid_rows": 0,
        "duplicate_id_rows": None,
        "duplicate_coordinate_rows": None,
        "statistic_consistency_max_abs_diff": None,
        "plink_test_non_add_rows": None,
        "plink_error_rows": None,
        "full_stream_parse": "PASS",
    }
    id_hashes = []
    coord_codes = []
    allele_re = re.compile(r"^[ACGT]+$")
    numeric_fields = [spec[k] for k in ("chrom", "pos", "freq", "beta", "se", "p", "n")]
    if spec["stat"]:
        numeric_fields.append(spec["stat"])
    if spec["name"] == "CPI_no_imaging":
        result["plink_test_non_add_rows"] = 0
        result["plink_error_rows"] = 0

    reader = pd.read_csv(
        path,
        sep=spec["sep"],
        compression="infer",
        chunksize=CHUNK_SIZE,
        dtype="string",
        keep_default_na=False,
        na_filter=False,
        engine="c",
    )

    for chunk in reader:
        result["chunks"] += 1
        if result["columns"] is None:
            result["columns"] = list(chunk.columns)
            result["columns_match_expected"] = list(chunk.columns) == spec["expected_columns"]
        result["rows"] += len(chunk)

        numeric = {}
        for col in numeric_fields:
            numeric[col] = pd.to_numeric(chunk[col], errors="coerce").to_numpy(dtype="float64")
            result["missing_or_non_numeric"][col] = result["missing_or_non_numeric"].get(col, 0) + int(np.isnan(numeric[col]).sum())

        chrom = numeric[spec["chrom"]]
        pos = numeric[spec["pos"]]
        valid_coord = np.isfinite(chrom) & np.isfinite(pos) & (chrom >= 0) & (pos >= 0)
        if valid_coord.any():
            chrom_int = chrom[valid_coord].astype(np.uint64)
            pos_int = pos[valid_coord].astype(np.uint64)
            coord_codes.append((chrom_int << np.uint64(32)) | pos_int)
            vals, counts = np.unique(chrom_int, return_counts=True)
            for val, count in zip(vals, counts):
                key = str(int(val))
                result["chromosomes"][key] = result["chromosomes"].get(key, 0) + int(count)
            result["position_min"] = finite_min(result["position_min"], pos[valid_coord])
            result["position_max"] = finite_max(result["position_max"], pos[valid_coord])

        ids = chunk[spec["id"]].astype("string")
        id_hashes.append(pd.util.hash_pandas_object(ids, index=False).to_numpy(dtype=np.uint64))
        is_rsid = ids.str.match(r"^rs\d+$", na=False).to_numpy()
        result["rsid_rows"] += int(is_rsid.sum())
        result["non_rsid_rows"] += int((~is_rsid).sum())

        ea = chunk[spec["ea"]].str.upper()
        oa = chunk[spec["oa"]].str.upper()
        valid_ea = ea.str.match(allele_re, na=False).to_numpy()
        valid_oa = oa.str.match(allele_re, na=False).to_numpy()
        result["invalid_allele_rows"] += int((~(valid_ea & valid_oa)).sum())
        result["identical_effect_other_allele_rows"] += int((ea == oa).sum())
        pal = (((ea == "A") & (oa == "T")) | ((ea == "T") & (oa == "A")) | ((ea == "C") & (oa == "G")) | ((ea == "G") & (oa == "C"))).to_numpy()
        result["palindromic_snp_rows"] += int(pal.sum())

        freq = numeric[spec["freq"]]
        result["freq_min"] = finite_min(result["freq_min"], freq)
        result["freq_max"] = finite_max(result["freq_max"], freq)
        result["freq_out_of_range"] += int((np.isfinite(freq) & ((freq < 0) | (freq > 1))).sum())
        result["palindromic_snp_maf_0_42_0_5_rows"] += int((pal & np.isfinite(freq) & (np.minimum(freq, 1 - freq) >= 0.42)).sum())

        p = numeric[spec["p"]]
        result["p_min"] = finite_min(result["p_min"], p)
        result["p_max"] = finite_max(result["p_max"], p)
        result["p_out_of_range"] += int((np.isfinite(p) & ((p < 0) | (p > 1))).sum())

        se = numeric[spec["se"]]
        result["se_min"] = finite_min(result["se_min"], se)
        result["se_max"] = finite_max(result["se_max"], se)
        result["se_nonpositive"] += int((np.isfinite(se) & (se <= 0)).sum())

        beta = numeric[spec["beta"]]
        result["beta_min"] = finite_min(result["beta_min"], beta)
        result["beta_max"] = finite_max(result["beta_max"], beta)
        result["beta_nonfinite"] += int((~np.isfinite(beta)).sum())

        n = numeric[spec["n"]]
        result["n_min"] = finite_min(result["n_min"], n)
        result["n_max"] = finite_max(result["n_max"], n)

        if spec["stat"]:
            stat = numeric[spec["stat"]]
            ok = np.isfinite(beta) & np.isfinite(se) & (se != 0) & np.isfinite(stat)
            if ok.any():
                diff = np.abs(beta[ok] / se[ok] - stat[ok])
                value = float(diff.max())
                old = result["statistic_consistency_max_abs_diff"]
                result["statistic_consistency_max_abs_diff"] = value if old is None else max(old, value)

        if spec["name"] == "CPI_no_imaging":
            result["plink_test_non_add_rows"] += int((chunk["TEST"] != "ADD").sum())
            result["plink_error_rows"] += int((chunk["ERRCODE"] != ".").sum())

    all_id_hashes = np.concatenate(id_hashes)
    all_coords = np.concatenate(coord_codes)
    result["duplicate_id_rows"] = int(all_id_hashes.size - np.unique(all_id_hashes).size)
    result["duplicate_coordinate_rows"] = int(all_coords.size - np.unique(all_coords).size)
    return result


def main():
    output = {
        "audit_version": "stage0_scan_v1",
        "chunk_size": CHUNK_SIZE,
        "datasets": [],
    }
    for spec in DATASETS:
        output["datasets"].append(audit_dataset(spec))
    OUT.write_text(json.dumps(output, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(output, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
