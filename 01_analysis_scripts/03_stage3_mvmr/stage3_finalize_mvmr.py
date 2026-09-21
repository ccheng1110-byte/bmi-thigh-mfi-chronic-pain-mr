import json
from math import erfc, exp, lgamma, log, sqrt
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"D:\SX\TMFI_CPI_MR")
RAW = ROOT / "00_stage0_audit" / "incoming"
OUT = ROOT / "03_stage3_mvmr"
CLUMPED = OUT / "plink_mvmr_union_r2_0p001_kb10000.clumped"
CHUNK = 500_000


def complement(a):
    return str(a).upper().translate(str.maketrans("ACGT", "TGCA"))


def orientation(target_a1, target_a2, source_a1, source_a2):
    ta1, ta2 = str(target_a1).upper(), str(target_a2).upper()
    sa1, sa2 = str(source_a1).upper(), str(source_a2).upper()
    if (sa1, sa2) == (ta1, ta2) or (complement(sa1), complement(sa2)) == (ta1, ta2):
        return 1
    if (sa1, sa2) == (ta2, ta1) or (complement(sa1), complement(sa2)) == (ta2, ta1):
        return -1
    return None


def palindromic(a1, a2):
    return {str(a1).upper(), str(a2).upper()} in ({"A", "T"}, {"C", "G"})


def sign_by_frequency(target_eaf, source_eaf):
    """Orientation for A/T and G/C variants, resolved from allele frequency.

    `orientation()` cannot decide these: both of its branches match, so it
    always returns +1 and the result then depends on nothing but the order in
    which the source file happens to list its alleles.  Here the source allele
    is assigned to whichever strand puts its frequency closest to the target
    frequency.  The target is TMFI (UKB imaging; mean |dEAF| against 1000G EUR
    = 0.017), which is far more precise than the GIANT BMI frequency column
    (multi-cohort average, deviations up to 0.12).
    """
    same = abs(source_eaf - target_eaf)
    opposite = abs(source_eaf - (1.0 - target_eaf))
    return -1 if opposite < same else 1


def read_index_snps():
    return pd.read_csv(CLUMPED, sep=r"\s+")["SNP"].astype(str).tolist()


def prepare_exposures():
    cross = pd.read_csv(OUT / "mvmr_candidate_exposure_associations.csv", dtype={"SNP": "string"})
    index_snps = read_index_snps()
    dat = cross.set_index("SNP").loc[index_snps].reset_index()
    signs, status = [], []
    for row in dat.itertuples(index=False):
        if palindromic(row.TMFI_A1, row.TMFI_A2) and min(row.TMFI_EAF, 1 - row.TMFI_EAF) >= 0.42:
            signs.append(np.nan)
            status.append("drop_palindromic_intermediate")
            continue
        sign = orientation(row.TMFI_A1, row.TMFI_A2, row.BMI_A1, row.BMI_A2)
        if sign is not None and palindromic(row.TMFI_A1, row.TMFI_A2):
            sign = sign_by_frequency(float(row.TMFI_EAF), float(row.BMI_EAF))
        signs.append(sign if sign is not None else np.nan)
        status.append("aligned" if sign == 1 else "swapped" if sign == -1 else "allele_mismatch")
    dat["BMI_ALIGN_SIGN"] = signs
    dat["EXPOSURE_HARMONISE_STATUS"] = status
    dat["BMI_BETA_HARMONISED"] = dat["BMI_BETA"] * dat["BMI_ALIGN_SIGN"]
    dat["BMI_EAF_HARMONISED"] = np.where(dat["BMI_ALIGN_SIGN"] == 1, dat["BMI_EAF"], 1 - dat["BMI_EAF"])
    dat.to_csv(OUT / "mvmr_union_exposures_all.csv", index=False)
    usable = dat.loc[dat["EXPOSURE_HARMONISE_STATUS"].isin(["aligned", "swapped"])].copy()
    usable.to_csv(OUT / "mvmr_union_exposures_usable.csv", index=False)
    return dat, usable


def match_outcome(exposures, name, filename, mapping):
    wanted = set(exposures["SNP"].astype(str))
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
    raw = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame()
    rows = []
    for erow in exposures.itertuples(index=False):
        subset = raw.loc[raw[mapping["snp"]].astype(str) == erow.SNP]
        chosen = None
        sign = None
        for _, orow in subset.iterrows():
            trial = orientation(erow.TMFI_A1, erow.TMFI_A2, orow[mapping["ea"]], orow[mapping["oa"]])
            if trial is not None:
                if palindromic(erow.TMFI_A1, erow.TMFI_A2):
                    trial = sign_by_frequency(float(erow.TMFI_EAF),
                                              float(orow[mapping["freq"]]))
                chosen, sign = orow, trial
                break
        base = erow._asdict()
        if chosen is None:
            base.update({"OUTCOME_STATUS": "missing_or_mismatch", "OUTCOME_BETA_H": np.nan})
        else:
            eaf = float(chosen[mapping["freq"]])
            maf = min(eaf, 1 - eaf)
            if palindromic(erow.TMFI_A1, erow.TMFI_A2) and maf >= 0.42:
                base.update({"OUTCOME_STATUS": "drop_palindromic_intermediate", "OUTCOME_BETA_H": np.nan})
            else:
                base.update(
                    {
                        "OUTCOME_STATUS": "aligned" if sign == 1 else "swapped",
                        "OUTCOME_A1": str(chosen[mapping["ea"]]).upper(),
                        "OUTCOME_A2": str(chosen[mapping["oa"]]).upper(),
                        "OUTCOME_EAF_H": eaf if sign == 1 else 1 - eaf,
                        "OUTCOME_BETA_H": float(chosen[mapping["beta"]]) * sign,
                        "OUTCOME_SE": float(chosen[mapping["se"]]),
                        "OUTCOME_P": float(chosen[mapping["p"]]),
                        "OUTCOME_N": int(chosen[mapping["n"]]),
                    }
                )
        rows.append(base)
    merged = pd.DataFrame(rows)
    merged.to_csv(OUT / f"{name}_mvmr_harmonised_all.csv", index=False)
    usable = merged.loc[merged["OUTCOME_STATUS"].isin(["aligned", "swapped"])].copy()
    usable.to_csv(OUT / f"{name}_mvmr_harmonised_usable.csv", index=False)
    return merged, usable


def conditional_f(target_beta, target_se, other_beta, other_se):
    # Match MVMR::strength_mvmr: unweighted no-intercept regression of the
    # target SNP effects on the other exposure, followed by Q/L.
    delta = float(np.sum(other_beta * target_beta) / np.sum(other_beta**2))
    variance = target_se**2 + (delta**2) * other_se**2
    q_strength = float(np.sum((target_beta - delta * other_beta) ** 2 / variance))
    return {"delta": delta, "q": q_strength, "denominator": len(target_beta), "conditional_f": q_strength / len(target_beta)}


def chi2_sf(x, df):
    """Upper-tail chi-square probability via a stable incomplete-gamma series."""
    a = df / 2.0
    z = x / 2.0
    if z <= 0:
        return 1.0
    if z < a + 1.0:
        term = 1.0 / a
        total = term
        ap = a
        for _ in range(10000):
            ap += 1.0
            term *= z / ap
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
        lower = total * exp(-z + a * log(z) - lgamma(a))
        return max(0.0, min(1.0, 1.0 - lower))
    b = z + 1.0 - a
    c = 1.0 / 1e-300
    d = 1.0 / b
    h = d
    for i in range(1, 10001):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < 1e-300:
            d = 1e-300
        c = b + an / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return max(0.0, min(1.0, exp(-z + a * log(z) - lgamma(a)) * h))


def weighted_model(y, se_y, x, intercept=False):
    design = np.column_stack([np.ones(len(y)), x]) if intercept else x
    w = 1 / se_y**2
    xtw = design.T * w
    bread = np.linalg.inv(xtw @ design)
    beta = bread @ (xtw @ y)
    residual = y - design @ beta
    q = float(np.sum(w * residual**2))
    df = len(y) - design.shape[1]
    phi = max(1.0, q / df)
    se_fixed = np.sqrt(np.diag(bread))
    se_mre = se_fixed * sqrt(phi)
    return beta, se_fixed, se_mre, q, df, phi


def analyze(dat, outcome_name):
    bx = dat[["TMFI_BETA", "BMI_BETA_HARMONISED"]].to_numpy(float)
    sx1 = dat["TMFI_SE"].to_numpy(float)
    sx2 = dat["BMI_SE"].to_numpy(float)
    by = dat["OUTCOME_BETA_H"].to_numpy(float)
    sy = dat["OUTCOME_SE"].to_numpy(float)
    beta, se_fixed, se_mre, q, df, phi = weighted_model(by, sy, bx, intercept=False)
    ebeta, ese_fixed, ese_mre, eq, edf, ephi = weighted_model(by, sy, bx, intercept=True)
    f_tmfi = conditional_f(bx[:, 0], sx1, bx[:, 1], sx2)
    f_bmi = conditional_f(bx[:, 1], sx2, bx[:, 0], sx1)

    methods = []
    for exposure, idx in [("TMFI", 0), ("BMI", 1)]:
        for method, se in [("MVMR_IVW_fixed", se_fixed[idx]), ("MVMR_IVW_multiplicative_random_effects", se_mre[idx])]:
            z = beta[idx] / se
            p = float(erfc(abs(z) / sqrt(2)))
            methods.append(
                {
                    "outcome": outcome_name,
                    "exposure": exposure,
                    "method": method,
                    "nsnp": len(dat),
                    "beta": beta[idx],
                    "se": se,
                    "p": p,
                    "ci_low": beta[idx] - 1.96 * se,
                    "ci_high": beta[idx] + 1.96 * se,
                }
            )
    validity_variance = sy**2 + beta[0] ** 2 * sx1**2 + beta[1] ** 2 * sx2**2
    validity_q = float(np.sum((by - bx @ beta) ** 2 / validity_variance))
    validity_df = len(dat) - bx.shape[1] - 1
    diagnostics = {
        "outcome": outcome_name,
        "nsnp": len(dat),
        "q_ivw_outcome_only": q,
        "q_ivw_df": df,
        "phi": phi,
        "q_validity_modified": validity_q,
        "q_validity_df": validity_df,
        "q_validity_p": chi2_sf(validity_q, validity_df),
        "egger_intercept": ebeta[0],
        "egger_intercept_se_mre": ese_mre[0],
        "egger_tmfi_beta": ebeta[1],
        "egger_tmfi_se_mre": ese_mre[1],
        "egger_bmi_beta": ebeta[2],
        "egger_bmi_se_mre": ese_mre[2],
        "egger_q": eq,
        "egger_df": edf,
        "tmfi_conditional_f": f_tmfi["conditional_f"],
        "bmi_conditional_f": f_bmi["conditional_f"],
        "tmfi_strength_q": f_tmfi["q"],
        "bmi_strength_q": f_bmi["q"],
    }
    return methods, diagnostics


def main():
    all_exp, exposures = prepare_exposures()
    full_all, full_use = match_outcome(
        exposures,
        "cpi_full",
        "Chronic_pain_intensity.txt.gz",
        {"snp": "SNP", "ea": "A1", "oa": "A2", "freq": "MAF", "beta": "BETA", "se": "SE", "p": "P", "n": "N"},
    )
    no_all, no_use = match_outcome(
        exposures,
        "cpi_no_imaging",
        "chronic_pain_intensity_no_imaging.txt.gz",
        {"snp": "ID", "ea": "A1", "oa": "OMITTED", "freq": "A1_FREQ", "beta": "BETA", "se": "SE", "p": "P", "n": "OBS_CT"},
    )
    methods1, diag1 = analyze(full_use, "CPI_full")
    methods2, diag2 = analyze(no_use, "CPI_no_imaging")
    methods = pd.DataFrame(methods1 + methods2)
    diagnostics = pd.DataFrame([diag1, diag2])
    methods.to_csv(OUT / "stage3_mvmr_methods.csv", index=False)
    diagnostics.to_csv(OUT / "stage3_mvmr_diagnostics.csv", index=False)
    summary = {
        "ld_clumped_union": int(len(all_exp)),
        "exposure_harmonised": int(len(exposures)),
        "exposure_status": {str(k): int(v) for k, v in all_exp["EXPOSURE_HARMONISE_STATUS"].value_counts().items()},
        "cpi_full_usable": int(len(full_use)),
        "cpi_no_imaging_usable": int(len(no_use)),
        "methods": methods.to_dict(orient="records"),
        "diagnostics": diagnostics.to_dict(orient="records"),
        "conditional_f_note": "Calculated to match MVMR::strength_mvmr: unweighted no-intercept exposure regression, first-order variance, Q divided by L, with zero covariance between exposure estimates.",
        "pleiotropy_q_note": "Calculated to match MVMR::pleiotropy_mvmr first-order modified Q with gencov=0 and df=L-K-1.",
    }
    (OUT / "stage3_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
