"""Sensitivity analysis after external phenotype review.

No R is executed. The two primary exclusions are variants with clear
non-BMI metabolic associations in the external review; BMI-locus variants are
reported separately and are not automatically treated as invalid instruments.
"""

import json
from pathlib import Path

import pandas as pd

from stage3b_robustness import fit_ivw, result_rows


OUT = Path(r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr")

SETS = {
    "external_metabolic_pleiotropy": ["rs731839", "rs13107325"],
    "all_five_review_targets": [
        "rs2867110",
        "rs731839",
        "rs543874",
        "rs11822234",
        "rs13107325",
    ],
}


def main():
    outputs = []
    summary = {}
    for filename, outcome in [
        ("cpi_full_mvmr_harmonised_usable.csv", "CPI_full"),
        ("cpi_no_imaging_mvmr_harmonised_usable.csv", "CPI_no_imaging"),
    ]:
        dat = pd.read_csv(OUT / filename)
        summary[outcome] = {}
        for label, snps in SETS.items():
            subset = dat.loc[~dat["SNP"].isin(snps)].copy()
            fit = fit_ivw(subset)
            outputs.extend(result_rows(outcome, f"exclude:{label}", fit, len(subset)))
            summary[outcome][label] = {
                "excluded": snps,
                "nsnp": int(len(subset)),
                "q": float(fit["q"]),
                "df": int(fit["df"]),
                "phi": float(fit["phi"]),
            }
    pd.DataFrame(outputs).to_csv(OUT / "stage3b_external_exclusion_results.csv", index=False)
    (OUT / "stage3b_external_exclusion_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
