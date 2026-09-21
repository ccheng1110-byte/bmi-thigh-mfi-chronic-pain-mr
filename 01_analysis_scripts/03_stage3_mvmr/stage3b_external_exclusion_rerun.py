"""Re-run the targeted-exclusion sensitivity analysis on the CURRENT (post-P0-1)
harmonised Stage 3 tables.

Why this script exists
----------------------
`stage3b_external_exclusion_results.csv` was written before the P0-1
palindromic-orientation fix, so it no longer matches the harmonised inputs that
the manuscript is built from. `stage3b_external_exclusion.py` needs pandas/numpy,
which are not available in the analysis sandbox, so the estimator is re-implemented
here with the standard library only. The maths is identical and small enough to
verify by eye:

    X  = [TMFI_BETA, BMI_BETA_HARMONISED]      (n x 2)
    y  = OUTCOME_BETA_H                        (n)
    W  = diag(1 / OUTCOME_SE^2)
    bread = (X' W X)^-1
    beta  = bread (X' W y)
    Q     = sum( w_i * r_i^2 ),  df = n - 2,  phi = max(1, Q/df)
    se_fixed = sqrt(diag(bread));  se_mre = se_fixed * sqrt(phi)

This is exactly `fit_ivw()` in `stage3b_robustness.py`. P values follow the
single declared convention now used throughout the manuscript: two-sided normal
on z = beta / se_mre, consistent with the +/- 1.96 * SE intervals. Student-t
P values with df = n - 2 are also written so the two rules can be compared.
"""

import csv
import json
import math
from pathlib import Path

OUT = Path(r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr")
Z = 1.959963984540054

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


def normal_p(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def student_t_two_sided_p(t, df):
    """Two-sided P for Student t via the regularised incomplete beta function."""
    return betai(df / 2.0, 0.5, df / (df + t * t))


def betacf(a, b, x):
    """Continued fraction for the incomplete beta function (Numerical Recipes)."""
    maxit, eps, fpmin = 300, 3.0e-16, 1.0e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, maxit + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < eps:
            break
    return h


def betai(a, b, x):
    """Regularised incomplete beta I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    bt = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x) + b * math.log(1.0 - x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * betacf(a, b, x) / a
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b


def load(filename):
    with open(OUT / filename, newline="", encoding="utf-8") as fh:
        return [r for r in csv.DictReader(fh)]


def fit_ivw(rows):
    n = len(rows)
    x1 = [float(r["TMFI_BETA"]) for r in rows]
    x2 = [float(r["BMI_BETA_HARMONISED"]) for r in rows]
    y = [float(r["OUTCOME_BETA_H"]) for r in rows]
    sy = [float(r["OUTCOME_SE"]) for r in rows]
    w = [1.0 / (s * s) for s in sy]

    a11 = sum(wi * a * a for wi, a in zip(w, x1))
    a12 = sum(wi * a * b for wi, a, b in zip(w, x1, x2))
    a22 = sum(wi * b * b for wi, b in zip(w, x2))
    b1 = sum(wi * a * yi for wi, a, yi in zip(w, x1, y))
    b2 = sum(wi * b * yi for wi, b, yi in zip(w, x2, y))

    det = a11 * a22 - a12 * a12
    i11, i12, i22 = a22 / det, -a12 / det, a11 / det
    beta1 = i11 * b1 + i12 * b2
    beta2 = i12 * b1 + i22 * b2

    q = 0.0
    for wi, a, b, yi in zip(w, x1, x2, y):
        r = yi - (a * beta1 + b * beta2)
        q += wi * r * r
    df = n - 2
    phi = max(1.0, q / df)
    se_fixed = (math.sqrt(i11), math.sqrt(i22))
    se_mre = (se_fixed[0] * math.sqrt(phi), se_fixed[1] * math.sqrt(phi))
    return {
        "beta": (beta1, beta2),
        "se_mre": se_mre,
        "se_fixed": se_fixed,
        "q": q,
        "df": df,
        "phi": phi,
        "n": n,
    }


def rows_for(outcome, method, fit):
    out = []
    for label, idx in (("TMFI", 0), ("BMI", 1)):
        b = fit["beta"][idx]
        se = fit["se_mre"][idx]
        z = b / se
        out.append(
            {
                "outcome": outcome,
                "method": method,
                "exposure": label,
                "nsnp": fit["n"],
                "beta": repr(b),
                "se": repr(se),
                "p_normal": repr(normal_p(z)),
                "p_t_df_n_minus_2": repr(student_t_two_sided_p(z, fit["df"])),
                "ci_low": repr(b - Z * se),
                "ci_high": repr(b + Z * se),
                "phi": repr(fit["phi"]),
                "q": repr(fit["q"]),
            }
        )
    return out


def fmt(v, nd=4):
    return f"{float(v):.{nd}f}"


def main():
    outputs = []
    summary = {}
    for filename, outcome in [
        ("cpi_full_mvmr_harmonised_usable.csv", "CPI_full"),
        ("cpi_no_imaging_mvmr_harmonised_usable.csv", "CPI_no_imaging"),
    ]:
        dat = load(filename)
        summary[outcome] = {"nsnp_total": len(dat)}
        # Primary IVW-MRE on the full set, for cross-checking against Table 3/4.
        outputs.extend(rows_for(outcome, "primary_IVW_MRE", fit_ivw(dat)))
        for label, snps in SETS.items():
            subset = [r for r in dat if r["SNP"] not in snps]
            fit = fit_ivw(subset)
            outputs.extend(rows_for(outcome, f"exclude:{label}", fit))
            summary[outcome][label] = {
                "excluded": snps,
                "nsnp": fit["n"],
                "q": fit["q"],
                "df": fit["df"],
                "phi": fit["phi"],
            }

    fields = [
        "outcome", "method", "exposure", "nsnp", "beta", "se",
        "p_normal", "p_t_df_n_minus_2", "ci_low", "ci_high", "phi", "q",
    ]
    with open(OUT / "stage3b_external_exclusion_results.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        wr.writerows(outputs)
    (OUT / "stage3b_external_exclusion_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print(f"{'outcome':<15}{'method':<38}{'exp':<6}{'n':>4} "
          f"{'beta':>10}{'se':>9}{'CI':>26}{'P(z)':>12}{'P(t)':>12}{'phi':>7}")
    for r in outputs:
        ci = f"[{fmt(r['ci_low'])}, {fmt(r['ci_high'])}]"
        print(f"{r['outcome']:<15}{r['method']:<38}{r['exposure']:<6}{r['nsnp']:>4} "
              f"{float(r['beta']):>10.4f}{float(r['se']):>9.4f}{ci:>26}"
              f"{float(r['p_normal']):>12.3g}{float(r['p_t_df_n_minus_2']):>12.3g}"
              f"{float(r['phi']):>7.3f}")


if __name__ == "__main__":
    main()
