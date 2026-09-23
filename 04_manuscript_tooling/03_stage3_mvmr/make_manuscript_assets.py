"""Create the manuscript table and publication figures for the Stage 3 paper.

Backend: Python/matplotlib only. All quantitative panels are generated from
the frozen Stage 3/3B CSV outputs; no mock data are used.
"""

from math import erfc, sqrt
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parent
FIG = BASE / "figures"
TABLE = BASE / "tables"
FIG.mkdir(exist_ok=True)
TABLE.mkdir(exist_ok=True)

COLORS = {"TMFI": "#3465A4", "BMI": "#D55E00"}
GREY = "#6E6E6E"
mpl = None
plt = None


def outcome_label(outcome):
    return {"CPI_full": "CPI full", "CPI_no_imaging": "CPI no imaging"}.get(outcome, outcome)


def sd_scale(outcome, exposure):
    """Multiplier taking a raw MVMR beta to SD-of-outcome per SD-of-exposure.

    Read from the frozen Stage 2B power outputs so the figure, Table 3 and the
    power section cannot drift apart.
    """
    import json

    payload = json.loads(
        (BASE / "stage2b_power" / "stage2b_mvmr_power_summary.json").read_text(encoding="utf-8")
    )
    record = payload["results"][f"{outcome}|{exposure}"]
    return record["implied_exposure_sd"] / record["implied_outcome_sd"]


def add_standardized(df):
    """Add *_std columns (SD per SD) alongside the raw-unit estimates."""
    df = df.copy()
    df["scale"] = [sd_scale(o, e) for o, e in zip(df["outcome"], df["exposure"])]
    for col in ("beta", "se", "ci_low", "ci_high"):
        df[col + "_std"] = df[col] * df["scale"]
    return df


def init_plotting():
    global mpl, plt
    try:
        import matplotlib as _mpl
        _mpl.use("Agg")
        import matplotlib.pyplot as _plt
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "matplotlib is required only for figure rendering. The result table "
            "was generated successfully; install matplotlib and rerun this script "
            "to render figures."
        ) from exc
    mpl, plt = _mpl, _plt
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "legend.frameon": False,
        }
    )


def save_pub(fig, stem):
    fig.savefig(FIG / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)


def primary_table():
    r = pd.read_csv(BASE / "R_validation_ivw.csv")
    r["ci_low"] = r["beta"] - 1.96 * r["se"]
    r["ci_high"] = r["beta"] + 1.96 * r["se"]
    r["analysis"] = "Primary MVMR-IVW (MRE)"
    r["estimate_unit"] = "CPI raw outcome units per one exposure unit"
    r["n_snp"] = 93
    # p is the R MVMR t-test value; the manuscript declares the standard-normal
    # rule (z = beta/SE), so both are exported to keep the package self-consistent.
    r["p_normal"] = [float(erfc(abs(float(b) / float(s)) / sqrt(2.0))) for b, s in zip(r["beta"], r["se"])]
    r = r[["analysis", "outcome", "exposure", "n_snp", "beta", "se", "ci_low", "ci_high", "p", "p_normal", "estimate_unit"]]

    ext = pd.read_csv(BASE / "stage3b_external_exclusion_results.csv")
    ext = ext.loc[ext["method"].isin(["exclude:external_metabolic_pleiotropy", "exclude:all_five_review_targets"])].copy()
    ext = ext.rename(columns={"method": "analysis", "p_normal": "p"})
    ext["p_normal"] = ext["p"]
    ext["analysis"] = ext["analysis"].map(
        {
            "exclude:external_metabolic_pleiotropy": "Sensitivity: exclude rs731839 + rs13107325",
            "exclude:all_five_review_targets": "Sensitivity: exclude five review targets",
        }
    )
    ext["n_snp"] = ext["nsnp"]
    ext["estimate_unit"] = "CPI raw outcome units per one exposure unit"
    ext = ext[["analysis", "outcome", "exposure", "n_snp", "beta", "se", "ci_low", "ci_high", "p", "p_normal", "estimate_unit"]]
    result = pd.concat([r, ext], ignore_index=True)
    result = add_standardized(result)
    result.to_csv(TABLE / "main_results_table.csv", index=False)

    lines = [
        "| Analysis | Outcome | Exposure | SNPs | beta | MRE SE | 95% CI | P (t, R MVMR) | P (normal, declared rule) |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in result.itertuples(index=False):
        lines.append(
            f"| {row.analysis} | {row.outcome} | {row.exposure} | {row.n_snp} | "
            f"{row.beta:.4f} | {row.se:.4f} | [{row.ci_low:.4f}, {row.ci_high:.4f}] | "
            f"{row.p:.3g} | {row.p_normal:.3g} |"
        )
    lines.extend(
        [
            "",
            "*Note.* Confidence intervals are estimate ± 1.96 × MRE SE. Two P values are shown because two conventions exist in the pipeline: `P (t, R MVMR)` is the Student-t value returned by the official R `MVMR` package, and `P (normal, declared rule)` is the two-sided standard-normal value applied to z = beta/SE, which is the rule declared in the manuscript's Statistical reporting section. The two differ only in the last one to two digits of significance for the BMI coefficients and change no conclusion. Estimates and confidence intervals retain the source exposure and CPI units.",
        ]
    )
    (TABLE / "main_results_table.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return add_standardized(r), add_standardized(ext)


def forest_panel(ax, data, title):
    data = data.copy().reset_index(drop=True)
    y = np.arange(len(data))[::-1]
    for yi, row in zip(y, data.itertuples(index=False)):
        # Standardised (SD per SD) so that the two exposures share one comparable axis.
        ax.errorbar(
            row.beta_std,
            yi,
            xerr=[[row.beta_std - row.ci_low_std], [row.ci_high_std - row.beta_std]],
            fmt="o",
            color=COLORS[row.exposure],
            ecolor=COLORS[row.exposure],
            markersize=4.2,
            capsize=2.0,
            linewidth=1.1,
        )
    ax.axvline(0, color="#333333", linewidth=0.8, linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{outcome_label(r.outcome)} · {r.exposure}" for r in data.itertuples(index=False)])
    ax.set_title(title, loc="left", fontsize=8, fontweight="bold")
    ax.set_xlabel("MVMR direct effect (SD of pain intensity per 1 SD of exposure)", fontsize=8)
    ax.grid(axis="x", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    ax.set_axisbelow(True)


def make_forest(primary, ext):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), sharex=True)
    forest_panel(axes[0], primary, "Primary MVMR-IVW (MRE)")
    sens = ext.loc[ext["analysis"] == "Sensitivity: exclude rs731839 + rs13107325"]
    forest_panel(axes[1], sens, "Sensitivity: metabolic pleiotropy targets")
    for ax in axes:
        ax.tick_params(axis="y", labelsize=8)
        ax.tick_params(axis="x", labelsize=8)
    fig.tight_layout()
    save_pub(fig, "Figure2_primary_and_sensitivity_forest")


def make_leave_one_out(primary):
    loo = pd.read_csv(BASE / "stage3b_leave_one_out.csv")
    loo["snp"] = loo["method"].str.replace("leave_one_out:", "", regex=False)
    full_order = pd.read_csv(BASE / "cpi_full_mvmr_harmonised_usable.csv")["SNP"].astype(str).tolist()
    rank = {snp: i for i, snp in enumerate(full_order)}
    loo["rank"] = loo["snp"].map(rank)
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.2), sharex=False)
    for i, outcome in enumerate(["CPI_full", "CPI_no_imaging"]):
        for j, exposure in enumerate(["TMFI", "BMI"]):
            ax = axes[i, j]
            d = loo.loc[(loo.outcome == outcome) & (loo.exposure == exposure)].sort_values("rank")
            ax.errorbar(
                d.beta,
                np.arange(len(d)),
                xerr=[d.beta - d.ci_low, d.ci_high - d.beta],
                fmt="o",
                markersize=2.4,
                linewidth=0.45,
                color=COLORS[exposure],
                ecolor=COLORS[exposure],
                alpha=0.78,
            )
            p = primary.loc[(primary.outcome == outcome) & (primary.exposure == exposure)].iloc[0]
            ax.axvline(p.beta, color="#222222", linewidth=0.9, label="Primary")
            ax.axvline(0, color=GREY, linewidth=0.7, linestyle="--")
            ax.set_title(f"{outcome_label(outcome)} · {exposure}", loc="left", fontsize=8, fontweight="bold")
            ax.set_xlabel("Leave-one-out beta")
            ax.set_ylabel("SNP index" if j == 0 else "")
            ax.set_yticks([])
            ax.grid(axis="x", color="#D9D9D9", linewidth=0.5, alpha=0.7)
            ax.set_axisbelow(True)
    fig.suptitle("The direct-effect pattern is not driven by a single SNP", fontsize=9, fontweight="bold", y=1.01)
    fig.tight_layout()
    save_pub(fig, "Figure2_leave_one_out")


def make_influence_sensitivity():
    infl = pd.read_csv(BASE / "stage3b_influence_diagnostics.csv")
    sens = pd.read_csv(BASE / "stage3b_robustness_results.csv")
    methods = [
        "primary_IVW_MRE",
        "exploratory_Huber_IRLS",
        "exclude_Cook_gt_4_over_n",
        "exclude_top_1_Cook",
        "exclude_top_3_Cook",
    ]
    short = {
        "primary_IVW_MRE": "Primary",
        "exploratory_Huber_IRLS": "Huber",
        "exclude_Cook_gt_4_over_n": "Cook threshold",
        "exclude_top_1_Cook": "Top 1 Cook",
        "exclude_top_3_Cook": "Top 3 Cook",
    }
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), gridspec_kw={"width_ratios": [1.05, 1.35]})
    for outcome, color in [("CPI_full", "#4C78A8"), ("CPI_no_imaging", "#F58518")]:
        d = infl.loc[infl.outcome == outcome].nlargest(10, "cook_distance").sort_values("cook_distance")
        axes[0].barh(d.SNP, d.cook_distance, color=color, alpha=0.72, label=outcome_label(outcome))
    axes[0].axvline(4 / 93, color="#333333", linestyle="--", linewidth=0.8, label="Cook threshold (4/n)")
    axes[0].set_xlabel("Cook distance")
    axes[0].set_title("Influence points", loc="left", fontsize=8, fontweight="bold")
    axes[0].legend(fontsize=6)
    axes[0].grid(axis="x", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    axes[0].set_axisbelow(True)

    d = sens.loc[(sens.method.isin(methods)) & (sens.exposure == "BMI")].copy()
    d["method_label"] = d.method.map(short)
    order = [short[m] for m in methods]
    ybase = np.arange(len(order))[::-1]
    for offset, outcome in enumerate(["CPI_full", "CPI_no_imaging"]):
        dd = d.loc[d.outcome == outcome].set_index("method_label").reindex(order).reset_index()
        y = ybase + (0.13 if offset == 0 else -0.13)
        axes[1].errorbar(
            dd.beta,
            y,
            xerr=[dd.beta - dd.ci_low, dd.ci_high - dd.beta],
            fmt="o",
            markersize=4,
            capsize=2,
            color="#4C78A8" if offset == 0 else "#F58518",
            label=outcome_label(outcome),
        )
    axes[1].axvline(0, color="#333333", linestyle="--", linewidth=0.8)
    axes[1].set_yticks(ybase)
    axes[1].set_yticklabels(order)
    axes[1].set_xlabel("BMI direct effect")
    axes[1].set_title("BMI signal survives influence-point checks", loc="left", fontsize=8, fontweight="bold")
    axes[1].legend(fontsize=6)
    axes[1].grid(axis="x", color="#D9D9D9", linewidth=0.5, alpha=0.7)
    axes[1].set_axisbelow(True)
    fig.suptitle("Robustness supports the same interpretation", fontsize=9, fontweight="bold", y=1.02)
    fig.tight_layout()
    save_pub(fig, "Figure4_influence_and_robustness")


def make_design_flowchart():
    """Figure 1: study design and analytic workflow.

    Layout is static (no data are plotted); the same figure is also committed as
    Figure1_study_design.svg so the vector version can be used directly.
    """
    fig, ax = plt.subplots(figsize=(6.8, 5.44))
    ax.set_xlim(0, 680)
    ax.set_ylim(544, 0)
    ax.axis("off")

    def box(x, y, w, h, fc, ec, lw=1.0):
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec,
                                   linewidth=lw, zorder=1))

    def txt(x, y, s, size=9, bold=False, color="#1A1A1A"):
        ax.text(x, y, s, ha="center", va="center", fontsize=size,
                fontweight="bold" if bold else "normal", color=color, zorder=2)

    def arrow(x1, y1, x2, y2):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1), zorder=0,
                    arrowprops=dict(arrowstyle="-|>", color="#6B7A8A",
                                    linewidth=1.2, shrinkA=0, shrinkB=0))

    # Row A - data sources
    for x in (25, 245, 465):
        box(x, 24, 190, 56, "#EDF2F7", "#8FA8BF")
    txt(120, 38, "UK Biobank", 11.5, True)
    txt(120, 52, "Thigh muscle fat infiltration", 9, False, "#4A4A4A")
    txt(120, 66, "N = 46,347", 9, False, "#4A4A4A")
    txt(340, 38, "GIANT 2015", 11.5, True)
    txt(340, 52, "Body mass index (BMI)", 9, False, "#4A4A4A")
    txt(340, 66, "N \u2264 339,224", 9, False, "#4A4A4A")
    txt(560, 38, "UK Biobank", 11.5, True)
    txt(560, 52, "Chronic pain intensity", 9, False, "#4A4A4A")
    txt(560, 66, "N = 134,627; 106,697", 9, False, "#4A4A4A")

    # Row B - instrument construction
    box(100, 98, 480, 58, "#F6F6F3", "#A8A89E")
    txt(340, 111, "Instrument construction", 11.5, True)
    txt(340, 136, "$P < 5\\times10^{-8}$; LD clumped r\u00b2 < 0.001 / 10 Mb; strand aligned on allele frequency",
        9, False, "#4A4A4A")

    # Row C - instrument sets
    box(115, 162, 210, 46, "#FFFFFF", "#8FA8BF", 1.1)
    box(355, 162, 210, 46, "#FFFFFF", "#8FA8BF", 1.1)
    txt(220, 176, "42 instruments", 11, True, "#3465A4")
    txt(220, 195, "TMFI (univariable)", 9, False, "#4A4A4A")
    txt(460, 176, "93 instruments", 11, True, "#D55E00")
    txt(460, 195, "TMFI + BMI (joint)", 9, False, "#4A4A4A")

    # Row D - analysis stages
    box(25, 222, 300, 62, "#F2F6FA", "#3465A4", 1.2)
    box(355, 222, 300, 62, "#FCF4EE", "#D55E00", 1.2)
    txt(175, 239, "Stage 1 \u2014 univariable MR", 11.5, True)
    txt(175, 253, "42 instruments", 9, False, "#4A4A4A")
    txt(175, 267, "TMFI \u2192 CPI (2 definitions)", 9, False, "#4A4A4A")
    txt(505, 239, "Stage 2 \u2014 multivariable MR", 11.5, True)
    txt(505, 253, "93 instruments", 9, False, "#4A4A4A")
    txt(505, 267, "TMFI + BMI \u2192 CPI (2 definitions)", 9, False, "#4A4A4A")

    # Row E - site-matched secondary analysis
    box(25, 306, 630, 62, "#F1F6F1", "#5F8A5F", 1.2)
    txt(340, 320, "Pre-specified site-matched secondary analysis", 11.5, True)
    txt(340, 335, "Knee pain (98,704 cases) \u00b7 back pain (118,471 cases)", 9, False, "#4A4A4A")
    txt(340, 348, "Second knee GWAS (N = 429,737): same-cohort Stage 1 concordance check",
        8.5, False, "#4A4A4A")
    txt(340, 360, "Frozen 42 / 93 instrument sets; BMI as internal comparator",
        8.5, False, "#4A4A4A")

    # Row F - validity audits
    for x in (25, 245, 465):
        box(x, 390, 190, 56, "#F5F5F5", "#909090")
    txt(120, 410, "Sample-overlap audit", 10.5, True)
    txt(120, 426, "LDSC intercept vs ceiling", 9, False, "#4A4A4A")
    txt(340, 410, "Power bound", 10.5, True)
    txt(340, 426, "MDE vs recovered BMI effect", 9, False, "#4A4A4A")
    txt(560, 410, "Reproducibility", 10.5, True)
    txt(560, 426, "Independent base-R re-run", 9, False, "#4A4A4A")

    # Row G - primary result
    box(25, 464, 630, 64, "#FBFBFB", "#4A4A4A", 1.2)
    txt(340, 477, "Primary result (CPI_full, Stage 2, per SD of exposure)", 10.5, True)
    txt(340, 497, "BMI: +0.090 SD (95% CI 0.038\u20130.143; $P = 7.5\\times10^{-4}$)", 9.5, True, "#D55E00")
    txt(340, 516, "TMFI: \u22120.002 SD (95% CI \u22120.036\u20130.031; P = 0.90)", 9.5, True, "#3465A4")

    for a in [(120, 80, 120, 98), (340, 80, 340, 98), (560, 80, 560, 98),
              (340, 156, 220, 162), (340, 156, 460, 162),
              (220, 208, 175, 222), (460, 208, 505, 222),
              (175, 284, 175, 306), (505, 284, 505, 306),
              (120, 368, 120, 390), (340, 368, 340, 390), (560, 368, 560, 390),
              (120, 446, 120, 464), (340, 446, 340, 464), (560, 446, 560, 464)]:
        arrow(*a)

    fig.tight_layout()
    save_pub(fig, "Figure1_study_design")


def main():
    primary, ext = primary_table()
    try:
        init_plotting()
    except RuntimeError as exc:
        print(str(exc))
        print(f"Wrote tables to {TABLE}")
        return
    make_design_flowchart()
    make_forest(primary, ext)
    make_leave_one_out(primary)
    make_influence_sensitivity()
    print(f"Wrote tables to {TABLE}")
    print(f"Wrote figures to {FIG}")


if __name__ == "__main__":
    main()
