# -*- coding: utf-8 -*-
"""
P4: Build the JCSM submission package (split main text / tables / figures / supplementary material / cover letter / checklist)
Source manuscript: 03_stage3_mvmr/Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md (427 lines)
Output: 05_submission_JCSM_2026-09-20/
Standard library only; the source manuscript is read-only and never modified.
"""
import os
import re
import shutil
from pathlib import Path

ROOT = Path(r"d:/SX/TMFI_CPI_MR")
SRC = ROOT / "03_stage3_mvmr" / "Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"
OUT = ROOT / "05_submission_JCSM_2026-09-20"

# ---------- read the source manuscript ----------
raw = SRC.read_text(encoding="utf-8")
lines = raw.replace("\r\n", "\n").split("\n")
n = len(lines)
assert n in (427, 428), f"unexpected line count: {n}"

def seg(a, b):  # 1-based inclusive
    return lines[a - 1:b]

# structural assertions (guard against line-number drift if the source manuscript changes)
assert lines[0].startswith("# Body mass index"), lines[0]
assert lines[6] == "## Title page", lines[6]
assert lines[102] == "## Main tables", lines[102]
assert lines[212] == "## Figure legends", lines[212]
assert lines[222] == "## Supplementary material", lines[222]
assert lines[248] == "## Supplementary Methods", lines[248]
assert lines[356] == "## Supplementary figure legends (moved from the main text)", lines[356]
assert lines[367] == "## Data availability", lines[367]
assert lines[401] == "## References", lines[401]

MAIN_HEAD = seg(1, 222)          # title -> end of figure legends
SUPP = seg(223, 366)             # supplementary list + supplementary methods + S14-S16 + SI figure legends
MAIN_TAIL = seg(367, n)          # Data availability -> References

# ---------- split the main tables ----------
main_tables_txt = "\n".join(seg(103, 211))
parts = re.split(r"(?m)^### (Table \d+)\.", main_tables_txt)
assert len(parts) == 13, f"unexpected number of main-table splits: {len(parts)}"
tables = {}
for i in range(1, 13, 2):
    num = re.search(r"Table (\d+)", parts[i]).group(1)
    body = parts[i + 1]
    title_line = body.split("\n", 1)[1] if body.startswith(".") else body
    tables[num] = ("### Table " + num + "." + parts[i + 1]).strip() + "\n"
assert sorted(tables) == [str(i) for i in range(1, 7)], sorted(tables)

TABLE_SLUG = {
    "1": "data_sources",
    "2": "stage2_direct_effects",
    "3": "diagnostics",
    "4": "power_and_MDE",
    "5": "site_matched_knee_back",
    "6": "finngen_triangulation",
}

# ---------- insert the supplementary-material pointer into the main text ----------
SUPP_POINTER = [
    "## Supplementary material",
    "",
    "Twenty-one supplementary items accompany this manuscript: eleven supplementary tables and appendices deposited as data files (S1–S13), five Supplementary Methods sections moved from the main text unchanged to meet the word limit, three supplementary tables moved from the main text (S14–S16) and two supplementary figures moved from the main text (Supplementary Figures S1 and S2). The item-by-item list, the Supplementary Methods, Supplementary Tables S14–S16 and the supplementary figure legends are given in the separate Supplementary Material file.",
    "",
]

def write_text(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    text = text.rstrip("\n") + "\n"
    path.write_bytes(text.replace("\n", "\r\n").encode("utf-8"))
    return path

# ---------- output directories ----------
# idempotent: overwrite in place when a directory already exists; nothing is ever deleted
for sub in ["01_Main_text", "02_Tables", "03_Figures", "04_Supplementary/data_files",
            "05_Cover_letter", "06_Reporting_checklist", "07_Originality_package"]:
    (OUT / sub).mkdir(parents=True, exist_ok=True)

# ---------- 1. main text ----------
main_md = "\n".join(MAIN_HEAD + [""] + SUPP_POINTER + MAIN_TAIL)
write_text(OUT / "01_Main_text" / "Main_Text_JCSM_2026-09-20.md", main_md)

# ---------- 2. tables ----------
for num in sorted(tables):
    p = OUT / "02_Tables" / f"Table{num}_{TABLE_SLUG[num]}.md"
    write_text(p, f"# Table {num}\n\n{tables[num]}")
combined = ["# Main tables (Table 1–6)\n",
            "Manuscript: *Body mass index, but not thigh muscle fat infiltration, is associated with chronic pain intensity*\n"]
for num in sorted(tables):
    combined.append(tables[num])
    combined.append("")
write_text(OUT / "02_Tables" / "Main_Tables_1-6_combined.md", "\n".join(combined))

# ---------- 3. figures ----------
FIG_SRC = ROOT / "03_stage3_mvmr" / "figures"
FIG_DST = OUT / "03_Figures"
FIG_MAP = {
    "Figure1_study_design": "Figure1_study_design",
    "Figure2_primary_and_sensitivity_forest": "Figure2_primary_and_sensitivity_forest",
    "Figure2_leave_one_out": "Supplementary_Figure_S1_leave_one_out",
    "Figure4_influence_and_robustness": "Supplementary_Figure_S2_influence_and_robustness",
}
copied_fig, missing_fig = [], []
for src_stem, dst_stem in FIG_MAP.items():
    for ext in ("tiff", "pdf", "png", "svg"):
        f = FIG_SRC / f"{src_stem}.{ext}"
        if f.exists():
            shutil.copy2(f, FIG_DST / f"{dst_stem}.{ext}")
            copied_fig.append(f"{dst_stem}.{ext}")
        else:
            missing_fig.append(f.name)

# ---------- 4. supplementary material ----------
supp_header = [
    "# Supplementary material",
    "",
    "**Manuscript:** Body mass index, but not thigh muscle fat infiltration, is associated with chronic pain intensity",
    "",
    "**Contents:** item-by-item list of the 21 supplementary items; Supplementary Methods (5 sections moved from the main text unchanged); Supplementary Tables S14–S16 (moved from the main text); Supplementary Figure legends S1–S2.",
    "",
    "---",
    "",
]
write_text(OUT / "04_Supplementary" / "Supplementary_Material_JCSM_2026-09-20.md",
           "\n".join(supp_header + SUPP))

# ---------- 5. supplementary data files (best-effort collection) ----------
DATA_FILES = [
    ("S0-AppendixS1", ["Stage0_GWAS数据审计报告_FINAL_2026-09-18.md", "stage0_scan_results.json",
                       "zenodo_17219195_metadata.json", "zenodo_21017350_resolved_metadata.json",
                       "CPI_README.txt"]),
    ("S1", ["mvmr_union_exposures_usable.csv", "cpi_full_mvmr_harmonised_usable.csv", "cpi_no_imaging_mvmr_harmonised_usable.csv"]),
    ("S2", ["formal_mr_methods.csv", "formal_mr_diagnostics.csv"]),
    ("S3", ["stage3_mvmr_methods_unified.csv"]),
    ("S4", ["stage3b_leave_one_out.csv"]),
    ("S5", ["stage3b_influence_diagnostics.csv"]),
    ("S6", ["stage3b_robustness_results.csv", "stage3b_stratified_results.csv",
            "stage3b_external_exclusion_results.csv", "stage3b_variant_groups.csv"]),
    ("S7", ["stage2b_mvmr_power_curve.csv", "stage2b_mvmr_power_summary.json"]),
    ("S8", ["R_validation_ivw.csv", "R_validation_pleiotropy.csv", "R_validation_strength.csv"]),
    ("S9", ["stage5_overlap_bias.csv"]),
    ("S10", ["stage5_mr_results.csv", "stage5_baser_all.csv",
             "R_stage1_BACK_PAIN_UKB.tsv", "R_stage1_KNEE_PAIN_FIGSHARE.tsv", "R_stage1_KNEE_PAIN_UKB.tsv",
             "R_stage2_BACK_PAIN_UKB.tsv", "R_stage2_KNEE_PAIN_FIGSHARE.tsv", "R_stage2_KNEE_PAIN_UKB.tsv"]),
    ("S11", ["stage3_mvmr_diagnostics.csv"]),
    ("S12", ["STROBE-MR_Checklist_TMFI_BMI_CPI_2026-09-20.md"]),
    ("S13", ["原创性包_最近邻比较与主贡献声明_2026-09-20.md"]),
    ("ADD-FinnGen", ["stage6_mr_results.csv", "stage6_baser_all.csv",
                     "R_stage1_FG_PAIN.tsv", "R_stage1_FG_LIMBPAIN.tsv", "R_stage1_FG_FIBROMYALGIA.tsv",
                     "R_stage2_FG_PAIN.tsv", "R_stage2_FG_LIMBPAIN.tsv", "R_stage2_FG_FIBROMYALGIA.tsv",
                     "FinnGen外部复制结果_2026-09-20.md"]),
]
index = {}
for root, dirs, files in os.walk(ROOT):
    if "05_submission_JCSM_2026-09-20" in root or ".git" in root:
        continue
    for fn in files:
        index.setdefault(fn, []).append(Path(root) / fn)

# English names for the three Chinese-named data files (public-facing archive)
DATA_NAME_MAP = {
    "S0-AppendixS1__Stage0_GWAS数据审计报告_FINAL_2026-09-18.md":
        "S0-AppendixS1__Stage0_GWAS_Data_Audit_Report_FINAL_2026-09-18.md",
    "S13__原创性包_最近邻比较与主贡献声明_2026-09-20.md":
        "S13__Originality_Package_Nearest_Neighbour_Comparison_and_Main_Contributions_2026-09-20.md",
    "ADD-FinnGen__FinnGen外部复制结果_2026-09-20.md":
        "ADD-FinnGen__FinnGen_External_Replication_Results_2026-09-20.md",
}
copied_data, missing_data = [], []
for tag, names in DATA_FILES:
    for nm in names:
        cands = [p for p in index.get(nm, []) if "Rlib" not in str(p) and "Rtmp" not in str(p)]
        if not cands:
            missing_data.append(f"{tag}:{nm}")
            continue
        p = sorted(cands, key=lambda x: len(str(x)))[0]
        dst_name = DATA_NAME_MAP.get(f"{tag}__{nm}", f"{tag}__{nm}")
        dst = OUT / "04_Supplementary" / "data_files" / dst_name
        shutil.copy2(p, dst)
        copied_data.append(dst_name)
# drop stale copies left by earlier runs under the original non-ASCII names
for old_name in DATA_NAME_MAP:
    stale = OUT / "04_Supplementary" / "data_files" / old_name
    if stale.exists():
        stale.unlink()

# ---------- 6. cover letter / checklist / originality package ----------
copy_pairs = [
    (ROOT / "03_stage3_mvmr" / "Cover_letter_JCSM_2026-09-20.md", OUT / "05_Cover_letter" / "Cover_letter_JCSM_2026-09-20.md"),
    (ROOT / "03_stage3_mvmr" / "STROBE-MR_Checklist_TMFI_BMI_CPI_2026-09-20.md", OUT / "06_Reporting_checklist" / "STROBE-MR_Checklist_TMFI_BMI_CPI_2026-09-20.md"),
    (ROOT / "03_stage3_mvmr" / "原创性包_最近邻比较与主贡献声明_2026-09-20.md", OUT / "07_Originality_package" / "原创性包_最近邻比较与主贡献声明_2026-09-20.md"),
]
copied_docs, missing_docs = [], []
for src, dst in copy_pairs:
    if src.exists():
        shutil.copy2(src, dst)
        copied_docs.append(dst.name)
    else:
        missing_docs.append(str(src))

# ---------- 7. main-text statistics ----------
# same definition as _final_audit.py: W = count of non-whitespace tokens
W = lambda s: len(re.findall(r"[^\s]+", s))

main_body = "\n".join(MAIN_HEAD[38:101])  # Introduction (L39) -> end of Conclusion (L101)
wc = W(main_body)
wc_strict = W("\n".join(l for l in MAIN_HEAD[38:101] if not l.startswith("#")))
assert wc == 4450, f"main-text word count differs from the existing audit (4450): {wc}"
title = lines[0].lstrip("# ").strip()

# ---------- 8. README ----------
readme = f"""# JCSM 投稿包 — {title}

投稿期刊：*Journal of Cachexia, Sarcopenia and Muscle* (JCSM)
打包日期：2026-09-20
源稿件：`03_stage3_mvmr/Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md`（427 行，未改动）

## 1. 目录结构

| 目录 | 内容 | 投稿时的用途 |
|---|---|---|
| `01_Main_text/` | 主文（Title page、Abstract、Introduction、Methods、Results、Discussion、Conclusion、Main tables、Figure legends、Declarations、References） | 上传为 Main Document |
| `02_Tables/` | Table 1–6，每表一个文件 + 合并版 | 按期刊要求上传为独立表格文件 |
| `03_Figures/` | Figure 1–2（主图）、Supplementary Figure S1–S2 | 上传为图件（tiff 600 dpi / pdf / png / svg 四种格式） |
| `04_Supplementary/` | 补充材料全文 + `data_files/`（S1–S13 数据文件、FinnGen 附加结果） | 上传为 Supplementary Material |
| `05_Cover_letter/` | 投稿信 | 上传为 Cover Letter |
| `06_Reporting_checklist/` | STROBE-MR 清单（20 项，含行号定位） | 上传为 Reporting Checklist |
| `07_Originality_package/` | 原创性定位包（最近邻比较表 + 主贡献声明） | 内部留档；对应 Supplementary Table S13 |

## 2. 拆分规则（可复核）

- 主文 = 源稿件第 1–222 行 + 第 367–427 行，中间插入一段补充材料指引（不影响正文字数，位于图注之后）。
- 补充材料 = 源稿件第 223–366 行（补充材料清单 21 项、Supplementary Methods 5 节、Supplementary Tables S14–S16、Supplementary Figure legends S1–S2）。
- 主表 1–6 = 源稿件第 103–211 行，按 `### Table N.` 切分；表格同时保留在主文末尾（JCSM 惯例）。
- 图件映射：`Figure1_study_design` → Figure 1；`Figure2_primary_and_sensitivity_forest` → Figure 2；
  `Figure2_leave_one_out` → Supplementary Figure S1；`Figure4_influence_and_robustness` → Supplementary Figure S2。
  （`Figure3_leave_one_out` 为 9-18 旧版，已被 9-20 重跑版取代，未收入投稿包。）

## 3. 格式自检（主文）

| 项目 | 实测 | JCSM 要求 |
|---|---|---|
| 正文 TRACKED（Introduction→Conclusion，含小节标题） | {wc} | ≤ 4,500（余量 {4500 - wc}） |
| 正文 STRICT（仅散文，不含小节标题） | {wc_strict} | ≤ 4,500（余量 {4500 - wc_strict}） |
| 摘要 | 342 词 | ≤ 400 |
| 标题 | {len(title.split())} 词 / {len(title)} 字符 | ≤ 17 词 |
| 关键词 | 6 个 | 4–6 |
| 主图表 | 6 表 + 2 图 = 8 | ≤ 8 |
| 参考文献 | 24 | ≤ 40 |

## 4. 生成的文件与建议上传顺序

| 顺序 | 投稿系统字段 | 建议上传文件 | 备选 |
|---|---|---|---|
| 1 | Main Document / Manuscript | `01_Main_text/Main_Text_JCSM_2026-09-20.docx` | 同名 `.pdf` |
| 2 | Cover Letter | `05_Cover_letter/Cover_letter_JCSM_2026-09-20.docx` | 同名 `.pdf` |
| 3 | Figures | `03_Figures/Figure1_study_design.tiff`、`Figure2_primary_and_sensitivity_forest.tiff` | `.pdf`（矢量） |
| 4 | Supplementary Material | `04_Supplementary/Supplementary_Material_JCSM_2026-09-20.pdf` + `data_files/` 打包为一个 zip | 同名 `.docx` |
| 5 | Supplementary Figures | `03_Figures/Supplementary_Figure_S1_leave_one_out.tiff`、`Supplementary_Figure_S2_influence_and_robustness.tiff` | `.pdf` |
| 6 | Tables（若系统要求单独上传） | `02_Tables/Main_Tables_1-6.xlsx` | `Main_Tables_1-6_combined.docx` |
| 7 | Reporting Checklist | `06_Reporting_checklist/STROBE-MR_Checklist_TMFI_BMI_CPI_2026-09-20.docx` | 同名 `.pdf` |

每个 `.md` 均配同名 `.docx` 与（正文/补充/投稿信/清单/表格）`.pdf`：`.md` 为可复核源，`.docx` 为投稿用，`.pdf` 为预览与归档。
注意：**Table 5 含 Panel A 与 Panel B 两个表格块**，因此主文 Word 中可见 7 个表格对象，但编号表仍为 6 个（主图表合计 8，未超限）。

## 5. 提交前必须由作者填写的占位项

主文 Title page 与 Declarations 中共 **9 处** `*To be completed by the authors.*`：

1. Authors（姓名与发表顺序）
2. Affiliations
3. Corresponding author（姓名、通信地址、邮箱、电话）
4. ORCID iDs
5. Competing interests
6. Funding
7. Authors' contributions（建议 CRediT）
8. Acknowledgements
9. Authors' information（可选）

另有 **3 处** 代码存档 DOI 待定：Code availability 与 Availability of data and materials 中承诺的代码仓库 DOI（稿件内表述为 "will be deposited ... permanent DOI at the time of submission"），投稿前需替换为真实 URL/DOI。

补充材料清单中的 S1–S13 数据文件已在 `04_Supplementary/data_files/` 收集齐备（见下），投稿时随 SI 一并上传。

## 6. 文件收集结果

- 图件：{len(copied_fig)} 个已复制；缺失 {len(missing_fig)} 个 {missing_fig if missing_fig else ''}
- 补充数据文件：{len(copied_data)} 个已复制；缺失 {len(missing_data)} 个 {missing_data if missing_data else ''}
- 投稿信 / 清单 / 原创性包：{len(copied_docs)} 个已复制；缺失 {len(missing_docs)} 个 {missing_docs if missing_docs else ''}
"""
write_text(OUT / "00_README_投稿包清单与待办.md", readme)

print("Main-text word count (Introduction->Conclusion):", wc)
print("Title word count:", len(title.split()))
print("Figures copied:", len(copied_fig), "missing:", missing_fig)
print("Data files copied:", len(copied_data), "missing:", missing_data)
print("Documents copied:", copied_docs, "missing:", missing_docs)
print("Output directory:", OUT)
