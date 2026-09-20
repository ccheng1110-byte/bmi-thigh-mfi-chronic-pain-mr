# -*- coding: utf-8 -*-
"""
P5: Build the delivery folder for the public code repository (Zenodo/GitHub + DOI).
Contents: all analysis scripts, frozen intermediate outputs (mapped one-to-one to supplementary items S0-S13),
      figure-generation code and figures, manuscript-number and manuscript-editing scripts, environment and logs, and a manifest (with MD5).
Copy-only plus generation of README/MANIFEST; no source file is modified.
"""
import hashlib
import os
import re
import shutil
from pathlib import Path

ROOT = Path(r"d:/SX/TMFI_CPI_MR")
SUB = ROOT / "05_submission_JCSM_2026-09-20"
OUT = ROOT / "06_code_repository_2026-09-20"

MAX_BYTES = 25 * 1024 * 1024          # per-file size cap; larger files are skipped and logged
SKIP_DIRS = {"__pycache__", ".git", "Rlib", "Rtmp", ".workbuddy", "incoming",
             "05_submission_JCSM_2026-09-20", "06_code_repository_2026-09-20",
             "reference", "tools", "_en_localized"}
SKIP_EXT = {".gz", ".bed", ".bim", ".fam", ".exe", ".dll", ".zip", ".bz2", ".tar", ".pyc"}
# one-off localisation migration tool: not needed to reproduce the analysis, and its
# Chinese->English mapping table would reintroduce Chinese into the public archive
SKIP_FILES = {"_translate_script_comments.py"}
# English localisation sources for the public archive (preferred over the Chinese originals)
EN_DIR = ROOT / "03_stage3_mvmr" / "_en_localized"

FIG_SRC = ROOT / "03_stage3_mvmr" / "figures"

for d in ["01_analysis_scripts", "02_frozen_outputs", "03_figure_code_and_outputs",
          "04_manuscript_tooling", "05_environment_and_logs"]:
    (OUT / d).mkdir(parents=True, exist_ok=True)


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def copy_to(src: Path, dst: Path, log: list, skipped: list):
    if not src.exists():
        skipped.append(f"MISSING {src}")
        return False
    if src.stat().st_size > MAX_BYTES:
        skipped.append(f"TOO LARGE {src} ({src.stat().st_size} bytes)")
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    log.append(dst)
    return True


# ---------- 1. analysis scripts (stages 00-06) ----------
scripts, tooling, skipped = [], [], []
for stage_dir in sorted(p for p in ROOT.iterdir() if p.is_dir() and re.match(r"^\d\d_stage", p.name)):
    for root, dirs, files in os.walk(stage_dir):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        rel_stage = Path(root).relative_to(ROOT)
        for fn in files:
            src = Path(root) / fn
            ext = src.suffix.lower()
            if ext in SKIP_EXT or fn in SKIP_FILES:
                continue
            # analysis scripts
            if ext in (".py", ".r", ".sh", ".m") or (ext in (".R",)):
                if fn.startswith("_") or fn in ("make_manuscript_assets.py", "manuscript_numbers.py"):
                    if copy_to(src, OUT / "04_manuscript_tooling" / rel_stage / fn, tooling, skipped):
                        pass
                else:
                    copy_to(src, OUT / "01_analysis_scripts" / rel_stage / fn, scripts, skipped)
            # stage logs / summaries / manifests (small text files)
            elif ext in (".json", ".log", ".txt", ".tsv", ".md") and "metadata" not in str(rel_stage).lower():
                if "sessionInfo" in fn or fn.endswith("_summary.json") or fn.endswith("_scan_results.json") \
                        or fn.endswith("manifest.tsv") or "plink" in fn.lower() or fn.endswith("_results.json"):
                    copy_to(src, OUT / "05_environment_and_logs" / rel_stage / fn, scripts, skipped)

# manuscript-level tooling (the _*.py files in 03_stage3 were already routed to 04 above, by stage)
extra_tooling = [
    ROOT / "03_stage3_mvmr" / "citation_dois.txt",
    ROOT / "03_stage3_mvmr" / "manuscript_numbers.py",
    ROOT / "03_stage3_mvmr" / "make_manuscript_assets.py",
]
for src in extra_tooling:
    copy_to(src, OUT / ("04_manuscript_tooling" if src.name != "make_manuscript_assets.py"
                        else "03_figure_code_and_outputs") / "03_stage3_mvmr" / src.name, tooling, skipped)

# ---------- 2. frozen intermediate outputs (mapped one-to-one to supplementary items S0-S13) ----------
frozen = []
df_src = SUB / "04_Supplementary" / "data_files"
if df_src.exists():
    # Public archive uses English filenames only; source files keep their original names
    NAME_MAP = {
        "S0-AppendixS1__Stage0_GWAS数据审计报告_FINAL_2026-09-18.md":
            "S0-AppendixS1__Stage0_GWAS_Data_Audit_Report_FINAL_2026-09-18.md",
        "S13__原创性包_最近邻比较与主贡献声明_2026-09-20.md":
            "S13__Originality_Package_Nearest_Neighbour_Comparison_and_Main_Contributions_2026-09-20.md",
        "ADD-FinnGen__FinnGen外部复制结果_2026-09-20.md":
            "ADD-FinnGen__FinnGen_External_Replication_Results_2026-09-20.md",
    }
    for p in sorted(df_src.iterdir()):
        if p.is_file():
            dst_name = NAME_MAP.get(p.name, p.name)
            if any(ord(c) > 127 for c in dst_name):
                skipped.append(f"NON-ASCII FILENAME (no mapping defined) {p.name}")
                continue
            # prefer an English-localised copy if one exists in _en_localized/
            src_file = EN_DIR / dst_name if (EN_DIR / dst_name).exists() else p
            copy_to(src_file, OUT / "02_frozen_outputs" / dst_name, frozen, skipped)
    # drop stale copies left by earlier runs under the original non-ASCII names
    for old_name in NAME_MAP:
        stale = OUT / "02_frozen_outputs" / old_name
        if stale.exists():
            stale.unlink()
else:
    skipped.append("MISSING data_files dir (run _build_submission_package.py first)")

# additional: key frozen inputs cited by the manuscript (not in the S list, but required for reproduction)
extra_frozen = [
    ("03_stage3_mvmr", "mvmr_union_exposures_all.csv"),
    ("03_stage3_mvmr", "mvmr_clump_input.tsv"),
    ("03_stage3_mvmr", "plink_mvmr_union_r2_0p001_kb10000.clumped"),
    ("03_stage3_mvmr", "stage3_mvmr_methods.csv"),
    ("03_stage3_mvmr", "stage3b_summary.json"),
    ("03_stage3_mvmr", "tables/main_results_table.csv"),
    ("01_stage1_instruments", "tmfi_gws_ldclumped_r2_0p001_kb10000.csv"),
    ("01_stage1_instruments", "plink_tmfi_r2_0p001_kb10000.clumped"),
    ("01_stage1_instruments", "p2_tagging_result.txt"),
    ("06_stage6_finngen", "instruments_master.csv"),
    ("04_stage4_sampling_cov", "stage4_results/stage4_cross_trait_ldsc.json"),
    ("04_stage4_sampling_cov", "stage4_results/stage4_bias_correction.csv"),
    ("04_stage4_sampling_cov", "stage4_results/stage5_cross_trait_ldsc.json"),
]
for rel, fn in extra_frozen:
    src = ROOT / rel / fn
    dst = OUT / "02_frozen_outputs" / "_extra" / f"{rel.replace('/', '-')}__{fn}"
    copy_to(src, dst, frozen, skipped)

# ---------- 3. figure-generation code and figures ----------
figs = []
# superseded 2026-09-18 leave-one-out figure; replaced by the 2026-09-20 rerun Figure2_leave_one_out, so excluded
LEGACY_FIG_PREFIXES = ("Figure3_leave_one_out",)
for p in sorted(FIG_SRC.glob("*")):
    if p.name.startswith(LEGACY_FIG_PREFIXES):
        continue
    if p.suffix.lower() in (".svg", ".pdf", ".png"):
        copy_to(p, OUT / "03_figure_code_and_outputs" / "figures" / p.name, figs, skipped)

# ---------- 4. environment and log additions ----------
for src in [ROOT / "06_stage6_finngen" / "results" / "R_validation_sessionInfo.txt"]:
    copy_to(src, OUT / "05_environment_and_logs" / "06_stage6_finngen" / src.name, scripts, skipped)

# ---------- 4b. licence files ----------
MIT = """MIT License

Copyright (c) 2026 Cheng Chen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

CCBY = """Data and text licence — Creative Commons Attribution 4.0 International (CC BY 4.0)

Copyright (c) 2026 Cheng Chen

All non-code content in this repository — frozen intermediate outputs (CSV/TSV/JSON),
tables, figures and documentation — is licensed under CC BY 4.0.

You are free to share and adapt this material for any purpose, including commercially,
provided that you give appropriate credit to the author, cite the associated article,
and indicate whether changes were made.

Human-readable summary : https://creativecommons.org/licenses/by/4.0/
Full legal code        : https://creativecommons.org/licenses/by/4.0/legalcode

Scope note: the source GWAS summary statistics from which these derivatives were
computed (TMFI; chronic pain intensity; GIANT BMI; OpenGWAS/figshare site-matched
pain; FinnGen R12) are NOT included in this repository and remain governed by the
terms of their respective providers. This licence covers only the authors' own code
and derived outputs.
"""

(OUT / "LICENSE").write_bytes(MIT.replace("\n", "\r\n").encode("utf-8"))
(OUT / "LICENSE_DATA").write_bytes(CCBY.replace("\n", "\r\n").encode("utf-8"))

# ---------- 4b. .gitignore ----------
(OUT / ".gitignore").write_bytes(
    b"__pycache__/\r\n*.pyc\r\n.DS_Store\r\nThumbs.db\r\n"
)

# ---------- 5. MANIFEST ----------
rows = []
for p in sorted(OUT.rglob("*")):
    if ".git" in p.relative_to(OUT).parts or "__pycache__" in p.relative_to(OUT).parts:
        continue
    if p.is_file() and p.name != "MANIFEST.tsv":
        rel = p.relative_to(OUT).as_posix()
        rows.append((rel, p.stat().st_size, md5(p)))
manifest = ["path\tsize_bytes\tmd5"]
manifest += [f"{r[0]}\t{r[1]}\t{r[2]}" for r in rows]
(OUT / "MANIFEST.tsv").write_bytes(("\r\n".join(manifest) + "\r\n").encode("utf-8"))

# ---------- 6. environment.md ----------
env = """# Environment

| Component | Version | Purpose |
|---|---|---|
| Python | 3.13.12 | Main analysis pipeline (harmonisation, MR, power analysis, figures, recomputation of all manuscript numbers) |
| R | 4.6.1 | Validation against the official `MVMR` package (v0.4.8): IVW, Egger, conditional F |
| R (second implementation) | 4.5.1 | Independent base-R re-implementation (no third-party packages), used to cross-validate Stage 5 / Stage 6 |
| PLINK | 1.9.0-rc2 | Instrument LD clumping (r²<0.001, 10,000 kb window) and the LD-panel coverage probe |
| Python packages | numpy / pandas / scipy / matplotlib | Numerical computation and plotting |

The binaries in `tools/plink19/` and `reference/1000G_EUR.*` (1000 Genomes EUR LD reference panel, GRCh37, 1,836,406 variants / 503 individuals) are **not distributed with this repository**:
- PLINK is available from https://www.cog-genomics.org/plink/;
- the 1000 Genomes reference panel is available from https://www.internationalgenome.org/.

Raw GWAS summary statistics are likewise excluded (GB scale). See the Data availability statement in the main text and the per-stage sources below for how to obtain them.

## Data sources (public)

| Data | Source | Notes |
|---|---|---|
| TMFI | Zenodo 21017351 (concept 21017350) | UK Biobank thigh MRI V2, n=46,347 |
| Chronic pain intensity (CPI) | Zenodo 17219195 | n=134,627 / 106,697 (two definitions) |
| BMI | GIANT 2015 (Locke et al.) | file `SNP_gwas_mc_merge_nogc.tbl.uniq.gz` |
| Knee pain / back pain | OpenGWAS ukb-b-16254 / ukb-b-9838; figshare 29377565 | BOLT-LMM linear-model coefficients (risk-difference scale) |
| FinnGen R12 | FinnGen-only files from `summary_stats/release/` | PAIN / M13_LIMBPAIN / M13_FIBROMYALGIA |

## Recommended execution order

The per-stage reports (`Stage*_*.md`) are the authoritative documentation; the list below is compiled from script names and stage reports as a guide for reproduction:

1. **Stage 0 data audit** — `00_stage0_audit/audit_gwas_stage0.py`
2. **Stage 1 instruments** — `stage1_extract_match.py` → PLINK clump → `stage1_formal_postclump.py` → `stage1_formal_mr.R`
3. **Stage 2 power** — `stage2_power_analysis.py`, `stage2b_mvmr_power.py`
4. **Stage 3 multivariable MR** — `stage3_prepare_candidates.py` → PLINK clump (union) → `stage3_finalize_mvmr.py` → `stage3b_robustness.py`, `stage3b_external_exclusion_rerun.py` → `stage3_official_MVMR_validation.R`; see `p2_probe_step1..7.py` for the LD-panel coverage probe
5. **Stage 4 sample overlap** — `harmonise_sumstats.py` → `ldscore_from_reference.py` → `cross_trait_ldsc.py` → `mrlap_bias_correction.py`
6. **Stage 5 site-matched outcomes** — `opengwas_vcf_to_tsv.py` → `patch_exposures.py` / `patch_palindromic.py` → `stage5_mr.py` → `stage5_validation_baser.R` → `stage5_overlap_bias.py`
7. **Stage 6 FinnGen external triangulation** — `download_finngen.py` → `extract_instruments.py` → `harmonise_finngen.py` → `build_instrument_master.py` → `stage6_mr.py` → `stage6_official_MVMR_validation.R` / `stage6_validation_baser.R`
8. **Figures** — `make_manuscript_assets.py` (writes to `figures/`)
9. **Manuscript numbers** — `manuscript_numbers.py` (recomputes every number reported in the manuscript from the frozen outputs)
"""
(OUT / "05_environment_and_logs" / "environment.md").write_bytes(env.replace("\n", "\r\n").encode("utf-8"))

# ---------- 7. README ----------
def cnt(path):
    return len([p for p in (OUT / path).rglob("*") if p.is_file()])

readme = f"""# Code and data archive — TMFI / BMI / chronic pain intensity Mendelian randomization study

**Associated manuscript:** *Body mass index, but not thigh muscle fat infiltration, is associated with chronic pain intensity*
**Target journal:** Journal of Cachexia, Sarcopenia and Muscle (JCSM)
**Archived:** 2026-09-20
**Purpose:** Public versioned repository (GitHub + Zenodo DOI) accompanying the manuscript, corresponding to its Code availability and Data availability statements.

## Repository contents

| Directory | Description | Files |
|---|---|---|
| `01_analysis_scripts/` | All Stage 0–6 analysis scripts (.py / .R), preserving the original per-stage directory structure | {cnt('01_analysis_scripts')} |
| `02_frozen_outputs/` | Frozen intermediate outputs; filename prefixes correspond one-to-one to supplementary items (`S1__`–`S13__`, `S0-AppendixS1__`, `ADD-FinnGen__`). `_extra/` holds additional frozen inputs required for reproduction | {cnt('02_frozen_outputs')} |
| `03_figure_code_and_outputs/` | Figure-generation code `make_manuscript_assets.py` and figures (svg / pdf / png; 600 dpi tiff files are in the submission package `03_Figures/`). The superseded 2026-09-18 leave-one-out plot `Figure3_leave_one_out.*` is excluded; the leave-one-out figure is the 2026-09-20 rerun `Figure2_leave_one_out.*` (Supplementary Figure S1 in the submission package) | {cnt('03_figure_code_and_outputs')} |
| `04_manuscript_tooling/` | Recomputation of manuscript numbers (`manuscript_numbers.py`), word-count and format audits (`_final_audit.py`, `_jcsm_audit.py`), and manuscript revision / submission-package build scripts (`_apply_*.py`, `_build_*.py`, `_md_to_docx.py`) | {cnt('04_manuscript_tooling')} |
| `05_environment_and_logs/` | `environment.md` (versions, data sources, execution order), PLINK logs, R sessionInfo, per-stage summary JSON | {cnt('05_environment_and_logs')} |
| `MANIFEST.tsv` | Full file inventory: relative path / bytes / MD5 | — |

## Mapping to supplementary material

Each filename prefix in `02_frozen_outputs/` corresponds to an item in the supplementary-material list:

- `S0-AppendixS1__` → Supplementary Appendix S1 (Stage 0 data audit)
- `S1__` → Supplementary Table S1 (93-SNP joint instrument set)
- `S2__`–`S8__` → Supplementary Tables S2–S8
- `S9__` → Supplementary Table S9 (overlap audit for site-matched outcomes)
- `S10__` → Supplementary Table S10 (Python / base-R dual implementations)
- `S11__` → Supplementary Table S11 (two estimators of conditional F)
- `S12__` → Supplementary Table S12 (STROBE-MR checklist)
- `S13__` → Supplementary Table S13 (originality positioning)
- `ADD-FinnGen__` → FinnGen external triangulation results (underlying output of Table 6, included as additional material)

## Reproducibility

Every number reported in the manuscript can be deterministically recomputed by `manuscript_numbers.py` from the frozen outputs in this repository. Two independent implementations (the Python pipeline and base-R / the `MVMR` package) agree line-by-line at Stage 5 and Stage 6. `MANIFEST.tsv` provides MD5 checksums so that files can be verified as unmodified.

## Not included in this repository

- Raw GWAS summary statistics (GB scale; sources are listed in `05_environment_and_logs/environment.md`);
- PLINK binaries and the 1000 Genomes LD reference panel (third-party distributions);
- 600 dpi tiff figures (large; submitted with the manuscript package).

## Licence

- **Code** (`01_analysis_scripts/`, scripts in `03_figure_code_and_outputs/`, `04_manuscript_tooling/`): **MIT**, see `LICENSE`.
- **Data and text** (`02_frozen_outputs/`, figures, documentation): **CC BY 4.0**, see `LICENSE_DATA`.
- Copyright held by the author **Cheng Chen** (consistent with the manuscript title page).
"""
(OUT / "README.md").write_bytes(readme.replace("\n", "\r\n").encode("utf-8"))

print("Analysis scripts:", len(scripts), "| manuscript tooling:", len(tooling), "| frozen outputs:", len(frozen), "| figures:", len(figs))
print("Skipped/missing:", len(skipped))
for s in skipped[:20]:
    print("   ", s)
print("MANIFEST lines:", len(rows))
print("Output:", OUT)
