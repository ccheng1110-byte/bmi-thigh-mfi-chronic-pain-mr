# Code and data archive — TMFI / BMI / chronic pain intensity Mendelian randomization study

**Associated manuscript:** *Body mass index, but not thigh muscle fat infiltration, is associated with chronic pain intensity*
**Target journal:** Journal of Cachexia, Sarcopenia and Muscle (JCSM)
**Archived:** 2026-09-20
**Purpose:** Public versioned repository (GitHub + Zenodo DOI) accompanying the manuscript, corresponding to its Code availability and Data availability statements.

## Repository contents

| Directory | Description | Files |
|---|---|---|
| `01_analysis_scripts/` | All Stage 0–6 analysis scripts (.py / .R), preserving the original per-stage directory structure | 61 |
| `02_frozen_outputs/` | Frozen intermediate outputs; filename prefixes correspond one-to-one to supplementary items (`S1__`–`S13__`, `S0-AppendixS1__`, `ADD-FinnGen__`). `_extra/` holds additional frozen inputs required for reproduction | 56 |
| `03_figure_code_and_outputs/` | Figure-generation code `make_manuscript_assets.py` and figures (svg / pdf / png; 600 dpi tiff files are in the submission package `03_Figures/`). The superseded 2026-09-18 leave-one-out plot `Figure3_leave_one_out.*` is excluded; the leave-one-out figure is the 2026-09-20 rerun `Figure2_leave_one_out.*` (Supplementary Figure S1 in the submission package) | 13 |
| `04_manuscript_tooling/` | Recomputation of manuscript numbers (`manuscript_numbers.py`), word-count and format audits (`_final_audit.py`, `_jcsm_audit.py`), and manuscript revision / submission-package build scripts (`_apply_*.py`, `_build_*.py`, `_md_to_docx.py`) | 40 |
| `05_environment_and_logs/` | `environment.md` (versions, data sources, execution order), PLINK logs, R sessionInfo, per-stage summary JSON | 19 |
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
