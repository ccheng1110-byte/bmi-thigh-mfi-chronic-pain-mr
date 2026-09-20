# Environment

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
