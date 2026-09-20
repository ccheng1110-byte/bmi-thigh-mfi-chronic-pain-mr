# Stage 0 GWAS data audit report (final)

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate (data integrity and reproducibility verification)
- Origin Date: 2026-09-18T13:35:00+08:00
- Verification Status: VERIFIED
- Version Label: stage0_gwas_audit_v2
- Upstream Dependency: `研究方案_TMFI_慢性疼痛MR_v2.0_2026-09-18.md` (study protocol; original filename retained verbatim)

## 1. Final verdict

**Stage 0: PASS WITH WARNINGS (proceed, with warnings that must be handled during preprocessing).**

All three target GWAS files are fully on disk; file sizes and MD5 checksums match the official Zenodo metadata exactly. Both gzip files decompress completely, and all three files were parsed stream-wise from first to last line with no truncation, download error pages or parse failures. P values, standard errors, allele frequencies and effect sizes in the core numeric fields show no out-of-range or illegal values.

The data can proceed to the next stage, "format standardisation and instrument extraction", but the three raw files must NOT be fed into MR directly: the TMFI duplicate loci, the CPI no-imaging file's special PLINK2 field mapping, chromosome X, and the CPI genome build not explicitly declared by the repository must all be handled first.

## 2. File identity and integrity

| Dataset | File | Bytes | Zenodo expected MD5 | Local MD5 | Verdict |
|---|---|---:|---|---|---|
| TMFI | `TMFI_GWAS.csv` | 948,739,649 | `b5ac4d2d1cb7269ca318c328d1ec0d5a` | `b5ac4d2d1cb7269ca318c328d1ec0d5a` | PASS |
| CPI full sample | `Chronic_pain_intensity.txt.gz` | 412,658,977 | `866f968c7eace91d281d9039f860a4b3` | `866f968c7eace91d281d9039f860a4b3` | PASS |
| CPI no-brain-imaging sample | `chronic_pain_intensity_no_imaging.txt.gz` | 420,494,187 | `47d3a1cb749c7d4c36817ded17bed1cd` | `47d3a1cb749c7d4c36817ded17bed1cd` | PASS |

Data sources:

- TMFI concept DOI: [10.5281/zenodo.21017350](https://doi.org/10.5281/zenodo.21017350); version DOI locked for this analysis: [10.5281/zenodo.21017351](https://doi.org/10.5281/zenodo.21017351).
- CPI version DOI: [10.5281/zenodo.17219195](https://doi.org/10.5281/zenodo.17219195).
- All three datasets are open access; the Zenodo metadata licence is CC BY 4.0.

## 3. Whole-file scan results

| Metric | TMFI | CPI full sample | CPI no-imaging sample |
|---|---:|---:|---:|
| Data rows (excluding header) | 12,128,084 | 12,971,234 | 12,685,460 |
| Chromosomes | 1–22 | 1–22, X | 1–22 |
| Min / max position | 302 / 249,230,914 | 302 / 249,230,914 | 302 / 249,230,914 |
| rsID rows | 11,659,539 | 12,386,849 | 12,114,898 |
| Non-rsID rows | 468,545 | 584,385 | 570,562 |
| Excess rows from duplicate IDs* | 11,821 | 0 | 0 |
| Excess rows from duplicate CHR:POS | 14,455 | 4,925 | 4,973 |
| Missing core numeric values | 0 | 0 | 0 |
| Illegal allele rows | 0 | 0 | 0 |
| Rows with A1 identical to the other allele | 0 | 0 | 0 |
| P out of range | 0 | 0 | 0 |
| SE ≤ 0 | 0 | 0 | 0 |
| Frequency out of range | 0 | 0 | 0 |
| Non-finite BETA | 0 | 0 | 0 |
| Full-stream parse | PASS | PASS | PASS |

\* ID de-duplication uses a deterministic 64-bit pandas hash; collision risk is negligible, but exact de-duplication by original string is still required when formally extracting instruments.

## 4. Fields and value ranges

### 4.1 TMFI

- Actual fields: `CHR, SNP, POS, A1, A2, N, AF1, BETA, SE, P`
- Rows: 12,128,084.
- P: `6.94489×10⁻79` to `1`.
- SE: `0.0132345` to `0.267076`.
- AF1: `0.000626147` to `0.999244`; this is the A1 frequency, not a MAF restricted to 0–0.5.
- BETA: `-0.907893` to `1.08565`.
- Per-SNP N: `43,880` to `46,347`.
- The paper reports European ancestry, UK Biobank V2 MRI, GRCh37, GCTA fastGWA-mlm, maximum sample size 46,347.

### 4.2 CPI full sample

- Actual fields: `CHR, BP, SNP, A1, A2, MAF, BETA, SE, Z, P, N`
- Rows: 12,971,234.
- P: `1.75506×10⁻9` to `1`.
- SE: `0.00926175` to `0.177022`.
- MAF column: `0.001` to `0.500863`. A small number of values slightly exceed 0.5, suggesting the column may be affected by floating-point or frequency-definition issues; MAF is subsequently computed uniformly as `min(f, 1-f)`.
- BETA: `-0.738241` to `0.835478`.
- Per-SNP N: `70,138` to `134,627`.
- Maximum absolute difference between `BETA/SE` and Z is `2.98×10⁻5`, consistent with rounding differences from output precision.
- 296,094 rows in the CHR column are `X`; the "non-numeric CHR" values seen during numeric conversion are all caused by legitimate chromosome X, not by missingness.

### 4.3 CPI no-brain-imaging sample

- Actual fields: `#CHROM, POS, ID, REF, ALT, PROVISIONAL_REF?, A1, OMITTED, A1_FREQ, TEST, OBS_CT, BETA, SE, T_STAT, P, ERRCODE`
- Rows: 12,685,460.
- P: `2.08442×10⁻9` to `1`.
- SE: `0.0125917` to `0.202278`.
- A1_FREQ: `0.001` to `0.5`.
- BETA: `-0.767891` to `0.961409`.
- Per-SNP OBS_CT: `55,638` to `106,697`.
- All `TEST` values are `ADD`; all `ERRCODE` values are `.`; no PLINK2 model-error rows were found.
- Maximum absolute difference between `BETA/SE` and T_STAT is `2.72×10⁻5`, consistent with rounding error.

## 5. Warnings that must be addressed

### W1. CPI no-imaging file does not match the README fields

The official README summarises CPI as `CHR/BP/SNP/A1/A2/MAF/BETA/SE/P/N`, but the no-imaging file actually retains the raw 16-column PLINK2 format. The following mapping should be used:

| Standard field | No-imaging file field |
|---|---|
| CHR | `#CHROM` |
| BP | `POS` |
| SNP | `ID` |
| effect allele | `A1` |
| other allele | `OMITTED` |
| effect allele frequency | `A1_FREQ` |
| N | `OBS_CT` |
| BETA / SE / P | same field names |

`REF/ALT` are used for coordinate and reference-allele checks, but the MR effect direction should be taken from `A1` and `OMITTED`.

### W2. TMFI contains duplicate IDs and duplicate coordinates

TMFI has 11,821 excess rows from duplicate IDs and 14,455 excess rows from duplicate coordinates. Duplicate coordinates also occur in both CPI files and typically arise from multi-allelic variants or different variant representations. The next stage must normalise by `CHR:POS:EA:OA` and perform exact de-duplication for instrument SNPs; simply retaining an arbitrary first row is not acceptable.

### W3. CPI genome build is not explicitly stated in the current Zenodo README

TMFI is explicitly GRCh37. The CPI coordinate range and large numbers of UK Biobank variant representations are highly consistent with TMFI, supporting compatibility with GRCh37, but this is an inference at the data level and does not amount to an explicit statement by the repository. Before formal harmonisation it should be re-confirmed via the original paper, author notes or reference variant coordinates; no undocumented liftOver may be performed before confirmation.

### W4. Sample overlap is still not excluded

TMFI comes from UK Biobank thigh MRI V2; the CPI no-imaging version only states that participants with brain imaging data were excluded. Brain imaging and thigh MRI are not the same definition, so the no-imaging CPI cannot be described as having "zero sample overlap" with TMFI. It is recommended to use the no-imaging CPI as the preferred primary-analysis candidate and the full-sample CPI as a sensitivity analysis, and to report potential overlap conservatively in the paper.

### W5. Palindromic SNPs require frequency-assisted handling

The three files contain 1,716,071, 1,808,929 and 1,768,670 rows of A/T or C/G palindromic SNPs respectively; of these, rows near intermediate frequency (MAF ≥ 0.42) number 110,027, 123,155 and 119,676 respectively. Palindromic instrument SNPs that cannot be reliably oriented should be removed during harmonisation.

## 6. Stage 0 gate verdict

| Gate item | Verdict |
|---|---|
| Source, version, licence | PASS |
| File on disk and MD5 | PASS |
| gzip / text integrity | PASS |
| Whole-file structural parse | PASS |
| Core field completeness | PASS WITH WARNING (no-imaging file needs a dedicated mapping) |
| Numeric validity | PASS |
| SNP uniqueness | PASS WITH WARNING (normalised de-duplication required) |
| Genome build | PASS WITH WARNING (CPI to be finally confirmed) |
| Sample overlap risk | WARNING |
| Proceed to next stage | **Permitted: format standardisation, exact de-duplication and instrument extraction only** |

## 7. Reproduction materials

- `audit_gwas_stage0.py`: chunked, whole-file streaming audit script.
- `stage0_scan_results.json`: all scan results in machine-readable form.
- `metadata/zenodo_21017350_resolved_metadata.json`: TMFI official metadata snapshot.
- `metadata/zenodo_17219195_metadata.json`: CPI official metadata snapshot.
- `metadata/CPI_README.txt`: CPI official field description; MD5 verified against Zenodo.

## 8. Recommended next step

Create a standardised intermediate file containing only the fields required for MR, fixed as:

`CHR, BP, SNP, EA, OA, EAF, BETA, SE, P, N, SOURCE`

Then extract candidate instrument SNPs from TMFI at `P < 5×10⁻8`, perform exact de-duplication, LD clumping, F-statistic filtering and outcome-data matching. Do not run the formal MR main model until instrument quality control is complete.
