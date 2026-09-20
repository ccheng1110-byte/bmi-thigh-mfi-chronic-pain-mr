# STROBE-MR checklist (remapped for the JCSM-adapted manuscript)

**Manuscript:** *Body mass index, but not thigh muscle fat infiltration, is associated with chronic pain intensity*
**Target journal:** Journal of Cachexia, Sarcopenia and Muscle (JCSM)

**Checklist source:** Skrivankova VW, Richmond RC, Woolf BAR, et al. Strengthening the Reporting of Observational Studies in Epidemiology Using Mendelian Randomization: The STROBE-MR Statement. *JAMA*. 2021;326(16):1614–1621. doi:10.1001/jama.2021.18236

**Line numbers refer to** `Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md` **as of 2026-09-20** (424 lines, JCSM-compressed version). The previous version of this checklist pointed at the pre-compression 531-line draft and is void. Every line reference below was re-derived from the current file; re-derive again if the manuscript is edited substantially.

---

## Title and abstract

| Item | Checklist item | Addressed at | Note |
|---|---|---|---|
| 1 | Indicate MR as the study's design in the title and/or the abstract | Title (L1); Abstract Methods (L27); Keywords (L37) | Title is a contrast statement; the design is named in the abstract Methods ("Mendelian randomization estimated…") and in two keywords |

## Introduction

| Item | Checklist item | Addressed at | Note |
|---|---|---|---|
| 2 | Background: rationale; what is the exposure; is a causal relationship plausible; why MR is helpful | L39–41 | Exposure = thigh muscle fat infiltration (TMFI); confounded observational literature and reverse causation motivate MR |
| 3 | Objectives, including pre-specified causal hypotheses; state that MR estimates causal effects only under assumptions | L41; secondary analyses L64, L66; preregistration status L68 | Hypothesis: TMFI has a direct effect on chronic pain intensity independent of BMI; BMI is the internal comparator. Study explicitly **not preregistered**; analysis choices fixed before Stage 2 estimates were inspected |

## Methods

| Item | Checklist item | Addressed at | Note |
|---|---|---|---|
| 4 | Study design and data sources | L44–47; Table 1 (L105–119) | Two analysis stages on the same summary data (a stage is an analysis step, not an independent sample); Table 1 lists every source with ancestry and per-SNP N |
| 4a | Setting, locations, relevant dates | L46; Table 1 | UK Biobank (TMFI, both CPI definitions, knee/back pain), GIANT 2015 (BMI), FinnGen R12 (three endpoints), figshare second knee-pain GWAS |
| 4b | Participants: eligibility, selection, sample size, prior power calculations | L46; L56; Table 1; Table 4 (L166–177) | Power reported as minimum detectable effect (MDE) rather than post hoc; MDE solved exactly from the observed SE |
| 4c | Measurement, quality control and selection of genetic variants | L48 | Genome-wide selection, LD clumping (r²<0.001 within 10 Mb, PLINK 1.9.0-rc2), allele-frequency strand alignment; palindromic instruments recovered rather than silently dropped |
| 4d | Methods of assessment and diagnostic criteria for each exposure, outcome and covariate | L46; L64; Supplementary Methods L277–279 | Phenotype identity verified independently of repository metadata (GDF5 for knee, SOX5-adjacent for back; allele-frequency agreement r≈0.998) |
| 4e | Ethics committee approval and participant consent | L376 (Ethics statement); L384 (Declarations) | Summary data only; underlying UKB approvals cited [14,15]; FinnGen obtained under its data-access policy |
| 5 | Assumptions: the three core IV assumptions plus assumptions of any additional or sensitivity analysis | L44–45; relevance L58; independence L62; exclusion restriction L54, L60 | Relevance by conditional F; independence by cross-trait LDSC overlap audit; exclusion restriction by MVMR-Egger, weighted median, Huber-IRLS and targeted exclusion |
| 6 | Statistical methods: main analysis | L50–53 | Univariable IVW-MRE (Stage 1); two-exposure no-intercept MVMR with MRE IVW weighting (Stage 2) |
| 6a | How quantitative variables were handled (scale, units, model) | L56; Supplementary Methods L253–256 (BMI scale), L259–262 (site-matched scale) | SD-standardised inference throughout; implied-SD estimator with three independent checks for BMI; BOLT-LMM risk differences converted by log OR = β/[μ(1−μ)] |
| 6b | How genetic variants were handled and how weights were selected | L48–49; L52 | Inverse-variance weighting; identical frozen instrument sets carried across every stage and outcome |
| 6c | The MR estimator and related statistics; covariates; whether the same covariate set was used in both samples | L50–55 | Summary-level data, so no covariate adjustment is possible in either sample; the joint model includes both exposures simultaneously |
| 6d | How missing data were addressed | L68 | Complete-case harmonisation: a variant enters an analysis only if exposure, outcome and LD panel all supply a usable allele pair, frequency, effect and SE; no imputation |
| 6e | How multiple testing was addressed | L68 | No correction applied; the primary contrast is a single comparison fixed in advance and the site-matched and FinnGen analyses are declared secondary and exploratory |
| 7 | Assessment of assumptions: methods or prior knowledge used | L58 (conditional F, Q, Cook); L50 (Steiger); L62 (cross-trait LDSC, structural ceiling); Supplementary Methods L265–270 (independent re-implementation) | — |
| 8 | Sensitivity analyses and additional analyses | L60 (leave-one-out, Huber-IRLS, Cook, targeted exclusion); L64 (site-matched); L66 (FinnGen) | — |
| 9a | Software and packages, with versions and settings | L68; L372 (Code availability) | R 4.6.1 with MVMR 0.4.8; PLINK 1.9.0-rc2; Python 3.13 with NumPy, pandas, SciPy; base-R re-implementation (R 4.5.1, no third-party packages) |
| 9b | Whether the protocol was preregistered, and when and where | L68 | **Not preregistered** — stated explicitly; instrument sets, outcome definitions, estimator hierarchy and sensitivity analyses were fixed before the Stage 2 estimates were inspected, and the frozen files are deposited |

## Results

| Item | Checklist item | Addressed at | Note |
|---|---|---|---|
| 10a | Numbers of individuals at each stage and reasons for exclusion; consider a flow diagram | L48; Table 1 (L105); Figure 1 (L215) | Candidates → 42 Stage 1 instruments → 93 joint instruments (59 BMI-only, 32 TMFI-only, 2 shared); FinnGen retains 108/109 |
| 10b | Summary statistics for exposures, outcomes and other relevant variables | Table 1; Table 4 (L166) | Per-SNP N ranges, ancestry, and the SDs used for standardisation |
| 10c | If meta-analyses contributed, assess heterogeneity across studies | Supplementary Methods L253; Limitations L99 | GIANT BMI is a meta-analysis; its genomic-control correction (λ=1.61) is documented and propagated |
| 10d(i) | Two-sample MR: justify similarity of variant–exposure associations across samples | L46; Table 1 | All samples European-ancestry; a single ancestry panel throughout (FinnGen is Finnish, disclosed) |
| 10d(ii) | Two-sample MR: information on the number of overlapping individuals | L62; L80–81; Supplementary Table S16 (L318–354) | Head-count not identifiable; the product N_shared·ρ_phen is bounded against its structural ceiling and reported as such |
| 11a | Variant–exposure and variant–outcome associations, on an interpretable scale | Supplementary Table S1 (listed L225) | Per-variant exposure and outcome coefficients for the full 93-SNP set |
| 11b | MR estimates with uncertainty, on an interpretable scale (e.g. OR per SD) | L71–75; Table 2 (L121), Table 3 (L142), Table 5 (L179) | SD of pain intensity per SD of exposure in the primary analysis; OR per SD for the site-matched and FinnGen outcomes |
| 11c | Translate relative measures into absolute risk where relevant | Table 5 note (L181) | Risk differences per SD reported alongside the odds ratios: +0.24 pp (TMFI, knee) vs +4.08 pp (BMI, knee) |
| 11d | Plots to visualise results | Figure 2 (L219); Supplementary Figure S1 (L359), S2 (L363) | Forest-style primary estimates, leave-one-out series, influence diagnostics |
| 12a | Report the assessment of the validity of the assumptions | L73; Table 3 (L142); Supplementary Table S11 (L238) | Conditional F 22.12 (TMFI) and 35.23 (BMI); Steiger 42/42; closed-form alternative 23.6 / 42.8 |
| 12b | Additional statistics: heterogeneity across variants (I², Q, E-value) | L73; L84; Table 3 | Modified and unmodified Q, φ dispersion, MVMR-Egger intercept |
| 13a | Sensitivity analyses assessing robustness to assumption violations | L78–79; Table 2 (L121); Supplementary Table S15 (L299), S16 (L318) | — |
| 13b | Other sensitivity or additional analyses | L60; L64; L82 | Targeted metabolic-pleiotropy exclusion; site-matched secondary outcomes |
| 13c | Any assessment of the direction of the causal relationship | L50; L71; Table 3 | Steiger directionality supports the exposure orientation for 42/42 Stage 1 instruments |
| 13d | Compare with estimates from non-MR analyses where relevant | L90–91 | BMI estimates agree in direction and magnitude with the adiposity-to-pain MR literature [19,20] |
| 13e | Additional plots (e.g. leave-one-out) | Supplementary Figure S1 (L359) | — |

## Discussion

| Item | Checklist item | Addressed at | Note |
|---|---|---|---|
| 14 | Key results with reference to the objectives | L87–89; Conclusion (L100) | — |
| 15 | Limitations, including validity of the IV assumptions, other sources of bias, imprecision; direction and magnitude of potential bias; efforts to address it | L98–99 (eight boundaries in three groups) | Source data (TMFI flag, LD panel sparsity, GIANT correction); model validity (residual heterogeneity, non-identifiable overlap head-count, BOLT-LMM approximation); external validity (reversed BMI coefficient in the second knee-pain specification, outcome-side-only replication, generalizability) |
| 16a | Meaning: cautious overall interpretation against limitations and other studies | L87–89; L90–91 | Null stated as bounded, not as an absent effect |
| 16b | Mechanism and whether the gene–environment equivalence assumption is reasonable | L92–93 | Causal language confined to the stated assumptions; the estimates describe lifelong genetically predicted differences, not interventions |
| 16c | Clinical or policy relevance, and what the results imply about intervention effect sizes | L96–97 | TMFI is not supported as an intervention target for chronic pain while adiposity itself is unaddressed |
| 17 | Generalizability: to other populations, other exposure periods/durations, other exposure levels | L99 | European ancestry only; lifetime-averaged genetic difference rather than a timed intervention; full genetic range rather than a clinically attainable change |

## Other information

| Item | Checklist item | Addressed at | Note |
|---|---|---|---|
| 18 | Funding and the role of funders, including funding of the source databases and original studies | L392 (Declarations, Funding) | *To be completed by the authors* |
| 19 | Data and data sharing; statistical code needed to reproduce the results, and where it can be obtained | L368 (Data availability); L372 (Code availability); L388 (Declarations) | All source GWAS publicly available; frozen CSV outputs and scripts to be deposited with a permanent DOI at submission — **the DOI must replace the three "will be deposited" placeholders before submission** |
| 20 | Conflicts of interest | L390 (Declarations) | *To be completed by the authors* |

---

## Two open items that are not author-administrative

1. **Item 19 is not yet fully satisfied.** Three places still say the code archive "will be deposited … at the time of submission" (L372, L388, and the Supplementary material preamble at L223). The repository must be public and the permanent DOI filled in before submission; JCSM's editorial office checks data and code availability statements.

2. **Item 15 carries one unresolved scientific caveat.** The TMFI source GWAS header records `--fastGWA-mlm-binary` (a logistic mixed model, Jiang et al. *Nat Genet* 2021;53:1616–1621) whereas the source publication describes a continuous-trait analysis. The manuscript now states the standardized bounds as conditional on the exposure scale (Limitations, L99; Supplementary Methods L273–274) rather than asserting that the conclusions are independent of the flag. Confirmation from the source authors would close this item.

## Changes since the 2026-09-19 version

- All line references remapped (the manuscript was compressed from 566 to 424 lines for JCSM and most derivations moved to Supplementary Methods).
- Item 1: the design name was dropped from the title during compression; it is now restored in the abstract Methods so the item is satisfied by the abstract rather than the title.
- Item 3 / 6e: "pre-specified" replaced throughout by "fixed in advance / fixed before the Stage 2 estimates were inspected", because the study is declared not preregistered; the site-matched and FinnGen analyses are now labelled secondary and exploratory.
- Item 15: limitations regrouped into three blocks (source data / model validity / external validity) and the claim that inference does not depend on re-clumping with a denser LD panel was removed.
