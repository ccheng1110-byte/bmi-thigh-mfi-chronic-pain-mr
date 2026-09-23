# -*- coding: utf-8 -*-
"""JCSM adaptation pass 2: Plan A + title candidate 1, trimmed to <=4500 body words."""
import re, sys

PATH = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"

with open(PATH, "rb") as f:
    raw = f.read()
text = raw.decode("utf-8")
NL = "\r\n"

def rep(old, new, n=1, what=""):
    global text
    c = text.count(old)
    assert c == n, "FAIL %s: expected %d, found %d" % (what or repr(old[:60]), n, c)
    text = text.replace(old, new)

# ---------------- 1. Title ----------------
OLD_TITLE = ("# Body mass index, but not thigh muscle fat infiltration, shows a direct genetic "
             "association with chronic pain intensity within the resolution of this design: a "
             "multivariable Mendelian randomization study with explicit power bounds")
NEW_TITLE = "# Body mass index, but not thigh muscle fat infiltration, is associated with chronic pain intensity"
rep(OLD_TITLE, NEW_TITLE, 1, "title")

# ---------------- 2. Keywords ----------------
rep("**Keywords:** chronic pain intensity; thigh muscle fat infiltration; body mass index; "
    "Mendelian randomization; multivariable Mendelian randomization; UK Biobank; pleiotropy",
    "**Keywords:** chronic pain intensity; thigh muscle fat infiltration; body mass index; "
    "Mendelian randomization; multivariable Mendelian randomization; UK Biobank", 1, "keywords")

# ---------------- 3. New body ----------------
BODY = NL.join([
"## Introduction",

"Chronic pain is a leading cause of disability; large-scale genetic studies have begun to define the architecture of pain phenotypes, including chronic pain intensity (CPI) in the UK Biobank [3]. Body composition is central to this literature because systemic adiposity may influence pain through mechanical loading, inflammation, metabolic dysfunction and altered pain processing, and genetically predicted adiposity is already an established cause of site-specific pain [19,20].",

"Muscle phenotypes have attracted parallel attention: genetically predicted lower-extremity muscle traits have been related to spinal diseases [7], sarcopenia has been reported to share genetic determinants with chronic pain [6], and the genetic and multi-omic architecture of thigh muscle fat infiltration (TMFI) \u2014 an MRI-derived measure of muscle quality \u2014 has recently been characterized [1], although with infiltration treated as the outcome rather than as an exposure on pain [21]. What neither literature provides is the joint test: a muscle-compartment exposure and an adiposity exposure estimated in one model, so that a muscle effect can be read against a bound rather than against a P value. Existing muscle\u2013pain genetic analyses share two weaknesses \u2014 they rarely include a comparator exposure able to demonstrate that the model detects an effect when one exists, and they rarely state the smallest effect the design could have detected \u2014 so a null reported without either is descriptively correct but interpretively uninformative.",

"We addressed both weaknesses in two analysis stages applied to the same summary data (a stage is an analysis step, not an independent sample; the design is not two-stage least squares). Stage 1 estimated the total effect of TMFI on CPI. Stage 2 fitted two-exposure multivariable Mendelian randomization (MVMR), estimating each exposure's effect conditional on the other [4,5], treating the BMI coefficient as an internal comparator and reporting study power as the minimum detectable effect on a standardized scale [8]. Because Stage 1 already showed no detectable TMFI total effect, the question was not whether a BMI signal could be separated from a TMFI signal, but how large a TMFI effect these data can exclude. We examined two CPI definitions, quantified exposure\u2013outcome sample overlap genetically, analysed two site-matched outcomes (knee and back pain), and replicated both stages in FinnGen, a cohort sharing no participants with UK Biobank.",

"## Methods",

"### Study design and estimand",

"Stage 1 estimated the univariable (total) effect of TMFI on CPI. Stage 2 estimated, for each exposure, its direct effect on CPI conditional on the other exposure using two-exposure summary-data MVMR. Direct effects are narrower than total effects: each coefficient represents the component of that exposure's genetic contrast not explained by the other. Because Stage 1 returned a null total effect, Stage 2 was not framed as separating one signal from another; BMI instead served as an internal comparator, demonstrating that the instrument set and outcome data can detect an effect of the size expected for a body-composition trait and permitting the null TMFI estimate to be read against the study's minimum detectable effect. The primary estimator was multiplicative-random-effects (MRE) inverse-variance-weighted (IVW) MVMR; robustness analyses tested dependence on single instruments, high-leverage observations and annotated metabolic-pleiotropy candidates.",

"### GWAS sources and phenotype definitions",

"Summary-level data were assembled for one imaging-derived exposure (TMFI, UK Biobank thigh MRI GWAS, n=46,347 [1]), one systemic adiposity exposure (BMI, GIANT 2015 [2]) and two UK Biobank CPI outcomes (n=134,627 and 106,697 [3]); sources, ancestries and effective sample sizes are given in Table 1. UK Biobank is a prospective cohort of over 500,000 individuals [14], with central genotyping, quality control and imputation [15]. All four source GWAS were restricted to European ancestry; the LD reference panel was 1000 Genomes phase 3 European on GRCh37 (503 unrelated individuals). The two exposure GWAS are treated as non-overlapping \u2014 verified rather than assumed, with a cross-trait LD score intercept of \u22120.018 (SE 0.009). Because the TMFI source reports no verifiable physical unit, all primary inference is standardized (SD of outcome per SD of exposure), with raw-scale coefficients reported alongside.",

"### Instrument selection and harmonization",

"Variants associated with either exposure at genome-wide significance (P<5\u00d710\u207b\u2078) were harmonized before clumping; incompatible and frequency-ambiguous palindromic variants were removed. LD clumping used PLINK v1.9.0-rc2 (r\u00b2<0.001 within 10,000 kb). Stage 1 applied this to the TMFI candidates alone, giving 42 instruments; the reference panel contained only 22.75% of TMFI candidates, so eligible signals outside it may have been missed \u2014 a disclosed limitation, since every retained variant was present. Stage 2 clumped the union of TMFI- and BMI-significant candidates, giving 93 variants (59 BMI-only, 32 TMFI-only, 2 shared); union-based clumping changes index SNPs within blocks, so the joint set is not a superset of the Stage 1 set (26 of 42 retained). The same 93-SNP definition was carried forward for both CPI outcomes. Five palindromic instruments were oriented on allele frequency against the TMFI anchor; the alternative orientation leaves both conclusions unchanged (per-variant decisions in Supplementary Table S1).",

"### Univariable MR stage (Stage 1)",

"For each CPI outcome we estimated the total TMFI effect with the 42 instruments: IVW-MRE as primary estimator, weighted median and MR-Egger as pleiotropy-robust alternatives, Steiger directionality testing [9], and heterogeneity quantified by Cochran's Q and the dispersion statistic \u03c6.",

"### Primary multivariable model (Stage 2)",

"We fitted the no-intercept joint model \u0393\u0302<sub>j,CPI</sub> = \u03b8<sub>TMFI</sub>\u03b3\u0302<sub>j,TMFI</sub> + \u03b8<sub>BMI</sub>\u03b3\u0302<sub>j,BMI</sub> + \u03b5<sub>j</sub>, whose coefficients are direct effects conditional on the other exposure [4,5], with MRE IVW weighting and standard errors from the R MVMR validation workflow.",

"### Pleiotropy-robust and alternative estimators",

"Fixed-effect IVW showed the effect of removing the heterogeneity correction; MVMR-Egger's intercept estimates the average direct pleiotropic effect. Influence diagnostics are described under Robustness analyses.",

"### Statistical power and minimum detectable effect",

"A null estimate is informative only if the design could have detected a relevant effect. We report the minimum detectable effect (MDE): the smallest absolute standardized effect detectable at two-sided \u03b1=0.05 with 80% (MDE<sub>80</sub>) or 90% (MDE<sub>90</sub>) power, solved exactly from the two-sided normal power function for the observed standard error [8]. Standardization used implied phenotypic standard deviations recovered as SE \u00d7 \u221a(2 \u00d7 N \u00d7 EAF \u00d7 (1 \u2212 EAF)): sd<sub>TMFI</sub>=2.03, sd<sub>CPI_full</sub>=2.89, sd<sub>CPI_no_imaging</sub>=2.93. For BMI, three independent checks \u2014 published per-allele effects at TMEM18, MC4R and SEC16B, variance explained across 71 independent loci, and the estimator's internal consistency for the other three sources \u2014 establish that the GIANT coefficients are already per SD of BMI, the implied-SD surplus reflecting genomic-control-corrected standard errors (single-trait intercept 0.62), so sd<sub>BMI</sub>=1.00 (derivation in Supplementary Methods). Effects are expressed as SD of pain intensity per SD of exposure.",

"### Strength, heterogeneity and influence diagnostics",

"Conditional instrument strength was assessed with multivariable conditional F statistics from the R `MVMR` package [5], the values reported throughout (Tables 3 and 5); a closed-form approximation computed independently returns 23.6 for TMFI and 42.8 for BMI against 22.12 and 35.23 from `MVMR` (Supplementary Table S11), both far above 10. The modified Q statistic quantified residual heterogeneity and pleiotropic misspecification [4,5]. Because the released BMI standard errors are inflated, the conditional F for BMI is conservative by a factor of about 1.6.",

"### Robustness analyses",

"Four layers were evaluated: leave-one-out refitting across all 93 variants; Cook-distance exclusions above the prespecified 4/n threshold and of the top one or three variants; exploratory Huber-IRLS reweighting; and targeted exclusion of the two highest-priority metabolic-pleiotropy candidates (rs731839, rs13107325) and, separately, all five review targets. Overlap between the Cook-threshold and review-target sets is reported explicitly and not counted as independent replication.",

"### Quantifying sample overlap by cross-trait LD score regression",

"In two-sample MR, overlap between the UK Biobank-derived exposure and outcome samples biases estimates toward the confounded observational association; we quantified it rather than asserting it. LD scores were computed from the 1000 Genomes European panel, and single-trait and cross-trait LD score regression [10] was fitted by two-step weighted least squares with a 200-block delete-one jackknife; under this model the cross-trait intercept equals N_shared\u00b7\u03c1_phen/\u221a(N_kN_m), exactly the sampling covariance that generates the MR bias. Two checks validated the pipeline: all three UK Biobank single-trait intercepts fell within 2% of 1 (0.983, 1.000, 1.010), and the nested-subset check \u2014 CPI_no_imaging is a subset of CPI_full and therefore maximally overlapping by design \u2014 returned an intercept of 0.910 (SE 0.010) at its structural ceiling of 0.890 within ~2 SE, with genetic correlation 1.0025. The intercept propagates to the estimates analytically as bias(\u03b8\u0302) = (X'WX)\u207b\u00b9u with u_i = intercept_i\u00b7\u03a3_j(se_Xij/se_Yj) \u2014 the analytic core of MRlap [11] re-implemented against the exact instrument sets. Because the bias enters through the noise covariance while the signal enters through \u03b2_X\u00b2, strong instruments dilute it, so a significant intercept need not imply consequential bias; the identified quantity is the product N_shared\u00b7\u03c1_phen, not a participant head-count. The GIANT release is genomic-control corrected (single-trait intercept 0.619), so intercepts involving BMI were rescaled by 1/\u221a(c_k c_m) and stress-tested over assumed non-negative values.",

"### Pre-specified site-matched secondary analysis",

"Because chronic pain intensity does not distinguish pain location, a site-specific effect could be diluted across locations; we therefore added two binary site-matched outcomes analysed with the identical frozen instrument sets: knee pain and back pain, self-reported for the previous month in the same UK Biobank cohort (461,857 participants; 98,704 and 118,471 cases; MRC IEU OpenGWAS ukb-b-16254 and ukb-b-9838 [12], fitted with BOLT-LMM, a linear mixed model, to the binary phenotypes). A second knee-pain GWAS from the same cohort (429,737 participants [13]) served as a same-cohort concordance check on the univariable TMFI estimate only. Phenotype identities were verified independently of repository metadata \u2014 the knee source leads at rs143384 in GDF5 (\u2212log\u2081\u2080P=21.0), absent from back pain, which leads near SOX5, and allele frequencies agree with the reference panel at r\u22480.998 (full checks in Supplementary Methods). Because BOLT-LMM applies a linear model to a 0/1 phenotype, the distributed coefficients are risk differences rather than log odds ratios; this was verified on the summary statistics themselves, and odds ratios were obtained by log OR = \u03b2/[\u03bc(1\u2212\u03bc)] with the standard error divided by the same constant \u2014 every z statistic, P value, F statistic and tipping point is exactly invariant to the conversion, which is an approximation only in the odds-ratio magnitudes (verification and derivation in Supplementary Methods).",

"### External replication in a cohort with no participant overlap",

"To test the two central contrasts where overlap bias cannot arise, we repeated both stages against three endpoints from the FinnGen R12 release [17] using FinnGen-only summary statistics; the FinnGen \u00d7 UK Biobank meta-analysis module was not used, to avoid reintroducing the overlap. Endpoints were chosen a priori: `PAIN` (237,944 cases, 261,418 controls), `M13_LIMBPAIN` (45,476 / 357,549) and `M13_FIBROMYALGIA` (3,623 / 357,549). The exposures are unchanged \u2014 TMFI has been genotyped in no other cohort \u2014 so the replication is outcome-side only. FinnGen fits a logistic mixed model (REGENIE), so its coefficients are log odds ratios; this was verified rather than assumed, the identity SE\u00b2 \u00d7 2Nf(1\u2212f) = 1/[\u03bc(1\u2212\u03bc)] holding to within 5% on non-significant instruments while the linear alternative is excluded by factors of 17 to 10,300, and no rescaling was applied. Variants were anchored on rsid (b37 coordinates, no liftover); 108 of the 109 union instruments were present, none ambiguous, median |allele-frequency difference| 0.035. Stage 1 ran on all 42 instruments and Stage 2 on 92 of 93, with the same estimators and the same double implementation as the preceding stages.",

"### Reproducibility and statistical reporting",

"All estimates, tables and figures were regenerated from the frozen stage CSV outputs by `manuscript_numbers.py`, the single source of truth for every number quoted. Software: PLINK v1.9.0-rc2; R 4.6.1 with MVMR 0.4.8 (Stage 2); Python 3.13 with NumPy, pandas and SciPy. Every estimator applied to the site-matched and FinnGen outcomes was computed twice by independently written implementations sharing only the frozen harmonized inputs, agreeing to at least six significant figures (Supplementary Table S10). Confidence intervals are estimate \u00b1 1.96 \u00d7 SE with two-sided normal P values, mutually consistent by construction. Missing data were handled by complete-case harmonization, with no imputation. No multiple-testing correction was applied, because the primary contrast is a single pre-specified comparison and the site-matched outcomes were pre-specified replications of it. The study was not preregistered; instrument sets, outcome definitions, estimator hierarchy and sensitivity analyses were fixed before the Stage 2 estimates were inspected (frozen files deposited with the code archive). Reporting follows the STROBE-MR statement [16]; the completed checklist is Supplementary Table S12.",

"## Results",

"### Stage 1: the total effect of TMFI on chronic pain intensity was null",

"Using 42 LD-clumped instruments, the IVW-MRE total effect was null for both outcomes (Supplementary Table S14): 0.008 SD per SD for CPI_full (95% CI \u22120.018 to 0.033, P=0.558) and 0.005 SD per SD for CPI_no_imaging (\u22120.022 to 0.032, P=0.721). Weighted-median and MR-Egger estimates agreed (all P\u22650.49), MR-Egger intercepts were small and non-significant (0.0022, P=0.344; 0.0025, P=0.324), heterogeneity was modest (Cochran Q=53.71 on 41 degrees of freedom, P=0.088; Q=48.75, P=0.189), Steiger directionality supported the TMFI\u2192CPI orientation for 42 of 42 instruments, and leave-one-out estimates stayed within \u22120.004 to 0.018. This stage was powered to detect 0.037 and 0.039 SD per SD at 80% power; the observed estimates were roughly five-fold smaller and their intervals exclude any total effect above 0.033 SD. Stage 2 was therefore needed to determine whether that null conceals an effect visible only once BMI is held constant, and to bound any such effect.",

"### Stage 2 primary direct effects: BMI positive, TMFI null",

"Joint selection produced 93 independent variants (59 BMI-only, 32 TMFI-only, 2 shared) with conditional F statistics of 22.12 and 35.23, identical across outcomes because the same 93-SNP definition was carried forward. In the primary CPI_full analysis the BMI direct effect was 0.090 SD of pain intensity per SD of BMI (95% CI 0.038\u20130.143; P=0.00075) against a null TMFI direct effect of \u22120.002 SD per SD (\u22120.036 to 0.031; P=0.898) (Table 2). CPI_no_imaging, a subset of CPI_full rather than an independent replication, gave the same pattern: BMI 0.101 SD per SD (0.044\u20130.159; P=0.00058), TMFI \u22120.008 SD per SD (\u22120.045 to 0.029; P=0.656). The two TMFI estimates carry signs opposite to the BMI estimates, so the null is not a down-weighted copy of the BMI signal.",

"Fixed-effect IVW returned nearly identical point estimates with smaller standard errors (BMI P=3.5\u00d710\u207b\u2076 and 3.6\u00d710\u207b\u2076), the roughly 27% reduction following from dispersion statistics \u03c6=1.90 and 1.81 rather than from any change in effect size; MVMR-Egger intercepts were indistinguishable from zero (\u22120.0013, P=0.481; 0.0004, P=0.863), with slopes preserving the contrast (BMI 0.272, P=0.00060 and 0.294, P=0.00091). The significant modified Q statistic therefore reflects balanced scatter rather than directional pleiotropy.",

"### Power bounds make the null informative",

"Stage 2 had 80% power to detect a TMFI direct effect of 0.048 SD per SD for CPI_full and 0.053 for CPI_no_imaging (Table 4); the loss of precision relative to Stage 1 (0.037, 0.039) is the variance inflation caused by conditioning on BMI, and is why a multivariable null cannot be interpreted without a stated bound. The bound is interpretable against the comparator estimated in the same model: the BMI effect was detected at approximately 92% and 93% power, so the minimum TMFI effect detectable at 80% power was about 1.9-fold smaller than the BMI effect the same model recovered. The upper confidence limits exclude any positive TMFI effect larger than 0.031 SD for CPI_full and 0.029 SD for CPI_no_imaging \u2014 35% and 28% of the respective BMI point estimates. Power reached 83% and 75% for a 0.05 SD TMFI effect and exceeded 99% for 0.075 SD. These data are inconsistent with a TMFI effect of the order of the adiposity effect; smaller effects remain compatible.",

"### The pattern was stable under influence-point, reweighting and targeted-exclusion analyses",

"Leave-one-out estimates preserved the pattern across all 93 variants, and no single variant determined the exposure contrast (Supplementary Figure S1). Across Huber-IRLS reweighting, Cook-distance exclusions and targeted exclusions, BMI remained positive with confidence intervals excluding zero, whereas TMFI estimates ranged from \u22120.032 to \u22120.001 in the source units (Supplementary Table S15; Supplementary Figure S2). For CPI_full, the Cook-threshold and five-target analyses used the same excluded set and therefore represent the same fitted model.",

"### Sample overlap was real but left both findings unchanged",

"The nested-subset check (CPI_no_imaging is by design a subset of CPI_full) returned a cross-trait intercept of 0.910 (SE 0.010) at its structural ceiling of 0.890 within ~2 SE, with genetic correlation 1.0025 (Supplementary Table S16). Against that calibration, TMFI and CPI_full showed genuine overlap (intercept 0.052, SE 0.008, P=1.2\u00d710\u207b\u00b9\u00b9; N_shared\u00b7\u03c1_phen = 4,112), reduced more than four-fold after excluding brain-imaging participants (0.012, P=0.097). No overlap was detectable between the GIANT BMI GWAS and either outcome (\u22120.015 and \u22120.009, both P>0.08) or between the two exposure GWAS (\u22120.018, SE 0.009). Propagating the measured intercepts changed nothing material: Stage 1 bias was +0.0005 and +0.0001 SD per SD (4% and 1% of the reported standard errors), and in Stage 2 the TMFI direct effect moved by at most 0.0014 SD per SD while the BMI direct effect moved by at most 0.004, remaining significant at P\u22649.6\u00d710\u207b\u2076 in every scenario, including one with the BMI\u2013outcome overlap assumed three times larger than measured. The intercept required to move the TMFI effect to significance, or the BMI effect to non-significance, was 1.27, 1.24, 1.03 and 1.28 \u2014 above the structural ceilings of 0.587\u20130.765 in every case \u2014 so under the stated assumptions no configuration of overlapping participants could reverse either Stage 2 finding.",

"### Site-matched outcomes: the TMFI null extends to knee and back pain",

"Neither site-matched outcome produced a TMFI effect (Table 5). Stage 1 gave odds ratios per SD of 1.0307 (P=0.20) for knee pain and 1.0135 (P=0.53) for back pain; Stage 2, conditioning on BMI, gave 1.0146 (P=0.58) and 1.0005 (P=0.98). The second knee-pain GWAS returned a concordant univariable null (0.9965, P=0.87). The internal comparator again behaved as expected and differentiated the sites: BMI gave 1.2751 (P=2.5\u00d710\u207b\u2079) for knee pain and 1.1344 (P=9.4\u00d710\u207b\u2074) for back pain \u2014 descriptively larger at the weight-bearing joint, consistent with mechanical loading, though the intervals overlap and the site contrast is not formally tested. Design sensitivity was tighter: 80% power was retained for TMFI effects of 1.0681 (knee) and 1.0621 (back) per SD, approximately 27% and 48% of the BMI effect the same model recovered; on the risk-difference scale, TMFI gave +0.24 and +0.01 percentage points against +4.08 and +2.40 for BMI. Overlap with these outcomes was smaller (intercepts 0.014\u20130.017 against ceilings of 0.316\u20130.328; induced bias at most 0.024 standard errors). Genetically predicted TMFI and knee pain are correlated (r_g=0.343), but conditioning on BMI (r_g=0.433 with knee pain) left no detectable direct association.",

"### External replication: both contrasts hold with no participant overlap",

"Repeating both stages against three FinnGen endpoints reproduced both central findings (Table 6). No TMFI effect emerged: multivariable estimates of 1.019 per SD (P=0.42) for `PAIN`, 0.983 (P=0.53) for limb pain and 1.001 (P=0.99) for fibromyalgia, with the univariable estimates agreeing. The comparator behaved as expected: BMI was associated with limb pain at 1.217 per SD (1.116\u20131.327, P=9.0\u00d710\u207b\u2076) and pointed in the same direction, non-significantly, for `PAIN` (1.055, P=0.15) and fibromyalgia (1.037, P=0.79). The limb-pain estimate lies between the UK Biobank estimates for back pain (1.134) and knee pain (1.275), so the independent cohort reproduces both the TMFI null and the adiposity gradient across weight-bearing sites. Conditional F statistics were 22.4 and 35.4 (against 22.1 and 35.2 in UK Biobank), Steiger directionality held for 42 of 42 instruments, and leave-one-out left the limb-pain BMI coefficient between 0.176 and 0.207 (full-sample 0.196). Heterogeneity was present (modified Q 3.03, 1.46, 1.20), so the random-effects specification is the relevant one; the Egger intercept was nominally significant for `PAIN` alone (P=0.041), within expectation for three endpoints under the declared no-correction rule. Minimum detectable effects were 1.07, 1.09 and 1.26 per SD, so the fibromyalgia endpoint contributes descriptive rather than inferential support. Because no participants are shared, the cross-trait intercept is zero by construction, so the concordance of the two designs indicates the bounded null is not an artefact of overlap.",

"## Discussion",

"### Principal finding",

"Across two stages, two outcome definitions and every robustness specification, genetically predicted BMI showed a positive association with chronic pain intensity while genetically predicted thigh muscle fat infiltration showed no evidence of an effect, either total or conditional on BMI. The contribution is the conversion of that null into a quantitative statement: the same model that returned a null TMFI estimate recovered a BMI effect of 0.090 SD per SD at approximately 92% power while retaining 80% power for a TMFI effect as small as 0.048 SD per SD, so the upper confidence limits exclude any positive TMFI effect larger than 35% (CPI_full) and 28% (CPI_no_imaging) of the BMI point estimate. Both contrasts were recovered in FinnGen, which shares no participants with UK Biobank (Table 6), so the bounded null is not an artefact of the exposure\u2013outcome overlap that the tipping-point arguments address.",

"Two features make this more than a failed test. The internal comparator shows that the instrument set, harmonization and outcome data could detect a body-composition effect in this specific contrast; a null obtained alongside such a comparator differs from a null in isolation. The comparator is not a positive control with a known answer \u2014 it establishes the model's sensitivity in this contrast, not that every excluded effect would have been recovered \u2014 and the opposite-signed estimates mean the TMFI result is not a weakened echo of the adiposity signal. The conclusion is bounded: thigh muscle fat infiltration is unlikely to drive chronic pain intensity at anything like the scale of systemic adiposity, but effects below about 0.05 SD per SD remain open, and muscle tissue may still matter through other phenotypes or pathways.",

"### Relation to existing evidence",

"Muscle\u2013pain Mendelian randomization has linked lower-extremity muscle traits to spinal diseases [7] and sarcopenia to chronic pain [6]; those analyses do not condition on adiposity, so they are not in conflict with ours \u2014 a muscle\u2013pain link may be carried by phenotypes other than intramuscular fat or by pathways shared with adiposity, and distinguishing these requires the joint conditioning reported here. Likewise, the adiposity-to-pain analyses established higher BMI as a cause of back, knee and hip pain [19,20], with which our BMI estimates agree in direction and magnitude, but they could not ask whether a muscle compartment contributes independently, because no comparator exposure was estimated in the same model. A structured comparison with the nearest published analyses is given in Supplementary Table S13.",

"### Why the joint model changes the interpretation",

"TMFI summarizes local muscle quality, adipose replacement and metabolic state; BMI captures systemic adiposity without resolving where fat is deposited or how muscle is remodelled. Conditioning on BMI asks what remains attributable to TMFI-specific variation once systemic adiposity is accounted for. The multivariable stage was less precise than the univariable stage \u2014 the 80%-power bound rose from 0.037 to 0.048 SD per SD \u2014 so a multivariable null reported without its bound cannot be distinguished from low power; conversely, the near-zero estimate sits against a comparator that was unambiguously detected, and that pairing is what licenses the interpretation.",

"### Robustness and interpretation boundary",

"The exposure contrast was stable under leave-one-out analysis, Huber-IRLS reweighting, Cook-distance exclusions and targeted removal of the five review-priority variants; precision changed modestly, the pattern did not. The diagnostics bound what can be claimed rather than merely qualifying it: the modified Q statistics were strongly positive, so absolute magnitudes should not be over-interpreted; the MVMR-Egger intercepts were indistinguishable from zero, providing no evidence of average directional pleiotropy although pleiotropic bias cannot be excluded; the fixed-effect fit reproduced both point estimates, so heterogeneity inflates variance rather than shifting estimates; and conditional F statistics of 22.12 and 35.23 reduce concern about severe conditional weak-instrument bias without excluding it.",

"### Implications for pain and muscle research",

"Methodologically, null results in multivariable MR should be reported with the minimum detectable effect and, where one exists, an internal comparator; absent both, a null cannot be distinguished from low power. Comparator-based power calibration of this kind is only beginning to be applied explicitly [22]. Substantively, systemic adiposity was associated with chronic pain intensity whereas thigh muscle composition was not; because the BMI instrument indexes adiposity as a whole, this localises the difference no further than systemic adiposity versus thigh muscle composition. These estimates describe lifelong genetically predicted differences and do not support advancing thigh muscle fat infiltration as an intervention target for chronic pain intensity while adiposity itself is unaddressed. The adjacent positive findings for other muscle phenotypes [6,7] suggest the next useful test: jointly conditioning several correlated muscle traits to establish whether any muscle compartment carries an adiposity-independent effect.",

"### Limitations",

"Eight boundaries qualify the estimates: (i) the TMFI source GWAS carries an unresolved preprocessing discrepancy affecting the absolute exposure scale; all inference is therefore standardized, and raw coefficients should not be translated into physical units until the source authors confirm the flag; (ii) the LD reference panel (503 individuals, 1.84 million variants) contained only 22.75% of TMFI candidate variants, so eligible instruments outside it may have been missed, although every retained variant was present and inference does not depend on re-clumping; (iii) residual heterogeneity indicates that some variants influence CPI through pathways outside both exposures; (iv) the overlap audit identifies the product N_shared\u00b7\u03c1_phen rather than a participant head-count, so the overlap figures should not be read as participant counts; (v) the GIANT release is genomic-control corrected (\u03bb=1.61 inferred from its single-trait intercept), so intercepts involving BMI required rescaling, its released standard errors are inflated by about 27% \u2014 making the BMI conditional F conservative by a factor of about 1.6 and sd<sub>BMI</sub> set to 1.00 by definition; (vi) the second knee-pain source returned a BMI association opposite to the epidemiological expectation (OR per SD 0.905, P=0.008), most plausibly because that GWAS adjusted for adiposity-related covariates, so it is used only to check the TMFI null and its BMI coefficient is not interpreted; (vii) the site-matched odds ratios rest on the approximate BOLT-LMM conversion log OR = \u03b2/[\u03bc(1\u2212\u03bc)], although every P value, F statistic and tipping point is invariant to it and the risk differences reported alongside require no approximation; and (viii) the external replication is outcome-side only \u2014 the exposures remain UK Biobank instruments, and the FinnGen endpoints are binary ICD-based diagnoses, so concordance is related-phenotype triangulation rather than same-phenotype replication, with the fibromyalgia endpoint (3,623 cases) able to exclude only effects larger than 1.26 per SD and therefore contributing descriptive rather than inferential support. Generalizability is correspondingly limited: all source GWAS were conducted in European-ancestry participants; the estimates describe lifetime-averaged genetically predicted differences rather than interventions at a particular age or duration; and they apply to a shift across the full genetic range of the exposure, not to clinically attainable changes.",

"## Conclusion",

"Genetically predicted BMI was positively associated with chronic pain intensity after conditioning on thigh muscle fat infiltration, whereas TMFI showed no evidence of a total or conditional direct effect. Because the same model recovered the BMI effect at approximately 92% power while retaining 80% power for a TMFI effect of 0.048 SD per SD, this is a bounded null: thigh muscle fat infiltration is unlikely to influence chronic pain intensity at anywhere near the scale of systemic adiposity, although effects below about 0.05 SD per SD remain open. The same pattern held for site-matched knee and back pain, with a descriptively larger BMI effect at the weight-bearing joint, and for three FinnGen endpoints with no participant overlap. The next decisive test is a replication that preserves this joint structure, uses a denser LD reference panel for instrument selection, and extends the conditioning to further correlated muscle traits.",
])

# ---------------- 3b. Final trim pass ----------------
TRIMS = [
    ("quantified exposure\u2013outcome sample overlap genetically, analysed two site-matched outcomes (knee and back pain), and replicated both stages in FinnGen, a cohort sharing no participants with UK Biobank.",
     "quantified sample overlap genetically, analysed site-matched knee and back pain, and replicated both stages in FinnGen, which shares no participants with UK Biobank."),
    ("has recently been characterized [1], although with infiltration treated as the outcome rather than as an exposure on pain [21]",
     "has recently been characterized [1], although with infiltration as the outcome rather than the exposure [21]"),
    ("A null estimate is informative only if", "A null is informative only if"),
    ("Conditional instrument strength was assessed with multivariable conditional F statistics from the R `MVMR` package [5], the values reported throughout (Tables 3 and 5); a closed-form approximation computed independently returns 23.6 for TMFI and 42.8 for BMI against 22.12 and 35.23 from `MVMR` (Supplementary Table S11), both far above 10.",
     "Conditional strength was assessed with multivariable conditional F statistics from the R `MVMR` package [5], reported throughout (Tables 3 and 5); a closed-form approximation computed independently returns 23.6 (TMFI) and 42.8 (BMI) against 22.12 and 35.23 (Supplementary Table S11), both far above 10."),
    ("Two checks validated the pipeline: all three UK Biobank single-trait intercepts fell within 2% of 1 (0.983, 1.000, 1.010), and the nested-subset check \u2014 CPI_no_imaging is a subset of CPI_full and therefore maximally overlapping by design \u2014 returned an intercept of 0.910 (SE 0.010) at its structural ceiling of 0.890 within ~2 SE, with genetic correlation 1.0025.",
     "Two checks validated the pipeline: the three UK Biobank single-trait intercepts fell within 2% of 1 (0.983, 1.000, 1.010), and the nested-subset check (CPI_no_imaging is a subset of CPI_full, maximally overlapping by design) recovered its structural ceiling of 0.890 within ~2 SE (intercept 0.910, SE 0.010; genetic correlation 1.0025)."),
    ("Phenotype identities were verified independently of repository metadata \u2014 the knee source leads at rs143384 in GDF5 (\u2212log\u2081\u2080P=21.0), absent from back pain, which leads near SOX5, and allele frequencies agree with the reference panel at r\u22480.998 (full checks in Supplementary Methods).",
     "Phenotype identities were verified independently of repository metadata (leading signals GDF5 for knee and SOX5-adjacent for back pain; allele-frequency agreement with the panel r\u22480.998; full checks in Supplementary Methods)."),
    ("so its coefficients are log odds ratios; this was verified rather than assumed, the identity SE\u00b2 \u00d7 2Nf(1\u2212f) = 1/[\u03bc(1\u2212\u03bc)] holding to within 5% on non-significant instruments while the linear alternative is excluded by factors of 17 to 10,300, and no rescaling was applied.",
     "so its coefficients are log odds ratios \u2014 verified rather than assumed: the identity SE\u00b2 \u00d7 2Nf(1\u2212f) = 1/[\u03bc(1\u2212\u03bc)] holds to within 5% on non-significant instruments, the linear alternative is excluded by factors of 17 to 10,300, and no rescaling was applied."),
    ("Missing data were handled by complete-case harmonization, with no imputation.",
     "Missing data were handled by complete-case harmonization."),
    ("Weighted-median and MR-Egger estimates agreed (all P\u22650.49), MR-Egger intercepts were small and non-significant (0.0022, P=0.344; 0.0025, P=0.324), heterogeneity was modest (Cochran Q=53.71 on 41 degrees of freedom, P=0.088; Q=48.75, P=0.189), Steiger directionality supported the TMFI\u2192CPI orientation for 42 of 42 instruments, and leave-one-out estimates stayed within \u22120.004 to 0.018.",
     "Weighted-median and MR-Egger estimates agreed (all P\u22650.49) with small, non-significant Egger intercepts (0.0022, P=0.344; 0.0025, P=0.324); heterogeneity was modest (Q=53.71 on 41 df, P=0.088; Q=48.75, P=0.189); Steiger directionality supported the TMFI\u2192CPI orientation for 42 of 42 instruments; and leave-one-out estimates stayed within \u22120.004 to 0.018."),
    ("identical across outcomes because the same 93-SNP definition was carried forward",
     "identical across outcomes because the same instrument set was carried forward"),
    (", the roughly 27% reduction following from dispersion statistics \u03c6=1.90 and 1.81 rather than from any change in effect size;",
     ", the roughly 27% reduction following from dispersion statistics \u03c6=1.90 and 1.81;"),
    ("Conditional F statistics were 22.4 and 35.4 (against 22.1 and 35.2 in UK Biobank), Steiger directionality held for 42 of 42 instruments, and leave-one-out left the limb-pain BMI coefficient between 0.176 and 0.207 (full-sample 0.196). Heterogeneity was present (modified Q 3.03, 1.46, 1.20), so the random-effects specification is the relevant one; the Egger intercept was nominally significant for `PAIN` alone (P=0.041), within expectation for three endpoints under the declared no-correction rule.",
     "Conditional F statistics were 22.4 and 35.4 (UK Biobank: 22.1 and 35.2), Steiger held for 42 of 42 instruments, and leave-one-out left the limb-pain BMI coefficient between 0.176 and 0.207 (full-sample 0.196). Heterogeneity was present (modified Q 3.03, 1.46, 1.20); the Egger intercept was nominally significant for `PAIN` alone (P=0.041), within expectation for three endpoints under the no-correction rule."),
    ("The bound is interpretable against the comparator estimated in the same model: the BMI effect was detected at approximately 92% and 93% power, so the minimum TMFI effect detectable at 80% power was about 1.9-fold smaller than the BMI effect the same model recovered.",
     "The bound is interpretable against the comparator in the same model: the BMI effect was detected at approximately 92% and 93% power, so the minimum TMFI effect detectable at 80% power was about 1.9-fold smaller than the BMI effect recovered."),
    ("BMI captures systemic adiposity without resolving where fat is deposited or how muscle is remodelled",
     "BMI captures systemic adiposity without resolving fat deposition or muscle remodelling"),
    ("The diagnostics bound what can be claimed rather than merely qualifying it:",
     "The diagnostics bound what can be claimed:"),
    ("suggest the next useful test: jointly conditioning several correlated muscle traits to establish whether any muscle compartment carries an adiposity-independent effect.",
     "suggest the next test: jointly conditioning several correlated muscle traits for an adiposity-independent effect."),
    ("but they could not ask whether a muscle compartment contributes independently, because no comparator exposure was estimated in the same model.",
     "but they could not ask whether a muscle compartment contributes independently, lacking a comparator in the same model."),
    ("most plausibly because that GWAS adjusted for adiposity-related covariates, so it is used only to check the TMFI null and its BMI coefficient is not interpreted",
     "plausibly because that GWAS adjusted for adiposity-related covariates; it is used only to check the TMFI null"),
    ("with the fibromyalgia endpoint (3,623 cases) able to exclude only effects larger than 1.26 per SD and therefore contributing descriptive rather than inferential support",
     "the fibromyalgia endpoint (3,623 cases) excluding only effects larger than 1.26 per SD"),
    ("Generalizability is correspondingly limited: all source GWAS were conducted in European-ancestry participants; the estimates describe lifetime-averaged genetically predicted differences rather than interventions at a particular age or duration; and they apply to a shift across the full genetic range of the exposure, not to clinically attainable changes.",
     "Generalizability is limited to European ancestry, to lifetime-averaged genetically predicted differences rather than interventions, and to shifts across the full genetic range of the exposure rather than clinically attainable changes."),
    ("uses a denser LD reference panel for instrument selection, and extends",
     "uses a denser LD panel, and extends"),
    ("Existing muscle\u2013pain genetic analyses share two weaknesses \u2014 they rarely include a comparator exposure able to demonstrate that the model detects an effect when one exists, and they rarely state the smallest effect the design could have detected \u2014 so a null reported without either is descriptively correct but interpretively uninformative.",
     "Existing muscle\u2013pain genetic analyses share two weaknesses \u2014 they rarely include a comparator exposure demonstrating that the model detects an effect when one exists, and they rarely state the smallest effect the design could detect \u2014 so a null without either is descriptively correct but interpretively uninformative."),
    ("The intercept required to move the TMFI effect to significance, or the BMI effect to non-significance, was 1.27, 1.24, 1.03 and 1.28 \u2014 above the structural ceilings of 0.587\u20130.765 in every case \u2014 so under the stated assumptions no configuration of overlapping participants could reverse either Stage 2 finding.",
     "The intercept required to move either Stage 2 estimate across the significance threshold was 1.27, 1.24, 1.03 and 1.28 \u2014 above the structural ceilings of 0.587\u20130.765 in every case \u2014 so under the stated assumptions no configuration of overlapping participants could reverse either finding."),
    ("approximately 27% and 48% of the BMI effect the same model recovered; on the risk-difference scale",
     "approximately 27% and 48% of the BMI effect recovered; on the risk-difference scale"),
    ("These estimates describe lifelong genetically predicted differences and do not support advancing thigh muscle fat infiltration as an intervention target for chronic pain intensity while adiposity itself is unaddressed.",
     "These estimates describe lifelong genetically predicted differences and do not support thigh muscle fat infiltration as an intervention target for chronic pain while adiposity itself is unaddressed."),
    ("sources, ancestries and effective sample sizes are given in Table 1",
     "sources and effective sample sizes are given in Table 1"),
    ("We addressed both weaknesses in two analysis stages applied to the same summary data",
     "We addressed both in two analysis stages applied to the same summary data"),
    ("and permitting the null TMFI estimate to be read against the study's minimum detectable effect",
     "and permitting the null TMFI estimate to be read against the minimum detectable effect"),
    ("the single source of truth for every number quoted",
     "the single source of truth for every number"),
    ("though the intervals overlap and the site contrast is not formally tested",
     "though the intervals overlap and the contrast is not formally tested"),
    ("smaller effects remain compatible", "smaller ones remain compatible"),
    ("approximately 27% and 48% of the BMI effect recovered",
     "27% and 48% of the BMI effect recovered"),
    ("and to bound any such effect.", "and to bound it."),
    ("Stage 1 estimated the univariable (total) effect of TMFI on CPI. Stage 2 estimated, for each exposure, its direct effect on CPI conditional on the other exposure using two-exposure summary-data MVMR.",
     "Stage 1 estimated the univariable (total) effect of TMFI on CPI; Stage 2 estimated each exposure's direct effect on CPI conditional on the other using two-exposure summary-data MVMR."),
    ("showed no evidence of an effect, either total or conditional on BMI. The contribution is",
     "showed no evidence of an effect, total or conditional. The contribution is"),
    ("The multivariable stage was less precise than the univariable stage \u2014 the 80%-power bound rose from 0.037 to 0.048 SD per SD \u2014 so a multivariable null reported without its bound cannot be distinguished from low power; conversely, the near-zero estimate sits against a comparator that was unambiguously detected, and that pairing is what licenses the interpretation.",
     "The multivariable stage was less precise than the univariable \u2014 the 80%-power bound rose from 0.037 to 0.048 SD per SD \u2014 so a multivariable null without its bound cannot be distinguished from low power; conversely, the near-zero estimate sits against a comparator that was unambiguously detected, and that pairing licenses the interpretation."),
    ("The internal comparator again behaved as expected and differentiated the sites:",
     "The comparator again behaved as expected and differentiated the sites:"),
    ("Body composition is central to this literature because systemic adiposity may influence pain through mechanical loading, inflammation, metabolic dysfunction and altered pain processing, and genetically predicted adiposity is already an established cause of site-specific pain [19,20].",
     "Systemic adiposity may influence pain through mechanical loading, inflammation, metabolic dysfunction and altered pain processing, and genetically predicted adiposity is an established cause of site-specific pain [19,20]."),
    ("is the variance inflation caused by conditioning on BMI, and is why",
     "is the variance inflation caused by conditioning, and is why"),
]
for old, new in TRIMS:
    assert BODY.count(old) == 1, "trim fail: %r" % old[:60]
    BODY = BODY.replace(old, new)

i_intro = text.index("## Introduction")
i_main = text.index("## Main tables")
text = text[:i_intro] + BODY + NL + NL + text[i_main:]

# ---------------- 4. Tables / figures / supplementary ----------------
i_mt = text.index("## Main tables")
i_fl = text.index("## Figure legends")
i_sm = text.index("## Supplementary material")
i_da = text.index("## Data availability")

mt_block = text[i_mt:i_fl]
fl_block = text[i_fl:i_sm]
sm_block = text[i_sm:i_da]

def split_tables(block):
    lines = block.split(NL)
    tables = {}
    cur_name, cur = None, []
    for ln in lines:
        m = re.match(r"^### Table (\d+)\.", ln)
        if m:
            if cur_name is not None:
                tables[cur_name] = cur
            cur_name, cur = "T" + m.group(1), [ln]
        elif cur_name is not None:
            cur.append(ln)
    if cur_name is not None:
        tables[cur_name] = cur
    return tables

tabs = split_tables(mt_block)
assert set(tabs) == {"T1","T2","T3","T4","T5","T6","T7","T8","T9"}, sorted(tabs)

def strip_block(lines):
    while lines and lines[-1].strip() == "":
        lines = lines[:-1]
    return NL.join(lines)

def patch(block, pairs, name):
    for old, new in pairs:
        assert block.count(old) == 1, "patch fail in %s: %r found %d" % (name, old[:50], block.count(old))
        block = block.replace(old, new)
    return block

t1 = strip_block(tabs["T1"])

t2 = patch(strip_block(tabs["T3"]), [("### Table 3.", "### Table 2."),
    ("listed in Table 6", "listed in Table 4")], "T3->T2")

t3 = patch(strip_block(tabs["T5"]), [("### Table 5.", "### Table 3."),
    ("Tables 5 and 8 therefore quote", "Tables 3 and 5 therefore quote"),
    ("see Figure 3 | see Figure 3", "see Supplementary Figure S1 | see Supplementary Figure S1")], "T5->T3")

t4 = patch(strip_block(tabs["T6"]), [("### Table 6.", "### Table 4.")], "T6->T4")

t5 = patch(strip_block(tabs["T8"]), [("### Table 8.", "### Table 5."),
    ("reported in Table 5 (R `MVMR`", "reported in Table 3 (R `MVMR`")], "T8->T5")

t6 = patch(strip_block(tabs["T9"]), [("### Table 9.", "### Table 6.")], "T9->T6")

s14 = patch(strip_block(tabs["T2"]), [("### Table 2.", "### Supplementary Table S14.")], "T2->S14")

s15 = patch(strip_block(tabs["T4"]), [("### Table 4.", "### Supplementary Table S15."),
    ("matching the convention in Table 3 and", "matching the convention in Table 2 and"),
    ("(see Transparency note in Results)", "(see Results)")], "T4->S15")

s16 = patch(strip_block(tabs["T7"]), [
    ("### Table 7.", "### Supplementary Table S16."),
    ("Intercept needed (Table 3 SE)", "Intercept needed (Table 2 SE)"),
    ("quoted against the standard errors reported in Table 2 and for Stage 2",
     "quoted against the standard errors reported in Supplementary Table S14 and for Stage 2"),
    ("multiplicative-random-effects errors of Table 3 gives",
     "multiplicative-random-effects errors of Table 2 gives")], "T7->S16")

fl_lines = fl_block.split(NL)
figs = {}
cur_name, cur = None, []
for ln in fl_lines:
    m = re.match(r"^### Figure (\d+)\.", ln)
    if m:
        if cur_name is not None:
            figs[cur_name] = cur
        cur_name, cur = "F" + m.group(1), [ln]
    elif cur_name is not None:
        cur.append(ln)
if cur_name is not None:
    figs[cur_name] = cur
assert set(figs) == {"F1","F2","F3","F4"}, sorted(figs)

f1 = patch(strip_block(figs["F1"]), [("are given in Table 3", "are given in Table 2"),
                                     ("(Results; Table 9)", "(Results; Table 6)")], "F1")
f2 = patch(strip_block(figs["F2"]), [("column of Table 3", "column of Table 2")], "F2")
sf1 = patch(strip_block(figs["F3"]), [("### Figure 3.", "### Supplementary Figure S1.")], "F3->SF1")
sf2 = patch(strip_block(figs["F4"]), [("### Figure 4.", "### Supplementary Figure S2.")], "F4->SF2")

with open(PATH + ".bak_jcsm2026-09-20", "rb") as f:
    old_text = f.read().decode("utf-8")

def grab_subsection(src, header):
    i = src.index(header)
    j = src.index(NL + "### ", i + len(header))
    return src[i:j].rstrip(NL)

sm_bmi = grab_subsection(old_text, "### Scale of the BMI summary statistics")
sm_scale = grab_subsection(old_text, "### Effect-size scale of the site-matched outcomes")
sm_reimpl = grab_subsection(old_text, "### Independent re-implementation of the site-matched estimators")
sm_integrity = grab_subsection(old_text, "### Source phenotype integrity and one unresolved data conflict")
sm_pheno = grab_subsection(old_text, "### Pre-specified site-matched secondary analysis")
ph_body = sm_pheno.split(NL)[2:].join([""]).split(NL + NL) if False else NL.join(sm_pheno.split(NL)[2:]).split(NL + NL)
sm_phenoid = "### Phenotype identity checks and harmonization of the site-matched instruments" + NL + NL + NL.join(ph_body[1:3])
assert "rs143384" in sm_phenoid and "0.998" in sm_phenoid and "palindromic" in sm_phenoid

ADD_LIST = NL.join([
"15. **Supplementary Methods** \u2014 Full derivations moved from the main text to meet the word limit, unchanged in content: the scale of the BMI summary statistics; the effect-size scale of the site-matched outcomes; the independent re-implementation of the site-matched estimators; source phenotype integrity and the unresolved data conflict; and phenotype identity checks and harmonization for the site-matched instruments.",
"16. **Supplementary Table S14** \u2014 Stage 1 univariable MR estimates of the total effect of TMFI on both CPI definitions, all estimators (moved from the main text; per-instrument data in Supplementary Table S2).",
"17. **Supplementary Table S15** \u2014 Stage 2 targeted-exclusion sensitivity analyses (moved from the main text).",
"18. **Supplementary Table S16** \u2014 Sample-overlap audit for the primary outcomes: cross-trait LD score intercepts, induced bias and tipping-point analysis (moved from the main text).",
"19. **Supplementary Figure S1** \u2014 Leave-one-out MVMR estimates for both outcomes and both exposures (moved from the main text).",
"20. **Supplementary Figure S2** \u2014 Influence-point diagnostics and the robustness envelope (moved from the main text).",
])
assert sm_block.rstrip(NL).endswith("no participant overlap."), "unexpected supp list tail"
sm_block_new = sm_block.rstrip(NL) + NL + ADD_LIST + NL + NL

new_mid = NL.join([x for x in [
    "## Main tables", "",
    t1, "", t2, "", t3, "", t4, "", t5, "", t6, "",
    "## Figure legends", "",
    f1, "", f2, "",
    sm_block_new.rstrip(NL), "",
    "## Supplementary Methods", "",
    "The following subsections were moved from the main text without change of content to meet the journal word limit; each is referenced from the main text where its result is summarized.", "",
    sm_bmi, "", sm_scale, "", sm_reimpl, "", sm_integrity, "", sm_phenoid, "",
    "## Supplementary tables (moved from the main text)", "",
    s14, "", s15, "", s16, "",
    "## Supplementary figure legends (moved from the main text)", "",
    sf1, "", sf2, "",
]]) + NL + NL

text = text[:i_mt] + new_mid + text[i_da:]

# ---------------- 5. Declarations: LLM statement ----------------
OLD_AI = '**Use of artificial intelligence tools.** Declared in the Methods, under "Use of large language models".'
assert text.count(OLD_AI) == 1
NEW_AI = ("**Use of artificial intelligence tools.** A large language model (WorkBuddy, an AI writing assistant) "
          "was used for three bounded tasks: checking grammar and wording, reorganising prose for concision, and "
          "cross-checking that every number quoted in the text, tables and figure legends matches the value in the "
          "underlying result files. It was not used to generate, clean or analyse any data, to select instruments or "
          "outcomes, to choose or implement statistical methods, or to draw scientific conclusions; no result, table "
          "or figure value originates from the model. Every numerical statement in this manuscript is regenerated "
          "deterministically by `manuscript_numbers.py` from the frozen Stage 1, Stage 2, Stage 2B, Stage 3 and "
          "Stage 5 CSV outputs, and the two independent implementations of the estimators (`stage5_mr.py` under "
          "Python and `stage5_validation_baser.R` under base R) were written and reviewed by the authors. The "
          "authors have verified all content and take full responsibility for it.")
text = text.replace(OLD_AI, NEW_AI)

# ---------------- 6. Assertions ----------------
checks = [
    "0.048", "0.090", "0.038\u20130.143", "0.00075", "\u22120.002", "0.898",
    "0.008 SD per SD", "P=0.558", "0.005 SD per SD", "P=0.721",
    "22.12", "35.23", "0.052", "1.2\u00d710\u207b\u00b9\u00b9", "0.0014",
    "1.0146", "P=0.58", "1.0005", "P=0.98", "1.2751", "1.1344",
    "0.983", "P=0.53", "1.217", "1.116\u20131.327", "9.0\u00d710\u207b\u2076",
    "46,347", "134,627", "106,697", "499,362", "80% power", "0.05 SD",
    "descriptively larger", "internal comparator",
    "at its structural ceiling of 0.890 within ~2 SE", "Supplementary Table S13",
    "Eight boundaries", "Supplementary Table S14", "Supplementary Table S15",
    "Supplementary Table S16", "Supplementary Figure S1", "Supplementary Figure S2",
    "Supplementary Methods", "bounded null",
]
for c in checks:
    assert c in text, "MISSING after edit: " + c

sup_start = text.index("## Supplementary Methods")
main_part, sup_part = text[:sup_start], text[sup_start:]
for bad in ["Table 7", "Table 8", "Table 9"]:
    assert bad not in main_part, "stale ref in main part: " + bad
assert "Figure 3" not in main_part and "Figure 4" not in main_part

n_main_tables = len(re.findall(r"^### Table \d+\.", main_part, flags=re.M))
n_main_figs = len(re.findall(r"^### Figure \d+\.", main_part, flags=re.M))
assert (n_main_tables, n_main_figs) == (6, 2), (n_main_tables, n_main_figs)

kw_line = [l for l in text.split(NL) if l.startswith("**Keywords:**")][0]
n_kw = kw_line.count(";") + 1
assert 4 <= n_kw <= 6, n_kw

title = text.split(NL, 1)[0].lstrip("# ").strip()
n_title = len(title.split())
assert n_title <= 17, n_title

refs = text[text.index("## References"):]
n_refs = len([l for l in refs.split(NL) if re.match(r"^\d+\. ", l)])
assert n_refs <= 40, n_refs

body = text[text.index("## Introduction"):text.index("## Main tables")]
tokens = len(re.findall(r"[^\s]+", body))
tokens_nohdr = len(re.findall(r"[^\s]+", NL.join(
    l for l in body.split(NL) if not l.startswith("#"))))
abstract = text.split("## Abstract",1)[1].split("## Introduction",1)[0].split("**Keywords:**")[0]
print("title words          :", n_title)
print("keywords             :", n_kw)
print("abstract (body)      :", len(re.findall(r"[^\s]+", abstract)))
print("body words (all)     :", tokens, " / 4500")
print("body words (no hdr)  :", tokens_nohdr)
print("main tables/figures  : %d / %d" % (n_main_tables, n_main_figs))
print("references           :", n_refs)
print("total lines          :", text.count(NL) + 1)
assert tokens <= 4500, "OVER LIMIT"

with open(PATH, "wb") as f:
    f.write(text.encode("utf-8"))
print("WRITTEN OK")
