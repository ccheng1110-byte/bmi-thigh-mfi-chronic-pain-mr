# Originality package: nearest-neighbour comparison table + main contribution statement

**Manuscript:** `Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md`
**Compiled:** 2026-09-20
**Purpose:** (i) to answer "where is the increment" before submission; (ii) for use in the opening paragraph of the cover letter and in positioning the Discussion; (iii) for submission as Supplementary Table S13.
**Verification:** all comparator studies were verified by online search on 2026-09-20 (PMID / DOI / volume-issue-pages given in the table); no literature data were taken from memory.

---

## 1. Main contribution statement (one sentence)

### Version A — for the opening paragraph of the cover letter (recommended, 48 words)

> Using body mass index as an internal comparator fitted in the same multivariable Mendelian randomisation model, this study is the first to bound the effect of a muscle-compartment exposure on chronic pain intensity: genetically predicted thigh muscle fat infiltration has no effect larger than 0.048 SD per SD — about one third of the adiposity effect the identical model recovered — and both the null and the adiposity gradient reproduce in FinnGen, a cohort sharing no participants with UK Biobank.

### Version B — academic / technical register (for the Abstract or the final sentence of the Introduction, 44 words)

> We show that, after conditioning on systemic adiposity, genetically predicted thigh muscle fat infiltration has no detectable effect on chronic pain intensity above 0.048 SD per SD (80% power), converting a muscle–pain null into a bounded exclusion rather than a non-significant result, with the null and the BMI gradient both reproduced in an independent, non-overlapping cohort.

### Version C — minimal version (for oral use or a response letter, 20 words)

> A bounded null: thigh muscle fat infiltration does not drive chronic pain intensity at more than one third of the scale on which adiposity does.

### Version D — three-part "Research in context" (if the target journal requires that section)

- **What is already known:** Genetically predicted adiposity causes site-specific pain (back, knee, hip) [19,20]; muscle phenotypes, including sarcopenia and lower-extremity muscle traits, show genetic overlap and MR associations with pain and spinal disease [6,7]; the genetic architecture of muscle fatty infiltration has recently been characterised, with muscle fatty infiltration treated as an outcome [21].
- **What this study adds:** A muscle-compartment exposure and an adiposity exposure estimated jointly, so that the null for the muscle compartment is read against (i) an internal comparator recovered in the same model, (ii) an explicit minimum detectable effect of 0.048 SD per SD, (iii) a genetically quantified overlap with a tipping-point bound, and (iv) replication in a cohort with zero participant overlap — none of which any of the nearest published analyses provides.
- **Clinical and research implication:** Thigh muscle fat infiltration should not be advanced as an intervention target for chronic pain intensity at anywhere near the scale of systemic adiposity; the next decisive test is joint conditioning of several correlated muscle compartments, not another univariable screen.

---

## 2. Nearest-neighbour comparison table (8 rows, including this manuscript's own positioning)

Column definitions:
- **Comparator in same model** = whether a comparator exposure "that should be detected" was estimated in the same (multivariable) model;
- **Power bound** = whether the minimum detectable effect (MDE) of the design at 80% power is reported;
- **Overlap quantified** = whether exposure–outcome sample overlap is quantified with a method such as cross-trait LDSC (rather than merely declaring "two-sample");
- **Independent replication** = whether there is independent outcome replication with **zero overlap** with the exposure-source cohort.

| # | Study (verified reference) | Exposure → Outcome | Design / estimator | Comparator in same model | Power bound (MDE) | Overlap quantified | Independent replication with zero overlap | Main result | Increment of this manuscript relative to it |
|---|---|---|---|---|---|---|---|---|---|
| 1 | Elgaeva EE, Tsepilov Y, Freidin MB, Williams FMK, Aulchenko Y, Suri P. *Eur Spine J*. 2020;29(4):686–691. doi:10.1007/s00586-019-06224-6. PMID 31797139 | BMI (GIANT, 60 SNP) → back pain / chronic back pain (UKB n=453,860) | Two-sample UVMR, IVW primary | No (single exposure) | No | **Avoided by design** (deliberately non-overlapping samples) | No | OR/SD 1.15 (1.06–1.25) back pain; 1.20 (1.09–1.32) chronic back pain | That paper removed overlap by "switching cohort"; this manuscript **cannot switch cohort** (TMFI and CPI are both in UKB), so it instead applies a triple treatment: quantification + tipping point + FinnGen replication. Its BMI effect sizes (back pain 1.15–1.20) are of the same order as this manuscript's back-pain estimate of 1.134, and the two corroborate each other |
| 2 | Chen X, Tang H, Lin J, Zeng R. *Front Endocrinol (Lausanne)*. 2022;13:971997. doi:10.3389/fendo.2022.971997. PMID 36082069 | BMI / WC / HC / WHR → knee, hip, back pain | Two-sample UVMR (69/42/52/29 SNP), IVW+WM+Egger | No (four exposures run **separately**, univariable, not jointly) | No | Not reported | No | BMI: knee 1.049, hip 1.034, back 1.022 | Four adiposity measures each run univariably cannot answer "local muscle vs systemic fat"; this manuscript puts both exposures into the same model, yielding a **conditional direct effect** rather than a total effect |
| 3 | Wang Q, Lu S. *Eur Spine J*. 2026. doi:10.1007/s00586-026-10184-z (ref 7 of this manuscript) | Lower-limb muscle traits → spinal disease | MR | No | No | Not reported | No | Reports positive genetically predicted association | This positive muscle–pain MR does **not adjust for adiposity**, so it cannot separate "the effect of muscle itself" from "correlation carried by body fat". That is exactly the step this manuscript performs, and it yields a null |
| 4 | Li J, Huang J, Han J. *Medicine*. 2026;105(20):e48819. doi:10.1097/MD.0000000000048819. PMID 42152289 (ref 6 of this manuscript) | Sarcopenia ↔ chronic pain (bidirectional) | Bidirectional UVMR + Bayesian colocalisation + drug repositioning | No | No | Not reported | No | Significant in both directions; 9 shared genes | **Positive** evidence for the same "muscle–pain" proposition. This manuscript is complementary and more conservative: after adjusting for adiposity, the thigh muscle fat compartment is undetectable; the two papers together show that "muscle" is not a single exposure |
| 5 | Thompson WD, Swain S, Zhao SS, Coupland C, Kuo CF, Doherty M, Zhang W. *Pain*. 2024;165(8):1882–1889. doi:10.1097/j.pain.0000000000003183. PMID 38358931 | BMI + multisite chronic pain (MCP) → knee OA | **MVMR** (501 SNP), modified IVW + MVMR-Egger | Partial: the two exposures act as each other's comparator (both significant) | No | **Not addressed** (the paper states "there is sample overlap between the 3 GWAS") | No | BMI 1.76 (1.64–1.88), MCP 1.83 (1.54–2.16), each independent | **Design closest to this manuscript** (musculoskeletal pain outcome + correlated dual-exposure MVMR). Difference: there both exposures are significant, so no power problem for a null arises; overlap is merely stated, not quantified. This manuscript adds exactly those two pieces: MDE + cross-trait LDSC tipping point |
| 6 | Norris T, Sanderson E, Cooper R, Garfield V, Pinto Pereira SM, et al. *Sci Rep*. 2023;13:16886. doi:10.1038/s41598-023-43908-y | BMI / WHR (+CRP mediation) → grip strength | MVMR (IVW/WM/Egger) | No | No | Mainly a single cohort | No | CRP mediation is **zero** | Methodological precedent: using MVMR to report a null of the form "this pathway does not mediate". But no MDE, no comparator exposure, no independent replication — the null in this manuscript adds three further layers of quantification |
| 7 | Tan J, Ding Z, Zheng J, Zhang J, Chen X, Li Z, Shi L, Chen J, Sun Y. *Int J Surg*. 2025. doi:10.1097/JS9.0000000000003523. PMID 41056018 | Metabolites / proteins / genes (eQTL, pQTL) → **muscle fat infiltration (MFI)** (33,300 UKB MRI) | GWAS (91 SNP / 26 locus) + MR + LDSC + wet-lab validation of KLF5 | No | No | Not applicable | No | Genetic architecture of MFI; KLF5 as a regulator | **Same source phenotype as this manuscript but opposite direction**: that paper treats MFI as the **outcome** and never as an **exposure** for pain. This is the first study to use a muscle fat compartment as a pain exposure with adiposity adjustment. **Cite both together, to avoid appearing to have scooped it** |
| 8 | Su C-Y, Butler-Laporte G. *medRxiv*. 2026. doi:10.64898/2026.05.29.26354472 (preprint) | BMI (**positive-control exposure**) → disease mortality | Two-sample MR with a positive-control exposure for **power calibration** | **Yes** (BMI explicitly used as positive-control exposure) | **Yes** (power calibration is the whole point of the paper) | Not applicable | No | Mortality GWAS underpowered, so "limited overlap" cannot be read biologically | Conceptual precedent, showing that "comparator exposure + power calibration" is a practice only beginning to be used explicitly in 2026. This manuscript goes further: it explicitly **denies** that the comparator exposure is a positive control with a known answer, claiming only that it is sensitivity evidence from the same model |
| **This manuscript** | TMFI + BMI → CPI (two UKB definitions) + site-specific pain + three FinnGen endpoints | **MVMR** (93 SNP) + UVMR (42 SNP) + LDSC + tipping point / structural ceiling + dual implementations | **Yes** (BMI as internal comparator, explicitly not claimed to be a positive control with a known answer) | **Yes** (0.048 SD/SD at 80% power, MVMR scale) | **Yes** (cross-trait LDSC intercept + tipping point vs structural ceiling) | **Yes** (FinnGen R12, zero UKB overlap) | TMFI null throughout; BMI 0.090 SD/SD (P=0.00075); BMI 1.217 for FinnGen limb pain | — | — |

---

## 3. Four pillars of originality (each mapped to a verifiable element of the manuscript)

| Pillar | Specific claim | Where it lands in the manuscript | Nearest-neighbour comparison |
|---|---|---|---|
| **P1 First test of a muscle compartment conditional on adiposity** | Published muscle–pain MR studies [6,7] are all single-exposure; their positive results cannot exclude "carried by body fat" | Table 3; Discussion "Relation to existing evidence" | rows 3, 4 |
| **P2 First use of a muscle fat compartment as a pain exposure** | The existing MFI GWAS [21] treats MFI as an outcome | Table 3, Table 8, Table 9 | row 7 |
| **P3 Bounded null (comparator + MDE)** | The null is accompanied by: a comparator exposure detected in the same model, plus a minimum detectable effect of 0.048 SD/SD at 80% power | Methods "Statistical power"; Table 7; Discussion "Principal finding" | rows 1–6 have no MDE; row 8 has a comparator but no compartment-specific test with a matching MDE |
| **P4 Overlap upgraded from argument to evidence** | Cross-trait LDSC tipping point (theoretical) + FinnGen zero-overlap replication (empirical) | Table 7 + Table 9 + two Methods subsections | row 1 avoids overlap by switching cohort; row 5 states it is unaddressed; this manuscript takes both routes |

---

## 4. Three "this is not novel" objections to disclose proactively, with responses

1. **Objection: "BMI → pain is already known; this paper merely repeats it."**
   → Response: yes, and that is precisely the design intent — BMI is the **internal comparator**, not the finding this paper claims. What the paper claims is that "in the same model that does detect the comparator effect, TMFI is excluded above 0.048 SD/SD". This is already written into the first paragraph of the Discussion and into Limitations (the comparator is not a positive control with a known answer).

2. **Objection: "A null result — no new locus, gene or target."**
   → Response: this paper is positioned as an **exclusion** study rather than a discovery study, so journals that require a new finding (e.g. PAIN) are a poor match; a methodologically or epidemiologically oriented journal, or one specialised in muscle / body composition, is more appropriate (see `投稿期刊选择与成功率建议_TMFI_BMI_CPI_2026-09-19.md`; original filename retained verbatim). The currency of an exclusion study is the **bound**, and this paper gives a verifiable numeric bound.

3. **Objection: "TMFI and CPI are both in UKB, so overlap bias cannot be excluded."**
   → Response: three independent treatments — (a) cross-trait LDSC intercept 0.052, induced bias ≤0.0015 SD/SD; (b) tipping point 1.03–1.28, all above the structural ceiling of 0.66–0.77; (c) **FinnGen R12 with zero participant overlap**, in which both comparators reproduce. The third is new in this round and is what upgrades (a) and (b) from argument to evidence.

---

## 5. Four changes recommended for the manuscript (already applied in the same batch)

1. Add a positioning sentence to the Introduction naming [19][20][21], stating that "what is missing is the joint test";
2. Add two sentences at the end of Discussion "Relation to existing evidence", fixing the division of labour with the adiposity–pain literature and pointing to Supplementary Table S13;
3. Add one sentence to Discussion "Implications" citing [22] (preprint), noting that comparator-exposure power calibration is an emergent practice in 2026;
4. Add references [19] Elgaeva 2020, [20] Chen 2022, [21] Tan 2025, [22] Su & Butler-Laporte 2026; add **S13** to the supplementary list (i.e. the table in Section 2 of this file).

---

## 6. Still open (requires an author decision)

- Whether to keep [22] (medRxiv preprint) as a reference: keeping it strengthens the argument that "comparator-exposure power calibration has precedent"; removing it leaves a reference list consisting entirely of peer-reviewed work.
- Whether to submit this table with the manuscript (as S13): submission is recommended — the leading reason for desk rejection at BMC Medicine is "insufficient sense of increment", and this table answers that question head-on.
- If the manuscript is ultimately redirected to IJE or JCSM, the "increment" wording of this table needs rewriting for that journal's positioning (IJE emphasises methodology, JCSM emphasises the body-composition theme).
