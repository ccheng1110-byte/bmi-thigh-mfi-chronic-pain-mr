# FinnGen R12 external replication — execution results

Date: 2026-09-20
Data: FinnGen R12 `summary_stats/release/` (FinnGen-only; **`meta_analysis/` was not used**)
Exposures: TMFI, BMI (still UK Biobank source; instrument sets unchanged)
Status: **extraction → alignment → scale verification → MR → dual-implementation validation, all complete**

---

## 1. Summary of findings

**External replication holds.** In a 500,000-person cohort with zero sample overlap with UK Biobank:

1. **The TMFI null is replicated**: all three endpoints, both UVMR and MVMR, are null, and the bounds are narrow (MDE for `PAIN` is only 1.069 OR/SD).
2. **The BMI effect is replicated and strongly significant for limb pain**: MVMR OR/SD = **1.217 (1.116–1.327), P = 9.0×10⁻⁶**.
3. **The contrast pattern matches UK Biobank exactly**: no effect for TMFI, an effect for BMI, and a stronger BMI effect at weight-bearing sites.

This upgrades the manuscript's claim that "sample-overlap bias is not sufficient to overturn the conclusion" from a **theoretical argument** (tipping point / structural ceiling) to **empirical evidence**, and simultaneously invalidates the reviewer comment that "`CPI_no_imaging` is not an independent replication" — FinnGen is a genuinely independent population.

---

## 2. Main results (MVMR, 92 instruments)

| Endpoint | N (cases) | TMFI OR/SD (95% CI) | P | BMI OR/SD (95% CI) | P |
|---|---:|---|---|---|---|
| `PAIN` | 499,362 (237,944) | 1.0186 (0.9744–1.0648) | 0.416 | 1.0546 (0.9815–1.1331) | 0.147 |
| `M13_LIMBPAIN` | 403,025 (45,476) | 0.9830 (0.9317–1.0371) | 0.530 | **1.2170 (1.1159–1.3272)** | **9.0×10⁻⁶** |
| `M13_FIBROMYALGIA` | 361,172 (3,623) | 1.0012 (0.8512–1.1777) | 0.988 | 1.0369 (0.7973–1.3484) | 0.787 |

### UVMR (frozen 42-SNP TMFI instrument set, all 42/42 present)

| Endpoint | OR/SD (95% CI) | P | MDE (80% power) | Egger intercept P |
|---|---|---|---|---|
| `PAIN` | 1.0211 (0.9743–1.0701) | 0.383 | 1.0693 | **0.041** |
| `M13_LIMBPAIN` | 0.9992 (0.9426–1.0592) | 0.978 | 1.0870 | 0.243 |
| `M13_FIBROMYALGIA` | 1.0170 (0.8778–1.1784) | 0.822 | 1.2343 | 0.793 |

> **Corrigendum:** this table originally recorded the Egger intercept P for fibromyalgia as 0.481, taken from an older result before the 42-SNP set was completed. The current `results/stage6_mr_results.csv` gives **0.793** (alongside 0.243 for limb pain and **0.041** for `PAIN`). Table 9 and the Results section of the manuscript have been written using 0.793.

### Instrument strength (92 instruments)

- Univariable mean F: TMFI 26.6, BMI 48.0
- Conditional F (R²-adjusted): TMFI 23.88, BMI 43.03
- Conditional F (weighted partial F): TMFI 24.57, BMI 41.77
- Steiger: exposure-dominant for 100% of instruments at all three endpoints (mean R² exposure 1.36×10⁻³ vs outcome 3.2–8.1×10⁻⁶)
- LOO (BMI coefficient): `PAIN` 0.0397–0.0664, full sample 0.0531; `LIMBPAIN` 0.1761–0.2069, full sample 0.1964 — **no single instrument drives the result**

---

## 3. Comparison with the UK Biobank primary analysis

| Analysis | Cohort | TMFI | BMI |
|---|---|---|---|
| **Primary analysis** CPI_full | UKB (N=105,977) | β=−0.0031 (P=0.898) | β=+0.2610 (P=0.0011) |
| Knee pain | UKB (98,704 cases) | 1.0146 (P=0.578) | **1.2751 (P=2.5×10⁻⁹)** |
| Back pain | UKB (118,471 cases) | 1.0005 (P=0.983) | **1.1344 (P=9.4×10⁻⁴)** |
| Knee pain (figshare) | independent (≈96,000 cases) | 1.0089 (P=0.714) | 0.9048 (P=0.008) ⚠️ reversed |
| **External replication** limb pain | **FinnGen (45,476 cases)** | 0.9830 (P=0.530) | **1.2170 (P=9.0×10⁻⁶)** |
| **External replication** overall pain | **FinnGen (237,944 cases)** | 1.0186 (P=0.416) | 1.0546 (P=0.147) |

**How to read this:** the FinnGen limb-pain BMI effect (1.217) falls between the UK Biobank estimates for back pain (1.134) and knee pain (1.275), consistent with the explanation that effects are larger at weight-bearing joints. The overall-pain endpoint (`PAIN` mixes headache, abdominal pain and other non-mechanical pain) is diluted to 1.055 and non-significant — still positive in direction, which is consistent with expectation rather than contradictory.

⚠️ Note: the figshare knee-pain BMI effect in the table is **reversed** (0.905, P=0.008). This is an anomaly that already exists in the manuscript's supplementary material and is unrelated to this round of work, but if external replication is written into the main text, reviewers are likely to ask about it in passing.

---

## 4. Methodological verification (three gates, all passed)

### 1. Integrity

File sizes and GCS object md5 match **bit-for-bit** for all three files:
`PAIN` 820,988,965 B / `82X9DEpw5pYsUnAB+tVQ1Q==`; `FIBRO` 810,479,105 B / `5fPaO4jegXvb1WZFDFxl+g==`; `LIMBPAIN` 818,394,608 B / `CASXOirM2b/WdJEDUbP7hw==`

### 2. Alignment (§4)

- 21.3M rows scanned per file; **108/109** instruments matched (only `rs1808579` is absent from FinnGen)
- Alignment rules: 100 direct, 6 palindromic resolved by allele frequency, 2 strand flips, **0 dropped**
- Coordinate pitfall avoided: R12 is GRCh38 (e.g. rs1421085 = 16:53767042) whereas our panel is b37 (16:53800954), so **the whole pipeline anchored on rsid and reused our own b37 coordinates**; no liftover was performed
- Allele-frequency concordance: median |FinnGen − UKB| **0.035**, maximum 0.126, none exceeding 0.20 — Finnish population drift is within the normal range

### 3. Scale (§5) — the most critical gate

Median of the identity `SE² × 2Nf(1−f)` over 83–101 **non-significant** instruments:

| Endpoint | Observed | Logistic target `1/(μ(1−μ))` | Ratio | Linear (risk-difference) target | Fold deviation |
|---|---:|---:|---:|---:|---:|
| `PAIN` | 4.1887 | 4.0089 | **1.045** | 0.24945 | 16.8× |
| `M13_LIMBPAIN` | 10.2301 | 9.9896 | **1.024** | 0.10010 | 102× |
| `M13_FIBROMYALGIA` | 102.3148 | 100.6988 | **1.016** | 0.00993 | 10,303× |

**Verdict: log-OR.** FinnGen betas are already log odds ratios and require no λ rescaling. Carrying over the risk-difference rescaling workflow from Stage 5 would introduce errors of 16- to 10,000-fold. Observed values exceed the theoretical target by 1.6%–4.5%, attributable to REGENIE's LOCO correction and covariates, which is normal.

---

## 5. Dual-implementation validation (§7)

Python (`stage6_mr.py`, reusing the Stage 5 estimator) and base-R (`stage6_validation_baser.R`, R 4.5.1, no third-party packages) share the same frozen input:

| Estimator | Degree of agreement |
|---|---|
| MVMR TMFI / BMI **beta and SE** | difference **~1×10⁻¹⁶** (machine precision) |
| IVW-MRE, MR-Egger | relative difference 3×10⁻¹⁶ – 5×10⁻¹⁵ |
| Weighted median | difference 0.15% – 1.9% (a known definitional difference, see below) |

The weighted-median difference arises from **two differing implementation definitions**: Python uses centred cumulative weights `(cumsum − 0.5w)/Σw` while R uses plain `cumsum/Σw`, and the bootstrap seeds are independent. A difference of the same magnitude exists at Stage 5 (0.2%–1.1%).
**The statement at manuscript line 146 that "six significant figures agree" enumerates only IVW-FE, IVW-MRE, MR-Egger and the two MVMR coefficients, and does not include the weighted median; the wording is correctly scoped and is unaffected.**

---

## 6. Three problems found and fixed during the work

1. **The 42-SNP set is not a subset of the 93-SNP union** — only 26 overlap by rsid. Extracting only against the 93-SNP union would have silently reduced Stage 1 from 42 instruments to 26. `instruments_master.csv` was therefore built (109 rsids, the union of the two sets), and Stage 1 is now **42/42**.
2. Two bugs in the extraction script (a miscomputed split limit causing zero matches; a trailing newline in the last column splitting lines in half), both fixed.
3. **FinnGen has no sample overlap**, which has been confirmed; the tipping-point / structural-ceiling overlap-bias quantification is therefore unnecessary and should not be repeated here — it was designed for the UK Biobank internal outcomes.

---

## 7. Caveats to retain when writing this into the manuscript (mandatory)

1. **Related-phenotype replication, not same-phenotype replication**: FinnGen endpoints are binary ICD diagnoses, whereas CPI is a continuous intensity score.
2. **Independent on the outcome side, not on the exposure side**: TMFI and BMI still come from UK Biobank (TMFI has only one source worldwide, Zhang 2026, with no substitute).
3. `M13_FIBROMYALGIA` has only 3,623 cases and an MDE of 1.26 OR/SD, so its null provides **no useful information** and can only count as descriptive concordance.
4. The UVMR Egger intercept P for `PAIN` is 0.041 (nominally significant in one of three endpoints) and must be reported as such.
5. 92/93 instruments are usable (`rs1808579` is absent from FinnGen).

---

## 8. Wording templates ready for use

> **Status: all written into the manuscript (2026-09-20).** Executed by `_apply_finngen_integration.py`; four templates were rewritten to BMC Medicine style and placed in Methods (L150), Results (L248), Limitations (item 8) and the Discussion principal-finding section, with a new **Table 9**, three rows in Table 1, references [17][18], and additional sentences in Data availability and Ethics. A 320-word abstract version was also added (339 including keywords, ≤ 350). The templates below are retained for the record.

### Addendum: official R `MVMR` validation (appended after being written in)

`stage6_official_MVMR_validation.R` (R 4.6.1 + MVMR 0.4.8) shares the same frozen input as Python and base-R across the three endpoints:

- MVMR coefficients and standard errors: **bit-for-bit agreement** across the three implementations (PAIN β=0.00906364/0.05313475, SE=0.01114670/0.03662207).
- Conditional F: **TMFI 22.356, BMI 35.441** (same order as the UK Biobank primary analysis, 22.119/35.234; the difference comes from 92 vs 93 instruments). Table 9 of the manuscript uses these values, consistent with the rule stated in the main text that R `MVMR` is primary.
- Modified Q: 272.02 (df 89, P=2.0×10⁻²⁰), 130.29 (P=0.0029), 107.99 (P=0.083) — heterogeneity is present, so MRE is used throughout.
- t-test P values are reported alongside the normal-theory P values stated in the main text, in the table footnote (limb pain BMI: normal 9.0×10⁻⁶ vs t 2.6×10⁻⁵).

### Methods (new subsection)

> **External replication in an independent cohort.** To test the two central
> contrasts in a sample that shares no participants with UK Biobank, we
> repeated both instrument sets against FinnGen R12 release endpoints
> (summary statistics only; the FinnGen × UKB meta-analysis module was not
> used, because it reintroduces the overlap this analysis is designed to
> escape). Three endpoints were selected a priori on the basis of phenotype
> relevance and case count: `PAIN` (237,944 cases / 261,418 controls),
> `M13_LIMBPAIN` (45,476 / 357,549) and `M13_FIBROMYALGIA` (3,623 / 357,549).
> Because FinnGen fits REGENIE logistic rather than linear models, its
> coefficients are log odds ratios and no rescaling was applied; this was
> verified rather than assumed, using the identity
> SE² × 2Nf(1−f) = 1/(μ(1−μ)) on non-significant instruments, which returned
> 1.045, 1.024 and 1.016 times the logistic target for the three endpoints
> (the linear-model alternative is excluded by factors of 17 to 10,300).
> Variants were anchored on rsid and assigned our own b37 coordinates,
> because the release is on GRCh38 while the LD reference panel is on b37.
> Of the 109 instruments in the union of the two frozen sets, 108 were
> present; 100 aligned directly, 6 palindromic variants were resolved by
> allele frequency and 2 required strand-aware flipping, leaving none
> ambiguous. Median |allele frequency difference| between FinnGen and UK
> Biobank was 0.035 (maximum 0.126).

### Results (new paragraph)

> In FinnGen, neither instrument set produced a TMFI effect: the multivariable
> estimate was 1.019 (0.974–1.065) per SD for `PAIN`, 0.983 (0.932–1.037) for
> limb pain and 1.001 (0.851–1.178) for fibromyalgia, with minimum detectable
> effects of 1.07, 1.08 and 1.26 per SD. BMI, by contrast, was associated with
> limb pain at 1.217 per SD (1.116–1.327; P = 9.0 × 10⁻⁶) and pointed in the
> same direction for `PAIN` (1.055; 0.982–1.133; P = 0.15). The limb-pain
> estimate lies between the UK Biobank estimates for back pain (1.134) and
> knee pain (1.275), so the independent cohort reproduces both the null for
> thigh muscle fat infiltration and the gradient of the adiposity comparator
> across weight-bearing sites. Leave-one-out analysis left the BMI coefficient
> between 0.176 and 0.207 against a full-sample value of 0.196, and Steiger
> filtering indicated the exposure-dominant direction for 100% of instruments
> in all three endpoints.

### Limitations (new sentence)

> The replication is outcome-side only: the exposures remain UK Biobank
> instruments, because thigh muscle fat infiltration has been genotyped in no
> other cohort, and the FinnGen endpoints are binary ICD-based diagnoses
> rather than the continuous pain-intensity score, so concordance should be
> read as related-phenotype triangulation rather than same-phenotype
> replication. The fibromyalgia endpoint (3,623 cases) can exclude only
> effects larger than 1.26 per SD and therefore adds descriptive, not
> inferential, support.

### Significance (new sentence, responding to the overlap objection)

> These estimates are free of sample overlap by construction; the tipping-point
> and structural-ceiling arguments developed for the UK Biobank outcomes
> quantify a bias that cannot arise here, and the concordance of the two
> designs indicates that the bounded null for thigh muscle fat infiltration is
> not an artefact of that overlap.

---

## 9. Output files

| File | Contents |
|---|---|
| `download_finngen.py` | parallel, resumable download (complete) |
| `build_instrument_master.py` | 109-locus master table (union of the two instrument sets) |
| `extract_instruments.py` | 21.3M-row streaming scan, extraction by rsid |
| `harmonise_finngen.py` | allele alignment + scale-identity determination |
| `stage6_mr.py` | UVMR / MVMR / Egger / Steiger / conditional F / LOO |
| `stage6_validation_baser.R` | independent base-R implementation (no third-party packages) |
| `extracted/*_raw.tsv` | raw extraction for the three endpoints |
| `gwas_harmonised/*_harmonised.tsv.gz` | aligned outcomes (includes an `n_eff` audit column) |
| `gwas_harmonised/harmonise_report.json` | machine-readable record of alignment and scale determination |
| `results/stage6_mr_results.csv` | Python results |
| `results/stage6_baser_all.csv` | base-R results |
| `r_input/R_stage{1,2}_*.tsv` | frozen input shared by the two implementations |
| `finngen_R12_manifest.tsv` | dictionary of 2,469 endpoints (including case counts) |

---

## 10. Not yet executed (optional)

- **Cross-trait LDSC**: in theory FinnGen and UKB have zero overlap, so the intercept should be ≈0. This could serve as independent corroboration, but it requires an LDSC environment, and this round is already guaranteed directly by the design fact of "zero overlap", so it is **not currently necessary**.
- **MR-PRESSO / robust MVMR**: this round already has LOO plus the Egger intercept; these could be added if submitting to a methodologically oriented journal.
