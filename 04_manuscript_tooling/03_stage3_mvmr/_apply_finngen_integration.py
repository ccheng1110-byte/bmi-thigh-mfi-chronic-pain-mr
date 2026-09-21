# -*- coding: utf-8 -*-
"""
Integrate the FinnGen R12 external replication into the manuscript.

Design rules (same convention as every previous edit script in this project):
  * every text substitution asserts that the anchor occurs exactly once;
  * insertions are anchored on an exact heading line, never on a line number;
  * the abstract word count is checked BEFORE the file is written;
  * a .bak copy is written first; line endings (CRLF) are restored on write.
"""
import io
import os
import re
import shutil
import sys

PATH = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"
BAK = PATH + ".bak_finngen2026-09-20"

raw = io.open(PATH, "rb").read()
text = raw.decode("utf-8")
LF = "\r\n" in text
lines = text.replace("\r\n", "\n").split("\n")
original_line_count = len(lines)

errors = []


def rep(old, new, expect=1, tag=""):
    """Replace an exact substring; assert the anchor occurs `expect` times."""
    n = text.count(old)
    if n != expect:
        errors.append("ANCHOR[%s] found %d times, expected %d: %r" % (tag, n, expect, old[:70]))
        return False
    text_local = text.replace(old, new)
    globals()["text"] = text_local
    return True


# ---------------------------------------------------------------- 1. Abstract
rep(
    "with standard diagnostics for instrument strength and pleiotropy. "
    "Because exposure and outcomes come from the same cohort, that overlap was "
    "quantified genetically and carried through every estimate; results are "
    "reported against the smallest effect detectable at 80% power rather than significance.",
    "with standard diagnostics. Exposure–outcome overlap was quantified genetically "
    "and carried through every estimate; results are reported against the smallest "
    "effect detectable at 80% power rather than significance. Both stages were then "
    "repeated against three FinnGen endpoints (up to 499,362 participants) with no UK "
    "Biobank overlap.",
    tag="abstract-methods",
)

rep(
    "Exposure and outcomes share participants (intercept 0.052, P=1.2×10⁻¹¹); this "
    "moved estimates by ≤0.0015 SD per SD, and no degree of sharing could reverse "
    "either finding. Two site-specific secondary outcomes also returned a TMFI null "
    "(knee pain OR 1.015 per SD, P=0.58; back pain 1.001, P=0.98), with a descriptively "
    "larger BMI effect at the weight-bearing joint (1.275 vs 1.134).",
    "Exposure–outcome overlap (intercept 0.052, P=1.2×10⁻¹¹) moved estimates by "
    "≤0.0015 SD per SD and could not reverse either finding. Two site-specific outcomes "
    "also returned a TMFI null (knee 1.015 per SD, P=0.58; back 1.001, P=0.98), with a "
    "descriptively larger BMI effect at the weight-bearing joint (1.275 vs 1.134). Both "
    "patterns held in FinnGen: TMFI null (limb pain OR 0.983 per SD, P=0.53), BMI "
    "associated with limb pain (1.217, 1.116–1.327, P=9×10⁻⁶).",
    tag="abstract-results",
)

rep(
    "Against an explicit power bound, an internal comparator that behaved as expected "
    "and quantified participant overlap, this is a bounded null:",
    "Against an explicit power bound, an internal comparator that behaved as expected, "
    "replication in an independent cohort and quantified overlap, this is a bounded null:",
    tag="abstract-conclusions",
)

# ------------------------------------------------------- 2. Table 1 data rows
rep(
    "| Stage 2 instruments | Same thresholds applied to the union of TMFI and BMI candidates | — | — | 93 SNPs | Multivariable MR |",
    "| Stage 2 instruments | Same thresholds applied to the union of TMFI and BMI candidates | — | — | 93 SNPs | Multivariable MR |\n"
    "| FinnGen `PAIN` | FinnGen R12 release, FinnGen-only summary statistics (`summary_stats/release/`) [17] | Finnish | 499,362 | 237,944 cases / 261,418 controls | External replication outcome |\n"
    "| FinnGen `M13_LIMBPAIN` | FinnGen R12 release, FinnGen-only summary statistics (`summary_stats/release/`) [17] | Finnish | 403,025 | 45,476 cases / 357,549 controls | External replication outcome |\n"
    "| FinnGen `M13_FIBROMYALGIA` | FinnGen R12 release, FinnGen-only summary statistics (`summary_stats/release/`) [17] | Finnish | 361,172 | 3,623 cases / 357,549 controls | External replication outcome |",
    tag="table1-rows",
)

rep(
    "All four source GWAS were restricted to European-ancestry participan",
    "All source GWAS were restricted to European-ancestry participan",
    tag="table1-note",
)

rep(
    "all four source GWAS were conducted in European-ancestry participants",
    "all source GWAS were conducted in European-ancestry participants",
    tag="limitations-ancestry",
)

# --------------------------------------------------------- 3. Limitations
rep(
    "Seven boundaries qualify the estimates.",
    "Eight boundaries qualify the estimates.",
    tag="limitations-count",
)

rep(
    "Each of these specifies what a replication study must fix, and none alters the primary contrast.",
    "Eighth, the external replication is outcome-side only: the exposures remain UK "
    "Biobank instruments, because thigh muscle fat infiltration has been genotyped in no "
    "other cohort, and the FinnGen endpoints are binary ICD-based diagnoses rather than "
    "the continuous pain-intensity score, so concordance between the two designs should be "
    "read as related-phenotype triangulation rather than same-phenotype replication; the "
    "fibromyalgia endpoint (3,623 cases) can exclude only effects larger than 1.26 per SD "
    "and therefore adds descriptive, not inferential, support. "
    "Each of these specifies what a replication study must fix, and none alters the primary contrast.",
    tag="limitations-eighth",
)

# ------------------------------------------------- 4. Discussion: significance
rep(
    "The upper confidence limits therefore exclude any positive TMFI effect larger than "
    "35% of the BMI point estimate for `CPI_full` and 28% for `CPI_no_imaging`.",
    "The upper confidence limits therefore exclude any positive TMFI effect larger than "
    "35% of the BMI point estimate for `CPI_full` and 28% for `CPI_no_imaging`. The same "
    "two contrasts were recovered in FinnGen, a cohort with no participant overlap at all "
    "(Table 9), so the bounded null is not an artefact of the exposure–outcome overlap "
    "that the tipping-point arguments address.",
    tag="discussion-significance",
)

# ----------------------------------------------- 5. Figure 1 legend / availability
rep(
    "the corresponding numbers for `CPI_no_imaging` are given in Table 3.",
    "the corresponding numbers for `CPI_no_imaging` are given in Table 3. The external "
    "replication stage (Results; Table 9) applies the identical frozen instrument sets to "
    "three FinnGen endpoints and is not shown separately in the figure.",
    tag="figure1-legend",
)

rep(
    "All datasets are publicly available; no restricted-access individual-level data were used.",
    "FinnGen R12 summary statistics were retrieved from the public FinnGen release [18]; "
    "the analysed files were `finngen_R12_PAIN.gz` (820,988,965 bytes; MD5 "
    "`82X9DEpw5pYsUnAB+tVQ1Q==`), `finngen_R12_M13_LIMBPAIN.gz` (818,394,608 bytes; MD5 "
    "`CASXOirM2b/WdJEDUbP7hw==`) and `finngen_R12_M13_FIBROMYALGIA.gz` (810,479,105 "
    "bytes; MD5 `5fPaO4jegXvb1WZFDFxl+g==`), all taken from `summary_stats/release/`; the "
    "FinnGen × UK Biobank meta-analysis module was not used, because it would reintroduce "
    "the participant overlap this replication is designed to exclude. All datasets are "
    "publicly available; no restricted-access individual-level data were used.",
    tag="data-availability",
)

rep(
    "No individual-level data were accessed, so no further ethical approval was required.",
    "No individual-level data were accessed, so no further ethical approval was required. "
    "The FinnGen R12 summary statistics were obtained from the publicly accessible "
    "FinnGen release under its data-access policy, and no individual-level Finnish data "
    "were accessed.",
    tag="ethics",
)

# ------------------------------------------------------------- 6. Insertions
# CRITICAL: every `rep()` above mutated `text` only. Re-sync the line list so the
# insertions below operate on the already-substituted text.
lines = text.replace("\r\n", "\n").split("\n")


def insert_before(anchor, block, tag):
    idx = None
    for i, l in enumerate(lines):
        if l.strip() == anchor:
            idx = i
            break
    if idx is None:
        errors.append("HEADING[%s] not found: %r" % (tag, anchor))
        return
    lines[idx:idx] = block.split("\n")


METHODS_BLOCK = """### External replication in a cohort with no participant overlap

Every argument above about sample overlap exists because the UK Biobank exposure and outcomes draw on the same participants. To test the two central contrasts where that bias cannot arise, we repeated both stages against three endpoints from the FinnGen R12 release [17], using the FinnGen-only summary statistics: the FinnGen × UK Biobank meta-analysis module was deliberately not used, because it reintroduces exactly the overlap this analysis is designed to escape. Endpoints were chosen a priori on phenotype relevance and case count — `PAIN` (237,944 cases, 261,418 controls), `M13_LIMBPAIN` (45,476 / 357,549) and `M13_FIBROMYALGIA` (3,623 / 357,549). The exposures are unchanged, because thigh muscle fat infiltration has been genotyped in no other cohort; the replication is therefore outcome-side only.

FinnGen fits a logistic mixed model (REGENIE), so its coefficients are log odds ratios, and — as with the site-matched outcomes — this was verified rather than assumed. The identity SE² × 2Nf(1−f) = 1/(μ(1−μ)) was evaluated on non-significant instruments only, because Firth-corrected standard errors at significant variants are back-derived from the P value and would confound the check: the median ratio to the logistic target was 1.045 for `PAIN` (83 variants), 1.024 for limb pain (93) and 1.016 for fibromyalgia (101), whereas the linear-model alternative is excluded by factors of 17 to 10,300. No rescaling was applied. This is the opposite of the site-matched conversion, where a linear model made rescaling necessary, and applying that correction here would have introduced errors of exactly that magnitude.

The release is on GRCh38 while the LD reference panel is on b37, so variants were anchored on rsid and assigned our own b37 coordinates rather than lifted over. Of the 109 instruments in the union of the two frozen sets, 108 were present (`rs1808579` is absent from FinnGen); 100 aligned directly, six palindromic variants were resolved on allele frequency and two required strand-aware flipping, leaving none ambiguous. Median |allele-frequency difference| between FinnGen and UK Biobank was 0.035 (maximum 0.126), consistent with ordinary Finnish drift. Stage 1 therefore ran on all 42 TMFI instruments and Stage 2 on 92 of the 93 joint instruments, and every estimator was computed with the same code as the preceding stages and, as before, twice — the Python pipeline and the base-R re-implementation sharing identical frozen input."""

RESULTS_BLOCK = """### External replication: both contrasts hold in a cohort with no participant overlap

Repeating Stage 1 and Stage 2 against three FinnGen endpoints reproduced both central findings (Table 9). No TMFI effect emerged: the multivariable estimate was 1.019 per SD (95% CI 0.974–1.065, P=0.42) for `PAIN`, 0.983 (0.932–1.037, P=0.53) for limb pain and 1.001 (0.851–1.178, P=0.99) for fibromyalgia, and the univariable estimates agreed (1.021, 0.999 and 1.017). The internal comparator again behaved as expected: BMI was associated with limb pain at 1.217 per SD (1.116–1.327; P=9.0×10⁻⁶) and pointed in the same direction, although not significantly, for `PAIN` (1.055; 0.982–1.133; P=0.15) and fibromyalgia (1.037; 0.797–1.348; P=0.79). The limb-pain estimate lies between the UK Biobank estimates for back pain (1.134) and knee pain (1.275), so the independent cohort reproduces not only the null for thigh muscle fat infiltration but also the gradient of the adiposity comparator across weight-bearing sites.

Instrument strength was comparable to the primary analysis — conditional F 22.4 for TMFI and 35.4 for BMI, against 22.1 and 35.2 for `CPI_full`, with univariable mean F of 26.6 and 48.0 — Steiger filtering indicated the exposure-dominant direction for 42 of 42 instruments in all three endpoints, and leave-one-out analysis left the limb-pain BMI coefficient between 0.176 and 0.207 against a full-sample value of 0.196. Heterogeneity was present (modified Q 3.03, 1.46 and 1.20), so the random-effects specification used throughout is the relevant one; the Egger intercept was nominally significant for `PAIN` alone (P=0.041, against 0.24 for limb pain and 0.79 for fibromyalgia), which is within expectation for three endpoints under the declared rule of no multiplicity correction. Minimum detectable effects for the univariable TMFI estimates were 1.07, 1.09 and 1.26 per SD, so the fibromyalgia endpoint, with 3,623 cases, contributes descriptive rather than inferential support.

Because no participants are shared, none of the overlap machinery of Table 7 applies to these estimates: the cross-trait intercept is zero by construction rather than by argument, and the concordance of the two designs indicates that the bounded null is not an artefact of that overlap."""

TABLE9_BLOCK = """### Table 9. External replication in FinnGen R12 (no participant overlap with UK Biobank)

| Endpoint | Cases / controls | Stage 1 TMFI OR per SD (95% CI) | P | MDE | Stage 2 TMFI OR per SD (95% CI) | P | Stage 2 BMI OR per SD (95% CI) | P |
|---|---|---|---|---|---|---|---|---|
| `PAIN` | 237,944 / 261,418 | 1.021 (0.974–1.070) | 0.383 | 1.069 | 1.019 (0.974–1.065) | 0.416 | 1.055 (0.982–1.133) | 0.147 |
| `M13_LIMBPAIN` | 45,476 / 357,549 | 0.999 (0.943–1.059) | 0.978 | 1.087 | 0.983 (0.932–1.037) | 0.530 | 1.217 (1.116–1.327) | 9.0×10⁻⁶ |
| `M13_FIBROMYALGIA` | 3,623 / 357,549 | 1.017 (0.878–1.178) | 0.822 | 1.234 | 1.001 (0.851–1.178) | 0.988 | 1.037 (0.797–1.348) | 0.787 |

**Table note:** Stage 1 uses the frozen 42-instrument TMFI set, all 42 of which were present in FinnGen; Stage 2 uses 92 of the 93 joint instruments (`rs1808579` is absent). Conditional F statistics from the R `MVMR` package were 22.36 (TMFI) and 35.44 (BMI) and univariable mean F was 26.6 and 48.0. P values follow the declared normal-approximation rule; the official `MVMR` t-test returns P = 0.418, 0.532 and 0.988 for the TMFI coefficients and 0.150, 2.6×10⁻⁵ and 0.788 for the BMI coefficients, and the multivariable coefficients and standard errors agree with the Python and base-R implementations to machine precision. Modified Q was 3.03, 1.46 and 1.20 and the Egger intercept P was 0.041, 0.243 and 0.793. FinnGen endpoints are binary ICD-based diagnoses, not the continuous pain-intensity score, so the concordance is related-phenotype triangulation rather than same-phenotype replication."""

REF_BLOCK = """17. Kurki MI, Karjalainen J, Palta P, et al. FinnGen provides genetic insights from a well-phenotyped isolated population. *Nature*. 2023;613:508–518. doi:[10.1038/s41586-022-05473-8](https://doi.org/10.1038/s41586-022-05473-8).
18. FinnGen. FinnGen documentation: data download (R12 release). 2025. Available from: https://finngen.gitbook.io/documentation/data-download. Accessed 20 September 2026."""

insert_before("### Reproducibility and statistical reporting", METHODS_BLOCK, "methods-block")
insert_before("## Discussion", RESULTS_BLOCK, "results-block")
insert_before("## Figure legends", TABLE9_BLOCK, "table9-block")

# references: append after reference 16
for i, l in enumerate(lines):
    if l.startswith("16. Skrivankova"):
        lines.insert(i + 1, REF_BLOCK)
        break
else:
    errors.append("HEADING[references] reference 16 not found")

if errors:
    sys.stdout.write("\n".join(errors) + "\n")
    sys.stdout.write("ABORTED — nothing was written.\n")
    sys.exit(1)

# ------------------------------------------------------------ 7. Verification
new_text = "\n".join(lines)


def wc(a, b):
    seg = "\n".join(new_text.split("\n")[a - 1 : b])
    seg = re.sub(r"^\|.*$", "", seg, flags=re.M)
    return len(re.findall(r"[A-Za-z][A-Za-z'\-–]*", seg))


# abstract = lines 9..25 (Background through the keywords line)
abs_i = next(i for i, l in enumerate(lines, 1) if l.strip() == "### Background")
kw_i = next(i for i, l in enumerate(lines, 1) if l.startswith("**Keywords:**"))
abs_words = wc(abs_i, kw_i)
body_words = wc(abs_i, kw_i - 1)
sys.stdout.write(
    "abstract body (Background..Conclusions) = %d words; including keywords = %d\n"
    % (body_words, abs_words)
)
if abs_words > 350:
    sys.stdout.write("ABORTED — abstract exceeds the 350-word limit.\n")
    sys.exit(1)

for pat in ["Table 9", "[17]", "[18]"]:
    sys.stdout.write("occurrences of %-9s : %d\n" % (pat, new_text.count(pat)))

if not LF:
    new_text = new_text
out = new_text.replace("\n", "\r\n") if LF else new_text
shutil.copy2(PATH, BAK)
io.open(PATH, "wb").write(out.encode("utf-8"))
sys.stdout.write("\nwritten: %d -> %d lines; backup %s\n" % (original_line_count, len(lines), os.path.basename(BAK)))
