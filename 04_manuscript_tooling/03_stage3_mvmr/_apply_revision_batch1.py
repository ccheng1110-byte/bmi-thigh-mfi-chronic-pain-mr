# -*- coding: utf-8 -*-
"""JCSM revision pass 1: batch 1-4 of the adjudicated revision list.

Batches applied here
  B1  stale table cross-references (5 sites)
  B2  logical / defensive sentences (L63, L73, L77, L81, L83)
  B3  FinnGen downgrade (L73, abstract, conclusion) + pre-specified -> fixed in advance
  B4  Limitations restructure, AI declaration, abstract MDE wording, L33 narrowing,
      Supplementary Methods flag wording
"""
import re, sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

P = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"
W = lambda s: len(re.findall(r"[^\s]+", s))

raw = open(P, 'rb').read().decode('utf-8')
text = raw.replace('\r\n', '\n')
before = text

# ---------------------------------------------------------------- counters
def section_counts(t):
    s0 = t.index('## Introduction'); s1 = t.index('## References')
    body = t[s0:s1]
    secs, cur, buf = [], None, []
    for ln in body.split('\n'):
        if ln.startswith('## '):
            if cur: secs.append((cur, '\n'.join(buf)))
            cur = ln[3:].strip(); buf = []
        else:
            buf.append(ln)
    if cur: secs.append((cur, '\n'.join(buf)))
    return secs

def counts(t):
    secs = section_counts(t)
    d = {}
    for name, cont in secs:
        c = '\n'.join(l for l in cont.split('\n') if not l.startswith('#'))
        d[name] = W(c)
    # strict: Introduction..Conclusion only (JCSM-countable prose)
    strict = 0
    started = False
    for name, cont in secs:
        if name == 'Introduction': started = True
        if not started: continue
        if name in ('Main tables', 'Figure legends'): continue
        if name in ('Supplementary material', 'Supplementary Methods',
                    'Supplementary tables (moved from the main text)',
                    'Supplementary figure legends (moved from the main text)',
                    'Data availability', 'Code availability', 'Ethics statement',
                    'Declarations'):
            break
        c = '\n'.join(l for l in cont.split('\n') if not l.startswith('#'))
        strict += W(c)
    # legacy: everything between Introduction and References minus tables/figs
    legacy = sum(v for k, v in d.items() if k not in ('Main tables', 'Figure legends'))
    # tracked metric: the one _apply_jcsm_compression.py used -> Introduction .. Main tables,
    # counting every word including ### subheadings
    a = t.index('## Introduction'); b = t.index('## Main tables')
    tracked = W(t[a:b])
    return d, strict, legacy, tracked

d0, strict0, legacy0, tracked0 = counts(text)
print("=== BEFORE ===")
for k, v in d0.items():
    print("  %-52s %6d" % (k[:52], v))
print("  TRACKED (Introduction..Main tables, incl. ### heads) =", tracked0)
print("  STRICT  (Introduction..Conclusion prose only)        =", strict0)
print("  LEGACY  (all sections minus tables/figs)             =", legacy0)

# ---------------------------------------------------------------- helpers
def sub1(old, new, tag):
    global text
    n = text.count(old)
    assert n == 1, "ABORT [%s]: %d occurrences (expect 1)" % (tag, n)
    text = text.replace(old, new)

def setline(idx, new, pre, suf, tag):
    global text
    lines = text.split('\n')
    old = lines[idx]
    assert old.startswith(pre), "ABORT [%s]: prefix mismatch -> %r" % (tag, old[:70])
    assert old.rstrip().endswith(suf), "ABORT [%s]: suffix mismatch -> %r" % (tag, old[-70:])
    lines[idx] = new
    text = '\n'.join(lines)
    print("  [L%d] %-22s %3d -> %3d words" % (idx + 1, tag, W(old), W(new)))

# ================================================================ B1: table numbers
print("\n--- B1 stale table cross-references ---")
sub1("reported in Tables 5 and 8,", "reported in Tables 3 and 5,", "S11 Tables 5 and 8")
sub1("Two consequences matter for reading Table 8.",
     "Two consequences matter for reading Table 5.", "L250a Table 8")
sub1("which require no approximation, in the note to Table 8.",
     "which require no approximation, in the note to Table 5.", "L250b Table 8 (2nd in same line)")
sub1("the confidence intervals and P values in Table 8 are implementation-independent",
     "the confidence intervals and P values in Table 5 are implementation-independent", "L256 Table 8")
sub1("narrower than the confidence intervals in Table 8",
     "narrower than the confidence intervals in Table 5", "L258 Table 8")
sub1("within the confidence intervals reported in Table 3.",
     "within the confidence intervals reported in Table 2.", "L348 Table 3 -> 2")
print("  5 sites patched")

# ================================================================ B2: logic / defensive
print("\n--- B2 logical and defensive sentences ---")
sub1("The significant modified Q statistic therefore reflects balanced scatter rather than directional pleiotropy.",
     "Modified Q was significantly elevated; the Egger intercepts provide no evidence that the residual scatter is directional.",
     "L63 balanced scatter")

setline(72,
    "Outcome-side triangulation in FinnGen supported the absence of a detectable TMFI association across three related pain endpoints (Table 6): multivariable estimates were 1.019 per SD (P=0.42) for `PAIN`, 0.983 (P=0.53) for limb pain and 1.001 (P=0.99) for fibromyalgia, with the univariable estimates agreeing. The BMI association was reproduced for limb pain (1.217 per SD, 1.116–1.327, P=9.0×10⁻⁶) and was directionally concordant but imprecise for `PAIN` (1.055, P=0.15) and fibromyalgia (1.037, P=0.79). The limb-pain estimate lies between the UK Biobank estimates for back pain (1.134) and knee pain (1.275). Conditional F statistics were 22.4 and 35.4 (UK Biobank: 22.1 and 35.2), Steiger held for 42 of 42 instruments, and leave-one-out left the limb-pain BMI coefficient between 0.176 and 0.207 (full-sample 0.196). Heterogeneity was present (modified Q 3.03, 1.46, 1.20); the Egger intercept was nominally significant for `PAIN` alone (P=0.041), an uncorrected nominal result that is not interpreted further. Minimum detectable effects were 1.07, 1.09 and 1.26 per SD, so the fibromyalgia endpoint contributes descriptive rather than inferential support. Because no participants are shared, no overlap bias was expected, so the concordance of the two designs indicates the bounded null is not an artefact of overlap.",
    "Repeating both stages against three FinnGen endpoints",
    "not an artefact of overlap.", "L73 FinnGen triangulation")

setline(76,
    "Two aspects give the null its interpretation. The internal comparator shows that the harmonization and outcome data could detect a body-composition effect of the size recovered for BMI in this contrast — not that the TMFI measurement is adequate — so a null obtained alongside it differs from a null in isolation. The comparator is not a positive control with a known answer: it establishes the model's sensitivity in this contrast, not that every excluded effect would have been recovered. The opposite-signed estimates mean the TMFI result is not a weakened echo of the adiposity signal. The conclusion is bounded: thigh muscle fat infiltration is unlikely to drive chronic pain intensity at anything like the scale of systemic adiposity, but effects below about 0.05 SD per SD remain open, and muscle tissue may still matter through other phenotypes or pathways.",
    "Two features make this more than a failed test",
    "other phenotypes or pathways.", "L77 principal finding")

sub1("conversely, the near-zero estimate sits against a comparator that was unambiguously detected, and that pairing licenses the interpretation.",
     "conversely, the near-zero estimate and its confidence interval sit against a comparator that was unambiguously detected.",
     "L81 licenses")

sub1("the fixed-effect fit reproduced both point estimates, so heterogeneity inflates variance rather than shifting estimates",
     "the fixed-effect fit returns the same point estimates by construction, so heterogeneity widens the intervals rather than shifting them",
     "L83 FE=MRE identity")

# ================================================================ B3: pre-specified / FinnGen elsewhere
print("\n--- B3 pre-specified and remaining FinnGen claims ---")
sub1("BMI was pre-specified as an internal comparator",
     "BMI was designated as an internal comparator", "abstract L11")
sub1("Cook-distance exclusions above the prespecified 4/n threshold",
     "Cook-distance exclusions above the 4/n threshold", "L49 prespecified")
sub1("### Pre-specified site-matched secondary analysis",
     "### Site-matched secondary analysis", "L52 heading")
sub1("### Table 5. Pre-specified site-matched secondary analysis: knee pain and back pain",
     "### Table 5. Site-matched secondary analysis: knee pain and back pain", "L167 table title")
sub1("the pre-specified site-matched secondary outcomes",
     "the site-matched secondary outcomes", "L205 fig legend")
sub1("No multiple-testing correction was applied, because the primary contrast is a single pre-specified comparison and the site-matched outcomes were pre-specified replications of it.",
     "No multiple-testing correction was applied, because the primary contrast is a single comparison fixed in advance; the site-matched and FinnGen analyses are secondary and exploratory.",
     "L57 multiplicity")
sub1("so eligible signals outside it may have been missed — a disclosed limitation, since every retained variant was present.",
     "so eligible signals outside it may have been missed — a disclosed limitation.", "L37 LD panel")

sub1("Both held in FinnGen: TMFI null for limb pain (OR 0.983 per SD, P=0.53), BMI associated with limb pain (1.217, 1.116–1.327, P=9×10⁻⁶).",
     "In FinnGen, TMFI was null for limb pain (OR 0.983 per SD, P=0.53) while BMI was associated with it (1.217, 1.116–1.327, P=9×10⁻⁶); estimates for the two further pain endpoints were directionally concordant but imprecise.",
     "abstract L19 FinnGen")

sub1("Against an explicit power bound, an internal comparator that behaved as expected, an independent replication and a quantified overlap, this is a bounded null: effects above about 0.05 SD per SD are unlikely; smaller ones remain open.",
     "Against an internal comparator, a quantified overlap and confidence intervals excluding effects above 0.031 SD per SD, this is a bounded null: effects of the order of the adiposity effect are excluded; smaller ones remain open.",
     "abstract L23 MDE threshold")

sub1("The same pattern held for site-matched knee and back pain, with a descriptively larger BMI effect at the weight-bearing joint, and for three FinnGen endpoints with no participant overlap.",
     "The TMFI null also held for site-matched knee and back pain, where BMI was descriptively larger at the weight-bearing joint, and across three related FinnGen endpoints, where the BMI association was reproduced for limb pain and imprecise for the other two.",
     "L89 conclusion FinnGen")

# ================================================================ B4: limitations, AI, L33, supp methods
print("\n--- B4 limitations, AI declaration, L33, supplementary flag ---")
setline(86,
    "Eight boundaries qualify the estimates, in three groups. **Source data.** The TMFI source GWAS carries an unresolved preprocessing discrepancy affecting the absolute exposure scale, so inference is standardized and raw coefficients should not be translated into physical units until the source authors confirm the flag, and the standardized bounds should be read as conditional on that confirmation. The LD reference panel (503 individuals, 1.84 million variants) contained only 22.75% of TMFI candidate variants, so eligible instruments outside it may have been missed. The GIANT release is genomic-control corrected (λ=1.61 inferred from its single-trait intercept), so intercepts involving BMI required rescaling and its released standard errors are inflated by about 27%, making the BMI conditional F conservative by a factor of about 1.6; sd<sub>BMI</sub> is 1.00 by definition. **Model validity.** Residual heterogeneity indicates that some variants influence CPI through pathways outside both exposures; the overlap audit identifies the product N_shared·ρ_phen rather than a participant head-count, so overlap figures should not be read as participant counts; and the site-matched odds ratios rest on the approximate BOLT-LMM conversion log OR = β/[μ(1−μ)], although every P value, F statistic and tipping point is invariant to it and the risk differences reported alongside require no approximation. **External validity.** The second knee-pain source returned a BMI association opposite to the epidemiological expectation (OR per SD 0.905, P=0.008), plausibly because that GWAS adjusted for adiposity-related covariates; it is used only to check the TMFI null. The replication is outcome-side only — the exposures remain UK Biobank instruments and the FinnGen endpoints are binary ICD-based diagnoses, so concordance is related-phenotype triangulation rather than same-phenotype replication, with the fibromyalgia endpoint (3,623 cases) excluding only effects larger than 1.26 per SD. Generalizability is limited to European ancestry, to lifetime-averaged genetically predicted differences rather than interventions, and to shifts across the full genetic range of the exposure rather than clinically attainable changes.",
    "Eight boundaries qualify the estimates",
    "clinically attainable changes.", "L87 limitations restructure")

sub1("BMI instead served as an internal comparator, demonstrating that the instrument set and outcome data can detect an effect of the size expected for a body-composition trait and permitting the null TMFI estimate to be read against the minimum detectable effect.",
     "BMI instead served as an internal comparator, showing that the outcome data and harmonization can recover an effect of the size observed for adiposity and permitting the null TMFI estimate to be read against the minimum detectable effect; it does not verify the TMFI measurement itself.",
     "L33 comparator narrowing")

sub1("Its consequence is confined to the vertical scale of the exposure coefficients: IVW-type estimators test whether the outcome associations are proportional to the exposure associations across instruments, and that proportionality is not created or destroyed by a monotone rescaling of the exposure betas. Because our conclusions are stated as standardized bounds and not as absolute physical effects, they do not depend on resolving the flag.",
     "To the extent that the two models return proportional exposure coefficients, which a liability-threshold model predicts only approximately, instrument selection and the z statistics and P values of the MR estimates are unchanged; the discrepancy then affects the vertical scale of the exposure coefficients and hence the standardized magnitudes. Because proportionality cannot be verified from the distributed file, the standardized bounds are reported as conditional on the exposure scale rather than as resolved physical effects.",
     "SuppMethods flag wording")

sub1("A large language model (WorkBuddy, an AI writing assistant) was used for three bounded tasks: checking grammar and wording, reorganising prose for concision, and cross-checking that every number quoted in the text, tables and figure legends matches the value in the underlying result files. It was not used to generate, clean or analyse any data, to select instruments or outcomes, to choose or implement statistical methods, or to draw scientific conclusions; no result, table or figure value originates from the model. Every numerical statement in this manuscript is regenerated deterministically by `manuscript_numbers.py` from the frozen Stage 1, Stage 2, Stage 2B, Stage 3 and Stage 5 CSV outputs, and the two independent implementations of the estimators (`stage5_mr.py` under Python and `stage5_validation_baser.R` under base R) were written and reviewed by the authors. The authors have verified all content and take full responsibility for it.",
     "A large language model (WorkBuddy, an AI writing assistant) was used for three bounded tasks: checking grammar and wording, reorganising prose for concision, and cross-checking that every number quoted in the text, tables and figure legends matches the underlying result files. It was not used to generate, clean or analyse data, to select instruments or outcomes, or to choose or implement statistical methods; no reported value originates from the model. Every number is regenerated deterministically by `manuscript_numbers.py` from the frozen stage outputs. The authors verified all content and take full responsibility for it.",
     "AI declaration")

# ================================================================ assertions
print("\n--- assertions ---")
MUST = [
    # key numbers (subset, full list re-checked in final audit)
    "0.008", "0.558", "0.005", "0.721", "0.090", "0.038–0.143", "0.00075",
    "−0.002", "0.898", "22.12", "35.23", "0.048", "0.052", "1.2×10⁻¹¹",
    "1.0146", "1.0005", "1.2751", "1.1344", "0.983", "1.217", "1.116–1.327",
    "9.0×10⁻⁶", "46,347", "134,627", "106,697", "499,362", "22.75%",
    "1.61", "0.905", "1.26", "3,623", "42/42", "0.031", "0.037",
    # reviewer-mandated wording
    "descriptively larger", "internal comparator", "structural ceiling of 0.890 within ~2 SE",
    "bounded null", "Supplementary Table S13", "Eight boundaries",
    # new wording present
    "Outcome-side triangulation", "directionally concordant but imprecise",
    "by construction", "conditional on the exposure scale", "fixed in advance",
    "secondary and exploratory", "no overlap bias was expected",
    "in three groups", "**Source data.**", "**Model validity.**", "**External validity.**",
]
for m in MUST:
    assert m in text, "ABORT: lost %r" % m
print("  %d must-keep strings present" % len(MUST))

GONE = [
    "Two features make this more than a failed test",
    "licenses the interpretation",
    "balanced scatter rather than directional pleiotropy",
    "within expectation for three endpoints under the no-correction rule",
    "reproduced both central findings",
    "adiposity gradient across weight-bearing sites",
    "cross-trait intercept is zero by construction",
    "inference does not depend on re-clumping",
    "pre-specified as an internal comparator",
    "pre-specified replications of it",
    "reproduced both point estimates, so heterogeneity inflates variance",
    "they do not depend on resolving the flag",
    "Table 8",
]
for g in GONE:
    if g in text:
        for m in re.finditer(re.escape(g), text):
            s = text[max(0, m.start() - 70):m.start() + 40].replace("\n", " ")
            print("   DEBUG leftover:", s)
    assert g not in text, "ABORT: still present %r" % g
print("  %d must-go strings absent" % len(GONE))
for bad in ("Table 7", "Table 9", "Figure 3", "Figure 4"):
    assert bad not in text, "ABORT: stale float ref %r" % bad
print("  no stale Table 7/9 or Figure 3/4 references")

# structural
assert text.count("### Table ") == 6, "main tables != 6"
assert len(re.findall(r"^### Figure \d", text, re.M)) == 2, "main figures != 2"
assert len(re.search(r"\*\*Keywords:\*\*", text).group(0)) > 0
kw = re.search(r"\*\*Keywords:\*\*\s*(.+)", text).group(1)
assert 4 <= len([k for k in kw.split(';') if k.strip()]) <= 6, "keywords out of 4-6"
print("  6 tables + 2 figures, keywords %d" % len([k for k in kw.split(';') if k.strip()]))

ab = text.split('## Abstract', 1)[1].split('\n## Introduction', 1)[0]
ab_body = '\n'.join(l for l in ab.split('\n') if not l.startswith('#') and '**Keywords:**' not in l)
print("  abstract words =", W(ab_body))
assert W(ab_body) <= 400, "ABORT: abstract > 400"

d1, strict1, legacy1, tracked1 = counts(text)
print("\n=== AFTER ===")
for k, v in d1.items():
    print("  %-52s %6d" % (k[:52], v))
print("  TRACKED =", tracked0, "->", tracked1)
print("  STRICT  =", strict0, "->", strict1)
print("  LEGACY  =", legacy0, "->", legacy1)

if tracked1 > 4500:
    sys.exit("ABORT: tracked body over 4500 (nothing written)")

# write with CRLF
out = text.replace('\n', '\r\n')
open(P, 'wb').write(out.encode('utf-8'))
print("\nWROTE", len(out.encode('utf-8')), "bytes")
