# -*- coding: utf-8 -*-
"""P1 revision pass: section 3.2 / 5.2 x4 / 5.3 / 5.4 of the targeted re-review,
plus the word-count reduction package. Every substitution asserts occurrence == 1."""
import re, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

P = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"
BAK = P + ".bak_p1_2026-09-20"

raw = open(P, 'rb').read()
text = raw.decode('utf-8')
assert text.count('\r\n') > 0 and len(re.findall(r'(?<!\r)\n', text)) == 0, "expected pure CRLF"
W = lambda s: len(re.findall(r"[^\s]+", s))

i0 = text.index("## Introduction")
i1 = text.index("## Main tables")
def tracked(t):
    return W(t[t.index("## Introduction"):t.index("## Main tables")])
def strict(t):
    seg = t[t.index("## Introduction"):t.index("## Main tables")]
    return W('\n'.join(l for l in seg.split('\n') if not l.startswith('#')))

print("BEFORE  tracked=%d  strict=%d" % (tracked(text), strict(text)))

def sub1(old, new, tag, expect=1):
    n = text.count(old)
    assert n == expect, "ABORT %s: found %d occurrences, expected %d\n  %r" % (tag, n, expect, old[:120])
    globals()['text'] = text.replace(old, new)
    return tag

# ---------------------------------------------------------------- 5.2  FinnGen downgrade (4 sites)
sub1("and replicated both stages in FinnGen, which shares no participants with UK Biobank.",
     "and repeated both contrasts against FinnGen endpoints, which share no participants with UK Biobank.",
     "5.2-L42 introduction")

sub1("### External replication in a cohort with no participant overlap",
     "### Outcome-side triangulation in a cohort with no participant overlap",
     "5.2-L66 methods heading")

sub1("### External replication: outcome-side triangulation with no participant overlap",
     "### Outcome-side triangulation in FinnGen: no participant overlap",
     "5.2-L84 results heading")

sub1("so the concordance of the two designs indicates the bounded null is not an artefact of overlap.",
     "so the concordance of the two designs is consistent with the bounded null not being an artefact of overlap.",
     "5.2-L85 results tail")

sub1("The TMFI null was reproduced across three related FinnGen endpoints that share no participants "
     "with UK Biobank (Table 6), and the BMI association was reproduced for limb pain, so the bounded "
     "null is not an artefact of the exposure–outcome overlap that the tipping-point arguments address.",
     "No detectable TMFI association was observed across three related FinnGen endpoints that share no "
     "participants with UK Biobank (Table 6), and the BMI association was reproduced for limb pain, so "
     "the concordance is consistent with the bounded null not being an artefact of the exposure–outcome "
     "overlap that the tipping-point arguments address.",
     "5.2-L88 discussion")

# consistency outside the tracked window
sub1("### Table 6. External replication in FinnGen R12 (no participant overlap with UK Biobank)",
     "### Table 6. Outcome-side triangulation in FinnGen R12 (no participant overlap with UK Biobank)",
     "Table 6 heading")
sub1("| External replication outcome |", "| External triangulation outcome |",
     "Table 1 role cells", expect=3)

# ---------------------------------------------------------------- 5.3  Conclusion: separate CI from MDE
sub1("Because the same model recovered the BMI effect at approximately 92% power while retaining 80% "
     "power for a TMFI effect of 0.048 SD per SD, this is a bounded null: thigh muscle fat infiltration "
     "is unlikely to influence chronic pain intensity at anywhere near the scale of systemic adiposity, "
     "although effects below about 0.05 SD per SD remain open.",
     "Because the same model recovered the BMI effect at approximately 92% power, this is a bounded "
     "null: the upper confidence limits exclude any positive TMFI effect larger than 0.031 SD per SD "
     "(0.029 for CPI_no_imaging), and the 80%-power minimum detectable effects were 0.048 and 0.053 SD "
     "per SD. Thigh muscle fat infiltration is unlikely to influence chronic pain intensity at anywhere "
     "near the scale of systemic adiposity, although effects below those confidence limits remain "
     "compatible.",
     "5.3-L101 conclusion")

# ---------------------------------------------------------------- 5.4  fixed-effect wording
sub1("the fixed-effect fit returns the same point estimates by construction, so heterogeneity widens "
     "the intervals rather than shifting them;",
     "the multiplicative random-effects specification retained the fixed-effect point estimate by "
     "construction while scaling its standard error for overdispersion;",
     "5.4-L95 robustness boundary")

# ---------------------------------------------------------------- 3.2  TMFI scale: keep the anchor, narrow the consequence
sub1("A public script nevertheless carries a `--fastGWA-mlm-binary` flag the text does not describe; "
     "if it reflected the released statistics the bounds would tighten, not loosen, and both contrasts "
     "are unaffected.",
     "A public script nevertheless carries a `--fastGWA-mlm-binary` flag the text does not describe. "
     "Because every z statistic, P value and conditional F statistic is invariant to a rescaling of the "
     "exposure, neither contrast's direction or significance depends on which specification the released "
     "statistics follow; the standardized magnitudes above are conditional on the published "
     "continuous-trait specification.",
     "3.2-L99 limitations")

# ---------------------------------------------------------------- word-count reduction package
sub1("For BMI, three independent checks — published per-allele effects at TMEM18, MC4R and SEC16B, "
     "variance explained across 71 independent loci, and the estimator's internal consistency for the "
     "other three sources — establish that the GIANT coefficients are already per SD of BMI, the "
     "implied-SD surplus reflecting genomic-control-corrected standard errors (single-trait intercept "
     "0.62), so sd<sub>BMI</sub>=1.00 (derivation in Supplementary Methods).",
     "Three independent checks establish that the GIANT coefficients are already per SD of BMI, the "
     "implied-SD surplus reflecting genomic-control-corrected standard errors (single-trait intercept "
     "0.62), so sd<sub>BMI</sub>=1.00 (Supplementary Methods).",
     "trim-L57 BMI scale derivation")

# power sentence leaves the body but goes into the Table 4 note (outside the tracked window)
sub1(" Power reached 83% and 75% for a 0.05 SD TMFI effect and exceeded 99% for 0.075 SD.", "",
     "trim-L77 power sentence")
sub1("The observed BMI direct effects correspond to approximately 92% power on CPI_full and 93% on "
     "CPI_no_imaging.",
     "The observed BMI direct effects correspond to approximately 92% power on CPI_full and 93% on "
     "CPI_no_imaging. For a 0.05 SD TMFI effect, power was 83% (CPI_full) and 75% (CPI_no_imaging), "
     "and exceeded 99% for 0.075 SD.",
     "trim-L77 -> Table 4 note")

# risk differences are already in the Table 5 note
sub1(" on the risk-difference scale, TMFI gave +0.24 and +0.01 percentage points against +4.08 and "
     "+2.40 for BMI.", "",
     "trim-L83 risk differences")

sub1(" Overlap between the Cook-threshold and review-target sets is reported explicitly and not counted "
     "as independent replication.", "",
     "trim-L61 duplicate of Results L79")

sub1("Systemic adiposity may influence pain through mechanical loading, inflammation, metabolic "
     "dysfunction and altered pain processing, and",
     "Systemic adiposity may influence pain through mechanical, inflammatory and metabolic pathways, and",
     "trim-L40 mechanism list")

sub1("MRC IEU OpenGWAS ukb-b-16254 and ukb-b-9838 [12], fitted with BOLT-LMM, a linear mixed model, "
     "to the binary phenotypes)",
     "MRC IEU OpenGWAS ukb-b-16254 and ukb-b-9838 [12], fitted with BOLT-LMM, a linear mixed model)",
     "trim-L65 BOLT-LMM gloss")

sub1(" Phenotype identities were verified independently of repository metadata (leading signals GDF5 "
     "for knee and SOX5-adjacent for back pain; allele-frequency agreement with the panel r≈0.998; "
     "full checks in Supplementary Methods).",
     " Phenotype identities were verified independently of repository metadata (leading signals GDF5 "
     "and SOX5-adjacent; allele-frequency agreement r≈0.998; Supplementary Methods).",
     "trim-L65 phenotype checks")

sub1("the identity SE² × 2Nf(1−f) = 1/[μ(1−μ)] holds to within 5% on non-significant instruments, "
     "the linear alternative is excluded by factors of 17 to 10,300, and no rescaling was applied.",
     "the identity SE² × 2Nf(1−f) = 1/[μ(1−μ)] holds to within 5% and excludes the linear alternative "
     "by factors of 17 to 10,300, and no rescaling was applied.",
     "trim-L67 scale verification")

sub1(", but effects below about 0.05 SD per SD remain open, and muscle tissue may still matter",
     "; muscle tissue may still matter",
     "trim-L89 duplicate of Conclusion")

# ---------------------------------------------------------------- content integrity
MUST = ["0.008", "0.558", "0.005", "0.721", "0.090", "0.038–0.143", "0.00075",
        "−0.002", "0.898", "22.12", "35.23", "0.048", "0.053", "0.052", "1.2×10⁻¹¹",
        "1.0146", "1.0005", "1.2751", "1.1344", "0.983", "1.217", "1.116–1.327",
        "9.0×10⁻⁶", "46,347", "134,627", "106,697", "499,362", "22.75%", "1.61",
        "0.905", "1.26", "3,623", "42/42", "0.031", "0.029", "0.037", "0.910", "0.890",
        "0.0014", "1.0681", "1.0621", "1.019", "1.001", "0.176", "0.196",
        "descriptively larger", "internal comparator", "bounded null",
        "structural ceiling of 0.890 within ~2 SE", "Supplementary Table S13",
        "Eight boundaries", "in three groups", "Outcome-side triangulation",
        "directionally concordant but imprecise", "conditional on the exposure scale",
        "fixed in advance", "secondary and exploratory", "no overlap bias was expected",
        "Mendelian randomization estimated", "MRlap", "STROBE-MR",
        # new in P1
        "no evidence of average directional pleiotropy",
        "is consistent with the bounded null not being an artefact of overlap",
        "No detectable TMFI association was observed",
        "invariant to a rescaling of the exposure",
        "conditional on the published continuous-trait specification",
        "scaling its standard error for overdispersion",
        "83%", "0.075 SD", "+0.24", "+4.08", "GDF5", "r≈0.998", "0.62"]
bad = [m for m in MUST if m not in text]
assert not bad, "ABORT must-keep missing: %s" % bad
print("  must-keep : %d checked, 0 missing" % len(MUST))

GONE = ["replicated both stages in FinnGen",
        "External replication in a cohort",
        "### External replication:",
        "indicates the bounded null is not an artefact",
        "The TMFI null was reproduced",
        "bounds would tighten, not loosen",
        "both contrasts are unaffected",
        "heterogeneity widens the intervals rather than shifting them",
        "effects below about 0.05 SD per SD remain open",
        "Power reached 83% and 75% for a 0.05 SD TMFI effect",
        "on the risk-difference scale, TMFI gave",
        "Overlap between the Cook-threshold and review-target sets is reported",
        "External replication outcome",
        "### Table 6. External replication"]
left = [g for g in GONE if g in text]
assert not left, "ABORT must-go remaining: %s" % left
print("  must-go   : %d checked, 0 remaining" % len(GONE))

# ---------------------------------------------------------------- limits
t_after = tracked(text)
s_after = strict(text)
print("\nAFTER   tracked=%d  strict=%d" % (t_after, s_after))
print("  delta tracked = %+d" % (t_after - tracked(raw.decode('utf-8'))))

ab = text.split('## Abstract', 1)[1].split('\n## Introduction', 1)[0]
ab_body = '\n'.join(l for l in ab.split('\n') if not l.startswith('#') and '**Keywords:**' not in l)
print("  abstract = %d / 400" % W(ab_body))

if t_after > 4400:
    sys.exit("ABORT: tracked body still %d, target <= 4400 (nothing written)" % t_after)

assert "## Introduction" in text and "## Main tables" in text
assert text.count('\r\n') > 0

with io.open(BAK, 'w', encoding='utf-8', newline='') as f:
    f.write(raw.decode('utf-8'))
with io.open(P, 'w', encoding='utf-8', newline='') as f:
    f.write(text)
print("\nWROTE %s (backup %s)" % (P, BAK))
