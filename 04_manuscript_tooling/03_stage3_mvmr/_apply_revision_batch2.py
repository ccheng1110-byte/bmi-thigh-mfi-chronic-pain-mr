# -*- coding: utf-8 -*-
"""JCSM revision pass 2: remaining FinnGen over-claims + word-margin compression."""
import re, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

P = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"
W = lambda s: len(re.findall(r"[^\s]+", s))
text = open(P, 'rb').read().decode('utf-8').replace('\r\n', '\n')

def tracked(t):
    a = t.index('## Introduction'); b = t.index('## Main tables')
    return W(t[a:b])

def strict(t):
    s0 = t.index('## Introduction')
    seg = t[s0:t.index('## Main tables')]
    return W('\n'.join(l for l in seg.split('\n') if not l.startswith('#')))

def sub1(old, new, tag):
    global text
    n = text.count(old)
    assert n == 1, "ABORT [%s]: %d occurrences" % (tag, n)
    text = text.replace(old, new)
    print("  %-38s %+d words" % (tag, W(new) - W(old)))

t0, s0 = tracked(text), strict(text)
print("BEFORE  tracked=%d  strict=%d" % (t0, s0))

# ---------------------------------------------------------------- remaining FinnGen over-claims
print("\n--- remaining FinnGen over-claims ---")
sub1("### External replication: both contrasts hold with no participant overlap",
     "### External replication: outcome-side triangulation with no participant overlap",
     "L72 heading")
sub1("Both contrasts were recovered in FinnGen, which shares no participants with UK Biobank (Table 6), so the bounded null is not an artefact of the exposure–outcome overlap that the tipping-point arguments address.",
     "The TMFI null was reproduced across three related FinnGen endpoints that share no participants with UK Biobank (Table 6), and the BMI association was reproduced for limb pain, so the bounded null is not an artefact of the exposure–outcome overlap that the tipping-point arguments address.",
     "L76 principal finding")

# ---------------------------------------------------------------- word-margin trims
print("\n--- word-margin trims (no numbers or reviewer wording touched) ---")
sub1("overlap between the UK Biobank-derived exposure and outcome samples biases estimates toward the confounded observational association; we quantified it rather than asserting it.",
     "overlap between the UK Biobank-derived exposure and outcome samples biases estimates toward the confounded observational association; we quantified it.",
     "L51 defensive clause")
sub1("so a significant intercept need not imply consequential bias; the identified quantity is the product N_shared·ρ_phen, not a participant head-count.",
     "so a significant intercept need not imply consequential bias.",
     "L51 duplicate of limitation (iv)")
sub1("R 4.6.1 with MVMR 0.4.8 (Stage 2);", "R 4.6.1 with MVMR 0.4.8;", "L57 redundant tag")
sub1("Stage 2 was therefore needed to determine whether that null conceals an effect visible only once BMI is held constant, and to bound it.",
     "Stage 2 was therefore needed to determine whether that null conceals an effect visible only once BMI is held constant.",
     "L60 redundant tail")
sub1(" For CPI_full the Cook-threshold set and the review-target set were the same five variants, so those two analyses are the same fitted model and are not independent replications; the convergence of a statistical influence criterion and an annotation-based criterion on the same variants is itself informative.",
     " For CPI_full the Cook-threshold set and the review-target set were the same five variants, so those two analyses are the same fitted model and are not independent replications.",
     "L67 self-commentary")
sub1("Influence diagnostics are described under Robustness analyses.", "", "L43 cross-ref")
sub1(" The adjacent positive findings for other muscle phenotypes [6,7] suggest the next test: jointly conditioning several correlated muscle traits for an adiposity-independent effect.",
     "",
     "L85 duplicated by Conclusion")

# ---------------------------------------------------------------- assertions
print("\n--- assertions ---")
MUST = ["0.008", "0.558", "0.005", "0.721", "0.090", "0.038–0.143", "0.00075",
        "−0.002", "0.898", "22.12", "35.23", "0.048", "0.052", "1.2×10⁻¹¹",
        "1.0146", "1.0005", "1.2751", "1.1344", "0.983", "1.217", "1.116–1.327",
        "9.0×10⁻⁶", "46,347", "134,627", "106,697", "499,362", "22.75%",
        "1.61", "0.905", "1.26", "3,623", "42/42", "0.031", "0.037", "0.910",
        "0.890", "0.0014", "1.0681", "1.0621",
        "descriptively larger", "internal comparator",
        "structural ceiling of 0.890 within ~2 SE", "bounded null",
        "Supplementary Table S13", "Eight boundaries", "in three groups",
        "Outcome-side triangulation", "directionally concordant but imprecise",
        "conditional on the exposure scale", "fixed in advance",
        "secondary and exploratory", "no overlap bias was expected",
        "not independent replications", "MRlap"]
for m in MUST:
    assert m in text, "ABORT: lost %r" % m
print("  %d must-keep strings present" % len(MUST))

GONE = ["Both contrasts were recovered in FinnGen", "both contrasts hold with no participant overlap",
        "Two features make this more than a failed test", "licenses the interpretation",
        "balanced scatter rather than directional pleiotropy",
        "within expectation for three endpoints under the no-correction rule",
        "reproduced both central findings", "adiposity gradient across weight-bearing sites",
        "cross-trait intercept is zero by construction", "inference does not depend on re-clumping",
        "pre-specified as an internal comparator", "pre-specified replications of it",
        "Table 8", "not a participant head-count. The GIANT"]
for g in GONE:
    assert g not in text, "ABORT: still present %r" % g
print("  %d must-go strings absent" % len(GONE))
assert "[6]" in text and "[7]" in text, "refs 6,7 must still be cited"
assert "[22]" in text, "ref 22 must still be cited"

t1, s1 = tracked(text), strict(text)
print("\nAFTER   tracked=%d (Δ%+d)  strict=%d (Δ%+d)" % (t1, t1 - t0, s1, s1 - s0))
assert t1 <= 4500, "ABORT: tracked over 4500"

ab = text.split('## Abstract', 1)[1].split('\n## Introduction', 1)[0]
ab_body = '\n'.join(l for l in ab.split('\n') if not l.startswith('#') and '**Keywords:**' not in l)
print("abstract words =", W(ab_body))
assert W(ab_body) <= 400
assert text.count("### Table ") == 6 and len(re.findall(r"^### Figure \d", text, re.M)) == 2

open(P, 'wb').write(text.replace('\n', '\r\n').encode('utf-8'))
print("WROTE", len(text.replace('\n', '\r\n').encode('utf-8')), "bytes")
