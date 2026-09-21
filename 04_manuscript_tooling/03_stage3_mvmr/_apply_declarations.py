"""Add the mandatory LLM-use statement and expand Declarations into the eight BMC
Medicine items. Uniqueness assertions on every anchor."""
from pathlib import Path

P = Path(r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md")
t = P.read_text(encoding="utf-8")

LLM_SUBSEC = """### Use of large language models

A large language model (WorkBuddy, an AI writing assistant) was used for three bounded tasks: checking grammar and wording, reorganising prose for concision, and cross-checking that every number quoted in the text, tables and figure legends matches the value in the underlying result files. It was not used to generate, clean or analyse any data, to select instruments or outcomes, to choose or implement statistical methods, or to draw scientific conclusions; no result, table or figure value originates from the model. Every numerical statement in this manuscript is regenerated deterministically by `manuscript_numbers.py` from the frozen Stage 1, Stage 2, Stage 2B, Stage 3 and Stage 5 CSV outputs, and the two independent implementations of the estimators (`stage5_mr.py` under Python and `stage5_validation_baser.R` under base R) were written and reviewed by the authors. The authors have verified all content and take full responsibility for it.

"""

OLD_DECL = ("Funding, competing interests, author contributions and correspondence details are to "
            "be completed by the authors before submission; none of these affects the analyses "
            "reported here.")

NEW_DECL = """All items below are stated explicitly. Items marked *to be completed by the authors* are journal-administrative and do not affect any analysis reported here.

**Ethics approval and consent to participate.** Not applicable as a separate approval: the study uses publicly available summary-level GWAS data and no individual-level data were accessed. The underlying UK Biobank studies were conducted under the relevant UK Biobank ethical and governance approvals as reported by the original data-generating studies [14,15], each of which reported participant consent and ethical review.

**Consent for publication.** Not applicable. No individual-level data, images or other potentially identifiable material are reported.

**Availability of data and materials.** See the Data availability and Code availability sections. All source GWAS summary statistics are publicly available from the repositories cited there; all frozen intermediate CSV outputs, analysis scripts and figure code will be deposited in a public version-controlled repository with a permanent DOI at the time of submission.

**Competing interests.** *To be completed by the authors.* [Declare any financial or non-financial competing interests, or state "The authors declare that they have no competing interests."]

**Funding.** *To be completed by the authors.* [Name all funding sources and grant numbers, or state "This research received no specific funding."]

**Authors' contributions.** *To be completed by the authors.* [Use CRediT-style statements, e.g. conceptualization, data curation, formal analysis, methodology, software, validation, writing – original draft, writing – review & editing.]

**Acknowledgements.** *To be completed by the authors.* [Name anyone who contributed but does not meet authorship criteria, and any technical assistance.]

**Authors' information.** *To be completed by the authors.* [Optional: affiliations, or a note on the corresponding author.]

**Use of artificial intelligence tools.** Declared in the Methods, under "Use of large language models"."""

REPL = [
    ("Raw-scale coefficients are reported alongside all standardized values.\n\n## Results\n",
     "Raw-scale coefficients are reported alongside all standardized values.\n\n"
     + LLM_SUBSEC + "## Results\n", 1),
    (OLD_DECL, NEW_DECL, 1),
]

for old, new, exp in REPL:
    c = t.count(old)
    if c != exp:
        raise SystemExit("ABORT (%d != %d) for:\n  %s" % (c, exp, old[:100]))
    t = t.replace(old, new, 1)

P.write_text(t, encoding="utf-8")
print("ok")
