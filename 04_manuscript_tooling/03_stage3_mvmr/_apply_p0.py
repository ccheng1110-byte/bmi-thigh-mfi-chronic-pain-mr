# -*- coding: utf-8 -*-
"""
P0 batch of the JCSM pre-submission fixes (2026-09-20).

1. Abstract: "with no directional pleiotropy" -> "with no evidence of average
   directional pleiotropy"  (an non-significant Egger intercept is absence of
   evidence, not evidence of absence; also makes the abstract agree with the
   Discussion, which already uses the hedged form).
2. Declarations: delete the internal working note "Items marked ... are
   journal-administrative and do not affect any analysis reported here."
3. Title page counts: refresh stale "22 references" / "4,417 words".
4. Ethics statement: add the JCSM-mandated certification sentence, with the
   required citation of the ethical-guidelines article -> new reference [24].

All four changes are outside the TRACKED word-count window
(Introduction .. Main tables), so the body count is unchanged.
"""
import io
import os
import re
import sys

P = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md")

with io.open(P, "r", encoding="utf-8", newline="") as f:
    text = f.read()

orig = text


def W(s):
    return len(re.findall(r"[^\s]+", s))


def sub1(old, new, tag):
    global text
    n = text.count(old)
    if n != 1:
        sys.exit("ABORT [%s]: expected 1 occurrence, found %d\n  %r"
                 % (tag, n, old[:130]))
    text = text.replace(old, new)
    print("  ok  %s" % tag)


# ---------------------------------------------------------------- 1. abstract
sub1("with no directional pleiotropy.",
     "with no evidence of average directional pleiotropy.",
     "Abstract: pleiotropy hedged")

# ------------------------------------------------- 2. declarations: drop note
sub1("All items below are stated explicitly. Items marked *to be completed by the "
     "authors* are journal-administrative and do not affect any analysis reported here.",
     "All items below are stated explicitly.",
     "Declarations: internal working note removed")

# ------------------------------------------------- 3. ethics: JCSM wording
sub1("and no individual-level Finnish data were accessed.",
     "and no individual-level Finnish data were accessed. All authors of this manuscript "
     "comply with the guidelines of ethical authorship and publishing in the Journal of "
     "Cachexia, Sarcopenia and Muscle [24]. The manuscript does not contain clinical "
     "studies or patient data.",
     "Ethics: JCSM certification + no-patient-data statement")

# ------------------------------------------------- 4. reference [24]
OLD_R = ("23. UK Biobank. Data-field 22435: muscle fat infiltration (anterior thigh). UK "
         "Biobank Showcase. Available from: "
         "https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=22435. "
         "Accessed 20 September 2026.")

NEW_R = (OLD_R +
         "\n24. von Haehling S, Coats AJS, Anker SD. Ethical guidelines for publishing in "
         "the Journal of Cachexia, Sarcopenia and Muscle: update 2021. "
         "*J Cachexia Sarcopenia Muscle*. 2021;12(6):2259-2261. "
         "doi:[10.1002/jcsm.12899](https://doi.org/10.1002/jcsm.12899).")

sub1(OLD_R, NEW_R, "Reference [24] JCSM ethical guidelines 2021")

# ------------------------------------------------- 5. title-page counts
def ab_words(t):
    a = t.split("## Abstract", 1)[1].split("\n## Introduction", 1)[0]
    a = re.sub(r"^###\s*.*$", "", a, flags=re.M)
    a = re.sub(r"^\*\*Keywords:\*\*.*$", "", a, flags=re.M)
    return W(a)


tracked = W(text[text.index("## Introduction"):text.index("## Main tables")])
abstract = ab_words(text)

sub1("**Counts.** Abstract 338 words (limit 400); main text 4,417 words from Introduction "
     "to Conclusion (limit 4,500); 6 main tables and 2 main figures (limit 8); "
     "22 references (limit 40); supplementary material unlimited.",
     "**Counts.** Abstract %d words (limit 400); main text %d words from Introduction to the "
     "end of the Conclusion (limit 4,500); 6 main tables and 2 main figures (limit 8); "
     "24 references (limit 40); supplementary material unlimited." % (abstract, tracked),
     "Title page: counts refreshed")

# ---------------------------------------------------------------- assertions
print()
print("  TRACKED body words = %d / 4500  (margin %d)" % (tracked, 4500 - tracked))
print("  abstract words     = %d / 400" % abstract)

if tracked > 4500:
    sys.exit("ABORT: body over 4500 words (nothing written)")
if abstract > 400:
    sys.exit("ABORT: abstract over 400 words (nothing written)")

KEEP = [
    "0.048", "0.031", "0.090", "0.038", "0.143", "0.00075",
    "22.12", "35.23", "2.03", "1.217", "1.116", "1.327", "0.983",
    "46,347", "134,627", "106,697", "499,362", "237,944", "42", "93", "80%",
    "descriptively larger", "internal comparator",
    "structural ceiling of 0.890 within ~2 SE", "Supplementary Table S13",
    "Eight boundaries", "bounded null", "Mendelian randomization",
    "raw-scale coefficients reported alongside",
    # new in this patch
    "no evidence of average directional pleiotropy",
    "guidelines of ethical authorship and publishing in the Journal of Cachexia, "
    "Sarcopenia and Muscle [24]",
    "does not contain clinical studies or patient data",
    "10.1002/jcsm.12899",
    "percentage points",
]
missing = [k for k in KEEP if k not in text]
print("  must-keep: %d/%d present" % (len(KEEP) - len(missing), len(KEEP)))
if missing:
    sys.exit("ABORT: lost %r" % missing)

GONE = [
    "with no directional pleiotropy",
    "journal-administrative and do not affect any analysis",
    "22 references (limit 40)",
    "main text 4,417 words",
]
left = [g for g in GONE if g in text]
print("  must-go  : %d/%d absent" % (len(GONE) - len(left), len(GONE)))
if left:
    sys.exit("ABORT: still present %r" % left)

refs = [int(x) for x in re.findall(r"(?m)^(\d+)\.\s", text[text.index("## References"):])]
print("  references         = %d (contiguous: %s)"
      % (len(refs), refs == list(range(1, len(refs) + 1))))
if refs != list(range(1, 25)):
    sys.exit("ABORT: reference numbering wrong (nothing written)")

title = text.split("\n")[0].lstrip("# ").strip()
kw = re.search(r"\*\*Keywords:\*\*(.+)", text).group(1)
print("  title words        = %d / 17" % W(title))
print("  keywords           = %d / 4-6" % len([k for k in kw.split(";") if k.strip()]))

tabs = re.findall(r"(?m)^### Table \d+\.",
                  text[text.index("## Main tables"):text.index("## Figure legends")])
figs = re.findall(r"(?m)^### Figure \d+\.",
                  text[text.index("## Figure legends"):text.index("## Supplementary material")])
print("  main floats        = %d (tables %d + figures %d) / 8"
      % (len(tabs) + len(figs), len(tabs), len(figs)))
if len(tabs) + len(figs) != 8:
    sys.exit("ABORT: main float count changed (nothing written)")

placeholders = text.count("To be completed by the authors")
print("  author placeholders= %d (title page 4 + declarations 5)" % placeholders)
if placeholders != 9:
    sys.exit("ABORT: placeholder count changed (nothing written)")

if text == orig:
    sys.exit("ABORT: no change")

# ---------------------------------------------------------------- write (CRLF)
data = text.replace("\r\n", "\n").replace("\n", "\r\n").encode("utf-8")
with io.open(P, "wb") as f:
    f.write(data)
print("\n  WRITTEN %d bytes" % len(data))
