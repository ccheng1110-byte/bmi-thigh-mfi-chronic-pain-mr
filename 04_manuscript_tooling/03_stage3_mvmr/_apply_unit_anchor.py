# -*- coding: utf-8 -*-
"""
Add the TMFI unit anchor (percentage points of fat fraction) to the manuscript.

Evidence gathered 2026-09-20 from the Aging Cell source paper (PMID 42675546,
PMC13530289, free full text) and UK Biobank data-field 22435:
  - Methods 4.6: continuous linear mixed model, GCTA v1.94.1 fastGWA-mlm,
    no standardisation / INT / rank transform / dichotomisation.
  - Methods 4.2: total TMFI = volume-weighted mean of four thigh-region fat
    fractions from AMRA neck-to-knee Dixon MRI.
  - UK Biobank field 22435 (anterior thigh MFI): units percent, mean 7.3361,
    SD 1.87868 -> implied sd_TMFI = 2.03 is ~2.0 percentage points.
  - Effect magnitudes: top instrument beta = 0.294 -> 0.145 SD per allele
    (plausible) but OR 1.34 per allele if read as log-OR (implausible).
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

# ---------------------------------------------------------------- helpers
def sub1(old, new, tag):
    global text
    n = text.count(old)
    if n != 1:
        sys.exit("ABORT [%s]: expected 1 occurrence, found %d\n  %r" % (tag, n, old[:120]))
    text = text.replace(old, new)
    print("  ok  %s" % tag)


# ---------------------------------------------------------------- 1. Methods
OLD_M = ("Because the TMFI source reports no verifiable physical unit, all primary "
         "inference is standardized (SD of outcome per SD of exposure), with raw-scale "
         "coefficients reported alongside.")

NEW_M = ("TMFI is a volume-weighted mean of fat fraction across four thigh regions and is "
         "expressed in percentage points, so its implied scale of 2.03 is about 2.0 percentage "
         "points; all primary inference is nevertheless standardized (SD of outcome per SD of "
         "exposure), with raw-scale coefficients reported alongside.")

sub1(OLD_M, NEW_M, "Methods: TMFI unit anchor")

# ---------------------------------------------------------------- 2. Limitations
OLD_L = ("The TMFI source GWAS carries an unresolved preprocessing discrepancy affecting the "
         "absolute exposure scale, so inference is standardized and raw coefficients should not "
         "be translated into physical units until the source authors confirm the flag, and the "
         "standardized bounds should be read as conditional on that confirmation.")

NEW_L = ("The TMFI source is a continuous linear mixed model reported in percentage points of "
         "fat fraction [1]; the implied scale of 2.03 is consistent with the UK Biobank "
         "anterior-thigh field produced by the same pipeline (mean 7.3%, SD 1.88% [23]). A "
         "public script nevertheless carries a `--fastGWA-mlm-binary` flag the text does not "
         "describe; if it reflected the released statistics the bounds would tighten, not "
         "loosen, and both contrasts are unaffected.")

sub1(OLD_L, NEW_L, "Limitations: TMFI scale reclassified")

# ---------------------------------------------------------------- 3. Reference [23]
OLD_R = ("22. Su C-Y, Butler-Laporte G. Positive-control Mendelian randomization highlights "
         "power constraints in disease-mortality GWAS. *medRxiv*. 2026. "
         "doi:[10.64898/2026.05.29.26354472](https://doi.org/10.64898/2026.05.29.26354472). "
         "Preprint.")

NEW_R = (OLD_R +
         "\n23. UK Biobank. Data-field 22435: muscle fat infiltration (anterior thigh). UK "
         "Biobank Showcase. Available from: "
         "https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=22435. "
         "Accessed 20 September 2026.")

sub1(OLD_R, NEW_R, "Reference [23] UK Biobank field 22435")

# ---------------------------------------------------------------- assertions
def W(s):
    return len(re.findall(r"[^\s]+", s))

# --- word count (TRACKED = the metric used throughout this project:
#     Introduction .. Main tables, counting every word incl. ### subheadings)
body = text[text.index("## Introduction"):text.index("## Main tables")]
tracked = W(body)

BAK = P + ".bak_unit2026-09-20"
with io.open(BAK, "r", encoding="utf-8", newline="") as f:
    btext = f.read()
bbody = btext[btext.index("## Introduction"):btext.index("## Main tables")]
btracked = W(bbody)
print("\n  TRACKED body words = %d -> %d / 4500  (margin %d, delta %+d)"
      % (btracked, tracked, 4500 - tracked, tracked - btracked))

# --- abstract unchanged (compare before/after with the same rule)
def ab_words(t):
    a = t.split("## Abstract", 1)[1].split("\n## Introduction", 1)[0]
    a = re.sub(r"^###\s*.*$", "", a, flags=re.M)
    a = re.sub(r"^\*\*Keywords:\*\*.*$", "", a, flags=re.M)
    return W(a)

print("  abstract words     = %d -> %d / 400" % (ab_words(btext), ab_words(text)))
if ab_words(text) != ab_words(btext):
    sys.exit("ABORT: abstract word count changed (nothing written)")

# --- must-keep
KEEP = [
    "0.048", "0.031", "0.090", "0.038", "0.143", "0.00075",
    "22.12", "35.23", "2.03", "1.217", "1.116", "1.327",
    "0.983", "46,347", "134,627", "106,697", "499,362", "237,944",
    "42", "93", "80%", "descriptively larger", "internal comparator",
    "structural ceiling of 0.890 within ~2 SE", "Supplementary Table S13",
    "Eight boundaries", "bounded null", "Mendelian randomization",
    "raw-scale coefficients reported alongside",
]
missing = [k for k in KEEP if k not in text]
print("  must-keep: %d/%d present" % (len(KEEP) - len(missing), len(KEEP)))
if missing:
    sys.exit("ABORT: lost %r" % missing)

# --- must-go
GONE = [
    "no verifiable physical unit",
    "unresolved preprocessing discrepancy",
    "should not be translated into physical units",
    "conditional on that confirmation",
]
left = [g for g in GONE if g in text]
print("  must-go  : %d/%d absent" % (len(GONE) - len(left), len(GONE)))
if left:
    sys.exit("ABORT: still present %r" % left)

# --- structure
refs = re.findall(r"(?m)^(\d+)\.\s", text[text.index("## References"):])
refs = [int(x) for x in refs]
print("  references         = %d (contiguous: %s)" % (len(refs), refs == list(range(1, len(refs) + 1))))
if len(refs) != 23 or refs != list(range(1, 24)):
    sys.exit("ABORT: reference numbering wrong")

# --- title / keywords / floats
title = text.split("\n")[0].lstrip("# ").strip()
print("  title words        = %d / 17" % W(title))
kw = re.search(r"\*\*Keywords:\*\*(.+)", text).group(1)
print("  keywords           = %d / 4-6" % len([k for k in kw.split(";") if k.strip()]))
tabs = re.findall(r"(?m)^### Table \d+\.", text[text.index("## Main tables"):text.index("## Figure legends")])
figs = re.findall(r"(?m)^### Figure \d+\.", text[text.index("## Figure legends"):text.index("## Supplementary material")])
print("  main floats        = %d (tables %d + figures %d) / 8" % (len(tabs) + len(figs), len(tabs), len(figs)))
if len(tabs) + len(figs) != 8:
    sys.exit("ABORT: main float count changed (nothing written)")

if tracked > 4500:
    sys.exit("ABORT: body over 4500 words (nothing written)")

# ---------------------------------------------------------------- write (CRLF)
if text == orig:
    sys.exit("ABORT: no change")

data = text.replace("\r\n", "\n").replace("\n", "\r\n").encode("utf-8")
with io.open(P, "wb") as f:
    f.write(data)
print("\n  WRITTEN %d bytes" % len(data))
