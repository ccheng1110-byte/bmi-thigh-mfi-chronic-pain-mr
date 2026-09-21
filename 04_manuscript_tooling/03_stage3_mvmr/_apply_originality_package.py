# -*- coding: utf-8 -*-
"""
Apply the originality package to the manuscript:
  1. Introduction positioning sentence (refs 19,20,21)
  2. Discussion "Relation to existing evidence" extension (+ pointer to S13)
  3. Discussion "Implications" clause (ref 22, comparator-based power calibration)
  4. Supplementary list: new item 14 = Supplementary Table S13
  5. References: append [19]-[22]

Every replacement is asserted to occur exactly once; the script aborts otherwise.
CRLF is restored after writing.
"""
import io, os, re, shutil, sys, unicodedata

MAN = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md")
BAK = MAN + ".bak_origpack2026-09-20"

if not os.path.exists(BAK):
    shutil.copy2(MAN, BAK)

with io.open(MAN, "r", encoding="utf-8") as f:
    text = f.read()

n0 = len(text)

# ---------------------------------------------------------------- replacements
REPL = []

# 1. Introduction positioning
REPL.append((
    "A null result reported without either is descriptively correct but interpretively uninformative.\n",
    "A null result reported without either is descriptively correct but interpretively uninformative. "
    "Genetically predicted adiposity is already an established cause of site-specific pain [19,20], and the "
    "genetic architecture of muscle fatty infiltration has recently been characterised, although with "
    "infiltration treated as the outcome rather than as an exposure on pain [21]. What neither literature "
    "provides is the joint test: a muscle-compartment exposure and an adiposity exposure estimated in one "
    "model, so that a muscle effect can be read against a bound rather than against a P value.\n"
))

# 2. Discussion: relation to existing evidence
REPL.append((
    "distinguishing those possibilities requires the joint conditioning reported here.\n",
    "distinguishing those possibilities requires the joint conditioning reported here. The same distinction "
    "separates this work from the adiposity-to-pain literature: those analyses established that higher BMI "
    "causes back, knee and hip pain [19,20], and the BMI estimates reported here agree with them in "
    "direction and magnitude, but they could not ask whether a muscle compartment contributes independently, "
    "because no comparator exposure was estimated in the same model. A structured comparison with the "
    "nearest published analyses is given in Supplementary Table S13.\n"
))

# 3. Discussion: implications
REPL.append((
    "absent both, a null cannot be distinguished from low power. Substantively,",
    "absent both, a null cannot be distinguished from low power. Comparator-based power calibration of this "
    "kind is only beginning to be applied explicitly [22]. Substantively,"
))

# 4. Supplementary item 14
REPL.append((
    "13. **Supplementary Table S12** — STROBE-MR reporting checklist: all 20 items with the manuscript "
    "section and line reference satisfying each one.\n",
    "13. **Supplementary Table S12** — STROBE-MR reporting checklist: all 20 items with the manuscript "
    "section and line reference satisfying each one.\n"
    "14. **Supplementary Table S13** — Originality positioning: the present design and results compared "
    "with the nearest published Mendelian randomization analyses of adiposity, muscle and pain, tabulated "
    "by estimator, presence of a comparator exposure in the same model, reported power bound, "
    "quantification of exposure–outcome overlap and availability of a replication cohort with no "
    "participant overlap.\n"
))

# 5. References 19-22
REPL.append((
    "18. FinnGen. FinnGen documentation: data download (R12 release). 2025. Available from: "
    "https://finngen.gitbook.io/documentation/data-download. Accessed 20 September 2026.\n",
    "18. FinnGen. FinnGen documentation: data download (R12 release). 2025. Available from: "
    "https://finngen.gitbook.io/documentation/data-download. Accessed 20 September 2026.\n"
    "19. Elgaeva EE, Tsepilov Y, Freidin MB, Williams FMK, Aulchenko Y, Suri P. ISSLS Prize in Clinical "
    "Science 2020. Examining causal effects of body mass index on back pain: a Mendelian randomization "
    "study. *European Spine Journal*. 2020;29(4):686–691. "
    "doi:[10.1007/s00586-019-06224-6](https://doi.org/10.1007/s00586-019-06224-6).\n"
    "20. Chen X, Tang H, Lin J, Zeng R. Causal relationships of obesity on musculoskeletal chronic pain: a "
    "two-sample Mendelian randomization study. *Frontiers in Endocrinology*. 2022;13:971997. "
    "doi:[10.3389/fendo.2022.971997](https://doi.org/10.3389/fendo.2022.971997).\n"
    "21. Tan J, Ding Z, Zheng J, Zhang J, Chen X, Li Z, Shi L, Chen J, Sun Y. Potential regulators and "
    "metabolic networks of muscle fatty infiltration: genomic and radiomics investigation based on 33 300 "
    "participants. *International Journal of Surgery*. 2025. Advance online publication. "
    "doi:[10.1097/JS9.0000000000003523](https://doi.org/10.1097/JS9.0000000000003523).\n"
    "22. Su C-Y, Butler-Laporte G. Positive-control Mendelian randomization highlights power constraints in "
    "disease-mortality GWAS. *medRxiv*. 2026. "
    "doi:[10.64898/2026.05.29.26354472](https://doi.org/10.64898/2026.05.29.26354472). Preprint.\n"
))

for old, new in REPL:
    k = text.count(old)
    if k != 1:
        sys.exit("ABORT: occurrence=%d for pattern starting: %r" % (k, old[:70]))
    text = text.replace(old, new, 1)

# ---------------------------------------------------------------- sanity checks
for tag in ["[19]", "[20]", "[21]", "[22]", "Supplementary Table S13"]:
    print("occurrences %-28s : %d" % (tag, text.count(tag)))

# abstract word count (unchanged by this patch, but verify the 350 cap).
# Count prose only: exclude the "###" section labels, the "**Keywords:**" label
# and the keyword list itself (journals count keywords separately).
ab = text.split("## Abstract", 1)[1].split("\n## Introduction", 1)[0]
tok = len(re.findall(r"[^\s]+", ab))
head = sum(len(re.findall(r"[^\s]+", l)) for l in ab.split("\n")
           if l.strip().startswith("#"))
mk = re.search(r"\*\*Keywords:\*\*(.*)$", ab, re.S)
kw = len(re.findall(r"[^\s]+", mk.group(1))) if mk else 0
words = tok - head - 1 - kw
print("abstract tokens total :", tok)
print("abstract prose words  :", words, "(cap 350)")
print("abstract keywords     :", kw)
if words > 350:
    sys.exit("ABORT: abstract prose over 350 words")

body = text.split("\n## Introduction", 1)[1].split("\n## Main tables", 1)[0]
print("body words (Intro..Conclusion):", len(re.findall(r"[^\s]+", body)))
print("total lines:", text.count("\n") + 1)
print("chars %d -> %d" % (n0, len(text)))

raw = text.encode("utf-8").replace(b"\n", b"\r\n")
with open(MAN, "wb") as f:
    f.write(raw)
print("written:", MAN)
