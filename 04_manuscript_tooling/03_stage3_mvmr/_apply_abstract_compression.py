# -*- coding: utf-8 -*-
"""
Compress the structured abstract to create headroom under the 350-word cap
(BMC Medicine: abstract <= 350 words, structured Background/Methods/Results/
Conclusions, no references, abbreviations minimised).

Numbers, P values, sample sizes and every substantive claim are preserved.
Every replacement is asserted to occur exactly once; abort otherwise.
CRLF restored after writing.
"""
import io, os, re, shutil, sys

D = os.path.dirname(os.path.abspath(__file__))
MAN = os.path.join(D, "Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md")
BAK = MAN + ".bak_abs2026-09-20"
if not os.path.exists(BAK):
    shutil.copy2(MAN, BAK)

with io.open(MAN, "r", encoding="utf-8") as f:
    text = f.read()

REPL = []

# ------------------------------------------------------------ Background 42 -> 37
REPL.append((
    "Thigh muscle fat infiltration (TMFI) and body mass index (BMI) may both influence chronic pain "
    "intensity (CPI) but are hard to separate. BMI was pre-specified as an internal comparator, so a null "
    "TMFI estimate could be read against an explicit power bound.",
    "Thigh muscle fat infiltration (TMFI) and body mass index (BMI) may both influence chronic pain "
    "intensity (CPI). BMI was pre-specified as an internal comparator, so a null TMFI estimate could be "
    "read against an explicit power bound."
))

# ------------------------------------------------------------ Methods 95 -> 85
REPL.append((
    "We used summary statistics for TMFI (n=46,347), BMI (GIANT 2015) and two UK Biobank chronic pain "
    "intensity outcomes (n=134,627 and 106,697). We first estimated the total effect of TMFI using 42 "
    "independent genetic variants, then estimated TMFI and BMI jointly using 93 so each effect is adjusted "
    "for the other, with standard diagnostics. Exposure\u2013outcome overlap was quantified genetically and "
    "carried through every estimate; results are reported against the smallest effect detectable at 80% "
    "power rather than significance. Both stages were then repeated against three FinnGen endpoints (up to "
    "499,362 participants) with no UK Biobank overlap.",
    "Summary statistics were used for TMFI (n=46,347), BMI (GIANT 2015) and two UK Biobank chronic pain "
    "intensity outcomes (n=134,627 and 106,697). The total TMFI effect was estimated with 42 variants, "
    "then TMFI and BMI jointly with 93 so each is adjusted for the other. Exposure\u2013outcome overlap was "
    "quantified genetically and carried through every estimate; results are reported against the smallest "
    "effect detectable at 80% power rather than significance. Both stages were repeated against three "
    "FinnGen endpoints (up to 499,362 participants) with no UK Biobank overlap."
))

# ------------------------------------------------------------ Results 159 -> 151
REPL.append((
    "instruments remained strong after conditioning (F 22.1 and 35.2) and showed no directional pleiotropy. "
    "The study could detect a TMFI effect of 0.048 SD at 80% power, 1.9-fold smaller than the BMI effect "
    "recovered.",
    "instruments remained strong after conditioning (F 22.1 and 35.2), with no directional pleiotropy. The "
    "design could detect 0.048 SD at 80% power, 1.9-fold smaller than the BMI effect recovered."
))
REPL.append((
    "Exposure\u2013outcome overlap (intercept 0.052, P=1.2\u00d710\u207b\u00b9\u00b9) moved estimates by "
    "\u22640.0015 SD per SD and could not reverse either finding.",
    "Exposure\u2013outcome overlap (intercept 0.052, P=1.2\u00d710\u207b\u00b9\u00b9) moved estimates by "
    "\u22640.0015 SD per SD and could not reverse either."
))
REPL.append((
    "Both patterns held in FinnGen: TMFI null (limb pain OR 0.983 per SD, P=0.53), BMI associated with limb "
    "pain (1.217, 1.116\u20131.327, P=9\u00d710\u207b\u2076).",
    "Both held in FinnGen: TMFI null for limb pain (OR 0.983 per SD, P=0.53), BMI associated with limb pain "
    "(1.217, 1.116\u20131.327, P=9\u00d710\u207b\u2076)."
))

# ------------------------------------------------------------ Conclusions 54 -> 52
REPL.append((
    "Against an explicit power bound, an internal comparator that behaved as expected, replication in an "
    "independent cohort and a quantified overlap, this is a bounded null:",
    "Against an explicit power bound, an internal comparator that behaved as expected, an independent "
    "replication and a quantified overlap, this is a bounded null:"
))

for old, new in REPL:
    k = text.count(old)
    if k != 1:
        sys.exit("ABORT: occurrence=%d for pattern: %r" % (k, old[:80]))
    text = text.replace(old, new, 1)

# ------------------------------------------------------------ word-count report
W = lambda s: len(re.findall(r"[^\s]+", s))
ab = text.split("## Abstract", 1)[1].split("\n## Introduction", 1)[0]
mk = re.search(r"\*\*Keywords:\*\*(.*)$", ab, re.S)
kw_line = mk.group(1).strip()
n_kw = [k.strip() for k in kw_line.split(";") if k.strip()]
print("keywords: %d entries, %d tokens  -> %s" % (len(n_kw), W(kw_line), "; ".join(n_kw)))
print("BMC Medicine requires 3-10 keywords:", "OK" if 3 <= len(n_kw) <= 10 else "VIOLATION")

total = 0
for s in re.split(r"\n### ", ab):
    if not s.strip():
        continue
    name = s.split("\n", 1)[0].strip()
    body = "\n".join(s.split("\n")[1:])
    body = re.sub(r"\*\*Keywords:\*\*.*$", "", body, flags=re.S)
    n = W(body)
    total += n
    print("  %-14s %4d" % (name, n))
print("ABSTRACT PROSE : %d words (cap 350 -> headroom %d)" % (total, 350 - total))
if total > 350:
    sys.exit("ABORT: abstract over 350")

# values that must survive compression
must = ["0.008", "0.558", "0.005", "0.721", "0.090", "0.038", "0.143", "0.00075",
        "\u22120.002", "0.898", "22.1", "35.2", "0.048", "0.052", "1.2", "0.0015",
        "1.015", "0.58", "1.001", "0.98", "1.275", "1.134", "0.983", "0.53",
        "1.217", "1.116", "1.327", "46,347", "134,627", "106,697", "499,362",
        "42", "93", "80%", "0.05"]
missing = [m for m in must if m not in ab]
print("missing key values:", missing if missing else "none")
if missing:
    sys.exit("ABORT: key value lost from abstract")

raw = text.encode("utf-8").replace(b"\n", b"\r\n")
with open(MAN, "wb") as f:
    f.write(raw)
print("written:", MAN)
