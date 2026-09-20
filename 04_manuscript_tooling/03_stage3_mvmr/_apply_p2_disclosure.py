# -*- coding: utf-8 -*-
"""
P2 disclosure upgrade.

1. Methods  L49 : "may have been missed -- a disclosed limitation"
                 -> bounded shortfall, cross-referenced to Limitations
2. Limit.   L99 : open-ended "eligible instruments outside it may have been missed"
                 -> bounded (49 loci / 42 visible / 7 blind / +8.4% R2 / MDE 0.048->0.046)
3. Concl.  L101 : drop "uses a denser LD panel" (contradicted by the P2 probe)
4. Suppl.       : add S17 = full LD coverage audit

Every substitution is asserted occurrence == 1 before writing.
"""
import io
import os
import re
import sys

P = r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md"
BAK = P + ".bak_p2_2026-09-20"

if not os.path.exists(BAK):
    with io.open(P, "r", encoding="utf-8", newline="") as f:
        raw = f.read()
    with io.open(BAK, "w", encoding="utf-8", newline="") as f:
        f.write(raw)
    print("backup written:", BAK)
else:
    with io.open(BAK, "r", encoding="utf-8", newline="") as f:
        raw = f.read()
    print("restored from backup:", BAK)

text = raw
assert "\r\n" in text, "ABORT: no CRLF in source"


def W(s):
    return len(re.findall(r"[^\s]+", s))


def sub1(old, new, tag):
    global text
    n = text.count(old)
    if n != 1:
        sys.exit("ABORT [%s]: occurrence = %d (expected 1)\n  %r" % (tag, n, old[:120]))
    text = text.replace(old, new)
    print("  ok  %-28s %+d words" % (tag, W(new) - W(old)))


# ---------------------------------------------------------------- 1. Methods
OLD_M = ("the reference panel contained only 22.75% of TMFI candidates, "
         "so eligible signals outside it may have been missed — a disclosed limitation.")
NEW_M = ("the reference panel contained only 22.75% of TMFI candidates, "
         "so eligible signals outside it may have been missed; the shortfall is "
         "bounded at seven instruments (Limitations).")
sub1(OLD_M, NEW_M, "L49 Methods")

# ------------------------------------------------------------- 2. Limitations
OLD_L = ("The LD reference panel (503 individuals, 1.84 million variants) contained "
         "only 22.75% of TMFI candidate variants, so eligible instruments outside it "
         "may have been missed.")
NEW_L = ("The LD reference panel (503 individuals, 1.84 million variants) contained "
         "only 22.75% of TMFI candidate variants, but these fall into 49 loci of "
         "which 42 contained a panel variant and yielded exactly the 42 instruments; "
         "the seven uncovered loci were weak (best P 2.3×10⁻⁹ to 4.8×10⁻⁸), so a "
         "denser panel could add at most seven instruments, raising the explained "
         "variance by 8.4% and moving the minimum detectable effect only from 0.048 "
         "to 0.046 SD per SD.")
sub1(OLD_L, NEW_L, "L99 Limitations")

# -------------------------------------------------------------- 3. Conclusion
OLD_C = ("The next decisive test is a replication that preserves this joint structure, "
         "uses a denser LD panel, and extends the conditioning to further correlated "
         "muscle traits.")
NEW_C = ("The next decisive test is a replication that preserves this joint structure "
         "and extends the conditioning to further correlated muscle traits.")
sub1(OLD_C, NEW_C, "L101 Conclusion")

# ------------------------------------------------------- 4. Supplementary S17
m = re.search(r"(?m)^16\. \*\*Supplementary [^\n]*\n", text)
if not m:
    sys.exit("ABORT: could not locate supplementary item 16")
S17 = ("17. **Supplementary Methods** — LD reference panel coverage audit: locus-level "
       "reconstruction of the 3,838 TMFI candidates (49 loci at a 500 kb gap), the tagging "
       "test of the 42 retained instruments against non-panel candidates together with its "
       "positive control, the r² noise floor of the 503-person panel, and the quantitative "
       "consequence of a denser panel (explained variance +8.4%; minimum detectable effect "
       "0.048 to 0.046 SD per SD).\r\n")
text = text[:m.end()] + S17 + text[m.end():]
print("  ok  %-28s S17 appended" % "supplementary list")

# ------------------------------------------------------------------ assertions
print("\n--- must-go (open-ended / contradicted wording) ---")
GONE = [
    "may have been missed.",            # Limitations: now bounded
    "a disclosed limitation",           # Methods: now bounded
    "uses a denser LD panel",           # Conclusion: contradicted by P2
]
for g in GONE:
    assert g not in text, "ABORT: still present %r" % g
print("  %d must-go strings absent" % len(GONE))

print("\n--- must-keep (P2 probe numbers now disclosed) ---")
KEEP = [
    "49 loci", "42 contained a panel variant", "exactly the 42 instruments",
    "seven uncovered loci", "2.3×10⁻⁹", "4.8×10⁻⁸",
    "at most seven instruments", "8.4%", "0.048 to 0.046 SD per SD",
    "22.75%", "503 individuals", "1.84 million variants",
    "bounded at seven instruments",
]
bad = [k for k in KEEP if k not in text]
if bad:
    sys.exit("ABORT: missing %r" % bad)
print("  %d must-keep strings present" % len(KEEP))

# --------------------------------------------------------------- word counts
i0 = text.index("## Introduction")
i1 = text.index("## Main tables")
tracked = W(text[i0:i1])

b = raw[raw.index("## Introduction"):raw.index("## Main tables")]
btracked = W(b)

print("\n--- word counts ---")
print("  TRACKED body (Introduction..Main tables) = %d -> %d  / 4500 (margin %d)"
      % (btracked, tracked, 4500 - tracked))

ab = text.split("## Abstract", 1)[1].split("\n## Introduction", 1)[0]
ab = re.sub(r"(?m)^###\s*.*$", "", ab)
print("  abstract = %d / 400" % W(ab))

if tracked > 4500:
    sys.exit("ABORT: TRACKED body over 4500 (nothing written)")

title = text.split("\n")[0].replace("# ", "").strip()
print("  title words = %d / 17" % W(title))

# structural
floats = len(re.findall(r"(?m)^### Table \d+\.", text)) + len(re.findall(r"(?m)^### Figure \d+\.", text))
refs = len(re.findall(r"(?m)^\d+\.\s", text[text.index("\n## References"):]))
kws = text.split("**Keywords:**", 1)[1].split("\n", 1)[0]
print("  main floats = %d / 8" % floats)
print("  references  = %d / 40" % refs)
print("  keywords    = %d / 4-6" % len([x for x in kws.split(",") if x.strip()]))

# ---------------------------------------------------------------------- write
with io.open(P, "w", encoding="utf-8", newline="") as f:
    f.write(text)

with io.open(P, "rb") as f:
    data = f.read()
crlf = data.count(b"\r\n")
bare = len(re.findall(rb"(?<!\r)\n", data))
print("\nwritten: %d bytes, CRLF=%d, bareLF=%d, lines=%d"
      % (len(data), crlf, bare, crlf + bare))
if bare:
    sys.exit("ABORT: bare LF introduced")
print("DONE")
