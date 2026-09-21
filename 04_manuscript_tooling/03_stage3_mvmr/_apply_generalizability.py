from pathlib import Path

P = Path(r"D:\SX\TMFI_CPI_MR\03_stage3_mvmr\Complete_Manuscript_Draft_TMFI_BMI_CPI_2026-09-18.md")
t = P.read_text(encoding="utf-8")

OLD = "Each of these specifies what a replication study must fix, and none alters the primary contrast."
NEW = ("Each of these specifies what a replication study must fix, and none alters the primary "
       "contrast. Generalizability is correspondingly limited: all four source GWAS were conducted "
       "in European-ancestry participants, so the estimates should not be extrapolated to other "
       "ancestries without replication; they describe the effect of a lifetime-averaged genetic "
       "difference in the exposure rather than an intervention applied at a particular age or over "
       "a particular duration; and they apply to a shift across the full genetic range of the "
       "exposure, which is not the same quantity as a clinically attainable change. None of these "
       "three extrapolations can be assessed from the present data.")

c = t.count(OLD)
if c != 1:
    raise SystemExit("ABORT: %d occurrences" % c)
t = t.replace(OLD, NEW, 1)
P.write_text(t, encoding="utf-8")
print("ok")
