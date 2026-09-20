# GitHub Upload Audit and Fix — bmi-thigh-mfi-chronic-pain-mr

**Repository:** https://github.com/ccheng1110-byte/bmi-thigh-mfi-chronic-pain-mr
**Local source:** `d:\SX\TMFI_CPI_MR\06_code_repository_2026-09-20\`
**Audited:** 2026-09-21
**Method:** remote file tree via GitHub API; per-file comparison by Git blob SHA-1 of the actual file bytes

---

## 1. Audit result

| Item | Result |
|---|---|
| Local files (expected) | **194** |
| Files on GitHub | **175** |
| **Missing on GitHub** | **19** |
| Extra files (on GitHub but not local) | **0** |
| Content byte-identical | **174 / 175** |
| Line-ending (CRLF/LF) mismatches | **0** |
| Content differences | **1** — `README.md` only |

The directories are all present and correctly named: `01_analysis_scripts` (61), `02_frozen_outputs` (39 of 56), `03_figure_code_and_outputs` (13), `04_manuscript_tooling` (40), `05_environment_and_logs` (19), plus root files.

**Bottom line: what you uploaded is correct — nothing was corrupted or misnamed. You are missing 19 files and one stray line in README.md.**

---

## 2. Issue A — README.md has one extra line

The file on GitHub begins with this line, which is NOT part of the README:

```
[README.md](https://github.com/user-attachments/files/32437795/README.md)
```

This is a GitHub-generated attachment link, inserted when the file is pasted/dropped into GitHub's online editor rather than uploaded as a file. It also converted the file from CRLF to LF. Everything else (the remaining 47 lines) is byte-identical to the local version.

**Fix:** the corrected `README.md` in this folder has that line removed. Upload it to overwrite the one on GitHub.
(Alternative: open README.md on GitHub, click the pencil icon, delete line 1, commit.)

---

## 3. Issue B — 19 missing files

Root cause: GitHub's web upload caps at 100 files per batch, so these were skipped when the batches were assembled.

| # | File | Note |
|---|---|---|
| 1 | `MANIFEST.tsv` | **Important** — README references it; it is the MD5 inventory reviewers use |
| 2 | `.gitignore` | Hidden file — verify it appears in the repo root afterwards |
| 3–5 | `02_frozen_outputs/S8__R_validation_ivw.csv`<br>`S8__R_validation_pleiotropy.csv`<br>`S8__R_validation_strength.csv` | Supplementary Table S8 |
| 6 | `02_frozen_outputs/S9__stage5_overlap_bias.csv` | Supplementary Table S9 |
| 7–19 | `02_frozen_outputs/_extra/` (13 files) | Additional frozen inputs needed to reproduce: PLINK clump results, `mvmr_union_exposures_all.csv`, stage3b summary, stage4 LDSC results, `instruments_master.csv`, etc. |

`_extra/` was skipped entirely as a subfolder — this is why the count is 13.

---

## 4. How to upload this folder

1. Open https://github.com/ccheng1110-byte/bmi-thigh-mfi-chronic-pain-mr
2. Click **Add file → Upload files**
3. Drag the **entire contents** of `07_missing_files_for_github_2026-09-21\` (all 20 items, well under the 100-file limit) into the upload area
4. Commit message: `Add missing MANIFEST, supplementary S8/S9, _extra inputs; fix README header`
5. Commit directly to `main`

**Two things to check afterwards:**
- `.gitignore` is visible in the repo root (Windows Explorer shows it normally; some upload methods skip dotfiles)
- README.md no longer starts with `[README.md](...)`

---

## 5. Expected final state

175 + 19 missing = **194 files**, matching the local repository exactly.
`MANIFEST.tsv` (193 entries + header) will then cover every file in the repository.

After that, create the **v1.0 release** so Zenodo can mint a DOI:
GitHub → Releases → Draft a new release → tag `v1.0`; then enable the repository in Zenodo and paste the DOI into the 3 placeholders in the manuscript's Code/Data availability section.
