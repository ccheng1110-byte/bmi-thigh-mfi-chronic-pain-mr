"""Stage 6 - Step 1: extract the 93 MVMR instrument SNPs from FinnGen R12 summary stats.

Design notes
------------
* FinnGen R12 is GRCh38; our LD panel / instrument coordinates are b37.
  Position joins would silently return zero rows, so we anchor on **rsid**
  (FinnGen `rsids` column) and carry our own b37 coordinates through.
* Files are ~800 MB gzipped (~20 M rows). We stream in binary, split on at most
  7 tabs per line, and only decode rows that match -- keeps the scan fast and
  memory flat.
* NEVER point this at `meta_analysis/` (FinnGen x UKB x Estonia meta); those
  reintroduce the sample overlap we are trying to escape.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import os
import sys
import time
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
INCOMING = os.path.join(HERE, "incoming")
EXTRACTED = os.path.join(HERE, "extracted")
# union of BOTH frozen instrument sets (42-SNP TMFI + 93-SNP MVMR = 109 rsids).
# See build_instrument_master.py: the 42-SNP set is NOT nested in the 93-SNP set.
UNION = os.path.join(HERE, "instruments_master.csv")

PHENOS = ["PAIN", "M13_FIBROMYALGIA", "M13_LIMBPAIN"]

# keep in sync with 06_stage6_finngen/FinnGen外部复制执行清单_2026-09-20.md
EXPECTED_N = {
    "PAIN": (237944, 261418),
    "M13_FIBROMYALGIA": (3623, 357549),
    "M13_LIMBPAIN": (45476, 357549),
}

RSID_SPLIT = b",;|"


def load_instruments(path):
    """Return (dict rsid -> row, list of rsids)."""
    with io.open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    by_rsid = {}
    for r in rows:
        by_rsid[r["rsid"].strip()] = r
    return by_rsid, list(by_rsid)


def split_any(raw: bytes):
    """FinnGen sometimes lists several rsids for one multiallelic row."""
    out = []
    cur = raw
    for sep in RSID_SPLIT:
        cur = cur.replace(bytes([sep]), b" ")
    for tok in cur.split():
        tok = tok.strip()
        if tok.startswith(b"rs"):
            out.append(tok.decode("ascii"))
    return out


def scan(pheno, wanted, path):
    """Stream one FinnGen file and return {rsid: [raw row dicts]}."""
    hits = defaultdict(list)
    n_lines = 0
    t0 = time.time()
    last = t0
    with gzip.open(path, "rb") as fh:
        header = fh.readline().decode("ascii").rstrip("\n").split("\t")
        idx = {c.lstrip("#"): i for i, c in enumerate(header)}
        # maxsplit must be large enough to expose the last column we read
        # (af_alt_controls) but no larger -- keeps the scan cheap
        last_idx = max(idx["rsids"], idx["af_alt_controls"])
        maxsplit = last_idx + 1
        n_rs = 0
        for raw in fh:
            n_lines += 1
            if n_lines % 4_000_000 == 0:
                now = time.time()
                sys.stdout.write(
                    "    %8d rows  %5.1fs  (+%.1fs)\n" % (n_lines, now - t0, now - last)
                )
                sys.stdout.flush()
                last = now
            parts = raw.split(b"\t", maxsplit)
            if len(parts) <= last_idx:
                continue
            rs = parts[idx["rsids"]]
            if not rs.startswith(b"rs"):
                continue
            n_rs += 1
            toks = split_any(rs)
            if not toks:
                continue
            for t in toks:
                if t in wanted:
                    hits[t].append(
                        {
                            "chrom": parts[idx["chrom"]].decode(),
                            "pos": parts[idx["pos"]].decode(),
                            "ref": parts[idx["ref"]].decode(),
                            "alt": parts[idx["alt"]].decode(),
                            "rsids": parts[idx["rsids"]].decode(),
                            "nearest_genes": parts[idx["nearest_genes"]].decode(),
                            "pval": parts[idx["pval"]].decode(),
                            "beta": parts[idx["beta"]].decode(),
                            "sebeta": parts[idx["sebeta"]].decode(),
                            "af_alt": parts[idx["af_alt"]].decode(),
                            "af_alt_cases": parts[idx["af_alt_cases"]].decode(),
                            # last column on the line -- strip the line terminator
                            "af_alt_controls": parts[idx["af_alt_controls"]].decode().strip(),
                        }
                    )
    print("    rows with an rs-prefixed rsid: %d" % n_rs)
    return hits, n_lines, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("phenos", nargs="*", default=PHENOS)
    args = ap.parse_args()

    by_rsid, rsids = load_instruments(UNION)
    wanted = set(rsids)
    print("instruments loaded: %d (union file %s)" % (len(rsids), os.path.basename(UNION)))
    os.makedirs(EXTRACTED, exist_ok=True)

    summary = {}
    for ph in args.phenos:
        path = os.path.join(INCOMING, "finngen_R12_%s.gz" % ph)
        if not os.path.isfile(path):
            print("\n[skip] %s -- file not found: %s" % (ph, path))
            continue
        print("\n=== %s  (%.0f MB)" % (ph, os.path.getsize(path) / 2**20))
        hits, n_lines, dt = scan(ph, wanted, path)
        missing = [r for r in rsids if r not in hits]
        multi = {k: len(v) for k, v in hits.items() if len(v) > 1}
        print("    scanned %d rows in %.1fs" % (n_lines, dt))
        print("    matched %d / %d ; missing %d ; multiallelic rows %d"
              % (len(hits), len(rsids), len(missing), len(multi)))
        if missing:
            print("    missing rsids: %s" % ", ".join(missing[:20]))

        out_tsv = os.path.join(EXTRACTED, "%s_raw.tsv" % ph)
        cols = ["rsid", "chrom", "pos", "ref", "alt", "beta", "sebeta", "pval",
                "af_alt", "af_alt_cases", "af_alt_controls", "nearest_genes", "n_rows"]
        with io.open(out_tsv, "w", encoding="utf-8", newline="") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rsids:  # keep instrument order stable
                for h in hits.get(r, []):
                    fh.write("\t".join([
                        r, h["chrom"], h["pos"], h["ref"], h["alt"],
                        h["beta"], h["sebeta"], h["pval"],
                        h["af_alt"], h["af_alt_cases"], h["af_alt_controls"],
                        h["nearest_genes"], str(len(hits.get(r, []))),
                    ]) + "\n")
        print("    wrote %s (%d data rows)" % (out_tsv, sum(len(v) for v in hits.values())))
        summary[ph] = {
            "n_scanned": n_lines,
            "n_matched": len(hits),
            "n_missing": len(missing),
            "missing_rsids": missing,
            "multiallelic": multi,
            "seconds": round(dt, 1),
        }

    with io.open(os.path.join(EXTRACTED, "extract_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    print("\ndone -> %s" % os.path.join(EXTRACTED, "extract_summary.json"))


if __name__ == "__main__":
    main()
