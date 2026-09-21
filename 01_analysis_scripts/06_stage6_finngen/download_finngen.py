#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""FinnGen R12 summary-statistics parallel downloader (resumable, md5-verified).

Usage
-----
    python download_finngen.py PAIN
    python download_finngen.py PAIN M13_FIBROMYALGIA M13_LIMBPAIN
    python download_finngen.py PAIN --threads 16

Files land in  06_stage6_finngen/incoming/finngen_R12_<PHENO>.gz

Why parallel: a single connection measured 0.23 MB/s from this network
(783 MB -> ~57 min); 8 concurrent Range requests gave 1.25 MB/s (~10 min).
Raise --threads if your link allows, lower it if the server throttles.

Resume: each chunk is written to its own .part file. Re-running skips chunks
that are already complete, so an interrupted download continues instead of
restarting.

Only ever fetch from summary_stats/release/ -- the meta_analysis/ module
contains FinnGen x UKB x Estonian meta-analysed results and would reintroduce
the UK Biobank sample overlap this replication exists to remove.
"""
import argparse
import base64
import hashlib
import json
import os
import socket
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request

BUCKET = "finngen-public-data-r12"
PREFIX = "summary_stats/release"
BASE_HTTP = f"https://storage.googleapis.com/{BUCKET}/{PREFIX}/"
API_OBJ = f"https://storage.googleapis.com/storage/v1/b/{BUCKET}/o/"
INCOMING = os.path.join(os.path.dirname(os.path.abspath(__file__)), "incoming")

socket.setdefaulttimeout(120)
UA = {"User-Agent": "Mozilla/5.0"}


def object_meta(pheno):
    """Exact size + md5 straight from the GCS object metadata."""
    name = f"{PREFIX}/finngen_R12_{pheno}.gz"
    url = API_OBJ + urllib.parse.quote(name, safe="")
    req = urllib.request.Request(url, headers=UA)
    return json.loads(urllib.request.urlopen(req).read())


def chunk_worker(url, tidx, start, end, part, lock, counter, t0):
    """Download [start, end] into `part`; skip if already complete."""
    want = end - start + 1
    if os.path.isfile(part) and os.path.getsize(part) == want:
        with lock:
            counter[0] += want
        return True
    tmp = part + ".dl"
    for attempt in range(5):
        try:
            req = urllib.request.Request(
                url, headers=dict(UA, Range=f"bytes={start}-{end}")
            )
            r = urllib.request.urlopen(req)
            got = 0
            with open(tmp, "wb") as fh:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    fh.write(b)
                    got += len(b)
                    with lock:
                        counter[0] += len(b)
            if got != want:
                raise IOError(f"short read {got} != {want}")
            os.replace(tmp, part)
            return True
        except Exception as exc:  # noqa: BLE001
            if os.path.isfile(tmp):
                os.remove(tmp)
            if attempt == 4:
                print(f"   [chunk {tidx}] FAILED: {exc}", flush=True)
                return False
            time.sleep(2 ** attempt)
    return False


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return base64.b64encode(h.digest()).decode()


def download(pheno, threads):
    meta = object_meta(pheno)
    size = int(meta["size"])
    want_md5 = meta.get("md5Hash")
    url = BASE_HTTP + f"finngen_R12_{pheno}.gz"
    os.makedirs(INCOMING, exist_ok=True)
    dst = os.path.join(INCOMING, f"finngen_R12_{pheno}.gz")

    if os.path.isfile(dst) and os.path.getsize(dst) == size:
        print(f"{pheno}: {dst} already complete ({size:,} bytes)", flush=True)
    else:
        chunk = max(4 << 20, (size + threads - 1) // threads)
        bounds = []
        off = 0
        while off < size:
            end = min(off + chunk - 1, size - 1)
            bounds.append((off, end))
            off = end + 1
        print(
            f"{pheno}: {size:,} bytes ({size/2**20:.1f} MB) in "
            f"{len(bounds)} chunks x {threads} threads",
            flush=True,
        )
        lock = threading.Lock()
        counter = [0]
        t0 = time.time()
        results = []
        slots = [None] * len(bounds)

        def run(i):
            s, e = bounds[i]
            slots[i] = chunk_worker(
                url, i, s, e, f"{dst}.part{i:04d}", lock, counter, t0
            )

        # bounded concurrency
        idx = 0
        live = []
        while idx < len(bounds) or live:
            while idx < len(bounds) and len(live) < threads:
                t = threading.Thread(target=run, args=(idx,))
                t.start()
                live.append(t)
                idx += 1
            time.sleep(0.4)
            live = [t for t in live if t.is_alive()]
            done = counter[0]
            el = max(time.time() - t0, 1e-9)
            sys.stdout.write(
                f"\r   {done/2**20:8.1f} / {size/2**20:.1f} MB "
                f"({100*done/size:5.1f}%)  {done/2**20/el:6.2f} MB/s"
            )
            sys.stdout.flush()
        print(flush=True)
        if not all(slots):
            print(f"{pheno}: INCOMPLETE - re-run to resume the missing chunks")
            return False

        with open(dst, "wb") as out:
            for i in range(len(bounds)):
                p = f"{dst}.part{i:04d}"
                with open(p, "rb") as fh:
                    for b in iter(lambda: fh.read(1 << 22), b""):
                        out.write(b)
                os.remove(p)
        print(f"{pheno}: assembled -> {dst}", flush=True)

    if want_md5:
        got = md5_of(dst)
        ok = got == want_md5
        print(f"{pheno}: md5 {'OK' if ok else 'MISMATCH'}  ({got} vs {want_md5})")
        return ok
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("phenos", nargs="+")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    allok = True
    for ph in a.phenos:
        allok &= download(ph, a.threads)
    print("\nALL OK" if allok else "\nSOME FAILED - see messages above")
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
