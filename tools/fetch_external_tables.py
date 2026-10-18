"""Download per-sample processed tables for the external NSCLC cohorts from EBI BioStudies.

Source: ArrayExpress/BioStudies mirrors of GEO series (NCBI is unreachable from this host):
  E-GEOD-68465 -> GSE68465 (Director's Challenge, LUAD, HG-U133A)
  E-GEOD-31210 -> GSE31210 (Japanese LUAD stage I-II, HG-U133 Plus 2)

Writes data/external/<gse>/tables/<sample>_sample_table.txt and a manifest with sizes + md5.
Idempotent: existing files with matching byte size are skipped.
"""
import concurrent.futures as cf
import hashlib
import json
import os
import sys
import time
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "work")
DEST_ROOT = os.path.join(ROOT, "data", "external")
BASE = "https://www.ebi.ac.uk/biostudies/files/"

STUDIES = {"E-GEOD-68465": "gse68465", "E-GEOD-31210": "gse31210"}


def file_list(acc):
    j = json.load(open(os.path.join(WORK, "_bs_%s.json" % acc), encoding="utf-8"))
    out = []

    def walk(n):
        if isinstance(n, dict):
            if "path" in n and "size" in n and str(n["path"]).endswith("sample_table.txt"):
                out.append((n["path"], int(n["size"])))
            for v in n.values():
                walk(v)
        elif isinstance(n, list):
            for v in n:
                walk(v)

    walk(j)
    return sorted(set(out))


def fetch(acc, dest_dir, item):
    path, size = item
    out = os.path.join(dest_dir, os.path.basename(path))
    if os.path.exists(out) and os.path.getsize(out) == size:
        return ("skip", path, size, None)
    url = BASE + acc + "/" + path
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                data = r.read()
            if len(data) != size:
                raise IOError("size mismatch %d != %d" % (len(data), size))
            with open(out, "wb") as fh:
                fh.write(data)
            return ("ok", path, size, hashlib.md5(data).hexdigest())
        except Exception as e:  # noqa: BLE001
            if attempt == 3:
                return ("fail", path, size, str(e))
            time.sleep(3 + 3 * attempt)


def main():
    all_rows = []
    for acc, sub in STUDIES.items():
        dest = os.path.join(DEST_ROOT, sub, "tables")
        os.makedirs(dest, exist_ok=True)
        items = file_list(acc)
        print("%s: %d tables -> %s" % (acc, len(items), dest), flush=True)
        done = 0
        with cf.ThreadPoolExecutor(max_workers=8) as ex:
            futs = {ex.submit(fetch, acc, dest, it): it for it in items}
            for f in cf.as_completed(futs):
                status, path, size, extra = f.result()
                done += 1
                if status == "fail":
                    print("FAIL %s %s" % (path, extra), flush=True)
                if done % 50 == 0:
                    print("  %s %d/%d" % (acc, done, len(items)), flush=True)
                all_rows.append((acc, os.path.basename(path), size, extra or "", status))
        ok = sum(1 for r in all_rows if r[0] == acc and r[4] in ("ok", "skip"))
        print("%s done: %d/%d present" % (acc, ok, len(items)), flush=True)

    man = os.path.join(DEST_ROOT, "download_manifest.tsv")
    with open(man, "w", encoding="utf-8", newline="") as fh:
        fh.write("study\tfile\tsize\tmd5\tstatus\n")
        for r in all_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print("manifest ->", man, flush=True)


if __name__ == "__main__":
    main()
    sys.exit(0)
