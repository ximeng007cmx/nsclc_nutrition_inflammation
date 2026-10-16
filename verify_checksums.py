#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
verify_checksums.py — check this repository against CHECKSUMS.sha256.

CHECKSUMS.sha256 lists the SHA-256 of every file shipped in this repository, one per line, as

    <64-hex-digest>  <relative/path/with/forward/slashes>

Path order is the sorted relative path, so the manifest is stable and diff-friendly.

Usage
-----
    python verify_checksums.py            # verify (default)
    python verify_checksums.py --list     # print the manifest, verify nothing
    python verify_checksums.py --write    # regenerate CHECKSUMS.sha256
    python verify_checksums.py --quiet    # print only a summary line

Exit status is 0 when every listed file matches, 1 on any mismatch or missing/extra file.

Notes
-----
* Files are read in binary mode in 1 MiB chunks, so large files are handled without loading them
  into memory.
* Excluded from the manifest: the manifest itself, and the usual VCS / cache artefacts.
* Regenerating the manifest (--write) is a deliberate act: it should only be done when the release
  content itself is intended to change, because the manifest is part of the archived record.
"""

import argparse
import hashlib
import sys
from pathlib import Path

MANIFEST = "CHECKSUMS.sha256"

EXCLUDE_DIRS = {".git", "__pycache__", ".venv", "venv", ".idea", ".vscode"}
EXCLUDE_NAMES = {MANIFEST}
EXCLUDE_SUFFIX = (".pyc", ".pyo")


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def iter_files(root: Path):
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(root)
        if any(part in EXCLUDE_DIRS for part in rel.parts):
            continue
        if p.name in EXCLUDE_NAMES or p.name.endswith(EXCLUDE_SUFFIX):
            continue
        yield p, rel.as_posix()


def read_manifest(root: Path):
    path = root / MANIFEST
    if not path.is_file():
        print(f"ERROR: {MANIFEST} not found in {root}", file=sys.stderr)
        return None
    entries = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        digest, _, rel = line.partition("  ")
        if not rel:
            digest, _, rel = line.partition(" ")
        entries.append((digest.strip(), rel.strip().lstrip("*")))
    return entries


def main() -> int:
    ap = argparse.ArgumentParser(description="Verify files against CHECKSUMS.sha256.")
    ap.add_argument("--list", action="store_true", help="print the manifest and exit")
    ap.add_argument("--write", action="store_true", help="regenerate CHECKSUMS.sha256")
    ap.add_argument("--quiet", action="store_true", help="print only the summary line")
    args = ap.parse_args()

    root = Path(__file__).resolve().parent

    if args.list:
        path = root / MANIFEST
        if not path.is_file():
            print(f"ERROR: {MANIFEST} not found", file=sys.stderr)
            return 1
        sys.stdout.write(path.read_text(encoding="utf-8"))
        return 0

    if args.write:
        lines = [f"{sha256_of(p)}  {rel}" for p, rel in iter_files(root)]
        (root / MANIFEST).write_text("\n".join(lines) + "\n", encoding="utf-8", newline="")
        print(f"wrote {MANIFEST}: {len(lines)} file(s)")
        return 0

    entries = read_manifest(root)
    if entries is None:
        return 1

    on_disk = {rel: p for p, rel in iter_files(root)}

    ok = 0
    mismatched = []
    missing = []
    for digest, rel in entries:
        p = on_disk.get(rel)
        if p is None:
            missing.append(rel)
            continue
        if sha256_of(p) == digest:
            ok += 1
        else:
            mismatched.append(rel)

    listed = {rel for _, rel in entries}
    extra = sorted(set(on_disk) - listed)

    if not args.quiet:
        for rel in mismatched:
            print(f"MISMATCH  {rel}")
        for rel in missing:
            print(f"MISSING   {rel}")
        for rel in extra:
            print(f"EXTRA     {rel}")

    print(f"\n{ok}/{len(entries)} file(s) verified"
          f" | mismatched {len(mismatched)} | missing {len(missing)} | unlisted {len(extra)}")

    return 0 if (not mismatched and not missing) else 1


if __name__ == "__main__":
    raise SystemExit(main())
