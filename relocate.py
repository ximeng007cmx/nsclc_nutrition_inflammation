#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
relocate.py — point the analysis scripts at *your* local checkout.

Background
----------
The scripts in this release were written for a fixed working directory. Every R script and several
Python tools contain hard-coded absolute paths, all of which are rooted at one project directory:

    D:/WorkBuddy_Y02/out_ali_directions/pathA

(The scripts are shipped unmodified so that the archived code is byte-identical to the code that
produced the published results.) Running this script once rewrites those paths to the location where
you cloned or extracted this repository.

Two replacements are performed
------------------------------
1. the GENCODE probe map, which in the original layout lived *next to* the project directory:

       D:/WorkBuddy_Y02/out_ali_directions/_gencode_probemap.gz
       ->  <root>/assets/_gencode_probemap.gz          (the copy shipped in this repository)

2. everything else under the project directory:

       D:/WorkBuddy_Y02/out_ali_directions/pathA[/...]
       ->  <root>[/...]

Both are applied in **a single left-to-right pass** over one alternation pattern, so a rewritten path
is never re-examined — cascading rewrites are impossible by construction.

The rules are anchored on the literal project path, not on a generic "drive letter" pattern. This is
deliberate: a broad pattern would also match the "s:" inside "https://" in the URLs that appear in
several scripts. Every other absolute-looking string in this repository is a URL and is left alone.

Usage
-----
    python relocate.py                    # target = this file's directory (recommended)
    python relocate.py /path/to/clone
    python relocate.py --dry-run          # show what would change, write nothing
    python relocate.py --restore          # revert to the original paths
    python relocate.py --list             # show the rules and exit

Only those paths are rewritten. Nothing else in any file is touched.
"""

import argparse
import re
import sys
from pathlib import Path

# ----------------------------------------------------------------------------- originals
ORIG_PARENT = "D:/WorkBuddy_Y02/out_ali_directions"          # parent directory
ORIG_ROOT = ORIG_PARENT + "/pathA"                            # the project directory
ORIG_PROBEMAP = ORIG_PARENT + "/_gencode_probemap.gz"         # lived beside the project directory

# Which files to process. Code only — the Markdown/JSON/CFF documentation deliberately keeps the
# original path as an illustrative example of what gets rewritten.
SUFFIXES = (".R", ".r", ".py", ".sh")

# Separator-agnostic matching. `_SEP` is a standalone one-character class; `_CLS_SEP` is the same
# two characters meant to be spliced *inside* another character class (splicing `_SEP` there would
# create a nested set, which Python's re rejects).
_SEP = r"[\\/]"
_CLS_SEP = r"\\/"

# What may legally follow a path literal: a separator, a closing quote, whitespace, a bracket, or EOL.
_BOUNDARY = r"(?=[" + _CLS_SEP + r"\"'\s,)\]}]|$)"

_ORIG_PARENT_RX = r"[A-Za-z]:" + _SEP + r"WorkBuddy_Y02" + _SEP + r"out_ali_directions"


def build_pattern(target, restoring):
    """Return (compiled_pattern, dispatch) for one direction.

    `dispatch(match) -> str` returns the replacement for whichever alternative matched, and raises
    for the impossible "no alternative matched" case.
    """
    if restoring:
        t = re.escape(target)
        rx = re.compile(
            r"(?P<asset>" + t + _SEP + r"assets" + _SEP + r"_gencode_probemap\.gz)"
            r"|(?P<root>" + t + r")" + _BOUNDARY
        )

        def dispatch(m):
            if m.group("asset") is not None:
                return ORIG_PROBEMAP
            return ORIG_ROOT

        return rx, dispatch

    rx = re.compile(
        r"(?P<probemap>" + _ORIG_PARENT_RX + _SEP + r"_gencode_probemap\.gz)"
        r"|(?P<root>" + _ORIG_PARENT_RX + _SEP + r"pathA)"
        r"|(?P<parent>" + _ORIG_PARENT_RX + r")" + _BOUNDARY
    )

    def dispatch(m):
        if m.group("probemap") is not None:
            return target + "/assets/_gencode_probemap.gz"
        # both `root` and `parent` map onto the checkout root
        return target

    return rx, dispatch


def iter_sources(here: Path):
    for p in sorted(here.rglob("*")):
        if p.is_file() and p.suffix in SUFFIXES and p.name != Path(__file__).name:
            yield p


def main() -> int:
    ap = argparse.ArgumentParser(description="Rewrite the hard-coded project paths.")
    ap.add_argument("target", nargs="?", default=None,
                    help="new project root (default: the directory containing this script)")
    ap.add_argument("--dry-run", action="store_true", help="show changes without writing")
    ap.add_argument("--restore", action="store_true", help="revert to the original paths")
    ap.add_argument("--list", action="store_true", help="print the rewrite rules and exit")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent

    if args.list:
        print("original parent   :", ORIG_PARENT)
        print("original project  :", ORIG_ROOT)
        print("original probemap :", ORIG_PROBEMAP)
        print("rule order        : probemap | root | parent   (single pass, no cascade)")
        return 0

    if args.restore:
        new_root = ORIG_ROOT
        target = str(here)
        if sys.platform == "win32":
            target = target.replace("\\", "/")
        print(f"restoring original paths (from {target})\n")
    else:
        new_root = str(Path(args.target).resolve()) if args.target else str(here)
        if sys.platform == "win32":
            new_root = new_root.replace("\\", "/")
        target = new_root
        if not Path(target).is_dir():
            print(f"ERROR: target directory does not exist: {target}", file=sys.stderr)
            return 2
        print(f"project root  ->  {target}\n")

    rx, dispatch = build_pattern(target, args.restore)

    changed_files = 0
    changed_hits = 0
    for src in iter_sources(here):
        # newline="" on BOTH read and write: no line-ending translation either way, so a file is
        # returned byte-identical after --restore. (Opening in the default text mode would silently
        # convert LF to CRLF on Windows and change the file's md5.)
        with open(src, "r", encoding="utf-8", errors="surrogateescape", newline="") as fh:
            text = fh.read()
        new_text, n = rx.subn(dispatch, text)
        if n:
            changed_files += 1
            changed_hits += n
            print(f"  [{n}] {src.relative_to(here)}")
            if not args.dry_run:
                with open(src, "w", encoding="utf-8", errors="surrogateescape", newline="") as fh:
                    fh.write(new_text)

    verb = "would rewrite" if args.dry_run else "rewrote"
    print(f"\n{verb} {changed_hits} path(s) across {changed_files} file(s).")
    if changed_files == 0:
        print("Nothing matched — the paths may already point elsewhere, or were edited.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
