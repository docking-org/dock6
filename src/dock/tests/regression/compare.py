#!/usr/bin/env python3
"""
compare.py — numeric-aware golden-master comparator for the DOCK_GA regression harness.

DOCK_GA's determinism target is same-seed bit-for-bit reproducibility, so by default
this compares files EXACTLY. The only reason it is "numeric-aware" is to let you apply
a *tight, explicit, justified* floating-point tolerance where (and only where) documented
output formatting genuinely requires it — never as a way to paper over real numerical
drift introduced by a refactor.

It works on any DOCK text output (mol2 poses/ensembles, .scores files, captured logs):
each line is split into whitespace tokens; tokens that parse as numbers are compared with
the chosen tolerance, all other tokens are compared exactly, and line/record ORDER must
match. That catches the things the refactor must preserve: molecule membership, ordering,
integer counts, and scores/coordinates.

Usage:
    compare.py GOLDEN CURRENT [--abs-tol A] [--rel-tol R] [--ignore-regex RE ...]
    compare.py --dir GOLDEN_DIR CURRENT_DIR [--manifest FILE] [...]

Exit status: 0 if equivalent within tolerance, 1 if differences found, 2 on usage error.

Examples:
    # strict exact match (default) of one ensemble file
    compare.py golden/restart0005.mol2 out/restart0005.mol2

    # allow a tight epsilon on printed floats (document WHY when you use this)
    compare.py --abs-tol 1e-6 golden/ga_output.scores out/ga_output.scores

    # compare a whole run directory against goldens, ignoring wall-clock log lines
    compare.py --dir golden/ out/ --ignore-regex 'seconds' --ignore-regex '^Elapsed'
"""

import argparse
import os
import re
import sys


_INT_RE = re.compile(r"^[+-]?\d+$")


def parse_number(tok):
    """Return float(tok) if tok is a number, else None. Handles ints, floats, sci-notation."""
    try:
        return float(tok)
    except (ValueError, TypeError):
        return None


def tokens_equal(a, b, abs_tol, rel_tol):
    """Compare two tokens: numeric-aware if both parse as numbers, else exact string match.

    Integers (atom/bond counts, molecule indices, ordering) are ALWAYS compared exactly,
    regardless of tolerance, per the determinism policy: tolerance is only ever for
    floating-point output formatting, never for counts or ordering.
    """
    if a == b:
        return True
    # Integer-looking tokens: exact only. A loose float tolerance must never mask a
    # changed count or index. If they differ as strings, they differ.
    if _INT_RE.match(a) and _INT_RE.match(b):
        return False
    na, nb = parse_number(a), parse_number(b)
    if na is None or nb is None:
        return False  # at least one is non-numeric and they differ as strings
    if na != na or nb != nb:  # NaN handling: NaN==NaN considered equal here
        return na != na and nb != nb
    diff = abs(na - nb)
    if diff <= abs_tol:
        return True
    scale = max(abs(na), abs(nb))
    return scale > 0 and diff <= rel_tol * scale


def line_matches_any(line, patterns):
    return any(p.search(line) for p in patterns)


def compare_lines(golden_lines, current_lines, abs_tol, rel_tol, ignore_patterns, label):
    """Return list of human-readable difference strings (empty == equivalent)."""
    diffs = []

    g = [ln for ln in golden_lines if not line_matches_any(ln, ignore_patterns)]
    c = [ln for ln in current_lines if not line_matches_any(ln, ignore_patterns)]

    if len(g) != len(c):
        diffs.append(
            "%s: line count differs (golden=%d, current=%d) after ignore filters"
            % (label, len(g), len(c))
        )

    for i, (gl, cl) in enumerate(zip(g, c), start=1):
        if gl == cl:
            continue
        gt, ct = gl.split(), cl.split()
        if len(gt) != len(ct):
            diffs.append("%s:%d token count differs\n  golden : %s\n  current: %s"
                         % (label, i, gl.rstrip(), cl.rstrip()))
            continue
        for j, (a, b) in enumerate(zip(gt, ct)):
            if not tokens_equal(a, b, abs_tol, rel_tol):
                diffs.append(
                    "%s:%d token %d differs: golden=%r current=%r\n  golden : %s\n  current: %s"
                    % (label, i, j + 1, a, b, gl.rstrip(), cl.rstrip()))
                break
    return diffs


def read_lines(path):
    with open(path, "r", errors="replace") as fh:
        return fh.readlines()


def compare_files(golden, current, abs_tol, rel_tol, ignore_patterns):
    if not os.path.isfile(golden):
        return ["MISSING golden file: %s" % golden]
    if not os.path.isfile(current):
        return ["MISSING current file: %s" % current]
    return compare_lines(read_lines(golden), read_lines(current),
                         abs_tol, rel_tol, ignore_patterns, os.path.basename(current))


def main(argv):
    ap = argparse.ArgumentParser(description="Numeric-aware golden-master comparator for DOCK_GA.")
    ap.add_argument("golden", help="golden file, or golden dir with --dir")
    ap.add_argument("current", help="current file, or current dir with --dir")
    ap.add_argument("--dir", action="store_true", help="treat args as directories")
    ap.add_argument("--manifest", help="file listing relative paths to compare (one per line); "
                                       "default with --dir is every file present in the golden dir")
    ap.add_argument("--abs-tol", type=float, default=0.0,
                    help="absolute float tolerance (default 0.0 = exact)")
    ap.add_argument("--rel-tol", type=float, default=0.0,
                    help="relative float tolerance (default 0.0 = exact)")
    ap.add_argument("--ignore-regex", action="append", default=[],
                    help="regex; lines matching it are ignored (repeatable). "
                         "Use for wall-clock/timestamp log lines.")
    args = ap.parse_args(argv)

    ignore_patterns = [re.compile(rx) for rx in args.ignore_regex]

    if args.abs_tol or args.rel_tol:
        sys.stderr.write(
            "NOTE: non-zero tolerance in use (abs=%g rel=%g). The determinism target is "
            "exact same-seed reproducibility; only relax tolerance where output formatting "
            "demonstrably requires it, and record why.\n" % (args.abs_tol, args.rel_tol))

    all_diffs = []
    if args.dir:
        if args.manifest:
            with open(args.manifest) as fh:
                rels = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
        else:
            rels = []
            for root, _dirs, files in os.walk(args.golden):
                for name in files:
                    full = os.path.join(root, name)
                    rels.append(os.path.relpath(full, args.golden))
            rels.sort()
        if not rels:
            sys.stderr.write("WARNING: no files to compare\n")
        for rel in rels:
            all_diffs += compare_files(os.path.join(args.golden, rel),
                                       os.path.join(args.current, rel),
                                       args.abs_tol, args.rel_tol, ignore_patterns)
    else:
        all_diffs = compare_files(args.golden, args.current,
                                  args.abs_tol, args.rel_tol, ignore_patterns)

    if all_diffs:
        print("REGRESSION: %d difference(s) found\n" % len(all_diffs))
        for d in all_diffs:
            print(d)
        return 1

    print("OK: outputs match within tolerance (abs=%g rel=%g)" % (args.abs_tol, args.rel_tol))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
