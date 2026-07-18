#!/usr/bin/env bash
#
# check_case.sh — the regression check. Run a DOCK_GA case with the CURRENT
# (possibly refactored) dock6 and diff its outputs against the stored goldens.
# Run this after EVERY refactoring step; it must stay green.
#
# Usage:
#     check_case.sh <case_dir> <dock6_binary> [-- <extra compare.py args>]
#
# By default the comparison is EXACT (same-seed bit-for-bit). Pass tolerance /
# ignore flags after `--` only with an explicit, recorded justification, e.g.:
#     check_case.sh cases/mini ./dock6 -- --ignore-regex 'seconds' --abs-tol 1e-6
#
# Exit status: 0 if outputs match the goldens, 1 if they diverge, 2 on setup error.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "$#" -lt 2 ]; then
    echo "usage: check_case.sh <case_dir> <dock6_binary> [-- <compare.py args>]" >&2
    exit 2
fi

CASE_DIR="$(cd "$1" && pwd)"; shift
DOCK6="$1"; shift
COMPARE_ARGS=()
if [ "${1:-}" = "--" ]; then
    shift
    COMPARE_ARGS=("$@")
fi

GOLDEN_DIR="$CASE_DIR/golden"
[ -d "$GOLDEN_DIR" ] || {
    echo "ERROR: no golden/ in $CASE_DIR — run capture_golden.sh on the UNMODIFIED build first." >&2
    exit 2; }

OUT_DIR="$("$HERE/run_case.sh" "$CASE_DIR" "$DOCK6" "$CASE_DIR/.check" | tail -1)"

MANIFEST_ARG=()
[ -f "$CASE_DIR/manifest.txt" ] && MANIFEST_ARG=(--manifest "$CASE_DIR/manifest.txt")

set +e
python3 "$HERE/compare.py" --dir "$GOLDEN_DIR" "$OUT_DIR" \
    "${MANIFEST_ARG[@]}" "${COMPARE_ARGS[@]}"
rc=$?
set -e

rm -rf "$CASE_DIR/.check"

if [ "$rc" -eq 0 ]; then
    echo "PASS: case '$(basename "$CASE_DIR")' matches golden."
else
    echo "FAIL: case '$(basename "$CASE_DIR")' diverged from golden (see diffs above)."
fi
exit "$rc"
