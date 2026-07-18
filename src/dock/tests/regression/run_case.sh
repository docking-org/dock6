#!/usr/bin/env bash
#
# run_case.sh — execute one DOCK_GA regression case in an isolated directory and
# collect its output files. Shared by capture_golden.sh and check_case.sh.
#
# A "case" is a directory laid out as:
#     cases/<name>/
#         dock.in          # a DOCK_GA input file with a FIXED random_seed
#         inputs/          # everything dock.in references (grids, fraglibs, ligands, ...)
#         manifest.txt     # (optional) output basenames to collect/compare; if absent,
#                          # every file created during the run is collected
#         golden/          # golden outputs (created by capture_golden.sh; compared by check_case.sh)
#
# DOCK_GA writes its output files (restart####.mol2, *_pruned, *_filtered, *_rejected, ...)
# into the current working directory, so this script runs each case in a fresh scratch dir
# to keep runs hermetic and reproducible.
#
# Usage:
#     run_case.sh <case_dir> <dock6_binary> [out_dir]
#
# On success prints the output directory path on stdout (last line) and exits 0.

set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "usage: run_case.sh <case_dir> <dock6_binary> [out_dir]" >&2
    exit 2
fi

CASE_DIR="$(cd "$1" && pwd)"
DOCK6="$2"
OUT_DIR="${3:-$CASE_DIR/outputs}"

[ -f "$CASE_DIR/dock.in" ] || { echo "ERROR: no dock.in in $CASE_DIR" >&2; exit 2; }
command -v "$DOCK6" >/dev/null 2>&1 || [ -x "$DOCK6" ] || {
    echo "ERROR: dock6 binary not found/executable: $DOCK6" >&2; exit 2; }

# Fresh hermetic run directory.
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"
cp "$CASE_DIR/dock.in" "$OUT_DIR/"
if [ -d "$CASE_DIR/inputs" ]; then
    cp -R "$CASE_DIR/inputs/." "$OUT_DIR/"
fi

echo "=== running case '$(basename "$CASE_DIR")' with $DOCK6 ===" >&2
(
    cd "$OUT_DIR"
    # Capture the full log; scores/counts printed to stdout are part of the golden.
    "$DOCK6" -i dock.in > run.log 2>&1
) || { echo "ERROR: dock6 exited non-zero; see $OUT_DIR/run.log" >&2; exit 1; }

echo "=== output collected in $OUT_DIR ===" >&2
echo "$OUT_DIR"
