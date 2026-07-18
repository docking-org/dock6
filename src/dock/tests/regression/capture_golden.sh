#!/usr/bin/env bash
#
# capture_golden.sh — capture golden-master outputs for a DOCK_GA regression case
# from the CURRENT dock6 binary. Run this against the UNMODIFIED build BEFORE any
# refactoring, so the goldens represent pre-refactor scientific behavior.
#
# Usage:
#     capture_golden.sh <case_dir> <dock6_binary>
#
# Writes the collected output files into <case_dir>/golden/. If a manifest.txt is
# present in the case, only the listed files are stored as golden; otherwise every
# file produced by the run is stored (run.log included).

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "$#" -ne 2 ]; then
    echo "usage: capture_golden.sh <case_dir> <dock6_binary>" >&2
    exit 2
fi

CASE_DIR="$(cd "$1" && pwd)"
DOCK6="$2"
GOLDEN_DIR="$CASE_DIR/golden"

OUT_DIR="$("$HERE/run_case.sh" "$CASE_DIR" "$DOCK6" "$CASE_DIR/.capture" | tail -1)"

rm -rf "$GOLDEN_DIR"
mkdir -p "$GOLDEN_DIR"

if [ -f "$CASE_DIR/manifest.txt" ]; then
    while IFS= read -r rel; do
        [ -z "$rel" ] && continue
        case "$rel" in \#*) continue ;; esac
        if [ -f "$OUT_DIR/$rel" ]; then
            mkdir -p "$GOLDEN_DIR/$(dirname "$rel")"
            cp "$OUT_DIR/$rel" "$GOLDEN_DIR/$rel"
        else
            echo "WARNING: manifest lists '$rel' but the run did not produce it" >&2
        fi
    done < "$CASE_DIR/manifest.txt"
else
    # No manifest: everything the run produced except copied-in inputs is golden.
    cp -R "$OUT_DIR/." "$GOLDEN_DIR/"
fi

rm -rf "$CASE_DIR/.capture"
echo "Golden outputs written to $GOLDEN_DIR"
echo "Review them, then commit the golden/ directory so the harness can diff against it."
