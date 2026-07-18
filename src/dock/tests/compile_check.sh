#!/usr/bin/env bash
#
# compile_check.sh — cheap, dependency-light regression gate for the DOCK_GA refactor.
#
# Why this exists: a full dock6 build (and therefore an end-to-end golden run) is not
# feasible in every dev environment — it needs a Fortran toolchain (nab, score modules)
# and, for RDKit configs, RDKit + Boost. But the ONE file being refactored,
# conf_gen_ga.cpp, compiles cleanly on its own with a C++ compiler. This script
# syntax-checks it (-fsyntax-only, no link) so that after every Phase 2 refactoring step
# you get an immediate "does it still compile?" signal without the full build.
#
# It does NOT replace the end-to-end golden-master regression harness
# (tests/regression/) or the unit suite (make test) — it is a fast first line of defense.
#
# Usage:
#   tests/compile_check.sh            # serial config (always runs)
#   tests/compile_check.sh --rdkit    # also check the RDKit config (needs $RDBASE, $BOOST)
#
# Exit status: 0 if all attempted configs compile, non-zero otherwise.

set -u

# Resolve paths relative to this script so it works from any CWD.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCK_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SRC="$DOCK_DIR/conf_gen_ga.cpp"

CXX="${CXX:-c++}"
STD="${STD:--std=c++11}"
INCLUDES=(-I"$DOCK_DIR" -I"$DOCK_DIR/nab")

want_rdkit=0
for arg in "$@"; do
    case "$arg" in
        --rdkit) want_rdkit=1 ;;
        *) echo "unknown arg: $arg" >&2; exit 2 ;;
    esac
done

if [ ! -f "$SRC" ]; then
    echo "ERROR: cannot find $SRC" >&2
    exit 2
fi

fail=0

run_check() {
    local label="$1"; shift
    echo "=== compile-check: $label ==="
    local log
    log="$(mktemp)"
    if "$CXX" $STD -fsyntax-only "${INCLUDES[@]}" "$@" "$SRC" 2>"$log"; then
        local warns
        warns="$(grep -c 'warning:' "$log" 2>/dev/null || echo 0)"
        echo "  PASS (0 errors, ${warns} warnings)"
    else
        echo "  FAIL — errors:"
        grep 'error:' "$log" | sed 's/^/    /' | head -40
        fail=1
    fi
    rm -f "$log"
}

# --- Serial configuration (always) ------------------------------------------
run_check "serial (no RDKit, no MPI)"

# --- RDKit configuration (optional) -----------------------------------------
# conf_gen_ga.cpp has no MPI code, so there is no separate MPI syntax config to
# check here — MPI only affects sibling files. RDKit, however, pulls rdtyper.h /
# RDKit headers, so only attempt it when an RDKit install is discoverable.
if [ "$want_rdkit" -eq 1 ]; then
    if [ -n "${RDBASE:-}" ] && [ -n "${BOOST:-}" ]; then
        run_check "rdkit (-DBUILD_DOCK_WITH_RDKIT)" \
            -DBUILD_DOCK_WITH_RDKIT -I"$RDBASE/Code" -I"$BOOST/include"
    else
        echo "=== compile-check: rdkit — SKIPPED (set \$RDBASE and \$BOOST to enable) ==="
    fi
fi

echo
if [ "$fail" -eq 0 ]; then
    echo "compile-check: ALL PASSED"
else
    echo "compile-check: FAILURES ABOVE"
fi
exit "$fail"
