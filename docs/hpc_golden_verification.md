# Verifying the DOCK_GA refactor on an HPC cluster (golden-master check)

This is the step-by-step for confirming, **on your cluster**, that the `conf_gen_ga`
refactor is **bit-for-bit behavior-preserving**: same inputs + same seed ⇒ identical
output. It is the one verification that cannot be done in the refactor dev environment
(no Fortran toolchain there, and no de-novo/GA case ships in `tutorials/`).

The idea: build the **baseline** (pre-refactor) `dock6`, capture its outputs as *golden*
references, build the **refactored** `dock6`, and diff. The harness under
`src/dock/tests/regression/` does the running and diffing (see its
[README](../src/dock/tests/regression/README.md) for the case layout and comparator
details). This doc is the HPC-specific how-to around it.

---

## 0. Prerequisites on the cluster

You need a toolchain that can build the full `dock6` (unlike the dev sandbox):

- A C/C++ compiler and a **Fortran** compiler (for `nab/`, `grid/`, `resp/`, score
  modules). GCC (`gcc`/`g++`/`gfortran`) or Intel (`icc`/`icpc`/`ifort`) both work.
- `make`, plus `python3` (for the comparator) and `git`.
- Optional: MPI (`-DBUILD_DOCK_WITH_MPI` profiles) and/or RDKit + Boost (`*.rdkit`
  profiles) if you want to verify those build configurations too.

Load them via your module system, e.g.:

```sh
module load gcc/11        # provides gcc/g++/gfortran
module load python/3.11
# module load openmpi/4    # only if verifying an MPI build
# module load rdkit boost  # only if verifying an RDKit build
```

Pick the `install/` profile that matches (`gnu` for GCC serial, `intel` for Intel,
`gnu.rdkit`, `intel.intelmpi.parallel`, …). The commands below use `gnu`.

---

## 1. Get two builds: baseline (pre-refactor) and refactored

The determinism check compares outputs from the **unmodified** code against the
**refactored** code. Keep two checkouts so both binaries exist side by side.

```sh
# --- baseline: the pre-refactor tag (no ga_*.cpp units, original conf_gen_ga.cpp) ---
git clone https://github.com/docking-org/dock6 dock6-baseline
cd dock6-baseline
git checkout v6.13.1              # the version the refactor started from
cd install && ./configure gnu && cd ../src/dock
make -j                          # produces ../../bin/dock6 after `make install`; or dock6 here
cd ../../.. 

# --- refactored: your working branch with the ga_* units ---
git clone <your-refactor-remote-or-path> dock6-refactor
cd dock6-refactor
# (check out the refactor branch)
cd install && ./configure gnu && cd ../src/dock
make -j
cd ../../..
```

Note the two binaries, e.g. `dock6-baseline/src/dock/dock6` and
`dock6-refactor/src/dock/dock6`. If `make install` moved them, they are under
`.../bin/dock6`.

> **Sanity check the refactored build first:** in the refactor checkout, the unit suite
> and compile gate should already be green:
> ```sh
> cd dock6-refactor/src/dock && make test && make check-compile
> ```

---

## 2. Create at least one regression case

A case is a directory the harness runs (layout in the regression README). The critical
requirements for a *reproducible* case:

- a fixed `random_seed` in `dock.in`;
- **minimization enabled** — DOCK_GA seeds its RNG via `srand(simplex.random_seed)` and,
  per a code comment, effectively only when the minimizer is on; without it the run may
  not be seed-reproducible;
- small and fast (a few generations, small ensemble) so it reruns quickly.

```
dock6-refactor/src/dock/tests/regression/cases/mini/
    dock.in          # your seeded DOCK_GA config (minimizer on)
    inputs/          # grids, fragment libraries, starting ligand mol2, parameter files
    manifest.txt     # output basenames to compare (restart####.mol2, *.scores, ...)
```

**Confirm the case is actually reproducible before trusting it** — run it twice with the
baseline binary and diff; if two baseline runs differ, the case is not seed-deterministic
(fix the seed / enable minimization) and is not a valid golden:

```sh
cd dock6-refactor/src/dock/tests/regression
BIN=../../../../../dock6-baseline/src/dock/dock6
./run_case.sh cases/mini "$BIN" /tmp/mini_runA >/dev/null
./run_case.sh cases/mini "$BIN" /tmp/mini_runB >/dev/null
python3 compare.py --dir /tmp/mini_runA /tmp/mini_runB --ignore-regex 'seconds' --ignore-regex 'Elapsed'
# -> "OK: outputs match" means the case is deterministic and usable
```

Aim for a small suite covering each surface the refactor touched: each selection method
(elitism/tournament/roulette), a representative of each mutation type
(addition/deletion/substitution/replacement), crossover (`breeding_rand` and
`breeding_exhaustive`), and each active filter incl. MW **hard and soft**. A ready-made
suite with the exact `dock.in` overrides for each is in
[`regression_case_matrix.md`](regression_case_matrix.md).

---

## 3. Capture goldens from the BASELINE build

```sh
cd dock6-refactor/src/dock/tests/regression
BASELINE=../../../../../dock6-baseline/src/dock/dock6

for c in cases/*/; do
    [ -d "$c" ] || continue
    ./capture_golden.sh "$c" "$BASELINE"
done
# golden outputs are written into each cases/*/golden/ ; review, then keep them
```

These goldens are the source of truth for "the science as it was before the refactor."

---

## 4. Check the REFACTORED build against the goldens

```sh
cd dock6-refactor/src/dock/tests/regression
REFACTORED=../../src/dock/dock6        # the refactor checkout's own binary
# (adjust the path to wherever your refactored dock6 landed)

rc=0
for c in cases/*/; do
    [ -d "$c" ] || continue
    ./check_case.sh "$c" "$REFACTORED" -- --ignore-regex 'seconds' --ignore-regex 'Elapsed' \
        || { echo "FAIL: $c"; rc=1; }
done
[ "$rc" -eq 0 ] && echo "ALL CASES MATCH GOLDEN — refactor is bit-for-bit." || echo "REGRESSIONS ABOVE"
```

- Exit 0 / "OK" per case ⇒ the refactored build reproduces the baseline exactly.
- Exit 1 ⇒ a difference; the comparator prints the file, line, and token that diverged.

The comparison is **exact** by default. Integers (atom/bond counts, indices, ordering)
are always exact; the only tolerance in play is what you pass. Use `--ignore-regex` to
skip volatile log lines (wall-clock seconds, timestamps). Only add a float tolerance
(`--abs-tol 1e-6`) if a case has a documented formatting reason, and record why — the
determinism target is exact same-seed reproducibility.

---

## 5. Running it as a batch job (SLURM example)

```sh
#!/bin/bash
#SBATCH --job-name=dockga-golden
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=dockga-golden-%j.log

module load gcc/11 python/3.11

cd "$SLURM_SUBMIT_DIR/dock6-refactor/src/dock/tests/regression"
REFACTORED=../../src/dock/dock6

rc=0
for c in cases/*/; do
    [ -d "$c" ] || continue
    ./check_case.sh "$c" "$REFACTORED" -- --ignore-regex 'seconds' --ignore-regex 'Elapsed' \
        || { echo "FAIL: $c"; rc=1; }
done
exit $rc
```

Submit with `sbatch`. (Goldens in step 3 can be captured in an interactive session or a
similar batch job pointed at the baseline binary.) If you verify an MPI build, run the
binary under `srun`/`mpirun` inside `run_case.sh` — this file's logic has no MPI code, so
serial and MPI runs of the *GA* should produce the same molecules; MPI verification is
mainly "does the `-DBUILD_DOCK_WITH_MPI` build compile, link, and still reproduce."

---

## 6. If a case fails

1. **Read the printed diff** — it names the file, line number, and exact tokens. A changed
   *score* or *coordinate* points at floating-point drift; a changed *count/ordering* or
   *molecule name* points at a logic/RNG-order change.
2. **Localize with the unit suite** — `cd src/dock && make test`. The extracted units have
   characterization tests; a failure there pinpoints the unit.
3. **Bisect the refactor** — the changelog (architecture doc §9) lists each structural
   change in order; check out before/after a given change and re-run the failing case to
   isolate which one moved the result.
4. **Expected zero drift.** Every change so far is designed to be bit-for-bit; a diff means
   a real regression to fix, not a tolerance to widen. Do not relax tolerance to make a
   case pass.

---

## Quick reference

| Step | Command |
|---|---|
| Build baseline | `cd dock6-baseline/install && ./configure gnu && cd ../src/dock && make -j` |
| Build refactored | `cd dock6-refactor/install && ./configure gnu && cd ../src/dock && make -j` |
| Unit tests (refactored) | `cd dock6-refactor/src/dock && make test` |
| Capture goldens (baseline) | `tests/regression/capture_golden.sh cases/<c> <baseline-dock6>` |
| Check refactored vs golden | `tests/regression/check_case.sh cases/<c> <refactored-dock6> -- --ignore-regex 'seconds'` |
