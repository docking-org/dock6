# DOCK_GA end-to-end regression harness

This is the **behavior-preservation safety net** for the `conf_gen_ga` refactor: it proves
that the *science does not move*. It runs a real DOCK_GA case with a **fixed seed** and diffs
the outputs (poses, ensembles, scores, ordering, counts, log lines) against **golden**
outputs captured from the **unmodified** build. Same inputs + same seed ⇒ identical output.

> **Determinism target: exact, same-seed, bit-for-bit** (decision A default). Comparisons are
> exact unless you pass an explicit, justified tolerance. See "Tolerances" below.

## ⚠️ This must be run on your cluster, not in the refactor dev environment

Full end-to-end runs need things that are **not present in the refactor sandbox**:

- A complete **dock6 build**, which needs a **Fortran toolchain** (the `nab/`, `grid/`, and
  score modules) plus the usual C/C++ compilers. The sandbox has no Fortran compiler, so it
  cannot build or run `dock6`.
- A **real DOCK_GA input case**: a receptor grid, fragment libraries, a seeded `dock.in`, etc.
  The `tutorials/` shipped in this repo do **not** include a de-novo / GA case.

Therefore, per the refactor spec §5.2, **capturing and running the end-to-end goldens is the
user's pre-merge responsibility.** Everything here is built, self-tested (with a stub binary),
and ready — you supply the case inputs and the built binary. The refactor sandbox instead
relies on the fast [`compile_check.sh`](../compile_check.sh) gate and the unit suite
([`make test`](../README.md)) after every step, and on this harness on your cluster before merge.

## Case layout

```
cases/<name>/
    dock.in          # DOCK_GA input with a FIXED random_seed (and minimizer on — see note)
    inputs/          # everything dock.in references: grids, fraglibs, ligands, parameter files
    manifest.txt     # (optional) output basenames to compare, one per line; '#' comments ok
    golden/          # golden outputs — created by capture_golden.sh, committed to the repo
```

If `manifest.txt` is omitted, every file the run produces is captured/compared (including
`run.log`). A manifest is recommended so you compare exactly the scientific outputs and can
`--ignore-regex` volatile log lines (timestamps, wall-clock seconds).

> **Seeding note (from Phase 0):** RNG is seeded once via `srand(simplex.random_seed)` and a
> code comment indicates seeding effectively happens only when the minimizer is enabled. Make
> your regression `dock.in` **enable minimization** and set a fixed `random_seed`, or the run
> may not be reproducible by seed. Verify a case reproduces (run it twice, diff) before trusting
> its golden.

## Workflow

1. **Build the unmodified `dock6`** (pre-refactor) on your cluster.
2. **Capture goldens** for each case, from that unmodified binary:
   ```sh
   tests/regression/capture_golden.sh tests/regression/cases/mini /path/to/dock6
   ```
   Review `cases/mini/golden/`, then commit it.
3. **After every refactoring step**, rebuild and check:
   ```sh
   tests/regression/check_case.sh tests/regression/cases/mini /path/to/dock6
   ```
   Exit 0 = outputs match golden. Exit 1 = regression (diffs printed). This must stay green.

Run all cases, e.g.:
```sh
for c in tests/regression/cases/*/; do
    tests/regression/check_case.sh "$c" /path/to/dock6 || echo "FAILED: $c"
done
```

## Coverage the cases should exercise (spec §5.2.3)

> A ready-made suite with the exact `dock.in` key overrides for each surface is in
> [`docs/regression_case_matrix.md`](../../../../docs/regression_case_matrix.md).

Aim for at least one small, fast, fixed-seed case per changed surface:

- each **selection** method reachable today: elitism, tournament, roulette
- a representative of each **mutation** type: addition, deletion, substitution, replacement
- **crossover**: both `breeding_rand` (random sampling) and `breeding_exhaustive`
- each active **filter/cutoff**: MW (hard *and* soft), rotatable bonds, H-acceptors/donors, formal charge
- serial build (required); MPI build only needs to *compile/link* (this file has no MPI logic)

Keep them tiny (few generations, small ensemble) so the whole suite runs after every step.

## Tolerances

`compare.py` defaults to **exact**. Integers (atom/bond counts, indices, ordering) are **always**
compared exactly — a float tolerance can never mask a count or ordering change. Use tolerance
only where output *formatting* demands it, and record why:

```sh
# ignore wall-clock log lines; allow a tight epsilon on printed score floats
tests/regression/check_case.sh cases/mini /path/to/dock6 -- \
    --ignore-regex 'seconds' --ignore-regex 'Elapsed' --abs-tol 1e-6
```

## What the sandbox already verified

The plumbing (`run_case.sh` → `capture_golden.sh` → `check_case.sh` → `compare.py`) was proven
with a deterministic **stub** binary: capture→check passes when outputs are identical, and fails
(exit 1) when a score drifts by 3e-4 or an atom count changes by 1 — i.e. it catches exactly the
"silent numerical drift" the spec is most worried about. `compare.py` has its own token-level
checks. You only need to supply real cases and a real binary.
