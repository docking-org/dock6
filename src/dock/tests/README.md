# DOCK_GA test suite (`src/dock/tests/`)

Test infrastructure for the `conf_gen_ga` refactor. It exists to hold two properties fixed
while the code's structure improves: **no silent numerical drift** and **preserved
determinism** (same seed ⇒ identical output). None of this is part of the normal DOCK build:
`make` / `make install` neither build nor require it, and it adds **no** dependency to DOCK.

## Three layers of safety net

| Layer | What it proves | Runs where | How to run |
|---|---|---|---|
| **Unit suite** (doctest) | Extracted pure logic behaves as characterized; new seams work incl. edge cases | anywhere with a C++ compiler | `make test` (or `cd tests && make run`) |
| **Compile-check gate** | `conf_gen_ga.cpp` still compiles after a refactor step | anywhere with a C++ compiler | `make check-compile` (or `tests/compile_check.sh`) |
| **End-to-end regression** | The *science* didn't move (poses/scores/ensembles/ordering) | your cluster (needs full build + real inputs) | `tests/regression/check_case.sh …` — see [regression/README.md](regression/README.md) |

The unit suite and compile-check run **in the refactor dev environment**; the end-to-end
harness must run **on your cluster** (a full `dock6` needs a Fortran toolchain absent here, and
no GA input case ships with the repo). See the regression README for why and how.

## Framework

- **doctest v2.4.11**, vendored as a single header ([`doctest.h`](doctest.h), MIT license) — nothing
  to install. This was decision **F** (default) in the Phase 0 doc.
- [`test_main.cpp`](test_main.cpp) is the single doctest entry point. Every other `test_*.cpp`
  just `#include "doctest.h"` and adds `TEST_CASE`s.
- The suite links only the standard library today ([`test_smoke.cpp`](test_smoke.cpp) proves the
  harness works and pins the RNG-determinism assumption). Real subject-code tests are added
  **as Phase 2 extracts low-dependency units** — see below.

## Running

```sh
# from src/dock (integrated: uses the compiler your ./configure selected)
make test            # build + run the unit suite
make check-compile   # fast syntax-only check of conf_gen_ga.cpp
make test-clean

# standalone (no ./configure needed — works in any checkout)
cd src/dock/tests
make run             # build + run unit suite
make check-compile   # compile gate
make clean
```

`make test` from `src/dock` requires a configured tree (`../../install/config.h`, produced by
`./configure`), like every other target in that Makefile. The **standalone** `cd tests && make run`
path deliberately needs no configure, so tests run even where DOCK can't be fully built.

## Adding tests as the refactor proceeds (Phase 2/3)

The rule from the spec: **every function extracted in Phase 2 gets characterization tests that
pin its current behavior before it moves**, and **every seam introduced in Phase 3 gets tests
for the interface plus each concrete implementation, including edge cases** (empty population,
boundary MW/charge, degenerate tournaments, fragment-incompatibility guards, …).

Mechanically, when Phase 2 extracts a low-dependency unit (e.g. `ga_filters.cpp`/`.h` with the
MW / rotatable-bond / HA-HD / charge predicates decoupled from `DOCKMol` and globals):

1. Add `tests/test_ga_filters.cpp` with `#include "doctest.h"` and the `TEST_CASE`s.
2. Add that file to `TEST_SRCS` in [`Makefile`](Makefile).
3. If the tests call into an extracted object file, add it to `UNIT_OBJS` in the Makefile (it's
   built by the normal DOCK build in `..`), and add the include/link wiring there.

RNG-dependent logic is tested deterministically by injecting a seeded/stub generator through
the centralized RNG interface introduced in Phase 2 (spec §6 step 5) — not by calling global
`rand()` in a test.

## Layout

```
tests/
    doctest.h            vendored framework (do not edit)
    test_main.cpp        doctest entry point (defines main once)
    test_smoke.cpp       framework + determinism smoke tests
    Makefile             standalone build of the unit suite + check-compile
    compile_check.sh     syntax-only regression gate for conf_gen_ga.cpp
    regression/          end-to-end golden-master harness (run on the cluster)
```
