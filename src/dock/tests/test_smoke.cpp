// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// test_smoke.cpp
//
// Framework smoke test. Proves the vendored doctest harness compiles, links, and
// runs green in the serial configuration with no DOCK dependencies. It intentionally
// links against nothing but the standard library, so `make test` works even before
// any GA logic has been extracted into independently-linkable translation units.
//
// Real characterization / seam unit tests are added alongside Phase 2/3 as the
// god-functions in conf_gen_ga.cpp are decoupled into low-dependency units
// (ga_filters.*, ga_selection.*, ...). Each such unit gets its own tests/test_*.cpp.
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "doctest.h"

#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>

// --- Sanity: the framework itself works -------------------------------------
TEST_CASE("harness: doctest is wired up and runs") {
    CHECK(1 + 1 == 2);
    CHECK(std::string("dock") + "_ga" == "dock_ga");
}

// --- Determinism guardrail --------------------------------------------------
// DOCK_GA reproducibility rests on the C stdlib rand() stream seeded once via
// srand(simplex.random_seed) (conf_gen_ga.cpp). This documents and pins the
// assumption the whole refactor depends on: for a fixed seed, rand() produces a
// fixed sequence. If this ever fails, the platform's rand() is not the stable
// stream the golden-master harness assumes, and same-seed reproducibility across
// machines cannot be guaranteed by seed alone.
TEST_CASE("determinism: srand(seed) makes rand() reproducible within a run") {
    srand(12345);
    std::vector<int> a;
    for (int i = 0; i < 8; ++i) a.push_back(rand());

    srand(12345);
    std::vector<int> b;
    for (int i = 0; i < 8; ++i) b.push_back(rand());

    CHECK(a == b);  // same seed -> identical draw sequence

    srand(54321);
    std::vector<int> c;
    for (int i = 0; i < 8; ++i) c.push_back(rand());
    CHECK(a != c);  // different seed -> (overwhelmingly likely) different sequence
}
