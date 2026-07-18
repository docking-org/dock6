// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// test_ga_filters.cpp
//
// Characterization tests for ga_filters — the candidate-filter predicates extracted
// from GA_Recomb (conf_gen_ga.cpp) in Phase 2: the MW cutoff (hard + soft) and the
// rotatable-bond / H-acceptor / H-donor / formal-charge checks.
//
// The soft MW path consumes RNG; these tests inject a deterministic stub generator so
// the probabilistic branch is exercised reproducibly, and assert the number and timing
// of draws (which must stay data-dependent exactly as in the original).
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "doctest.h"
#include "ga_filters.h"

#include <cmath>
#include <functional>

using namespace ga_filters;

// A stub "rand()" that returns a fixed value and counts how many times it is called.
struct StubRng {
    int value;
    int calls = 0;
    explicit StubRng(int v) : value(v) {}
    int operator()() { ++calls; return value; }
};

TEST_CASE("passes_mw_cutoff hard: accept inside [lower, upper], reject outside, no RNG") {
    StubRng rng(0);
    std::function<int()> gen = std::ref(rng);

    CHECK(passes_mw_cutoff(false, 300.0f, 100.0f, 500.0f, 35.0f, gen));   // inside -> pass
    CHECK_FALSE(passes_mw_cutoff(false, 600.0f, 100.0f, 500.0f, 35.0f, gen)); // above -> reject
    CHECK_FALSE(passes_mw_cutoff(false, 50.0f, 100.0f, 500.0f, 35.0f, gen));  // below -> reject
    CHECK(passes_mw_cutoff(false, 500.0f, 100.0f, 500.0f, 35.0f, gen));   // ==upper -> pass (not > )
    CHECK(passes_mw_cutoff(false, 100.0f, 100.0f, 500.0f, 35.0f, gen));   // ==lower -> pass (not < )

    CHECK(rng.calls == 0);  // hard path never draws
}

TEST_CASE("passes_mw_cutoff soft: no draw while inside bounds") {
    StubRng rng(50);
    std::function<int()> gen = std::ref(rng);

    CHECK(passes_mw_cutoff(true, 300.0f, 100.0f, 500.0f, 35.0f, gen));
    CHECK(rng.calls == 0);  // in-bounds -> no probabilistic draw
}

TEST_CASE("passes_mw_cutoff soft: one draw when over upper; decision matches exp(-Z^2)") {
    // mol_wt 535, upper 500, std_dev 35 -> excess 35, Z=1, acceptRate=exp(-1)=0.3679
    const double acceptRate = std::exp(-1.0);
    REQUIRE(acceptRate == doctest::Approx(0.367879));

    // rng=49 -> rand_num = 49%100+1 = 50 -> rand_num_dec = 0.50; 0.3679 < 0.50 -> REJECT
    {
        StubRng rng(49);
        std::function<int()> gen = std::ref(rng);
        CHECK_FALSE(passes_mw_cutoff(true, 535.0f, 100.0f, 500.0f, 35.0f, gen));
        CHECK(rng.calls == 1);  // exactly one draw, and only because upper was exceeded
    }
    // rng=19 -> rand_num = 20 -> rand_num_dec = 0.20; 0.3679 < 0.20 is false -> ACCEPT
    {
        StubRng rng(19);
        std::function<int()> gen = std::ref(rng);
        CHECK(passes_mw_cutoff(true, 535.0f, 100.0f, 500.0f, 35.0f, gen));
        CHECK(rng.calls == 1);
    }
}

TEST_CASE("passes_mw_cutoff soft: one draw when under lower; symmetric Z^2") {
    // mol_wt 65, lower 100, std_dev 35 -> excess 35, Z=1, acceptRate=exp(-1)=0.3679
    StubRng rng(49);  // rand_num_dec = 0.50 -> reject
    std::function<int()> gen = std::ref(rng);
    CHECK_FALSE(passes_mw_cutoff(true, 65.0f, 100.0f, 500.0f, 35.0f, gen));
    CHECK(rng.calls == 1);  // draw happens for the lower-bound branch only
}

TEST_CASE("passes_mw_cutoff soft: near-boundary excess almost always accepted") {
    // tiny excess -> acceptRate ~ 1, so accepted unless the draw is essentially 1.0
    StubRng rng(50);  // rand_num = 51 -> dec = 0.51
    std::function<int()> gen = std::ref(rng);
    CHECK(passes_mw_cutoff(true, 501.0f, 100.0f, 500.0f, 35.0f, gen)); // excess 1, acceptRate~0.999
    CHECK(rng.calls == 1);
}

TEST_CASE("exceeds_count_limit: strict greater-than (rot bonds / HA / HD)") {
    CHECK_FALSE(exceeds_count_limit(5u, 10));
    CHECK_FALSE(exceeds_count_limit(10u, 10));  // equal is allowed
    CHECK(exceeds_count_limit(11u, 10));
    CHECK_FALSE(exceeds_count_limit(0u, 0));
    CHECK(exceeds_count_limit(1u, 0));
}

TEST_CASE("outside_charge_range: symmetric [-limit, +limit], boundaries inclusive") {
    CHECK_FALSE(outside_charge_range(0.0f, 2.0f));
    CHECK_FALSE(outside_charge_range(2.0f, 2.0f));   // == +limit is inside
    CHECK_FALSE(outside_charge_range(-2.0f, 2.0f));  // == -limit is inside
    CHECK(outside_charge_range(2.5f, 2.0f));
    CHECK(outside_charge_range(-2.5f, 2.0f));
    // the real code uses ga_constraint_formal_charge = user + 0.1, e.g. 2.1
    CHECK_FALSE(outside_charge_range(2.0f, 2.1f));
    CHECK(outside_charge_range(2.2f, 2.1f));
}
