// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// test_ga_mutation.cpp
//
// Characterization tests for ga_mutation::build_weighted_type_pool — the probability
// weighting behind STEP 1 of GA_Recomb::mutation_selection. Pins the push order and
// per-type multiplicity that determine how often each mutation operator is chosen.
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "doctest.h"
#include "ga_mutation.h"

#include <vector>

using ga_mutation::TypeWeight;
using ga_mutation::build_weighted_type_pool;

// Type codes mirror conf_gen_ga.h: DELETION 0, ADDITION 1, SUBSTITUTION 2, REPLACEMENT 3.
static const int DEL = 0, ADD = 1, SUB = 2, REP = 3;

TEST_CASE("weighted pool: each enabled type appears coefficient times, in order") {
    TypeWeight w[] = {
        {true, 5, DEL},
        {true, 1, ADD},
        {true, 1, SUB},
        {true, 5, REP},
    };
    std::vector<int> pool = build_weighted_type_pool(w, 4);

    // 5 + 1 + 1 + 5 = 12 entries, in deletion, addition, substitution, replacement order
    std::vector<int> expected = {DEL,DEL,DEL,DEL,DEL, ADD, SUB, REP,REP,REP,REP,REP};
    CHECK(pool == expected);
}

TEST_CASE("weighted pool: disabled types contribute nothing") {
    TypeWeight w[] = {
        {false, 5, DEL},
        {true,  3, ADD},
        {false, 9, SUB},
        {true,  2, REP},
    };
    std::vector<int> pool = build_weighted_type_pool(w, 4);
    std::vector<int> expected = {ADD, ADD, ADD, REP, REP};
    CHECK(pool == expected);
}

TEST_CASE("weighted pool: a zero coefficient adds nothing even if enabled") {
    TypeWeight w[] = {
        {true, 0, DEL},
        {true, 2, ADD},
    };
    std::vector<int> pool = build_weighted_type_pool(w, 2);
    std::vector<int> expected = {ADD, ADD};
    CHECK(pool == expected);
}

TEST_CASE("weighted pool: all disabled yields empty pool (caller then exits)") {
    TypeWeight w[] = {
        {false, 1, DEL},
        {false, 1, ADD},
        {false, 1, SUB},
        {false, 1, REP},
    };
    std::vector<int> pool = build_weighted_type_pool(w, 4);
    CHECK(pool.empty());  // mutation_selection treats this as "no types enabled" -> exit(1)
}

TEST_CASE("weighted pool: equal coefficients give uniform representation") {
    TypeWeight w[] = {
        {true, 1, DEL},
        {true, 1, ADD},
        {true, 1, SUB},
        {true, 1, REP},
    };
    std::vector<int> pool = build_weighted_type_pool(w, 4);
    std::vector<int> expected = {DEL, ADD, SUB, REP};
    CHECK(pool == expected);
}
