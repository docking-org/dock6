// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// test_ga_selection.cpp
//
// Tests for ga_selection::first_enabled — the pure decision behind the now
// table-driven selection dispatcher in GA_Recomb::selection_method. Pins the
// "first enabled wins, none => no selection" contract that replaced the original
// if/else-if ladder over the five selection-method booleans.
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "doctest.h"
#include "ga_selection.h"

using ga_selection::first_enabled;

TEST_CASE("first_enabled: picks the first true flag (ladder order preserved)") {
    // Order mirrors the registry: elitism, tournament, roulette, sus, metropolis.
    bool elitism_only[]    = {true,  false, false, false, false};
    bool tournament_only[] = {false, true,  false, false, false};
    bool roulette_only[]   = {false, false, true,  false, false};
    bool sus_only[]        = {false, false, false, true,  false};
    bool metropolis_only[] = {false, false, false, false, true };

    CHECK(first_enabled(elitism_only, 5)    == 0);
    CHECK(first_enabled(tournament_only, 5) == 1);
    CHECK(first_enabled(roulette_only, 5)   == 2);
    CHECK(first_enabled(sus_only, 5)        == 3);
    CHECK(first_enabled(metropolis_only, 5) == 4);
}

TEST_CASE("first_enabled: first match wins when several are set") {
    bool multi[] = {false, true, true, false, true};
    CHECK(first_enabled(multi, 5) == 1);  // like if/else-if: earliest true dispatches

    bool all[] = {true, true, true, true, true};
    CHECK(first_enabled(all, 5) == 0);
}

TEST_CASE("first_enabled: none set => -1 (dispatcher no-ops, as the original did)") {
    bool none[] = {false, false, false, false, false};
    CHECK(first_enabled(none, 5) == -1);
    CHECK(first_enabled(none, 0) == -1);   // empty range
}
