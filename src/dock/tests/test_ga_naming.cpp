// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// test_ga_naming.cpp
//
// Characterization tests for ga_naming::build_molecule_title — the molecule-title
// formatter extracted from GA_Recomb::naming_function. Pins the exact zero-padding,
// including the deliberate gen-vs-loc asymmetry (a quirk preserved, not fixed).
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "doctest.h"
#include "ga_naming.h"

#include <string>

using ga_naming::build_molecule_title;

TEST_CASE("title: generation is zero-padded to width 4") {
    CHECK(build_molecule_title("", "ga", 0, 0)    == "ga_g0000_i0000");
    CHECK(build_molecule_title("", "ga", 5, 0)    == "ga_g0005_i0000");
    CHECK(build_molecule_title("", "ga", 42, 0)   == "ga_g0042_i0000");
    CHECK(build_molecule_title("", "ga", 300, 0)  == "ga_g0300_i0000");
    CHECK(build_molecule_title("", "ga", 1234, 0) == "ga_g1234_i0000");
}

TEST_CASE("title: prefix and identifier are prepended verbatim") {
    CHECK(build_molecule_title("run7.", "ga", 1, 2) == "run7.ga_g0001_i0002");
    CHECK(build_molecule_title("", "myid", 1, 2)    == "myid_g0001_i0002");
}

TEST_CASE("title: location padding is asymmetric with generation (preserved quirk)") {
    // loc < 10 and < 100 ARE padded...
    CHECK(build_molecule_title("", "ga", 0, 7)  == "ga_g0000_i0007");
    CHECK(build_molecule_title("", "ga", 0, 42) == "ga_g0000_i0042");
    // ...but loc >= 100 is NOT padded to width 4 (unlike gen). This is the quirk.
    CHECK(build_molecule_title("", "ga", 0, 100)  == "ga_g0000_i100");
    CHECK(build_molecule_title("", "ga", 0, 500)  == "ga_g0000_i500");
    CHECK(build_molecule_title("", "ga", 0, 1234) == "ga_g0000_i1234");
    // side-by-side: same numeric value, different width for g vs i
    CHECK(build_molecule_title("", "ga", 300, 300) == "ga_g0300_i300");
}
