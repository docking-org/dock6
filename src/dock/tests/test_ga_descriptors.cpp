// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// test_ga_descriptors.cpp
//
// Characterization tests for ga_descriptors — the pure molecular-descriptor
// computations extracted from GA_Recomb (conf_gen_ga.cpp) in Phase 2.
//
// These pin the CURRENT behavior of calc_mol_wt / calc_rot_bonds / num_HA_HD /
// calc_formal_charge so any future change that moves the numbers is caught here,
// at unit granularity, without needing a full DOCK build.
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "doctest.h"
#include "ga_descriptors.h"

#include <string>

using namespace ga_descriptors;

TEST_CASE("atomic_weight: recognized Sybyl types return exact literals") {
    bool ok = false;
    CHECK(atomic_weight("H", ok)   == doctest::Approx(1.00794));   CHECK(ok);
    CHECK(atomic_weight("C.3", ok) == doctest::Approx(12.011));    CHECK(ok);
    CHECK(atomic_weight("C.ar", ok)== doctest::Approx(12.011));    CHECK(ok);
    CHECK(atomic_weight("N.am", ok)== doctest::Approx(14.00674));  CHECK(ok);
    CHECK(atomic_weight("O.co2",ok)== doctest::Approx(15.9994));   CHECK(ok);
    CHECK(atomic_weight("S.o2", ok)== doctest::Approx(32.066));    CHECK(ok);
    CHECK(atomic_weight("P.3", ok) == doctest::Approx(30.973762)); CHECK(ok);
    CHECK(atomic_weight("F", ok)   == doctest::Approx(18.9984032));CHECK(ok);
    CHECK(atomic_weight("Cl", ok)  == doctest::Approx(35.4527));   CHECK(ok);
    CHECK(atomic_weight("Br", ok)  == doctest::Approx(79.904));    CHECK(ok);
    CHECK(atomic_weight("I", ok)   == doctest::Approx(126.90447)); CHECK(ok);
}

TEST_CASE("atomic_weight: dummy atom Du is recognized and weighs 0") {
    bool ok = false;
    CHECK(atomic_weight("Du", ok) == 0.0);
    CHECK(ok);  // Du is a RECOGNIZED type (no warning), distinct from unknown
}

TEST_CASE("atomic_weight: unknown types are flagged and contribute 0") {
    bool ok = true;
    CHECK(atomic_weight("Xx", ok) == 0.0);
    CHECK_FALSE(ok);
    ok = true;
    CHECK(atomic_weight("", ok) == 0.0);
    CHECK_FALSE(ok);
    ok = true;
    CHECK(atomic_weight("c.3", ok) == 0.0);  // case-sensitive: lowercase is unknown
    CHECK_FALSE(ok);
}

TEST_CASE("molecular_weight: sum over atoms, honoring num_atoms bound") {
    // methane-ish: C.3 + 4 H
    std::string atoms[] = {"C.3", "H", "H", "H", "H"};
    CHECK(molecular_weight(atoms, 5) == doctest::Approx(12.011 + 4 * 1.00794));

    // num_atoms bound is honored: only the first 1 atom counted
    CHECK(molecular_weight(atoms, 1) == doctest::Approx(12.011));

    // empty molecule
    CHECK(molecular_weight(atoms, 0) == 0.0f);

    // dummy atoms add nothing
    std::string withDu[] = {"C.3", "Du", "Du"};
    CHECK(molecular_weight(withDu, 3) == doctest::Approx(12.011));
}

TEST_CASE("molecular_weight: float+double accumulation matches the original") {
    // Reproduce the original accumulation exactly: float mw, double literals.
    std::string atoms[] = {"C.3", "N.am", "O.3", "S.3", "H", "H"};
    const int n = 6;

    float expected = 0.0f;   // same type as the original `float mw = 0.0;`
    expected += 12.011;      // C.3   (double literal, as in the source)
    expected += 14.00674;    // N.am
    expected += 15.9994;     // O.3
    expected += 32.066;      // S.3
    expected += 1.00794;     // H
    expected += 1.00794;     // H

    // Bit-for-bit equality (not Approx): the extraction must not change rounding.
    CHECK(molecular_weight(atoms, n) == expected);
}

TEST_CASE("count_rotatable_bonds: counts amber_bt_id != -1") {
    int ids[] = {-1, 0, 5, -1, 3};
    CHECK(count_rotatable_bonds(ids, 5) == 3);
    CHECK(count_rotatable_bonds(ids, 0) == 0);
    int none[] = {-1, -1, -1};
    CHECK(count_rotatable_bonds(none, 3) == 0);
    int all[] = {0, 1, 2};
    CHECK(count_rotatable_bonds(all, 3) == 3);
}

TEST_CASE("count_true_flags: counts true entries (HA/HD)") {
    bool f[] = {true, false, true, true, false};
    CHECK(count_true_flags(f, 5) == 3);
    CHECK(count_true_flags(f, 0) == 0);
    bool none[] = {false, false};
    CHECK(count_true_flags(none, 2) == 0);
}

TEST_CASE("sum_charges: float accumulation of partial charges") {
    float q[] = {0.5f, -0.25f, -0.25f, 1.0f};
    CHECK(sum_charges(q, 4) == doctest::Approx(1.0f));
    CHECK(sum_charges(q, 0) == 0.0f);

    // net-neutral molecule
    float neutral[] = {0.33f, -0.33f};
    CHECK(sum_charges(neutral, 2) == doctest::Approx(0.0f));
}

TEST_CASE("covalent_radius: recognized Sybyl types return the CRC-generalized radii") {
    bool ok = false;
    CHECK(covalent_radius("H", ok)   == doctest::Approx(0.32)); CHECK(ok);
    CHECK(covalent_radius("C.ar", ok)== doctest::Approx(0.75)); CHECK(ok);
    CHECK(covalent_radius("N.am", ok)== doctest::Approx(0.71)); CHECK(ok);
    CHECK(covalent_radius("O.co2",ok)== doctest::Approx(0.64)); CHECK(ok);
    CHECK(covalent_radius("S.o2", ok)== doctest::Approx(1.04)); CHECK(ok);
    CHECK(covalent_radius("P.3", ok) == doctest::Approx(1.09)); CHECK(ok);
    CHECK(covalent_radius("F", ok)   == doctest::Approx(0.60)); CHECK(ok);
    CHECK(covalent_radius("Cl", ok)  == doctest::Approx(1.00)); CHECK(ok);
    CHECK(covalent_radius("Br", ok)  == doctest::Approx(1.17)); CHECK(ok);
}

TEST_CASE("covalent_radius: unknown types warn-and-fallback to 0.71, flagged unrecognized") {
    bool ok = true;
    CHECK(covalent_radius("I", ok) == doctest::Approx(0.71));  // iodine has no entry here
    CHECK_FALSE(ok);
    ok = true;
    CHECK(covalent_radius("Xx", ok) == doctest::Approx(0.71));
    CHECK_FALSE(ok);
    // N atoms also happen to be 0.71, but they are RECOGNIZED (distinct from fallback)
    bool ok2 = false;
    CHECK(covalent_radius("N.3", ok2) == doctest::Approx(0.71));
    CHECK(ok2);
}
