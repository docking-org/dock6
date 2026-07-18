// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ga_descriptors.cpp
//
// Implementation of the pure molecular-descriptor computations for DOCK_GA.
// Extracted verbatim (same values, same order, same accumulation types) from the
// GA_Recomb methods calc_mol_wt / calc_rot_bonds / num_HA_HD / calc_formal_charge
// in conf_gen_ga.cpp. See ga_descriptors.h for the license header.
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ga_descriptors.h"

#include <iostream>

namespace ga_descriptors {

// Atomic weights from General Chemistry 3rd Edition, Darrell D. Ebbing (per the
// original calc_mol_wt comment). Double literals preserved exactly so that the
// float accumulation in molecular_weight() is bit-for-bit identical to the original.
double atomic_weight( const std::string & atom, bool & recognized )
{
    recognized = true;

    if ( atom == "H" )
        return 1.00794;

    else if ( atom == "C.3" || atom == "C.2" || atom == "C.1" || atom == "C.ar" || atom == "C.cat" )
        return 12.011;

    else if ( atom == "N.4" || atom == "N.3" || atom == "N.2" || atom == "N.1" || atom == "N.ar" ||
              atom == "N.am" || atom == "N.pl3" )
        return 14.00674;

    else if ( atom == "O.3" || atom == "O.2" || atom == "O.co2" )
        return 15.9994;

    else if ( atom == "S.3" || atom == "S.2" || atom == "S.O" || atom == "S.o" || atom == "S.O2" ||
              atom == "S.o2" )
        return 32.066;

    else if ( atom == "P.3" )
        return 30.973762;

    else if ( atom == "F" )
        return 18.9984032;

    else if ( atom == "Cl" )
        return 35.4527;

    else if ( atom == "Br" )
        return 79.904;

    else if ( atom == "I" )
        return 126.90447;

    else if ( atom == "Du" )
        return 0;

    recognized = false;
    return 0;
}

float molecular_weight( const std::string * atom_types, int num_atoms )
{
    float mw = 0.0;

    for ( int i = 0; i < num_atoms; i++ ) {
        bool recognized;
        double w = atomic_weight( atom_types[i], recognized );
        if ( recognized ) {
            mw += w;   // float += double, matching the original `mw += <literal>`
        } else {
            // Preserved verbatim from the original (note: the legacy message names
            // DN_GA_Build, a pre-existing copy-paste artifact kept for output parity).
            std::cout << "WARNING: Did not recognize the atom_type " << atom_types[i]
                      << " in DN_GA_Build::calc_mol_wt()\n";
        }
    }

    return mw;
}

int count_rotatable_bonds( const int * amber_bt_id, int num_bonds )
{
    int counter = 0;
    for ( int i = 0; i < num_bonds; i++ ) {
        if ( amber_bt_id[i] != -1 ) {
            counter++;
        }
    }
    return counter;
}

int count_true_flags( const bool * flags, int n )
{
    int counter = 0;
    for ( int i = 0; i < n; i++ ) {
        if ( flags[i] == true ) {
            counter++;
        }
    }
    return counter;
}

float sum_charges( const float * charges, int num_atoms )
{
    float charge = 0.0;
    for ( int i = 0; i < num_atoms; i++ ) {
        charge += charges[i];
    }
    return charge;
}

// Covalent radii (Angstrom), generalized by element from the CRC handbook. Values are
// the COV_RADII_* constants from conf_gen_ga.h, kept identical here.
float covalent_radius( const std::string & atom, bool & recognized )
{
    recognized = true;

    if ( atom == "H" )
        return 0.32;   // COV_RADII_H

    else if ( atom == "C.3" || atom == "C.2" || atom == "C.1" || atom == "C.ar" || atom == "C.cat" )
        return 0.75;   // COV_RADII_C

    else if ( atom == "N.4" || atom == "N.3" || atom == "N.2" || atom == "N.1" || atom == "N.ar" ||
              atom == "N.am" || atom == "N.pl3" )
        return 0.71;   // COV_RADII_N

    else if ( atom == "O.3" || atom == "O.2" || atom == "O.co2" )
        return 0.64;   // COV_RADII_O

    else if ( atom == "S.3" || atom == "S.2" || atom == "S.O" || atom == "S.o" || atom == "S.O2" ||
              atom == "S.o2" )
        return 1.04;   // COV_RADII_S

    else if ( atom == "P.3" )
        return 1.09;   // COV_RADII_P

    else if ( atom == "F" )
        return 0.60;   // COV_RADII_F

    else if ( atom == "Cl" )
        return 1.00;   // COV_RADII_CL

    else if ( atom == "Br" )
        return 1.17;   // COV_RADII_BR

    recognized = false;
    return 0.71;   // original fallback for an unrecognized type
}

} // namespace ga_descriptors
