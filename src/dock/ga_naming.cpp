// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ga_naming.cpp
//
// Implementation of DOCK_GA molecule-title formatting. Extracted verbatim (same
// padding branches, same output) from GA_Recomb::naming_function in conf_gen_ga.cpp.
// See ga_naming.h for the license header and the padding-quirk note.
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ga_naming.h"

#include <sstream>

namespace ga_naming {

std::string build_molecule_title( const std::string & prefix,
                                  const std::string & identifier,
                                  int gen, int loc )
{
    std::ostringstream new_title;

    new_title << prefix;

    // Generation field, zero-padded to width 4.
    if ( gen < 10 ) {
        new_title << identifier << "_g000" << gen;
    } else if ( gen < 100 ) {
        new_title << identifier << "_g00" << gen;
    } else if ( gen < 1000 ) {
        new_title << identifier << "_g0" << gen;
    } else {
        new_title << identifier << "_g" << gen;
    }

    // Location field — NOTE the padding is intentionally NOT symmetric with gen:
    // loc >= 100 gets no extra padding. Preserved exactly from the original.
    if ( loc < 10 ) {
        new_title << "_i000" << loc;
    } else if ( loc < 100 ) {
        new_title << "_i00" << loc;
    } else if ( loc < 1000 ) {
        new_title << "_i" << loc;
    } else {
        new_title << "_i" << loc;
    }

    return new_title.str();
}

} // namespace ga_naming
