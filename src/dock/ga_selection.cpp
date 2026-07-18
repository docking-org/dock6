// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ga_selection.cpp
//
// Pure selection-seam helpers for DOCK_GA. See ga_selection.h for the license header
// and the rationale (table-driven selection dispatch).
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ga_selection.h"

namespace ga_selection {

int first_enabled( const bool * flags, int count )
{
    for ( int i = 0; i < count; i++ ) {
        if ( flags[i] ) {
            return i;
        }
    }
    return -1;
}

} // namespace ga_selection
