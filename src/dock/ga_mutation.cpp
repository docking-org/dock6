// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ga_mutation.cpp
//
// Implementation of the pure mutation-support helpers for DOCK_GA. Extracted verbatim
// (same push order, same weighting) from STEP 1 of GA_Recomb::mutation_selection in
// conf_gen_ga.cpp. See ga_mutation.h for the license header and scope rationale.
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ga_mutation.h"

namespace ga_mutation {

std::vector<int> build_weighted_type_pool( const TypeWeight * weights, int count )
{
    std::vector<int> pool;
    for ( int t = 0; t < count; t++ ) {
        if ( weights[t].enabled ) {
            for ( int i = 0; i < weights[t].coefficient; i++ ) {
                pool.push_back( weights[t].type_code );
            }
        }
    }
    return pool;
}

} // namespace ga_mutation
