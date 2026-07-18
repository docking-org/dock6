// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ga_filters.cpp
//
// Implementation of the DOCK_GA candidate-filter predicates. Extracted verbatim
// (same comparisons, same float/double arithmetic, same RNG usage points) from the
// GA_Recomb methods mw_cutoff / hard_filter* in conf_gen_ga.cpp. See ga_filters.h
// for the license header and the RNG-injection rationale.
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ga_filters.h"

#include <cmath>

namespace ga_filters {

bool passes_mw_cutoff( bool soft, float mol_wt, float lower_bound, float upper_bound,
                       float std_dev, const std::function<int()> & rng )
{
    // Variable types mirror the original mw_cutoff exactly: rand_num / rand_num_dec
    // are float; excessMW / Z_scoreExcess / acceptRate are double.
    float  rand_num{};
    float  rand_num_dec{};
    double excessMW{};
    double Z_scoreExcess{};
    double acceptRate{};
    bool   result = true;

    if ( soft ) {
        // if exceeds upper boundary accept or reject with some probability
        if ( mol_wt > upper_bound ) {
            rand_num = ( rng() % 100 + 1 );            // random number from 1 to 100
            rand_num_dec = rand_num / 100;
            excessMW = mol_wt - upper_bound;            // how far above the cutoff
            Z_scoreExcess = excessMW / std_dev;         // ~ how many std devs over
            acceptRate = exp( -1 * Z_scoreExcess * Z_scoreExcess );  // Metropolis-like
            if ( acceptRate < rand_num_dec ) {
                result = false;
            }
        }

        // if exceeds lower boundary accept or reject with some probability
        if ( mol_wt < lower_bound ) {
            rand_num = ( rng() % 100 + 1 );
            rand_num_dec = rand_num / 100;
            excessMW = lower_bound - mol_wt;            // how far below the cutoff
            Z_scoreExcess = excessMW / std_dev;
            acceptRate = exp( -1 * Z_scoreExcess * Z_scoreExcess );
            if ( acceptRate < rand_num_dec ) {
                result = false;
            }
        }
    } else {
        if ( mol_wt > upper_bound ) {
            result = false;
        }
        if ( mol_wt < lower_bound ) {
            result = false;
        }
    }

    return result;
}

bool exceeds_count_limit( unsigned int value, int limit )
{
    return value > limit;
}

bool outside_charge_range( float formal_charge, float limit )
{
    return ( formal_charge > limit ) || ( formal_charge < -limit );
}

} // namespace ga_filters
