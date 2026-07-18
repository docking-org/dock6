// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ga_filters.h
//
// Candidate-filter predicates for DOCK_GA, extracted from the GA_Recomb god-class
// (conf_gen_ga.cpp) during the structural refactor. These are the property cutoffs
// applied to evolved molecules: molecular weight (hard and the 6.13 "soft" cutoff),
// rotatable bonds, H-bond acceptors/donors, and formal-charge range.
//
// Each predicate is a free function over primitive descriptor values and constraint
// values, so it carries no dependency on DOCKMol or globals and is independently
// unit-testable. Behavior is bit-for-bit identical to the original inline checks,
// including exact float/double arithmetic. The soft molecular-weight path consumes a
// random draw; to preserve the exact RNG order/count of the original while keeping
// the logic deterministically testable, the generator is INJECTED (the caller passes
// a callback that returns values like C rand()); tests inject a stub.
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//
// This software is copyrighted, 2004-2026,
// by the DOCK Developers.
//
// The authors hereby grant permission to use, copy, modify, and re-distribute
// this software and its documentation for any purpose, provided
// that existing copyright notices are retained in all copies and that this
// notice is included verbatim in any distributions. No written agreement,
// license, or royalty fee is required for any of the authorized uses.
// Modifications to this software may be distributed provided that
// the nature of the modifications are clearly indicated.
//
// IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY
// FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
// ARISING OUT OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY
// DERIVATIVES THEREOF, EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES,
// INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE
// IS PROVIDED ON AN "AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE
// NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
// MODIFICATIONS.
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef GA_FILTERS_H
#define GA_FILTERS_H

#include <functional>

namespace ga_filters {

// Molecular-weight cutoff decision (extracted from GA_Recomb::mw_cutoff).
// Returns true if the molecule PASSES (is accepted), false if rejected.
//   soft == true  : outside a bound is rejected only probabilistically, with
//                   acceptRate = exp(-((excess/std_dev)^2)) compared against a draw
//                   rand()%100+1 scaled to (0,1]. Draws occur ONLY when a bound is
//                   exceeded, exactly as in the original (data-dependent RNG use).
//   soft == false : hard reject outside [lower, upper].
// The float/double types mirror the original exactly (mol_wt/bounds/std_dev are the
// float DOCKMol/constraint fields). `rng` must behave like C rand(): each call is one
// draw. Types match: value passed to rng()%100+1 is the raw rand() int.
bool passes_mw_cutoff( bool soft, float mol_wt, float lower_bound, float upper_bound,
                       float std_dev, const std::function<int()> & rng );

// True if an integer count strictly exceeds its upper limit (the "invalid" condition
// in the original filters). Used for rotatable bonds, H-bond acceptors, and donors.
// Argument types (unsigned value, signed limit) reproduce the original comparison's
// implicit conversion exactly.
bool exceeds_count_limit( unsigned int value, int limit );

// True if a formal charge is outside the symmetric range [-limit, +limit] (the
// "invalid" condition in the original charge filter).
bool outside_charge_range( float formal_charge, float limit );

} // namespace ga_filters

#endif // GA_FILTERS_H
