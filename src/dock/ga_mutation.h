// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ga_mutation.h
//
// Pure mutation-operator support for DOCK_GA, extracted from the GA_Recomb god-class
// (conf_gen_ga.cpp) during the structural refactor.
//
// Only the genuinely pure, testable pieces of the mutation machinery live here. The
// operator dispatch in mutation_selection is intentionally NOT turned into a uniform
// registry: the four operators (deletion / addition / substitution / replacement) have
// materially different signatures, tagging, and pre/post processing (substitution is
// literally deletion-then-addition), so a common interface would be forced abstraction.
// It stays an explicit if/else chain; what moves out is the probability-weighting math.
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

#ifndef GA_MUTATION_H
#define GA_MUTATION_H

#include <vector>

namespace ga_mutation {

// One mutation type's weighting: whether it is enabled, how many times it enters the
// pool (its probability coefficient), and the type code pushed into the pool.
struct TypeWeight {
    bool enabled;
    int  coefficient;
    int  type_code;
};

// Build the weighted pool of mutation-type codes: for each enabled type, its type_code
// is appended `coefficient` times, in the order given. The caller then picks uniformly
// at random from the returned pool (pool[rand() % pool.size()]), so the relative
// frequency of a type equals its coefficient share. Mirrors STEP 1 of mutation_selection
// exactly (order: deletion, addition, substitution, replacement).
std::vector<int> build_weighted_type_pool( const TypeWeight * weights, int count );

} // namespace ga_mutation

#endif // GA_MUTATION_H
