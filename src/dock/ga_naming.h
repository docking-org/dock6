// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ga_naming.h
//
// Molecule-title formatting for DOCK_GA, extracted from the GA_Recomb god-class
// (conf_gen_ga.cpp) during the structural refactor. Pure string building over
// primitive inputs, so it is independently unit-testable and carries no dependency
// on DOCKMol or globals. Behavior is identical to GA_Recomb::naming_function.
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

#ifndef GA_NAMING_H
#define GA_NAMING_H

#include <string>

namespace ga_naming {

// Build a molecule title of the form "<prefix><identifier>_g<gen>_i<loc>", with the
// generation zero-padded to width 4 ("_g0005", "_g0042", "_g0300", "_g1234").
//
// NOTE (behavior preserved exactly, quirk included): the location field is padded
// LESS than the generation field. For loc < 10 it is "_i000<loc>" and for loc < 100
// it is "_i00<loc>", but for loc >= 100 it is just "_i<loc>" (no further padding) —
// i.e. loc is NOT consistently 4-wide the way gen is. This matches the original
// naming_function and must not be "fixed" here, as it would change output titles.
std::string build_molecule_title( const std::string & prefix,
                                  const std::string & identifier,
                                  int gen, int loc );

} // namespace ga_naming

#endif // GA_NAMING_H
