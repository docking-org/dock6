// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ga_selection.h
//
// Selection-strategy support for DOCK_GA, introduced during the structural refactor
// to turn selection into an additive seam. Historically GA_Recomb::selection_method
// dispatched to one of five strategies through an if/else-if ladder over five parallel
// boolean members. That ladder is now table-driven (a registry of {label, enabling
// flag, handler} rows in conf_gen_ga.cpp), so a new strategy is added by registering a
// row rather than editing the dispatcher.
//
// This header holds the small, pure, dependency-free part of that seam so it can be
// unit-tested in isolation. The member-pointer registry itself lives with GA_Recomb
// (it necessarily depends on the class).
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

#ifndef GA_SELECTION_H
#define GA_SELECTION_H

namespace ga_selection {

// Index of the first enabled strategy in flags[0, count), or -1 if none is enabled.
// This reproduces the original if/else-if ladder's semantics exactly: first match
// wins, and if no flag is set nothing is selected (the dispatcher then no-ops, as
// the original did when none of the five booleans was true).
int first_enabled( const bool * flags, int count );

} // namespace ga_selection

#endif // GA_SELECTION_H
