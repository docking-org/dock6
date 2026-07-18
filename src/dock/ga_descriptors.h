// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// ga_descriptors.h
//
// Pure molecular-descriptor computations for DOCK_GA, extracted from the
// GA_Recomb god-class (conf_gen_ga.cpp) during the structural refactor.
//
// These operate on the primitive fields of a DOCKMol (raw arrays + counts) rather
// than on DOCKMol itself, so they carry no dependency on DOCKMol, globals, or the
// scoring/typing machinery. That keeps them independently unit-testable. Behavior
// is bit-for-bit identical to the original in-class methods, including the exact
// floating-point accumulation order/type and the legacy warning text.
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

#ifndef GA_DESCRIPTORS_H
#define GA_DESCRIPTORS_H

#include <string>

namespace ga_descriptors {

// Atomic weight (g/mol) for a Sybyl atom type, returned as a double so that callers
// reproduce the original `float += <double literal>` accumulation exactly. Unknown
// types set recognized=false and return 0.0 (the original loop added nothing for them).
double atomic_weight( const std::string & atom_type, bool & recognized );

// Sum of atomic weights over atom_types[0, num_atoms). Accumulates into a float using
// double-precision per-term addition, matching GA_Recomb::calc_mol_wt exactly. For each
// unrecognized atom type it emits the original warning to std::cout, verbatim.
float molecular_weight( const std::string * atom_types, int num_atoms );

// Number of rotatable bonds: count of amber_bt_id[i] != -1 over [0, num_bonds).
int count_rotatable_bonds( const int * amber_bt_id, int num_bonds );

// Count of true entries over flags[0, n). Used for both H-bond acceptors and donors.
int count_true_flags( const bool * flags, int n );

// Sum of partial charges over charges[0, num_atoms), accumulated in float
// (matching GA_Recomb::calc_formal_charge). Molecule must be pre-charged (e.g. Gasteiger).
float sum_charges( const float * charges, int num_atoms );

// Covalent radius (Angstrom) for a Sybyl atom type (CRC handbook, generalized by
// element; extracted from GA_Recomb::calc_cov_radius). Unknown types set
// recognized=false and return 0.71 — the fallback the original used while warning.
float covalent_radius( const std::string & atom_type, bool & recognized );

} // namespace ga_descriptors

#endif // GA_DESCRIPTORS_H
