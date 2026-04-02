/* Parts of this code use RDKit and we are obliged to include the 
 * COPYRIGHT NOTICE below
 
Copyright (c) 2006-2015, Rational Discovery LLC, Greg Landrum, and Julie Penzotti and others
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <iostream>
#include <sstream>
#include <string.h>
#include "dockmol.h"

#ifdef BUILD_DOCK_WITH_RDKIT
#include <fstream>
// RDKit stuff
#include <boost/range/iterator_range.hpp>
#include <Geometry/point.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Conformer.h>

// Important definitions
// Necessary for RDKit::RWMol object
#define ATOMIC_NUMBER_H 1
#define ATOMIC_NUMBER_C 6
#define ATOMIC_NUMBER_N 7
#define ATOMIC_NUMBER_O 8
#define ATOMIC_NUMBER_P 15
#define ATOMIC_NUMBER_S 16
#define ATOMIC_NUMBER_F 9
#define ATOMIC_NUMBER_CL 17
#define ATOMIC_NUMBER_BR 35
#define ATOMIC_NUMBER_I 53
#endif

using namespace std;

#ifdef BUILD_DOCK_WITH_RDKIT
// GDRM Feb 19, 2020
int find_atomic_number( std::string atom_type ){

    // This function takes the SYBYL atom type (or its lower/upper case variant)
    // and returns the atomic number (Z), i.e., the number of protons in a nucleus.

    if ( atom_type == "H" || atom_type == "Du" ){
        return ATOMIC_NUMBER_H;
    } else if ( atom_type == "C.3" || atom_type == "C.2" || atom_type == "C.1" || atom_type == "C.ar" || atom_type == "C.cat" ){
        return ATOMIC_NUMBER_C;
    } else if ( atom_type == "N.4" || atom_type == "N.3" || atom_type == "N.2" || atom_type == "N.1" || atom_type == "N.ar" ||
                atom_type == "N.am" || atom_type == "N.pl3" ){
        return ATOMIC_NUMBER_N;
    } else if ( atom_type == "O.3" || atom_type == "O.2" || atom_type == "O.co2" || atom_type == "O.CO2" || atom_type == "O.Co2" ){
        return ATOMIC_NUMBER_O;
    } else if ( atom_type == "S.3" || atom_type == "S.2" || atom_type == "S.O" || atom_type == "S.o" || atom_type == "S.O2" ||
                atom_type == "S.o2" ){
        return ATOMIC_NUMBER_S;
    } else if ( atom_type == "P.3" ){
        return ATOMIC_NUMBER_P;
    } else if ( atom_type == "F" ){
        return ATOMIC_NUMBER_F;
    } else if ( atom_type == "Cl" || atom_type == "CL" ){
        return ATOMIC_NUMBER_CL;
    } else if ( atom_type == "Br" || atom_type == "BR" ){
        return ATOMIC_NUMBER_BR;
    } else if ( atom_type == "I") {
        return ATOMIC_NUMBER_I;
    } else {
        cout << "======================================================================" << endl;
        cout << "WARNING: Did not recognize " << atom_type << " in find_atomic_number()" << endl;
        cout << "======================================================================" << endl;
        return ATOMIC_NUMBER_N; //GA_Recomb::calc_cov_radius() returns N's radius.
    }
}
// GDRM Feb 19, 2020
RDKit::Bond::BondType select_bond_type( std::string bond_order ){

    // DOCK6 stores bond order information as strings. RDKit has its own BondType object
    // and requires conversion to a numeric type (double or int)

    RDKit::Bond::BondType bond_type;

    if ( bond_order == "1" || bond_order == " 1" || bond_order == "1 " || bond_order == "am" ){
        bond_type = RDKit::Bond::SINGLE;
    } else if ( bond_order == "2" || bond_order == " 2" || bond_order == "2 " ){
        bond_type = RDKit::Bond::DOUBLE;
    } else if ( bond_order == "3" || bond_order == " 3" || bond_order == "3 " ){
        bond_type = RDKit::Bond::TRIPLE;
    } else if ( bond_order == "ar" ){
        bond_type = RDKit::Bond::AROMATIC;
    } else {
        cout << "======================================================================" << endl;
        cout << "    WARNING: DOCK6 does not recognize a bond order of " << bond_order << endl;
        cout << "======================================================================" << endl;
    }
    return bond_type;
}
#endif

/***********************************************************************/
bool
Read_Mol2(DOCKMol & mol, istream & ifs, bool read_color, bool read_solvation,
          bool read_amber)
{
    char            line[1000];
    int             count;
    int             i;
    int             n1;
    bool            atom_line;
    char            typ[100],
                    col[100];

    char            tmp1[100],
                    tmp2[100],
                    jwu_tmp4[100],
		            jwu_tmp3[100],
                    subst_name[100];
    int             natoms,
                    nbonds,
                    nresidues; //defined by jwu
    string          l1,
                    l2,
                    l3,
                    l4,
                    l5,
                    l6;
    float           f1,
                    f2,
                    f3,
                    f4,
                    f5;

    bool            found_solvation = false;
    bool            found_color = false;
    bool            found_amber = false;

    bool            check_cols = false;
    bool            seven_columns = false; // default @<TRIPOS>ATOM line has 9 columns. By GDRM
    //int             bad_line{0}; // number of bad atom lines. By GDRM

    // init state vars
    atom_line = false;

    // read forward until the tripos molecule tag is reached
    for (;;) {
        if (!ifs.getline(line, 1000)) {
            mol.clear_molecule();
            // cout << endl << "ERROR: Ligand file empty.  Program will
            // terminate." << endl;
            return false;
        }

        if (!strncmp(line, "@<TRIPOS>MOLECULE", 17))
            break;
    }

    // loop over the header info
    for (count = 0;; count++) {

        if (!ifs.getline(line, 1000)) {
            mol.clear_molecule();
            cout << endl <<
                "ERROR:  Ligand file empty.  Program will terminate." << endl;
            return false;
        }

        if (!strncmp(line, "@<TRIPOS>ATOM", 13)) {
            atom_line = true;
            break;
        }
        // assign the first 5 header lines to the proper fields
        switch (count) {

        case 0:
            l1 = line;
            break;

        case 1:
            l2 = line;
            break;

        case 2:
            l3 = line;
            break;

        case 3:
            l4 = line;
            break;

        case 4:
            l5 = line;
            break;

        case 5:
            l6 = line;
            break;
        }

    }

    // if there are no atoms, throw error and return false
    if (!atom_line) {
        mol.clear_molecule();
        cout <<
            "ERROR: @<TRIPOS>ATOM indicator missing from ligand file.  Program will terminate."
            << endl;
        return false;
    }
    // get # of atoms and bonds from mol info line
    // sscanf(l2.c_str(), "%d %d", &natoms, &nbonds);
    sscanf(l2.c_str(), "%d %d %d", &natoms, &nbonds, &nresidues);

    // If nresidues is zero, DOCK6 needs to be told to read the @<TRIPOS>ATOM lines
    // correctly. By GDRM
    if (nresidues == 0){
        check_cols = true;
        nresidues = 1;
    }

    // initialize molecule vectors
    // mol.allocate_arrays(natoms, nbonds);
    mol.allocate_arrays(natoms, nbonds, nresidues);

    mol.title = l1;
    mol.mol_info_line = l2;
    mol.comment1 = l3;
    mol.comment2 = l4;
    mol.energy = l5;
    mol.comment3 = l6;

    // loop over atoms and read in atom info
    for (i = 0; i < mol.num_atoms; i++) {
        if (!ifs.getline(line, 1000)) {
            mol.clear_molecule();
            cout <<
                "ERROR:  Atom information missing from ligand file.  Program will terminate."
                << endl;
            return false;
        }
        //jwu comment
        // 1     N       15.2586  -59.3416   35.3528 N.4     1 PRO1  -0.2020
        // GDRM comment:
        // when nres = 0, @<TRIPOS>ATOM might contain 7 columns instead of 9:
        // 1     N       15.2586  -59.3416   35.3528 N.4  -0.2020
        if (check_cols){
            // Test number of columns per line
            stringstream string_test;
            string_test << line;
            int countCols{0};
            double value{};
            while(string_test >> value){
                ++countCols;
            }
            if (countCols == 7){
                seven_columns == true;
            }
        }
        //  Read data from mol2 as a function of the number of columns
        if (seven_columns){
            sscanf(line, "%s %s %f %f %f %s %f", jwu_tmp4, tmp1, &mol.x[i], &mol.y[i],
                   &mol.z[i], tmp2, &mol.charges[i]);
            mol.atom_names[i] = tmp1;
            mol.atom_types[i] = tmp2;
            mol.atom_number[i] = jwu_tmp4;
            mol.atom_residue_numbers[i] = "1";
            mol.subst_names[i] = subst_name;//mol.title
        } else {
            sscanf(line, "%s %s %f %f %f %s %s %s %f", jwu_tmp4, tmp1, &mol.x[i], &mol.y[i],
                   &mol.z[i], tmp2, jwu_tmp3, subst_name, &mol.charges[i]);
            mol.atom_names[i] = tmp1;
            mol.atom_types[i] = tmp2;
	        mol.atom_number[i] = jwu_tmp4;
	        mol.atom_residue_numbers[i] = jwu_tmp3;
            mol.subst_names[i] = subst_name;//mol.title
        }
    }
    // If all partial charges are zero, it is a bad molecule
    int num_zeros = 0;
    for (int i = 0; i < mol.num_atoms; ++i){
        if (mol.charges[i] == 0.0){
            num_zeros += 1;
        }
    }
    // Is it a bad molecule?
    if (num_zeros == mol.num_atoms){
        mol.bad_molecule = true;
    } else {
        mol.bad_molecule = false;
    }

    // skip down to the bond section
    for (;;) {
        if (!ifs.getline(line, 1000)) {
            mol.clear_molecule();
            cout <<
                "ERROR: @<TRIPOS>BOND indicator missing from ligand file.  Program will terminate."
                << endl;
            return false;
        }

        if (!strncmp(line, "@<TRIPOS>BOND", 13))
            break;
    }

    // loop over bonds and add them
    for (i = 0; i < mol.num_bonds; i++) {
        if (!ifs.getline(line, 1000)) {
            mol.clear_molecule();
            cout <<
                "ERROR: Bond information missing from ligand file.  Program will terminate."
                << endl;
            return false;
        }


        sscanf(line, "%*d %d %d %s", &mol.bonds_origin_atom[i],
               &mol.bonds_target_atom[i], tmp1);

        // adjust bond atom #'s to start at 0
        mol.bonds_origin_atom[i]--;
        mol.bonds_target_atom[i]--;

        mol.bond_types[i] = tmp1;

    }


    // ID ring atoms/bonds
    mol.id_ring_atoms_bonds();


    // Read Atom Color from a pre-colored mol2 file kxr 0206
    // skip to the color section

    if (read_color) {
        for (;;) {

            if (!ifs.getline(line, 1000)) {
                mol.clear_molecule();
                cout <<
                    "ERROR: @<TRIPOS>COLOR indicator missing from ligand file.  Program will terminate."
                    << endl;
                return false;
            }

            if (!strncmp(line, "@<TRIPOS>COLOR", 14)) {
                found_color = true;
                break;
            }


        }


        if (found_color) {

            for (i = 0; i < mol.num_atoms; i++) {

                if (!ifs.getline(line, 1000)) {
                    mol.clear_molecule();
                    cout <<
                        "ERROR: Coloring information missing from ligand file.  Program will terminate."
                        << endl;
                    return false;
                }

                sscanf(line, "%d %99s %99s", &n1, typ, col);
                mol.atom_color[i] = col;
            }
        }
    }
    // Read Total and atomic desolvation numbers kxr
    // skip to the solvation info section
    if (read_solvation) {

        for (;;) {
            if (!ifs.getline(line, 1000)) {
                mol.clear_molecule();
                cout << endl <<
                    "ERROR: @<TRIPOS>SOLVATION indicator missing from ligand file.  Program will terminate."
                    << endl;
                return false;
            }
            // if(!found_solvation) {
            if (!strncmp(line, "@<TRIPOS>SOLVATION", 18)) {
                found_solvation = true;
                break;
            }                   // else { cout << endl << "ERROR: Desolvation
                                // input absent.  Program will terminate." <<
                                // endl; return false; } } else break;
        }

        if (found_solvation) {

            if (!ifs.getline(line, 1000)) {
                mol.clear_molecule();
                cout <<
                    "ERROR: Solvation information missing from ligand file.  Program will terminate."
                    << endl;
                return false;
            }

            mol.total_dsol = 0.0f;
            int             start_atm = 0;

            if (ifs.getline(line, 1000)) {
                sscanf(line, "%f %f %f %f %f", &f1, &f2, &f3, &f4, &f5);
                mol.total_dsol += f5;
                mol.atom_psol[start_atm] = f2;
                mol.atom_apsol[start_atm] = f4;
                start_atm = 1;
            }

            for (i = start_atm; i < mol.num_atoms; i++) {

                if (ifs.getline(line, 1000)) {
                    sscanf(line, "%f %f %f %f %f", &f1, &f2, &f3, &f4, &f5);
                    mol.total_dsol += f5;
                    mol.atom_psol[i] = f2;
                    mol.atom_apsol[i] = f4;
                }
            }
        }
    }

    if (read_amber) {
        for (;;) {
            // check if AMBER_SCORE_ID is even present
            if (!ifs.getline(line, 1000)) {
                mol.clear_molecule();
                cout << endl <<
                    "ERROR: @<TRIPOS>AMBER_SCORE_ID indicator missing from ligand file.  Program will terminate."
                    << endl;
                return false;
            }
            // identify that we are in the amber_score_id portion of the file
            if (!strncmp(line, "@<TRIPOS>AMBER_SCORE_ID", 23)) {
                found_amber = true;
                break;
            }
        }

        if (found_amber) {
            // check if there is any information in the AMBER_SCORE_ID field
            if (!ifs.getline(line, 1000)) {
                mol.clear_molecule();
                cout <<
                    "ERROR: No AMBER Score location information available.  Program will terminate."
                    << endl;
                return false;
            }

            sscanf(line, "%99s", tmp1);
            mol.amber_score_ligand_id = tmp1;
        }
    }


    return true;
}


/***********************************************************************/
bool Write_Mol2(DOCKMol & mol, ostream & ofs)
{
    char            line[1000];
    int             i;
    vector < int   >renumber;
    int             current_atom,
                    current_bond;

    // init atom/bond renumbering data
    renumber.resize(mol.num_atoms + 1);
    current_atom = 1;
    current_bond = 1;

    // write out header information
    ofs << "@<TRIPOS>MOLECULE" << endl;

    if (mol.title.empty())
        ofs << "*****" << endl;
    else
        ofs << mol.title << endl;

    sprintf(line, " %d %d 1 0 0", mol.num_active_atoms, mol.num_active_bonds);
    ofs << line << endl;

    ofs << mol.comment1 << endl;
    ofs << mol.comment2 << endl;
    ofs << mol.energy << endl;
    ofs << mol.comment3 << endl;

    // write out atom lines
    ofs << "@<TRIPOS>ATOM" << endl;

    for (i = 0; i < mol.num_atoms; i++) {
        if (mol.atom_active_flags[i]) {

            sprintf(line,
                    "%7d%1s%-6s%12.4f%10.4f%10.4f%1s%-5s%4s%1s %-8s%10.4f",
                    current_atom, "", mol.atom_names[i].c_str(), mol.x[i],
                    mol.y[i], mol.z[i], "", mol.atom_types[i].c_str(),
                    mol.atom_residue_numbers[i].c_str(), "",
                    mol.subst_names[i].c_str(), mol.charges[i]);

            ofs << line << endl;

            renumber[i] = current_atom;
            current_atom++;
        }
    }

    // write out bond lines
    ofs << "@<TRIPOS>BOND" << endl;

    for (i = 0; i < mol.num_bonds; i++) {
        if (mol.bond_active_flags[i]) {

            sprintf(line, "%6d%6d%6d%3s%2s", current_bond,
                    renumber[mol.bonds_origin_atom[i]],
                    renumber[mol.bonds_target_atom[i]], "",
                    mol.bond_types[i].c_str());

            ofs << line << endl;
            current_bond++;
        }
    }

    ofs << "@<TRIPOS>SUBSTRUCTURE" << endl;
    sprintf(line,
            "     1 %-4s        1 TEMP              0 ****  ****    0 ROOT",
            mol.subst_names[0].c_str());
    ofs << line << endl << endl;


    return true;
}

/***********************************************************************/
//bool
//DOCKMol::Write_Footprint(DOCKMol & mol, ostream & ofs)
//{
//return true;
//}
/***********************************************************************/
void
copy_molecule(DOCKMol & target, const DOCKMol & original)
{
    int             i;

    target.allocate_arrays(original.num_atoms, original.num_bonds, original.num_residues);

    // copy scalar data
    target.title = original.title;
    target.mol_info_line = original.mol_info_line;
    target.comment1 = original.comment1;
    target.comment2 = original.comment2;
    target.comment3 = original.comment3;
    target.simplex_text = original.simplex_text;
    target.energy = original.energy;
    target.parent = original.parent; // for GA parent status

    target.mol_data = original.mol_data;

    target.num_atoms = original.num_atoms;
    target.num_bonds = original.num_bonds;
    target.num_residues = original.num_residues; // this is the number in this molecule
    target.heavy_atoms = original.heavy_atoms;

    target.num_active_atoms = original.num_active_atoms;
    target.num_active_bonds = original.num_active_bonds;
    target.score_text_data = original.score_text_data;

    target.amber_at_assigned = original.amber_at_assigned;
    target.amber_bt_assigned = original.amber_bt_assigned;
    target.chem_types_assigned = original.chem_types_assigned;
    target.ph4_types_assigned = original.ph4_types_assigned;

    target.current_score = original.current_score;
    // Score components - CS - 06-05-16
    target.score_nrg = original.score_nrg;          // Energy Score
    target.score_mg_nrg = original.score_mg_nrg;    // MG Score
    target.score_cont_nrg = original.score_cont_nrg;// Continuous Score
    target.score_fps = original.score_fps;          // Footprint Score
    target.score_ph4 = original.score_ph4;          // Ph4 Score
    target.score_tan = original.score_tan;        // Fingerprint Score
    target.score_hun = original.score_hun;          // Hungarian Score
    target.score_vol = original.score_vol;          // Volume Score
    target.internal_energy = original.internal_energy;
    target.current_data = original.current_data;
    target.primary_data = original.primary_data;
    target.hbond_text_data = original.hbond_text_data;

    target.rank = original.rank;
    target.fitness = original.fitness;
    target.crowding_dist = original.crowding_dist;


    target.grid_num = original.grid_num;
    target.amber_score_ligand_id = original.amber_score_ligand_id;
    target.total_dsol = original.total_dsol;

    target.rot_bonds = original.rot_bonds;    //YUCHEN
    target.mol_wt = original.mol_wt;
    target.formal_charge = original.formal_charge;
    target.atom_keys_0 = original.atom_keys_0;
    target.atom_keys_1 = original.atom_keys_1;
    target.atom_keys_2 = original.atom_keys_2;
    target.hb_donors = original.hb_donors;    //CS-09-15-2017
    target.hb_acceptors = original.hb_acceptors;
    target.mol_dnm_flag = original.mol_dnm_flag; // JDB
    ////  xlogp code : commented out sudipto and trent
    //target.number_of_H_Donors = original.number_of_H_Donors;
    //target.number_of_H_Acceptors = original.number_of_H_Acceptors;

    // GDRM, 2020_10_13
    #ifdef BUILD_DOCK_WITH_RDKIT
    target.num_stereocenters = original.num_stereocenters;
    target.num_spiro_atoms = original.num_spiro_atoms;
    target.clogp = original.clogp;
    target.tpsa = original.tpsa;
    target.num_arom_rings = original.num_arom_rings;
    target.num_alip_rings = original.num_alip_rings;
    target.num_sat_rings = original.num_sat_rings;
    target.smiles = original.smiles;
    target.qed_score = original.qed_score;
    target.sa_score = original.sa_score;
    target.esol = original.esol;
    target.fail_clogp = original.fail_clogp;
    target.fail_esol = original.fail_esol;
    target.fail_qed = original.fail_qed;
    target.fail_sa = original.fail_sa;
    target.fail_stereo = original.fail_stereo;
    target.pns = original.pns;
    target.pns_name = original.pns_name;
    target.MACCS = original.MACCS;
    #endif

    // copy arrays
    for (i = 0; i < original.num_atoms; i++) {

        target.atom_data[i] = original.atom_data[i];

        target.x[i] = original.x[i];
        target.y[i] = original.y[i];
        target.z[i] = original.z[i];

        //hbond
        target.flag_acceptor[i] = original.flag_acceptor[i]; // this is a h-bond acceptor
        target.flag_donator[i]  = original.flag_donator[i]; // this is true for polar hydrogen
        target.acc_heavy_atomid[i] = original.acc_heavy_atomid[i]; // atom id of acceptor that is connected to polar h
                           // equal to zero if not polar hydrogen

        target.charges[i] = original.charges[i];
        target.atom_types[i] = original.atom_types[i];
        target.atom_names[i] = original.atom_names[i];

        //footprint info.
        target.atom_number[i] = original.atom_number[i];
        target.atom_residue_numbers[i] = original.atom_residue_numbers[i];
        target.subst_names[i] = original.subst_names[i];

        target.atom_color[i] = original.atom_color[i];  // kxr205
        target.atom_psol[i] = original.atom_psol[i];    // kxr205
        target.atom_apsol[i] = original.atom_apsol[i];  // kxr205
        target.atom_dnm_flag[i] = original.atom_dnm_flag[i]; // JDB
        //  xlogp
        // commmented out: sudipto and trent
        /* target.H_bond_donor[i]    = original.H_bond_donor[i];
        target.H_bond_acceptor[i] = original.H_bond_acceptor[i]; */
        // uncommmented out because we're using XLogP // GDRM (check if really needed)
        //target.H_bond_donor[i]    = original.H_bond_donor[i];
        //target.H_bond_acceptor[i] = original.H_bond_acceptor[i];


        target.atom_ring_flags[i] = original.atom_ring_flags[i];
        target.atom_active_flags[i] = original.atom_active_flags[i];

        target.amber_at_id[i] = original.amber_at_id[i];
        target.amber_at_radius[i] = original.amber_at_radius[i];
        target.amber_at_well_depth[i] = original.amber_at_well_depth[i];
        target.amber_at_heavy_flag[i] = original.amber_at_heavy_flag[i];
        target.amber_at_valence[i] = original.amber_at_valence[i];
        target.amber_at_bump_id[i] = original.amber_at_bump_id[i];

        target.atom_envs[i] = original.atom_envs[i];
        target.chem_types[i] = original.chem_types[i];
        target.ph4_types[i] = original.ph4_types[i];

        target.gb_hawkins_radius[i] = original.gb_hawkins_radius[i];
        target.gb_hawkins_scale[i] = original.gb_hawkins_scale[i];

        target.neighbor_list[i] = original.neighbor_list[i];

        // DTM - 11-12-08 - copy the new atom_segment_ids array
        target.atom_segment_ids[i] = original.atom_segment_ids[i];
    }

    
    target.dummy1  = original.dummy1;
    target.dummy2  = original.dummy2;
    target.atomtag = original.atomtag;

    for (i = 0; i < original.num_bonds; i++) {

        target.bonds_origin_atom[i] = original.bonds_origin_atom[i];
        target.bonds_target_atom[i] = original.bonds_target_atom[i];
        target.bond_types[i] = original.bond_types[i];
        target.bond_segment_ids[i] = original.bond_segment_ids[i]; // CDS - 04/09/15

        target.bond_ring_flags[i] = original.bond_ring_flags[i];
        target.bond_active_flags[i] = original.bond_active_flags[i];
        target.bond_keep_flags[i] = original.bond_keep_flags[i];

        target.amber_bt_id[i] = original.amber_bt_id[i];
        target.amber_bt_minimize[i] = original.amber_bt_minimize[i];
        target.amber_bt_torsion_total[i] = original.amber_bt_torsion_total[i];
        target.amber_bt_torsions[i] = original.amber_bt_torsions[i];

    }
    // footprints
    target.footprints.clear();
    FOOTPRINT_ELEMENT temp_footprint_ele;
    int num_residues_in_rec = original.footprints.size();
    // loops over all residues and appends the values of the orginal mol
    // on to the new mol
    for (i = 0; i < num_residues_in_rec; i++) {
        temp_footprint_ele.resname = original.footprints[i].resname;
        temp_footprint_ele.resid   = original.footprints[i].resid;
        temp_footprint_ele.vdw     = original.footprints[i].vdw;
        temp_footprint_ele.es      = original.footprints[i].es;
        temp_footprint_ele.hb      = original.footprints[i].hb;
        target.footprints.push_back(temp_footprint_ele);
    }

/*  //pharmacophore
    target.ph4.clear(); //LINGLING
    PH4_ELEMENT temp_ph4_ele;
    //initialize temp_ph4_ele using original pharmacophore (?)
*/

    // test_child_list = new bool[2*num_bonds*num_atoms];

    // for(i=0;i<original.num_atoms*original.num_atoms;i++)
    // target.ie_neighbor_list[i] = original.ie_neighbor_list[i];

    for (i = 0; i < 2 * original.num_bonds; i++)
        target.atom_child_list[i] = original.atom_child_list[i];

}

//DOCKMol::initialize_from_mol2(std::ifstream &)
//{
//
//
//}
/**********************************************************************/
void
copy_molecule_shallow(DOCKMol & target, const DOCKMol & original)
{
    int             i;

    target.allocate_arrays(original.num_atoms, original.num_bonds, original.num_residues);
    target.title = original.title;
    target.mol_info_line = original.mol_info_line;
    target.energy = original.energy;
    target.parent = original.parent; // for GA parent status
    target.mol_data = original.mol_data;
    target.num_atoms = original.num_atoms;
    target.num_bonds = original.num_bonds;
    target.num_residues = original.num_residues; // this is the number in this molecule
    target.heavy_atoms = original.heavy_atoms;
    target.num_active_atoms = original.num_active_atoms;
    target.num_active_bonds = original.num_active_bonds;
    target.amber_at_assigned = original.amber_at_assigned;
    target.amber_bt_assigned = original.amber_bt_assigned;
    target.rot_bonds = original.rot_bonds;    //YUCHEN
    target.mol_wt = original.mol_wt;
    target.atom_keys_0 = original.atom_keys_0;
    target.atom_keys_1 = original.atom_keys_1;
    target.atom_keys_2 = original.atom_keys_2;
    target.mol_dnm_flag = original.mol_dnm_flag; // JDB

    for (i = 0; i < original.num_atoms; i++) {

        target.atom_data[i] = original.atom_data[i];

        target.x[i] = original.x[i];
        target.y[i] = original.y[i];
        target.z[i] = original.z[i];

        target.acc_heavy_atomid[i] = original.acc_heavy_atomid[i]; // atom id of acceptor that is connected to polar h
        target.atom_types[i] = original.atom_types[i];
        target.atom_names[i] = original.atom_names[i];
        target.atom_number[i] = original.atom_number[i];
        target.atom_active_flags[i] = original.atom_active_flags[i];
        target.subst_names[i] = original.subst_names[i];
        target.atom_dnm_flag[i] = original.atom_dnm_flag[i]; // JDB

        target.amber_at_id[i] = original.amber_at_id[i];
        target.atom_envs[i] = original.atom_envs[i];
        target.neighbor_list[i] = original.neighbor_list[i];
        target.atom_segment_ids[i] = original.atom_segment_ids[i];
    }

    for (i = 0; i < original.num_bonds; i++) {

        target.bonds_origin_atom[i] = original.bonds_origin_atom[i];
        target.bonds_target_atom[i] = original.bonds_target_atom[i];
        target.bond_types[i] = original.bond_types[i];
        target.bond_segment_ids[i] = original.bond_segment_ids[i]; // CDS - 04/09/15

        target.bond_active_flags[i] = original.bond_active_flags[i];
        target.bond_keep_flags[i] = original.bond_keep_flags[i];

        target.amber_bt_id[i] = original.amber_bt_id[i];
        target.amber_bt_minimize[i] = original.amber_bt_minimize[i];
        target.amber_bt_torsion_total[i] = original.amber_bt_torsion_total[i];
        target.amber_bt_torsions[i] = original.amber_bt_torsions[i];

    }

    for (i = 0; i < 2 * original.num_bonds; i++) {
        target.atom_child_list[i] = original.atom_child_list[i];

    }

}



/***********************************************************************/
void
copy_crds(DOCKMol & target, DOCKMol & original)
{
    int             i;

    for (i = 0; i < original.num_atoms; i++) {
        target.x[i] = original.x[i];
        target.y[i] = original.y[i];
        target.z[i] = original.z[i];
    }

}

/*************************************/
DOCKMol::DOCKMol()
{

    initialize();

}

/*************************************/
DOCKMol::DOCKMol(const DOCKMol & original)
{

    initialize();
    copy_molecule(*this, original);

}

/*************************************/
DOCKMol::~DOCKMol()
{

    clear_molecule();

}

/*************************************/
void
DOCKMol::operator=(const DOCKMol & original)
{

    copy_molecule(*this, original);

}


/*************************************/
void
DOCKMol::input_parameters(Parameter_Reader & parm)
{


    // read_solvation = (parm.query_param("read_mol_solvation","no","yes no")
    // == "yes")?true:false;

    // read_color = (parm.query_param("read_mol_color","no","yes no") ==
    // "yes")?true:false;


  /***
  string tmp;

  cout << "\nMolecule Parameters" << endl;
  cout << "------------------------------------------------------------------------------------------" << endl;

  tmp = parm.query_param("calc_internal_energy", "no", "yes no");
  if(tmp == "yes")
	  use_internal_energy = true;
  else
	  use_internal_energy = false;

  if(use_internal_energy) {
	  att_exp = atoi(parm.query_param("internal_energy_att_exp", "6").c_str());
	  rep_exp = atoi(parm.query_param("internal_energy_rep_exp", "12").c_str());
	  dielectric = atof(parm.query_param("internal_energy_dielectric", "4.0").c_str());

	  // internal energy hardcoded params
	  dielectric = 332.0 / dielectric;
  }
  ***/
    // use_internal_energy = false;

}

/*************************************/
void
DOCKMol::initialize()
{

    arrays_allocated = false;
    // use_internal_energy = false;
    clear_molecule();           // /////////////////////

}

/*************************************/
void
DOCKMol::allocate_arrays(int natoms, int nbonds, int nresidues) //added jwu
{
    int             i;

    clear_molecule();

    num_atoms    = natoms;
    num_bonds    = nbonds;
    num_residues = nresidues; // number in current mol

//    num_residues_in_rec = 1000; // number in current rec


    atom_data = new string[num_atoms];

    x = new float[num_atoms];
    y = new float[num_atoms];
    z = new float[num_atoms];

    charges = new float[num_atoms];
    atom_types = new string[num_atoms];
    // hbond
    flag_acceptor    = new bool[num_atoms];
    flag_donator     = new bool[num_atoms];
    acc_heavy_atomid = new int[num_atoms];

    //jwu code
    atom_number=new string[num_atoms];
    atom_residue_numbers=new string[num_atoms];
    atom_names = new string[num_atoms];
    subst_names = new string[num_atoms];
    atom_color = new string[num_atoms]; // kxr205
    atom_psol = new float[num_atoms];   // kxr205
    atom_apsol = new float[num_atoms];  // kxr205

    // DTM - 11-12-08 - allocate the atom_segment_ids array
    atom_segment_ids = new int[num_atoms];
    bond_segment_ids = new int[num_bonds]; // CDS - 04/09/15
    atom_dnm_flag = new bool[num_atoms];   // LEP
    mol_dnm_flag = false; // JDB

    //// added for inter molecular H-Bonds in xlogp
    //// commented out by sudipto and trent
    //// uncommented out by GDRM
    //number_of_H_Donors = 0;
    //number_of_H_Acceptors = 0;
    //H_bond_donor = new int[num_atoms];
    //H_bond_acceptor = new int[num_atoms];
    //xlogp = 0.0;

    atom_ring_flags = new bool[num_atoms];
    bond_ring_flags = new bool[num_bonds];

    bonds_origin_atom = new int[num_bonds];
    bonds_target_atom = new int[num_bonds];
    bond_types = new string[num_bonds];

    atom_active_flags = new bool[num_atoms];
    bond_active_flags = new bool[num_bonds];
    bond_keep_flags = new bool[num_bonds];

    num_active_atoms = num_atoms;
    num_active_bonds = num_bonds;

    amber_at_bump_id = new int[num_atoms];
    amber_at_heavy_flag = new int[num_atoms];
    amber_at_id = new int[num_atoms];
    amber_at_valence = new int[num_atoms];
    amber_at_radius = new float[num_atoms];
    amber_at_well_depth = new float[num_atoms];

    atom_envs  = new string[num_atoms]; //BCF fingerprint for gasteiger
    chem_types = new string[num_atoms];
    ph4_types  = new string[num_atoms];

    amber_bt_id = new int[num_bonds];
    amber_bt_minimize = new int[num_bonds];
    amber_bt_torsion_total = new int[num_bonds];
    amber_bt_torsions = new TOR_LIST[num_bonds];

    gb_hawkins_radius = new float[num_atoms];
    gb_hawkins_scale = new float[num_atoms];

    // footprints.
    //footprints.clear(); // performed in clear_molecule()

    neighbor_list = new INTVec[num_atoms];
    // ie_neighbor_list = new bool[num_atoms*num_atoms];

    atom_child_list = new INTVec[num_bonds * 2];

    arrays_allocated = true;

    // assign init values
    for (i = 0; i < num_atoms; i++) {

        atom_data[i] = "";
        atom_names[i] = "";

        //footprint
        atom_residue_numbers[i] ="";
        atom_names[i] = "";

        subst_names[i] = "";
        atom_types[i] = "";
        atom_color[i] = "null"; // kxr205
        atom_psol[i] = 0.0;     // kxr205
        atom_apsol[i] = 0.0;
        charges[i] = 0.0;

        // DTM - 11-12-08 - initialize the atom_segment_ids array
       	atom_segment_ids[i] = -1;

        // hbond
        flag_acceptor[i]    = false;
        flag_donator[i]     = false;
        acc_heavy_atomid[i] = 0;

        atom_ring_flags[i] = false;
        atom_active_flags[i] = true;

        amber_at_bump_id[i] = 0;
        amber_at_heavy_flag[i] = 0;
        amber_at_id[i] = 0;
        amber_at_valence[i] = 0;
        amber_at_radius[i] = 0.0;
        amber_at_well_depth[i] = 0.0;

        atom_envs[i] = "";
        chem_types[i] = "";
        ph4_types[i] = "";

        gb_hawkins_radius[i] = 0.0;
        gb_hawkins_scale[i] = 0.0;

    }

    for (i = 0; i < num_bonds; i++) {

        bonds_origin_atom[i] = 0;
        bonds_target_atom[i] = 0;
        bond_types[i] = "";

        bond_segment_ids[i] = -1;

        bond_ring_flags[i] = false;
        bond_active_flags[i] = true;
        bond_keep_flags[i] = false;

        amber_bt_id[i] = 0;
        amber_bt_minimize[i] = 0;
        amber_bt_torsion_total[i] = 0;
    }
    // test_child_list = new bool[2*num_bonds*num_atoms];

    #ifdef BUILD_DOCK_WITH_RDKIT
    num_stereocenters = 0;
    num_spiro_atoms = 0;
    clogp = 0.0;
    tpsa = 0.0;
    num_arom_rings = 0;
    num_alip_rings = 0;
    num_sat_rings = 0;
    smiles = "";
    qed_score = 0.0;
    sa_score = 0.0;
    esol = 0.0;
    fail_clogp = false;
    fail_esol = false;
    fail_qed = false;
    fail_sa = false;
    fail_stereo = false;
    pns = 0.0;
    pns_name = {};
    MACCS = {};
    #endif

}

/*************************************/
void
DOCKMol::clear_molecule()
{
    int             i;

    amber_at_assigned = false;
    amber_bt_assigned = false;


    if (arrays_allocated) {

        for (i = 0; i < num_bonds; i++)
            amber_bt_torsions[i].clear();

        for (i = 0; i < num_atoms; i++)
            neighbor_list[i].clear();

        for (i = 0; i < 2 * num_bonds; i++)
            atom_child_list[i].clear();

        delete[]atom_data;
        delete[]x;
        delete[]y;
        delete[]z;

        // hbond
        delete[]flag_acceptor;
        delete[]flag_donator;
        delete[]acc_heavy_atomid;

        //// Arrays obtained from within xlogp.cpp
        //// Commented out since we are not defining this any more: sudipto and trent
        //delete[]H_bond_donor;
        //delete[]H_bond_acceptor;
        //xlogp = 0.0;

        delete[]charges;
        delete[]atom_types;
        delete[]atom_names;

        //footprint info.
        delete[]atom_number;
        delete[]atom_residue_numbers;
        delete[]subst_names; // residue name

        delete[]atom_color;     // kxr205
        delete[]atom_psol;      // kxr205
        delete[]atom_apsol;

	// DTM - 11-12-08 - delete the atom_segment_ids array
        delete[]atom_segment_ids;
        delete[]bond_segment_ids;

        delete[]bonds_origin_atom;
        delete[]bonds_target_atom;
        delete[]bond_types;

        delete[]atom_ring_flags;
        delete[]bond_ring_flags;
        delete[]bond_keep_flags;

        delete[]atom_active_flags;
        delete[]bond_active_flags;
        delete[]atom_dnm_flag;

        delete[]amber_at_id;
        delete[]amber_at_radius;
        delete[]amber_at_well_depth;
        delete[]amber_at_heavy_flag;
        delete[]amber_at_valence;
        delete[]amber_at_bump_id;

        delete[]atom_envs;
        delete[]chem_types;
        delete[]ph4_types;

        delete[]amber_bt_id;
        delete[]amber_bt_minimize;
        delete[]amber_bt_torsion_total;
        delete[]amber_bt_torsions;

        delete[]gb_hawkins_radius;
        delete[]gb_hawkins_scale;

        // footprints
        footprints.clear();

        delete[]neighbor_list;
        // delete [] ie_neighbor_list;
        // delete [] test_child_list;
        delete[]atom_child_list;

        arrays_allocated = false;
    }
    // clear MOL2 header info
    title = "";
    mol_info_line = "";
    comment1 = "";
    comment2 = "";
    comment3 = "";
    simplex_text = "";
    energy = "";
    mol_data = "";

    // clear scalar data
    num_atoms = 0;
    num_bonds = 0;
    num_residues = 0;
    //num_residues_in_rec = 0;
    num_active_atoms = 0;
    num_active_bonds = 0;
    score_text_data = "";
    current_score = 0.0;
    rot_bonds = 0;
    heavy_atoms = 0;
    hb_donors = 0;
    hb_acceptors = 0;
    mol_wt = 0.0;
    formal_charge = 0.0;
    xlogp = 0.0;
    // Score components for Desc Score - CS 06-05-16
    score_nrg = 0.0;     // Energy Score
    score_mg_nrg = 0.0;  // MG Score
    score_cont_nrg = 0.0;// Continuous Score
    score_fps = 0.0;     // Footprint Score
    score_ph4 = 0.0;     // Ph4 Score
    score_tan = 0.0;     // Fingerprint Score
    score_hun = 0.0;     // Hungarian Score
    score_vol = 0.0;     // Volume Score
    // GA Rank & Crowding - CS 06-06-16
    rank = 0;         // Rank between 0 - n molecules
    fitness = 0.0;    // Individual molecule fitness score
    crowding_dist = 0.0;// Distance from 0 - 9999
    internal_energy = 0.0;
    current_data = "ERROR: Conformation could not be scored by DOCK.\n";
    primary_data = "ERROR: Conformation could not be scored by DOCK.\n";
    hbond_text_data = "";
    amber_score_ligand_id = "";
    grid_num = 0; // This defaults to the first grid read in
    total_dsol = 0.0;
    mol_dnm_flag=false; //JDB

    #ifdef BUILD_DOCK_WITH_RDKIT
    // GDRM, Feb 26, 2020
    // Clear RDKit-related data
    num_stereocenters = 0;
    num_spiro_atoms = 0;
    clogp = 0.0;
    tpsa = 0.0;
    num_arom_rings = 0;
    num_alip_rings = 0;
    num_sat_rings = 0;
    smiles = "";
    qed_score = 0.0;
    sa_score = 0.0;
    esol = 0.0;
    pns = 0.0;
    pns_name = {};
    MACCS = {};
    #endif

}

/*************************************/
vector < int   >
DOCKMol::get_atom_neighbors(int index)
{
    vector < int   >nbrs;
    int             i;

    nbrs.clear();

    // loop over bonds and find any that include the index atom
    for (i = 0; i < num_bonds; i++) {

        if (bonds_origin_atom[i] == index)
            nbrs.push_back(bonds_target_atom[i]);

        if (bonds_target_atom[i] == index)
            nbrs.push_back(bonds_origin_atom[i]);

    }

    return nbrs;
}

/*************************************/
vector < int   >
DOCKMol::get_bond_neighbors(int index)
{
    vector < int   >nbrs;
    int             i;

    nbrs.clear();

    // loop over bonds- id any that contain index atom
    for (i = 0; i < num_bonds; i++) {

        if ((bonds_origin_atom[i] == index) || (bonds_target_atom[i] == index))
            nbrs.push_back(i);

    }
    return nbrs;
}

/*************************************/
int
DOCKMol::get_bond(int a1, int a2)
{
    int             bond_id;
    int             i;

    bond_id = -1;

    for (i = 0; i < num_bonds; i++) {

        if ((bonds_origin_atom[i] == a1) && (bonds_target_atom[i] == a2))
            bond_id = i;

        if ((bonds_origin_atom[i] == a2) && (bonds_target_atom[i] == a1))
            bond_id = i;

    }

    return bond_id;
}

/*************************************/
void
DOCKMol::id_ring_atoms_bonds()
{
    vector < bool > atoms_visited,
                    bonds_visited;
    vector < int  > atom_path,
                    bond_path;
    int             i;

    atoms_visited.clear();
    bonds_visited.clear();

    atoms_visited.resize(num_atoms, false);
    bonds_visited.resize(num_bonds, false);
    atom_path.clear();
    bond_path.clear();

    for (i = 0; i < num_atoms; i++)
        atoms_visited[i] = false;

    for (i = 0; i < num_bonds; i++)
        bonds_visited[i] = false;

    // loop over all atoms and find rings
    for (i = 0; i < num_atoms; i++)
        if (!atoms_visited[i])
            find_rings(atom_path, bond_path, atoms_visited, bonds_visited, i);

}

/*************************************/
void
DOCKMol::find_rings(vector < int >apath, vector < int >bpath,
                    vector < bool > &atoms, vector < bool > &bonds, int atnum)
{
    int             i,
                    j,
                    nbr_bond;
    vector < int  > nbrs;

    if (atoms[atnum]) {

        i = apath.size() - 1;
        j = bpath.size() - 1;

        while ((i >= 0) && (j >= 0)) {  // ///// changed from > to >= fixes the
                                        // first atom bug
            atom_ring_flags[apath[i--]] = true;
            bond_ring_flags[bpath[j--]] = true;

            if (i == -1)        // added fix to address if i = 0
                break;
            else if (apath[i] == atnum)
                break;

        }

    } else {

        atoms[atnum] = true;
        nbrs = get_atom_neighbors(atnum);

        for (long unsigned int i = 0; i < nbrs.size(); i++) {

            nbr_bond = get_bond(nbrs[i], atnum);

            if (!bonds[nbr_bond]) {

                bonds[nbr_bond] = true;

                apath.push_back(nbrs[i]);
                bpath.push_back(nbr_bond);

                find_rings(apath, bpath, atoms, bonds, nbrs[i]);
                apath.pop_back();
                bpath.pop_back();
            }

        }

    }

}

/*************************************/
vector < int   >
DOCKMol::get_atom_children(int a1, int a2)
{
    vector < int   >children;
    bool           *visited;
    vector < int   >nbrs,
                    new_nbrs;

    if (get_bond(a1, a2) != -1) {
        visited = new bool[num_atoms];
        memset(visited, 0, num_atoms * sizeof(bool));
        // visited.resize(num_atoms, false);

        visited[a1] = true;
        visited[a2] = true;

        // nbrs = get_atom_neighbors(a2);
        nbrs = neighbor_list[a2];

        while (nbrs.size() > 0) {

            if (!visited[nbrs[nbrs.size() - 1]]) {

                children.push_back(nbrs[nbrs.size() - 1]);
                new_nbrs.clear();
                // new_nbrs = get_atom_neighbors(nbrs[nbrs.size()-1]);
                new_nbrs = neighbor_list[nbrs[nbrs.size() - 1]];
                visited[nbrs[nbrs.size() - 1]] = true;
                nbrs.pop_back();
                for (long unsigned int i = 0; i < new_nbrs.size(); i++)
                    nbrs.push_back(new_nbrs[i]);

            } else {
                nbrs.pop_back();
            }

        }

        delete[]visited;
    }

    return children;
}

/*************************************/
float
DOCKMol::get_torsion(int a1, int a2, int a3, int a4)
{
    float           torsion;
    DOCKVector      v1,
                    v2,
                    v3,
                    v4;

    v1.x = x[a1];
    v1.y = y[a1];
    v1.z = z[a1];

    v2.x = x[a2];
    v2.y = y[a2];
    v2.z = z[a2];

    v3.x = x[a3];
    v3.y = y[a3];
    v3.z = z[a3];

    v4.x = x[a4];
    v4.y = y[a4];
    v4.z = z[a4];

    torsion = get_torsion_angle(v1, v2, v3, v4);

    return torsion;
}

/*************************************/
// This function seems to be using simple rotation matrices
// Why not use the fancy quternion stuff? :sudipto
void
DOCKMol::set_torsion(int a1, int a2, int a3, int a4, float angle)
{
    int             tor[4];
    vector < int   >atoms;
    float           v1x,
                    v1y,
                    v1z,
                    v2x,
                    v2y,
                    v2z,
                    v3x,
                    v3y,
                    v3z;
    float           c1x,
                    c1y,
                    c1z,
                    c2x,
                    c2y,
                    c2z,
                    c3x,
                    c3y,
                    c3z;
    float           c1mag,
                    c2mag,
                    radang,
                    costheta,
                    m[9];
    float           nx,
                    ny,
                    nz,
                    mag,
                    rotang,
                    sn,
                    cs,
                    t,
                    tx,
                    ty,
                    tz;
    int             j,
                    idx;

    tor[0] = a1;
    tor[1] = a2;
    tor[2] = a3;
    tor[3] = a4;

    idx = get_bond(a2, a3);

    if (a2 < a3)
        idx = 2 * idx;
    else
        idx = 2 * idx + 1;

    atoms = atom_child_list[idx];
    // atoms = get_atom_children(a2, a3);
    // atoms = child_list[a2][a3]; //////////////////

    // calculate the torsion angle
    v1x = x[tor[0]] - x[tor[1]];
    v2x = x[tor[1]] - x[tor[2]];
    v1y = y[tor[0]] - y[tor[1]];
    v2y = y[tor[1]] - y[tor[2]];
    v1z = z[tor[0]] - z[tor[1]];
    v2z = z[tor[1]] - z[tor[2]];
    v3x = x[tor[2]] - x[tor[3]];
    v3y = y[tor[2]] - y[tor[3]];
    v3z = z[tor[2]] - z[tor[3]];

    c1x = v1y * v2z - v1z * v2y;
    c2x = v2y * v3z - v2z * v3y;
    c1y = -v1x * v2z + v1z * v2x;
    c2y = -v2x * v3z + v2z * v3x;
    c1z = v1x * v2y - v1y * v2x;
    c2z = v2x * v3y - v2y * v3x;
    c3x = c1y * c2z - c1z * c2y;
    c3y = -c1x * c2z + c1z * c2x;
    c3z = c1x * c2y - c1y * c2x;

    c1mag = pow(c1x, 2) + pow(c1y, 2) + pow(c1z, 2);
    c2mag = pow(c2x, 2) + pow(c2y, 2) + pow(c2z, 2);

    if (c1mag * c2mag < 0.01)
        costheta = 1.0;         // avoid div by zero error
    else
        costheta = (c1x * c2x + c1y * c2y + c1z * c2z) / (sqrt(c1mag * c2mag));

    if (costheta < -0.999999)
        costheta = -0.999999f;
    if (costheta > 0.999999)
        costheta = 0.999999f;

    if ((v2x * c3x + v2y * c3y + v2z * c3z) > 0.0)
        radang = -acos(costheta);
    else
        radang = acos(costheta);

    //
    // now we have the torsion angle (radang) - set up the rot matrix
    //

    // find the difference between current and requested
    rotang = angle - radang;

    sn = sin(rotang);
    cs = cos(rotang);
    t = 1 - cs;

    // normalize the rotation vector
    mag = sqrt(pow(v2x, 2) + pow(v2y, 2) + pow(v2z, 2));
    nx = v2x / mag;
    ny = v2y / mag;
    nz = v2z / mag;

    // set up the rotation matrix
    m[0] = t * nx * nx + cs;
    m[1] = t * nx * ny + sn * nz;
    m[2] = t * nx * nz - sn * ny;
    m[3] = t * nx * ny - sn * nz;
    m[4] = t * ny * ny + cs;
    m[5] = t * ny * nz + sn * nx;
    m[6] = t * nx * nz + sn * ny;
    m[7] = t * ny * nz - sn * nx;
    m[8] = t * nz * nz + cs;

    //
    // now the matrix is set - time to rotate the atoms
    //
    tx = x[tor[1]];
    ty = y[tor[1]];
    tz = z[tor[1]];

    for (long unsigned int i = 0; i < atoms.size(); i++) {
        j = atoms[i];

        // for(i=0;i<num_atoms;i++) {
        // j = i;

        // if(child_list[a2*num_atoms + i] == a3) { //////////////

        x[j] -= tx;
        y[j] -= ty;
        z[j] -= tz;

        nx = x[j] * m[0] + y[j] * m[1] + z[j] * m[2];
        ny = x[j] * m[3] + y[j] * m[4] + z[j] * m[5];
        nz = x[j] * m[6] + y[j] * m[7] + z[j] * m[8];

        x[j] = nx;
        y[j] = ny;
        z[j] = nz;
        x[j] += tx;
        y[j] += ty;
        z[j] += tz;

        // } /////////////////////////
    }

}

/*********************************************/
void
DOCKMol::translate_mol(const DOCKVector & vec)
{
    int             i;

    for (i = 0; i < num_atoms; i++) {
        x[i] += vec.x;
        y[i] += vec.y;
        z[i] += vec.z;
    }

}

/*********************************************/
void
DOCKMol::translate_mol(const DOCKVector & vec, bool flag_all)
{
    int             i;

    for (i = 0; i < num_atoms; i++) {
        if (atom_active_flags[i] || flag_all){ // if flag_all is true then always enter regardless of active state
          x[i] += vec.x;
          y[i] += vec.y;
          z[i] += vec.z;
        }
    }

}

/*********************************************/
void
DOCKMol::rotate_mol(double mat[9],bool flag_all)
{
    int             i;
    double          nx,
                    ny,
                    nz;

    for (i = 0; i < num_atoms; i++) {
        if (atom_active_flags[i] || flag_all){
            nx = x[i];
            ny = y[i];
            nz = z[i];
           
            x[i] = mat[0] * nx + mat[1] * ny + mat[2] * nz;
            y[i] = mat[3] * nx + mat[4] * ny + mat[5] * nz;
            z[i] = mat[6] * nx + mat[7] * ny + mat[8] * nz;
        }
    }

}

/*********************************************/
void
DOCKMol::rotate_mol(double mat[3][3], bool flag_all)
{
    double          new_mat[9];
    int             i,
                    j,
                    k;

    k = 0;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            new_mat[k++] = mat[i][j];

    rotate_mol(new_mat, flag_all);

}

/*********************************************/
void
DOCKMol::rotate_mol(double mat[3][3])
{
    double          new_mat[9];
    int             i,
                    j,
                    k;

    k = 0;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            new_mat[k++] = mat[i][j];

    rotate_mol(new_mat);

}

/*********************************************/
void
DOCKMol::rotate_mol(double mat[9])
{
    int             i;
    double          nx,
                    ny,
                    nz;

    for (i = 0; i < num_atoms; i++) {
        nx = x[i];
        ny = y[i];
        nz = z[i];

        x[i] = mat[0] * nx + mat[1] * ny + mat[2] * nz;
        y[i] = mat[3] * nx + mat[4] * ny + mat[5] * nz;
        z[i] = mat[6] * nx + mat[7] * ny + mat[8] * nz;
    }

}
/*********************************************/
bool
DOCKMol::atoms_are_one_three(int a1, int a2)
{
    vector < int   >nbrs1,
                    nbrs2;

    nbrs1 = get_atom_neighbors(a1);
    nbrs2 = get_atom_neighbors(a2);

    for (long unsigned int i = 0; i < nbrs1.size(); i++) {
        for (long unsigned int j = 0; j < nbrs2.size(); j++) {
            if (nbrs1[i] == nbrs2[j])
                return true;
        }
    }

    return false;
}

/*********************************************/
bool
DOCKMol::atoms_are_one_four(int a1, int a2)
{
    vector < int   >nbrs1,
                    nbrs2;

    nbrs1 = get_atom_neighbors(a1);
    nbrs2 = get_atom_neighbors(a2);

    for (long unsigned int i = 0; i < nbrs1.size(); i++) {
        for (long unsigned int j = 0; j < nbrs2.size(); j++) {
            if (get_bond(nbrs1[i], nbrs2[j]) != -1)
                return true;
        }
    }

    return false;
}

/*********************************************
float DOCKMol::compute_internal_energy() {
        int             i,j;
        float   total_energy;
        float   total_vdw;
        float   total_es;
        float   distance;

        total_energy = total_vdw = total_es = 0.0;

        if(use_internal_energy) {

                for(i=0;i<num_atoms;i++) {

                        for(j=i+1;j<num_atoms;j++) {

                                if(atom_active_flags[i] && atom_active_flags[j]) {

                                        if( (get_bond(i, j)==-1) && (!atoms_are_one_three(i, j)) && (!atoms_are_one_four(i, j)) ) {
//                                      if(!ie_neighbor_list[i*num_atoms + j]) {

                                                distance = pow((x[i]-x[j]),2) + pow((y[i]-y[j]),2) + pow((z[i]-z[j]),2);
//                                              distance = sqrt(distance);

                                                total_vdw += (vdwA[i]*vdwA[j])/pow(distance, (rep_exp)) - (vdwB[i]*vdwB[j])/pow(distance, (att_exp)); // removed /2 term
                                                total_es  += (charges[i] * charges[j] * dielectric) / pow(distance, 1/2);//sqrt(distance);

                                        }

                                }
                        }

                }

        }

        total_energy = total_vdw + total_es;

        return total_energy;
}
**/

/*************************************/
void
DOCKMol::prepare_molecule()
{
    int             i,
                    j,
                    idx;

    // pre-cache neighbor list
    for (i = 0; i < num_atoms; i++)
        neighbor_list[i] = get_atom_neighbors(i);

        /**
	// pre-cache ie_neighbor_list
	for(i=0;i<num_atoms;i++) {
		for(j=0;j<num_atoms;j++) {
			if((get_bond(i, j)==-1)&&(!atoms_are_one_three(i, j))&&(!atoms_are_one_four(i, j))&&(i!=j))
				ie_neighbor_list[i*num_atoms + j] = false;
			else
				ie_neighbor_list[i*num_atoms + j] = true;
		}
	}
	**/

    // pre-cache atom child list
    for (i = 0; i < num_atoms; i++) {
        for (j = 0; j < num_atoms; j++) {
            if (i != j) {
                idx = get_bond(i, j);
                if (idx != -1) {
                    if (i < j)
                        idx = 2 * idx;
                    else
                        idx = 2 * idx + 1;

                    atom_child_list[idx] = get_atom_children(i, j);
                }
            }
        }
    }


}

/*********************************************/
#ifdef BUILD_DOCK_WITH_RDKIT

void fixNitroSubstructureAndCharge(RDKit::RWMol &res, unsigned int atIdx) {
    unsigned int noODblNeighbors = 0;
    RDKit::ROMol::ADJ_ITER nbrIdxIt, nbrEndIdxIt;
    vector<unsigned int> toModIdx {};
    boost::tie(nbrIdxIt, nbrEndIdxIt) =
        res.getAtomNeighbors(res.getAtomWithIdx(atIdx));
    while (nbrIdxIt != nbrEndIdxIt) {
        RDKit::Bond *curBond = res.getBondBetweenAtoms(atIdx, *nbrIdxIt);
        if (res.getAtomWithIdx(*nbrIdxIt)->getAtomicNum() == 8 &&
            (curBond->getBondType() == RDKit::Bond::DOUBLE ||
             curBond->getBondType() == RDKit::Bond::AROMATIC)) {
            ++noODblNeighbors;
            toModIdx.push_back(*nbrIdxIt); // save the indices if condition is met
        }
        ++nbrIdxIt;
    }
    if (noODblNeighbors == 2) {
        unsigned int count = 0;
        RDKit::Atom *nbr = res.getAtomWithIdx(atIdx);
        for (vector<unsigned int>::iterator it = toModIdx.begin(); it < toModIdx.end(); ++it) {
            RDKit::Bond *curBond = res.getBondBetweenAtoms(atIdx, *it);
            if (count == (noODblNeighbors - 1)) {
                curBond->setBondType(RDKit::Bond::SINGLE);
                res.getAtomWithIdx(atIdx)->setFormalCharge(1);
                res.getAtomWithIdx(*it)->setFormalCharge(-1);
                res.getAtomWithIdx(*it)->setIsAromatic(false);
                nbr->setIsAromatic(false);
            } else {
                curBond->setBondType(RDKit::Bond::DOUBLE);
                res.getAtomWithIdx(*it)->setIsAromatic(false);
            }
            ++count;
        }
    }
}

void fixProtonatedAmidineSubstructureAndCharge(RDKit::RWMol &res, unsigned int atIdx){
    unsigned int noNArNeighbors = 0;
    unsigned int noCNeighbors = 0;
    RDKit::ROMol::ADJ_ITER nbrIdxIt, nbrEndIdxIt;
    vector<unsigned int> toModIdx {};
    boost::tie(nbrIdxIt, nbrEndIdxIt) =
        res.getAtomNeighbors(res.getAtomWithIdx(atIdx));
    while (nbrIdxIt != nbrEndIdxIt) {
        RDKit::Bond *curBond = res.getBondBetweenAtoms(atIdx, *nbrIdxIt);
        if (res.getAtomWithIdx(*nbrIdxIt)->getAtomicNum() == 7 &&
            curBond->getBondType() == RDKit::Bond::AROMATIC) {
            ++noNArNeighbors;
            toModIdx.push_back(*nbrIdxIt);
        } else if (res.getAtomWithIdx(*nbrIdxIt)->getAtomicNum() == 6) {
            ++noCNeighbors;
        }
        ++nbrIdxIt;
    }
    if (noNArNeighbors == 2 && noCNeighbors == 1) {
        unsigned int count = 0;
        RDKit::Atom *nbr = res.getAtomWithIdx(atIdx);
        for (vector<unsigned int>::iterator it = toModIdx.begin(); it < toModIdx.end(); ++it){
            RDKit::Bond *curBond = res.getBondBetweenAtoms(atIdx, *it);
            if (count == (noNArNeighbors - 1)) {
                string tATT;
                res.getAtomWithIdx(*it)->getProp(RDKit::common_properties::_TriposAtomType, tATT);
                if (tATT == "N.3"){
                    curBond->setBondType(RDKit::Bond::DOUBLE);
                    // set formal charge of +1 over the nitrogen
                    res.getAtomWithIdx(*it)->setFormalCharge(1);
                    res.getAtomWithIdx(*it)->setIsAromatic(false);
                    nbr->setIsAromatic(false);
                }
            } else {
                string tATT;
                res.getAtomWithIdx(*it)->getProp(RDKit::common_properties::_TriposAtomType, tATT);
                if (tATT == "N.3"){
                    curBond->setBondType(RDKit::Bond::SINGLE);
                    res.getAtomWithIdx(*it)->setIsAromatic(false);
                }
            }
            ++count;
        }
    }
}

void fixPhosphateSubstructureAndCharge(RDKit::RWMol &res, unsigned int atIdx) {
    unsigned int noOArNeighbors = 0;
    RDKit::ROMol::ADJ_ITER nbrIdxIt, nbrEndIdxIt;
    vector<unsigned int> toModIdx {};
    boost::tie(nbrIdxIt, nbrEndIdxIt) =
        res.getAtomNeighbors(res.getAtomWithIdx(atIdx));
    while (nbrIdxIt != nbrEndIdxIt) {
        RDKit::Bond *curBond = res.getBondBetweenAtoms(atIdx, *nbrIdxIt);
        if (res.getAtomWithIdx(*nbrIdxIt)->getAtomicNum() == 8 &&
            (curBond->getBondType() == RDKit::Bond::DOUBLE ||
            curBond->getBondType() == RDKit::Bond::AROMATIC)) {
            ++noOArNeighbors;
            toModIdx.push_back(*nbrIdxIt);
        }
        ++nbrIdxIt;
    }
    if (noOArNeighbors > 1) { // I might need to change it to consider non-terminal phosphates
        unsigned int count = 0;
        RDKit::Atom *nbr = res.getAtomWithIdx(atIdx);
        for (vector<unsigned int>::iterator it = toModIdx.begin(); it < toModIdx.end(); ++it) {
            RDKit::Bond *curBond = res.getBondBetweenAtoms(atIdx, *it);
            if (count == (noOArNeighbors - 1)) {
                curBond->setBondType(RDKit::Bond::SINGLE);
                res.getAtomWithIdx(atIdx)->setFormalCharge(0);
                res.getAtomWithIdx(atIdx)->setIsAromatic(false);
                res.getAtomWithIdx(*it)->setFormalCharge(-1);
                res.getAtomWithIdx(*it)->setIsAromatic(false);
                nbr->setIsAromatic(false);
            } else {
                curBond->setBondType(RDKit::Bond::DOUBLE);
                res.getAtomWithIdx(*it)->setIsAromatic(false);
            }
            ++count;
        }
    }
}

void guessFormalCharges(RDKit::RWMol &res) {
  for (RDKit::RWMol::AtomIterator atomIt = res.beginAtoms();
       atomIt != res.endAtoms(); ++atomIt) {
    RDKit::Atom *at = (*atomIt);

    if (at->getFormalCharge() == 0 && at->getSymbol() != "C" &&
        !(at->hasQuery())) {
      int noAromBonds = 0;
      double accum = 0;
      RDKit::ROMol::OEDGE_ITER beg, end;
      boost::tie(beg, end) = res.getAtomBonds(at);
      while (beg != end) {
        accum += res[*beg]->getValenceContrib(at); //(*res)[*beg]
        if (res[*beg]->getBondType() == RDKit::Bond::AROMATIC) {
          ++noAromBonds;
        }
        ++beg;
      }


      // Assumption: if there is an aromatic bridge atom the accum will be 4.5
      //(three aromatic bonds), e.g. naphthalenes. However those are not charged
      // so we can stop here
      if (noAromBonds > 2 && at->getSymbol() == "C") {
        continue;
      }

      // dbtranslate problems - for mols with no UNITY_ATOM_ATTR set - we won't
      // be able to guess
      // if this mol needs to be guessed or not ... hence, we don't guess on
      // atoms with
      // ar specification and not atom type X.ar and in 5-membered ring (there
      // is stuff like
      // c1cc[o+]cc1 that should be charged
      // e.g. string with N.pl3 as NH atom or other atoms without ar
      // specification in aromatic ring
      // FIX: do we need make sure this only happens for atoms in ring?
      std::string tATT;
      at->getProp(RDKit::common_properties::_TriposAtomType, tATT);
      RDKit::MolOps::findSSSR(res); //(*res) if method is passed a pointer
      if (tATT.find("ar") == std::string::npos && at->getIsAromatic() &&
          res.getRingInfo()->isAtomInRingOfSize(at->getIdx(), 5)) {
        continue;
      }

      // for dbtranslate a problem will also occur for N.ar with 3 aromatic
      // bonds and no charge
      // assigned for these we don't assign charges (corina will create
      // kekulized input for this
      //(at least in most cases) - anyway, throw a warning!
      if (noAromBonds == 3 && tATT == "N.ar") {
        std::string nm;
        res.getProp(RDKit::common_properties::_Name, nm);
        BOOST_LOG(rdWarningLog)
            << nm
            << ": warning - aromatic N with 3 aromatic bonds - "
               "skipping charge guess for this atom"
            << std::endl;
        continue;
      }

      // sometimes things like benzimidazoles can have only one bond of the
      // imidazole ring as aromatic and the other one as a single bond ...
      // catch that this way - see also the trick from GL
      int expVal = static_cast<int>(round(accum + 0.1));
      const RDKit::INT_VECT &valens =
          RDKit::PeriodicTable::getTable()->getValenceList(at->getAtomicNum());
      RDKit::INT_VECT_CI vi;

      // check default valence and compare to expVal - chg
      // the hypothesis is that we prefer positively charged atoms over
      // negatively charged ones
      // for multi default valence atoms (e.g. CS(O)(O) should end up being
      // C[S+]([O-])[O-] rather
      // than C[S-][O-][O-] but that might change based no different examples
      int nElectrons =
          RDKit::PeriodicTable::getTable()->getNouterElecs(at->getAtomicNum());
      int assignChg;
      if (nElectrons >= 4)
        assignChg = expVal - (*valens.begin());
      else
        assignChg = (*valens.begin()) - expVal;
      if (assignChg > 0 && nElectrons >= 4) {
        for (vi = valens.begin(); vi != valens.end(); ++vi) {
          // Since we do this only for nocharged atoms we can get away without
          // including the
          // charge into this
          // apart from that we do not assign charges higher than +/- 1 for
          // atoms with multiple valence states
          // otherwise the early break would have to go away which in turn would
          // result in things like [S+4] for
          // sulfonamides
          assignChg = expVal - (*vi);
          //// DEBUG LINES BELOW
          //// They will help you dealing with valence and formal charges.
          //cout << "Assigned charge of atom Z = " << at->getAtomicNum();
          //cout << " equals to " << assignChg << "\n" << endl; //debug line
          if ((*vi) <= expVal && abs(assignChg) < 2) {
            break;
          }
        }
      }
      if (assignChg) {
        // no aromatic atom will get aa abs(charge) > 1
        if (at->getIsAromatic() && abs(assignChg) > 1) {
          at->setFormalCharge((assignChg > 0) -
                              (assignChg < 0));  // this results in -1 or +1
        } else {
          at->setFormalCharge(assignChg);
        }
        // corina will create strange nitro groups which look like N(=O)(=O)
        // which in turn will result in
        //[N+2](=O)(=O) this needs to be fixed at this stage since otherwise we
        // would have to check on all
        // N.pl3 during the cleanup substructures step
        if (assignChg == 2 && expVal == 5 && at->getSymbol() == "N") {
          fixNitroSubstructureAndCharge(res, at->getIdx());
        } else if (assignChg > 0 && tATT == "P.3") {
          at->setFormalCharge(0);
        }
        // what if expVal != at->calcExplicitValence();?
        // cannot imagine a case where that will happen now.
      }
    }
  }
}

unsigned int chkNoHNeighbNOx(RDKit::RWMol &res, RDKit::ROMol::ADJ_ITER atIdxIt,
                             int &toModIdx) {
  RDKit::Atom *at = res.getAtomWithIdx(*atIdxIt);
  unsigned int noHNbrs = 0;
  RDKit::ROMol::ADJ_ITER nbrIdxIt, nbrEndIdxIt;
  boost::tie(nbrIdxIt, nbrEndIdxIt) = res.getAtomNeighbors(at);
  while (nbrIdxIt != nbrEndIdxIt) {
    if (res.getAtomWithIdx(*nbrIdxIt)->getAtomicNum() == 1) {
      ++noHNbrs;
    } else if (res.getAtomWithIdx(*nbrIdxIt)->getAtomicNum() == 8 &&
               res.getAtomDegree(res.getAtomWithIdx(*nbrIdxIt)) == 1) {
      // this is a N in an N-oxide constellation
      // we can do the above if clause since mol2 have explicit hydrogens
      toModIdx = *atIdxIt;
    }
    ++nbrIdxIt;
  }
  return noHNbrs;
}

bool cleanUpSubstructures(RDKit::RWMol &res) {
    // NOTE: check the nitro fix in guess formal charges!
    boost::dynamic_bitset<> isFixed(res.getNumAtoms());
    for (RDKit::ROMol::AtomIterator atIt = res.beginAtoms(); atIt != res.endAtoms();
        ++atIt) {
    std::string tAT;
    RDKit::Atom *at = *atIt;
    unsigned int idx = at->getIdx();
    at->getProp(RDKit::common_properties::_TriposAtomType, tAT);

    if (tAT == "N.4") {
        at->setFormalCharge(1);
    } else if (tAT == "N.2") { // sulfonamides can be tricky
        RDKit::ROMol::ADJ_ITER nbrIdxIt, endNbrsIdxIt;
        boost::tie(nbrIdxIt, endNbrsIdxIt) = res.getAtomNeighbors(at);
        RDKit::Atom *nbr = res.getAtomWithIdx(*nbrIdxIt);
        string tATT;
        nbr->getProp(RDKit::common_properties::_TriposAtomType, tATT);
        // check if neighbor is S.3 and bond order is 2
        if (tATT == "S.3") {
            RDKit::Bond *b =  res.getBondBetweenAtoms(idx, *nbrIdxIt);
            if (b->getBondType() == RDKit::Bond::DOUBLE) {
                b->setBondType(RDKit::Bond::SINGLE);
            }
        }
    } else if (tAT == "N.pl3") {
        fixNitroSubstructureAndCharge(res, idx);
    } else if (tAT == "C.2") { // Discover amidines
        fixProtonatedAmidineSubstructureAndCharge(res, idx);
    } else if (tAT == "O.2") {
        // Sulfonamines frequently have S.3 and O.2 instead of S.o2 and O.co2 atom types.
        if (at->getDegree() != 1) { // O.2 can only have one neighbor
            BOOST_LOG(rdWarningLog) << "Warning - O.2 with degree >1." << endl;
            return false;
        }
        RDKit::ROMol::ADJ_ITER nbrIdxIt, endNbrsIdxIt;
        boost::tie(nbrIdxIt, endNbrsIdxIt) = res.getAtomNeighbors(at);
        RDKit::Atom *nbr = res.getAtomWithIdx(*nbrIdxIt); // We need to get the neighbor's info
        string tATT;
        nbr->getProp(RDKit::common_properties::_TriposAtomType, tATT);
        // bizzarely typed sulfonamide
        if (tATT == "S.3") {
            // We have to make sure that bond orders surrounding this sulfur atom are correct 
            RDKit::Bond *b =  res.getBondBetweenAtoms(idx, *nbrIdxIt);
            b->setBondType(RDKit::Bond::DOUBLE);
            b->setIsAromatic(false);
            at->setIsAromatic(false);
            nbr->setIsAromatic(false);
        }
    } else if (tAT == "O.co2") {
        // negatively charged carboxylates with O.co2
        // according to Tripos, those should only appear in carboxylates and
        // phosphates,
        if (at->getDegree() != 1) { // degrees are the #neighbors in the graph
            BOOST_LOG(rdWarningLog) << "Warning - O.co2 with degree >1." << endl;
            return false;
        }
        RDKit::ROMol::ADJ_ITER nbrIdxIt, endNbrsIdxIt;
        // getAtomNeighbors returns 2 adjacency iterators (index iterators)
        boost::tie(nbrIdxIt, endNbrsIdxIt) = res.getAtomNeighbors(at);
        // this should return only the C.2
        RDKit::Atom *nbr = res.getAtomWithIdx(*nbrIdxIt);
        string tATT;
        nbr->getProp(RDKit::common_properties::_TriposAtomType, tATT);
        // carboxylates
        if (tATT == "C.2" || tATT == "S.o2") {
        //if (tATT == "C.2" || tATT == "S.o2") {
            // this should return only the bond between C.2 and O.co2
            RDKit::Bond *b = res.getBondBetweenAtoms(idx, *nbrIdxIt);
            if (!isFixed[*nbrIdxIt]) {
                // the first occurrence is negatively charged and has a single bond
                b->setBondType(RDKit::Bond::SINGLE);
                b->setIsAromatic(false);
                at->setFormalCharge(-1);
                at->setIsAromatic(false);
                nbr->setIsAromatic(false);
                isFixed[idx] = 1;
                isFixed[*nbrIdxIt] = 1;
            } else {
                // the other occurrences are not charged and have a double bond
                b->setBondType(RDKit::Bond::DOUBLE);
                b->setIsAromatic(false);
                at->setIsAromatic(false);
                isFixed[idx] = 1;
            }
        // GDRM 2020-11-03 patch
        } else if (tATT == "P.3") {
            fixPhosphateSubstructureAndCharge(res, *nbrIdxIt);
        } else {
            std::string nm;
            res.getProp(RDKit::common_properties::_Name, nm);
            BOOST_LOG(rdWarningLog)
            << nm << ": warning - O.co2 with non C.2 or S.o2 neighbor."
            << std::endl;
            return false;
        }
    } else if (tAT == "C.cat") {
        // positively charged guanidinium groups with C.cat
        // according to Tripos these should only appear in guanidinium groups
        // for the structural fix - the last nitrogen with the least number of
        // heavy atoms will get the double bond and the positive charge.
        // remember : this is not canonical!
        // first - set the C.cat as fixed
        isFixed[idx] = 1;
        RDKit::ROMol::ADJ_ITER nbrIdxIt, endNbrsIdxIt, tmpIdxIt;
        unsigned int lowestDeg = 100;
        boost::tie(nbrIdxIt, endNbrsIdxIt) = res.getAtomNeighbors(at);
        // one problem of programs like Corina is, that they will create also
        // C.cat
        // for groups that are not guanidinium. We cannot fix all, but the charged
        // amidine
        // in a ring is taken care of too.
        tmpIdxIt = nbrIdxIt;
        // declare and initialise toModIdx
        int toModIdx = -1;
        unsigned int noNNeighbors = 0;
        while (tmpIdxIt != endNbrsIdxIt) {
            if (res.getAtomWithIdx(*tmpIdxIt)->getSymbol() == "N") {
                ++noNNeighbors;
            }
            ++tmpIdxIt;
        }
        if (noNNeighbors < 2 || noNNeighbors > 3) {
            std::string nm;
            res.getProp(RDKit::common_properties::_Name, nm);
            BOOST_LOG(rdWarningLog)
            << nm << ": Error - C.Cat with bad number of N neighbors."
            << std::endl;
            return false;
        } else if (noNNeighbors == 2) {
            // the idea is that we assign the positive charge according to the
            // following precedence:
            // 1. is part of N-oxide
            // 2. atom with highest number of hydrogen atoms
            // 3. atom in ring
            // 4. random
            // first we identify the N atoms
            RDKit::ROMol::ADJ_ITER idxIt1 = nbrIdxIt, idxIt2 = nbrIdxIt;
            bool firstIdent = false;
            while (nbrIdxIt != endNbrsIdxIt) {
                if (res.getAtomWithIdx(*nbrIdxIt)->getSymbol() == "N") {
                    // fix the bond to one - only the modified N will have a double bond
                    // to C.cat
                    res.getBondBetweenAtoms(idx, *nbrIdxIt)->setBondType(RDKit::Bond::SINGLE);
                    res.getBondBetweenAtoms(idx, *nbrIdxIt)->setIsAromatic(false);
                    res.getAtomWithIdx(*nbrIdxIt)->setIsAromatic(false);
                    // FIX: what is happening if we hit an atom that was fixed before -
                    // propably nothing.
                    // since I cannot think of a case where this is a problem - throw a
                    // warning
                    if (isFixed[*nbrIdxIt]) {
                        std::string nm;
                        res.getProp(RDKit::common_properties::_Name, nm);
                        BOOST_LOG(rdWarningLog)
                        << nm << ": warning - charged amidine and isFixed atom."
                        << std::endl;
                    }
                    isFixed[*nbrIdxIt] = 1;
                    if (firstIdent) {
                        idxIt2 = nbrIdxIt;
                    } else {
                        idxIt1 = nbrIdxIt;
                        firstIdent = true;
                    }
                }
                ++nbrIdxIt;
            }
            // now that we know which are the relevant atoms we check the above
            // features
            // is part of N-oxide?
            // number of hydrogens on each neighbour
            unsigned int noHNbrs1 = chkNoHNeighbNOx(res, idxIt1, toModIdx);
            unsigned int noHNbrs2 = chkNoHNeighbNOx(res, idxIt2, toModIdx);
            if (toModIdx < 0) {
                // no N-oxide
                if (noHNbrs1 != noHNbrs2) {
                    if (noHNbrs1 > noHNbrs2) {
                        toModIdx = *idxIt1;
                    } else {
                        toModIdx = *idxIt2; // this is random if both have the same
                                            // number of atoms
                            }
                    } else {
                        // perceive the rings
                        RDKit::MolOps::findSSSR(res); //(*res) if passed a pointer
                        // then we check if both atoms are in a ring
                        unsigned int rIdx1 = res.getRingInfo()->numAtomRings((*idxIt1));
                        unsigned int rIdx2 = res.getRingInfo()->numAtomRings((*idxIt2));
                        if (rIdx1 > rIdx2) {
                            toModIdx = *idxIt1;
                        } else {
                            toModIdx = *idxIt2;
                        }
                    }
                }
                res.getBondBetweenAtoms(idx, toModIdx)->setBondType(RDKit::Bond::DOUBLE);
                res.getBondBetweenAtoms(idx, toModIdx)->setIsAromatic(false);
                res.getAtomWithIdx(toModIdx)->setFormalCharge(1);
                at->setIsAromatic(false);
            } else {
                while (nbrIdxIt != endNbrsIdxIt) {
                    if (!isFixed[*nbrIdxIt]) {
                        // we get in here if this N.pl3 was not seen / fixed before
                        RDKit::Atom *nbr = res.getAtomWithIdx(*nbrIdxIt);
                        // get the number of heavy atoms connected to this atom
                        RDKit::ROMol::ADJ_ITER nbrNbrIdxIt, nbrEndNbrsIdxIt;
                        unsigned int hvyAtDeg = 0;
                        boost::tie(nbrNbrIdxIt, nbrEndNbrsIdxIt) =
                        res.getAtomNeighbors(nbr);
                        while (nbrNbrIdxIt != nbrEndNbrsIdxIt) {
                            if (res.getAtomWithIdx(*nbrNbrIdxIt)->getAtomicNum() > 1) {
                                std::string nbrAT;
                                res.getAtomWithIdx(*nbrNbrIdxIt)
                                ->getProp(RDKit::common_properties::_TriposAtomType, nbrAT);
                                if (nbrAT == "C.cat") {
                                    hvyAtDeg += 2;  // that way we reduce the risk of ionising the
                                                    // N attached to another C.cat ...
                                } else {
                                    ++hvyAtDeg;
                                }
                            }
                            ++nbrNbrIdxIt;
                        }
                        // now check for lowest heavy atom degree
                        if (hvyAtDeg < lowestDeg) {
                            toModIdx = *nbrIdxIt;
                            lowestDeg = hvyAtDeg;
                        }
                        // modify the bond between C.Cat and the N.pl3
                        RDKit::Bond *b = res.getBondBetweenAtoms(idx, *nbrIdxIt);
                        b->setBondType(RDKit::Bond::SINGLE);
                        b->setIsAromatic(false);
                        nbr->setIsAromatic(false);
                        // set N.pl3 as fixed
                        isFixed[*nbrIdxIt] = 1;
                    } else {
                        // the N is allready fixed - since we don't touch this atom make the
                        // bond to single
                        // FIX: check on 3-way symmetric guanidinium mol -
                        //     this could produce a only single bonded C.cat for bad H mols
                        res.getBondBetweenAtoms(idx, *nbrIdxIt)->setBondType(RDKit::Bond::SINGLE);
                        res.getBondBetweenAtoms(idx, *nbrIdxIt)->setIsAromatic(false);
                    }
                    ++nbrIdxIt;
                }
                // now modify the respective N and the C.cat
                RDKit::Bond *b = res.getBondBetweenAtoms(idx, toModIdx);
                b->setBondType(RDKit::Bond::DOUBLE);
                b->setIsAromatic(false);
                res.getAtomWithIdx(toModIdx)->setFormalCharge(1);
                at->setIsAromatic(false);
            }
        }
        idx++;
    }
    return true;
}

/////////

bool DOCKMol::isNitro( string atom_type, int idx ){
    /*
    We need to identify nitro groups if we are to have them properly read
    by RDKit from our DOCKMol objects.
    */
    if (atom_type == "N.pl3"){
        int numONeighbors = 0;
        vector<int> neighbor_atoms = get_atom_neighbors(idx);
        for (vector<int>::iterator ptr = neighbor_atoms.begin();
            ptr < neighbor_atoms.end(); ++ptr) {
            if (atom_types[*ptr] == "O.2") { ++numONeighbors; }
        }
        if (numONeighbors == 2) { //Nitros have 2 O.2 atoms bonded to N.pl3
            return true;
        } else { return false; }
    } else { return false; }
}

// bool DOCKMol::isPhosphate( string atom_type, int idx ){
//     /*
//     We need to identify phosphate groups if we are to have them properly read
//     by RDKit from our DOCKMol objects.
//     */
//     if (atom_type == "P.3"){
//         int numONeighbors = 0;
//         vector<int> neighbor_atoms = get_atom_neighbors(idx);
//         for (vector<int>::iterator ptr = neighbor_atoms.begin();
//             ptr < neighbor_atoms.end(); ++ptr) {
//             if (atom_types[*ptr] == "O.co2") { ++numONeighbors; }
//         }
//         if (numONeighbors == 3) { //Phosphates have 3 O.co2 atoms bonded to P.3
//             return true;
//         } else {
//             return false;
//         }
//     } else {
//         return false;
//     }
// }

// June 01, 2020
RDKit::Atom *DOCKMol::createAtom( unsigned int atomIdx, RDGeom::Point3D & pos, bool noImplicit=true ){
    auto *rdatom = new RDKit::Atom();

    // Add positions
    pos.x = x[atomIdx];
    pos.y = y[atomIdx];
    pos.z = z[atomIdx];

    // Figure out and assign atomic number
    int atomic_num = find_atomic_number( atom_types[atomIdx] );
    rdatom->setAtomicNum( atomic_num );
    bool npl3 = isNitro( atom_types[atomIdx], atomIdx );

    // Set properties
    rdatom->setProp( RDKit::common_properties::_TriposAtomType, atom_types[atomIdx] );
    rdatom->setProp( "atomLabel", atom_names[atomIdx] );
    rdatom->setProp( "_TriposPartialCharge", charges[atomIdx] );
    rdatom->setNoImplicit( noImplicit ); // no implicit atoms
    if (npl3) {
        rdatom->setFormalCharge(2); // Future RDKit-related methods fix this.
    }

    return rdatom;
}


// June 01, 2020
RDKit::Bond *DOCKMol::createBond( unsigned int bondIdx ){

    // // Debug
    // ofstream outfile;
    // outfile.open( "DEBUG_CREATE_BOND.txt", ios::app );

    // Figure out bond type
    // outfile << "Select bond type.\n";
    RDKit::Bond::BondType bond_type = select_bond_type( bond_types[bondIdx] );

    //// Remove aromatic bonds in phosphate and nitro groups
    //// Identify atom identities and get index of N.pl3 (or P.3)
    //string tAT1 { atom_types[ bonds_origin_atom[bondIdx] ] };
    //string tAT2 { atom_types[ bonds_target_atom[bondIdx] ] };
    //if (tAT1 == "N.pl3") {
    //    int atIdx = bonds_origin_atom[bondIdx];
    //} else if (tAT2 == "N.pl3") {
    //    int atIdx = bonds_target_atom[bondIdx];
    //} else if (tAT1 == "P.3") {
    //    int atIdx = bonds_origin_atom[bondIdx];
    //} else if (tAT2 == "P.3") {
    //    int atIdx = bonds_target_atom[bondIdx];
    //} else {
    //}

    // // Debug
    // outfile << "Atom types of interest: " << tAT1 << " and " << tAT2 << endl;
    //
    // if (bond_type == RDKit::Bond::AROMATIC && (tAT1 == "N.pl3" || tAT2 == "N.pl3")){
    //
    //     // Debug
    //     outfile << "Potentially problematic N=O bond between atoms ";
    //     outfile << bonds_origin_atom[bondIdx] << " and ";
    //     outfile << bonds_target_atom[bondIdx] << endl;
    //
    //     bond_type = RDKit::Bond::DOUBLE;
    //}

    // Create bond and assign its ends
    auto *bond = new RDKit::Bond( bond_type );
    bond->setBeginAtomIdx( bonds_origin_atom[bondIdx] );
    bond->setEndAtomIdx( bonds_target_atom[bondIdx] );

    // Debug
    //outfile.close();

    return bond;
}

// June 01, 2020
void DOCKMol::addAtomsToRWMol( RDKit::RWMol & tmp ){
    // Create uninitialized set of coordinates
    vector<RDGeom::Point3D> allCoords;
    allCoords.reserve(num_atoms);

    //// Debug
    //ofstream outfile;
    //outfile.open("DEBUG_ATOM.txt");

    // Loop over all atoms and add them to mol object
    for (auto atomIdx = 0; atomIdx < num_atoms; ++atomIdx){
        // Create atom
        RDGeom::Point3D pos;
        //outfile << "Creating atom idx = " << atomIdx << endl;
        RDKit::Atom *rdatom = createAtom( atomIdx, pos );
        //outfile << "Atomic number: " << rdatom->getAtomicNum() << endl;
        tmp.updatePropertyCache();
        // Add atom to mol object
        //int idx = tmp.addAtom(rdatom, false, true); // idx can be useful for debugging
        tmp.addAtom(rdatom, false, true);
        allCoords.push_back(pos);
    }

    // Create conformer based on allCoords
    // and add to RWMol
    //// Debug
    //outfile << "Create conformation" << endl;

    auto *conf = new RDKit::Conformer(allCoords.size());
    vector<RDGeom::Point3D>::const_iterator allCoordsIt = allCoords.begin();
    for (unsigned int i = 0; i < allCoords.size(); ++i){
        conf->setAtomPos(i, *allCoordsIt);
        ++allCoordsIt;
    }

    //// Debug
    //outfile << "Add conformer" << endl;

    tmp.addConformer(conf, true);

    // Debug
    //outfile << "Closing debug file." << endl;
    //outfile.close();
}

// June 01, 2020
void DOCKMol::addBondsToRWMol( RDKit::RWMol & tmp ){

    //// Debug
    //ofstream outfile;
    //outfile.open("DEBUG_BONDS.txt");

    // Loop over all bond indices and create bonds
    for (auto bondIdx = 0; bondIdx < num_bonds; ++bondIdx){

        //// Debug
        //outfile << "Create bond " << bondIdx << endl;

        RDKit::Bond *bond = createBond( bondIdx );
        if (bond->getBondType() == RDKit::Bond::AROMATIC){

            //// Debug
            //outfile << "Bond between atoms " << bond->getBeginAtomIdx();
            //outfile << " and " << bond->getEndAtomIdx() << " is aromatic.\n";

            bond->setIsAromatic(true);
            tmp.getAtomWithIdx(bond->getBeginAtomIdx())->setIsAromatic(true);
            tmp.getAtomWithIdx(bond->getBeginAtomIdx())->updatePropertyCache();
            tmp.getAtomWithIdx(bond->getEndAtomIdx())->setIsAromatic(true);
            tmp.getAtomWithIdx(bond->getEndAtomIdx())->updatePropertyCache();
        } else {
            tmp.getAtomWithIdx(bond->getBeginAtomIdx())->updatePropertyCache();
            tmp.getAtomWithIdx(bond->getEndAtomIdx())->updatePropertyCache();
        }
        tmp.addBond(bond, true);
    }

    //// Debug
    //outfile.close();
}

// GDRM June 01, 2020
RDKit::ROMol  DOCKMol::DOCKMol_to_ROMol( bool create_smiles ){
    // create temporary RWMol object
    RDKit::RWMol tmp;

    // assign parameters from DOCKMol to RWMol
    // odbfile << "Adding atoms and bonds\n";
    addAtomsToRWMol( tmp );
    // odbfile << "Atoms added to RWMol\n";
    addBondsToRWMol( tmp );
    // odbfile << "Bonds added to RWMol\n";

    // Make molecule Lewis structure-compliant
    bool fix = cleanUpSubstructures( tmp );
    // odbfile << "Clean up done!\n";
    if (fix){
        guessFormalCharges( tmp );
        // odbfile << "Charging done!\n";
    } else {
        cerr << "Error: " << title << " has issues." << endl;
        auto tmp = RDKit::RWMol();
        return tmp;
    }

    // assign chirality
    RDKit::MolOps::assignChiralTypesFrom3D( tmp );

    // clean up
    RDKit::MolOps::cleanUp( tmp );

    // sanitize
    // odbfile << "Starting sanitization\n";

    RDKit::MolOps::sanitizeMol( tmp );
    // odbfile << "Sanitization done\n";

    // detect bond stereochemistry
    RDKit::MolOps::detectBondStereochemistry( tmp );

    // Update
    tmp.updatePropertyCache( false );

    // assign stereochemistry
    RDKit::MolOps::assignStereochemistryFrom3D( tmp );

    if (create_smiles) {
        RDKit::RWMol smi_gen {tmp};
        RDKit::MolOps::removeHs(smi_gen);
        smiles = RDKit::MolToSmiles(smi_gen, true, false, -1, true, false, false, false);
    }

    // Copy all assigned attributes to an immutable RDKit::ROMol object.
    RDKit::ROMol rdmol{tmp}; // The RDKit::ROMol object cannot be modified.

    return rdmol;
}

//void DOCKMol::createRDKitDescriptors( bool create_smiles ){
//    // Create temporary rdmol object
//    // When given "true" it creates and assigns SMILES attribute.
//    try {
//        RDKit::ROMol rdmol = DOCKMol_to_ROMol( create_smiles );
//
//        // Create descriptors
//        num_arom_rings = RDKit::Descriptors::calcNumAromaticRings(rdmol);
//        num_alip_rings = RDKit::Descriptors::calcNumAliphaticRings(rdmol);
//        num_sat_rings = RDKit::Descriptors::calcNumSaturatedRings(rdmol);
//        num_stereocenters = RDKit::Descriptors::numAtomStereoCenters(rdmol);
//        num_spiro_atoms = RDKit::Descriptors::calcNumSpiroAtoms(rdmol);
//        clogp = RDKit::Descriptors::calcClogP(rdmol);
//        tpsa = RDKit::Descriptors::calcTPSA(rdmol);
//    } catch (...) {
//        //RDKit::ROMol rdmol {};
//        //cout << "RDKit couldn't interpret molecule ";
//        //cout << title << ". Investigate structure for details." << endl;
//        smiles = title + ": RDKIT ERROR";
//    }
//
//}

#endif

/*********************************************/
int
less_than_pair(SCOREMol a, SCOREMol b)
{
    return (a.first < b.first);
}

/************************************************/
int
BREADTH_SEARCH::get_search_radius(DOCKMol & mol, int root_atom, int avoid_atom)
{
    int             root,
                    max_radius;
    int             nbr_atom;
    vector < int   >atom_nbrs;

    atoms.clear();
    nbrs.clear();
    nbrs_next.clear();

    atoms.resize(mol.num_atoms, -1);
    nbrs.push_back(root_atom);
    atoms[root_atom] = 0;
    atoms[avoid_atom] = -2;

    while (nbrs.size() > 0) {
        root = nbrs[nbrs.size() - 1];
        nbrs.pop_back();

        // atom_nbrs = mol.get_atom_neighbors(root);
        atom_nbrs = mol.neighbor_list[root];

        for (long unsigned int i = 0; i < atom_nbrs.size(); i++) {
            nbr_atom = atom_nbrs[i];

            if (atoms[nbr_atom] == -1) {
                atoms[nbr_atom] = atoms[root] + 1;
                nbrs.push_back(nbr_atom);
            }
        }

        if ((nbrs.size() == 0) && (nbrs_next.size() > 0)) {
            nbrs = nbrs_next;
            nbrs_next.clear();
        }

    }

    max_radius = 0;
    for (long unsigned int x = 0; x < atoms.size(); x++) {
        if (atoms[x] > max_radius)
            max_radius = atoms[x];
    }

    return max_radius;
}

/******************************************************/
void
transform(DOCKMol & mol, float rmat[3][3], DOCKVector trans, DOCKVector com)
{
    int             i;
    float           temp1,
                    temp2,
                    temp3,
                    temp4,
                    temp5,
                    temp6;

    for (i = 0; i < mol.num_atoms; i++) {
        temp1 = mol.x[i] - com.x;
        temp2 = mol.y[i] - com.y;
        temp3 = mol.z[i] - com.z;

        temp4 = rmat[0][0] * temp1;
        temp5 = rmat[1][0] * temp2;
        temp6 = rmat[2][0] * temp3;

        mol.x[i] = temp4 + temp5 + temp6 + com.x + trans.x;

        temp4 = rmat[0][1] * temp1;
        temp5 = rmat[1][1] * temp2;
        temp6 = rmat[2][1] * temp3;

        mol.y[i] = temp4 + temp5 + temp6 + com.y + trans.y;

        temp4 = rmat[0][2] * temp1;
        temp5 = rmat[1][2] * temp2;
        temp6 = rmat[2][2] * temp3;

        mol.z[i] = temp4 + temp5 + temp6 + com.z + trans.z;

    }

}

// +++++++++++++++++++++++++++++++++++++++++
// Copy the Cartesian coordinates from the argument to this DOCKMol.
// The linear memory representation of argument xyz, starting from *xyz, is
// first atom x, first atom y, first atom z, second atom x, ...
void
DOCKMol::setxyz( const double * xyz )
{
    for (int i = 0; i < num_atoms; ++ i) {
        x[i] = xyz[ 3 * i + 0 ];
        y[i] = xyz[ 3 * i + 1 ];
        z[i] = xyz[ 3 * i + 2 ];
    }
}

// +++++++++++++++++++++++++++++++++++++++++
// Copy the charges from the argument to this DOCKMol
// multiplying them by the units conversion factor.
void
DOCKMol::setcharges(const double* chrgs, double units_factor )
{
    for (int i = 0; i < num_atoms; ++ i) {
        charges[i] = units_factor * chrgs[i];
    }
}

FOOTPRINT_ELEMENT::FOOTPRINT_ELEMENT(){
   resname = " ";
   resid   = 0;
   vdw     = 0;
   es      = 0;
   hb      = 0;
}

PH4_ELEMENT::PH4_ELEMENT(){
   ph4_type   = " ";
   ph4_id     = 0;
   x_crd      = 0.0;
   y_crd      = 0.0;
   z_crd      = 0.0;
   v_x        = 0.0;
   v_y        = 0.0;
   v_z        = 0.0;
   ph4_radius = 0.0;
}

/*********************************************/
HDB_atom::HDB_atom(){
}
/*********************************************/
HDB_atom::~HDB_atom(){
}
//void HDB_atom::initialize(int,    std::string,  std::string, int,    int,     float,     float,   float,    float)
void HDB_atom::initialize(int anum, string aname, string ast, int asn, int adt, float crg, float p, float ap, float sa){
//void HDB_atom::initialize(int anum, char aname[], char ast[], int asn, int adt, float crg, float p, float ap, float sa){
  atom_num       = anum;
  //atom_name      = string(aname);
  //atom_syb_type  = string(ast);
  atom_name      = aname;
  atom_syb_type  = ast;
  atom_seg_num   = asn;
  atom_dock_type = adt;
  q              = crg;
  p_desolv       = p;
  a_desolv       = ap;
  sasa           = sa;
}
void HDB_atom::print_atom(){
  cout << "atom_num       " <<  atom_num       << endl; 
  cout << "atom_name      " <<  atom_name      << endl; 
  cout << "atom_syb_type  " <<  atom_syb_type  << endl; 
  cout << "atom_seg_num   " <<  atom_seg_num   << endl; 
  cout << "atom_dock_type " <<  atom_dock_type << endl; 
  cout << "q              " <<  q              << endl; 
  cout << "p_desolv       " <<  p_desolv       << endl; 
  cout << "a_desolv       " <<  a_desolv       << endl; 
  cout << "sasa           " <<  sasa           << endl; 
}
/*********************************************/
HDB_bond::HDB_bond(){
}
/*********************************************/
HDB_bond::~HDB_bond(){
}
/*********************************************/
void HDB_bond::initialize(int b, int a1, int a2, string t){
  bond_num  = b;
  atom1_num = a1;
  atom2_num = a2;
  bond_name = t; // sybyl bond type
}

/*********************************************/
HDB_coordinate::HDB_coordinate(){
}
/*********************************************/
HDB_coordinate::~HDB_coordinate(){
}
/*********************************************/
void HDB_coordinate::initialize(int c, int a, int s, float tx, float ty, float tz){
  coord_num = c;
  atom_num  = a;
  seg_num   = s;
  x         = tx;
  y         = ty;
  z         = tz;
}
void HDB_coordinate::print_coord(){
  cout << "coord_num " << coord_num << endl; 
  cout << "atom_num  " << atom_num  << endl;
  cout << "seg_num   " << seg_num   << endl;
  cout << "x         " << x         << endl;
  cout << "y         " << y         << endl;
  cout << "z         " << z         << endl;
}
/*********************************************/
HDB_segment::HDB_segment(){
}
/*********************************************/
HDB_segment::~HDB_segment(){
}
/*********************************************/
void HDB_segment::initialize(int num, int start, int stop){
  seg_num     = num;
  start_coord = start;
  stop_coord  = stop;
}
/*********************************************/
HDB_conformer::HDB_conformer(){
  list_of_seg = NULL;
}
/*********************************************/
HDB_conformer::~HDB_conformer(){
  delete[]list_of_seg;
}
/*********************************************/
//void HDB_conformer::initialize(int num, int size, float ie, int * sl){
void HDB_conformer::initialize(int num, int size, float ie){
    conf_num = num;
    num_of_seg = size;
    internal_energy = ie;
    list_of_seg = new int[size];
}

/*********************************************/
HDB_Mol::HDB_Mol(){
name = new char [100];
atoms  = NULL;
bonds  = NULL;
coords = NULL;
rigid  = NULL;
segs   = NULL;
confs  = NULL;
}

/*********************************************/
HDB_Mol::~HDB_Mol(){
  delete[]name;  
  delete[]atoms;
  delete[]bonds;
  delete[]coords;
  delete[]rigid;
  delete[]segs;
  delete[]confs;
}

