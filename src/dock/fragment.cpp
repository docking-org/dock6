#include <iostream>
#include <string.h>
#include "dockmol.h"
#include "fragment.h"

using namespace std;



// +++++++++++++++++++++++++++++++++++++++++
// Some constructors and destructors

AttPoint::AttPoint(){
               heavy_atom = -1;
               dummy_atom = -1;
       };
AttPoint::~AttPoint(){};


Fragment::Fragment(){
               mol.initialize();
               used = false;
               tanimoto = 0.0;
               last_ap_heavy = -1;
               scaffolds_this_layer = 0;
               size = 0;
               ring_size = 0;
               mut_type = 500; // mutations start at 0, so initalizing beyond
        
       };
Fragment::~Fragment(){
            mol.clear_molecule();
            aps.clear();
            torenv_recheck_indices.clear();
            frag_growth_tree.clear();
            aps_bonds_type.clear();
            aps_cos.clear();
            aps_dist.clear();
       };


// +++++++++++++++++++++++++++++++++++++++++
// Part of the initialize function to prepare the fragment object for a molecule
// explicitly called in conf_gen_ga.cpp, (need to update in conf_gen_dn.cpp)
Fragment
Fragment::read_mol ( Fragment & tmp_frag )
{
    // Create temporary vectors of atom indices
    vector <int> tmp_dummy_atoms;
    vector <int> tmp_heavy_atoms;
    
    // Create temp pair for bond information - CS: 09/19/16
    std::pair <int, std::string> bond_info;

    // Loop over every bond in the dockmol object
    for (int i=0; i<tmp_frag.mol.num_bonds; i++){
    
        // Check if the bond origin is a dummy atom
        if (tmp_frag.mol.atom_types[tmp_frag.mol.bonds_origin_atom[i]].compare("Du") == 0){
    
            // if so, save origin as dummy_atom and target as heavy atom
            tmp_dummy_atoms.push_back(tmp_frag.mol.bonds_origin_atom[i]);
            tmp_heavy_atoms.push_back(tmp_frag.mol.bonds_target_atom[i]);
 
            // Update the bond number and type - CS: 09/19/16
            bond_info.first = i;
            bond_info.second = tmp_frag.mol.bond_types[i];
            tmp_frag.aps_bonds_type.push_back(bond_info);                  
        }
   
        // Check if the bond target is a dummy atom
        if (tmp_frag.mol.atom_types[tmp_frag.mol.bonds_target_atom[i]].compare("Du") == 0){
    
            // if so, save target as dummy atom and origin as heavy atom
            tmp_dummy_atoms.push_back(tmp_frag.mol.bonds_target_atom[i]);
            tmp_heavy_atoms.push_back(tmp_frag.mol.bonds_origin_atom[i]);

            // Update the bond number and type - CS: 09/19/16
            bond_info.first = i;
            bond_info.second = tmp_frag.mol.bond_types[i];
            tmp_frag.aps_bonds_type.push_back(bond_info);                  
        }
    }
    
    // Create temporary attachment point object
    AttPoint tmp_ap;
    // Loop over all of the dummy_atom / heavy_atom pairs
    for (int i=0; i<tmp_dummy_atoms.size(); i++){
    
        // Copy the heavy atom / dummy atom pairs onto the attachment point class
        tmp_ap.heavy_atom = tmp_heavy_atoms[i];
        tmp_ap.dummy_atom = tmp_dummy_atoms[i];
    
        // And add that information to the attachment points (aps) vector of the scaffold
        tmp_frag.aps.push_back(tmp_ap);
    } 
    // Clear temporary vectors
    tmp_dummy_atoms.clear();
    tmp_heavy_atoms.clear();
    
    
    return tmp_frag;
}// Fragment::read_mol() 
