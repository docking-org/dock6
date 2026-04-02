#include <iostream>
#include <string.h>
#include "dockmol.h"
#include "fragment.h"
#include "iso_align.h"

using namespace std;

class Iso_Align;


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
               num_du = 0; 
       };
Fragment::~Fragment(){
            mol.clear_molecule();
            aps.clear();
            torenv_recheck_indices.clear();
            frag_growth_tree.clear();
            aps_bonds_type.clear();
            aps_cos.clear();
            aps_dist.clear();
            num_du = 0;
       };


// methods that sets radial distanace for iso_align protocol
void Fragment::set_radial_dist_distri(int atom_num, std::vector<float> radial){

    Iso_Align iso_ali;
    (*this).radial_dist_distri[atom_num] = radial;
}

void Fragment::calc_radial_dist_distri(){

    Iso_Align iso_ali;
    iso_ali.get_all_radial_dist(*this);

}

void Fragment::print_radial_dist_distri(){
    for (unsigned int i; i< this->mol.num_atoms; i++){
        int obj_radial_vec_size = this->radial_dist_distri[i].size();
        for (unsigned int j=0; j<obj_radial_vec_size; j++){
            std::cout << this->radial_dist_distri[i][j] << " " ;
            if ((j+1)==obj_radial_vec_size){
                std::cout << std::endl;
            }
        }
    }

}

std::vector<float> Fragment::get_radial_dist_distri(int atom_num){
    return this->radial_dist_distri[atom_num];
}

void Fragment::allocate_radial_dist_distri(){

    radial_dist_distri    = new      std::vector<float> [this->mol.num_atoms];
}

void Fragment::clear_radial_dist_distri(){

    delete[]radial_dist_distri; 

}

void Fragment::calc_num_du(){
    int counter = 0;
    for (unsigned int i =0; i< this->mol.num_atoms; i++){
        if (this->mol.atom_types[i] == "Du"){counter++;}
    }
    this->num_du = counter; 
}

int Fragment::get_num_du(){

    return this->num_du;
}

void Fragment::is_iso_aligned(){
    this->iso_aligned=true;
}
void Fragment::is_not_iso_aligned(){
    this->iso_aligned=false;
}
bool Fragment::is_it_iso_aligned(){

    return this->iso_aligned;
}

// 



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
