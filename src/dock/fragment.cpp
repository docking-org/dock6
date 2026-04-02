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
               frag_name = "XXX";
       };
AttPoint::~AttPoint(){};

AttPoint::AttPoint(const AttPoint &original){
    (*this).heavy_atom = original.heavy_atom;
    (*this).dummy_atom = original.dummy_atom;
    (*this).frag_name = original.frag_name;
       };

void AttPoint::operator=(const AttPoint &original){
    (*this).heavy_atom = original.heavy_atom;
    (*this).dummy_atom = original.dummy_atom;
    (*this).frag_name = original.frag_name;

       };


Fragment::Fragment(){
            mol.initialize();
            used = false;
            tanimoto = 0.0;
            last_ap_heavy = -1;
            scaffolds_this_layer = 0;
            size = 0;
            ring_size = 0;
            mut_type = 500; // mutations start at 0, so initalizing beyond
            iso_score = 9999;
            radial_set = false;
            alloc_set = false;
            num_du = 0;
            freq_num = 0;
            aps = {};

            iso_ori_mat[0][0]=1; iso_ori_mat[0][1]=0; iso_ori_mat[0][2]=0;
            iso_ori_mat[1][0]=0; iso_ori_mat[1][1]=1; iso_ori_mat[1][2]=0;
            iso_ori_mat[2][0]=0; iso_ori_mat[2][1]=0; iso_ori_mat[2][2]=1;

            iso_tors_turned = 0.0;
            iso_targeted_AP = -1;
            iso_frag_att_num_APs = -1;
            
            iso_centers_com.x = 0;
            iso_centers_com.y = 0;
            iso_centers_com.z = 0;

            iso_spheres_com.x = 0;
            iso_spheres_com.y = 0;
            iso_spheres_com.z = 0;
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
            radial_set = false;
            clear_radial_dist_distri();
            alloc_set = false;
            freq_num = 0;
            aps = {};

            iso_ori_mat[0][0]=1; iso_ori_mat[0][1]=0; iso_ori_mat[0][2]=0;
            iso_ori_mat[1][0]=0; iso_ori_mat[1][1]=1; iso_ori_mat[1][2]=0;
            iso_ori_mat[2][0]=0; iso_ori_mat[2][1]=0; iso_ori_mat[2][2]=1;

            iso_tors_turned = 0.0;
            iso_targeted_AP = -1;
            iso_frag_att_num_APs = -1;

            iso_centers_com.x = 0;
            iso_centers_com.y = 0;
            iso_centers_com.z = 0;

            iso_spheres_com.x = 0;
            iso_spheres_com.y = 0;
            iso_spheres_com.z = 0;
       };


//// Fragment class methods for isoalign
//// methods that sets radial distanace for iso_align protocol START
void Fragment::operator=(const Fragment & original){

    num_du                 =   original.num_du;
    iso_aligned            =   original.iso_aligned;
    iso_score              =   original.iso_score;
    radial_set             =   original.radial_set;
    
    copy_molecule((*this).mol, original.mol);
    for (unsigned int i =0; i<original.best_tri.three_pairs.size(); i++){
        (*this).best_tri.three_pairs.push_back(original.best_tri.three_pairs[i]);
    }
    (*this).alloc_set = original.alloc_set;
    if ((*this).alloc_set){
        (*this).allocate_radial_dist_distri();

        for (unsigned int i =0;i<(*this).mol.num_atoms;i++){
            for (unsigned int j=0; j<original.radial_dist_distri[i].size(); j++){
                (*this).radial_dist_distri[i].push_back(original.radial_dist_distri[i][j]);
            }
        }

    }


    for (unsigned int i =0; i<original.aps.size(); i++){
        aps.push_back(original.aps[i]);
    }
   
    used                   =   original.used;
    tanimoto               =   original.tanimoto;
    last_ap_heavy          =   original.last_ap_heavy;
    scaffolds_this_layer   = original.scaffolds_this_layer;

    for (unsigned int i =0; i<original.torenv_recheck_indices.size(); i++){
        torenv_recheck_indices.push_back(original.torenv_recheck_indices[i]);
    }
    
    for (unsigned int i =0; i< original.frag_growth_tree.size(); i++){
        DOCKMol   tmp_mol;
        copy_molecule(tmp_mol,original.frag_growth_tree[i]);
        (*this).frag_growth_tree.push_back(tmp_mol);
    }
    aps_bonds_type         = original.aps_bonds_type;
    size                   = original.size;
    ring_size              = original.ring_size;
   
    for (unsigned int i = 0; i< original.aps_cos.size(); i++){
        aps_cos.push_back(original.aps_cos[i]);
    }

    for (unsigned int i = 0; i<original.aps_dist.size(); i++){
        aps_dist.push_back(original.aps_dist[i]);
    }
    
    mut_type              = original.mut_type;
    freq_num              = original.freq_num;

    for(int i=0;i<3;i++){
        for( int j=0;j<3;j++){
            iso_ori_mat[i][j] = original.iso_ori_mat[i][j];
        }
    }
 
    iso_tors_turned = original.iso_tors_turned;
    iso_targeted_AP = original.iso_targeted_AP;
    iso_frag_att_num_APs = original.iso_frag_att_num_APs;
    iso_head_attached_ind = original.iso_head_attached_ind;

    iso_centers_com = original.iso_centers_com;
    iso_spheres_com = original.iso_spheres_com;
}


Fragment::Fragment(const Fragment & original){

    num_du                 =   original.num_du;
    iso_aligned            =   original.iso_aligned;
    iso_score              =   original.iso_score;
    radial_set             =   original.radial_set;
    
    copy_molecule((*this).mol, original.mol);
    for (unsigned int i =0; i<original.best_tri.three_pairs.size(); i++){
        (*this).best_tri.three_pairs.push_back(original.best_tri.three_pairs[i]);
    }
     (*this).alloc_set = original.alloc_set;
    if ((*this).alloc_set){
        (*this).allocate_radial_dist_distri();

        for (unsigned int i =0;i<(*this).mol.num_atoms;i++){
            for (unsigned int j=0; j<original.radial_dist_distri[i].size(); j++){
                (*this).radial_dist_distri[i].push_back(original.radial_dist_distri[i][j]);
            }
        }
    }

    for (unsigned int i =0; i<original.aps.size(); i++){
        aps.push_back(original.aps[i]);
    }
   
    used                   =   original.used;
    tanimoto               =   original.tanimoto;
    last_ap_heavy          =   original.last_ap_heavy;
    scaffolds_this_layer   =   original.scaffolds_this_layer;

    for (unsigned int i =0; i<original.torenv_recheck_indices.size(); i++){
        torenv_recheck_indices.push_back(original.torenv_recheck_indices[i]);
    }
    
    for (unsigned int i =0; i< original.frag_growth_tree.size(); i++){
        DOCKMol   tmp_mol;
        copy_molecule(tmp_mol,original.frag_growth_tree[i]);
        (*this).frag_growth_tree.push_back(tmp_mol);
    }
    aps_bonds_type         = original.aps_bonds_type;
    size                   = original.size;
    ring_size              = original.ring_size;
   
    for (unsigned int i = 0; i< original.aps_cos.size(); i++){
        aps_cos.push_back(original.aps_cos[i]);
    }

    for (unsigned int i = 0; i<original.aps_dist.size(); i++){
        aps_dist.push_back(original.aps_dist[i]);
    }
    
    mut_type              = original.mut_type;
    freq_num              = original.freq_num;

    for(int i=0;i<3;i++){
        for( int j=0;j<3;j++){
            iso_ori_mat[i][j] = original.iso_ori_mat[i][j];
        }
    }
 
    iso_tors_turned = original.iso_tors_turned;
    iso_targeted_AP = original.iso_targeted_AP;
    iso_frag_att_num_APs = original.iso_frag_att_num_APs;
    iso_head_attached_ind = original.iso_head_attached_ind;

    iso_centers_com = original.iso_centers_com;
    iso_spheres_com = original.iso_spheres_com;
}

void Fragment::set_radial_dist_distri(int atom_num, std::vector<float> radial){
    for (unsigned int i=0; i<radial.size(); i++){
        (*this).radial_dist_distri[atom_num].push_back(radial[i]);
    }
}

void Fragment::calc_radial_dist_distri(){

    Iso_Align iso_ali;
    iso_ali.get_all_radial_dist(*this);
    (*this).radial_set = true;
}

void Fragment::print_radial_dist_distri(){
    for (unsigned int i; i< (*this).mol.num_atoms; i++){
        for (unsigned int j=0; j<(*this).radial_dist_distri[i].size(); j++){
            std::cout << (*this).radial_dist_distri[i][j] << " " ;
            if ((j+1)==(*this).radial_dist_distri[i].size()){
                std::cout << std::endl;
            }
        }
    }

}

std::vector<float> Fragment::get_radial_dist_distri(int atom_num){
    return (*this).radial_dist_distri[atom_num];
}

void Fragment::allocate_radial_dist_distri(){
    radial_dist_distri    = new      std::vector<float> [(*this).mol.num_atoms];
    alloc_set = true;
}

bool Fragment::is_it_radial_set(){return (*this).radial_set;}

void Fragment::clear_radial_dist_distri(){

    if (alloc_set){
        delete[]radial_dist_distri; 
        alloc_set = false;
    }
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
void Fragment::set_iso_score(float score){
    this->iso_score = score;
}

float Fragment::get_iso_score(){
	
    if (this->iso_score== 9999){
	std::cout << "you are trying to access an unscored iso frag" << std::endl;
	exit(0);
    }
    return this->iso_score;
}
// Fragment class methods for isoalign END

void
Fragment::calc_mol_wt( )
{
    // atomic weights from General Chemistry 3rd Edition Darrell D. Ebbing
    float mw = 0.0;

    // For every atom in the dockmol object
    for (int i=0; i<mol.num_atoms; i++) {

        string atom = mol.atom_types[i];

        if ( atom == "H")
            { mw += ATOMIC_WEIGHT_H; }

        else if ( atom == "C.3" || atom == "C.2" || atom == "C.1" || atom == "C.ar" ||
                  atom  == "C.cat" )
            { mw += ATOMIC_WEIGHT_C; }

        else if ( atom == "N.4" || atom == "N.3" || atom == "N.2" || atom == "N.1" ||
                  atom  == "N.ar" || atom == "N.am" || atom == "N.pl3" )
            { mw += ATOMIC_WEIGHT_N; }

        else if ( atom == "O.3" || atom == "O.2" || atom == "O.co2" )
            { mw += ATOMIC_WEIGHT_O; }

        else if ( atom == "S.3" || atom == "S.2" || atom == "S.O" || atom == "S.o" ||
                  atom == "S.O2" || atom == "S.o2" )
            { mw += ATOMIC_WEIGHT_S; }

        else if ( atom == "P.3" )
            { mw += ATOMIC_WEIGHT_P; }

        else if ( atom == "F" )
            { mw += ATOMIC_WEIGHT_F; }

        else if ( atom == "Cl" )
            { mw += ATOMIC_WEIGHT_Cl; }

        else if ( atom == "Br" )
            { mw += ATOMIC_WEIGHT_Br; }

        else if ( atom == "I" )
            { mw += ATOMIC_WEIGHT_I; }

        else if ( atom == "Du" )
            { mw += 0.0; }

        else
            { cout <<"Warning: Did not recognize the atom_type " <<atom 
                   <<" in Fragment::calc_mol_wt()" <<endl; }

    }    

    (*this).mol.mol_wt =  mw;
}


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

//printing function for debug
void
Fragment::print(int index, std::string label){
    //cout << "num_du: " <<  num_du << endl << endl;
    //cout << "iso_aligned: " <<  iso_aligned << endl;
    //cout << "iso_score: " <<  iso_score << endl;
    //cout << "radial_set:" << radial_set << endl;
    //print_radial_dist_distri();
    string PRINT_DELIM = "##############";
    cout << PRINT_DELIM << endl;
    cout << "printing Fragment " << mol.energy <<  " index: " << index << " with label: " << label << endl;
    cout << "frag mw: " << mol.mol_wt << endl;
    cout << "mol # atoms: " << mol.num_atoms << endl;
    cout << "ring size: " << ring_size << endl;
    

    int print_vecs = 1;
    std::stringstream ss;
    for (int i = 0; i < aps.size(); i++){
        ss << "{ aps " << i << "th Du value: " << aps[i].dummy_atom << " Hvy Atom #:  " << aps[i].heavy_atom << "}; ";
    }
    if (print_vecs){
        //cout << ss.str() << endl;

        ss.str("");
    }

    for (int i = 0; i < aps_bonds_type.size(); i++){
        //ss << "{ aps_bonds_type " << i << "th value: " << aps_bonds_type[i].first;
        ss << "," << aps_bonds_type[i].second << "}; ";

    }
    if (print_vecs){
        ss << "aps size: " << aps.size() << " ## aps bond size: " << aps_bonds_type.size();
        cout << ss.str() << endl;
    }
    cout << PRINT_DELIM << endl;

}

