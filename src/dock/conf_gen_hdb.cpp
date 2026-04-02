/** HDB_IS_NOT_YET_COMPLETELY_IMPLEMENTED **/

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
//#include <string>
#include <limits.h>
#include <time.h>

#include "conf_gen_hdb.h"
#include "orient.h"
#include "master_score.h"
#include "simplex.h"
#include "trace.h"
#include "math.h"

class Bump_Filter;
class Parameter_Reader;

using namespace std;


// +++++++++++++++++++++++++++++++++++++++++
 HDB_Conformer_Search::HDB_Conformer_Search()
 {
 }
// +++++++++++++++++++++++++++++++++++++++++

HDB_Conformer_Search::~HDB_Conformer_Search()
 {
  all_poses.clear(); 
 }
// +++++++++++++++++++++++++++++++++++++++++

void 
HDB_Conformer_Search::input_parameters(Parameter_Reader &parm) 
{
        string        tmp;

        cout << "\nHierarchy DB Search Parameters" << endl;
        cout << "------------------------------------------------------------------------------------------" << endl;

        //tmp = parm.query_param("hdb_conf_search","no","yes no");
        //num_per_hierarchy = atoi(parm.query_param("num_per_hierarchy", "10").c_str());
        num_per_search = atoi(parm.query_param("num_per_search", "10").c_str());

        //if(tmp == "yes")
        //        flexible_ligand = true;
        //else
                flexible_ligand = false;
        // this will eventually read in an index file containing a list of gz index files
        //db2filename = parm.query_param("hdb_db2_input_file","hdb_db2_file.db2"); 
        db2filename = parm.query_param("hdb_db2_input_file","hdb_db2_file.db2.gz"); 
        score_thres = atof(parm.query_param("hdb_db2_search_score_threshold", "10.0").c_str());

}


// +++++++++++++++++++++++++++++++++++++++++
void 
HDB_Conformer_Search::initialize() 
{
  // store the fist branch in the molecule.  
  //
  //total_mol.allocate_arrays(num_of_atoms, num_of_bonds, 0)

}


// +++++++++++++++++++++++++++++++++++++++++
void 
HDB_Conformer_Search::prepare_molecule(HDB_Mol &mol) 
{



}

// +++++++++++++++++++++++++++++++++++++++++
void 
HDB_Conformer_Search::create_mol(DOCKMol &mol, DOCKMol & mol_ac, HDB_Mol &hmol, int confnum) 
{
    int i,j,k,l; 
    cout << "HDB_Conformer_Search::create_mol" << endl;
    mol.allocate_arrays(hmol.num_of_atoms, hmol.num_of_bonds,0);

    // copy atoms
    //
    mol.title     = hmol.name;
    for (i=0; i < hmol.num_of_atoms; i++){
        // mol.= hmol
        mol.charges[i]     = hmol.atoms[i].q;
        mol.atom_types[i]  = hmol.atoms[i].atom_syb_type;
        mol.atom_names[i]  = hmol.atoms[i].atom_name;
        mol.atom_number[i]          = hmol.atoms[i].atom_num;
        mol.atom_residue_numbers[i] = "";
        mol.atom_psol[i]            = hmol.atoms[i].p_desolv;    // kxr205
        mol.atom_apsol[i]           = hmol.atoms[i].a_desolv;  // kxr205
        //mol.subst_names[i]          = original.subst_names[i];
        //mol.atom_color[i]           = hmol.atoms[i].atom_color;  // kxr205
    }
    // copy bonds
    //
    for (j=0; j < hmol.num_of_bonds; j++){
       //cout << "temp" <<endl;
    //int            *bonds_origin_atom;
    //int            *bonds_target_atom;
    //std::string    *bond_types;
    //int            bond_num;
    //int            atom1_num;
    //int            atom2_num;
    //std::string    bond_name; // sybyl bond type
      //cout << j << " " << hmol.bonds[j].atom1_num << " " << hmol.bonds[j].atom2_num << " " << hmol.bonds[j].bond_name << endl;
      mol.bonds_origin_atom[j] = hmol.bonds[j].atom1_num-1;
      mol.bonds_target_atom[j] = hmol.bonds[j].atom2_num-1;
      mol.bond_types[j]        = hmol.bonds[j].bond_name;
    }

    // look up coordenates for a branch.  confs -> seg -> coord -> atom . 
    //cout << "number of segments:" << "  " << hmol.confs[confnum].num_of_seg << endl;
    //cout << "conf list seg: = " <<endl;
    // write out segment numbers: 
    //for (k=0;k<hmol.confs[confnum].num_of_seg;k++){ 
    //   cout << "  " << hmol.confs[confnum].list_of_seg[k] ; 
    //}
    //cout << "\nstop list seg" << endl;
    for (k=0;k<hmol.confs[confnum].num_of_seg;k++){ 
       //cout << "seg_num" << "  " << hmol.segs[hmol.confs[confnum].list_of_seg[k]-1].seg_num << endl;
       //cout << "coord_start"<< "  " << hmol.segs[hmol.confs[confnum].list_of_seg[k]-1].start_coord <<endl;
       //cout << "coord_stop" << "  " << hmol.segs[hmol.confs[confnum].list_of_seg[k]-1].stop_coord << endl;

       int start = hmol.segs[hmol.confs[confnum].list_of_seg[k]-1].start_coord;
       int stop  = hmol.segs[hmol.confs[confnum].list_of_seg[k]-1].stop_coord;
       for (l=start-1; l<stop;l++){
           //hmol.coords[l].print_coord();
           int an = hmol.coords[l].atom_num-1;
           //hmol.atoms[an].print_atom();
           mol.x[an] = hmol.coords[l].x; 
           mol.y[an] = hmol.coords[l].y;
           mol.z[an] = hmol.coords[l].z; 
       }
       
    }
    //ofstream myfile;
    //myfile.open ("debug.mol2");
    //Write_Mol2(mol, myfile);

    //  store all coords in a dockmol object.  
    //DOCKMol mol_ac; //# all_coords_rigid_seg;
    mol_ac.allocate_arrays(hmol.num_of_coords+hmol.num_of_rigid,0,0); 
    for (l=0; l<hmol.num_of_coords;l++){
//        mol_ac.atom_types[l] = "C"; 
//        mol_ac.atom_names[l] = "C"; 
        mol_ac.x[l] = hmol.coords[l].x; 
        mol_ac.y[l] = hmol.coords[l].y;
        mol_ac.z[l] = hmol.coords[l].z; 
        mol_ac.atom_active_flags[l] = false;

        int an = hmol.coords[l].atom_num-1;
        //cout << "atom_num == " << an << endl;
        mol_ac.charges[l]              = hmol.atoms[an].q;
        mol_ac.atom_types[l]           = hmol.atoms[an].atom_syb_type;
        mol_ac.atom_names[l]           = hmol.atoms[an].atom_name;
        mol_ac.atom_number[l]          = hmol.atoms[an].atom_num;
        mol_ac.atom_residue_numbers[l] = "";
        mol_ac.atom_psol[l]            = hmol.atoms[an].p_desolv;  // kxr205
        mol_ac.atom_apsol[l]           = hmol.atoms[an].a_desolv;  // kxr205
    }
    //exit (0);
    for (l=0; l<hmol.num_of_rigid;l++){
        mol_ac.atom_types[hmol.num_of_coords+l] = "C.3"; 
        mol_ac.atom_names[hmol.num_of_coords+l] = "C.3"; 
        mol_ac.x[hmol.num_of_coords+l] = hmol.rigid[l].x; 
        mol_ac.y[hmol.num_of_coords+l] = hmol.rigid[l].y;
        mol_ac.z[hmol.num_of_coords+l] = hmol.rigid[l].z; 
        mol_ac.atom_active_flags[hmol.num_of_coords+l] = true;
        mol_ac.amber_at_heavy_flag[hmol.num_of_coords+l] = true;
    }

    
    //myfile.close();

    //ofstream myfile;
    //myfile.open ("debug2.mol2");
    //Write_Mol2(mol_ac, myfile);
    //myfile.close();

//   for (i=0,i<hmol.num_of_atoms

}

// +++++++++++++++++++++++++++++++++++++++++
void 
HDB_Conformer_Search::set_branch_mol(DOCKMol &mol, DOCKMol &mol_ac, HDB_Mol &hmol, int confnum_in) 
{
    //int confnum = confnum_in-1;
    int confnum = confnum_in;
    //mol.current_data =  mol_ac.current_data;
    //mol.current_score = mol_ac.current_score;
  
    for (int k=0;k<hmol.confs[confnum].num_of_seg;k++){ 
       //cout << "seg_num" << "  " << hmol.segs[hmol.confs[confnum].list_of_seg[k]-1].seg_num << endl;
       //cout << "coord_start"<< "  " << hmol.segs[hmol.confs[confnum].list_of_seg[k]-1].start_coord <<endl;
       //cout << "coord_stop" << "  " << hmol.segs[hmol.confs[confnum].list_of_seg[k]-1].stop_coord << endl;

       int start = hmol.segs[hmol.confs[confnum].list_of_seg[k]-1].start_coord;
       int stop  = hmol.segs[hmol.confs[confnum].list_of_seg[k]-1].stop_coord;
       for (int l=start-1; l<stop;l++){
           //hmol.coords[l].print_coord();
           int an = hmol.coords[l].atom_num-1;
           //hmol.atoms[an].print_atom();
           //mol.x[an] = hmol.coords[l].x; 
           //mol.y[an] = hmol.coords[l].y;
           //mol.z[an] = hmol.coords[l].z; 
           mol.x[an] = mol_ac.x[l]; 
           mol.y[an] = mol_ac.y[l];
           mol.z[an] = mol_ac.z[l]; 
           mol.atom_active_flags[an] = true;
       }
       
    }


}

// +++++++++++++++++++++++++++++++++++++++++
// there is only one anchor specified per hiarchy 
bool 
HDB_Conformer_Search::next_anchor(DOCKMol &mol) 
{

        return true;

}


// +++++++++++++++++++++++++++++++++++++++++
bool 
HDB_Conformer_Search::submit_anchor_orientation(DOCKMol &mol, bool more_orients ) 
{
                  return false;

}



// +++++++++++++++++++++++++++++++++++++++++
void 
HDB_Conformer_Search::activate_anchor(DOCKMol &mol) 
{

}

/*
// +++++++++++++++++++++++++++++++++++++++++
bool 
HDB_Conformer_Search::generate_mols_from_confs(HDB_MULTICONF &conf, DOCKMol &ret_mol) 
{
        return false;
}
*/

float score_segment( Master_Score &score, Orient &orient, DOCKMol &mol, DOCKMol & mol_ac, HDB_Mol &hmol, int i ){
       float seg_score;
       //activate the segment atoms
       int atom_start = hmol.segs[i].start_coord;
       int atom_stop = hmol.segs[i].stop_coord;
       //cout << "  start = " << atom_start ;
       //cout << "  stop  = " << atom_stop << endl;
       //cout << "atom_nums = ";
       for(int j = atom_start-1; j<atom_stop; j++){
          int atom_num = hmol.coords[j].atom_num-1;
          //cout << " " << atom_num;
          //mol_ac.atom_active_flags[j] = true;
          mol.x[atom_num] = mol_ac.x[j];
          mol.y[atom_num] = mol_ac.y[j];
          mol.z[atom_num] = mol_ac.z[j];
          mol.atom_active_flags[atom_num] = true;
          //cout << mol.atom_types[atom_num] << " " << mol_ac.atom_types[j] << endl;
       }
       //cout << endl;
       //Write_Mol2(mol_ac, myfile);

       //if (score.compute_primary_score(mol)){
       //float s_flag = score.compute_primary_score(mol_ac);
       //cout << "I AM HERE" << endl;
       orient.orientation_HDB(mol);
       bool seg_flag = score.compute_primary_score(mol);
       //float s_flag = score.compute_primary_score(mol);
       //seg_flag[i] = score.compute_primary_score(mol_ac);
       if (seg_flag){
           float s = mol.current_score;
           seg_score = s;
           //cout << "seg_score, seg = " << i << " score = " << s << endl;
       }
       else{  
           seg_score = 100000.0;
       }
       //myfile << "###### score1:   " << seg_scores[i] << endl;
       //Write_Mol2(mol_ac, myfile);
       //Write_Mol2(mol, myfile);
       // deactivate (turn off) just atoms that are currently turned on
       for(int j = atom_start-1; j<atom_stop; j++){
          int atom_num = hmol.coords[j].atom_num-1;
          mol.atom_active_flags[atom_num] = false;
       }
       return seg_score;
}

// +++++++++++++++++++++++++++++++++++++++++
// Trent Balius, FNLCR, 2020/10/01
// The db2 search is fast, 
// we start with each conformer (branch), 
// we then will orient and score each segement in the branch, 
// if one segment clashes then go to the next branch, 
// if the anchor (rigid segment) clashes go to next orient.
void 
HDB_Conformer_Search::search(Master_Score &score, Orient &orient, Bump_Filter &bump, DOCKMol &mol, DOCKMol & mol_ac, HDB_Mol &hmol, int num ) 
{

   //ofstream myfile;
   //myfile.open ("debug3.mol2");
   //Write_Mol2(mol_ac, myfile);
   //myfile.close();

   char filename[50];
   //sprintf(filename,"debug_confs_%d.mol2",num);
   //myfile.open ("debug_confs.mol2");
   //myfile.open (filename);

   //cout << filename << endl;

   float *seg_scores;
   bool *seg_flag;

   seg_scores = new float[hmol.num_of_segs];
   seg_flag = new bool[hmol.num_of_segs];

   for (int i=0; i< hmol.num_of_segs; i++){
        seg_scores[i] = 10000.0;
   }

   float *conf_scores;
   bool *conf_flag;

   conf_scores = new float[hmol.num_of_confs];
   conf_flag = new bool[hmol.num_of_confs];

   for (int i=0; i< hmol.num_of_confs; i++){
        conf_scores[i] = 0.0;
   }
   // deactivate (turn off) all atoms
   for(int j = 0; j<mol_ac.num_atoms;j++){
      mol_ac.atom_active_flags[j] = false;
   } 

   for(int j = 0; j<mol.num_atoms;j++){
       mol.atom_active_flags[j] = false;
   } 

   // make amber parms agree.  
   /*
   for(int j = 0; j<mol_ac.num_atoms;j++){
          int atom_num = hmol.coords[j].atom_num-1;
          mol_ac.amber_at_id[j]         = mol.amber_at_id[atom_num]        ;
          mol_ac.amber_at_radius[j]     = mol.amber_at_radius[atom_num]    ;
          mol_ac.amber_at_well_depth[j] = mol.amber_at_well_depth[atom_num];
          mol_ac.amber_at_heavy_flag[j] = mol.amber_at_heavy_flag[atom_num];
          mol_ac.amber_at_valence[j]    = mol.amber_at_valence[atom_num]   ;
          mol_ac.amber_at_bump_id[j]    = mol.amber_at_bump_id[atom_num]   ;
          mol_ac.gb_hawkins_radius[j]   = mol.gb_hawkins_radius[atom_num]  ;
          mol_ac.gb_hawkins_scale[j]    = mol.gb_hawkins_scale[atom_num]   ;
   }
   */

   //cout << "I AM HERE (1)" << endl;
   //bool seg_flag[<hmol.num_of_segs];
   
   for (int i=0;i<hmol.num_of_segs;i++){
       seg_flag[i] = true; // has not been seen yet. 
   }

/*
   for (int i=0;i<hmol.num_of_segs;i++){
       seg_scores[i] = score_segment( score, mol, mol_ac, hmol, i );
       //cout << endl;
   } 
   //myfile.close();
*/

   //float score_thres = 10.0;
   // loop over all conformations.  
   float s = 0.0;
   for (int i=0; i< hmol.num_of_confs; i++){
      conf_flag[i] = true;
      //cout << i << " num of seg " << hmol.confs[i].num_of_seg<< endl;
      for (int k=0;k<hmol.confs[i].num_of_seg;k++){ // loop of each seg in conformation.  
           //hmol.segs[hmol.confs[i].list_of_seg[k]-1]
           //if (seg_flag[hmol.confs[i].list_of_seg[k]-1]){ 
           if (seg_flag[hmol.confs[i].list_of_seg[k]-1]) { // if segment score has not been calculated;
              //seg_scores[hmol.confs[i].list_of_seg[k]-1] < score_thres){ 
              s = score_segment( score, orient, mol, mol_ac, hmol, (hmol.confs[i].list_of_seg[k]-1) );
              //cout << "seg " << hmol.confs[i].list_of_seg[k] << ": "<< seg_scores[hmol.confs[i].list_of_seg[k]-1] << endl;
              seg_scores[hmol.confs[i].list_of_seg[k]-1] = s;
              seg_flag[hmol.confs[i].list_of_seg[k]-1] = false; // already seen once
           }else{ 
              s = seg_scores[hmol.confs[i].list_of_seg[k]-1];
           }
           conf_scores[i] = conf_scores[i] + seg_scores[hmol.confs[i].list_of_seg[k]-1];

           if (s > score_thres) { 
              if (i == 0) {// this is the anchor or rigid segment. if it is clashing then skip the whole orient.  
                  //cout << " anchor or rigid segment is clashing, so skip orient" << endl;
                  return;
              }
              conf_scores[i] = 10000.000;
              conf_flag[i] = false;
              //cout << "break out of loop" << endl;
              break;
           }
      }
      //cout << "conf " << i << ": " << "score " << conf_scores[i] << endl;
   }

   // sort molecules
   //
   vector <INT_FLOAT_Pair> conf_score_index ;
   conf_score_index.clear();

   for (int i=0; i< hmol.num_of_confs; i++){
      if (conf_flag[i]){
         INT_FLOAT_Pair val;
         val.first = i;
         val.second = conf_scores[i];
         conf_score_index.push_back(val);
      }
   } 

   //int num_to_add = 10;
   int num_to_add = num_per_search;
   if (conf_score_index.size() > 1) sort(conf_score_index.begin(), conf_score_index.end(), int_float_pair_less_than);
   if (conf_score_index.size() < num_to_add) num_to_add = conf_score_index.size();


   // get best scoring molecules
   //myfile.open ("debug_best_conf.mol2","app");
   //char filename[50];
   //sprintf(filename,"debug_best_conf_%d.mol2",num);
   //myfile.open (filename);
   SCOREMol smol;
   for (int i=0; i< num_to_add; i++){
       //cout << conf_score_index[i].first << " " << conf_score_index[i].second << endl;
       set_branch_mol(mol, mol_ac, hmol, conf_score_index[i].first) ;
       orient.orientation_HDB(mol);
       score.compute_primary_score(mol);
       //myfile << "###### score1:   " << conf_score_index[i].second << endl;
       //cout << "###### score1:   " << conf_score_index[i].second << endl;
       //myfile << "###### score2:   " << mol.current_score << endl;
       smol.first = conf_score_index[i].second;
       copy_molecule(smol.second, mol); 
       //copy_crds(smol.second, mol); 
       all_poses.push_back(smol);
       //Write_Mol2(mol, myfile);
   }
   //cout <<  "all_poses size = " << all_poses.size()<<endl;
   //myfile.close();

}


// +++++++++++++++++++++++++++++++++++++++++
bool 
HDB_Conformer_Search::next_conformer(DOCKMol &mol) 
{

                        return false;
}

