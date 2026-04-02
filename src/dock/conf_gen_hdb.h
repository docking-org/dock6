#include <string>
#include <vector>

#include "dockmol.h"
#include "orient.h"
#include "utils.h"  // INTVec

class Bump_Filter;
class Master_Score;
class Parameter_Reader;
class Simplex_Minimizer;


/********************************************************************
class           HDB_MULTICONF {

  public:
    INTVecVecVec branch_confs;
    FLOATVecVec     branch_conf_scores;
    DOCKMol         mol;

    void            clear() {
        branch_confs.clear();
        branch_conf_scores.clear();
        mol.clear_molecule();
    };

};
*/

/********************************************************************
 * db2 file has atom section that defines the atom type its charge and desolvation value
class           HDB_atom {
     int            atom_num;
     int            atom_type;
     std::string    atom_name;
     float          q  ;
     float          p_desolv  ;
     float          a_desolv  ;
     float          sasa  ;

     void           initialize();
     HDB_atom();
     ~HDB_atom();
};
*/

/********************************************************************
 * db2 file has bond section this defineds conectivity 
class           HDB_bond {
     int            bond_num;
     int            atom1_num;
     int            atom2_num;
     std::string    bond_name; // sybyl bond type

     void           initialize();
     HDB_bond();
     ~HDB_bond();
};
*/


/********************************************************************
class           HDB_coordinate {
     int            coord_num;
     float          x  ;
     float          y  ;
     float          z  ;
     int            atom_num;

     void           initialize();
     HDB_coordinate();
     ~HDB_coordinate();
};
*/

/********************************************************************
class           HDB_segment {
      int seg_num;    
      int start_coord;    
      int stop_coord;    
  
     void           initialize();
     HDB_segment();
     ~HDB_segment();
};
*/

/********************************************************************
class           HDB_conformer {

    int conf_num;
    std::vector <int> list_of_seg;

     void           initialize();
     HDB_conformer();
     ~HDB_conformer();
};
*/


/********************************************************************/
class           HDB_Conformer_Search {

  public:

    bool flexible_ligand;
    std::string     db2filename;

/*
    // heirarchy data
    // db2 info when reading in
    HDB_atom *      atoms;
    HDB_bond *      bonds;
    HDB_segment *   segs;
    HDB_conformer * conf; 
*/

    HDB_Mol         total_hdb_mol;
    // store db2 info in total_mol for scoring.
    DOCKMol         total_mol;
    float           score_thres;
    //int             num_per_hierarchy;
    int             num_per_search;

    // Here is what TEB thinks we need, 2019. 
    // atoms ^A -- type number charge
    // bonds ^B -- 2 atom numbers
    // rigid segement ^R 
    // coordate list ^X -- this is all of the coordenates for all of the atoms in the database of conformations overlaping regons only appear once. 
    // segment list ^C -- conf in the db2 formate -- this is on segment of the molecule in a specific conformation
    // conformer list ^S -- set in the db2 formate -- this discribes a branch or a fully grown molecule

    // family level data

    // atom hierarchy assignments

    int             num_anchor_poses;
    std::vector < SCOREMol > anchor_positions;
//    std::vector < HDB_MULTICONF > docked_confs;
    std::vector < DOCKMol > final_confs;
    std::vector < SCOREMol > all_poses;

    INTVec          conf_state;
    INTVec          conf_state_max;

    bool            return_anchor_value;        // used when no flex ligand is
                                                // requested
    bool            return_periph_value;        // used when no flex ligand is
                                                // requested

    HDB_Conformer_Search();
    virtual ~ HDB_Conformer_Search();
    void            initialize();
    void            input_parameters(Parameter_Reader &);
    void            prepare_molecule(HDB_Mol &);
    void            create_mol(DOCKMol &,DOCKMol &,HDB_Mol &, int); // this will put the conf in a dockmol and activate the atoms
    void            set_conf_mol(DOCKMol &,HDB_Mol &); // this will put the conf in a dockmol and activate the atoms
    void            set_branch_mol(DOCKMol &, DOCKMol &, HDB_Mol &, int); // this will put the confs from a branch in dockmol and activate all the atoms
    bool            next_anchor(DOCKMol &);
    bool            submit_anchor_orientation(DOCKMol &, bool ); // 
//    void            grow_periphery(Master_Score &, Bump_Filter & );
    void            search(Master_Score &, Orient &, Bump_Filter &, DOCKMol &,DOCKMol &,HDB_Mol &,int);
    void            deactivate_molecule(DOCKMol &);
    void            activate_anchor(DOCKMol &);
    bool            next_conformer(DOCKMol &);
    //bool            generate_mols_from_confs(HDB_MULTICONF &, DOCKMol &);

};


typedef         std::pair < int, float > INT_FLOAT_Pair;

static int      int_float_pair_less_than(INT_FLOAT_Pair a, INT_FLOAT_Pair b) {
        return a.second < b.second;
        //return a.second > b.second;
};
