#include <string>
#include <vector>
#include <map>

#include "dockmol.h"
#include "utils.h"  // INTVec

#define ATOMIC_WEIGHT_H 1.00794
#define ATOMIC_WEIGHT_C 12.011
#define ATOMIC_WEIGHT_N 14.00647
#define ATOMIC_WEIGHT_O 15.9994
#define ATOMIC_WEIGHT_S 32.066
#define ATOMIC_WEIGHT_P 30.973762
#define ATOMIC_WEIGHT_F 18.9984032
#define ATOMIC_WEIGHT_Cl 35.4527
#define ATOMIC_WEIGHT_Br 79.904
#define ATOMIC_WEIGHT_I 126.90447

class Bump_Filter;
class Master_Score;
class Parameter_Reader;
class Simplex_Minimizer;


/********************************************************************/
class           SEGMENT {

  public:

    SEGMENT() {
        num_hvy_atoms = 0;
    };
    // //////////////////////////////////////////////////
    INTVec          atoms;      // list of atoms in the segment
    INTVec          bonds;      // list of bonds completely contained in the
                                // segment
    int             num_hvy_atoms;      // number of atoms in the segment
    INTVec          neighbors;  // neighboring segments
    INTVec          neighbor_bonds;     // rot bonds btwn segments
    INTVec          neighbor_atoms;     // test for segment at a time conf gen

    friend int      operator<(SEGMENT s1, SEGMENT s2) {
        return (s1.atoms.size() > s2.atoms.size());
    };
};


/********************************************************************/
class           LAYER_SEGMENT {

  public:

    LAYER_SEGMENT() {
        num_hvy_atoms = 0;
        rot_bond = -1;
        origin_segment = -1;
    };
    // //////////////////////////////////////////////////
    INTVec          atoms;      // list of atoms in the segment
    INTVec          bonds;      // list of bonds completely contained in the
                                // segment
    int             num_hvy_atoms;      // number of atoms in the segment
    int             rot_bond;   // rot bond to this segment from previous layer
    int             origin_segment;     // which segment this one is derived
                                        // from

    friend int      operator<(LAYER_SEGMENT s1, LAYER_SEGMENT s2) {
        return (s1.atoms.size() > s2.atoms.size());
    };
};


/********************************************************************/
class           LAYER {

  public:

    // //////////////////////////////////////////////////
    INTVec segments;            // List of segments in this layer
    int             num_segments;       // Number of segments in this layer

};


/********************************************************************/
class           ROT_BOND {

  public:

    // //////////////////////////////////////////////////
    int             atom1;      // first atom bound to atom 2 (to define angle)
    int             atom2;      // first atom of rotatable bond
    int             atom3;      // second atom of rotatable bond
    int             atom4;      // first atom bound to atom 3 (to define angle)
    int             seg1;       // Segment containing atom #2
    int             seg2;       // Segment containing atom #3
    int             bond_num;   // bond number in Mol
    float           initial_angle;      // initial angle from database
    float           current_angle;      // current angle

};


/********************************************************************/
class           CONFORMER {

  public:

    // General info
    float           score;      // score of current partial conf
    int             layer_num;  // current layer
    bool            used;       // flag to delete conformer once its used
    int             anchor_num;  // each conformer originated from anchor position.
    int             conformer_num; // independent number;
    int             parent_num;   // number;
    std::string     header;       // text with energy score, rmsd etc. for branch_*.mol2
 
    DOCKMol         structure;  // structure of conformer

};


/********************************************************************/
class           BRANCH {

  public:

    int             beatm,
                    bmatm,
                    bnhvy,
                    bnhyd,
                    confnum,
                    bi,
                    iconf;
    float           solvat,
                    apol;
    INTVec          level_values,
                    level_counts;

};


/********************************************************************/
class           ATOM_INFO {

  public:

    int             level,
                    vdwtype,
                    flagat,
                    lcolor;
    float           charge,
                    polsolv,
                    apolsolv;
    std::string     type;

    int             branch_num;
    INTVec          level_confs;

};


/********************************************************************/
class           AG_Conformer_Search {

  public:

    // growth tree confromers 
    int             count_conf_num;                //counter - trent balius Dec 13, 2008

    //bool flexible_ligand;                        // no
    bool            verbose;
    int             anchor_size;                   // 10

    // TEB ADD 2010-01-23
    bool            user_specified_anchor;
    std::string     atom_in_anchor;
    bool            limit_max_anchors;
    int             max_anchor_num;

    int             num_anchor_poses;
    float           anchor_score_cutoff;           // energy cutoff value for next layer

    int             num_growth_poses;
    //bool            growth_cutoff;
    float           growth_score_scaling_factor;
    float           growth_score_cutoff;           // energy cutoff value for growth
    float           growth_score_cutoff_begin;

    bool            cluster;
    int             pruning_clustering_cutoff;     // number of confs/layer

    bool            use_clash_penalty;             // flag to use clash penalty or not
    float           clash_penalty;                 // atom clash threshold value

    // move to Master_Conformer_Search 
    bool            use_internal_energy;           // int nrg function superseded by func in base_score
    int             ie_att_exp;
    int             ie_rep_exp;
    float           ie_diel;
    float          *ie_vdwA;
    float          *ie_vdwB; 
    float           internal_energy_cutoff;        //BCF internal energy cutoff

    bool            print_growth_tree;             // sudipto & trent 30-01-09 

    DOCKMol         orig;                          // original molecule struct

    std::vector < SEGMENT > orig_segments;         // 
    std::vector < LAYER_SEGMENT > layer_segments;  // 
    std::vector < LAYER > layers;                  // 
    std::vector < INTPair > anchors;               // pair list of anchor segments and sizes
    std::vector < SCOREMol > anchor_positions;     // 
    std::vector < DOCKMol > anchor_confs;          // 
    std::vector < ROT_BOND > bond_list;            // 
    std::vector < SCOREMol > pruned_confs;         //

    // list of anchors for current molecule after minimization and pruning
    std::vector < CONFORMER > conf_anchors;        // sudipto & trent - Dec 08, 2008
    std::vector < bool > assigned_atoms;           // flags to track when
                                                   // an atom is assigned
                                                   // to a layer
    bool            last_conformer;                // 

    int             current_anchor;                // current anchor (anchors.size() -> 0)
    bool            return_anchor_value;           // used when no flex ligand is
                                                   // requested
    bool            return_periph_value;           // used when no flex ligand is
                                                   // requested

    INTVec          atom_seg_ids;                  // id of which segment each atom
                                                   // belongs to
    INTVec          bond_seg_ids;                  // id of which segment each bond
                                                   // belongs to (flex bonds assigned -1)
    INTVec          bond_tors_vectors;             // id of "start" atom for each bond WRT 
                                                   // layers- for flexible minimization
                                                   // during growth
    //PAK vects
    std::vector < int > next_nbrs;
    std::vector < int > tmp_nextnbrs;
    std::vector < std::string > tmp_vec_atom_strings;
    std::vector < float > tmp_vec_atom_wt;
    // Functions in AG_Conformer_Search //////////////////////////////////////

    AG_Conformer_Search();
    virtual ~ AG_Conformer_Search();
    void            initialize();       // 
    void            initialize_internal_energy_parms(bool uie, int rep_exp, int att_exp, float diel, float iec); 
    void            initialize_internal_energy_null(bool uie); 
                    // get values from Master_Conformer_Search called in input_parameters of Master
    void            input_parameters(Parameter_Reader & parm);  // 
    void            prepare_molecule(DOCKMol &);        // 
    void            identify_rigid_segments(DOCKMol &); // 
    void            extend_segments(int, int, DOCKMol &);       // 
    void            id_anchor_segments();       // 
    bool            next_anchor(DOCKMol &);     // 
    void            extend_layers(int, int, int);       // 
    bool            submit_anchor_orientation(DOCKMol &, bool); // 
    float           calc_layer_rmsd(CONFORMER &, CONFORMER &);  // 
    float           calc_active_rmsd(CONFORMER &, CONFORMER &);  // 2008-11-17 trent balius add  
    void            grow_periphery(Master_Score &, Simplex_Minimizer &, Bump_Filter &);
    void            conf_header(CONFORMER &, std::string, Master_Score &);
    void            segment_torsion_drive(CONFORMER &, int, std::vector < CONFORMER > &, int);
    void            activate_layer_segment(DOCKMol &, int, int);        // 
    void            print_branch(std::vector < CONFORMER > &, std::vector < CONFORMER > &,const CONFORMER &,Master_Score & ); // trent balius 2008-12-03
    void            print_conformer(CONFORMER &); // trent balius 2008-12-08
    void            reset_active_lists(DOCKMol &);      // 
    bool            segment_clash_check(DOCKMol &, int, int);   // 
    bool            next_conformer(DOCKMol &);  // 
    static int      conformer_less_than(CONFORMER a, CONFORMER b) {
        return a.score < b.score;
    };
    bool            atom_in_anchor_segments(SEGMENT); // TEB ADD 2010-01-23

    void            print_atom_in_anchor_segments(SEGMENT); // TEB ADD 2019-06-17

    // For generating fragment libraries and torsion environments for de novo growth
    bool            write_fragment_libraries;		// if true, libraries are written, most of dock.cpp is skipped
    std::string     fragment_library_prefix;		// prefix for libraries (four files total)
    int             fragment_library_freq_cutoff;	// frequency cutoff for writing fragments to library
    std::string     fragment_library_sort_method;       // method for sorting libraries, freq or fingerprint
    bool            fragment_library_trans_origin;	// if true, trans frags to origin before writing libs

    std::map < std::string, std::pair <DOCKMol, int> > segment_fingerprints; // hash for writing fragments
    std::map < std::string, int >                      torsions_map;         // hash for writing tor envs
    std::map < std::string, std::string >              torsions_map_ref;     // hash for writing tor envs

    int             global_frag_index;

    void            count_fragments(DOCKMol &);		// writes fragments to libraries depending on # of attachment points
    void            activate_fragment(DOCKMol &, int);	// activates atoms/bonds of a segment + neighbors
    void            write_unique_fragments();           // rewrites unique fragment libraries
    void            calc_mol_wt(DOCKMol &);
    std::vector <float>                               calc_atoms_wt(std::vector<std::string>, bool); 


};


/********************************************************************/
// Sort Functions
int       frequency_sort(std::pair <std::string, int>, std::pair <std::string, int> );
int       fingerprint_sort(std::pair <std::string, int>, std::pair <std::string, int> );
