#include <string>
#include <vector>
#include <map>

#include "dockmol.h"
#include "utils.h"  // INTVec
#include "conf_gen_ag.h"

class Bump_Filter;
class Master_Score;
class Parameter_Reader;
class Simplex_Minimizer;
/*
class           AG_Conformer_Search ;
class           SEGMENT ;
class           LAYER_SEGMENT ;
class           LAYER ;
class           ROT_BOND ;
class           CONFORMER ;
class           BRANCH ;
class           ATOM_INFO ;
*/

/********************************************************************/
class           CG_Conformer_Search {

  public:

    // growth tree confromers 
    int             count_conf_num;                //counter - trent balius Dec 13, 2008

    //bool flexible_ligand;                        // no
    bool            verbose;
    int             anchor_size;                   // 10

    // covalent bond
    int             dummy1 ;                    // attachment point 
    int             dummy2 ;                    // attachment point , vector
    int             atomtag ;
    //float           angle;                         // angle of attachment   
    //float           bond1len;
    //float           bond2;


    // TEB ADD 2010-01-23
    bool            user_specified_anchor;
    std::string     atom_in_anchor;
    bool            limit_max_anchors;
    int             max_anchor_num;

    int             num_anchor_poses;
    float           anchor_score_cutoff;           // energy cutoff value for next layer

    int             num_growth_poses;
    bool            growth_cutoff;
    float           growth_score_cutoff;           // energy cutoff value for growth

    bool            cluster;
    int             pruning_clustering_cutoff;     // number of confs/layer

    float           bondlenth;                     // bond lenth of covalent angle.
    //int             num_sample_angles;           // number of angles to sample about covalent bond.
    float           dihideral_step;                // dihideral step size about the covalent attachment.

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

    // Functions in CG_Conformer_Search //////////////////////////////////////

    CG_Conformer_Search();
    virtual ~ CG_Conformer_Search();
    void            initialize();       // 
    void            initialize_internal_energy_parms(bool uie, int rep_exp, int att_exp, float diel, float iec); 
                    // get values from Master_Conformer_Search called in input_parameters of Master
    void            input_parameters(Parameter_Reader & parm);  // 
    void            prepare_molecule(DOCKMol &);        // 
    void            identify_rigid_segments(DOCKMol &); // 
    void            extend_segments(int, int, DOCKMol &);       // 
    void            id_anchor_segments(DOCKMol &);       // 
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

    void            set_torsion(float x0, float x1, float x2, float x3,
                                float y0 , float y1, float y2, float y3,
                                float z0 , float z1, float z2, float z3,
                                float angle, DOCKMol & mol);

    static int      conformer_less_than(CONFORMER a, CONFORMER b) {
        return a.score < b.score;
    };
    bool            atom_in_anchor_segments(SEGMENT, DOCKMol &); // TEB ADD 2010-01-23


    // For generating fragment libraries and torsion environments for de novo growth
    bool            write_fragment_libraries;		// if true, libraries are written, most of dock.cpp is skipped
    std::string     fragment_library_prefix;		// prefix for libraries (four files total)
    int             fragment_library_freq_cutoff;	// frequency cutoff for writing fragments to library
    std::string     fragment_library_sort_method;       // method for sorting libraries, freq or fingerprint
    bool            fragment_library_trans_origin;	// if true, trans frags to origin before writing libs

    std::map < std::string, std::pair <DOCKMol, int> > segment_fingerprints; // hash for writing fragments
    std::map < std::string, int >                      torsions_map;         // hash for writing tor envs
    int             global_frag_index;

//    void            count_fragments(DOCKMol &);		// writes fragments to libraries depending on # of attachment points
//    void            activate_fragment(DOCKMol &, int);	// activates atoms/bonds of a segment + neighbors
//    void            write_unique_fragments();           // rewrites unique fragment libraries


};


