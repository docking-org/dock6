#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
//#include <cstdlib>

#include "amber_typer.h"
#include "dockmol.h"
#include "fragment.h"
#include "fraggraph.h"
#include "master_score.h"
#include "orient.h"
#include "simplex.h"
#include "utils.h"

// These covalent radii are defined according to the CRC Handbook
// of Chemistry and Physics, 94th Ed., 2013-2014, pp 9-49 to 9-50
#define COV_RADII_H 0.32
#define COV_RADII_C 0.75
#define COV_RADII_N 0.71
#define COV_RADII_O 0.64
#define COV_RADII_P 1.09
#define COV_RADII_S 1.04
#define COV_RADII_F 0.60
#define COV_RADII_CL 1.00
#define COV_RADII_BR 1.17

// Atomic weights from General Chemistry 3rd Edition Darrell D. Ebbing
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

/* LEP - test then delete this later
// +++++++++++++++++++++++++++++++++++++++++
// Attachment points are the dummy atom + heavy atom it is connected to
class           AttPoint {
  public:
    int    dummy_atom;
    int    heavy_atom;

    AttPoint();
    ~AttPoint();
};


// +++++++++++++++++++++++++++++++++++++++++
// Scaffolds, Linkers, Sidechains, and Rigids are condensed into one class called Fragment
class           Fragment {
  public:
    DOCKMol                                 mol;
    std::vector <AttPoint>                  aps;
    bool                                    used;
    float                                   tanimoto;
    int                                     last_ap_heavy;
    int                                     scaffolds_this_layer;
    std::vector < std::pair < int, int > >  torenv_recheck_indices;
    std::vector < DOCKMol >                 frag_growth_tree;

    Fragment();
    ~Fragment();
};
*/

// +++++++++++++++++++++++++++++++++++++++++
// Torsion Environments: A data structure that remembers an atom enviroment (as string) and a vector
// of all the things it can be connected to, as well as corresponding frequencies
class           TorEnv {
  public:
    std::string                  origin_env;
    std::vector <std::string>    target_envs;
    std::vector <int>            target_freqs;
    // use this for storing the environment frequency as well
    // std::vector < std::pair <std::string, int> >  target_envs;

    TorEnv();
    ~TorEnv();
};

/*
// +++++++++++++++++++++++++++++++++++++++++
// Fragment Graph: A data structure that stores a multi-directional graph of similarities for 
// fragment libraries based on Tanimoto (or something similar).
class           FragGraph {
  public:

    std::vector < std::pair <float, int> >   tanvec;
    std::vector <int>                        rankvec;
    bool                                     visited;

    FragGraph();
    ~FragGraph();
};
*/

// +++++++++++++++++++++++++++++++++++++++++
// De novo build molecules from fragment libraries
class           DN_Build {

  public:

    /** Variables **/

    // Format for output components the same as Base_Score
    static const std::string DELIMITER;
    static const int         FLOAT_WIDTH;
    static const int         STRING_WIDTH;

    // Input parameters: files / filenames
    std::string        dn_fraglib_scaffold_file;
    std::string        dn_fraglib_linker_file;
    std::string        dn_fraglib_sidechain_file;
    //std::string      dn_fraglib_rigid_file;
    bool               dn_user_specified_anchor;    // user can specify anchor(s) (good for lead opt)
    std::string        dn_fraglib_anchor_file;
    bool               dn_use_torenv_table;         // use torsion environment table?
    bool               dn_use_roulette;             //
    std::string        dn_torenv_table;
    std::string        denovo_name;                 // LEP- unique denovo name

    // Input parameters: sampling method
    bool               dn_refinement_random;        // LEP put in Du's for H's
    int                dn_refinement_random_picks;  // LEP
    std::string        dn_sampling_method;          // type of sampling (ex, rand, graph)
    bool               dn_sampling_method_ex;
    bool               dn_sampling_method_rand;
    bool               dn_sampling_method_graph;
    int                dn_num_random_picks;         // number of random picks for samp method rand
    int                dn_graph_max_picks;
    int                dn_graph_breadth;
    int                dn_graph_depth;
    float              dn_graph_temperature;
    float              dn_temp_begin;

    // Input parameters: pruning
    float              dn_pruning_conformer_score_cutoff;  // hard upper cutoff for the score during torsion sampling
    float              dn_pruning_conformer_score_cutoff_begin;  // hard upper cutoff for the score during torsion sampling
    float              dn_pruning_conformer_score_scaling_factor;
    float              dn_pruning_clustering_cutoff;
    float              dn_upper_constraint_mol_wt;       // upper bound for mol wt
    float              dn_lower_constraint_mol_wt;       // lower bound for mol wt
    std::string        dn_MW_cutoff_type;           // user specified cutoff type
    bool               dn_MW_cutoff_type_hard;      //  specified type == hard eval true and type == soft eval false
    float              dn_MW_std_dev;              //user defined standard deviation when using a soft MW cutoff

    // Directed de novo parameters
    #ifdef BUILD_DOCK_WITH_RDKIT

    // paths for RDKit-related descriptors
    std::string        sa_fraglib_path;
    std::string        PAINS_path;
    // SA fraglib 
    std::map<unsigned int, double> fragMap;
    // PAINSlib
    std::map<std::string, std::string> PAINSMap;
    // Reading atributes
    bool               dn_normal; // normal de novo
    bool               dn_drive_clogp;
    bool               dn_drive_esol;
    bool               dn_drive_qed;
    bool               dn_drive_sa;
    bool               dn_drive_stereocenters;
    bool               dn_calc_pains;
    // Driving attributes
    double             dn_lower_clogp;
    double             dn_upper_clogp;
    double             dn_clogp_std_dev;
    double             dn_lower_esol;
    double             dn_upper_esol;
    double             dn_esol_std_dev;
    double             dn_lower_qed;
    double             dn_qed_std_dev;
    double             dn_upper_sa;
    double             dn_sa_std_dev;
    int                dn_upper_stereocenter;
    // Attributes for debugging
    bool               dn_drive_verbose;
    bool               dn_save_all_molecules;
    int                dn_start_at_layer;
    #endif


    int                dn_constraint_rot_bon;       // upper bound for # rot bonds
    float              dn_constraint_formal_charge; // maximum absolute charge (- or +)
    //int                dn_constraint_charge_groups; // number of charged groups (- or +)
    int                dn_heur_unmatched_num;       // num of unmatched atoms for H. RMSD pruning
    float              dn_heur_matched_rmsd;        // rmsd of matched atoms for H. RMSD pruning

    // Input parameters: growth control
    int                dn_unique_anchors;           // number of anchors to seed growth
    int                dn_max_grow_layers;          // maximum number of layers
    int                dn_max_root_size;            // control over combinatorics
    int                dn_max_layer_size;           // control over combinatorics
    int                dn_max_current_aps;          // maximum number of aps one frag can have
    int                dn_max_scaffolds_per_layer;  // maximum number of scaffolds you can add at each layer of growth

    // Input parameters: output
    bool               dn_write_checkpoints;        // write checkpoint files at each layer
    bool               dn_write_prune_dump;         // write molecules pruned at each layer
    bool               dn_write_orients;            // write anchor orients from beginning of growth
    bool               dn_write_growth_trees;       // write growth trees
    std::string        dn_output_prefix;            // output filename prefix
    bool               verbose;                     // output verbose statistics
    //GA parameters
    bool               dn_ga_flag;                  // for mutation, make sure 1 valid molecule is produced

    // Vectors of fragments and other important things
    std::vector <Fragment>     scaffolds;
    std::vector <Fragment>     linkers;
    std::vector <Fragment>     sidechains;
    std::vector <Fragment>     scaf_link_sid;       // a combination of the first three
    //std::vector <Fragment>   rigids;
    std::vector <Fragment>     anchors;
    std::vector <TorEnv>       torenv_vector;
    std::vector <FragGraph>    scaf_link_sid_graph; // for the combined fragment vector
    std::vector <Fragment>     tmp_mutants;         // save mutation results to pass to GA

    // Variables for calculating internal energy (same as in conf_gen_ag.h)
    bool               use_internal_energy;    // int energy function superseded by funct in base_score
    int                ie_att_exp;             // attractive VDW exponent (6 by default)
    int                ie_rep_exp;             // repulsive VDW exponent (12 by default, can change in input file)
    float              ie_diel;                // dielectric constant (4.0 by default)
    //float              internal_energy_cutoff;        //BCF internal energy cutoff
    float              ie_cutoff;              //BCF internal energy cutoff

    // Miscellaneous
    int                molecule_counter;       // for counting molecules
    int                growth_tree_index;      // for naming growth trees
    std::vector <int>  bond_tors_vectors;      // for the minimizer

    //JDB - roulette params
    std::vector<double> roulette_vect;
    std::vector<string> torsion_vect;
    std::vector<double> roulette_bin_vect;
    std::string roulette_dummy_replace_1;
    std::string roulette_dummy_replace_2;
    bool tmp_first_check;
    bool tmp_second_check;
    /** Functions **/

    // Read parameters from file, initialize stuff, prepare vectors and molecules
    void            input_parameters( Parameter_Reader & parm );
    void            initialize();
    void            initialize_internal_energy_parms( bool, int, int, float, float );
    void            read_library_anchor( std::vector <Fragment> &, std::string );
    void            read_library( std::vector <Fragment> &, std::string );
    void            read_torenv_table( std::string );
    void            prepare_fragment_graph( std::vector<Fragment> &, std::vector <FragGraph> & );
    void            read_roulette( string roulette_table);
    void            generate_roulette();
    // Main build functions
    void            build_molecules( Master_Score &, Simplex_Minimizer &, AMBER_TYPER &, Orient & );
    void            orient_fragments( std::vector <Fragment> &, Fragment &, Master_Score &, Simplex_Minimizer &,
                                      AMBER_TYPER &, Orient & );
    void            simple_build( Master_Score &, Simplex_Minimizer &, AMBER_TYPER & );
    void            sample_fraglib_exhaustive( Fragment &, int, std::vector <Fragment> &, std::vector <Fragment> &,
                                               Master_Score &, Simplex_Minimizer &, AMBER_TYPER &, bool );
    void            sample_fraglib_rand( Fragment &, int, std::vector <Fragment> &, std::vector <Fragment> &,
                                         Master_Score &, Simplex_Minimizer &, AMBER_TYPER &, bool );
    void            sample_fraglib_graph( Fragment &, int, std::vector <Fragment> &, std::vector <FragGraph> &,
                                          std::vector <Fragment> &,
                                          Master_Score &, Simplex_Minimizer &, AMBER_TYPER &, bool );

    // Functions for attaching two fragments together
    Fragment        combine_fragments( Fragment &, int, int, Fragment, int, int );
    float           calc_cov_radius( std::string atom );
    Fragment        attach( Fragment &, int, int, Fragment &, int, int );
    bool            valid_torenv( Fragment & );
    bool            valid_torenv_multi( Fragment & );
    bool            compare_atom_environments( std::string, std::string );
    bool            compare_dummy_bonds( Fragment, int, int, Fragment &, int, int );
    bool            roulette_valid_torenv ( Fragment & );
    // Functions for sampling torsions, computing energy, and minimizing
    void            sample_minimized_torsions( Fragment &, std::vector <Fragment> &, Master_Score &,
                                               Simplex_Minimizer &, AMBER_TYPER & );
    float           calc_fragment_rmsd( Fragment &, Fragment & );
    void            frag_torsion_drive( Fragment &, std::vector <Fragment> & );
    void            prepare_internal_energy( Fragment &, Master_Score & );

    // Horizontal pruning functions specific to denovo
    bool            mw_cutoff( Fragment &, int ); 
    void            prune_h_rmsd( std::vector <Fragment> & );
    void            calc_mol_wt( DOCKMol & );
    void            calc_rot_bonds( DOCKMol & );
    void            calc_formal_charge( DOCKMol & );
    void            calc_num_HA_HD( DOCKMol & );

    // Miscellaneous functions
    void            activate_mol( DOCKMol & );
    bool            dummy_in_mol( DOCKMol & );
    void            dummy_to_H( Fragment &, int, int );
    void            print_torenv( std::vector <TorEnv> );
    void            print_fraggraph();
    void            convert_H_to_Du( DOCKMol &, int );

    // RDKit-related functions for descriptor-driven processes
    #ifdef BUILD_DOCK_WITH_RDKIT
    bool            clogp_cutoff( Fragment & ); 
    bool            esol_cutoff( Fragment & ); 
    bool            qed_cutoff( Fragment & ); 
    bool            sa_cutoff( Fragment & ); 
    bool            stereocenter_cutoff( Fragment & ); 
    bool            drive_growth( Fragment & );
    #endif

    // Functions that are turned off right now
    //bool            prune_molecular_weight( DOCKMol & );
    //bool            prune_rotatable_bonds( DOCKMol & );


    /** Constructor and Destructor **/

    DN_Build();
    ~DN_Build();



    /** Old Functions **/
/*
    // The functions below here were made obsolete by the new combine_fragments() function
    // (They still work, just slightly slower)
    void                      orient_fragments_old( vector <Fragment> &, Master_Score &,
                                                    Simplex_Minimizer &, AMBER_TYPER &, Orient & );
    Fragment                  combine_fragments_old( Fragment, int, int, Fragment &, int, int );
    DOCKVector                translate_frag( DOCKMol &, int );
    std::pair <float, float>  rotate_frag( DOCKMol &, int );
    void                      unrotate_frag( DOCKMol &, std::pair <float, float> );
    void                      untranslate_frag( DOCKMol &, DOCKVector );
    void                      Rx_func( float, double Rx[3][3] );
    void                      Ry_func( float, double Ry[3][3] );
    void                      Rz_func( float, double Rz[3][3] );
    void                      fix_bond_length( Fragment &, int, Fragment &, int );
    void                      prune_tanimoto( std::vector <Fragment> & );
*/

}; 


// +++++++++++++++++++++++++++++++++++++++++
// Sort Functions
// Used with 'sort()' to sort frag vectors by a given descriptor
bool            fragment_sort(const Fragment &, const Fragment &);
bool            fgpair_sort(const std::pair<float,int> & a, const std::pair<float,int> & b);
bool            size_sort(const Fragment &, const Fragment &);
//Comparator functions JDB
bool            roulette_sort_order_frequencies( const vector<double>& vect1, const vector<double>& vect2);
bool            roulette_sort_order_torsions(const vector<double>& vect1, const vector<double>& vect2);    

double time_seconds();
