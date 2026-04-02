#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "amber_typer.h"
#include "dockmol.h"
#include "fragment.h"
#include "orient.h"
#include "master_score.h"
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


// +++++++++++++++++++++++++++++++++++++++++
// Torsion Environments: A data structure that remembers an atom enviroment (as string) and a vector
// of all the things it can be connected to, as well as corresponding frequencies
class           Tor_Env {
  public:
    std::string                  origin_env;
    std::vector <std::string>    target_envs;
    std::vector <int>            target_freqs;

    // std::vector < std::pair <std::string, int> >  target_envs;

    Tor_Env();
    ~Tor_Env();
};




// +++++++++++++++++++++++++++++++++++++++++
// Class that contains all information about segments needed for mutations
class           SEGMENTS {
  public:

    SEGMENTS() {
        num_aps = 0;
        num_atoms = 0;
        num_bonds = 0;
    };

    int       num_aps;                    // number of attachments for a segmen
    int       num_atoms;                  // number of atoms for a segment
    int       num_bonds;                  // number of atoms for a segment
    INTVec    atoms;                      // atoms in a segment
    INTVec    bonds;                      // bonds in a segment
    INTVec    interseg_bonds;             // bond index associated with aps
};



// +++++++++++++++++++++++++++++++++++++++++
// Genetic algorithm function
class           GA_Recomb {

  public:

    // use these for consistent formatting in output_score_summary
    static const std::string DELIMITER;
    static const int         FLOAT_WIDTH;
    static const int         STRING_WIDTH;                         // standard output header


    /** Variables **/

    // Initial ensemble taken from the dock.in file
    std::string              ga_molecule_file;     
    bool                     ga_utilities;                        // LEP - fight me
	    
    // Filenames of fragment libraries are taken from the dock.in file
    std::string              ga_fraglib_scaffold_file;
    std::string              ga_fraglib_linker_file;
    std::string              ga_fraglib_sidechain_file;
    std::string              ga_fraglib_rigid_file;
    std::string              ga_torenv_table; 

    // User defined input parameters
    // Input parameters: crossover
    bool                     ga_charge_parent_gasteiger;         // gasteiger charger
    bool                     ga_calc_parent_pairwise_tan;        // pairwise tam
    bool                     ga_tan_similarity_prune;            // pairwise tan prune
    bool                     ga_tan_best_first_clustering;	 // tan clsutering unfinished
    bool                     ga_calc_parent_pairwise_hms;        // pairwise hms

    int                      max_generations;                    // max number of generations  
    bool                     ga_xover_on;                        // turn on/off crossover
    bool                     ga_xover_sampling_method_rand;      // crossover sampling method
    int                      ga_xover_max;                       // max number of offspring generated
    float                    ga_bond_tolerance;                  // user-specified cutoff for allowable atom sq dist
    float                    ga_angle_cutoff;                    // user-specified cutoff for bond angle
    bool                     ga_check_overlap;                   // check the molecules that overlap
    bool                     ga_check_only;                      // no molecule generation

    // Input parameters: GA mutation
    bool                     ga_mutations;                       // turn on mutations?
    bool                     ga_mutate_addition;                 // turn on add
    bool                     ga_mutate_deletion;                 // turn on del
    bool                     ga_mutate_substitution;             // turn on sub
    bool                     ga_mutate_replacement;              // turn on repl 
    bool                     ga_mutate_parents;                  // mutate parents bool 
    float                    ga_pmut_rate;                       // parent mutation rate - for random mutation sampling
    float                    ga_omut_rate;                       // offspring mutation rate - for random mutation sampling
    int                      ga_max_mut_cycles;                  // max mutation attemps - for random mutation sampling

    // Input parameters: DN growth
    bool                     ga_use_dn_roulette;                 // turn dn roulette on or off  
    std::string              ga_mut_sampling_method;             // type of sampling (ex, rand, graph)
    bool                     ga_mut_sampling_method_rand;   
    bool                     ga_mut_sampling_method_graph;
    int                      ga_num_random_picks;                // number of random picks for random sampling
    int                      ga_max_root_size;                   // max root size
    int                      ga_graph_max_picks;                 // graph max picks
    int                      ga_graph_breadth;              
    int                      ga_graph_depth;
    float                    ga_graph_temperature;               // initial temp for graph

    // Input parameters: Pruning cutoffs
    bool                     ga_use_torenv_table;                // not used anymore
    //std::string              ga_torenv_table                     // the path to the torenv table
    int                      ga_energy_cutoff;                   // upper bounds for energy pruning
    int                      ga_heur_unmatched_num;              // num of unmatched atoms for H. RMSD pruning
    float                    ga_heur_matched_rmsd;               // rmsd of matched atoms for H. RMSD pruning
    int                      ga_constraint_mol_wt;               // upper bound for mol wt
    int                      ga_constraint_rot_bon;              // upper bound for # rot bonds
    int                      ga_constraint_H_accept;             // upper bound for # hydrogen acceptors
    int                      ga_constraint_H_donor;              // upper bound for # hydrogen donors
    float                    ga_constraint_formal_charge;        // upper and lower bounds for formal charge with 0.1
    int                      ga_user_constraint_formal_charge;   // upper and lower bounds for formal charge 
    // Input parameters: Selection Methods
    int                      ga_ensemble_size;                   // number of survivors to carry to next generation
    std::string              ga_selection_method;                // type of selection ( elitism, tournament, roulette, metropolis )
    bool                     ga_selection_method_elitism;        // carry only the top X% of parents to next generation
    bool                     ga_selection_method_tournament;     // pick the top scoring molecule of the pair
    bool                     ga_selection_method_roulette;       // pick molecules based on their fitness score
    bool                     ga_selection_method_sus;            // pick molecules based on their fitness score - TODO
    bool                     ga_selection_method_metropolis;     // pick based on MCM - TODO
    bool                     ga_elitism_combined;                // if parents + offspring ensembles are combined
    std::string              ga_elitism_option;                  // type of selection ( percent or number )
    float                    ga_elitism_percent;                 // top percent of ensemble to be carried
    int                      ga_elitism_number;                  // top scored parent to next generation
    bool                     ga_tournament_p_vs_c;               // Select the top scored mol between parents and offspring
    bool                     ga_roulette_separate;               // Select the top scored mol between parents and offspring
    bool                     ga_sus_separate;                    // Select the top scored mol between parents and offspring
    float                    ga_mc_initial_temp;                 // Metropolis initial temperature
    float                    ga_mc_rate_constant;                // Temperature modulating constant  Tf = Ti(1-R)^(generation) 
    float                    ga_mc_k;                            // Temperature constant [ P=exp(-(Ff-Fi)/(k*T)) ]
    bool                     ga_metropolis_separate;             
    bool                     ga_niching;                         // Selection niching method
    std::string              ga_niching_method  ;                // type of selection ( elitism, tournament, roulette, metropolis )
    bool                     ga_niche_sharing;                   // Selection niching method - fitness sharing
    //bool                     ga_niche_crowding;                  // Selection niching method - crowding distance
    bool                     ga_selection_extinction;            // Turn on or off extinction
    int                      ga_extinction_start;                // Starting gen to turn on
    int                      ga_extinction_frequency;            // Number of Generations to turn on
    int                      ga_extinction_duration;             // Number of Generations to remain on
    int                      ga_extinction_keep;                 // Number of Generations to remain on
    bool                     ga_extinction_on;                   // Signifies if extinction is currently on
    bool                     ga_extinction_switch_selection;     // Change mutation type
    std::string              ga_secondary_selection_method;      // type of selection ( elitism, tournament, roulette, metropolis )
    bool                     ga_secondary_selection_method_elitism;        // carry only the top X% of parents to next generation
    bool                     ga_secondary_selection_method_tournament;     // pick the top scoring molecule of the pair
    bool                     ga_secondary_selection_method_roulette;       // pick molecules based on their fitness score
    bool                     ga_secondary_selection_method_sus;            // pick molecules based on their fitness score
    bool                     ga_secondary_selection_method_metropolis;     // pick based on MCM

    std::string              ga_name_identifier; // 2021.01.11 - LEP naming prefix for each molecule
    std::string              ga_output_prefix; //2018.01.23 - LEP - user defined output name
    std::string              prefix; //2018.06.20 - LEP - var to name initial generation
    int                      ga_no_xover_exit_gen; //2018.02.24 - YZ - user defined parameter for when to kill a ga
                                                   //run when there is no crossover
    bool                     verbose;              //output verbose statistics - LEP 2018.02.27

    // Global vectors for fragment libraries
    std::vector  <Fragment>  tmp_scaffolds;
    std::vector  <Fragment>  tmp_linkers;
    std::vector  <Fragment>  tmp_sidechains;
    //std::vector  <Fragment>  rigids;
    
    //LEP stuff for ga debug
    std:: vector <DOCKMol>   shitty_molecules_go_here;           // array filled in distance calculator
    std:: vector <DOCKMol>   good_molecules_go_here;
    std:: vector <DOCKMol>   prunemol;

    // Global vectors for breeding 
    std::vector  <DOCKMol>   parents;                            // DOCKMol vector holds the initial ensemble
    std::vector  <DOCKMol>   xover_parents;                      // DOCKMol vector holds unique parents that overlap
    std::vector  <DOCKMol>   tmp_parents;                        // DOCKMol placeholder for selection
    std::vector  <DOCKMol>   saved_parents;                      // DOCKMol placeholder of unmanipulated molecules
    std::vector  <DOCKMol>   mutated_parents;                    // DOCKMol vector of mutated parent molecules
    std::vector  <DOCKMol>   children;                           // DOCKMol vector of crossover results
    std::vector  <DOCKMol>   tmp_children;                       // DOCKMol placeholder for random crossover
    std::vector  <DOCKMol>   pruned_children;                    // DOCKMol vector of energy scored conformers to be pruned
    std::vector  <DOCKMol>   scored_generation;                  // DOCKMol vector of final scored, pruned ensemble 
    std::vector  <DOCKMol>   mutants;                            // DOCKMol vector of mutated offspring 

    std::vector  <SEGMENTS>  orig_segments;                      // Hold segment index of each molecule

    std::vector  <Tor_Env>   torenv_vector;                      // Torenv - should switch to DN version 
    std::vector  <int>       bond_tors_vectors;                  // holds mol torsions for minimization

    DOCKMol tmp_parent1;                                         // DOCKMols holders
    DOCKMol tmp_parent2;
    DOCKMol new_tmp_parent1;                                     // LEP - fix
    DOCKMol new_tmp_parent2;                                     // LEP - fix

    std::vector <int> seg_exclude_indices;                       // JDB - tracking excluded segments for DNM
    bool              dnm_encountered;                           // JDB boolean for tracking encountering DNM in deletion
    bool              dnm_enabled;                               // JDB boolean for tracking if initial ensemble had DNM
    // Variables for calculating internal energy (same as in conf_gen_ag.h)
    bool                     use_internal_energy;                // int energy function superseded by funct in base_score
    int                      ie_att_exp;                         // attractive VDW exponent (6 by default)
    int                      ie_rep_exp;                         // repulsive VDW exponent (12 by default, can change in input file)
    float                    ie_diel;                            // dielectric constant (4.0 by default)
    float                    ie_cutoff;                          // BCF internal energy cutoff


    // Variables for mutations
    int                      num_segments;                       // number of segments per molecule
    int                      num_segs_removed;                   // used to determine max layers ofr DN growth
    INTVec                   no_seg_bonds;                       // list of intersegment bonds
    std::vector <INTPair>    bond_btwn_seg;                      // segment index associated with the bonds
    INTVec                   terminal_seg;                       // terminal segment index


    // Vectors for success rates
    INTVec                   success_add;                        // bool vectors for successful mutations
    INTVec                   success_del;
    INTVec                   success_sub;
    INTVec                   success_replace;


    // Variable for selection
    float                    new_temp;                           // updated Metropolis temp; Tf = Ti(1-R)^(generation)


    // Miscellaneous variables
    bool                     used;                               // used was added to dockmol.cpp 
    


    /** Functions **/

    // Read parameters from file, initialize stuff, prepare vectors and molecules
    void                     input_parameters( Parameter_Reader & parm );                          // initialize GA input parameters
    void                     input_parameters_selection( Parameter_Reader & parm );                // initialize GA selection method params
    void                     initialize();                                                         // initialize all DOCK params
    void                     initialize_fraglib( std::vector <Fragment> & );                       // initialize frag library for GA 
    void                     initialize_internal_energy_parms( bool, int, int, float, float );     // initialize internal energy
    void                     read_library( std::vector <DOCKMol> &, std::string );                 // read library
    void                     activate_half_mol( DOCKMol &, int, int, DOCKMol &, int, int, bool );  // turn on subset of atom/bond active flags
    Fragment                 mol_to_frag ( DOCKMol &);                                             // convert dockmol to fragment type
    void                     prepare_mol_torenv ( DOCKMol &);                                      // prepare mol for torenv check
    void                     prepare_mut_torenv ( DOCKMol &, int );                                // prepare mol by num for toenv check
    void                     prepare_torenv_indices ( Fragment &);                                 // prepare toenv indicies 
    void                     deactivate_mol ( DOCKMol &);                                          // turn off active atom/bond flags
    void                     activate_mol ( DOCKMol &);                                            // activate atom/bond flags
    void                     activate_vector( std::vector <DOCKMol> & );                           // activate a vector of mols

    // Main functions for breeding
    void                     max_breeding( Master_Score &, Simplex_Minimizer &, AMBER_TYPER &, Orient & );        // main function
    void                     check_exhaustive( Master_Score &, Simplex_Minimizer &, AMBER_TYPER &, Orient & );    // check overlapping bonds only
    void                     breeding_exhaustive( Master_Score &, Simplex_Minimizer &, AMBER_TYPER &, Orient &, int ); // exhaustive xover & mut
    void                     breeding_rand( Master_Score &, Simplex_Minimizer &, AMBER_TYPER &, Orient &, int );       // rand xover & mut
 
  
    // Functions for combining two fragments for crossover
    bool                     similarity_compare( DOCKMol &, int, int, DOCKMol &, int, int, 
                                                 bool, Master_Score &, Simplex_Minimizer &, AMBER_TYPER & );  // compare similarity of two segments
    void                     recomb_mols( DOCKMol &, int, int, DOCKMol &, int, int, std::vector <DOCKMol> & );// main xover function
    void                     recomb_mols_xover( DOCKMol &, int, int, DOCKMol &, int, int, std::vector <DOCKMol> & ); // main xover function, w/ DNA JDB
    void                     rotate( float, float, DOCKVector, DOCKMol &, DOCKVector);                        // rotate/align segment for attachment
    float                    calc_cov_radius( std::string atom );                                             // calc covalent radius 
    DOCKMol                  attach( DOCKMol & , int, DOCKMol &, int);                                        // create bond
    void                     switch_active_halves( DOCKMol &, int, int, DOCKMol &, int, int, 
                                                   std::vector <DOCKMol> & );                                 // switch active atoms/bonds
    std::vector <std::pair <int, int> > replace_combine_mols( DOCKMol &, int, int, 
                                                              DOCKMol &, int, int,
                                                              std::vector < std::pair <int,int> > &, 
                                                              std::vector <DOCKMol> &);                       // translate & rotate seg

    // Functions for mutations
    void                    master_mut_exhaustive( std::vector <DOCKMol> &, std::vector <DOCKMol> &, 
                                                   std::vector <DOCKMol> &, bool, Master_Score &, 
                                                   Simplex_Minimizer &, AMBER_TYPER &, Orient &, int );            // exhaustive mutations
    void                    master_mut_rand( std::vector <DOCKMol> &, std::vector <DOCKMol> &, 
                                             std::vector <DOCKMol> &, bool, Master_Score &, 
                                             Simplex_Minimizer &, AMBER_TYPER &, Orient &, int );                  // random mutations
    void                    segment_id(DOCKMol &, AMBER_TYPER &);                                             // get segment information (ag)
    void                    mutation_selection(DOCKMol &, std::vector <DOCKMol> &, Master_Score &, 
                                               Simplex_Minimizer &, AMBER_TYPER &, Orient &, int, int );                // select mut type
    void                    deletion_mutation(DOCKMol &, INTVec, INTVec, INTVec, bool, 
                                              Master_Score &, Simplex_Minimizer &, AMBER_TYPER &);            // deletion
    void                    additions(DOCKMol &, int, Master_Score &, Simplex_Minimizer &, 
                                      AMBER_TYPER &, Orient &, bool, int, int);                                         // additions
    void                    replacement (DOCKMol &, int, INTVec, INTVec, INTVec, Master_Score &, 
                                         Simplex_Minimizer &, AMBER_TYPER &, Orient &);                       // main replacements
    void                    prepare_replacement_segment (Fragment &, int, INTVec, INTVec, 
                                                         INTVec, AMBER_TYPER & );                             // prepare for repl
    void                    rand_replacement(Fragment &, int, std::vector <Fragment> &, 
                                             Master_Score &, Simplex_Minimizer &, AMBER_TYPER &, Orient &);   // select & attach new frag 
    std::vector <std::vector <int> > compare_bond_types ( Fragment &, Fragment & );                           // compare bond types
    std::vector <int>       allowed_bond_combos ( Fragment &, Fragment &, std::vector < std::vector <int> > & );// save same bond types
    std::vector <DOCKMol>   orient_frag_to_ref ( Fragment &, Fragment &, int, 
                                                 Master_Score &, AMBER_TYPER &, Orient &);                    // orient frag to ligand
    std::vector <std::pair <int, float> > overlapping_aps ( Fragment &, int, Fragment &, 
                                                            DOCKMol &, AMBER_TYPER &);                        // list of overlapping aps
    std::vector <DOCKMol>   replace_scaffold ( Fragment &, int, int, bool, bool, Fragment &, 
                                               std::vector <DOCKMol> &, AMBER_TYPER &, Master_Score &);       // 
    std::vector <float>     calc_norm_ref ( Fragment &, int); 
    bool                    compare_norm ( Fragment &, int, std::vector <float>, Fragment &, DOCKMol & ); 
    std::vector <float>     calc_norm ( DOCKMol &, std::vector <int> ); 
    float                   exit_vector(DOCKMol &, int, int, int); 
    int                     compare_rmsd_bt ( Fragment &, Fragment &, bool, int, std::vector <std::pair <int, float> > & );
    bool                    rand_H_to_Du (DOCKMol &, int );
    void                    H_to_Du (DOCKMol &, int );
    void                    add_H(DOCKMol &, int, int, bool, Master_Score &, Simplex_Minimizer &, AMBER_TYPER & );

    // Functions for sampling torsions, computing energy, and minimizing
    void                    score_parents( std::vector <DOCKMol> &, Master_Score &, AMBER_TYPER & );
    void                    unique_parents( std::vector <DOCKMol> &, Master_Score &, Simplex_Minimizer &, AMBER_TYPER & ); // Eliminate redundancy based on hungarian
    void                    score_children(int, Master_Score &, Simplex_Minimizer &, AMBER_TYPER & ); //REMOVE
    void                    minimize_children( DOCKMol &, std::vector <DOCKMol> &, Master_Score &, Simplex_Minimizer &, AMBER_TYPER & ); //
    void                    prepare_internal_energy( DOCKMol &, Master_Score & );

   
    // Pruning and fitness functions
    void                    uniqueness_prune( std::vector <DOCKMol> &, std::vector <DOCKMol> & , Master_Score &, AMBER_TYPER & );
    void                    uniqueness_prune_mut( std::vector <DOCKMol> &, std::vector <DOCKMol> & , Master_Score &, AMBER_TYPER & );
    void                    fitness_pruning( std::vector <DOCKMol> &, std::vector <DOCKMol> &, Master_Score &, Simplex_Minimizer &, AMBER_TYPER & );
    void                    hard_filter( std::vector <DOCKMol> & ); // COMMENT
    bool                    hard_filter_mol( DOCKMol & ); // COMMENT
    void                    calc_descriptors( DOCKMol & ); // COMMENT
    void                    calc_mol_wt( DOCKMol & );
    void                    calc_rot_bonds ( DOCKMol & );
    void                    num_HA_HD( DOCKMol & );
    void                    calc_formal_charge( DOCKMol & );


    // Selection functions
    void                   selection_method( Master_Score &, AMBER_TYPER & );
    void                   selection_elite( std::vector <DOCKMol> &, Master_Score &, AMBER_TYPER & );
    void                   selection_tournament( std::vector <DOCKMol> &, Master_Score &, AMBER_TYPER & );
    void                   tournament_trad( std::vector <DOCKMol> &, int, std::vector <DOCKMol> &, int, std::vector <DOCKMol> & );
    void                   tournament_niche( std::vector <DOCKMol> &, int, std::vector <DOCKMol> &, int, std::vector <DOCKMol> & );
    void                   selection_roulette( std::vector <DOCKMol> &, Master_Score &, AMBER_TYPER & );
    void                   selection_sus( std::vector <DOCKMol> &, Master_Score &, AMBER_TYPER & );
    void                   roulette( std::vector <DOCKMol> &, int, bool );
    void                   selection_metropolis( std::vector <DOCKMol> &, Master_Score &, AMBER_TYPER & );
    void                   mc_metropolis( std::vector <DOCKMol> &, int );
    void                   true_mc_metropolis( std::vector <DOCKMol> &, int );


    // Multiobjective function
    void                   niche_sharing( std::vector <DOCKMol> &, Master_Score &);
    void                   niche_crowding( std::vector <DOCKMol> &, Master_Score &);
    void                   score_scaling( std::vector <DOCKMol> &, int, Master_Score &);
    std::vector <std::vector <float> > distance_matrix( std::vector <DOCKMol> &, int, Master_Score &);
    float                  niche_radius( std::vector <DOCKMol> & mol, std::vector <std::vector <float> >);
    void                   sharing_function( std::vector <DOCKMol> &, int, std::vector <std::vector <float> >, float, Master_Score &);
    void                   increment_rank( std::vector <DOCKMol> & );
    void                   crowding_dist( std::vector <DOCKMol> &, int, Master_Score & );


    // +++++++++++++++++++++++++++++++++++
    // Check functions
    void                    is_active( DOCKMol & );
    void                    is_active_vec( std::vector <DOCKMol> & );
    void                    print_torenv( std::vector <Tor_Env>  );


    //Counters
    double                  time_seconds();
    /** Constructor and Destructor **/

    // Naming function
    void                   naming_function ( DOCKMol & , int, int);
    void                   renumber_atom_numbering( DOCKMol & , bool );    
    //Utilities
    void                   calc_pairwise_distance ( DOCKMol & ); //JDB
    void		   pairwise_tanimoto_calc ( std::vector <DOCKMol> & );
    void                   pairwise_hms_calc ( std::vector <DOCKMol> & );

    // Because Rob wanted it LEP
    void                   set_DNM_bools(DOCKMol &);
    void                   set_DNM_bools_start( std::vector <DOCKMol> & ); //Start of run DNM check
    /** Archived Functions and Related Variables **/
    //void                     read_torenv_table( std::string ); //
    //int                      molecule_counter;                   // for counting molecules
    //bool                     valid_torenv( DOCKMol  );
    //bool                     compare_atom_environments( std::string, std::string );
    // bool                    compare_dummy_bonds( Fragment, int, int, Fragment &, int, int );
    //void                    prune_children( Master_Score &, Simplex_Minimizer &, AMBER_TYPER & );//REMOVE
    //float                   calc_rmsd( DOCKMol &, DOCKMol & ); //REMOVE
    //void                    child_torsion_drive( int,std::vector <DOCKMol> & );
    //void                    rmsd_prune( Master_Score &, Simplex_Minimizer &, AMBER_TYPER & );
    #ifdef BUILD_DOCK_WITH_RDKIT
    std::string        ga_sa_fraglib_path;
    std::string        ga_PAINS_path;
    std::map<unsigned int, double> ga_fragMap;
    std::map<std::string, std::string> ga_PAINSMap;
    #endif
    GA_Recomb();
    ~GA_Recomb();

}; 

// ++++++++++++++++++++++++++++++++++++
// Sort Functions
// Used with 'sort()' to sort frag vectors by 'current_score'
    bool                   mol_sort(DOCKMol , DOCKMol );
    //bool            fragment_sort(const Fragment &, const Fragment &);
    bool                   grid_sort(DOCKMol , DOCKMol );
    bool                   mg_sort(DOCKMol , DOCKMol );
    bool                   cont_sort(DOCKMol , DOCKMol );
    bool                   fps_sort(DOCKMol , DOCKMol );
    bool                   ph4_sort(DOCKMol , DOCKMol );
    bool                   tan_sort(DOCKMol , DOCKMol );
    bool                   hun_sort(DOCKMol , DOCKMol );
    bool                   vol_sort(DOCKMol , DOCKMol );
    bool                   rank_sort(DOCKMol , DOCKMol );
    bool                   fitness_sort(DOCKMol , DOCKMol );
    bool                   crowd_sort(DOCKMol , DOCKMol );
    bool                   compare_rank(DOCKMol , DOCKMol );
    
//
