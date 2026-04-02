//
#ifndef SIMPLEX_H
#define SIMPLEX_H 

#include <vector>
#include "utils.h"  // INTVec, TORSION
class AMBER_TYPER;
class Base_Score;
class DOCKMol;
class Master_Score;
class Parameter_Reader;
class           Simplex_Minimizer {

  public:
    bool            minimize_ligand;
    std::vector<TORSION> torsions;
    INTVec          torsion_scale_factors;
    INTVec          bond_vectors;

    // state variables
    int             ga_gen;
    bool            simplex_ga_flag;
    int             current_cycle;
    int             random_seed;
    int 	        current_layer;
    int             num_layers;	
    // tethering parameters 
    bool            restrained_min; 
    float           coefficient_restraint;

    // simplex parameters
    bool            advanced_min_params;
    int             max_iterations;
    int             torsion_iterations;
    int             max_cycles;
    float           initial_score_converge;
    float           score_converge;
    float           cycle_converge;
    float           trans_step_size;
    float           rot_step_size;
    float           tors_step_size;

    // anchor minimization parameters
    bool            use_min_rigid_anchor;
    int             anchor_min_max_iterations;
    int             anchor_min_max_cycles;
    float           anchor_min_score_converge;
    float           anchor_min_cycle_converge;
    float           anchor_min_trans_step_size;
    float           anchor_min_rot_step_size;
    float           anchor_min_tors_step_size;

    // flex search minimization parameters
    bool            use_min_flex_growth_ramp;
    bool            use_min_flex_growth;
    int             flex_min_max_iterations;
    int             flex_min_torsion_iterations;   // # of steps of torsion pre-min during growth
    int             flex_min_max_cycles;
    float           flex_min_score_converge;
    float           flex_min_cycle_converge;
    float           flex_min_trans_step_size;
    float           flex_min_rot_step_size;
    float           flex_min_tors_step_size;
    //bool            flex_min_add_internal;     // no longer used. internal nrg function in base_score 
                                                 // returns zero when internal energy is switched off
    // flex mine ramp parameters
    bool            use_min_ramp_flex_growth;
    int             flex_min_ramp_max_iterations;
    int             flex_min_ramp_torsion_iterations;   // # of steps of torsion pre-min during growth
    int             flex_min_ramp_max_cycles;
    float           flex_min_ramp_score_converge;
    float           flex_min_ramp_cycle_converge;
    float           flex_min_ramp_trans_step_size;
    float           flex_min_ramp_rot_step_size;
    float           flex_min_ramp_tors_step_size;

/**/    // final minimization parameters
    bool            final_min;
    int             final_min_max_iterations;
    int             final_min_max_cycles;
    float           final_min_score_converge;
    float           final_min_cycle_converge;
    float           final_min_trans_step_size;
    float           final_min_rot_step_size;
    float           final_min_tors_step_size;
    float	    final_min_rep_radius_scale;
/**/

    // parameters for minimization with secondary score
    bool            secondary_min_pose;
    int             secondary_min_max_iterations;
    int             secondary_min_max_cycles;
    float           secondary_min_score_converge;
    float           secondary_min_cycle_converge;
    float           secondary_min_trans_step_size;
    float           secondary_min_rot_step_size;
    float           secondary_min_tors_step_size;
    bool            secondary_advanced_min_params;

    // functions
    void            input_parameters(Parameter_Reader & parm,
                                     bool flexible_ligand, bool genetic_algorithm, 
                                     bool denovo_design, Master_Score &);
    void            initialize();
    // minimize iterates over cycles and calls simplex_minimize
    void            minimize(Base_Score &, DOCKMol &, FLOATVec &, int, float, int, float, float, float, float);
    // simplex_minimize perform actual minimization.
    float           simplex_minimize(Base_Score &, DOCKMol &, FLOATVec &, int, float, float, float, float);
    bool            simplex_score(Base_Score &, DOCKMol &, DOCKMol &, FLOATVec &, float, float, float);
    void            scale_simplex_vector(FLOATVec &, FLOATVec &, float, float, float);
    void            vector_to_dockmol(DOCKMol &, FLOATVec &);

    void            id_torsions(DOCKMol &, FLOATVec &);
    //Wrapper functions that call minimize.
    void            minimize_rigid_anchor(DOCKMol &, Master_Score &);  //only for ag anchors
    void            minimize_flexible_growth(DOCKMol &, Master_Score &, INTVec & ); // for growth
    void            minimize_flexible_ramp_growth(DOCKMol &, Master_Score &, INTVec &, int current_layer, int num_layers ); // for rampgrowth 	
    void            minimize_final_pose(DOCKMol &, Master_Score &, AMBER_TYPER &);  
    void            minimize_pose_final_min(DOCKMol &, Master_Score &);  
    void            secondary_minimize_pose(DOCKMol &, Master_Score &);
    float           calc_active_rmsd2(DOCKMol &, DOCKMol &);


};

#endif  // SIMPLEX_H

