#include <iostream>
#include <string.h>
#include <stdlib.h>
#include "master_score.h"
#include "simplex.h"
#include "conf_gen_ag.h"

using namespace std;


/******************************************************/
void
Simplex_Minimizer::input_parameters(Parameter_Reader & parm,
                                    bool flexible_ligand, bool genetic_algorithm, bool denovo_design, Master_Score & score)
{

    advanced_min_params = false;
    use_min_rigid_anchor = false;
    use_min_flex_growth = false;
    use_min_flex_growth_ramp = false;
    final_min = false;
    secondary_min_pose = false;

    minimize_ligand = false;
    random_seed = 0;

    // restrained minimization parameters
    restrained_min = false;
    coefficient_restraint = 0.0;


    if (score.primary_min) {

        cout << "\nSimplex Minimization Parameters" << endl;
        cout <<
            "------------------------------------------------------------------------------------------"
            << endl;

        minimize_ligand =
            (parm.query_param("minimize_ligand", "yes", "yes no") ==
             "yes") ? true : false;

        if (minimize_ligand) {

            // if anchor and grow is flagged, check which portions of molecule
            // to minimize
            // PAK added logic and queries so user can choose if they want the simplex ramp
            if (flexible_ligand) {
                use_min_rigid_anchor = (parm.query_param("minimize_anchor", "yes", "yes no") == "yes") ? true : false;
                use_min_flex_growth = (parm.  query_param("minimize_flexible_growth", "yes", "yes no") == "yes") ? true : false;
		        advanced_min_params = (parm.  query_param("use_advanced_simplex_parameters", "no", "yes no") == "yes") ? true : false;
                // LEP - ga and dn not ready for simplex ramp
                if (!genetic_algorithm && !denovo_design){
                    use_min_flex_growth_ramp = (parm.  query_param("minimize_flexible_growth_ramp", "yes", "yes no") == "yes") ? true : false;
                }

            }
            // if docking is rigid or if the same parameters will be used for
            // all levels of primary minimization
            if (!advanced_min_params) {
                if (!flexible_ligand || (!use_min_rigid_anchor && !use_min_flex_growth)) {
                    const char simplex_max_it[] = "simplex_max_iterations";
                    max_iterations = atoi(parm.query_param(simplex_max_it, "1000").c_str());
                    if (max_iterations < 0) { // this is so that we can run only tors_premin
                        cout << "ERROR:  Parameter \"" << simplex_max_it
                            << "\" must be an integer greater than or equal to zero."
                            << endl
                            << "Program will terminate."
                            << endl;
                        exit(0);
                    }
                   
                    // ask for # of iterations of the torsion premin 
                    torsion_iterations = atoi(parm.query_param("simplex_tors_premin_iterations", "0").  c_str());
                    if (torsion_iterations < 0) {
                        cout << "ERROR:  simplex_tors_premin_iterations cannot be negative. Program will terminate." << endl;
                        exit(0);
                    }

                }
                max_cycles =
                    atoi(parm.query_param("simplex_max_cycles", "1").c_str());
                if (max_cycles <= 0) {
                    cout <<
                        "ERROR:  Parameter must be an integer greater than zero.  Program will terminate."
                        << endl;
                    exit(0);
                }
                /// PAK
                if (!use_min_flex_growth_ramp){
                    score_converge =
                        atof(parm.query_param("simplex_score_converge", "0.1").
                             c_str());
                    if (score_converge <= 0.0) {
                        cout <<
                            "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                            << endl;
                        exit(0);
                    }
                }
                if (use_min_flex_growth_ramp){
                    score_converge =
                           atof(parm.query_param("simplex_score_converge", "0.1").
                                c_str());
                       if (score_converge <= 0.0) {
                           cout <<
                               "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                               << endl;
                           exit(0);
                       }
                   initial_score_converge = atof(parm.query_param("simplex_initial_score_coverge","5").c_str());
                       if (initial_score_converge <= score_converge) {
                           cout <<
                               "ERROR:  Parameter must be larger than score converge value.  Program will terminate."
                              << endl;
                           exit(0);
                       }
                 }
                cycle_converge =
                    atof(parm.query_param("simplex_cycle_converge", "1.0").
                         c_str());
                if (cycle_converge <= 0.0) {
                    cout <<
                        "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                        << endl;
                    exit(0);
                }
                trans_step_size =
                    atof(parm.query_param("simplex_trans_step", "1.0").c_str());
                if (trans_step_size <= 0.0) {
                    cout <<
                        "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                        << endl;
                    exit(0);
                }
                rot_step_size =
                    atof(parm.query_param("simplex_rot_step", "0.1").c_str());
                if (rot_step_size <= 0.0) {
                    cout <<
                        "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                        << endl;
                    exit(0);
                }
                tors_step_size =
                    atof(parm.query_param("simplex_tors_step", "10.0").c_str());
                if (tors_step_size <= 0.0) {
                    cout <<
                        "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                        << endl;
                    exit(0);
                }
            }
            // parameters for anchor minimization
            if (use_min_rigid_anchor) {
                anchor_min_max_iterations = atoi(parm.  query_param("simplex_anchor_max_iterations", "500").c_str());
                if (anchor_min_max_iterations <= 0) {
                    cout <<
                        "ERROR:  Parameter must be an integer greater than zero.  Program will terminate."
                        << endl;
                    exit(0);
                }

                if (advanced_min_params) {
                    anchor_min_max_cycles =
                        atoi(parm.query_param("simplex_anchor_max_cycles", "1").
                             c_str());
                    if (anchor_min_max_cycles <= 0) {
                        cout <<
                            "ERROR:  Parameter must be an integer greater than zero.  Program will terminate."
                            << endl;
                        exit(0);
                    }
                    anchor_min_score_converge =
                        atof(parm.
                             query_param("simplex_anchor_score_converge",
                                         "0.1").c_str());
                    if (anchor_min_score_converge <= 0.0) {
                        cout <<
                            "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                            << endl;
                        exit(0);
                    }
                    anchor_min_cycle_converge =
                        atof(parm.
                             query_param("simplex_anchor_cycle_converge",
                                         "1.0").c_str());
                    if (anchor_min_cycle_converge <= 0.0) {
                        cout <<
                            "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                            << endl;
                        exit(0);
                    }
                    anchor_min_trans_step_size =
                        atof(parm.
                             query_param("simplex_anchor_trans_step",
                                         "1.0").c_str());
                    if (anchor_min_trans_step_size <= 0.0) {
                        cout <<
                            "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                            << endl;
                        exit(0);
                    }
                    anchor_min_rot_step_size =
                        atof(parm.query_param("simplex_anchor_rot_step", "0.1").
                             c_str());
                    if (anchor_min_rot_step_size <= 0.0) {
                        cout <<
                            "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                            << endl;
                        exit(0);
                    }
                    anchor_min_tors_step_size =
                        atof(parm.
                             query_param("simplex_anchor_tors_step",
                                         "10.0").c_str());
                    if (anchor_min_tors_step_size <= 0.0) {
                        cout <<
                            "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                            << endl;
                        exit(0);
                    }
                } else {
                    anchor_min_max_cycles = max_cycles;
                    anchor_min_score_converge = score_converge;
                    anchor_min_cycle_converge = cycle_converge;
                    anchor_min_trans_step_size = trans_step_size;
                    anchor_min_rot_step_size = rot_step_size;
                    anchor_min_tors_step_size = tors_step_size;
                }
            }
            // parameters for flexible minimization
            if (use_min_flex_growth) {

                // If ramp is on simplex_grow_max_iterations should be 250, otherwise 500
                if (use_min_flex_growth_ramp) {
                   flex_min_max_iterations = atoi(parm.query_param("simplex_grow_max_iterations", "250").  c_str());
                }
                else {
                   flex_min_max_iterations = atoi(parm.query_param("simplex_grow_max_iterations", "500").  c_str());
                }

                if (flex_min_max_iterations < 0) {
                    cout <<
                        "ERROR:  simplex_grow_max_iterations cannot be negative.  Program will terminate."
                        << endl;
                    exit(0);
                }

                flex_min_torsion_iterations = atoi(parm.query_param("simplex_grow_tors_premin_iterations", "0").  c_str());
                if (flex_min_torsion_iterations < 0) {
                    cout << "ERROR:  simplex_grow_tors_premin_iterations cannot be negative. Program will terminate." << endl;
                    exit(0);
                }

                if (advanced_min_params) {
                    flex_min_max_cycles =
                        atoi(parm.query_param("simplex_grow_max_cycles", "1").
                             c_str());
                    if (flex_min_max_cycles <= 0) {
                        cout <<
                            "ERROR:  Parameters must be an integer greater than zero.  Program will terminate."
                            << endl;
                        exit(0);
                    }
                    flex_min_score_converge =
                        atof(parm.
                             query_param("simplex_grow_score_converge",
                                         "0.1").c_str());
                    if (flex_min_score_converge <= 0.0) {
                        cout <<
                            "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                            << endl;
                        exit(0);
                    }
                    flex_min_cycle_converge =
                        atof(parm.
                             query_param("simplex_grow_cycle_converge",
                                         "1.0").c_str());
                    if (flex_min_cycle_converge <= 0.0) {
                        cout <<
                            "ERROR:  Parameters must be a float greater than zero.  Program will terminate."
                            << endl;
                        exit(0);
                    }
                    flex_min_trans_step_size =
                        atof(parm.query_param("simplex_grow_trans_step", "1.0").
                             c_str());
                    if (flex_min_trans_step_size <= 0.0) {
                        cout <<
                            "ERROR:  Parameters must be a float greater than zero.  Program will terminate."
                            << endl;
                        exit(0);
                    }
                    flex_min_rot_step_size =
                        atof(parm.query_param("simplex_grow_rot_step", "0.1").
                             c_str());
                    if (flex_min_rot_step_size <= 0.0) {
                        cout <<
                            "ERROR:  Parameters must be a float greater than zero.  Program will terminate."
                            << endl;
                        exit(0);
                    }
                    flex_min_tors_step_size =
                        atof(parm.query_param("simplex_grow_tors_step", "10.0").
                             c_str());
                    if (flex_min_tors_step_size <= 0.0) {
                        cout <<
                            "ERROR:  Parameters must be a float greater than zero.  Program will terminate."
                            << endl;
                        exit(0);
                    }
                } else {
                    flex_min_max_cycles = max_cycles;
                    flex_min_score_converge = score_converge;
                    flex_min_cycle_converge = cycle_converge;
                    flex_min_trans_step_size = trans_step_size;
                    flex_min_rot_step_size = rot_step_size;
                    flex_min_tors_step_size = tors_step_size;
                }
            }

        }
        // final_min has been superseded by the option to perform anchor and grow 
        // docking with internal energy at every level of growth
        // option to perform one more round of minimization 
        // TEB put final min back, 2023
        //
        final_min = (parm.query_param("simplex_final_min", "no", "yes no") == "yes") ? true : false;
        if (final_min) {
            final_min_rep_radius_scale = atof(parm.query_param("simplex_final_min_rep_rad_scale", "1").c_str());
            if (final_min_rep_radius_scale <= 0.0) {
                cout <<
                    "ERROR:  Parameter must be a float greater than zero. Program will terminate."
                    << endl;
                exit(0);
            }
            final_min_max_iterations =
                atoi(parm.query_param("simplex_final_max_iterations", "500").
                     c_str());
            if (final_min_max_iterations <= 0) {
                cout <<
                    "ERROR:  Parameters must be an integer greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }
            //if (advanced_min_params) {
                final_min_max_cycles =
                    atoi(parm.query_param("simplex_final_max_cycles", "1").
                         c_str());
                if (final_min_max_cycles <= 0) {
                    cout <<
                        "ERROR:  Parameters must be an integer greater than zero.  Program will terminate."
                        << endl;
                    exit(0);
                }
                final_min_score_converge =
                    atof(parm.
                         query_param("simplex_final_score_converge",
                                     "0.1").c_str());
                if (final_min_score_converge <= 0.0) {
                    cout <<
                        "ERROR:  Parameters must be a float greater than zero.  Program will terminate."
                        << endl;
                    exit(0);
                }
                final_min_cycle_converge =
                    atof(parm.
                         query_param("simplex_final_cycle_converge",
                                     "1.0").c_str());
                if (final_min_cycle_converge <= 0.0) {
                    cout <<
                        "ERROR:  Parameters must be a float greater than zero.  Program will terminate."
                        << endl;
                    exit(0);
                }
                final_min_trans_step_size =
                    atof(parm.
                         query_param("simplex_final_trans_step",
                                     "1.0").c_str());
                if (final_min_trans_step_size <= 0.0) {
                    cout <<
                        "ERROR:  Parameters must be a float greater than zero.  Program will terminate."
                        << endl;
                    exit(0);
                }
                final_min_rot_step_size =
                    atof(parm.query_param("simplex_final_rot_step", "0.1").
                         c_str());
                if (final_min_rot_step_size <= 0.0) {
                    cout <<
                        "ERROR:  Parameters must be a float greater than zero.  Program will terminate."
                        << endl;
                    exit(0);
                }
                final_min_tors_step_size =
                    atof(parm.
                         query_param("simplex_final_tors_step",
                                     "10.0").c_str());
                if (final_min_tors_step_size <= 0.0) {
                    cout <<
                        "ERROR:  Parameters must be a float greater than zero.  Program will terminate."
                        << endl;
                    exit(0);
                }
            //} else {
            //    final_min_max_cycles = max_cycles;
            //    final_min_score_converge = score_converge;
            //    final_min_cycle_converge = cycle_converge;
            //    final_min_trans_step_size = trans_step_size;
            //    final_min_rot_step_size = rot_step_size;
            //    final_min_tors_step_size = tors_step_size;
            //}
        }
      
        if (final_min or minimize_ligand) {  

           random_seed =
               atoi(parm.query_param("simplex_random_seed", "0").c_str());
           //trent balius 2009/08/27
           restrained_min = (parm.query_param("simplex_restraint_min", "no", "yes no") == "yes") ? true : false;
           if (restrained_min) {
               coefficient_restraint =
                     atof(parm.query_param("simplex_coefficient_restraint", "10.0").c_str());
           }
        }
    }
}

/******************************************************/
void
Simplex_Minimizer::initialize()
{
    //cout << "Initializing simplex" << endl;
    srand(random_seed);

}

/******************************************************/
// this function is used only for flexible docking anchor minimization
// not for rigid docking
void
Simplex_Minimizer::minimize_rigid_anchor(DOCKMol & mol, Master_Score & score)
{
    int             i;
    FLOATVec        vertex;

    if (use_min_rigid_anchor) {
        // initialize degrees of freedom as all zeros (rigid DOF only)
        vertex.clear();

        // rigid DOF
        for (i = 0; i < 6; i++)
            vertex.push_back(0.000);

        bond_vectors.clear();
        bond_vectors.resize(mol.num_bonds, -1);

        // mc_premin_override = false;
        minimize(*score.primary_score, mol, vertex, anchor_min_max_cycles,
                 anchor_min_cycle_converge, anchor_min_max_iterations,
                 anchor_min_score_converge, anchor_min_trans_step_size,
                 anchor_min_rot_step_size, anchor_min_tors_step_size);

    }

}

//PAK new function for negative  ramping. this was duplicated from the function from line 566 and then editted. 
void
Simplex_Minimizer::minimize_flexible_ramp_growth(DOCKMol & mol, Master_Score & score, INTVec & bvectors, int current_layer, int num_layers) 
{
     float 	     adjustable_score_converge;
     int             i;
     FLOATVec        vertex;
     //5.0 was used for the b due to Guilherme's preliminary data.  	
     float diff_interval  = initial_score_converge - flex_min_score_converge;
     float adj_unit = diff_interval / (num_layers-1); 
     adj_unit = adj_unit * (float)current_layer;  	
     adjustable_score_converge = initial_score_converge - adj_unit;
 
     if (use_min_flex_growth) {
        vertex.clear();

	for (i = 0; i < 6; i++)
		vertex.push_back(0.000);

  	id_torsions(mol, vertex);
	
  	torsion_scale_factors.resize(torsions.size(), 1);

	bond_vectors.clear();
        bond_vectors = bvectors;
        
        
	minimize(*score.primary_score, mol, vertex, flex_min_max_cycles,
                 flex_min_cycle_converge, flex_min_torsion_iterations,
                 adjustable_score_converge, 0, 0, flex_min_tors_step_size);

	minimize(*score.primary_score, mol, vertex, flex_min_max_cycles,
                 flex_min_cycle_converge, flex_min_max_iterations,
                 adjustable_score_converge, flex_min_trans_step_size,
                 flex_min_rot_step_size, flex_min_tors_step_size);

   }
}
/******************************************************/
void
Simplex_Minimizer::minimize_flexible_growth(DOCKMol & mol, Master_Score & score,
                                            INTVec & bvectors)
{
    int             i;
    FLOATVec        vertex;
    if (use_min_flex_growth) {
        // initialize degrees of freedom as all zeros (all DOF)
        vertex.clear();

        // rigid DOF
        for (i = 0; i < 6; i++)
            vertex.push_back(0.000);

        // flex DOF
        id_torsions(mol, vertex);

        torsion_scale_factors.resize(torsions.size(), 1);

        bond_vectors.clear();
        bond_vectors = bvectors;

        // trent & sudipto 01-05-2009
        // just minimize the torsional degrees of freedom
	minimize(*score.primary_score, mol, vertex, flex_min_max_cycles,
                 flex_min_cycle_converge, flex_min_torsion_iterations,
                 flex_min_score_converge, 0, 0, flex_min_tors_step_size);
	
        // NOTE: same pre-min is done in minimize_final_pose
        // sudipto believes that setting trans_step_size=0 && rot_step_size=0
        // in minimize might be using minimization iterations doing 0 size
        // steps of trans or rotational minimization


        // minimize all degrees of freedom for the remainder of the time
        // ceil means we run at least 1 step unless frac_time = 1
        minimize(*score.primary_score, mol, vertex, flex_min_max_cycles,
                 flex_min_cycle_converge, flex_min_max_iterations,
                 flex_min_cycle_converge, flex_min_trans_step_size,
                 flex_min_rot_step_size, flex_min_tors_step_size);	
		 
    }

}



/******************************************************/
// rigid docking uses this function exclusively for minimization
// no longer used for flexible growth now
// this function is used for all minimization except flexible docking
void
Simplex_Minimizer::minimize_final_pose(DOCKMol & mol, Master_Score & score, AMBER_TYPER &typer)
{
    int             i;
    FLOATVec        vertex;

    //cout << "Entering minimize_final_pose" << endl;
    //cout << minimize_ligand << " " << use_min_rigid_anchor << " " << use_min_flex_growth << endl;

    if (minimize_ligand) {
        // perform minimization if anchor and grow minimization was not called
        if (!use_min_rigid_anchor && !use_min_flex_growth) {
            // initialize degrees of freedom as all zeros (all DOF)
            vertex.clear();

            // rigid DOF
            for (i = 0; i < 6; i++)
                vertex.push_back(0.000);

            // flex DOF
            id_torsions(mol, vertex);

            torsion_scale_factors.resize(torsions.size(), 1);

            bond_vectors.clear();
            bond_vectors.resize(mol.num_bonds, -1);

            // trent & sudipto 03-10-2010
            // if (num of iteration) is zero or num_cycles is zero , 
            // the function minimize() will not call simplex_minimize()

            // trent & sudipto 01-05-2009
            // just minimize the torsional degrees of freedom
            minimize(*score.primary_score, mol, vertex, max_cycles,
                     cycle_converge, torsion_iterations,
                     score_converge, 0, 0, tors_step_size);

            // NOTE: same pre-min is done in minimize_flexible_growth
            // sudipto believes that setting trans_step_size=0 && rot_step_size=0
            // in minimize might be using minimization iterations doing 0 size 
            // steps of trans or rotational minimization

            minimize(*score.primary_score, mol, vertex, max_cycles,
                     cycle_converge, max_iterations, score_converge,
                     trans_step_size, rot_step_size, tors_step_size);

            // compute the final score for the molecule
            score.compute_primary_score(mol);

        }

/*      // this code has been superseded by global use of internal energy minimization
        if (final_min) {
            cout << "In Simplex_Minimizer::minimize_final_pose" << endl;
        // initialize degrees of freedom as all zeros (all DOF)
            vertex.clear();

            // rigid DOF
            for (i = 0; i < 6; i++)
                vertex.push_back(0.000);

            // flex DOF
            id_torsions(mol, vertex);

            torsion_scale_factors.resize(torsions.size(), 1);

            bond_vectors.clear();
            bond_vectors.resize(mol.num_bonds, -1);

         // DTM 11-12-08 - fix bug with bad energies in final min - due to re-init vdw energy
            //score.primary_score->rep_radius_scale = final_min_rep_radius_scale;
            //score.primary_score->init_vdw_energy(typer, 6, 12);

            minimize(*score.primary_score, mol, vertex, final_min_max_cycles,
                     final_min_cycle_converge, final_min_max_iterations,
                     final_min_score_converge, final_min_trans_step_size,
                     final_min_rot_step_size, final_min_tors_step_size);

            // compute the final score for the molecule
            score.compute_primary_score(mol);
        }

*/

    }
}

void
Simplex_Minimizer::minimize_pose_final_min(DOCKMol & mol, Master_Score & score)
{
    int             i;
    FLOATVec        vertex;

    //cout << "Entering minimize_final_pose" << endl;
    //cout << minimize_ligand << " " << use_min_rigid_anchor << " " << use_min_flex_growth << endl;
    //
        if (final_min) {
            cout << "In Simplex_Minimizer::minimize_final_pose" << endl;
        // initialize degrees of freedom as all zeros (all DOF)
            vertex.clear();

            // rigid DOF
            for (i = 0; i < 6; i++)
                vertex.push_back(0.000);

            // flex DOF
            id_torsions(mol, vertex);

            torsion_scale_factors.resize(torsions.size(), 1);

            bond_vectors.clear();
            bond_vectors.resize(mol.num_bonds, -1);

         // DTM 11-12-08 - fix bug with bad energies in final min - due to re-init vdw energy
            //score.primary_score->rep_radius_scale = final_min_rep_radius_scale;
            //score.primary_score->init_vdw_energy(typer, 6, 12);

            minimize(*score.primary_score, mol, vertex, final_min_max_cycles,
                     final_min_cycle_converge, final_min_max_iterations,
                     final_min_score_converge, final_min_trans_step_size,
                     final_min_rot_step_size, final_min_tors_step_size);

            // compute the final score for the molecule
            score.compute_primary_score(mol);
        }


}

/******************************************************/
void
Simplex_Minimizer::secondary_minimize_pose(DOCKMol & mol, Master_Score & score)
{
    int             i;
    FLOATVec        vertex;
        // perform additional round of minimization using secondary scoring
        // function
    if(minimize_ligand){
        if (secondary_min_pose) {

            // initialize degrees of freedom as all zeros (all DOF)
            vertex.clear();

            // rigid DOF
            for (i = 0; i < 6; i++)
                vertex.push_back(0.000);

            // flex DOF
            id_torsions(mol, vertex);

            torsion_scale_factors.resize(torsions.size(), 1);

            bond_vectors.clear();
            bond_vectors.resize(mol.num_bonds, -1);

            minimize(*score.secondary_score, mol, vertex,
                     secondary_min_max_cycles, secondary_min_cycle_converge,
                     secondary_min_max_iterations, secondary_min_score_converge,
                     secondary_min_trans_step_size, secondary_min_rot_step_size,
                     secondary_min_tors_step_size);

            // compute the final score for the molecule
            score.compute_secondary_score(mol);

        }
    }
}
/******************************************************/
void
Simplex_Minimizer::id_torsions(DOCKMol & mol, FLOATVec & vertex)
{
    //cout << "Entering id_torsions ..." << endl;
    int             i,
                    j,
                    max_central,
                    central;
    int             a1,
                    a2,
                    a3,
                    a4;
    int             nbr;
    vector < int   >nbrs;
    TORSION         tmp_torsion;

    BREADTH_SEARCH  bfs;

    torsions.clear();

    // loop over bonds- add flex bonds to torsion list
    for (i = 0; i < mol.num_bonds; i++) {
        if (mol.bond_active_flags[i]) {
            // Fochmod Sep 29, 2014 
            // altered to make minimization decisions only if bond is not a
            // rotor and if the mol2 says that it is a single bond
            //cout<< i << "  "<< mol.amber_bt_id[i] << endl;
            if (mol.bond_is_rotor(i)) {
                //cout<< i << "  "<< mol.bond_types[i] << endl;
                if (mol.bond_types[i] == "1") {
                    //cout<< i << "  "<< mol.bond_types[i] << endl;
                    torsions.push_back(tmp_torsion);
                    torsions[torsions.size() - 1].bond_num = i;
                    vertex.push_back(0.0);
                }
            }
        }
    }

    //cout << "torsions.size() = " << torsions.size() <<endl;
    // ID the inter-segment rot-bonds
    for (i = 0; i < torsions.size(); i++) {

        a2 = mol.bonds_origin_atom[torsions[i].bond_num];
        a3 = mol.bonds_target_atom[torsions[i].bond_num];

        max_central = -1;
        // nbrs = mol.get_atom_neighbors(a2);
        nbrs = mol.neighbor_list[a2];

        for (j = 0; j < nbrs.size(); j++) {
            nbr = nbrs[j];
            if (nbr != a3) {
                central = bfs.get_search_radius(mol, nbr, a2);
                if (central > max_central) {
                    a1 = nbr;
                    max_central = central;
                }
            }
        }

        max_central = -1;
        // nbrs = mol.get_atom_neighbors(a3);
        nbrs = mol.neighbor_list[a3];

        for (j = 0; j < nbrs.size(); j++) {
            nbr = nbrs[j];
            if (nbr != a2) {
                central = bfs.get_search_radius(mol, nbr, a3);
                if (central > max_central) {
                    a4 = nbr;
                    max_central = central;
                }
            }
        }

        torsions[i].atom1 = a1;
        torsions[i].atom2 = a2;
        torsions[i].atom3 = a3;
        torsions[i].atom4 = a4;
    }

}

/******************************************************/
void
Simplex_Minimizer::minimize(Base_Score & score, DOCKMol & mol,
                            FLOATVec & vertex, int max_cycles,
                            float cycle_converge, int max_iterations,
                            float score_converge, float trans_step_size,
                            float rot_step_size, float tors_step_size)
{
    // TEB 2010-03-10
    // only minimize if both max_iterations and max_cycles are more than 0.
    if (max_iterations > 0 && max_cycles > 0 ){

        int             i;
        float           distance;
 
        // initialize the minimization structures
        current_cycle = 0;
        distance = 0.0;
        // loop over simplex cycles
        while ((current_cycle < max_cycles)
               && ((distance > cycle_converge) || (current_cycle == 0))) {
 
            // call simplex minimizer
            simplex_minimize(score, mol, vertex, max_iterations, score_converge,
                             trans_step_size, rot_step_size, tors_step_size);
 
            // compute the distance moved, and re-zero the vertex vector
            distance = 0;
            for (i = 0; i < vertex.size(); i++) {
                distance += vertex[i] * vertex[i];
                vertex[i] = 0.0;
            }
            distance = sqrt(distance) / (float) (current_cycle + 1);
 
            current_cycle++;
        }
    }

}

/******************************************************/
float
Simplex_Minimizer::simplex_minimize(Base_Score & score, DOCKMol & mol,
                                    FLOATVec & vertex, int max_iterations,
                                    float score_converge, float trans_step_size,
                                    float rot_step_size, float tors_step_size)
{
    // This is the function the does the work!! not a wrapper function.


    int             iteration;
    int             size;

    float           delta = 0.0;

    // old variables
    int             i,
                    j,
                    x;
    int             ihi,
                    inhi;
    int             ilo = 0;
    float         **p;
    float          *pr;
    float          *prr;
    float          *pbar;
    float          *y;

    float           ypr = 0;
    float           yprr = 0;
    float           alpha = 1.0;        /* range: 0=no extrap, 1=unit step
                                         * extrap, higher OK */
    float           beta = 0.5; /* range: 0=no contraction, 1=full contraction */
    float           optimum;
    int             replace_flag;       /* flag for whether bad vertex replaced 
                                         */
    // end old vars

    FLOATVec        new_vec;
    DOCKMol         tmp_mol,
                    ref_mol,
                    min_mol;
    DOCKMol         rmsd_ref; // this is used to restrain min to starting position.

    double          Econstraint;

    float          *old_vertex;
    float           temp1,
                    temp2;

    float           diff;
    int             step_count;
    int             fail_count;

    // move to main dock loop
    //srand(random_seed); // reset the seed so that molecule order in the file does not mater. 

    size = vertex.size();

    // allocate arrays
    old_vertex = new float[size];
    memset(old_vertex, '\0', sizeof(float) * size);

    p = new float  *[size + 1];
    memset(p, '\0', sizeof(float *) * (size + 1));
    for (i = 0; i < size + 1; i++) {
        p[i] = new float[size];
        memset(p[i], '\0', sizeof(float) * size);
    }
    y = new float[size + 1];
    memset(y, '\0', sizeof(float) * (size + 1));
    pr = new float[size];
    memset(pr, '\0', sizeof(float) * size);
    prr = new float[size];
    memset(prr, '\0', sizeof(float) * size);
    pbar = new float[size];
    memset(pbar, '\0', sizeof(float) * size);

    // End allocation

    // copy molecules
    copy_molecule(ref_mol, mol);
    copy_molecule(min_mol, mol);
    copy_molecule(tmp_mol, mol);
    copy_molecule(rmsd_ref, mol);

    iteration = 0;

    do {

        // initialize all the simplex points
        if (iteration == 0) {

/***
            if((use_mc_premin)&&(!mc_premin_override)) {

                // Monte Carlo Preminimizer
                // Used to generate initial simplex points

                for(i=0;i<size;i++)
                    old_vertex[i] = vertex[i];

                temp2 = simplex_score(score, ref_mol, tmp_mol, vertex, trans_step_size, rot_step_size, tors_step_size);
                int count_v = 0;
                fail_count = 0;
                step_count = 0;
                bool accept_flag = false;

                do {

                    temp1 = simplex_score(score, ref_mol, tmp_mol, vertex, trans_step_size, rot_step_size, tors_step_size);
                    diff = temp1-temp2;
                    step_count++;

                    if((diff < 0)||( exp(-diff) > r_value )||(step_count > max_steps)) {
                        for(i=0;i<size;i++)
                            old_vertex[i] = vertex[i];

                        temp2 = temp1;

                        if((fail_count > fail_threshold)||(step_count > max_steps))
                            accept_flag = true;
                        fail_count = 0;

                        if(accept_flag) {
                            for(i=0;i<size;i++)
                                p[count_v][i] = old_vertex[i];

                            count_v++;
                        }

                    } else
                        fail_count++;

                    for(i=0;i<size;i++)
                        vertex[i] = old_vertex[i] + 0.5*(((float)rand()/(float)RAND_MAX) - 0.5); // 1.0 rather than 2.0 works better

                } while( count_v < (size+1));

            } else {

                // generate random initial simplex points

                for(i=0;i<size;i++)
                    p[0][i] = vertex[i];

                for(i=1;i<size+1;i++) {
                    for(j=0;j<size;j++) {
                        p[i][j] = p[0][j] + 2.0*(((float)rand()/(float)RAND_MAX) - 0.5);
                    }
                }

            }
***/

            // generate random initial simplex points
            for (i = 0; i < size; i++) {
                p[0][i] = vertex[i];
                //p[0][i] = 0.0;
            }

            for (i = 1; i < size + 1; i++) {
                for (j = 0; j < size; j++) {
                    p[i][j] =
                        //p[0][j] + 2.0 * (((float) rand() / (float) RAND_MAX) -
                        vertex[j] + 2.0 * (((float) rand() / (float) RAND_MAX) -
                                         0.5);
                }
            }

            // score initial simplex points
            for (i = 0; i < size + 1; i++) {
                for (j = 0; j < size; j++){
                    vertex[j] = p[i][j];
                }
                // check if simplex point is valid
                if (simplex_score (score, ref_mol, tmp_mol, vertex, trans_step_size,
                     rot_step_size, tors_step_size)) {
                     
                     if (restrained_min){
                         Econstraint = coefficient_restraint * calc_active_rmsd2(rmsd_ref,tmp_mol);
                         //cout << "rmsd2: " << calc_active_rmsd2(rmsd_ref,tmp_mol) << endl;
                         //cout << "Econstraint: " << Econstraint << endl;
                     } else{ 
                         Econstraint = 0.0;
                     }
                     y[i] = tmp_mol.current_score + tmp_mol.internal_energy
                              + Econstraint;

                } else {

                    // store best scoring vertex
                    for (x = 0; x < size; x++){
                        vertex[x] = p[ilo][x];
                    }
                    //optimum = ref_mol.current_score;
                    if (restrained_min){
                         Econstraint = coefficient_restraint * calc_active_rmsd2(rmsd_ref,ref_mol);
                    } else{ 
                         Econstraint = 0.0;
                    }
                    optimum = ref_mol.current_score + ref_mol.internal_energy //trent balius 2008-12-05
                              + Econstraint;

                    // copy best mol to min_mol, and generate min structure
                    copy_crds(min_mol, ref_mol);
                    scale_simplex_vector(new_vec, vertex, trans_step_size,
                                         rot_step_size, tors_step_size);
                    vector_to_dockmol(min_mol, new_vec);
                    copy_crds(mol, min_mol);

                    // free arrays
                    for (x = 0; x < size + 1; x++) {
                        delete[]p[x];
                        p[x] = NULL;
                    }

                    delete[]p;
                    p = NULL;

                    delete[]y;
                    y = NULL;

                    delete[]pr;
                    pr = NULL;

                    delete[]prr;
                    prr = NULL;

                    delete[]pbar;
                    pbar = NULL;

                    delete[]old_vertex;
                    old_vertex = NULL;

                    delta -= optimum;
                    //cout << "Optimum= " << optimum << endl;
                    return optimum;
                }
            }

            delta = y[0];
        
        } else {

            // Begin a new iteration
            for (i = 0; i < size; i++){
                pbar[i] = 0.0;
            }
            // compute vector ave. of all points except the highest
            for (i = 0; i < size + 1; i++){
                if (i != ihi){
                    for (j = 0; j < size; j++){
                        pbar[j] += p[i][j];
                    }
                }
            }
            // extrapolate by a factor alpha through the face
            for (i = 0; i < size; i++) {
                pbar[i] /= (float) size;
                vertex[i] = pr[i] = (1.0 + alpha) * pbar[i] - alpha * p[ihi][i];
            }

            // evaluate the fxn at the reflected point

            if (simplex_score (score, ref_mol, tmp_mol, vertex, trans_step_size,
                 rot_step_size, tors_step_size)) {

                 if (restrained_min){
                     Econstraint = coefficient_restraint * calc_active_rmsd2(rmsd_ref,tmp_mol);
                 } else{
                     Econstraint = 0.0;
                 }
                 ypr = tmp_mol.current_score + tmp_mol.internal_energy
                              + Econstraint;

            } else {
                // store best scoring vertex
                for (x = 0; x < size; x++){
                    vertex[x] = p[ilo][x];
                }
                //optimum = ref_mol.current_score;
                if (restrained_min){
                     Econstraint = coefficient_restraint * calc_active_rmsd2(rmsd_ref,ref_mol);
                } else{ 
                     Econstraint = 0.0;
                }
                optimum = ref_mol.current_score + ref_mol.internal_energy //trent balius 2008-12-05
                              + Econstraint;

                // copy best mol to min_mol, and generate min structure
                copy_crds(min_mol, ref_mol);
                scale_simplex_vector(new_vec, vertex, trans_step_size,
                                     rot_step_size, tors_step_size);
                vector_to_dockmol(min_mol, new_vec);
                copy_crds(mol, min_mol);

                // free arrays
                for (x = 0; x < size + 1; x++) {
                    delete[]p[x];
                    p[x] = NULL;
                }

                delete[]p;
                p = NULL;

                delete[]y;
                y = NULL;

                delete[]pr;
                pr = NULL;

                delete[]prr;
                prr = NULL;

                delete[]pbar;
                pbar = NULL;

                delete[]pbar;
                pbar = NULL;

                delete[]old_vertex;
                old_vertex = NULL;

                delta -= optimum;
                //cout << "Optimum= " << optimum << endl;
                return optimum;
            }


            if (ypr <= y[ilo]) {

                // Gives a better result than the best point, so try
                // extrapolation by alpha
                for (i = 0; i < size; i++){
                    vertex[i] = prr[i] =
                        (1.0 + alpha) * pr[i] - alpha * pbar[i];
                }
                // check if new point is valid
                if (simplex_score
                    (score, ref_mol, tmp_mol, vertex, trans_step_size,
                     rot_step_size, tors_step_size)) {
                     if (restrained_min){
                          Econstraint = coefficient_restraint * calc_active_rmsd2(rmsd_ref,tmp_mol);
                     } else{ 
                          Econstraint = 0.0;
                     }
                     yprr = tmp_mol.current_score + tmp_mol.internal_energy
                              + Econstraint;

                } else {

                    // store best scoring vertex
                    for (x = 0; x < size; x++){
                        vertex[x] = p[ilo][x];
                    }
                    //optimum = ref_mol.current_score;
                    if (restrained_min){
                         Econstraint =  coefficient_restraint * calc_active_rmsd2(rmsd_ref,ref_mol);
                    } else{ 
                         Econstraint = 0.0;
                    }
                    optimum = ref_mol.current_score + ref_mol.internal_energy //trent balius 2008-12-05
                              + Econstraint;

                    // copy best mol to min_mol, and generate min structure
                    copy_crds(min_mol, ref_mol);
                    scale_simplex_vector(new_vec, vertex, trans_step_size,
                                         rot_step_size, tors_step_size);
                    vector_to_dockmol(min_mol, new_vec);
                    copy_crds(mol, min_mol);

                    // free arrays
                    for (x = 0; x < size + 1; x++) {
                        delete[]p[x];
                        p[x] = NULL;
                    }

                    delete[]p;
                    p = NULL;

                    delete[]y;
                    y = NULL;

                    delete[]pr;
                    pr = NULL;

                    delete[]prr;
                    prr = NULL;

                    delete[]pbar;
                    pbar = NULL;

                    delete[]old_vertex;
                    old_vertex = NULL;

                    delta -= optimum;
                    //cout << "Optimum= " << optimum << endl;
                    return optimum;
                }


                if (yprr < y[ilo]) {
                    // The additional extrap succeeded, and replaces the high
                    // point
                    for (i = 0; i < size; i++){
                        p[ihi][i] = prr[i];
                    }
                    y[ihi] = yprr;

                } else {
                    // The additional extrap failed, but still use the
                    // reflected point
                    for (i = 0; i < size; i++){
                        p[ihi][i] = pr[i];
                    }
                    y[ihi] = ypr;
                }

            } else if (ypr >= y[inhi]) {

                // the reflected point is worse than the 2nd highest.  If
                // better than the highest, replace the highest
                replace_flag = false;

                if (ypr < y[ihi]) {
                    for (i = 0; i < size; i++){
                        p[ihi][i] = pr[i];
                    }
                    y[ihi] = ypr;
                    replace_flag = true;
                }
                // contract simplex in 1-D, then eval fxn
                for (i = 0; i < size; i++){
                    vertex[i] = prr[i] =
                        beta * p[ihi][i] + (1.0 - beta) * pbar[i];
                }
                if (simplex_score (score, ref_mol, tmp_mol, vertex, trans_step_size,
                     rot_step_size, tors_step_size)) {
                     if (restrained_min){
                          Econstraint = coefficient_restraint * calc_active_rmsd2(rmsd_ref,tmp_mol);
                     } else {
                          Econstraint = 0.0;
                     }
                     yprr = tmp_mol.current_score + tmp_mol.internal_energy
                              + Econstraint;

                } else {
                    // store best scoring vertex
                    for (x = 0; x < size; x++){
                        vertex[x] = p[ilo][x];
                    }
                    //optimum = ref_mol.current_score;
                    if (restrained_min){
                        Econstraint = coefficient_restraint * calc_active_rmsd2(rmsd_ref,ref_mol);
                    } else{ 
                        Econstraint = 0.0;
                    }
                    optimum = ref_mol.current_score + ref_mol.internal_energy //trent balius 2008-12-05
                              + Econstraint;

                    // copy best mol to min_mol, and generate min structure
                    copy_crds(min_mol, ref_mol);
                    scale_simplex_vector(new_vec, vertex, trans_step_size,
                                         rot_step_size, tors_step_size);
                    vector_to_dockmol(min_mol, new_vec);
                    copy_crds(mol, min_mol);

                    // free arrays
                    for (x = 0; x < size + 1; x++) {
                        delete[]p[x];
                        p[x] = NULL;
                    }

                    delete[]p;
                    p = NULL;

                    delete[]y;
                    y = NULL;

                    delete[]pr;
                    pr = NULL;

                    delete[]prr;
                    prr = NULL;

                    delete[]pbar;
                    pbar = NULL;

                    delete[]old_vertex;
                    old_vertex = NULL;

                    delta -= optimum;
                    //cout << "Optimum= " << optimum << endl;
                    return optimum;
                }

                if (yprr < y[ihi]) {
                    // contraction is an improvement, so accept it

                    for (i = 0; i < size; i++){
                        p[ihi][i] = prr[i];
                    }
                    y[ihi] = yprr;
                    replace_flag = true;
                }

                if (replace_flag == false) {
                    // can't elim high point.  contract about low point
                    for (i = 0; i < size + 1; i++){
                        if (i != ilo) {
                            for (j = 0; j < size; j++) {
                                vertex[j] = p[i][j] = pr[j] =
                                    0.5 * (p[i][j] + p[ilo][j]);
                            }

                            // check if low point is valid
                            if (simplex_score (score, ref_mol, tmp_mol, vertex,
                                 trans_step_size, rot_step_size, tors_step_size)) {
                                if (restrained_min){
                                      Econstraint = coefficient_restraint * calc_active_rmsd2(rmsd_ref,tmp_mol);
                                } else{ 
                                      Econstraint = 0.0;
                                }
                                y[i] = tmp_mol.current_score + tmp_mol.internal_energy
                                        + Econstraint;
                            }
                            else {
                                // store best scoring vertex
                                for (x = 0; x < size; x++)
                                    vertex[x] = p[ilo][x];

                                //optimum = ref_mol.current_score;
                                if (restrained_min){
                                     Econstraint = coefficient_restraint * calc_active_rmsd2(rmsd_ref,ref_mol);
                                }else{ 
                                     Econstraint = 0.0;
                                }
                                
                                optimum = ref_mol.current_score + ref_mol.internal_energy //trent balius 2008-12-05
                                                  + Econstraint;

                                // copy best mol to min_mol, and generate min
                                // structure
                                copy_crds(min_mol, ref_mol);
                                scale_simplex_vector(new_vec, vertex,
                                                     trans_step_size,
                                                     rot_step_size,
                                                     tors_step_size);
                                vector_to_dockmol(min_mol, new_vec);
                                copy_crds(mol, min_mol);

                                // free arrays
                                for (x = 0; x < size + 1; x++) {
                                    delete[]p[x];
                                    p[x] = NULL;
                                }

                                delete[]p;
                                p = NULL;

                                delete[]y;
                                y = NULL;

                                delete[]pr;
                                pr = NULL;

                                delete[]prr;
                                prr = NULL;

                                delete[]pbar;
                                pbar = NULL;

                                delete[]old_vertex;
                                old_vertex = NULL;

                                delta -= optimum;
                                //cout << "Optimum= " << optimum << endl;
                                return optimum;
                            }
                        }
                     
               }
            }
            } else {
                // orig reflection gives a middling point.  Replace high point
                // & move on

                for (i = 0; i < size; i++){
                    p[ihi][i] = pr[i];
                }
                y[ihi] = ypr;
            }

        }

        // ID Best & Worst vertices in current simplex

        if (y[0] > y[1]) {
            ihi = 0;
            inhi = 1;
        } else {
            ihi = 1;
            inhi = 0;
        }

        // loop over simplex points
        for (i = 0; i < size + 1; i++) {

            if (y[i] < y[ilo]){
                ilo = i;
            }
            if (y[i] > y[ihi]) {
                inhi = ihi;
                ihi = i;
            } else if (y[i] > y[inhi]) {
                if (i != ihi)
                    inhi = i;
            }
        }

        
/**
        // print out simplex trajectory
        for(i=0;i<size;i++)
            vertex[i] = prr[i] = p[ilo][i];

        copy_molecule(min_mol, ref_mol);
        scale_simplex_vector(new_vec, vertex, trans_step_size, rot_step_size, tors_step_size);
        vector_to_dockmol(min_mol, new_vec);
        //Write_Mol2(min_mol, cout);

        // DTM - 11-12-08 output simplex scores!
        score.compute_score(min_mol);
        // end print out of trajectory
**/

    } while ((iteration++ < max_iterations)
             && (fabs(y[ihi] - y[ilo]) > score_converge));

    // store best scoring vertex
    for (i = 0; i < size; i++){
        vertex[i] = p[ilo][i];
    }
    optimum = y[ilo];

/***
    // Monte Carlo Post-minimizer
    if(use_mc_postmin) {

        for(i=0;i<size;i++)
            old_vertex[i] = vertex[i];

        temp2 = simplex_score(score, ref_mol, tmp_mol, vertex, trans_step_size, rot_step_size, tors_step_size);
        fail_count = 0;
        step_count = 0;

        float ave_score = 0.0;
        int ave_count = 0;

        do {

            temp1 = simplex_score(score, ref_mol, tmp_mol, vertex, trans_step_size, rot_step_size, tors_step_size);
            diff = temp1-temp2;
            step_count++;

            if((diff < 0)||( exp(-diff) > r_value )) {
                for(i=0;i<size;i++)
                    old_vertex[i] = vertex[i];

                ave_score += temp1;
                ave_count++;

                temp2 = temp1;

            } else
                fail_count++;

            for(i=0;i<size;i++)
                vertex[i] = old_vertex[i] + 0.5*(((float)rand()/(float)RAND_MAX) - 0.5);

        } while((step_count < max_steps)&&(fail_count < fail_threshold));

        for(i=0;i<size;i++)
                vertex[i] = old_vertex[i];

        optimum = temp2;
    }
***/

    // copy best mol to min_mol, and generate min structure
    copy_crds(min_mol, ref_mol);
    scale_simplex_vector(new_vec, vertex, trans_step_size, rot_step_size,
                         tors_step_size);
    vector_to_dockmol(min_mol, new_vec);
    copy_crds(mol, min_mol);

    // free arrays
    for (i = 0; i < size + 1; i++) {
        delete[]p[i];
        p[i] = NULL;
    }

    delete[]p;
    p = NULL;

    delete[]y;
    y = NULL;

    delete[]pr;
    pr = NULL;

    delete[]prr;
    prr = NULL;

    delete[]pbar;
    pbar = NULL;

    delete[]old_vertex;
    old_vertex = NULL;

    delta -= optimum;
    //cout << "Optimum= " << optimum << endl;
    return optimum;

}

/******************************************************/
void
Simplex_Minimizer::scale_simplex_vector(FLOATVec & new_vec, FLOATVec & vertex,
                                        float trans_step_size,
                                        float rot_step_size,
                                        float tors_step_size)
{
    int             i;

    new_vec.resize(vertex.size(), 0);

    for (i = 0; i < 3; i++) {
        new_vec[i] =
            (vertex[i] * trans_step_size) / (float) (current_cycle + 1);
        new_vec[i + 3] =
            (vertex[i + 3] * rot_step_size) / (float) (current_cycle + 1);
    }

    for (i = 6; i < vertex.size(); i++) {
        new_vec[i] =
            (vertex[i] * tors_step_size) / ((float) (current_cycle + 1) *
                                            (float) (torsion_scale_factors
                                                     [i - 6]));
    }

}

/******************************************************/
bool
Simplex_Minimizer::simplex_score(Base_Score & score, DOCKMol & ref_mol,
                                 DOCKMol & tmp_mol, FLOATVec & vertex,
                                 float trans_step_size, float rot_step_size,
                                 float tors_step_size)
{
    FLOATVec        new_vec;
    bool            return_val;

    copy_crds(tmp_mol, ref_mol);
    scale_simplex_vector(new_vec, vertex, trans_step_size, rot_step_size,
                         tors_step_size);
    vector_to_dockmol(tmp_mol, new_vec);

    // compute internal energy as well
    score.compute_ligand_internal_energy(tmp_mol);

    return_val = score.compute_score(tmp_mol);

    //cout << "In Simplex_Minimizer::simplex_score: "
    //     << "score = " << tmp_mol.current_score
    //     << ";int  = " << tmp_mol.internal_energy << endl;

    return return_val;
}



/******************************************************/
void
Simplex_Minimizer::vector_to_dockmol(DOCKMol & mol, FLOATVec & v)
{
    DOCKVector      com, // Centre of Mass
                    dv;  // Translation Vector
    int             i;
    float           rmat[3][3];  // Rotational Matrix
    float           quat[3];
    float           current_angle,
                    new_angle;

    // calc COM of active atoms
    com.x = 0;
    com.y = 0;
    com.z = 0;

    for (i = 0; i < mol.num_atoms; i++) {
        if (mol.atom_active_flags[i]) {
            com.x += mol.x[i];
            com.y += mol.y[i];
            com.z += mol.z[i];
        }
    }

    com.x = com.x / mol.num_active_atoms;
    com.y = com.y / mol.num_active_atoms;
    com.z = com.z / mol.num_active_atoms;

    // build a rotation matrix
    quat[0] = v[3];
    quat[1] = v[4];
    quat[2] = v[5];

    get_matrix_from_quaternion(rmat, quat);

    // build translation vector
    dv.x = v[0];
    dv.y = v[1];
    dv.z = v[2];

    // transform mol
    transform(mol, rmat, dv, com);

    // set new torsion angles
    for (i = 6; i < v.size(); i++) {

        if (bond_vectors[torsions[i - 6].bond_num] == -1) {     // if bond
                                                                // directions
                                                                // don't matter

            current_angle =
                mol.get_torsion(torsions[i - 6].atom1, torsions[i - 6].atom2,
                                torsions[i - 6].atom3, torsions[i - 6].atom4);
            new_angle = (PI / 180.0) * (current_angle + v[i]);
            mol.set_torsion(torsions[i - 6].atom1, torsions[i - 6].atom2,
                            torsions[i - 6].atom3, torsions[i - 6].atom4,
                            new_angle);

        } else {                // if bond directions do matter (during flex
                                // growth)

            if (torsions[i - 6].atom2 == bond_vectors[torsions[i - 6].bond_num]) {
                current_angle =
                    mol.get_torsion(torsions[i - 6].atom1,
                                    torsions[i - 6].atom2,
                                    torsions[i - 6].atom3,
                                    torsions[i - 6].atom4);
                new_angle = (PI / 180.0) * (current_angle + v[i]);
                mol.set_torsion(torsions[i - 6].atom1, torsions[i - 6].atom2,
                                torsions[i - 6].atom3, torsions[i - 6].atom4,
                                new_angle);
            }

            if (torsions[i - 6].atom3 == bond_vectors[torsions[i - 6].bond_num]) {
                current_angle =
                    mol.get_torsion(torsions[i - 6].atom4,
                                    torsions[i - 6].atom3,
                                    torsions[i - 6].atom2,
                                    torsions[i - 6].atom1);
                new_angle = (PI / 180.0) * (current_angle + v[i]);
                mol.set_torsion(torsions[i - 6].atom4, torsions[i - 6].atom3,
                                torsions[i - 6].atom2, torsions[i - 6].atom1,
                                new_angle);
            }

        }

        //cout << "torsion:: " << (PI / 180.0) * current_angle << " " << new_angle << endl; 

    }

}

/******************************************************/
float
Simplex_Minimizer::calc_active_rmsd2(DOCKMol & ref, DOCKMol & conf)
{
 // This function is used to tether the molecule to prevent the previous growth step.
 // this function calculates rmsd2 between the active atoms of the ref structure and the same atoms of conf.

 // The rmsd2 can be thought of as the mean of the squared distances.

 // only heavy atom rmsd2 is reported
 // the rmsd2 of the active atoms in the reference is reported

//    if (! restrained_min) { // if restrained minimum is not used do not compute the rmsd2
//        return 0;
//    }


    int    i;
    float  rmsd2 = 0.0;
    int    atom_num_total = 0;

    
    for (i = 0; i < ref.num_atoms; i++) {
          if (ref.atom_active_flags[i] && ref.amber_at_heavy_flag[i]){
                    rmsd2 +=
                        ((ref.x[i] - conf.x[i]) * (ref.x[i] - conf.x[i]) +
                         (ref.y[i] - conf.y[i]) * (ref.y[i] - conf.y[i]) +
                         (ref.z[i] - conf.z[i]) * (ref.z[i] - conf.z[i]));

                    atom_num_total += 1;
          }
    }

    if (atom_num_total > 0)
        rmsd2 = rmsd2 / (float) atom_num_total;
    else
        rmsd2 = 0.0;

    return rmsd2;
}

