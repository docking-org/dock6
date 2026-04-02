#include <iostream>
#include "master_score.h"

using namespace std;


/******************************************************/
Master_Score::Master_Score(){
    primary_score_found = false;
    secondary_score_found = false;
    use_primary_score = false;
    use_secondary_score = false;
    read_gb_parm = false;
    primary_min = false;
    secondary_min = false;
    read_vdw = false;
    use_chem = false;
    amber = false;
    use_cnt_grid = false;
    use_nrg_grid = false;

    ir_ensemble  = false;
}

/******************************************************/
Master_Score::~Master_Score(){
   //close_all();
   //delete primary_score;
   //delete secondary_score;
}

/******************************************************/
void
Master_Score::input_parameters(Parameter_Reader & parm)
{
    string          tmp;

    primary_score_found = false;
    secondary_score_found = false;
    use_primary_score = false;
    use_secondary_score = false;

    // flags for whether other portions of code should be
    // activated depending on scoring function
    read_gb_parm = false;
    primary_min = false;
    secondary_min = false;
    read_vdw = false;
    use_ph4 = false;//LINGLING
    use_volume = false;//YUCHEN
    amber = false;
    use_cnt_grid = false;
    use_nrg_grid = false;

    cout << "\nMaster Score Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    tmp = parm.query_param("score_molecules", "yes", "yes no");

    if (tmp == "yes")
        use_score = true;
    else
        use_score = false;

    if (use_score) {

        c_cnt.use_score = false;
        c_cnt.use_primary_score = false;
        c_cnt.use_secondary_score = false;
        c_nrg.use_score = false;
        c_nrg.use_primary_score = false;
        c_nrg.use_secondary_score = false;
        c_gist.use_score = false;
        c_gist.use_primary_score = false;
        c_gist.use_secondary_score = false;
        c_mg_nrg.use_score = false;
        c_mg_nrg.use_primary_score = false;
        c_mg_nrg.use_secondary_score = false;
        c_cmg.use_score = false;
        c_cmg.use_primary_score = false;
        c_cmg.use_secondary_score = false;
        c_cont_nrg.use_score = false;
        c_cont_nrg.use_primary_score = false;
        c_cont_nrg.use_secondary_score = false;
        c_fps.use_score = false;
        c_fps.use_primary_score = false;
        c_fps.use_secondary_score = false;
        c_ph4.use_score = false;
        c_ph4.use_primary_score = false;
        c_ph4.use_secondary_score = false;//LINGLING
        //c_volume.use_score = false;
        //c_volume.use_primary_score = false;
        //c_volume.use_secondary_score = false;//YUCHEN
        c_hbond.use_score = false;
        c_hbond.use_primary_score = false;
        c_hbond.use_secondary_score = false;
        c_int.use_score = false;
        c_int.use_primary_score = false;
        c_int.use_secondary_score = false;
        c_desc.use_score = false;
        c_desc.use_primary_score = false;
        c_desc.use_secondary_score = false;
        c_gbsa.use_score = false;
        c_gbsa.use_primary_score = false;
        c_gbsa.use_secondary_score = false;
        c_gbsa_hawkins.use_score = false;
        c_gbsa_hawkins.use_primary_score = false;
        c_gbsa_hawkins.use_secondary_score = false;
        c_sasa.use_score = false;
        c_sasa.use_primary_score = false;
        c_sasa.use_secondary_score = false;
        c_amber.use_score = false;
        c_amber.use_primary_score = false;
        c_amber.use_secondary_score = false;


        if ((!primary_score_found) || (!secondary_score_found)) {
            c_cnt.input_parameters(parm, primary_score_found,
                                   secondary_score_found);

            if (c_cnt.use_primary_score) {
                primary_score = &c_cnt;
                primary_min = true;
                read_vdw = true;
                use_cnt_grid = true;
            }
            if (c_cnt.use_secondary_score) {
                secondary_score = &c_cnt;
                secondary_min = true;
                read_vdw = true;
                use_cnt_grid = true;
            }
        }

        if ((!primary_score_found) || (!secondary_score_found)) {
            c_nrg.input_parameters(parm, primary_score_found,
                                   secondary_score_found);

            if (c_nrg.use_primary_score) {
                primary_score = &c_nrg;
                primary_min = true;
                read_vdw = true;
                use_nrg_grid = true;
            }
            if (c_nrg.use_secondary_score) {
                secondary_score = &c_nrg;
                secondary_min = true;
                read_vdw = true;
                use_nrg_grid = true;
            }
        }

        if ((!primary_score_found) || (!secondary_score_found)) {
            c_gist.input_parameters(parm, primary_score_found,
                                   secondary_score_found);

            if (c_gist.use_primary_score) {
                primary_score = &c_gist;
                primary_min = true;
                read_vdw = true;
                use_nrg_grid = true;
            }
            if (c_gist.use_secondary_score) {
                secondary_score = &c_gist;
                secondary_min = true;
                read_vdw = true;
                use_nrg_grid = true;
            }
        }

        if ((!primary_score_found) || (!secondary_score_found)) {
            c_mg_nrg.input_parameters(parm, primary_score_found,
                                   secondary_score_found);

            if (c_mg_nrg.use_primary_score) {
                primary_score = &c_mg_nrg;
                primary_min = true;
                read_vdw = true;
                use_nrg_grid = true;

                if (c_mg_nrg.ir_ensemble)
                    ir_ensemble = true;
            }
            if (c_mg_nrg.use_secondary_score) {
                secondary_score = &c_mg_nrg;
                secondary_min = true;
                read_vdw = true;
                use_nrg_grid = true;
                if (c_mg_nrg.ir_ensemble) {
                    cout << "IR ensemble does not suport secondary score" << endl;
                    exit(0);
                }
            }
        }

        if ((!primary_score_found) || (!secondary_score_found)) {
            c_cmg.input_parameters(parm, primary_score_found,
                                   secondary_score_found);

            if (c_cmg.use_primary_score) {
                primary_score = &c_cmg;
                primary_min = true;
                read_vdw = true;
            }
            if (c_cmg.use_secondary_score) {
                secondary_score = &c_cmg;
                secondary_min = true;
                read_vdw = true;
            }
        }

        if ((!primary_score_found) || (!secondary_score_found)) {
            c_cont_nrg.input_parameters(parm, primary_score_found,
                                        secondary_score_found);

            if (c_cont_nrg.use_primary_score) {
                primary_score = &c_cont_nrg;
                primary_min = true;
                read_vdw = true;
            }
            if (c_cont_nrg.use_secondary_score) {
                secondary_score = &c_cont_nrg;
                secondary_min = true;
                read_vdw = true;
            }
        }

        if ((!primary_score_found) || (!secondary_score_found)) {
            c_fps.input_parameters(parm, primary_score_found,
                                   secondary_score_found);

            if (c_fps.use_primary_score) {
                primary_score = &c_fps;
                primary_min = true;
                read_vdw = true;
            }
            if (c_fps.use_secondary_score) {
                secondary_score = &c_fps;
                secondary_min = true;
                read_vdw = true;
            }
        }

        if ((!primary_score_found) || (!secondary_score_found)) {
            c_ph4.input_parameters(parm, primary_score_found,
                                        secondary_score_found);
            // LINGLING
            if (c_ph4.use_primary_score) {
                primary_score = &c_ph4;
                primary_min = true;
                read_vdw = true;
                use_ph4 = true;
            }
            if (c_ph4.use_secondary_score) {
                secondary_score = &c_ph4;
                secondary_min = true;
                read_vdw = true;
                use_ph4 = true;
            }
        }

/*
        if ((!primary_score_found) || (!secondary_score_found)) {
            c_volume.input_parameters(parm, primary_score_found,
                                   secondary_score_found);

            if (c_volume.use_primary_score) {
                primary_score = &c_volume;
                primary_min = true;
                read_vdw = true;
                use_volume = true;
            }
            if (c_volume.use_secondary_score) {
                secondary_score = &c_volume;
                secondary_min = true;
                read_vdw = true;
                use_volume = true;
            }
        }
*/

        if ((!primary_score_found) || (!secondary_score_found)) {
            c_hbond.input_parameters(parm, primary_score_found,
                                            secondary_score_found);

            if (c_hbond.use_primary_score) {
                primary_score = &c_hbond;
                primary_min = true;
                read_vdw = true;
                read_gb_parm = true;
            }
            if (c_hbond.use_secondary_score) {
                secondary_score = &c_hbond;
                secondary_min = true;
                read_vdw = true;
                read_gb_parm = true;
            }
        }

        if ((!primary_score_found) || (!secondary_score_found)) {
            c_int.input_parameters(parm, primary_score_found,
                                            secondary_score_found);

            if (c_int.use_primary_score) {
                primary_score = &c_int;
                primary_min = true;
                read_vdw = true;
                read_gb_parm = true;
            }
            if (c_int.use_secondary_score) {
                secondary_score = &c_int;
                secondary_min = true;
                read_vdw = true;
                read_gb_parm = true;
            }
        }


        if ((!primary_score_found) || (!secondary_score_found)) {
            c_desc.input_parameters(parm, primary_score_found,
                                        secondary_score_found);

            if (c_desc.use_primary_score) {
                primary_score = &c_desc;
                primary_min = true;
                read_vdw = true;
                if(c_desc.desc_use_ph4){// input ph4_defn via use_ph4 when ph4 score is used for desc score calc
                   use_ph4 = true;
                }
                if(c_desc.desc_use_volume) {
                    use_volume = true;
                }
            }
            if (c_desc.use_secondary_score) {
                secondary_score = &c_desc;
                secondary_min = true;
                read_vdw = true;
                if(c_desc.desc_use_ph4){// input ph4_defn via use_ph4 when ph4 score is used for desc score calc
                   use_ph4 = true;
                }
                if(c_desc.desc_use_volume) {
                    use_volume = true;
                }
            }
        }

        if ((!primary_score_found) || (!secondary_score_found)) {
            c_gbsa.input_parameters(parm, primary_score_found,
                                    secondary_score_found);

            if (c_gbsa.use_primary_score) {
                primary_score = &c_gbsa;
                primary_min = true;
                read_vdw = true;

            }
            if (c_gbsa.use_secondary_score) {
                secondary_score = &c_gbsa;
                secondary_min = true;
                read_vdw = true;
            }
        }

        if ((!primary_score_found) || (!secondary_score_found)) {
            c_gbsa_hawkins.input_parameters(parm, primary_score_found,
                                            secondary_score_found);

            if (c_gbsa_hawkins.use_primary_score) {
                primary_score = &c_gbsa_hawkins;
                primary_min = true;
                read_vdw = true;
                read_gb_parm = true;
                if (!c_gbsa_hawkins.cont_vdw)
                    use_nrg_grid = true;
            }
            if (c_gbsa_hawkins.use_secondary_score) {
                secondary_score = &c_gbsa_hawkins;
                secondary_min = true;
                read_vdw = true;
                read_gb_parm = true;
                if (!c_gbsa_hawkins.cont_vdw)
                    use_nrg_grid = true;
            }
        }
        if ((!primary_score_found) || (!secondary_score_found)) {
            c_sasa.input_parameters(parm, primary_score_found,
                                            secondary_score_found);

            if (c_sasa.use_primary_score) {
                primary_score = &c_sasa;
                primary_min = true;
                read_vdw = true;
                read_gb_parm = true;
            }
            if (c_sasa.use_secondary_score) {
                secondary_score = &c_sasa;
                secondary_min = true;
                read_vdw = true;
                read_gb_parm = true;
            }
        }

        if ((!primary_score_found) || (!secondary_score_found)) {
            c_amber.input_parameters(parm, primary_score_found,
                                     secondary_score_found);
            if (c_amber.use_primary_score) {
                primary_score = &c_amber;
                amber = true;
            }
            if (c_amber.use_secondary_score) {
                secondary_score = &c_amber;
                amber = true;
            }
        }

        if (!primary_score_found) {
            cout <<
                "\nError:  No scoring function selected for the primary score.  Program will now terminate.\n"
                << endl;
            exit(0);
        } else
            use_primary_score = true;

        if (!secondary_score_found) {
            //cout << "\nNote: No secondary scoring function selected.\n" <<
            //    endl;
        } else
            use_secondary_score = true;

    } else {
        use_primary_score = false;
        use_secondary_score = false;
        // TEB 2010-03-31
        // if no no scoring function then use base_score
        // this is so we still have access to the internal energy.
        c_bas.use_score = false;
        c_bas.use_primary_score = true;
        c_bas.use_secondary_score = false;
        primary_score = &c_bas;
        use_primary_score = false;
        use_secondary_score = false;

    }

}

/******************************************************/
void
Master_Score::initialize_all(AMBER_TYPER & typer, int argc, char **argv)  //internal_min no longer used
// internal_min was not used in this function. this is likely why dock was not using internal energy correctly earlier
{

    if (use_primary_score) {

        if (c_cnt.use_score) {
            c_cnt.initialize(typer);
        }

        if (c_nrg.use_score) {
            c_nrg.initialize(typer);
        }

        if (c_gist.use_score) {
            c_gist.initialize(typer);
        }

        if (c_mg_nrg.use_score) {
            c_mg_nrg.initialize(typer);
        }

        if (c_cmg.use_score) {
            c_cmg.initialize(typer);
        }

        if (c_cont_nrg.use_score) {
            c_cont_nrg.initialize(typer);
        }

        if (c_fps.use_score) {
            c_fps.initialize(typer);
            c_fps.submit_footprint_reference(typer);
        }

        if (c_ph4.use_score) {
            c_ph4.initialize(typer); //LINGLING
            // Ph4 citation
            cout << "To cite Pharmacophore Score use: \n Jiang, L.; Rizzo, R. C. Pharmacophore-Based Similarity Scoring for DOCK, J. Phys. Chem. B., 2015, 119, 1083-1102. \n" << endl;
/*
            c_ph4.submit_ph4_reference(typer); // pharmacophore
*/
        }

/*
        if (c_volume.use_score) {
            c_volume.initialize(typer);
        }    //YUCHEN
*/

        if (c_hbond.use_score == true) {
            c_hbond.initialize(typer);
        }

        if (c_int.use_score == true) {
            c_int.initialize(typer);
        }

        if (c_desc.use_score) {
            c_desc.initialize(typer);
        }

        if (c_gbsa.use_score) {
            c_gbsa.initialize(typer);
        }

        if (c_gbsa_hawkins.use_score) {
            c_gbsa_hawkins.initialize(typer);
        }
 
        if (c_sasa.use_score) {
            c_sasa.initialize(typer);
        }

        if (c_amber.use_score) {
            c_amber.initialize(typer);
        }

    } 

}

/******************************************************/
void
Master_Score::close_all()  //
{
    if (use_primary_score) {

        if (c_mg_nrg.use_score) {
            c_mg_nrg.close();
        }
        if (c_fps.use_score) {
            c_fps.close();
        }
        if (c_ph4.use_score) {
            c_ph4.close();
        }
        if (c_desc.use_score) {
            c_desc.close();
        }
    }
}

/******************************************************/
bool
Master_Score::compute_primary_score(DOCKMol & mol)
{

    if (use_primary_score){
        primary_score->compute_ligand_internal_energy(mol); // trent 2009-02-13
        return primary_score->compute_score(mol);
    } else
        return true;
}

/******************************************************/
bool
Master_Score::compute_secondary_score(DOCKMol & mol)
{

    if (use_secondary_score){
        secondary_score->compute_ligand_internal_energy(mol); // trent 2009-02-13
        return secondary_score->compute_score(mol);
    } else
        return true;
}

