#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <sstream>
#include <math.h>

#include "amber_typer.h"
#include "base_score.h"
#include "dockmol.h"
#include "grid.h" 
#include "score_descriptor.h"
#include "utils.h"

using namespace std;


// +++++++++++++++++++++++++++++++++++++++++
// Constructor and destructor
Descriptor_Energy_Score::Descriptor_Energy_Score()
{
}
Descriptor_Energy_Score::~Descriptor_Energy_Score()
{
}



// +++++++++++++++++++++++++++++++++++++++++
// Clear contents of C++ strings and arrays
void
Descriptor_Energy_Score::close()
{

    // close each of the classes here?

}



// +++++++++++++++++++++++++++++++++++++++++
// Get parameters from input file
void
Descriptor_Energy_Score::input_parameters(Parameter_Reader & parm, bool & primary_score, bool & secondary_score)
{
    string tmp;
    use_primary_score = false;
    use_secondary_score = false;
    desc_use_nrg = false;
    desc_use_mg_nrg = false;
    desc_use_cont_nrg = false;
    desc_use_fps = false;
    desc_use_ph4 = false;
    desc_use_tan = false;
    desc_use_hun = false;
    desc_use_volume = false;
    desc_use_gist = false;
    desc_use_cmg = false;

    desc_c_nrg.use_score = false;
    desc_c_mg_nrg.use_score = false;
    desc_c_cont_nrg.use_score = false;
    desc_c_fps.use_score = false;
    desc_c_ph4.use_score = false;


    cout << "\nDescriptor Score Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;


    if (!primary_score) {
        tmp = parm.query_param("descriptor_score_primary", "no", "yes no");
        use_primary_score = (tmp == "yes");
        primary_score = use_primary_score;
    }

    if (!secondary_score) {
     //tmp = parm.query_param("descriptor_score_secondary", "no", "yes no");
       tmp = "no";
       use_secondary_score = (tmp == "yes");
       secondary_score = use_secondary_score;
    }

    use_score = (use_primary_score || use_secondary_score);

    if (use_score) {

        // Ask the user which scoring functions they would like to use
        //
        // Note: There is a very specific hierarchy of what can and cannot be used in combination.
        // For example, the user should either turn on continuous_energy AND footprint_similarity,
        // OR they can turn on grid_score, OR they can turn on multigrid_score. (I.e. it does not 
        // make sense to turn on grid_score AND continuous_energy, so they should be mutually 
        // exclusive)
        //
        // -> Functions including tanimoto, pharmacophore_score, volume_overlap, mol_properties,
        //    and SASA are not mutually exclusive and can be combined with anything.

        desc_use_nrg = parm.query_param("descriptor_use_grid_score", "yes", "yes no") == "yes";

        if (!desc_use_nrg){
            desc_use_mg_nrg = parm.query_param("descriptor_use_multigrid_score", "yes", "yes no") == "yes";
        }

        if (!desc_use_nrg && !desc_use_mg_nrg){
            desc_use_cont_nrg = parm.query_param("descriptor_use_continuous_score", "yes", "yes no") == "yes";
            desc_use_fps = parm.query_param("descriptor_use_footprint_similarity", "yes", "yes no") == "yes";
        }

        desc_use_ph4 = parm.query_param("descriptor_use_pharmacophore_score", "no", "yes no") == "yes";
        desc_use_tan = parm.query_param("descriptor_use_tanimoto", "no", "yes no") == "yes";
	desc_use_hun = parm.query_param("descriptor_use_hungarian", "no", "yes no") == "yes";
        desc_use_volume = parm.query_param("descriptor_use_volume_overlap", "no", "yes no") == "yes";
        //desc_use_molprop = parm.query_param("descriptor_use_molecular_properties", "no", "yes no") == "yes";
//        desc_use_gist = parm.query_param("descriptor_use_gist", "no", "yes no") == "yes";
//        desc_use_cmg = parm.query_param("descriptor_use_dock3.5", "no", "yes no") == "yes";


        // If using grid_score, read in the parameters for grid_score. These lines are esentially
        // duplicated from the Energy_Score::input_parameters function, except the input names in
        // the param query are prefixed with 'descriptor', and the options are stored in an object
        // of the Energy_Score class.
        if (desc_use_nrg){

            cout <<"\n--- Descriptor Score Parameters: Grid Score ---" <<endl;
          
            // This must be set to true for initialize and compute score functions
            desc_c_nrg.use_score = true;

            desc_c_nrg.rep_radius_scale = atof(parm.query_param("descriptor_grid_score_rep_rad_scale", "1").c_str());
            if (desc_c_nrg.rep_radius_scale <= 0.0) {
                cout << "ERROR: Parameter must be a float greater than zero. Program will terminate."
                     << endl;
                exit(0);
            }

            desc_c_nrg.vdw_scale = atof(parm.query_param("descriptor_grid_score_vdw_scale", "1").c_str());
            if (desc_c_nrg.vdw_scale <= 0.0) {
                bool off = parm.query_param("descriptor_grid_score_turn_off_vdw", "yes", "yes no") == "yes";
                if (!off) {
                    cout << "ERROR:  Parameter must be a float greater than zero. Program will terminate."
                         << endl;
                    exit(0);
                }
            }

            desc_c_nrg.es_scale = atof(parm.query_param("descriptor_grid_score_es_scale", "1").c_str());
            if (desc_c_nrg.es_scale <= 0.0) {
                bool off = parm.query_param("descriptor_grid_score_turn_off_es", "yes", "yes no") == "yes";
                if (!off) {
                    cout << "ERROR:  Parameter must be a float greater than zero. Program will terminate."
                         << endl;
                    exit(0);
                }
            }

            desc_c_nrg.grid_file_name = parm.query_param("descriptor_grid_score_grid_prefix", "grid" );
        }



        // If using multigrid_score, read in the parameters for multigrid_score. These lines are
        // esentially duplicated from the Multigrid_Energy_Score::input_parameters_main_function, except
        // the input names in the param query are prefixed with 'descriptor', and the options are
        // stored in an object of the Multigrid_Energy_Score class.


        if (desc_use_mg_nrg){

            cout <<"\n--- Descriptor Score Parameters: Multigrid Score ---" <<endl;

            // This must be set to true for initialize and compute score functions
            desc_c_mg_nrg.use_score = true;
            // Call the main input paramter function in multigrid
            desc_c_mg_nrg.input_parameters_main(parm,"descriptor_multigrid_score");
        }



        // If using continuous_score, read in the parameters for continuous_score. These lines are
        // esentially duplicated from the Continuous_Energy_Score::input_parameters function, except
        // the input names in the param query are prefixed with 'descriptor', and the options are
        // stored in an object of the Continuous_Energy_Score class.
        if (desc_use_cont_nrg){

            cout <<"\n--- Descriptor Score Parameters: Continuous Energy Score ---" <<endl;

            // This must be set to true for initialize and compute score functions
            desc_c_cont_nrg.use_score = true;

            desc_c_cont_nrg.receptor_filename = parm.query_param("descriptor_cont_score_rec_filename", "receptor.mol2");

            desc_c_cont_nrg.att_exp = atoi(parm.query_param("descriptor_cont_score_att_exp", "6").c_str());
            if (desc_c_cont_nrg.att_exp <= 0) {
                cout << "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                     << endl;
                exit(0);
            }

            desc_c_cont_nrg.rep_exp = atoi(parm.query_param("descriptor_cont_score_rep_exp", "12").c_str());
            if (desc_c_cont_nrg.rep_exp <= 0) {
                cout << "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                     << endl;
                exit(0);
            }

            desc_c_cont_nrg.rep_radius_scale = atof(parm.query_param("descriptor_cont_score_rep_rad_scale", "1").c_str());
            if (desc_c_cont_nrg.rep_radius_scale <= 0.0) {
                cout << "ERROR:  Parameter must be a float greater than zero. Program will terminate."
                     << endl;
                exit(0);
            }

            // Ask whether we should use a distance dependent dieleectric
            desc_c_cont_nrg.use_ddd = parm.query_param("descriptor_cont_score_use_dist_dep_dielectric", "yes", "yes no") == "yes";

            // Use default dielectric = 4.0 for distance dependent dielectric
            // else use d = 1.0 assuming gas phase
            char default_diel[4];
            if (desc_c_cont_nrg.use_ddd) strcpy(default_diel, "4.0");
            else strcpy(default_diel,"1.0");

            desc_c_cont_nrg.diel_screen = atof(parm.query_param("descriptor_cont_score_dielectric", default_diel).c_str());
            if (desc_c_cont_nrg.diel_screen <= 0.0) {
                cout << "ERROR: Parameter must be a float greater than zero.  Program will terminate."
                     << endl;
                exit(0);
            }

            desc_c_cont_nrg.vdw_scale = atoi(parm.query_param("descriptor_cont_score_vdw_scale", "1").c_str());
            if (desc_c_cont_nrg.vdw_scale <= 0) {
                bool off;
                off = (parm.query_param("descriptor_cont_score_turn_off_vdw", "yes", "yes no") == "yes") ? true : false;
                if (!off) {
                    cout << "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                         << endl;
                    exit(0);
                }
            }

            desc_c_cont_nrg.es_scale = atoi(parm.query_param("descriptor_cont_score_es_scale", "1").c_str());
            if (desc_c_cont_nrg.es_scale <= 0) {
                bool off;
                off = (parm.query_param("descriptor_cont_score_turn_off_es", "yes", "yes no") == "yes") ? true : false;
                if (!off) {
                    cout << "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                         << endl;
                    exit(0);
                }
            }
        }



        // If using fps_score, read in the parameters for fps_score. These lines are esentially 
        // duplicated from the Footprint_Similarity_Score::input_parameters_main function, except the
        // input names in the param query are prefixed with 'descriptor', and the options are 
        // stored in an object of the Footprint_Similarity_Score class.
        if (desc_use_fps){

            cout <<"\n--- Descriptor Score Parameters: Footprint Similarity Score ---" <<endl;

            // This must be set to true for initialize and compute score functions
            desc_c_fps.use_score = true;


            // initialize some parameters
            string tmp;

            desc_c_fps.fps_foot_specify_a_range = false;
            desc_c_fps.fps_foot_specify_a_threshold = false;
            desc_c_fps.fps_normalize_foot = false;
            desc_c_fps.fps_use_remainder = false;     // remainder is an element of the footprint vector of
                                                      // the remaining residues
          
            desc_c_fps.input_parameters_main(parm,"descriptor_fps_score");
        }

         // Read in parameters to compute FMS score. Note: copied from score_ph4.cpp by LINGLING
        if (desc_use_ph4){

            cout <<"\n--- Descriptor Score Parameters: Pharmacophore ---" <<endl;

            // This must be set to true for initialize and compute score functions
            desc_c_ph4.use_score = true;

            desc_c_ph4.bool_ph4_ref_mol2 = false;
            desc_c_ph4.bool_ph4_ref_txt = false;
            // input_parameters_main is used here and in ph4_score. 
            desc_c_ph4.input_parameters_main(parm,"descriptor_fms_score");
        }


        // Read in parameters to compute the fingerprint / Tanimoto. Note: this is not really a
        // scoring function, it fits better in the category of molecular property.
        if (desc_use_tan){

            cout <<"\n--- Descriptor Score Parameters: Fingerprint / Tanimoto ---" <<endl;
            desc_fing_ref_filename = parm.query_param("descriptor_fingerprint_ref_filename", "fing_ref.mol2");
        }

        // Read in parameters to compute the hungarian distance according to Allen and Rizzo, JCIM, 2014
        // The functional form is: Score = C1 ( (#refatoms - #unmatched) / #refatoms) + C2 (RMSD matched)
        if (desc_use_hun){

            cout <<"\n--- Descriptor Score Parameters: Hungarian Matching Similarity ---" <<endl;
            desc_hun_ref_filename = parm.query_param("descriptor_hms_score_ref_filename", "hun_ref.mol2");
            desc_hun_matching_coeff = atof(parm.query_param("descriptor_hms_score_matching_coeff", "-5").c_str());
            desc_hun_rmsd_coeff = atof(parm.query_param("descriptor_hms_score_rmsd_coeff", "1").c_str());
        }

        if (desc_use_volume){
            cout <<"\n--- Descriptor Score Parameters: Volume Overlap Score ---" <<endl;

            // This must be set to true for initialize and compute score functions
            desc_c_vol.use_score = true;
            desc_c_vol.input_parameters_main(parm,"descriptor_volume_score");
        }

        if (desc_use_gist){
            cout <<"\n--- Descriptor Score Parameters: GIST Score ---" <<endl;

            // This must be set to true for initialize and compute score functions
            desc_c_gist.use_score = true;
            desc_c_gist.input_parameters_main(parm,"descriptor_gist_score_");
        }
/*
        if (desc_use_cmg){
            cout <<"\n--- Descriptor Score Parameters: DOCK3.5 Score ---" <<endl;

            // This must be set to true for initialize and compute score functions
            desc_c_cmg.use_score = true;
            //desc_c_cmg.input_parameters_main(parm,"descriptor_dock3.5_score_");
            desc_c_cmg.input_parameters_main(parm,"descriptor_dock3.5_");
        }
*/
        // These parameters are used to toggle the weights of each function
        cout <<"\n--- Descriptor Score Weights ---" <<endl;

        if (desc_use_nrg){
            desc_weight_nrg = atoi(parm.query_param("descriptor_weight_grid_score", "1").c_str());
        }

        if (desc_use_mg_nrg){
            desc_weight_mg_nrg = atoi(parm.query_param("descriptor_weight_multigrid_score", "1").c_str());
        }

        if (desc_use_cont_nrg){
            desc_weight_cont_nrg = atoi(parm.query_param("descriptor_weight_cont_score", "1").c_str());
        }

        if (desc_use_fps){
            desc_weight_fps = atoi(parm.query_param("descriptor_weight_fps_score", "1").c_str());
        }
        
        if (desc_use_ph4){
            desc_weight_ph4 = atoi(parm.query_param("descriptor_weight_pharmacophore_score", "1").c_str());
        }

        if (desc_use_tan){
            desc_weight_tan = atoi(parm.query_param("descriptor_weight_fingerprint_tanimoto", "-1").c_str());
        }

        if (desc_use_hun){
            desc_weight_hun = atoi(parm.query_param("descriptor_weight_hms_score", "1").c_str());
        }
        if (desc_use_volume){
            desc_weight_volume = atoi(parm.query_param("descriptor_weight_volume_overlap_score", "-1").c_str());
        }
        if (desc_use_gist){
            desc_weight_gist = atoi(parm.query_param("descriptor_weight_gist_score", "-1").c_str());
        }
/*
        if (desc_use_cmg){
            desc_weight_cmg = atoi(parm.query_param("descriptor_weight_dock3.5_score", "1").c_str());
        }
*/

/*
        rot_bonds_scale  = atoi(parm.query_param("descriptor_score_rot_bonds_scale", "0").c_str());
        heavyatoms_scale = atoi(parm.query_param("descriptor_score_Heavy_Atoms_scale", "0").c_str());
        hb_don_scale = atoi(parm.query_param("descriptor_score_HB_Donors_scale", "0").c_str());
        hb_acc_scale = atoi(parm.query_param("descriptor_score_HB_Acceptors_scale", "0").c_str());
        molecular_wt_scale = atoi(parm.query_param("descriptor_score_Molecular_Wt_scale", "0").c_str());
        formal_charge_scale = atoi(parm.query_param("descriptor_score_Formal_Charge_scale", "0").c_str());
        XlogP_scale = atoi(parm.query_param("descriptor_score_XLogP_scale", "0").c_str());
*/


        // Count the number of functions that are being used
        int num=0;
        if (desc_use_nrg){ num++; }
        if (desc_use_mg_nrg){ num++; }
        if (desc_use_cont_nrg){ num++; }
        if (desc_use_fps){ num++; }
        if (desc_use_ph4){ num++; }
        if (desc_use_tan){ num++; }
        if (desc_use_hun){ num++; }
        if (desc_use_volume){ num++; }
        if (desc_use_gist){ num++; }
//        if (desc_use_cmg){ num++; }
        int orig_num = num;

        // Print the functional form for the user
        ostringstream function;
        function <<"Descriptor_Score = ";
        if (desc_use_nrg){
            function <<"(" <<desc_weight_nrg <<" * Grid_Score)";
            num--; 
            if (num>0 && num!=orig_num){function <<" + "; }
        }

        if (desc_use_mg_nrg){
            function <<"(" <<desc_weight_mg_nrg <<" * Multigrid_Score)";
            num--; 
            if (num>0 && num!=orig_num){function <<" + "; }
        }

        if (desc_use_cont_nrg){
            function <<"(" <<desc_weight_cont_nrg <<" * Continuous_Score)";
            num--;
            if (num>0 && num!=orig_num){function <<" + "; }
        }

        if (desc_use_fps){
            function <<"(" <<desc_weight_fps <<" * Footprint_Similarity_Score)";
            num--;
            if (num>0 && num!=orig_num){function <<" + "; }
        }

        if (desc_use_ph4){
            function <<"(" <<desc_weight_ph4 <<" * Pharmacophore_Similarity_Score)";
            num--;
            if (num>0 && num!=orig_num){function <<" + "; }
        }

        if (desc_use_tan){
            function <<"(" <<desc_weight_tan <<" * Tanimoto)"; 
            num--;
            if (num>0 && num!=orig_num){function <<" + "; }
        }

        if (desc_use_hun){ 
            function <<"(" <<desc_weight_hun <<" * Hungarian_Matching_Similarity_Score)"; 
            num--; 
            if (num>0 && num!=orig_num){function <<" + "; }
        }

        if (desc_use_volume){ function <<"(" <<desc_weight_volume <<" * Volume_Score)"; num--; }

        if (desc_use_gist){ function <<"(" <<desc_weight_gist <<" * GIST_Score)"; num--; }
//        if (desc_use_cmg){ function <<"(" <<desc_weight_cmg <<" * DOCK3.5_Score)"; num--; }

        //if (desc_use_molprop){ function <<8 * molecular_property"; }

        cout <<endl;
        cout <<"Using these settings, the Descriptor Score function takes the form:" <<endl;
        cout <<function.str();
        cout <<endl;
  
        function.clear();


    } // end if (use_score)

    return;

} // end Descriptor_Energy_Score::input_parameters()



// +++++++++++++++++++++++++++++++++++++++++
// Intialize some variables for the Amber typer
// (this is done when you use continuous / footprint / multigrid score,
// so we probably don't have to do it here)
void
Descriptor_Energy_Score::initialize(AMBER_TYPER & typer)
{

    if (desc_use_nrg){
        cout << "Initializing typer for Grid Score within Descriptor Score..."
             << endl;
        desc_c_nrg.initialize(typer);
    }

    if (desc_use_mg_nrg){
        cout << "Initializing typer for Multigrid Score within Descriptor Score..."
             << endl;
        desc_c_mg_nrg.initialize(typer);
    }

    if (desc_use_cont_nrg){
        cout << "Initializing typer for Continuous Score within Descriptor Score..."
             << endl;
        desc_c_cont_nrg.initialize(typer);
    }

    if (desc_use_fps){
        cout << "Initializing typer for Footprint Similarity Score within Descriptor Score..."
             << endl;
        desc_c_fps.initialize(typer);
        desc_c_fps.submit_footprint_reference(typer);
    }

    if (desc_use_ph4){
        cout << "Initializing typer for Pharmacophore Similarity Score within Descriptor Score..."
             << endl;
        desc_c_ph4.initialize(typer);
    }

    if (desc_use_tan){
        cout << "Initializing Fingerprint computation within Descriptor Score..."
             << endl;

        ifstream fin;
        fin.open(desc_fing_ref_filename.c_str());
        if (fin.fail()){
            cout <<"Error: Could not open " <<desc_fing_ref_filename <<endl;
            exit (0);
        } else {
            Read_Mol2(fing_ref_mol, fin, false, false, false);
            typer.prepare_molecule(fing_ref_mol, true, false, false, false);
        }
    }

    if (desc_use_hun){
        cout << "Initializing Hungarian computation within Descriptor Score..."
             << endl;

        ifstream fin;
        fin.open(desc_hun_ref_filename.c_str());
        if (fin.fail()){
            cout <<"Error: Could not open " <<desc_hun_ref_filename <<endl;
            exit (0);
        } else {
            Read_Mol2(hun_ref_mol, fin, false, false, false);
        }

    }

    if (desc_use_volume){
        cout << "Initializing typer for Volume Overlap Score within Descriptor Score..."
             << endl;
        desc_c_vol.initialize(typer);
    }

    if (desc_use_gist){
        cout << "Initializing typer for GIST Score within Descriptor Score..."
             << endl;
        desc_c_gist.initialize(typer);
    }
/*
    if (desc_use_cmg){
        cout << "Initializing typer for DOCK3.5 Score within Descriptor Score..."
             << endl;
        desc_c_cmg.initialize(typer);
    }
*/
    return;

} // end Descriptor_Energy_Score::initialize()


// +++++++++++++++++++++++++++++++++++++++++
// Compute the score and add it to mol.current_score
bool
Descriptor_Energy_Score::compute_score(DOCKMol & mol)
{
    // Function here to compute components of score
    temp_nrg_score = 0;
    temp_mg_nrg_score = 0;
    temp_cont_nrg_score = 0;
    temp_fps_score = 0;
    temp_ph4_score = 0;
    temp_tan_score = 0;
    temp_hun_score = 0;
    temp_vol_score = 0;
    temp_gist_score = 0;
    temp_cmg_score = 0;
    temp_desc_score = 0;

    // check if all scoring fuctions are sucessful, if not then return false 
    // bool nrg_bool_suc, mg_nrg_bool_suc, fps_bool_suc, ph4_bool_suc, tan_bool_suc, hun_bool_suc, vol_bool_suc, gist_bool_suc ;

    if (desc_use_nrg){
        // Returns valid_orient value from grid or multi-grid through 
        // descriptor score wapper function.
        if (!desc_c_nrg.compute_score(mol)) return false;
        temp_nrg_score = mol.current_score;
        mol.score_nrg = mol.current_score; 
        temp_desc_score += temp_nrg_score * desc_weight_nrg;
    }

    if (desc_use_mg_nrg){
        if (!desc_c_mg_nrg.compute_score(mol)) return false;
        temp_mg_nrg_score = mol.current_score;
        mol.score_mg_nrg = mol.current_score; 
        temp_desc_score += temp_mg_nrg_score * desc_weight_mg_nrg;
    }

    if (desc_use_cont_nrg){
        if (!desc_c_cont_nrg.compute_score(mol)) return false;
        temp_cont_nrg_score = mol.current_score;
        mol.score_cont_nrg = mol.current_score; 
        temp_desc_score += temp_cont_nrg_score * desc_weight_cont_nrg;
    }

    if (desc_use_fps){
        if (!desc_c_fps.compute_score(mol)) return false;
        temp_fps_score = mol.current_score;
        mol.score_fps = mol.current_score; 
        temp_desc_score += temp_fps_score * desc_weight_fps;
    }

    if (desc_use_ph4){
        if (!desc_c_ph4.compute_score(mol)) return false;
        temp_ph4_score = mol.current_score;
        mol.score_ph4 = mol.current_score;
        temp_desc_score += temp_ph4_score * desc_weight_ph4;
    }

    if (desc_use_tan){
        temp_tan_score = desc_c_fing.compute_tanimoto(mol, fing_ref_mol);
        mol.score_tan = temp_tan_score; 
        temp_desc_score += temp_tan_score * desc_weight_tan;
    }

    if (desc_use_hun){
        hun_ref_heavy_atoms=0;
        desc_hun_result = desc_c_hun.calc_Hungarian_RMSD_dissimilar(mol, hun_ref_mol);
        for (int i=0; i<hun_ref_mol.num_atoms; i++){
            if (hun_ref_mol.atom_types[i].compare("H") != 0 && hun_ref_mol.atom_types[i].compare("Du") != 0){
                hun_ref_heavy_atoms++;
            }
        }
        temp_hun_score = desc_hun_matching_coeff * (hun_ref_heavy_atoms - desc_hun_result.second) / hun_ref_heavy_atoms +
                         desc_hun_rmsd_coeff * (desc_hun_result.first);
        mol.score_hun = temp_hun_score; 
        temp_desc_score += temp_hun_score * desc_weight_hun;
    }

    if (desc_use_volume){
        if (!desc_c_vol.compute_score(mol)) return false;
        temp_vol_score = mol.current_score;
        mol.score_vol = mol.current_score; 
        temp_desc_score += temp_vol_score * desc_weight_volume;
    }

    if (desc_use_gist){
        //bool gist_bool_suc = desc_c_gist.compute_score(mol);
        //if (!gist_bool_suc){ return false;} 
        // if gist score was unsuccessful return false. 
        if (!desc_c_gist.compute_score(mol)) return false;
        temp_gist_score = mol.current_score;
        //cout << " temp_gist_score =" << temp_gist_score <<endl;
        temp_desc_score += temp_gist_score * desc_weight_gist;
    }
/*
    if (desc_use_cmg){
        if (!desc_c_cmg.compute_score(mol)) return false;
        temp_cmg_score = mol.current_score;
        temp_desc_score += temp_cmg_score * desc_weight_cmg;
    }
*/

    mol.current_score = temp_desc_score;
    mol.current_data = output_score_summary(mol);

    return true;

} // end Descriptor_Energy_Score::compute_score()


// +++++++++++++++++++++++++++++++++++++++++
// Return a bunch of useful info for mol.current_data
string
Descriptor_Energy_Score::output_score_summary(DOCKMol & mol)
{
    ostringstream text;

    if (use_score) {

        //changed output score component names to manage length and be consistent.  Yuchen 10/24/2016
        text << DELIMITER << setw(STRING_WIDTH) << "Descriptor_Score:"
             << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl;

 
        if (desc_use_nrg){
            text << DELIMITER << setw(STRING_WIDTH) << "Grid_Score:"
                 << setw(FLOAT_WIDTH) << fixed << temp_nrg_score << endl
                 << DELIMITER << setw(STRING_WIDTH) << "desc_Grid_vdw_energy:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_nrg.vdw_component << endl
                 << DELIMITER << setw(STRING_WIDTH) << "desc_Grid_es_energy:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_nrg.es_component  << endl;
        }


        if (desc_use_mg_nrg){
            text << DELIMITER << setw(STRING_WIDTH) << "MultiGrid_Score:"
                 << setw(FLOAT_WIDTH) << fixed << temp_mg_nrg_score << endl
                 << DELIMITER << setw(STRING_WIDTH) << "desc_MGS_vdw_energy:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_mg_nrg.vdw_sum << endl
                 << DELIMITER << setw(STRING_WIDTH) << "desc_MGS_es_energy:"  
                 << setw(FLOAT_WIDTH) << fixed << desc_c_mg_nrg.es_sum << endl
                 << DELIMITER << setw(STRING_WIDTH) << "desc_MGS_vdw+es_energy:"  
                 << setw(FLOAT_WIDTH) << fixed << desc_c_mg_nrg.vdw_sum+desc_c_mg_nrg.es_sum  << endl;

            if ( desc_c_mg_nrg.fp_mol || desc_c_mg_nrg.fp_txt ) {
               if (desc_c_mg_nrg.use_cor){
                   text << DELIMITER << setw(STRING_WIDTH) << "desc_MGS_vdw_fps:"
                        << setw(FLOAT_WIDTH) << fixed << desc_c_mg_nrg.vdw_cor << endl;
                   text << DELIMITER << setw(STRING_WIDTH) << "desc_MGS_es_fps:"
                        << setw(FLOAT_WIDTH) << fixed << desc_c_mg_nrg.es_cor << endl;
                   text << DELIMITER << setw(STRING_WIDTH) << "desc_MGS_vdw+es_fps:"
                        << setw(FLOAT_WIDTH) << fixed << desc_c_mg_nrg.es_cor + desc_c_mg_nrg.vdw_cor << endl;
               }
               else if (desc_c_mg_nrg.use_euc){
                   text << DELIMITER << setw(STRING_WIDTH) << "desc_MGS_vdw_fps:"
                        << setw(FLOAT_WIDTH) << fixed << desc_c_mg_nrg.vdw_euc << endl;
                   text << DELIMITER << setw(STRING_WIDTH) << "desc_MGS_es_fps:"
                        << setw(FLOAT_WIDTH) << fixed << desc_c_mg_nrg.es_euc << endl;
                   text << DELIMITER << setw(STRING_WIDTH) << "desc_MGS_vdw+es_fps:"
                        << setw(FLOAT_WIDTH) << fixed << desc_c_mg_nrg.vdw_euc + desc_c_mg_nrg.es_euc << endl;
               }
               else if (desc_c_mg_nrg.use_norm){
                   text << DELIMITER << setw(STRING_WIDTH) << "desc_MGS_vdw_fps:"
                        << setw(FLOAT_WIDTH) << fixed << desc_c_mg_nrg.vdw_norm << endl;
                   text << DELIMITER << setw(STRING_WIDTH) << "desc_MGS_es_fps:"
                        << setw(FLOAT_WIDTH) << fixed << desc_c_mg_nrg.es_norm << endl;
                   text << DELIMITER << setw(STRING_WIDTH) << "desc_MGS_vdw+es_fps:"
                        << setw(FLOAT_WIDTH) << fixed << desc_c_mg_nrg.es_norm + desc_c_mg_nrg.vdw_norm << endl;
               }
            }
            if ( desc_c_mg_nrg.fp_mol || desc_c_mg_nrg.fp_txt ){
                // Fix this later WJA
                // printarrays(desc_c_mg_nrg.vdw_pose_array, desc_c_mg_nrg.es_pose_array,
                //             desc_c_mg_nrg.numgrids, desc_c_mg_nrg.text);
            }
        }

        if (desc_use_cont_nrg){
        text << DELIMITER << setw(STRING_WIDTH) << "Continuous_Score:"
             << setw(FLOAT_WIDTH) << fixed << temp_cont_nrg_score << endl
             << DELIMITER << setw(STRING_WIDTH) << "desc_Continuous_vdw_energy:"
             << setw(FLOAT_WIDTH) << fixed << desc_c_cont_nrg.vdw_component << endl
             << DELIMITER << setw(STRING_WIDTH) << "desc_Continuous_es_energy:"
             << setw(FLOAT_WIDTH) << fixed << desc_c_cont_nrg.es_component << endl;

        }

        if (desc_use_fps){
            text << DELIMITER << setw(STRING_WIDTH) << "Footprint_Similarity_Score:"
                 << setw(FLOAT_WIDTH) << fixed << temp_fps_score << endl
                 << DELIMITER << setw(STRING_WIDTH) << "desc_FPS_vdw_energy:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_fps.vdw_component << endl
                 << DELIMITER << setw(STRING_WIDTH) << "desc_FPS_es_energy:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_fps.es_component << endl
                 << DELIMITER << setw(STRING_WIDTH) << "desc_FPS_vdw+es_energy:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_fps.vdw_component+desc_c_fps.es_component << endl
                 << DELIMITER << setw(STRING_WIDTH) << "desc_FPS_num_hbond:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_fps.hbond  << endl;
    
            text << DELIMITER << setw(STRING_WIDTH) << "desc_FPS_vdw_fps:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_fps.vdw_foot_dist << endl
                 << DELIMITER << setw(STRING_WIDTH) << "desc_FPS_es_fps:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_fps.es_foot_dist << endl
                 << DELIMITER << setw(STRING_WIDTH) << "desc_FPS_hb_fps:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_fps.hbond_foot_dist << endl
                 << DELIMITER << setw(STRING_WIDTH) << "desc_FPS_vdw_fp_numres:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_fps.fps_vdw_num_resid << endl
                 << DELIMITER << setw(STRING_WIDTH) << "desc_FPS_es_fp_numres:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_fps.fps_es_num_resid << endl
                 << DELIMITER << setw(STRING_WIDTH) << "desc_FPS_hb_fp_numres:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_fps.fps_hb_num_resid  << endl;
             if (desc_c_fps.fps_fp_info != "")
                 text << desc_c_fps.fps_fp_info << endl;
        }

        if (desc_use_ph4){
            text << DELIMITER << setw(STRING_WIDTH) << "Pharmacophore_Score:"
                 << setw(FLOAT_WIDTH) << fixed << temp_ph4_score << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_FMS_num_match_tot:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_ph4.match_num[0] << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_FMS_max_match_ref:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_ph4.match_num[1] << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_FMS_max_match_mol:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_ph4.match_num[2] << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_FMS_match_rate:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_ph4.match_term[0] << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_FMS_match_resid:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_ph4.match_term[1] << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_FMS_num_hydrophobic_matched:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_ph4.match_comp[0] << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_FMS_num_donor_matched:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_ph4.match_comp[1] << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_FMS_num_acceptor_matched:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_ph4.match_comp[2] << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_FMS_num_aromatic_matched:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_ph4.match_comp[3] << endl;
            //text << DELIMITER << setw(STRING_WIDTH) << "desc_FMS_num_aroAcc_matched:"
            //     << setw(FLOAT_WIDTH) << fixed << desc_c_ph4.match_comp[4] << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_FMS_num_positive_matched:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_ph4.match_comp[5] << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_FMS_num_negative_matched:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_ph4.match_comp[6] << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_FMS_num_ring_matched:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_ph4.match_comp[7] << endl;
            //Numer of matched features changed from a loop to hardcoded for formatting purpose Yuchen 10/24/2016
            //Each matched term has a component score that is not printed here in this version Yuchen 10/24/2016
        }

        if (desc_use_tan){
            text << DELIMITER << setw(STRING_WIDTH) << "Tanimoto_Score:"
                 << setw(FLOAT_WIDTH) << fixed << temp_tan_score << endl;
        }

        if (desc_use_hun){
            text << DELIMITER << setw(STRING_WIDTH) << "Hungarian_Matching_Similarity_Score:"
                 << setw(FLOAT_WIDTH) << fixed << temp_hun_score << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_HMS_num_unmtchd_hvy_atms:"
                 << setw(FLOAT_WIDTH) << fixed << desc_hun_result.second << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_HMS_num_ref_hvy_atms:"
                 << setw(FLOAT_WIDTH) << fixed << hun_ref_heavy_atoms << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_HMS_rmsd_mtchd_hvy_atms:"
                 << setw(FLOAT_WIDTH) << fixed << desc_hun_result.first << endl;
        }

        if (desc_use_volume){
            text << DELIMITER << setw(STRING_WIDTH) << "Property_Volume_Score:"
                 << setw(FLOAT_WIDTH) << fixed << temp_vol_score << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_VOS_geometric_vo:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_vol.total_component << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_VOS_hvy_atom_vo:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_vol.heavy_atom_component << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_VOS_pos_chrg_atm_vo:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_vol.positive_component << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_VOS_neg_chrg_atm_vo:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_vol.negative_component << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_VOS_hydrophobic_atm_vo:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_vol.hydrophobic_component << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_VOS_hydrophilic_atm_vo:"
                 << setw(FLOAT_WIDTH) << fixed << desc_c_vol.hydrophilic_component << endl;
        }

        if (desc_use_gist){
            text << DELIMITER << setw(STRING_WIDTH) << "Property_GIST_Score:"
                 << setw(FLOAT_WIDTH) << fixed << temp_gist_score << endl;
        }
/*
        if (desc_use_cmg){
            //text << DELIMITER << setw(STRING_WIDTH) << "Property_dock3.5_Score:"     << setw(FLOAT_WIDTH) << fixed << temp_cmg_score                    << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_Chemgrid_Score:"        << setw(FLOAT_WIDTH) << fixed << temp_cmg_score                    << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_Chemgrid_vdw_energy:"   << setw(FLOAT_WIDTH) << fixed << desc_c_cmg.vdw_component          << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_Chemgrid_es_energy:"    << setw(FLOAT_WIDTH) << fixed << desc_c_cmg.es_component           << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_Chemgrid_rec_desolv:"   << setw(FLOAT_WIDTH) << fixed << desc_c_cmg.rdsol_component        << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_Chemgrid_lig_polsol:"   << setw(FLOAT_WIDTH) << fixed << desc_c_cmg.polsolv_component      << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_Chemgrid_lig_apolsol:"  << setw(FLOAT_WIDTH) << fixed << desc_c_cmg.apolsolv_component     << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "desc_Chemgrid_conf_entropy:" << setw(FLOAT_WIDTH) << fixed << desc_c_cmg.conf_entropy_component << endl;

//            sprintf(line,
//                    "########## Chemgrid_Score:%20f\n##########%14s%20f\n##########%14s%20f\n##########%14s%20f\n##########%14s%20f\n##########%14s%20f##########%14s%20f\n##########\n",
//                    score, "Chemgrid_vdw_energy:", vdw_component, "Chemgrid_es_energy:", es_component,
//                    "Chemgrid_rec_desol:", rdsol_component,
//                    "Chemgrid_polsol:", polsolv_component, "Chemgrid_apolsol:",
//                    apolsolv_component, "Chemgrid_conf_entropy:",
//                    conf_entropy_component);
//
        }
*/

        if (use_internal_energy){
            text << DELIMITER << setw(STRING_WIDTH) << "Internal_energy_repulsive:"
                 << setw(FLOAT_WIDTH) << fixed << mol.internal_energy << endl;
        }
    }

/*
    if (use_score) {
        text << DELIMITER << setw(STRING_WIDTH) << "Rot Bonds:"
             << setw(FLOAT_WIDTH) << fixed << mol.rot_bonds << endl
             << DELIMITER << setw(STRING_WIDTH) << "Heavy Atoms:"
             << setw(FLOAT_WIDTH) << fixed << mol.heavy_atoms << endl
             << DELIMITER << setw(STRING_WIDTH) << "HB Donors:"
             << setw(FLOAT_WIDTH) << fixed << mol.hb_donors << endl
             << DELIMITER << setw(STRING_WIDTH) << "HB Acceptors:"
             << setw(FLOAT_WIDTH) << fixed << mol.hb_acceptors << endl
             << DELIMITER << setw(STRING_WIDTH) << "Molecular Wt:"
             << setw(FLOAT_WIDTH) << fixed << mol.mol_wt << endl
             << DELIMITER << setw(STRING_WIDTH) << "Formal Charge:"
             << setw(FLOAT_WIDTH) << fixed << mol.formal_charge << endl
             << DELIMITER << setw(STRING_WIDTH) << "XLogP:"
             << setw(FLOAT_WIDTH) << fixed << mol.xlogp << endl
             << DELIMITER << setw(STRING_WIDTH) << "Total Score:"
             << setw(FLOAT_WIDTH) << fixed << mol.current_score

             << endl;
    }
*/

    return text.str();

} // end Descriptor_Energy_Score::output_score_summary()


