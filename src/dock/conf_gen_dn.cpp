#include "conf_gen_dn.h"
#include "conf_gen_ga.h"
#include "fingerprint.h"
#include "hungarian.h"
#include "gasteiger.h"
#include "iso_align.h"
#include "utils.h"
#include "score_volume.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <map>
#include <list>
#include <functional>
#include <unistd.h>
#include <cerrno>
#include <errno.h>
#include <string.h>
#include <random>
#include <sys/stat.h>
#include <chrono>

#include "trace.h"
//#include <cstdlib>
#include <time.h>

#ifdef BUILD_DOCK_WITH_RDKIT
#include "rdtyper.h"
#endif

#ifndef __APPLE__
#include "mem_info.h"
#endif

using namespace std;

class Bump_Filter;
class Fingerprint;
class AMBER_TYPER;
class Master_Score;

// static member initializers
// These parameters are consistent with those in Base_Score()
const string DN_Build::DELIMITER    = "########## ";
const int    DN_Build::FLOAT_WIDTH  = 20;
const int    DN_Build::STRING_WIDTH = 17 + 19;



//
// +++++++++++++++++++++++++++++++++++++++++
// Some constructors and destructors

TorEnv::TorEnv(){
    origin_env = "";
}
TorEnv::~TorEnv(){
    target_envs.clear();
    target_freqs.clear();
}


/*FragGraph::FragGraph(){
    visited = false;
}
FragGraph::~FragGraph(){
    tanvec.clear();
    rankvec.clear();
}
*/

DN_Build::DN_Build(){
    // Flag used to save & print mutated molecules for GA - CDS
    dn_ga_flag      = false;
    dn_iso_write_it = 0;

}
DN_Build::~DN_Build(){
    scaffolds.clear();
    linkers.clear();
    sidechains.clear();
    anchors.clear();
    scaf_link_sid.clear();
}



// +++++++++++++++++++++++++++++++++++++++++
// Read parameters from the dock.in file (called in master_conf.cpp)
void
DN_Build::input_parameters( Parameter_Reader & parm )
{
    cout <<endl <<"De Novo Build Parameters" <<endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    dn_mode = "sample";
    //dn_mode = parm.query_param("dn_mode", "sample", "sample | utilities");    

    //if (dn_mode.compare("utilities") == 0){ 

    //    dn_ult_file   = parm.query_param( "dn_utl_file",
    //                                                  "input.mol2" );
    //    dn_ult_method = parm.query_param("dn_ult_method", "best_first", "best_first");
  
    //    if ( dn_ult_method.compare("best_first") == 0 ){

    //        dn_ult_bf_anc_diverse = parm.query_param("dn_ult_bf_anc_diverse", "no", "yes no") == "yes";
    //        dn_ult_score = 
    //            parm.query_param("dn_ult_score", "hungarian", "hungarian | volume");
    //  
    //        if (dn_ult_score.compare("hungarian") == 0){

    //            dn_ult_hun_matching_coeff = 
    //                atof(parm.query_param("dn_ult_hun_matching_coeff", "-5").c_str());
    //            dn_ult_hun_rmsd_coeff = 
    //                atof(parm.query_param("dn_ult_hun_rmsd_coeff", "1").c_str());
    //            dn_ult_clust_tol = 
    //                atof(parm.query_param("dn_ult_clust_tol", "-3.0").c_str());
    //        } 
    //      
    //        if (dn_ult_score.compare("volume") == 0) {
    //            dn_ult_vol_type  = 
    //                parm.query_param("dn_ult_vol_type", "volume_heavy", "volume_heavy");
    //            dn_ult_clust_tol = 
    //                atof(parm.query_param("dn_ult_clust_tol", "0.80").c_str());
    //        }

    //        dn_ult_clust_size_lim = 
    //            atoi(parm.query_param("dn_ult_clust_size_lim", "2").c_str());
    //        dn_ult_num_clusters = 
    //            atoi(parm.query_param("dn_ult_num_clusters", "10000").c_str()); 
    //    }

    //    dn_sampling_method_ex = false;
    //    dn_sampling_method_rand = false;
    //    dn_sampling_method_graph = false;
    //    dn_sampling_method_matrix = false;
    //    dn_sampling_method_isoswap = false;
    //    dn_unique_anchors = 0; 
    //    verbose = false;
    //    dn_user_specified_anchor = false; 
    //    dn_max_grow_layers = 0;
 
    //    return;
    //}

    dn_fraglib_scaffold_file  = parm.query_param( "dn_fraglib_scaffold_file",
                                                  "fraglib_scaffold.mol2" );
    dn_fraglib_linker_file    = parm.query_param( "dn_fraglib_linker_file",
                                                  "fraglib_linker.mol2" );
    dn_fraglib_sidechain_file = parm.query_param( "dn_fraglib_sidechain_file",
                                                  "fraglib_sidechain.mol2" );
    // Uncomment if you want to include rigids in de novo
    // dn_fraglib_rigid_file     = parm.query_param( "dn_fraglib_rigid_file",
    //                                               "fraglib_rigid.mol2" );

    // User has the option to specify anchor(s)
    dn_user_specified_anchor = (parm.query_param
                               ("dn_user_specified_anchor", "yes", "yes no") == "yes");


    if (dn_user_specified_anchor){
        dn_fraglib_anchor_file = parm.query_param("dn_fraglib_anchor_file", "fraglib_anchor.mol2");
    }

    // denovo random refinement option LEP
    //dn_refinement_random = (parm.query_param
    //                       ("dn_refinement_random", "yes", "yes no") == "yes");
    //if (dn_refinement_random){
    //    dn_refinement_random_picks = atoi(parm.query_param("dn_refinement_random_picks", "3").c_str());
    //}

    // Ask the user for a torsion environment table of allowable environments
    // (The only time you might NOT want this is if you are doing simple_build. Otherwise it really
    // is useful in preventing some garbage from being created.)
    //dn_use_torenv_table = (parm.query_param("dn_use_torenv_table", "yes", "yes no") == "yes");
    // LEP 2019.03.13 hardcode torsion table on 
    dn_use_torenv_table = true;

    if (dn_use_torenv_table){
        dn_torenv_table = parm.query_param("dn_torenv_table", "fraglib_torenv.dat");

        //dn_use_roulette = parm.query_param("dn_use_roulette", "no", "yes no") == "yes";
        dn_use_roulette = false;
    }
    
    denovo_name = parm.query_param( "dn_name_identifier", "denovo" );

    // Specify a specific dn sampling method
    // ex = exhaustive (slow)
    // rand = random (user specified limit)
    // graph = sample from graph of related fragments
    //dn_sampling_method = parm.query_param("dn_sampling_method", "graph", "ex | rand | graph | matrix | isoswap");
    //JDB removed matrix being displayed here for 6.12
    //dn_sampling_method = parm.query_param("dn_sampling_method", "graph", "ex | rand | graph | matrix ");
    dn_sampling_method = parm.query_param("dn_sampling_method", "graph", "ex | rand | graph");
    dn_sampling_method_ex = false;
    dn_sampling_method_rand = false;
    dn_sampling_method_graph = false;
    dn_sampling_method_matrix = false;
    dn_sampling_method_isoswap = false;

    if (dn_sampling_method.compare("ex") == 0){
        dn_sampling_method_ex = true;
    } else if (dn_sampling_method.compare("rand") == 0){
        dn_sampling_method_rand = true;
    } else if (dn_sampling_method.compare("graph") == 0){
        dn_sampling_method_graph = true;
    //JDB Disabled matrix completely - defaults to random.
    } else if (dn_sampling_method.compare("matrix") == 0){
        cout <<  "Matrix disabled in this version - reverting to RAND" << endl;
        dn_sampling_method_matrix = false;
        dn_sampling_method_rand = true;
        //dn_sampling_method_matrix = true;
        //dn_frequencies_matrix_file = parm.query_param("dn_frequencies_matrix_file", "fraglib_matrix.dat");
    } else if (dn_sampling_method.compare("isoswap") == 0){
        dn_sampling_method_isoswap = true;
    } else {
        cout <<"You chose...poorly." <<endl;
        exit(0);
    }

    // If you are using the random sampling method, limit the number of picks at each connection
    // point
    if (dn_sampling_method_rand){
        dn_num_random_picks = atoi(parm.query_param("dn_num_random_picks", "20").c_str());
    }

    // Lots of parameters if you are doing the graph sampling method
    if (dn_sampling_method_graph){
        // The number of random starting points in the graph
        dn_graph_max_picks = atoi(parm.query_param("dn_graph_max_picks", "30").c_str());

        // The number of children explored relative to a starting point
        dn_graph_breadth = atoi(parm.query_param("dn_graph_breadth", "3").c_str());

        // The number of generations, e.g. 2=explore children, then children's children
        dn_graph_depth = atoi(parm.query_param("dn_graph_depth", "2").c_str());

        // A starting value for the graph temp [ P=exp(-(Ef-Ei)/T) ]
        dn_graph_temperature = atof(parm.query_param("dn_graph_temperature", "100.0").c_str());
        dn_temp_begin = dn_graph_temperature;
    }
    //Exhaustive method uses the random sampling method with a high number of picks
    if (dn_sampling_method_ex){
        dn_num_random_picks = atoi(parm.query_param("dn_num_ex_picks", "200").c_str());
    }

    if (!dn_sampling_method_isoswap){
        dn_bias_with_fraglib = false;
        //dn_bias_with_fraglib = parm.query_param("dn_bias_with_fragments", "no", "yes | no") == "yes";
        if ( dn_bias_with_fraglib && !dn_sampling_method_isoswap) {
            dn_sel_frag_by_freq_bool = parm.query_param("dn_sel_frag_by_freq", "no", "yes | no") == "yes";
            dn_acc_frag_by_freq_bool = parm.query_param("dn_accept_frag_by_freq", "no", "yes | no") == "yes";
            if ((dn_sel_frag_by_freq_bool) || (dn_acc_frag_by_freq_bool)) {
                dn_frag_frequency_file = parm.query_param("dn_frag_frequency_file", "frag_freqs.dat");
            } else {
                cout << "Warning: Biasing with fragment selected, but no method specified.\nFile not read in and no biasing will be performed." << std::endl;
            }
        } else {
            dn_sel_frag_by_freq_bool = false;
            dn_acc_frag_by_freq_bool = false;
        }

    // Isoswap doesn't use biasing
    } else { 
        dn_sel_frag_by_freq_bool = false;
        dn_acc_frag_by_freq_bool = false;
    }

    if (dn_sampling_method_isoswap){	
        dn_iso_write_libraries        = parm.query_param("dn_iso_write_libraries", "no", "yes no").c_str();
      
        if (dn_iso_write_libraries.compare("yes")==0) {
	    dn_iso_rank           = parm.query_param("dn_iso_rank", "yes", "yes no").c_str();
            if (dn_iso_rank.compare("yes")==0) {
	        dn_iso_reverse        = parm.query_param("dn_iso_reverse", "no", "yes no").c_str();
                dn_iso_rank_score_sel = parm.query_param("dn_iso_rank_score_sel", "heavy_atom_component").c_str();
                dn_iso_num_top        = atoi(parm.query_param("dn_iso_num_top_rank", "100").c_str());
            } else {
                dn_iso_reverse        = "";
                dn_iso_rank_score_sel = "";
                dn_iso_num_top        = atoi(parm.query_param("dn_iso_num_top","100").c_str());
            }
            dn_iso_score_sel      = parm.query_param("dn_iso_score_sel", "heavy_atom_component",
                                                     "heavy_atom_component | hms_score").c_str();
            dn_iso_output_path_libraries  = parm.query_param("dn_iso_output_path_libraries", "./").c_str();

            dn_fraglib_iso_scaffold_file  = parm.query_param("dn_fraglib_iso_scaffold_file",
                                                                "fraglib_iso_iso_scaffold.mol2");
            dn_fraglib_iso_linker_file    = parm.query_param("dn_fraglib_iso_linker_file",
                                                                "fraglib_iso_linker.mol2");
            dn_fraglib_iso_sidechain_file = parm.query_param("dn_fraglib_iso_sidechain_file",
                                                                "fraglib_iso_sidechain.mol2");

            dn_iso_bond_angle_tol_sid     = atof(parm.query_param("dn_iso_bond_angle_tol_sid",
            								                        "5.0" ).c_str());
            dn_iso_bond_angle_tol_lnk     = atof(parm.query_param("dn_iso_bond_angle_tol_lnk",
                                                                    "5.0" ).c_str());
            dn_iso_bond_angle_tol_scf     = atof(parm.query_param("dn_iso_bond_angle_tol_scf",
                                                                    "5.0" ).c_str());
            
            dn_iso_dist_du_du_inter       = atof(parm.query_param( "dn_iso_dist_du_du_inter",
                                                                    "0.25" ).c_str());

            dn_iso_dist_tol_sid           = atof(parm.query_param("dn_iso_dist_tol_sid",
                                                                    "1.0" ).c_str());
            dn_iso_dist_tol_lnk           = atof(parm.query_param("dn_iso_dist_tol_lnk",
                                                                    "1.0" ).c_str());
            dn_iso_dist_du_du_lnk         = atof(parm.query_param("dn_iso_dist_du_du_lnk",
                                                                    "1.0" ).c_str());
            dn_iso_dist_du_du_scf         = atof(parm.query_param("dn_iso_dist_du_du_scf",
                                                                    "1.0" ).c_str());
            dn_iso_diff_num_atoms         = atoi(parm.query_param( "dn_iso_diff_num_atoms",
                                                                    "5" ).c_str());
            dn_iso_freq_cutoff            = atoi(parm.query_param("dn_iso_freq_cutoff", "130").c_str());
            dn_iso_cos_score_cutoff       = atof(parm.query_param("dn_iso_cos_score_cutoff", "0.0").c_str());
        } else if (dn_iso_write_libraries.compare("no") == 0){
            dn_num_rand_head_picks        = atoi(parm.query_param("dn_num_rand_head_picks", "20").c_str());
            dn_iso_fraglib                = parm.query_param( "dn_iso_fraglib", "almade", "almade | onthefly ").c_str();
            dn_iso_pick_meth              = parm.query_param("dn_iso_pick_meth", "top_rank", "top_rank | rand | worst_rank").c_str();
            dn_iso_skip                   = parm.query_param("dn_iso_skip", "yes", "yes | no").c_str();
            dn_num_iso_picks              = atoi(parm.query_param("dn_num_iso_picks", "20").c_str());
            dn_iso_print_out              = parm.query_param("dn_iso_print_out", "no", "yes | no").c_str();
            if (dn_iso_print_out.compare("yes") == 0) { 
                dn_iso_output_verbose_path =  parm.query_param("dn_iso_output_verbose_path", "./").c_str();
            }

            if (dn_iso_fraglib.compare("almade")==0){
                dn_iso_fraglib_dir    = parm.query_param("dn_iso_fraglib_dir", "isofrags").c_str();
                dn_iso_num_gets       = atoi(parm.query_param("dn_iso_num_gets","100").c_str());

            } else if (dn_iso_fraglib.compare("onthefly")==0){
                dn_iso_rank           = "yes";

                dn_iso_reverse        = parm.query_param("dn_iso_reverse", "no", "yes no").c_str();
                dn_iso_rank_score_sel = parm.query_param("dn_iso_rank_score_sel","heavy_atom_component").c_str();
                dn_iso_score_sel      = parm.query_param("dn_iso_score_sel","heavy_atom_component", 
                                                         "heavy_atom_component | hms_score").c_str();
                dn_iso_num_top        = atoi(parm.query_param("dn_iso_num_top_rank","100").c_str());

                dn_fraglib_iso_scaffold_file  = parm.query_param("dn_fraglib_iso_scaffold_file",
                                                                  "fraglib_iso_iso_scaffold.mol2");
                dn_fraglib_iso_linker_file    = parm.query_param("dn_fraglib_iso_linker_file",
                                                                  "fraglib_iso_linker.mol2");
                dn_fraglib_iso_sidechain_file = parm.query_param("dn_fraglib_iso_sidechain_file",
                                                                      "fraglib_iso_sidechain.mol2");
                dn_iso_bond_angle_tol_sid     = atof(parm.query_param("dn_iso_bond_angle_tol_sid",
                                                                        "5.0" ).c_str());
                dn_iso_bond_angle_tol_lnk     = atof(parm.query_param("dn_iso_bond_angle_tol_lnk",
                                                                        "5.0" ).c_str());
                dn_iso_bond_angle_tol_scf     = atof(parm.query_param("dn_iso_bond_angle_tol_scf",
                                                                        "5.0" ).c_str());
                dn_iso_dist_du_du_inter       = atof(parm.query_param("dn_iso_dist_du_du_inter",
                                                                        "0.25" ).c_str());
                dn_iso_dist_tol_sid           = atof(parm.query_param("dn_iso_dist_tol_sid",
                                                                        "1.0" ).c_str());
                dn_iso_dist_tol_lnk           = atof(parm.query_param("dn_iso_dist_tol_lnk",
                                                                        "1.0" ).c_str());
                dn_iso_dist_du_du_lnk         = atof(parm.query_param("dn_iso_dist_du_du_lnk",
                                                                        "1.0" ).c_str());
                dn_iso_dist_du_du_scf         = atof(parm.query_param("dn_iso_dist_du_du_scf",
                                                                        "1.0" ).c_str());
                dn_iso_diff_num_atoms         = atoi(parm.query_param("dn_iso_diff_num_atoms",
             							                                "5" ).c_str());
                dn_iso_freq_cutoff            = atoi(parm.query_param("dn_iso_freq_cutoff", "130").c_str());
                dn_iso_cos_score_cutoff       = atof(parm.query_param("dn_iso_cos_score_cutoff", "0.0").c_str());
            }
        }
    }

    // User specifies a hard upper score cutoff, similar to pruning_conformer_score_cutoff in 
    // anchor-and-grow
    dn_pruning_conformer_score_cutoff = atof(parm.query_param("dn_pruning_conformer_score_cutoff",
                                                              "100.0").c_str());
    dn_pruning_conformer_score_cutoff_begin = dn_pruning_conformer_score_cutoff;
    dn_pruning_conformer_score_scaling_factor = atof(parm.query_param("dn_pruning_conformer_score_scaling_factor", "2.0").c_str());

    // User specifies a cutoff for pruning torsions - equivalent to the anchor and grow heuristic
    dn_pruning_clustering_cutoff = atof(parm.query_param("dn_pruning_clustering_cutoff",
                                                         "100.0").c_str());
    //2022.08.05 Brock Boysan
    //User specifies whether to prune identical molecules from the ensemble
    dn_make_unique = (parm.query_param("dn_remove_duplicates","yes", " no | yes") == "yes");

    //if make unique specify number of duplicates allowed, duplicates per molecule, and whether to dump them out to a duplicate prune file
    if (dn_make_unique) {
        dn_max_duplicates_per_molecule = atoi(parm.query_param("dn_max_duplicates_per_mol", "0").c_str());
        dn_write_out_duplicates = (parm.query_param("dn_write_pruned_duplicates", "no", " no | yes ") == "yes");
    }

    dn_advanced_pruning = parm.query_param("dn_advanced_pruning", "yes", "yes | no") == "yes";
    if (dn_advanced_pruning) {
        if (!dn_sampling_method_isoswap) {
            dn_prune_initial_sample = parm.query_param("dn_prune_initial_sample", "yes", "yes | no") == "yes";
        } else {
            dn_prune_initial_sample = false;
        }
        dn_legacy_prune_hrmsd = false;
        //dn_legacy_prune_hrmsd = parm.query_param("dn_legacy_prune_hrmsd", "no", "yes | no") == "yes";
        if (!dn_sampling_method_isoswap) {
            //dn_legacy_prune_complete = parm.query_param("dn_legacy_prune_complete", "no", "yes | no") == "yes";
            dn_legacy_prune_complete = false;
        } else {
            dn_legacy_prune_complete = true; 
        }

        dn_sample_torsions = parm.query_param("dn_sample_torsions", "yes", "yes | no") == "yes";
        if (dn_sample_torsions && !dn_sampling_method_isoswap) {
            dn_prune_individual_torsions = parm.query_param("dn_prune_individual_torsions", "yes", "yes | no") == "yes";
            dn_prune_combined_torsions = parm.query_param("dn_prune_combined_torsions", "yes", "yes | no") == "yes";
        } else if (!dn_sample_torsions) {
            cout << "WARNING: Torsional sampling has been DISABLED. This should only be done during the SIMPLE BUILD protocol!" << endl;
            dn_prune_individual_torsions = false;
            dn_prune_combined_torsions = false;            
        } else {
            dn_prune_individual_torsions = false;
            dn_prune_combined_torsions = false;
        }
        dn_random_root_selection = parm.query_param("dn_random_root_selection", "no", "yes | no") == "yes";
    } else {
        dn_sample_torsions = true;
        dn_prune_initial_sample = false;
        dn_prune_individual_torsions = false;
        dn_prune_combined_torsions = false;
        dn_random_root_selection = false;
        dn_legacy_prune_hrmsd = true;
        dn_legacy_prune_complete = true;
    }
    // User can specify whether they want a hard or soft MW cutoff
    // In both cases the upper bound is assessed during growth, while the lower bound is assesed
    // before molecules are written to file so as not to remove viable molecules from the growing ensemble
    dn_MW_cutoff_type_hard = false;
    dn_MW_cutoff_type = parm.query_param("dn_mol_wt_cutoff_type", "soft", "hard | soft");
    dn_upper_constraint_mol_wt = atof(parm.query_param("dn_upper_constraint_mol_wt", "550.0").c_str());
    dn_lower_constraint_mol_wt = atof(parm.query_param("dn_lower_constraint_mol_wt", "0.0").c_str());

    // User can specify an upper bound/lower bound for mw and if a soft cutoff also standard deviation
    if (dn_MW_cutoff_type.compare("hard") == 0) {
        dn_MW_cutoff_type_hard = true;
    } else if (dn_MW_cutoff_type.compare("soft") == 0) {
        dn_MW_std_dev = atof(parm.query_param("dn_mol_wt_std_dev", "35.0").c_str());
    } else {
        cout <<"You chose...poorly." <<endl;
        exit(0);
    }

    #ifdef BUILD_DOCK_WITH_RDKIT

    dn_drive_verbose = (parm.query_param("dn_drive_verbose","no") == "yes");
    if (dn_drive_verbose){
        dn_save_all_molecules = (parm.query_param("dn_save_all_mols","no") == "yes");
    } else {
        dn_save_all_molecules = false;
    }

    dn_drive_clogp = (parm.query_param("dn_drive_clogp","no") == "yes");
    if (dn_drive_clogp){
        dn_lower_clogp = atof(parm.query_param("dn_lower_clogp","-0.30").c_str());
        dn_upper_clogp = atof(parm.query_param("dn_upper_clogp","3.75").c_str());
        dn_clogp_std_dev = atof(parm.query_param("dn_clogp_std_dev","2.02").c_str());
    }
    
    dn_drive_esol = (parm.query_param("dn_drive_esol","no") == "yes");
    if (dn_drive_esol){
        dn_lower_esol = atof(parm.query_param("dn_lower_esol","-5.23").c_str());
        dn_upper_esol = atof(parm.query_param("dn_upper_esol","-1.35").c_str());
        dn_esol_std_dev = atof(parm.query_param("dn_esol_std_dev","1.94").c_str());
    }

    //TPSA 
    dn_drive_tpsa = (parm.query_param("dn_drive_tpsa","no") == "yes");
    if (dn_drive_tpsa){
        dn_lower_tpsa = atof(parm.query_param("dn_lower_tpsa","28.53").c_str());
        dn_upper_tpsa = atof(parm.query_param("dn_upper_tpsa","113.20").c_str());
        dn_tpsa_std_dev = atof(parm.query_param("dn_tpsa_std_dev","42.33").c_str());
    }
 
    dn_drive_qed = (parm.query_param("dn_drive_qed","no") == "yes");
    if (dn_drive_qed){
        dn_lower_qed = atof(parm.query_param("dn_lower_qed","0.61").c_str());
        dn_qed_std_dev = atof(parm.query_param("dn_qed_std_dev","0.19").c_str());
    }
    
    dn_drive_sa = (parm.query_param("dn_drive_sa","no") == "yes");
    if (dn_drive_sa){
        dn_upper_sa = atof(parm.query_param("dn_upper_sa","3.34").c_str());
        dn_sa_std_dev = atof(parm.query_param("dn_sa_std_dev","0.9").c_str());
        if (dn_upper_sa == 0.0) {
            std::cout << "ERROR: Lowest score for SynthA is not 0.0, please try again"<< std::endl;
            exit(0);
        }
    }
    
    dn_drive_stereocenters = (parm.query_param("dn_drive_stereocenters","no") == "yes");
    if (dn_drive_stereocenters){
        dn_upper_stereocenter = atoi(parm.query_param("dn_upper_stereocenter","2").c_str());
        if (dn_upper_stereocenter == 0) {
            std::cout << "ERROR: Upperbound for stereocenters during DN_drive cannot be 0, please try again" <<std::endl;
            exit(0);
        }
    }
    
    dn_drive_pains = (parm.query_param("dn_drive_pains","no") == "yes");
    if (dn_drive_pains){
        dn_upper_pains = atoi(parm.query_param("dn_upper_pains","1").c_str());
        if (dn_upper_pains == 0) {
            std::cout << "ERROR: Upperbound for number of pains during DN_drive cannot be 0, please try again" << std::endl;
            exit(0);
        }
    }    
    

    if (dn_drive_clogp || dn_drive_esol || dn_drive_qed ||
        dn_drive_sa || dn_drive_stereocenters || dn_drive_pains || dn_drive_tpsa){
        dn_start_at_layer = atoi(parm.query_param("dn_start_at_layer","1").c_str());

        // Add fragMap to memory once
        sa_fraglib_path = parm.query_param("sa_fraglib_path","sa_fraglib.dat");

        if (fragMap.empty() == true) {
            std::ifstream fin(sa_fraglib_path);
            double key;
            double val;
            while (fin >> key >> val) {
                fragMap[key] = val;
            }
            //the fragMap is empty so quit 
            if (fragMap.empty() == true){
                std::cout << "sa_fraglib.dat is empty. Please give the correct fraglib data table...";
                exit(0);
            }
            fin.close();
        }

        // Add PAINSMap to memory once
        PAINS_path = parm.query_param("PAINS_path","pains_table_2019_09_01.dat");

        std::vector<std::string> PAINStmp;
        if (PAINStmp.empty() == true) {
            std::ifstream fin(PAINS_path);
            std::string tmp_string;
            while (fin >> tmp_string) {
                PAINStmp.push_back(tmp_string);
            }
            fin.close();
        }
        //if the PAINStmp is still empty, quit 
        if (PAINStmp.empty() == true){
            std::string tmp_string;
            std::cout << "pains_table.dat is empty. Please give the correct PAINS table...";
            exit(0);
        }

        //processing the pains maps.
        int numofvectors = PAINStmp.size() - 1;
        for (int i{0}; i < numofvectors; ) {
            PAINSMap[PAINStmp[i]] = " " + PAINStmp[i + 1];
            i += 2;
        }
        PAINStmp.clear();
    }

    #endif

    // User can specify an upper bound for rotatable bonds
    dn_constraint_rot_bon = atoi(parm.query_param("dn_constraint_rot_bon", "15").c_str());

    // User can specify a bound for total net charge (+ or -)
    dn_constraint_formal_charge = atof(parm.query_param("dn_constraint_formal_charge", "2.0").c_str());

    // User can specify a maximum number of charged groups (- and + treated the same)
    //dn_constraint_charge_groups = atoi(parm.query_param("dn_constraint_charge_groups", "10").c_str());

    // These two are for the RMSD horizontal pruning heuristic using Hungarian algorithm
    dn_heur_unmatched_num = atoi(parm.query_param("dn_heur_unmatched_num", "1").c_str());
    dn_heur_matched_rmsd = atof(parm.query_param("dn_heur_matched_rmsd", "2.0").c_str());

    // This is the number of unique anchors used to seed dn growth
    dn_unique_anchors = atoi(parm.query_param("dn_unique_anchors", "1").c_str());

    // Maximum number of layers to grow outward
    dn_max_grow_layers = atoi(parm.query_param("dn_max_grow_layers", "9").c_str());

    // Control over combinatorics, see "root" vector below
    dn_max_root_size = atoi(parm.query_param("dn_max_root_size", "25").c_str());

    // Control over combinatorics, see "layer" vector below
    dn_max_layer_size = atoi(parm.query_param("dn_max_layer_size", "25").c_str());

    // Maximum attachment points the molecule can have at any one time (in other words, stop adding
    // new scaffolds when the current fragment has this many open attachment points)
    dn_max_current_aps = atoi(parm.query_param("dn_max_current_aps", "5").c_str());

    // Control number of scaffolds that can be added per layer
    dn_max_scaffolds_per_layer = atoi(parm.query_param("dn_max_scaffolds_per_layer", "1").c_str());

    // Control the number of successful attachments per ROOT
    dn_max_sccful_att_per_root = atoi(parm.query_param("dn_max_successful_att_per_root", "50000").c_str());

    // Write checkpoint files at each layer, useful to monitor progress and/or restart jobs
    dn_write_checkpoints = (parm.query_param("dn_write_checkpoints", "yes", "yes no") == "yes");

    // Write the molecules that are pruned at each layer, useful to see what we are throwing away
    dn_write_prune_dump = (parm.query_param("dn_write_prune_dump", "yes", "yes no") == "yes");

    // Write the original orients to file
    dn_write_orients = (parm.query_param("dn_write_orients", "no", "yes no") == "yes");

    // Track molecule provenance by storing vectors of dockmols (growth trees)
    dn_write_growth_trees = (parm.query_param("dn_write_growth_trees", "no", "yes no") == "yes");

    if (dn_write_growth_trees){
        cout <<endl <<"Warning: Writing growth trees uses a lot of memory and writes a lot of"
             <<" files to disk.\n         We highly recommend carefully monitoring jobs that"
             <<" have this parameter turned on." <<endl <<endl;
    }

    // Output filename prefix
    dn_output_prefix = parm.query_param("dn_output_prefix", "output");

    // For verbose statistics
    verbose = 0 != parm.verbosity_level();

    return;

} // end DN_Build::input_parameters()

std::pair<std::vector<Iso_Acessory::Scored_Triangle>, std::vector<Fragment>>
DN_Build::frag_iso_align(Fragment & reffrag){

    class Iso_Align    iso_ali;

    //std::fstream frefout_molecules;
    //frefout_molecules.open("refsource.mol2",std::ios_base::out);
    ////for the testing fragments
    //std::fstream testfout_molecules;
    //testfout_molecules.open("output.mol2",std::ios_base::out);

    std::vector<Fragment> return_vec {};
    std::vector<Iso_Acessory::Scored_Triangle> return_tri {};
    std::vector<Iso_Acessory::Scored_Triangle> best_triangles {};
    std::pair<std::vector<Iso_Acessory::Scored_Triangle>, std::vector<Fragment>> return_pair;

    reffrag.calc_num_du();

    if (reffrag.get_num_du()==1){

	iso_ali.align(reffrag, isosidechains, iso_parm, best_triangles);
        if (best_triangles.size() != isosidechains.size()) {
            std::cout << "Exiting... triangles and the isosidechains are not same size" 
                      << std::endl;
            exit(0);
        }
 
       	//activate each heavy atom 
       	for (unsigned int i =0; i<reffrag.mol.num_atoms;i++){
       	    reffrag.mol.atom_active_flags[i] = true;
       	}

       	//Write_Mol2(anchors[0].mol,frefout_molecules);
       	for (unsigned int i=0;i<isosidechains.size();i++){
       	    for (unsigned int j =0; j<isosidechains[i].mol.num_atoms;j++){
       	        isosidechains[i].mol.atom_active_flags[j] = true;
       	    }
       	    if (isosidechains[i].is_it_iso_aligned()){
       	 	return_vec.push_back(isosidechains[i]);
                return_tri.push_back(best_triangles[i]);
       	 	isosidechains[i].is_not_iso_aligned();
       	    }
       	}
        return_pair = std::make_pair(return_tri,return_vec);

    }else if (reffrag.get_num_du()==2){

        iso_ali.align(reffrag, isolinkers, iso_parm, best_triangles);
        if (best_triangles.size() != isolinkers.size()) {
            std::cout << "Exiting... triangles and the isolinkers are not same size" 
                      << std::endl;
            exit(0);
        }

        //activate each heavy atom 
        for (unsigned int i =0; i<reffrag.mol.num_atoms;i++){
            reffrag.mol.atom_active_flags[i] = true;
        }

        //Write_Mol2(anchors[0].mol,frefout_molecules);
        for (unsigned int i=0;i<isolinkers.size();i++){
    

            for (unsigned int j =0; j<isolinkers[i].mol.num_atoms;j++){
                isolinkers[i].mol.atom_active_flags[j] = true;
            }
            if (isolinkers[i].is_it_iso_aligned()){
		return_vec.push_back(isolinkers[i]);
                return_tri.push_back(best_triangles[i]);
		isolinkers[i].is_not_iso_aligned();
                //Write_Mol2(isolinkers[i].mol,testfout_molecules);
            }
        }

        return_pair = std::make_pair(return_tri,return_vec);
    }else if (reffrag.get_num_du()>=3){

        iso_ali.align(reffrag, isoscaffolds,iso_parm,best_triangles);
        if (best_triangles.size() != isoscaffolds.size()) {
            std::cout << "Exiting... triangles and the isoscaffolds are not same size" 
                      << std::endl;
            exit(0);
        }
        //activate each heavy atom 
        for (unsigned int i =0; i<reffrag.mol.num_atoms;i++){
            reffrag.mol.atom_active_flags[i] = true;
        }

        //Write_Mol2(anchors[0].mol,frefout_molecules);
        for (unsigned int i=0;i<isoscaffolds.size();i++){
    

            for (unsigned int j =0; j<isoscaffolds[i].mol.num_atoms;j++){
                isoscaffolds[i].mol.atom_active_flags[j] = true;
            }
            if (isoscaffolds[i].is_it_iso_aligned()){
                Fragment tmp_frag;
                tmp_frag= isoscaffolds[i];
		return_vec.push_back(tmp_frag);
                return_tri.push_back(best_triangles[i]);

		isoscaffolds[i].is_not_iso_aligned();
		
                //Write_Mol2(isoscaffolds[i].mol,testfout_molecules);
            }
        }

        return_pair = std::make_pair(return_tri,return_vec);
    }

    return return_pair;

}


// +++++++++++++++++++++++++++++++++++++++++
// The initialize function tells the user the filenames it is reading from, and it also populates 
// the fragment library vectors. (called in master_conf.cpp)
void
DN_Build::initialize()
{
    Trace trace("DN_Build::initialize()");


    if ( dn_mode.compare("utilities") == 0 ){
        cout <<endl<<"Initializing De Novo Utilities Routines..." << endl;
   
        cout <<" Reading the utilties input library from " <<dn_ult_file <<"...";
        read_library(vec_utilities, dn_ult_file);
        cout <<"Done (#=" <<vec_utilities.size() <<")" <<endl;

        return;
    }




    cout <<endl <<"Initializing De Novo Growth Routines..." <<endl;

    cout <<" Reading the scaffold library from " <<dn_fraglib_scaffold_file <<"...";
    read_frag_library( scaffolds, dn_fraglib_scaffold_file );
    cout <<"Done (#=" <<scaffolds.size() <<")" <<endl;

    cout <<" Reading the linker library from " <<dn_fraglib_linker_file <<"...";
    read_frag_library( linkers, dn_fraglib_linker_file );
    cout <<"Done (#=" <<linkers.size() <<")" <<endl;
       
    cout <<" Reading the sidechain library from " <<dn_fraglib_sidechain_file <<"...";
    read_frag_library( sidechains, dn_fraglib_sidechain_file );
    cout <<"Done (#=" <<sidechains.size() <<")" <<endl;

    //set iso parameters
    //and load up mol2 files
    if (dn_sampling_method_isoswap){


	if (dn_iso_rank.compare("yes") == 0){
	    iso_parm.set_rank(true);
	}else{ 
	    iso_parm.set_rank(false);
	}

        if (dn_iso_write_libraries.compare("yes") == 0){
        }else{
            iso_parm.set_write_libraries(false);
        }

        if (dn_iso_write_libraries.compare("yes")==0) {
            iso_parm.set_write_libraries(true);
            iso_parm.set_iso_num_top(dn_iso_num_top);
            iso_parm.set_iso_score_sel(dn_iso_score_sel);
            iso_parm.set_bond_angle_tol_sid(dn_iso_bond_angle_tol_sid);
    	    iso_parm.set_bond_angle_tol_lnk(dn_iso_bond_angle_tol_lnk);
    	    iso_parm.set_bond_angle_tol_scf(dn_iso_bond_angle_tol_scf);
    	    iso_parm.set_dist_du_du_inter(dn_iso_dist_du_du_inter);
    	    iso_parm.set_dist_tol_sid(dn_iso_dist_tol_sid);
    	    iso_parm.set_dist_tol_lnk(dn_iso_dist_tol_lnk);
    	    iso_parm.set_dist_du_du_lnk(dn_iso_dist_du_du_lnk);
    	    iso_parm.set_dist_du_du_scf(dn_iso_dist_du_du_scf);
            iso_parm.set_diff_num_atoms(dn_iso_diff_num_atoms);
            iso_parm.set_iso_rank_score_sel(dn_iso_rank_score_sel);
            iso_parm.set_iso_rank_reverse(dn_iso_reverse.compare("yes")==0);
            iso_parm.set_iso_write_freq_cutoff(dn_iso_freq_cutoff);
            iso_parm.set_iso_cos_score_cutoff(dn_iso_cos_score_cutoff);

            cout <<" Reading the iso_scaffold library from " <<dn_fraglib_iso_scaffold_file <<"...";
    	    read_frag_library( isoscaffolds, dn_fraglib_iso_scaffold_file );
    	    cout <<"Done (#=" <<isoscaffolds.size() <<")" <<endl;

    	    cout <<" Reading the iso_linker library from " <<dn_fraglib_iso_linker_file <<"...";
    	    read_frag_library( isolinkers, dn_fraglib_iso_linker_file );
    	    cout <<"Done (#=" <<isolinkers.size() <<")" <<endl;

    	    cout <<" Reading the iso_sidechain library from " <<dn_fraglib_iso_sidechain_file <<"...";
    	    read_frag_library( isosidechains, dn_fraglib_iso_sidechain_file );
    	    cout <<"Done (#=" <<isosidechains.size() <<")" <<endl;
        } else if (dn_iso_write_libraries.compare("no")==0){
            iso_parm.set_write_libraries(false);

            if (dn_iso_fraglib.compare("almade")==0){
                iso_parm.set_iso_fraglib(true,dn_iso_fraglib_dir);
            } else if (dn_iso_fraglib.compare("onthefly")==0){
                std::string iso_fl_path(dn_iso_fraglib_dir);
                iso_parm.set_iso_fraglib(false,iso_fl_path);
                iso_parm.set_iso_num_top(dn_iso_num_top);
                iso_parm.set_iso_score_sel(dn_iso_score_sel);
	        iso_parm.set_bond_angle_tol_sid(dn_iso_bond_angle_tol_sid);
    	        iso_parm.set_bond_angle_tol_lnk(dn_iso_bond_angle_tol_lnk);
    	        iso_parm.set_bond_angle_tol_scf(dn_iso_bond_angle_tol_scf);
    	        iso_parm.set_dist_du_du_inter(dn_iso_dist_du_du_inter);
    	        iso_parm.set_dist_tol_sid(dn_iso_dist_tol_sid);
    	        iso_parm.set_dist_tol_lnk(dn_iso_dist_tol_lnk);
    	        iso_parm.set_dist_du_du_lnk(dn_iso_dist_du_du_lnk);
    	        iso_parm.set_dist_du_du_scf(dn_iso_dist_du_du_scf);
	        iso_parm.set_diff_num_atoms(dn_iso_diff_num_atoms);
                iso_parm.set_iso_rank_score_sel(dn_iso_rank_score_sel);
                iso_parm.set_iso_rank_reverse(dn_iso_reverse.compare("yes")==0);
                iso_parm.set_iso_write_freq_cutoff(dn_iso_freq_cutoff);
                iso_parm.set_iso_cos_score_cutoff(dn_iso_cos_score_cutoff);

	        cout <<" Reading the iso_scaffold library from " <<dn_fraglib_iso_scaffold_file <<"...";
    	        read_frag_library( isoscaffolds, dn_fraglib_iso_scaffold_file );
    	        cout <<"Done (#=" <<isoscaffolds.size() <<")" <<endl;

    	        cout <<" Reading the iso_linker library from " <<dn_fraglib_iso_linker_file <<"...";
    	        read_frag_library( isolinkers, dn_fraglib_iso_linker_file );
    	        cout <<"Done (#=" <<isolinkers.size() <<")" <<endl;

    	        cout <<" Reading the iso_sidechain library from " <<dn_fraglib_iso_sidechain_file <<"...";
    	        read_frag_library( isosidechains, dn_fraglib_iso_sidechain_file );
    	        cout <<"Done (#=" <<isosidechains.size() <<")" <<endl;
            }
        }

    }

    // We don't need to include rigids in de novo growth, but this can be uncommented for testing
    /*
    cout <<" Reading the rigid library from " <<dn_fraglib_rigid_file <<"...";
    read_library( rigids, dn_fraglib_rigid_file );
    cout <<"Done (#=" <<rigids.size() <<")" <<endl;
    */

    if (dn_user_specified_anchor) {
        cout <<" Reading the anchor library from " <<dn_fraglib_anchor_file <<"...";
        read_frag_library( anchors, dn_fraglib_anchor_file );
        cout <<"Done (#=" <<anchors.size() <<")" <<endl;
    }
    if (dn_use_torenv_table) {
        cout <<" Reading the torenv table from " <<dn_torenv_table <<"...";
        read_torenv_table(dn_torenv_table);
    }
    if (dn_use_roulette) {
        cout << " Generating the roulette wheel from the torenv table ... ";
        generate_roulette();
    }
    //JDB - reads in the fragment attachment matrix
    //also generates a relative matrix
    if(dn_sampling_method_matrix){
        cout << "Reading in the fragment library matrix from " << dn_frequencies_matrix_file <<"...";
        read_frequencies_matrix(dn_frequencies_matrix_file);


        //JDB Debugging statement to see if the fragment matrix was read in correctly
        /*
        for (int i=0; i<fragment_frequency_matrix.size(); ++i){
            for(int x=0; x<fragment_frequency_matrix[i].size(); ++x){
                cout << fragment_frequency_matrix[i][x] << ",";
            }
            cout << endl << endl;
        }
        */

    }
    // We used to do this for just graph, but now we do it for all sampling functions for simplicity

    if (dn_sampling_method_graph || dn_sampling_method_rand || 
        dn_sampling_method_ex    || dn_sampling_method_matrix || dn_sampling_method_isoswap){
        // put scaffolds, linkers, and sidechains in the same vector
        for (int i=0; i<scaffolds.size(); i++){ scaf_link_sid.push_back(scaffolds[i]); }
        for (int i=0; i<linkers.size(); i++){ scaf_link_sid.push_back(linkers[i]); }
        for (int i=0; i<sidechains.size(); i++){ scaf_link_sid.push_back(sidechains[i]); }
    }

    // BCF compute MW and formal charge of each fragment (stored in frag.mol) for easy pruning later
    for (int i=0; i<scaf_link_sid.size(); i++) {
          // Declare a temporary fingerprint object and compute atom environments (necessary for Gasteiger)
          activate_mol(scaf_link_sid[i].mol);
          Fingerprint temp_finger;
          for (int j=0; j<scaf_link_sid[i].mol.num_atoms; j++){
               scaf_link_sid[i].mol.atom_envs[j] = temp_finger.return_environment(scaf_link_sid[i].mol, j);
          }

          scaf_link_sid[i].mol.prepare_molecule();
          // Compute the charges, saving them on the mol object
          float total_charges = 0;
          total_charges = compute_gast_charges(scaf_link_sid[i].mol);
          calc_mol_wt(scaf_link_sid[i].mol);
          calc_formal_charge(scaf_link_sid[i].mol);
          calc_rot_bonds(scaf_link_sid[i].mol);
          scaf_link_sid[i].mol.rot_bonds = 0;
    }

    if (scaf_link_sid.size() > 0) {
        frag_sort(scaf_link_sid, fragment_sort_by_name);
    }

    // This function will find what fragments are similar to what other fragments in each of the
    // libraries, storing that information as a FragGraph object
    if (dn_sampling_method_graph){
        prepare_fragment_graph( scaf_link_sid, scaf_link_sid_graph );
        // To visualize or inspect graph, uncomment this. Note: exits program on finishing.
        //print_fraggraph();
    }

    if ((dn_sel_frag_by_freq_bool) || (dn_acc_frag_by_freq_bool)) {
        cout << "Reading in fragment frequencies..." << endl;
        read_frag_frequencies(dn_frag_frequency_file);
        cout << "Fragment frequencies read in!" << endl;
        //The check is (-1) because of placeholder value @ beginning
        if ((ordered_fragments.size() - 1) < (scaffolds.size() + linkers.size() + sidechains.size())) {
            cout << "WARNING: Number of fragments != size of fragment number line! May be using incorrect library. Exiting!" << endl;
            exit(0);
        }
    }
    return;

} // end DN_Build::initialize()



// +++++++++++++++++++++++++++++++++++++++++
// Get the internal energy parameters from Master_Conformer_Search
void
DN_Build::initialize_internal_energy_parms( bool uie, int rep_exp, int att_exp, float diel, float iec )
{
    Trace trace("DN_Build::initialize_internal_energy_parms()");
    // These just need to be defined before we try to call the minimizer

    use_internal_energy = uie;         // boolean, use internal energy in ligand minimization?
    ie_rep_exp          = rep_exp;     // int. energy vdw repulsive exponent (default 12)
    ie_att_exp          = att_exp;     // int. energy vdw attractive exponent (6)
    ie_diel             = diel;        // int. energy dielectric (4.0)
  //internal_energy_cutoff = iec;
    ie_cutoff           = iec;         // BCF internal energy cutoff 
    return;

} // end DN_Build::initialize_internal_energy_parms()

// +++++++++++++++++++++++++++++++++++++++++
// // Read anchor libraries from file and save in appropriate vector of <Fragments>
// // LEP function for anchors only for random refinement functionality
void
DN_Build::read_library_anchor( vector <Fragment> & frag_vec, string filename )
{
    //Create temporary DOCKMol obj
    DOCKMol tmp_frag_mol;

    //Create filehandle for reading mol2 library
    ifstream fin_frags;

    //Open filehandle for reading an anchor
    fin_frags.open(filename.c_str());
    if (fin_frags.fail()){
        cout <<"Could not open " <<filename <<". Program will terminate." <<endl;
        cout << strerror(errno);
        exit(0);
    }

    //Read the mol2 file and copy DOCKMol objs into vec Frag obj
    while ( Read_Mol2(tmp_frag_mol, fin_frags, false, false, false) ){

        //Create temp Frag obj
        Fragment tmp_frag;

        //Copy the DOCKMol obj into DOCKMol member of Frag obj
        tmp_frag.mol = tmp_frag_mol;

        //Create temp vecs of atom indices
        vector <int> tmp_dummy_atoms;
        vector <int> tmp_heavy_atoms;

        //Create temp pair for bond info for incase of GA
        std::pair <int, std::string> bond_info;


        // LEP exchanging H for Du functionality
        INTVec h_index;

        // if random refinement is on
        if (dn_refinement_random){
            cout << "dn refinement true" << endl;
            //for each random pick per mol
            for (int k=0; k < dn_refinement_random_picks; k++){
                cout << "heres a random pick" << endl;
                //Get number of H's
                for (int j=0; j< tmp_frag.mol.num_atoms; j++){
                    if ( tmp_frag.mol.atom_types[j] == "H" ){
                        h_index.push_back(j);
                        cout << "atom num that is H:" << j << endl;
                    }
                }
                //There needs to be H's for this function to work
                if ( h_index.size() > 0 ){
                    cout << "H-index greater than 0" << endl;
                    int sel_index=rand() % h_index.size();
                    convert_H_to_Du( tmp_frag.mol, h_index[sel_index]);
                    for (int j=0; j< h_index.size(); j++){
                            convert_H_to_Du( tmp_frag.mol, h_index[j]);
                    }
                }
                h_index.clear();
            }
        }
        //Loop over every bond in dockmol obj
         for (int i=0; i<tmp_frag.mol.num_bonds; i++){
            cout << "heres a bond" << endl;
            //CHeck if the bond origin is a dummy atom
            if (tmp_frag.mol.atom_types[tmp_frag.mol.bonds_origin_atom[i]].compare("Du") == 0){
                cout << "here's a Du" << endl;
                tmp_dummy_atoms.push_back(tmp_frag.mol.bonds_origin_atom[i]);
                tmp_heavy_atoms.push_back(tmp_frag.mol.bonds_target_atom[i]);

                //update some stuff for GA
                bond_info.first = i;
                bond_info.second = tmp_frag.mol.bond_types[i];
                tmp_frag.aps_bonds_type.push_back(bond_info);
            }

            // CHeck if the bond target is a Du
            if (tmp_frag.mol.atom_types[tmp_frag.mol.bonds_target_atom[i]].compare("Du") == 0){

                tmp_dummy_atoms.push_back(tmp_frag.mol.bonds_target_atom[i]);
                tmp_heavy_atoms.push_back(tmp_frag.mol.bonds_origin_atom[i]);

                //GA stuff
                bond_info.first = i;
                bond_info.second = tmp_frag.mol.bond_types[i];
                tmp_frag.aps_bonds_type.push_back(bond_info);
            }
        }
        // Create temp aps object
        AttPoint tmp_ap;

        //Loop over all Du/heavy atom pairs
        for (int i=0; i<tmp_dummy_atoms.size(); i++){

            // Copy the heavy atom/du pairs onto aps class
            tmp_ap.heavy_atom = tmp_heavy_atoms[i];
            tmp_ap.dummy_atom = tmp_dummy_atoms[i];
            tmp_ap.frag_name = tmp_frag.mol.energy;
            // Add that info to aps vector of the scaffold
            tmp_frag.aps.push_back(tmp_ap);
        }

        // Add complete Frag obj to the vector
        frag_vec.push_back(tmp_frag);
        cout << "anchor aps" << tmp_frag.aps[0].frag_name << endl;
        //Clear temp vecs
        tmp_dummy_atoms.clear();
        tmp_heavy_atoms.clear();

    }

    //Close filehandle
    fin_frags.close();

    return;

} // end DN_Build::read_library_anchor()

// DN_Build::read_library
void
DN_Build::read_library( vector <DOCKMol> & mol_vec, string filename ){

    // Create temporary DOCKMol object
    DOCKMol tmp_mol;
    
    // Create filehandle for reading mol2 library
    ifstream fin_frags;

    // Open the filehandle for reading a fragment library
    fin_frags.open(filename.c_str());
    if (fin_frags.fail()){
        cout <<"Could not open " <<filename <<". Program will terminate." <<endl;
        cout << strerror(errno);
        exit(0);
    }

    // Read the mol2 file and copy DOCKMol objects into a vector Fragment objects
    while ( Read_Mol2_retain(tmp_mol, fin_frags, false, false, false) ){
        mol_vec.push_back(tmp_mol);
    }
    

    return;
} // end DN_Build::read_library()


// +++++++++++++++++++++++++++++++++++++++++
// Read fragment libraries from file and save in appropriate vector of <Fragments>
void
DN_Build::read_library( vector <Fragment> & frag_vec, string filename )
{
    // Create temporary DOCKMol object
    DOCKMol tmp_frag_mol;

    // Create filehandle for reading mol2 library
    ifstream fin_frags;

    // Open the filehandle for reading a fragment library
    fin_frags.open(filename.c_str());
    if (fin_frags.fail()){
        cout <<"Could not open " <<filename <<". Program will terminate." <<endl;
        cout << strerror(errno);
        exit(0);
    }

    // Read the mol2 file and copy DOCKMol objects into a vector Fragment objects
    while ( Read_Mol2(tmp_frag_mol, fin_frags, false, false, false) ){

        // Create a temporary Fragment object
        Fragment tmp_frag;

        // Copy the DOCKMol object into the DOCKMol member of Fragment object
        tmp_frag.mol = tmp_frag_mol;

        // Create temporary vectors of atom indices
        vector <int> tmp_dummy_atoms;
        vector <int> tmp_heavy_atoms;
       
        // Create temp pair for bond information - CS: 09/19/16
        std::pair <int, std::string> bond_info;

        // Loop over every bond in the dockmol object
        for (int i=0; i<tmp_frag.mol.num_bonds; i++){

            // Check if the bond origin is a dummy atom
            if (tmp_frag.mol.atom_types[tmp_frag.mol.bonds_origin_atom[i]].compare("Du") == 0){

                // if so, save origin as dummy_atom and target as heavy atom
                tmp_dummy_atoms.push_back(tmp_frag.mol.bonds_origin_atom[i]);
                tmp_heavy_atoms.push_back(tmp_frag.mol.bonds_target_atom[i]);

                // Update the bond number and type - CS: 09/19/16
                // TODO: update sampling methods to use stored bond info for fraf selection (LEP)
                bond_info.first = i;
                bond_info.second = tmp_frag.mol.bond_types[i];
                tmp_frag.aps_bonds_type.push_back(bond_info);
            }

            // Check if the bond target is a dummy atom
            if (tmp_frag.mol.atom_types[tmp_frag.mol.bonds_target_atom[i]].compare("Du") == 0){

                // if so, save target as dummy atom and origin as heavy atom
                tmp_dummy_atoms.push_back(tmp_frag.mol.bonds_target_atom[i]);
                tmp_heavy_atoms.push_back(tmp_frag.mol.bonds_origin_atom[i]);

                // Update the bond number and type - CS: 09/19/16
                bond_info.first = i;
                bond_info.second = tmp_frag.mol.bond_types[i];
                tmp_frag.aps_bonds_type.push_back(bond_info);
            }
        }


        // Create temporary attachment point object 
        AttPoint tmp_ap;

        // Loop over all of the dummy_atom / heavy_atom pairs
        for (int i=0; i<tmp_dummy_atoms.size(); i++){

             // Copy the heavy atom / dummy atom pairs onto the attachment point class
             tmp_ap.heavy_atom = tmp_heavy_atoms[i];
             tmp_ap.dummy_atom = tmp_dummy_atoms[i];
             // And add that information to the attachment points (aps) vector of the scaffold
             tmp_frag.aps.push_back(tmp_ap);
        } 

        // Add the complete Fragment object to the appropriate vector
        frag_vec.push_back(tmp_frag); 

        // Clear temporary vectors
        tmp_dummy_atoms.clear();
        tmp_heavy_atoms.clear();

    } // end while( ReadMol2() ) loop

    // Close the filehandle
    fin_frags.close();

    return;

} // end DN_Build::read_library()



// +++++++++++++++++++++++++++++++++++++++++
// Updated Read fragment libraries from file and save in appropriate vector of <Fragments>
// planned to be the main function
void
DN_Build::read_frag_library( vector <Fragment> & frag_vec, string filename )
{

    // Create filehandle for reading mol2 library
    ifstream fin_frags;

    // Open the filehandle for reading a fragment library
    fin_frags.open(filename.c_str());
    if (fin_frags.fail()){
        cout <<"Could not open " <<filename <<". Program will terminate." <<endl;
        cout << strerror(errno);
        exit(0);
    }

    bool stop_or_no = true;

    DOCKMol mol;
    while (stop_or_no==true){
        char            line[1000];
        int             count;
        int             i;
        int             n1;
        bool            atom_line;
        char            typ[100],
                        col[100];

        char            tmp1[100],
                        tmp2[100],
                        jwu_tmp4[100],
                                jwu_tmp3[100],
                        subst_name[100];
        int             natoms,
                        nbonds,
                        nresidues; //defined by jwu
        string          l1,
                        l2,
                        l3,
                        l4,
                        l5,
                        l6;
        float           f1,
                        f2,
                        f3,
                        f4,
                        f5;
 
        int             freq        = 0;
        float           iso_score   = 0.0;
        std::string     three_pairs ="";

        bool            found_solvation = false;
        bool            found_color = false;
        bool            found_amber = false;

        bool            check_cols = false;
        bool            seven_columns = false; // default @<TRIPOS>ATOM line has 9 columns. By GDRM
        //int             bad_line{0}; // number of bad atom lines. By GDRM


        // init state vars
        atom_line = false;


        // read forward until the tripos molecule tag is reached
        for (;;) {
            if (!fin_frags.getline(line, 1000)) {
                mol.clear_molecule();
                // cout << endl << "ERROR: Ligand file empty.  Program will
                // terminate." << endl;
                stop_or_no = false;
                return;
            }

            if (!strncmp(line, "@<TRIPOS>MOLECULE", 17))
                break;
             

            std::string tmp_string = line;
            if (tmp_string.find("Iso_Score:")!=string::npos){
                std::istringstream iss (tmp_string);
                std::vector<std::string> string_vec{};
                std::string tmp_str;
                while(iss >> tmp_str){
                    string_vec.push_back(tmp_str);
                }
                if (string_vec.size() == 3) {
                    iso_score = stof(string_vec[2]);
                }
            }

            if (tmp_string.find("FREQ:")!=string::npos){
                std::istringstream iss (tmp_string);
                std::vector<std::string> string_vec{};
                std::string tmp_str;
                while(iss >> tmp_str){
                    string_vec.push_back(tmp_str);
                }
       
                if (string_vec.size() == 3) {
                    freq = stoi(string_vec[2]);
                }

            }

            if (tmp_string.find("Three_pairs:")!=string::npos){
                std::istringstream iss (tmp_string);
                std::vector<std::string> string_vec{};
                std::string tmp_str;
                while(iss >> tmp_str){
                    string_vec.push_back(tmp_str);
                }
       
                if (string_vec.size() == 3) {
                    three_pairs = string_vec[2];
                }

            }
        }
   

        // loop over the header info
        for (count = 0;; count++) {

            if (!fin_frags.getline(line, 1000)) {
                mol.clear_molecule();
                cout << endl <<
                    "ERROR:  Ligand file empty.  Program will terminate." << endl;
                stop_or_no = false;
                return;
            }

            if (!strncmp(line, "@<TRIPOS>ATOM", 13)) {
                atom_line = true;
                break;
            }
            // assign the first 5 header lines to the proper fields
            switch (count) {

            case 0:
                l1 = line;
                break;

            case 1:
                l2 = line;
                break;

            case 2:
                l3 = line;
                break;

            case 3:
                l4 = line;
                break;

            case 4:
                l5 = line;
                break;

            case 5:
                l6 = line;
                break;
            }

        }

        // if there are no atoms, throw error and return false
        if (!atom_line) {
            mol.clear_molecule();
            cout <<
                "ERROR: @<TRIPOS>ATOM indicator missing from ligand file.  Program will terminate."
                << endl;
            stop_or_no = false;
            return;
        }
        // get # of atoms and bonds from mol info line
        // sscanf(l2.c_str(), "%d %d", &natoms, &nbonds);
        sscanf(l2.c_str(), "%d %d %d", &natoms, &nbonds, &nresidues);

        // If nresidues is zero, DOCK6 needs to be told to read the @<TRIPOS>ATOM lines
        // correctly. By GDRM
        if (nresidues == 0){
            check_cols = true;
            nresidues = 1;
        }

        // initialize molecule vectors
        // mol.allocate_arrays(natoms, nbonds);
        mol.allocate_arrays(natoms, nbonds, nresidues);

        mol.title = l1;
        mol.mol_info_line = l2;
        mol.comment1 = l3;
        mol.comment2 = l4;
        mol.energy = l5;
        mol.comment3 = l6;
        // loop over atoms and read in atom info
        for (i = 0; i < mol.num_atoms; i++) {
            if (!fin_frags.getline(line, 1000)) {
                mol.clear_molecule();
                cout <<
                    "ERROR:  Atom information missing from ligand file.  Program will terminate."
                    << endl;
                stop_or_no = false;
                return;
            }
            //jwu comment
            // 1     N       15.2586  -59.3416   35.3528 N.4     1 PRO1  -0.2020
            // GDRM comment:
            // when nres = 0, @<TRIPOS>ATOM might contain 7 columns instead of 9:
            // 1     N       15.2586  -59.3416   35.3528 N.4  -0.2020
            if (check_cols){
                // Test number of columns per line
                stringstream string_test;
                string_test << line;
                int countCols{0};
                double value{};
                while(string_test >> value){
                    ++countCols;
                }
                if (countCols == 7){
                    seven_columns == true;
                }
            }
            //  Read data from mol2 as a function of the number of columns
            if (seven_columns){
                sscanf(line, "%s %s %f %f %f %s %f", jwu_tmp4, tmp1, &mol.x[i], &mol.y[i],
                       &mol.z[i], tmp2, &mol.charges[i]);
                mol.atom_names[i] = tmp1;
                mol.atom_types[i] = tmp2;
                mol.atom_number[i] = jwu_tmp4;
                mol.atom_residue_numbers[i] = "1";
                mol.subst_names[i] = subst_name;//mol.title
            } else {
                sscanf(line, "%s %s %f %f %f %s %s %s %f", jwu_tmp4, tmp1, &mol.x[i], &mol.y[i],
                       &mol.z[i], tmp2, jwu_tmp3, subst_name, &mol.charges[i]);
                mol.atom_names[i] = tmp1;
                mol.atom_types[i] = tmp2;
                    mol.atom_number[i] = jwu_tmp4;
                    mol.atom_residue_numbers[i] = jwu_tmp3;
                mol.subst_names[i] = subst_name;//mol.title
            }
        }
        // If all partial charges are zero, it is a bad molecule
        int num_zeros = 0;
        for (int i = 0; i < mol.num_atoms; ++i){
            if (mol.charges[i] == 0.0){
                num_zeros += 1;
            }
        }
        // Is it a bad molecule?
        if (num_zeros == mol.num_atoms){
            mol.bad_molecule = true;
        } else {
            mol.bad_molecule = false;
        }
        // skip down to the bond section
        for (;;) {
            if (!fin_frags.getline(line, 1000)) {
                mol.clear_molecule();
                cout <<
                    "ERROR: @<TRIPOS>BOND indicator missing from ligand file.  Program will terminate."
                    << endl;
                stop_or_no = false;
                return;
            }

            if (!strncmp(line, "@<TRIPOS>BOND", 13))
                break;
        }

        // loop over bonds and add them
        for (i = 0; i < mol.num_bonds; i++) {
            if (!fin_frags.getline(line, 1000)) {
                mol.clear_molecule();
                cout <<
                    "ERROR: Bond information missing from ligand file.  Program will terminate."
                    << endl;
                stop_or_no = false;
                return;
            }


            sscanf(line, "%*d %d %d %s", &mol.bonds_origin_atom[i],
                   &mol.bonds_target_atom[i], tmp1);

            // adjust bond atom #'s to start at 0
            mol.bonds_origin_atom[i]--;
            mol.bonds_target_atom[i]--;

            mol.bond_types[i] = tmp1;

        }


        // ID ring atoms/bonds
        mol.id_ring_atoms_bonds();

        // Create a temporary Fragment object
        Fragment tmp_frag;

        // Copy the DOCKMol object into the DOCKMol member of Fragment object
        tmp_frag.mol = mol;

        // Create temporary vectors of atom indices
        vector <int> tmp_dummy_atoms;
        vector <int> tmp_heavy_atoms;
       
        // Create temp pair for bond information - CS: 09/19/16
        std::pair <int, std::string> bond_info;

        // Loop over every bond in the dockmol object
        for (int i=0; i<tmp_frag.mol.num_bonds; i++){

            // Check if the bond origin is a dummy atom
            if (tmp_frag.mol.atom_types[tmp_frag.mol.bonds_origin_atom[i]].compare("Du") == 0){

                // if so, save origin as dummy_atom and target as heavy atom
                tmp_dummy_atoms.push_back(tmp_frag.mol.bonds_origin_atom[i]);
                tmp_heavy_atoms.push_back(tmp_frag.mol.bonds_target_atom[i]);


                // Update the bond number and type - CS: 09/19/16
                // TODO: update sampling methods to use stored bond info for fraf selection (LEP)
                bond_info.first = i;
                bond_info.second = tmp_frag.mol.bond_types[i];
                tmp_frag.aps_bonds_type.push_back(bond_info);
            }

            // Check if the bond target is a dummy atom
            if (tmp_frag.mol.atom_types[tmp_frag.mol.bonds_target_atom[i]].compare("Du") == 0){

                // if so, save target as dummy atom and origin as heavy atom
                tmp_dummy_atoms.push_back(tmp_frag.mol.bonds_target_atom[i]);
                tmp_heavy_atoms.push_back(tmp_frag.mol.bonds_origin_atom[i]);

                // Update the bond number and type - CS: 09/19/16
                bond_info.first = i;
                bond_info.second = tmp_frag.mol.bond_types[i];
                tmp_frag.aps_bonds_type.push_back(bond_info);
            }
        }


        // Create temporary attachment point object 
        AttPoint tmp_ap;

        // Loop over all of the dummy_atom / heavy_atom pairs
        for (int i=0; i<tmp_dummy_atoms.size(); i++){

             // Copy the heavy atom / dummy atom pairs onto the attachment point class
             tmp_ap.heavy_atom = tmp_heavy_atoms[i];
             tmp_ap.dummy_atom = tmp_dummy_atoms[i];
             // And add that information to the attachment points (aps) vector of the scaffold
             tmp_frag.aps.push_back(tmp_ap);
        }
        // JDB JDB JDB
        for(int i=0; i<tmp_frag.aps.size(); i++){
            tmp_frag.aps[i].frag_name = tmp_frag.mol.energy;
        }
        //fill in all information as much as possible for the frag clas
        tmp_frag.set_iso_score(iso_score);
        tmp_frag.freq_num    = freq;


        if (!three_pairs.empty()){
            //three_pair handling
            std::replace(three_pairs.begin(), three_pairs.end(), '/', ' ');
            vector<std::string> three_pairs_vec;
            std::stringstream ss(three_pairs);
            std::string temp;
            while (ss >> temp) { three_pairs_vec.push_back(temp); }
            //save the three pairs into the best_tri attribute
            Iso_Acessory::A_Pair first_pair;
            Iso_Acessory::A_Pair sec_pair;
            Iso_Acessory::A_Pair third_pair;

            first_pair.set_ref(std::stoi(three_pairs_vec[0]));
            first_pair.set_test(std::stoi(three_pairs_vec[1]));

            sec_pair.set_ref(std::stoi(three_pairs_vec[2]));
            sec_pair.set_test(std::stoi(three_pairs_vec[3]));

            third_pair.set_ref(std::stoi(three_pairs_vec[4]));
            third_pair.set_test(std::stoi(three_pairs_vec[5]));

            tmp_frag.best_tri.three_pairs.push_back(first_pair);
            tmp_frag.best_tri.three_pairs.push_back(sec_pair);
            tmp_frag.best_tri.three_pairs.push_back(third_pair);
        }
      
        // Add the complete Fragment object to the appropriate vector
        frag_vec.push_back(tmp_frag); 

        // Clear temporary vectors
        tmp_dummy_atoms.clear();
        tmp_heavy_atoms.clear();

        // Clear Molecules
        mol.clear_molecule();
    }

    // Close the filehandle
    fin_frags.close();

    return;

} // end DN_Build::read_frag_library()


// +++++++++++++++++++++++++++++++++++++++++
// This function reads the table of allowable environments and stores them in a data structure 
void
DN_Build::read_torenv_table( string torenv_table ) 
{
    // Open input filestream given the filename from parameter file 
    ifstream fin_torenv_table;
    fin_torenv_table.open(torenv_table.c_str());

    if (fin_torenv_table.fail()){
        cout <<"Could not open " <<torenv_table <<". Program will terminate." <<endl;
        exit(0);
    }

    // Declare some temporary variables
    TorEnv temp_torenv;
    string line;
    string current_origin = "";
    vector <string> temp_tokens;
    int counter = 0;
    
    // For each line in the torenv data file...
    while (!fin_torenv_table.eof()){

        // Get the next line
        getline (fin_torenv_table, line);

        // And split on the '-' character (this Tokenizer function is in utils.cpp)
        Tokenizer(line, temp_tokens, '-');

        if (line.size() == 0){
           continue;
        }

        // If the next token has previously been seen, add its pair to the growing vector of targets
        if (current_origin == temp_tokens[0]){
            temp_torenv.target_envs.push_back(temp_tokens[1]);
            temp_torenv.target_freqs.push_back(atoi(temp_tokens[2].c_str()));

        } else {  // If it hasn't been seen,

            // And if the data structure is not currently empty,
            if (temp_torenv.target_envs.size() > 0){

                // Add it to the vector of TorEnv data structures
                counter += temp_torenv.target_envs.size();
                torenv_vector.push_back(temp_torenv);
                temp_torenv.origin_env.clear();
                temp_torenv.target_envs.clear();
            }

            // Then begin populating a new TorEnv
            current_origin = temp_tokens[0];
            temp_torenv.origin_env = current_origin;
            temp_torenv.target_envs.push_back(temp_tokens[1]);
            temp_torenv.target_freqs.push_back(atoi(temp_tokens[2].c_str()));
        }
    }

    // Add the last data structure to the vector
    if (temp_torenv.target_envs.size() > 0){
        counter += temp_torenv.target_envs.size();
        torenv_vector.push_back(temp_torenv);
    }
    if (dn_ga_flag == false){
       cout <<"Done: (#=" <<counter <<")" <<endl;
    }

    return;

} // end DN_Build::read_torenv_table()
// +++++++++++++++++++++++++++++++++++++++++
// This function reads the frequencies of an input library and generates a 
// number line for selection by frequency during growth.
//JDB READ_THE_FRAGMENT_FREQUENCIES
void
DN_Build::read_frag_frequencies( string frequencies_file ){
    //open some strings for recordkeeping
    string filename, line, line_save_freq, line_save_frag, temp_line_save, freq_str;
    long double total_values = 0; // the running total for this line
    double freq; // the raw frequency
    size_t found; //last position of "-" character
    int line_number=1; // record the line number
    ifstream fragment_freq_file; //ifstream for the file
    fragment_freq_file.open(frequencies_file.c_str()); //freq file

    //Put in placeholders at the start for proper logic later
    ordered_fragment_frequencies.push_back(0.0);
    ordered_fragments.push_back("PLACEHOLDER");

    //Pulls all the lines from an input file
    while( getline (fragment_freq_file,line)){
        stringstream freq_convert;
        found = line.find_last_of("-");
        line_save_freq = line.substr(found+1);
        freq_convert << line_save_freq;
        freq_convert >> freq;

        line_save_frag = line.substr(0,found);
        //collects the 'running total' to generate bin right walls
        total_values = total_values + freq;
        
        //save things
        ordered_fragments.push_back(line_save_frag);
        ordered_fragment_frequencies.push_back(total_values);
        line_number++;
    }

    //uncomment below to get output for the entire number line
    /*
    for(int i=0; i<ordered_fragments.size(); i++){
        cout << "Fragment:" << ordered_fragments[i] << endl;
        cout << "   Frequency:" << ordered_fragment_frequencies[i] << endl;
    }
    cout << ordered_fragments.size() << endl;
    cout << ordered_fragment_frequencies.size() << endl;
    */ 
} // end DN_Build::read_frag_frequencies



//JDB
//Comparator functions for vector of vector sorting
//order_torsions sorts least to greatest by the second position of each vector
bool roulette_sort_order_torsions(const vector<double>& vect1, const vector<double>& vect2)
{
    if(vect1[1] > vect2[1])
    {
        return false;
    } else if(vect1[1] < vect2[1])
    {
        return true;
    }
    return false;
}
//order_frequencies sorts greatest to least by the first position of each vector
bool roulette_sort_order_frequencies(const vector<double>& vect1, const vector<double>& vect2)
{
    if(vect1[0] > vect2[0])
    {
        return true;
    } else if(vect1[0] < vect2[0])
    {
        return false;
    }
    return false;
}

bool fragment_sort_by_name(const Fragment & frag1, const Fragment & frag2){
    std::string name1 = frag1.mol.energy;
    std::string name2 = frag2.mol.energy;
    int frag1_pos, frag2_pos;
    std::string::const_iterator it1 = name1.begin(), it2 = name2.begin();
    std::stringstream ss;
    while ((*it1) != '.'){
        it1++;
    }

    it1++;
    while (it1 != name1.end()){
        ss << (*it1);
        it1++;
    }
    ss >> frag1_pos;
    ss.clear();

    while ((*it2) != '.'){
        it2++;
    }

    it2++;
    while (it2 != name2.end()){
        ss << (*it2);
        it2++;
    }
    ss >> frag2_pos;
    if (frag1_pos > frag2_pos) {
        return false;
    } else { 
        return true;
    }

    return true;
}
//JDB - function to generate a roulette wheel from the provided torsion table
void
DN_Build::generate_roulette(){
    //initialization of variables to extract info from the torsion table
    string filename, line, line_save_freq, line_save_tors, temp_line_save, freq_str;
    long double total_values = 0;
    double freq;
    size_t found;
    double line_number=1;
    ifstream roulette_file;
    
    //initialization of vectors to hold roulette information for later construction and reordering
    std::vector<double> freq_raw;
    std::vector<vector<double> > freq_ordering;
    roulette_file.open(dn_torenv_table.c_str());
    
    //extracts each line from the torsion table
    while( getline (roulette_file,line)){
        //saves the frequency string, pushes to double(freq)
        stringstream freq_convert;
        found = line.find_last_of("-");
        line_save_freq = line.substr(found+1);
        freq_convert << line_save_freq;
        freq_convert >> freq;
        
        //saves the torsion
        line_save_tors = line.substr(0,found);
        
        //gets a running total of each frequency
        total_values = total_values + freq;
        
        //saves the torsion to a vector
        torsion_vect.push_back(line_save_tors);
        /*saves the frequency and line number to a vector, then pushes that vector to a vector of vectors
          this is done for reordering later, as the torsion table is order-specific*/
        freq_raw.push_back(freq);
        freq_raw.push_back(line_number);
        freq_ordering.push_back( freq_raw );
        
        //clears frequency and prepares line number for next line
        freq_raw.clear();
        line_number++;
    }
    //sorts vector of vector greatest to least by its frequencies
    //this does not need to be frag_sorted, it a vector of vectors of doubles
    sort(freq_ordering.begin(), freq_ordering.end(), roulette_sort_order_frequencies);
    
    //initializes collection vectors for bookkeeping and frequency modifying
    std::vector<double> freq_percent;
    std::vector<int> duplicate_pos_collector;
    
    for(int i=0; i<freq_ordering.size();i++)
    {
        freq_percent.push_back(freq_ordering[i][0]/total_values);
    }
    
    //loops through the entire vector, up to the N-1 position
    for (int i=0; i<freq_percent.size(); i++){
        
        //checks at each step if loop is at the N-1 position - the and structure
        //will overload as false and not check the second half of the statement
        //if N-1 position is reached.
        if ((i != freq_percent.size()-1) && (freq_percent[i] == freq_percent[i+1])){
            //test for all equivalent positions up to N-2; loop stops at N-2
            for (int y=i; y<freq_percent.size()-1; y++){
                //checks to see if each position is N-2
                if (y == freq_percent.size()-2) {
                    //if y = N-2, check if last value in vector (N-1) is
                    //a duplicate
                    if (freq_percent[y] == freq_percent[y+1]) {
                        
                        //if duplicate, record current position and +1
                        duplicate_pos_collector.push_back(y);
                        duplicate_pos_collector.push_back(y+1);
                        
                    } else {
                        //if not duplicate, record current position only.
                        duplicate_pos_collector.push_back(y);
                        break;
                    }
                    
                    //if not at N-1, record each individual position as duplicate and continue.
                } else {
                    if (freq_percent[y] == freq_percent[y+1]){
                        duplicate_pos_collector.push_back(y);
                        
                        //if next not duplicate, record current and break loop entirely.
                    } else { //break out if not equal
                        duplicate_pos_collector.push_back(y);
                        break;
                    }
                }
                
                
            }
            
            //combine all frequencies that are duplicates into a single bin.
            double duplicate_frequency_collector;
            if(i==0){
                duplicate_frequency_collector = (double(duplicate_pos_collector.size()) * freq_percent[i]);
            } else {
                duplicate_frequency_collector = double(freq_percent[i-1]) + (double(duplicate_pos_collector.size()) * freq_percent[i]);
            }
            
            
            //loops through each position as indexed above and changes them all
            //to the same total frequency.
            for(int z=0; z<duplicate_pos_collector.size(); z++)
            {
                freq_percent[duplicate_pos_collector[z]] = duplicate_frequency_collector;
            }
            
            //sets the position of the outer loop to the last one in the duplicate bin
            //and clears the duplicate collector for the next set of duplicates
            i=duplicate_pos_collector.back();
            duplicate_pos_collector.clear();
        } else { // if next value is not a duplicate
            //if it's the very first position, do nothing
            if (i==0) {
                continue;
            }
            
            //otherwise, add previous value to the current position
            freq_percent[i] = freq_percent[i] + freq_percent[i-1];
        }
    }
    //places the ordered freq_percents, now in roulette form, in the original vector
    for(int i=0; i<freq_ordering.size(); i++)
    {
        freq_ordering[i][0] = freq_percent[i];
    }
    freq_percent.clear();
    //sorts the roulette by the order found in the torsion list
    //this does not need to be frag_sorted, it a vector of vectors of doubles
    sort(freq_ordering.begin(), freq_ordering.end(), roulette_sort_order_torsions);
    
    //
    for(int i=0; i<freq_ordering.size(); i++)
    {
        roulette_vect.push_back(freq_ordering[i][0]);
    }
    //initialize variable for the roulette bins- this will be used for
    //keeping the roulette in an ordered bin format for comparison later
    roulette_bin_vect.push_back(0);
    for(int i=0; i<roulette_vect.size(); i++)
    {
        double roulette_bin_val = roulette_vect[i];
        roulette_bin_vect.push_back(roulette_bin_val);
    }
    
    //sorts the vector to be smallest to largest
    //This doesn't need to be frag_sorted because it is vec of doubles
    sort(roulette_bin_vect.begin(), roulette_bin_vect.end());
    
    //iterates over roulette_bin_vect and collects the unqiue values, then resizes to only unique values
    
    vector<double>::iterator bin_it;
    bin_it = unique(roulette_bin_vect.begin(),roulette_bin_vect.end());
    roulette_bin_vect.resize(distance(roulette_bin_vect.begin(),bin_it));
    
    
    //VERBOSE OUTPUT
    if (Parameter_Reader::verbosity_level() > 1)
    {
        //Torsion with its associated roulette frequency
        for(int i=0; i<roulette_vect.size(); i++)
        {
            cout << torsion_vect[i] << " " << roulette_vect[i] << endl;
        }
        //the list of bins
        for(int i=0; i<roulette_bin_vect.size(); i++)
        {
            cout << roulette_bin_vect[i] << endl;
        }
    }
    cout << "Done (# Roulette = " << freq_ordering.size() << ")";
    cout << " (# Bins = " << roulette_bin_vect.size() <<")" << endl << endl;
} // end DN_Build::generate_roulette()

//JDB - used to read in a matrix of pairwise frequencies of fragments
//made during fragment library generation
void
DN_Build::read_frequencies_matrix( std::string freq_matrix_filename){
    //initalize some variables, open file
    std::string line, collect_string;
    ifstream freq_matrix_file;
    freq_matrix_file.open(freq_matrix_filename.c_str());

    //extracts first line of the file
    if (getline(freq_matrix_file,line)) {
        //ensures proper formatting (lnk-0,lnk-1)
        if (line.find(".") != std::string::npos){
            //for all characters in the string
            for (int i=0; i<line.size(); ++i){
                //append to vector if comma is reached
                if (line[i] == ','){
                    ordered_frag_list.push_back(collect_string);
                    collect_string.clear(); //clear for next loop
                    continue; //skip next collection step
                } //otherwise collect character to string
                collect_string+=line[i];
            }  
        } else {
            //line must have dashes to be properly formatted - ie lnk-0
            cout << "IMPROPER INPUT FORMAT FOR MATRIX" << endl;
        }

    }
    //prepare the 2D vectors that contains the entire frequency connection matrices (raw and relative)
    std::vector<std::vector<int> > frag_matrix(ordered_frag_list.size(), vector<int> (ordered_frag_list.size(), 0));
    std::vector<std::vector<float> > rel_frag_matrix(ordered_frag_list.size(), vector<float> (ordered_frag_list.size(), 0.0));
    

    //a couple integers to keep track of matrix position while in the loop
    int current_mat_position_v=-1; // -1 so it equals 0 at the start of the loop
    int current_mat_position_h=0;

    //loops through second+ line in the input matrix file
    while( getline (freq_matrix_file,line)){
        current_mat_position_v+=1;
        current_mat_position_h=0;
        for(int i =0; i<line.size(); ++i){
            if (line[i] == ','){
                frag_matrix[current_mat_position_v][current_mat_position_h] = stoi(collect_string);
                collect_string.clear();
                current_mat_position_h+=1;
                continue;
            }
            collect_string += line[i];
        }
    }

    for (int i=0; i<frag_matrix.size(); ++i){
        int row_sum=0;
        for (int x=0; x<frag_matrix[i].size(); ++x){
            row_sum+=frag_matrix[i][x];
        }

        for (int x=0; x<rel_frag_matrix[i].size(); ++x){
            if(row_sum > 0){
                rel_frag_matrix[i][x] = static_cast<float>(frag_matrix[i][x]) / row_sum;                
            }
        }
    }

    /*
    //Uncomment to print out the matrix that's read in
    for (int i=0; i<frag_matrix.size(); ++i){
        for(int x=0; x<frag_matrix[i].size(); ++x){
            cout << frag_matrix[i][x] << ",";
        }
        cout << endl << endl;
    }
    */
    
    /*
    //Uncomment to print out the generated relative values matrix
    for (int i=0; i<rel_frag_matrix.size(); ++i){
        for(int x=0; x<rel_frag_matrix[i].size();++x){
            cout << rel_frag_matrix[i][x] << ",";
        }
        cout << endl << endl;
    }
    */

    fragment_frequency_matrix = frag_matrix;
    rel_fragment_frequency_matrix = rel_frag_matrix;
    cout << " Done!" << endl;
    return;
} //end DN_Build::read_frequencies_matrix

// +++++++++++++++++++++++++++++++++++++++++
// Iterate through a vector of fragments, compute pair-wise Tanimotos, and remember rank ordered
// lists of similarity using the FragGraph data structure.
void
DN_Build::prepare_fragment_graph( vector <Fragment> & frag_vec, vector <FragGraph> & fraggraph_vec )
{

    // Declare a fingerprint object for computing Tanimotos
    Fingerprint finger;
    float tanimoto;

    // For everything in the fragment vector
    for (int i=0; i<frag_vec.size(); i++){

        // Declare some temporary variables to assemble the FragGraph
        FragGraph temp_fraggraph;
        pair <float, int> temp_pair;

        // First compute pair-wise Tanimotos
        for (int j=0; j<frag_vec.size(); j++){
            tanimoto = finger.compute_tanimoto_fragment(frag_vec[i].mol, frag_vec[j].mol);
            //cout << tanimoto <<endl;

            // Assemble the Tanimoto and index into a pair
            temp_pair.first = tanimoto;
            temp_pair.second = j;

            // Push that pair onto the temp_fraggraph
            temp_fraggraph.tanvec.push_back(temp_pair);
        }        

        // Sort the Tanimoto vector in the temp_fraggraph

        // This does not need to be frag_sorted. it is a vector a pairs of float and ints
        sort(temp_fraggraph.tanvec.begin(), temp_fraggraph.tanvec.end(), fgpair_sort);

        // Record the ranks in the rank vector
        // (start at 1 to avoid recording the self-self Tanimoto)
        for (int j=1; j<temp_fraggraph.tanvec.size(); j++){
            temp_fraggraph.rankvec.push_back( temp_fraggraph.tanvec[j].second );
        }

        // Add this object to the vector of FragGraph
        fraggraph_vec.push_back(temp_fraggraph);
    }

    return;

} // end DN_Build::prepare_fragment_graph()


// 
void
DN_Build::utilities_methods( AMBER_TYPER & typer ){

    if ( dn_ult_score.compare("volume")==0 ){
        char result[ 256];
        ssize_t count = readlink( "/proc/self/exe", result, 256);

        std::string EXE_PATH( result, (count > 0) ? count : 0 ) ;
        std::size_t botDirPos = EXE_PATH.find_last_of("/");
        std::string BINDIR    = EXE_PATH.substr(0, botDirPos); 
        botDirPos             = BINDIR.find_last_of("/");
        std::string BASEDIR   = BINDIR.substr(0, botDirPos); 

        std::ostringstream oss;
        oss << BASEDIR << "/parameters/chem.defn";
        std::string paramater_path = oss.str();
    
        CHEM_TYPER chem_typer;
        chem_typer.get_chem_labels(paramater_path);
 
        for (DOCKMol & mol: vec_utilities){
            typer.prepare_molecule( mol, true, false, false, false ); 
            chem_typer.apply_chem_labels( mol );
        }
    }

 
    float clust_tol = dn_ult_clust_tol;


    if ( dn_ult_method.compare("best_first") == 0 ) {

        // Sort by grid score
        sort(vec_utilities.begin(),vec_utilities.end(),[](DOCKMol a, DOCKMol b)
                              {
                                  return a.current_score < b.current_score;
                              }); 

        // Make directories called best first
        std::stringstream    base;
        base << "./best_first";
        string base_str = base.str();
        mkdir(base_str.c_str(), 0744);

        bool done = true;

        // declaring a vector with size 5
        vector<int> vec_index(vec_utilities.size());
 
        // initializing using iota() to make a vector index
        iota( vec_index.begin(), vec_index.end(), 0);

        vector<int> vec_index_after;  

        int clust_num = 1;
        int num_clusters = 0;


        // Start clustering
        while( done ){
 
            cout << "#######################" << endl;
            cout << "Cluster: " << clust_num << endl;
            cout << "index start size: " << vec_index.size() << endl;

            num_clusters++; 

            DOCKMol refmol;
 
            refmol = vec_utilities[vec_index[0]];

            vector<DOCKMol> cluster_vec;
            Volume_Score_Comp final_comp_perfect;
            final_comp_perfect = final_comp_perfect.compute_score_analytical( refmol, refmol );

            float hms_score_perfect  = dn_ult_hun_matching_coeff * 1 + 
                                       dn_ult_hun_rmsd_coeff * 0.0;

            for (const int & i :vec_index) {
                vector <string> temp_tokens_ref {};
                vector <string> temp_tokens_test {};
                if ( dn_ult_bf_anc_diverse ) {
  
                    Tokenizer(refmol.energy, temp_tokens_ref, *"--");
                    Tokenizer(vec_utilities[i].energy, temp_tokens_test, *"--");
                }

                /* 
                    HUNGARIAN BEST FIRST SCORE
                */

                if ( dn_ult_score.compare("hungarian")==0 ){
                  
                    Hungarian_RMSD h;
                    int reffrag_heavy_atoms = 0;

                    for (int j=0; j<refmol.num_atoms; j++){
                        if (refmol.atom_types[j].compare("H") != 0 
                            && refmol.atom_types[j].compare("Du") != 0){
                            reffrag_heavy_atoms++;
                        }
                    }
            
                    if (reffrag_heavy_atoms == 0){
                        cout << "EXIT: some of molecules have no heavy atoms" << endl;
                        exit(2);
                    }
 
                    std::pair <float, int> h_result = 
                                 h.calc_Hungarian_RMSD_dissimilar( refmol, vec_utilities[i] );

                    float hms_score  = dn_ult_hun_matching_coeff * 
                                       (reffrag_heavy_atoms - h_result.second) / 
                                       reffrag_heavy_atoms + 
                                       dn_ult_hun_rmsd_coeff * (h_result.first);


                    if (dn_ult_bf_anc_diverse){

                        if ( (hms_score_perfect == hms_score ) || ( hms_score < clust_tol  && 
                               temp_tokens_ref[0] != temp_tokens_test[0] )) {
                            vec_utilities[i].score_hun = hms_score;
                            cluster_vec.push_back(vec_utilities[i]);

                        } else if( hms_score < clust_tol && 
                                   temp_tokens_ref[0] == temp_tokens_test[0] ) {
                            continue;
                        } else {
                            vec_index_after.push_back(i);
                        }
                    } else {
                       if ( hms_score < clust_tol ){
                            vec_utilities[i].score_hun = hms_score;
                            cluster_vec.push_back(vec_utilities[i]);
                       } else {
                            vec_index_after.push_back(i);
                       }
                    } 
                } 


                /* 
                    VOLUME BEST FIRST SCORE
                */

                if ( dn_ult_score.compare("volume")==0 ){

                    Volume_Score_Comp final_comp;
                    
                    final_comp = final_comp.compute_score_analytical( refmol, vec_utilities[i] );
 
                    float vol_score = -1; 
                    float vol_score_perfect = -1; 
                    if ( dn_ult_vol_type.compare("volume_heavy")==0) {
                        vol_score = final_comp.heavy_atom_component;
                        vol_score_perfect = final_comp_perfect.heavy_atom_component;
                    }

                    if (vol_score == -1) { 
                        cout << "A Volume score is negative. Not possible" << std::endl;
                        exit(2);
                    }
     
                    if (dn_ult_bf_anc_diverse){
                        if ( ( vol_score_perfect == vol_score )|| ( vol_score > clust_tol  && 
                               temp_tokens_ref[0] != temp_tokens_test[0] ) ){
                            vec_utilities[i].score_vol = vol_score;
                            cluster_vec.push_back( vec_utilities[i] );

                        } else if ( vol_score > clust_tol  && 
                                   temp_tokens_ref[0] == temp_tokens_test[0] ) {
                            continue;
                        } else {
                            vec_index_after.push_back( i );
                        }
                    } else {
                       if ( vol_score > clust_tol  ){
                            vec_utilities[i].score_vol = vol_score;
                            cluster_vec.push_back(vec_utilities[i]);
                       } else {
                            vec_index_after.push_back(i);
                       }
                    } 
                } 

            }
         
            if ( cluster_vec.size() > dn_ult_clust_size_lim) {

                ostringstream fout_cluster_name;

                fout_cluster_name << base_str << "/cluster_" << clust_num << ".mol2";
                fstream fout_cluster;
                fout_cluster.open( fout_cluster_name.str().c_str(), 
                                   fstream::out|fstream::app);
                fout_cluster_name.clear();
                // Sort and Write mol2 for each cluster
                if ( dn_ult_score.compare("hungarian")==0 ){
                    sort(cluster_vec.begin(),cluster_vec.end(),[](DOCKMol a, DOCKMol b)
                                          {
                                              return a.score_hun < b.score_hun;
                                          }); 
                }

                if ( dn_ult_score.compare("volume")==0 ){

                    sort(cluster_vec.begin(),cluster_vec.end(),[](DOCKMol a, DOCKMol b)
                                          {
                                              return a.score_vol >  b.score_vol;
                                          }); 
                }
                int rank = 0;

                for (DOCKMol & mol: cluster_vec){


                    fout_cluster << DELIMITER << setw(STRING_WIDTH)
                                   << "Grid_Score:" << setw(FLOAT_WIDTH) 
                                   << mol.current_score<<endl;

                    fout_cluster << DELIMITER << setw(STRING_WIDTH)
                                   << "Grid_vdw_energy:" << setw(FLOAT_WIDTH) 
                                   << mol.vdw_comp <<endl;

                    fout_cluster << DELIMITER << setw(STRING_WIDTH)
                                   << "Grid_es_energy:" << setw(FLOAT_WIDTH) 
                                   << mol.es_comp <<endl;

                    fout_cluster << DELIMITER << setw(STRING_WIDTH) 
                                   << "Internal_energy_repulsive:" << setw(FLOAT_WIDTH) 
                                   << mol.internal_energy << endl;

                    fout_cluster << DELIMITER << setw(STRING_WIDTH)
                                   << "Name:" << setw(FLOAT_WIDTH) 
                                   << mol.title <<endl;

                    fout_cluster << DELIMITER << setw(STRING_WIDTH)
                                   << "Molecular_Weight:" << setw(FLOAT_WIDTH) 
                                   << mol.mol_wt<<endl;

                    fout_cluster << DELIMITER << setw(STRING_WIDTH)
                                   << "DOCK_Rotatable_Bonds:" << setw(FLOAT_WIDTH) 
                                   << mol.rot_bonds<<endl;

                    fout_cluster << fixed << DELIMITER << setw(STRING_WIDTH)
                                   << "Formal_Charge:" << setw(FLOAT_WIDTH) 
                                   << mol.formal_charge<<endl;

                    fout_cluster << DELIMITER << setw(STRING_WIDTH)
                                   << "HBond_Acceptors:" << setw(FLOAT_WIDTH) 
                                   << mol.hb_acceptors<<endl;

                    fout_cluster << DELIMITER << setw(STRING_WIDTH)
                                   << "HBond_Donors:" << setw(FLOAT_WIDTH) 
                                   << mol.hb_donors<<endl;

                    fout_cluster << DELIMITER << setw(STRING_WIDTH)
                                   << "Frag_String:" << setw(FLOAT_WIDTH) 
                                   << mol.energy<<endl;
                    
                    fout_cluster << DELIMITER << setw(STRING_WIDTH)
                                   << "Cluster:" << setw(FLOAT_WIDTH) 
                                   << clust_num <<endl;

                    fout_cluster << DELIMITER << setw(STRING_WIDTH)
                                   << "Rank:" << setw(FLOAT_WIDTH) 
                                   << rank <<endl;

                    if ( dn_ult_score.compare("hungarian")==0){
                        fout_cluster << DELIMITER << setw(STRING_WIDTH)
                                       << "Cluster_Score:" << setw(FLOAT_WIDTH) 
                                       << mol.score_hun <<endl;
                    } 
              
                    if ( dn_ult_score.compare("volume")==0){
                        fout_cluster << DELIMITER << setw(STRING_WIDTH)
                                       << "Cluster_Score:" << setw(FLOAT_WIDTH) 
                                       << mol.score_vol <<endl;
                    }
                    Write_Mol2(mol, fout_cluster);
                    rank ++;
                }
            }

            cout << "index end size: " << vec_index_after.size() << endl;
            cout << "num of mols: " 
                 << vec_index.size() - vec_index_after.size() << endl;
            cout << "num of write out mols: "  << cluster_vec.size() << std::endl;
 
            cout << "#######################" << endl;
            refmol.clear_molecule();
            clust_num++;


            if ( vec_index_after.size() == 0 || 
                 num_clusters >= dn_ult_num_clusters ){
                done = false;
            }

            // update the vec of index after collecting non tol checks
            vec_index = vec_index_after;
            vec_index_after.clear();
            cluster_vec.clear();
        }
     

    }
    return;
}


// +++++++++++++++++++++++++++++++++++++++++
// This function is the main de novo engine for building molecules. (called in dock.cpp)
void
DN_Build::build_molecules( Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & typer,
                           Orient & orient )
{   
    Trace trace("DN_Build::build_molecules()");

    if (dn_mode.compare("utilities") == 0){
       
        cout << endl << endl <<"##### Entering the main de novo utiliies engine #####" << endl << endl;



        utilities_methods( typer); 
    }


    cout << endl << endl <<"##### Entering the main de novo build engine #####" << endl << endl;
    

    Iso_Table::Iso_Tab iso_table;
    
    if ( dn_sampling_method_isoswap ){ //START iso_table preparation
        char result[ 256];
        ssize_t count = readlink( "/proc/self/exe", result, 256);

        std::string EXE_PATH( result, (count > 0) ? count : 0 ) ;
        std::size_t botDirPos = EXE_PATH.find_last_of("/");
        std::string BINDIR    = EXE_PATH.substr(0, botDirPos); 
        botDirPos             = BINDIR.find_last_of("/");
        std::string BASEDIR   = BINDIR.substr(0, botDirPos); 

        std::ostringstream oss;
        oss << BASEDIR << "/parameters/chem.defn";
        std::string paramater_path = oss.str();
       
        CHEM_TYPER chem_typer;
        chem_typer.get_chem_labels(paramater_path);
        std::cout << "amber typing isofragment input library..."<< std::endl;
        for(unsigned int i = 0; i < isosidechains.size(); i++){
           typer.prepare_molecule(isosidechains[i].mol, true, false, false, false); 
           chem_typer.apply_chem_labels(isosidechains[i].mol);
        }
        for(unsigned int i = 0; i < isolinkers.size(); i++){
           typer.prepare_molecule(isolinkers[i].mol, true, false, false, false); 
           chem_typer.apply_chem_labels(isolinkers[i].mol);
        }
        for(unsigned int i = 0; i < isoscaffolds.size(); i++){
           typer.prepare_molecule(isoscaffolds[i].mol, true, false, false, false); 
           chem_typer.apply_chem_labels(isoscaffolds[i].mol);
        }
        for (unsigned int i =0;i<scaf_link_sid.size(); i++){
           typer.prepare_molecule(scaf_link_sid[i].mol, true, false, false, false);
           chem_typer.apply_chem_labels(scaf_link_sid[i].mol);
        }
        //This is where we write stuff
        if (iso_parm.get_write_libraries()){
            cout <<endl <<endl <<"##### Entering writing the Iso_Table #####" <<endl <<endl;
            std::cout << "aligning fragments..." << std::endl;

            int progress = 0;
            for (unsigned int i =0;i<scaf_link_sid.size(); i++){
                std::pair<std::vector<Iso_Acessory::Scored_Triangle>, std::vector<Fragment>> vec_iso_frags {};
                vec_iso_frags =  frag_iso_align(scaf_link_sid[i]);
                iso_table.set(scaf_link_sid[i],vec_iso_frags.second);
                vec_iso_frags.first.clear();
                vec_iso_frags.second.clear();


                 
                std::cout << "\r" << "[ "; 
                float percentage = (float)(i+1)/(float)scaf_link_sid.size();
                int val = (int) (percentage * 100);
                int lpad = (int) (percentage * 50);
                
                for (int z=0; z<lpad; z++){
                    std::cout << "|" ;
                }
                std::cout << " "<<std::fixed<<std::setprecision(2)<<(float)(percentage * 100) << "% ]";

                if (i == scaf_link_sid.size()-1){
                    std::cout << "\nFinished setting isotable"<<std::endl; 
                    std::cout << "Iso_Table All Heads Size: " << iso_table.get_size() << std::endl;
                    std::cout << "Iso_Table All Tails Size: " << iso_table.get_size_total() << std::endl;
                }

                std::cout.flush();
            }

            
            std::cout << "Writing mol2s" << std::endl;
            iso_table.write_table(dn_iso_output_path_libraries,iso_parm.get_iso_num_top());
            std::cout << "Finished on writing molecules" << std::endl;
        } else if (iso_parm.get_iso_fraglib().first){
            std::cout << "Reading in iso tables\n Loading..."<< std::endl;
            std::vector<std::string> filenames = iso_parm.get_iso_fraglib().second;
            //This does not need to be fragsorted because it is vec of strings
            std::sort(filenames.begin(),filenames.end(),[](std::string a, std::string b){

                size_t pos_a = a.find( ".mol2" );
                if ( pos_a != string::npos ) {
                   a.replace( pos_a, 4, "" );
                }
                size_t pos_b = b.find( ".mol2" );
                if ( pos_b != string::npos ) {
                   b.replace( pos_b, 4, "" );
                }
                size_t last_index_a = a.find_last_not_of("0123456789");
                int result_a = atoi(a.substr(last_index_a + 1).c_str());

                size_t last_index_b = b.find_last_not_of("0123456789");
                int result_b = atoi(b.substr(last_index_b + 1).c_str());
                return result_a<result_b;
            });

            int progress = 0;
            std::cout << "[";
            for (unsigned int i = 0; i < filenames.size(); i++){ 
                if (filenames[i] == "." || filenames[i] == ".."){
                    continue;
                }
	            //cout <<" Reading library from " <<filenames[i]<<"..."<<std::endl;
                std::stringstream path;
                std::vector<Fragment> tmp_vec {};
                path <<dn_iso_fraglib_dir << "/"<<filenames[i];
    	        read_frag_library(tmp_vec, path.str());


                //This is where we prepare the read in iso tail molecules
                for (int j=0; j<tmp_vec.size(); j++) {
                    activate_mol(tmp_vec[j].mol);
                    Fingerprint temp_finger;
                    for (int z=0; z<tmp_vec[j].mol.num_atoms; z++){
                         tmp_vec[j].mol.atom_envs[z] = temp_finger.return_environment(tmp_vec[j].mol, z);
                    }
                    tmp_vec[j].mol.prepare_molecule();
                    // Compute the charges, saving them on the mol object
                    float total_charges = 0;
                    total_charges = compute_gast_charges(tmp_vec[j].mol);
                    calc_mol_wt(tmp_vec[j].mol);
                    calc_formal_charge(tmp_vec[j].mol);
                    calc_rot_bonds(tmp_vec[j].mol);
                    tmp_vec[j].mol.rot_bonds = 0;

                }
            
                ifstream fin_tmp;
                fin_tmp.open(path.str());
                if (fin_tmp.fail()){
                    cout <<"Could not open " <<path.str()<<". Program will terminate." <<endl;
                    cout << strerror(errno);
                    exit(0);
                }    
                std::string line; 
                std::string string_head;

                while(getline(fin_tmp,line)){
                    if (line.find("#HEAD") != string::npos){
                        istringstream ss(line);
                        std::vector<std::string> head_vec{};
                        std::string word;
                        while(ss >> word){
                            head_vec.push_back(word);
                        }
                        if (head_vec.size() == 2){ 
                            string_head = head_vec[1];
                            string_head.insert(0,"iso_");
                        }else{
                            std::cout << "You are using isoswap incompatible files. exiting..." << std::endl;
                            exit(1); 
                        }
                    }
                }
                Fragment rep_head_frag; 
                for (unsigned int i=0; i<tmp_vec.size();i++){
                    if (tmp_vec[i].mol.energy == string_head){
                        rep_head_frag = tmp_vec[i];
                        break;
                    }
                    if (i==tmp_vec.size()-1){
                        std::cout << "Your isoswap files don't have the head molecular structure. exiting..." << std::endl;
                        exit(1); 
                    }
                }

                iso_table.set(rep_head_frag,tmp_vec);

                std::cout << "\r" << "[ "; 
                float percentage = (float)(i+1)/(float)filenames.size();
                int val = (int) (percentage * 100);
                int lpad = (int) (percentage * 50);
                
                for (int z=0; z<lpad; z++){
                    std::cout << "|" ;
                }
                std::cout << " "<<std::fixed<<std::setprecision(2)<<(float)(percentage * 100) << "% ]";

                if (i == filenames.size()-1){
                    std::cout << "\nFinished reading and loading up isotable"<<std::endl; 
                    std::cout << "Iso_Table All Heads Size: " << iso_table.get_size() << std::endl;
                    std::cout << "Iso_Table All Tails Size: " << iso_table.get_size_total() << std::endl;
                }

                std::cout.flush();
            }
    
        }
    }//END preparation of isotable

    molecule_counter = 0;    // this is an estimate of how many molecules are sampled
    growth_tree_index = 0;   // in case we are writing growth trees
    num_all_frags_layer = 0;
    num_prune_frags_layer = 0;

    if (verbose){ cout <<"verbose is turned on" <<endl; }

    if (orient.orient_ligand == false && dn_user_specified_anchor == false){
        cout << endl
             << "# Note: No anchors were specified AND orient is turned off." << endl
             << "# DOCK now assumes that all scaffolds, linkers, and sidechains" << endl
             << "# provided are already oriented to the binding site." << endl << endl;
    }

    // fout_molecules is the filestream for complete molecules
    ostringstream fout_molecules_name;
    fout_molecules_name << dn_output_prefix << ".denovo_build.mol2";
    fstream fout_molecules;

    //LEP - 2018.06.20 
    fout_molecules.open(fout_molecules_name.str().c_str(), fstream::out|fstream::app);
    fout_molecules_name.clear();

    #ifdef BUILD_DOCK_WITH_RDKIT
    // fout_rejected is the filestream for complete molecules
    ostringstream fout_rejected_name;
    fout_rejected_name << dn_output_prefix << ".denovo_rejected.mol2";
    fstream fout_rejected;
    if (dn_save_all_molecules) {
        fout_rejected.open(fout_rejected_name.str().c_str(), fstream::out|fstream::app);
    }    
    fout_rejected_name.clear();
    #endif

    // Declare some fragment vectors for growing
    // UTILITY VECTORS
    vector <Fragment> root;          // starting fragments, returned each layer
    vector <Fragment> layer;         // vector for completing a layer
    vector <Fragment> next_layer;    // to be moved to the next layer
    vector <Fragment> growing;       // temporary vector for torsions / minimizing

    // PRUNING VECTORS
    // IF NOT WRITING OUT, USE A COUNTER TO DISPLAY NUMBER
    vector <Fragment> prune_rmsd_mw;    // not as bad as it sounds
    prune_rmsd_mw_counter_prune = 0;
    prune_root_counter_prune = 0;
    cand_root_ign_counter_prune = 0;
    filtered_comp_counter_prune = 0;

    failed_cap_H_prune_rmsd_mw = 0;
    failed_cap_H_prune_root = 0;
    failed_cap_H_cand_root_ign = 0;
    failed_cap_H_filtered_comp = 0;
    failed_cap_H_duplicate_dump = 0;

    failed_vtm_prune_rmsd_mw = 0;
    failed_vtm_prune_root = 0;
    failed_vtm_cand_root_ign = 0;
    failed_vtm_filtered_comp = 0;
    failed_vtm_duplicate_dump = 0;

    vector <Fragment> duplicate_dump;
    duplicate_dump_counter_prune = 0;

    vector <int> available_fragment_indices;
    vector <int> selected_fragments;

    #ifdef BUILD_DOCK_WITH_RDKIT
    vector <Fragment> rejected;      // rejected molecules
    #endif


    // If orient_ligand is yes, then either choose up to dn_unique_anchors fragments from the
    // anchor file, or if the anchor file is not provided, choose them from the other fraglibs.
    // Then, one at a time, generate up to max_orientations orients for an anchor, put them all in
    // root, and go through de novo growth. Repeat the process for all other anchors. Each starting
    // orient should have its own set of checkpoint files and orient files, but they can share the
    // final output file.
    //
    // If orient_ligand is no, then do the exact same thing above except do not enter the
    // orient_fragments function. Thus, each new instance of de novo growth will start with exactly
    // one pose of one anchor in root.


    // Pick the anchors (either for orienting or fixed-anchor de novo) here
    if (dn_user_specified_anchor == false){

        // TODO anchors is similar to scaf_link_sid could probably save effort 

        // for every fragment in all three libs, push them onto anchors
        if (scaffolds.size() != 0){
            for (int i=0; i<scaffolds.size(); i++){
                anchors.push_back(scaffolds[i]);
            }
        }
        if (linkers.size() != 0){
            for (int i=0; i<linkers.size(); i++){
                anchors.push_back(linkers[i]);
            }
        }
        if (sidechains.size() != 0){
            for (int i=0; i<sidechains.size(); i++){
                anchors.push_back(sidechains[i]);
            }
        }
    }

    // Prepare each fragment now to simplify future things
    for (int i=0; i<anchors.size(); i++){
        anchors[i].mol.prepare_molecule();
        stringstream anchor_info;
        anchor_info <<DELIMITER<<setw(STRING_WIDTH)<<"Anchor:"<<setw(FLOAT_WIDTH)<<i <<"\n";
        anchors[i].mol.current_data = anchor_info.str();
        anchor_info.clear();

        //trace.boolean( "DN:score.use_chem", score.use_chem);
        trace.boolean( "DN:score.use_ph4", score.use_ph4);
        trace.boolean( "DN:score.use_volume", score.use_volume);
        typer.prepare_molecule(anchors[i].mol, true, score.use_chem, score.use_ph4, score.use_volume);
    }

    // This will check if the dummy atoms have vdw param = 0 and will exit if so.
    if (typer.atom_typer.atom_types.size() > 0) {
        if ( (typer.atom_typer.types[typer.atom_typer.atom_types[0]].radius == 0.0) || (typer.atom_typer.types[typer.atom_typer.atom_types[0]].well_depth == 0.0) ){
            cout <<"Dummy atoms must have well depth > 0 and radius > 0 . Use vdw_de_novo.defn in parameters folder. Program will terminate." <<endl;
            exit(0);
        } 
    }
   // ADDITION: CS preparing descriptors, which are called in sampling classes
    if (dn_user_specified_anchor == true){
       // For each use specified anchor
       for (int i=0; i<anchors.size(); i++){
           activate_mol(anchors[i].mol);
           calc_mol_wt(anchors[i].mol);
           calc_rot_bonds(anchors[i].mol);
           calc_formal_charge(anchors[i].mol);
           //calc_num_HA_HD(anchors[i].mol);
       }
    }

    // Sort anchors by # of heavy atoms, in case of a tie, by # aps
    frag_sort(anchors, size_sort);

    // Resize dn_unique_anchors in case the user requests more anchors than are provided
    if (dn_unique_anchors > anchors.size()){
        dn_unique_anchors = anchors.size();
    }

    // Resize anchors to keep only dn_unique_anchors, and warn the user if they provided more
    // anchors than we are keeping here
    if (dn_user_specified_anchor && (dn_unique_anchors < anchors.size())){
        cout <<endl
             <<"# Note: " <<anchors.size() <<" anchors were provided in the anchor file, but the" <<endl
             <<"# parameter 'dn_unique_anchors' is set to " <<dn_unique_anchors <<". Thus, de novo" <<endl
             <<"# growth will only be seeded by " <<dn_unique_anchors <<" anchors." <<endl <<endl;
    }
    anchors.erase(anchors.begin()+dn_unique_anchors, anchors.end());


    ////////////////////////////////////////////////////////////////////////////
    // Begin Main Loop                                                        //
    ////////////////////////////////////////////////////////////////////////////


    // Initiate a new round of de novo growth for each anchor
    for (int a=0; a<anchors.size(); a++){


        // Reset some values here that may have been scaled for the previous anchor
        if (dn_sampling_method_graph){
           dn_graph_temperature = dn_temp_begin;
        }

        dn_pruning_conformer_score_cutoff = dn_pruning_conformer_score_cutoff_begin;
        
        // Reset the random seed so that it is the same for every anchor
        // This behavior is quasi-consistent with A&G
        simplex.initialize();

        // Clear the non-bonded pairlist() from the previous anchor
        score.primary_score->nb_int.clear();
        trace.note("nb array clear");
    
        // If orient_ligand is set to yes
        if (orient.orient_ligand){
    
            // Orient one anchor and return the results on the vector root
            cout <<"##### Now orienting anchor #" <<(a+1) <<endl;

            if (verbose){
                // This outputs what anchor is being chosen CPC
                cout<<"##### Name of anchor chosen is " <<anchors[a].mol.energy<<endl;
            }    

            orient_fragments(root, anchors[a], score, simplex, typer, orient);
            //cout << "Orient Anchor Charge Debug" <<  
        } else {

            // Just put one anchor onto root - it should already be in the binding site
            cout <<"##### Now growing from anchor #" <<(a+1) <<endl;
            root.push_back(anchors[a]);
        }
    
        // counter for controlling the next while loop
        int counter = 0;
    
        // last_layer feature to help control sampling
        bool last_layer = false;
        int growing_size_tracker = 0;
    
        // while there are still molecules in root && for as many layers as you want to grow out
        // (each iteration of this while loop is equivalent to growing out one layer)
        while (counter < dn_max_grow_layers){

            //Start the clock now for per layer
            std::chrono::time_point<std::chrono::system_clock> start_TIME, end_TIME;
            start_TIME = std::chrono::system_clock::now();
   
            num_all_frags_layer = 0;
            num_prune_frags_layer = 0;
            num_att_layer = 0;
            num_att_root = 0;

            if (counter == (dn_max_grow_layers-1)) {last_layer = true;}
            if (root.size() == 0){ cout <<"Root is empty, continuing." <<endl <<endl; break;} 
            cout <<"  ##### Entering layer of growth #" <<(counter+1) <<endl;
            cout <<"  Root size at the beginning of this layer = " <<root.size() <<endl;
            vector< Fragment > completed_molecules; // JDB JDB JDB
            //PAK
            #ifdef BUILD_DOCK_WITH_RDKIT
            int total_rej_molecules = 0;
            #endif
            // for each fragment in root
            for (int i=0; i<root.size(); i++){
                std::chrono::time_point<std::chrono::system_clock> root_TIME_START, root_TIME_END ;
                root_TIME_START = std::chrono::system_clock::now();

                num_pruned_roots = 0;
                num_filtered_comp = 0;
                num_completed = 0;
                int root_mol_counter = 0;
                int root_frag_counter = 0;
                num_att_root = 0;

                // In case the fragment has no attachment points and it is the beginning of growth,
                // prepare it now just so there will actually be a score
                if (counter == 0 && root[i].aps.size() == 0){ 
                    activate_mol(root[i].mol);
                    root[i].mol.prepare_molecule();

                    typer.prepare_molecule(root[i].mol, true, score.use_chem, score.use_ph4, score.use_volume);
                    score.compute_primary_score(root[i].mol);
                    root[i].scaffolds_this_layer = 0;
                }
    
                // Set scaffolds this layer to 0 before we start growth
                root[i].scaffolds_this_layer = 0;
              
                // Copy the current frag from root into layer
                layer.push_back(root[i]);

                // Count the number of attachment points in layer[0]
                int num_root_att = root[i].aps.size();
                // For each of those attachment points
                // (do the last attachment point first, then increment down. this will have the effect
                // of filling every attachment point on the fragment before going to the next layer)
                for (int j=num_root_att-1; j>-1; j--){
                    // For all of the fragments in layer
                    for (int k=0; k<layer.size(); k++){
                        // Copy on the old dockmol for the growth tree
                        if (dn_write_growth_trees){
                            layer[k].frag_growth_tree.push_back(layer[k].mol);
                        }
//TODO
// in some of these     different sampling methods, add additional print statements that are activated by a verbose flag on the command line
                        // remember how many things are on growing before sampling
                        int growing_size_before = growing.size();
    
                        // Sampling method = exhaustive
                        if (dn_sampling_method_ex){
                            sample_fraglib_exhaustive( layer[k], j, scaf_link_sid, growing, 
                                                       score, simplex, typer, last_layer );
    
                        // Sampling method = random
                        } else if (dn_sampling_method_rand){
                            //select_frags_from_fraglib(scaf_link_sid, available_fragment_indices, selected_fragments,
                            //                          dn_num_random_picks, last_layer);
                            //attach_selected_frags( layer[k], j, scaf_link_sid, scaf_link_sid_graph,
                            //                      growing, score, simplex, typer, selected_fragments, last_layer );
                            sample_fraglib_rand( layer[k], j, scaf_link_sid, growing,
                                                 score, simplex, typer, available_fragment_indices,
                                                 selected_fragments, last_layer );
    
                        // Sampling method = graph
                        } else if (dn_sampling_method_graph){
                            sample_fraglib_graph( layer[k], j, scaf_link_sid, scaf_link_sid_graph,
                                                  growing, score, simplex, typer, available_fragment_indices, 
                                                  selected_fragments, last_layer  );
                        // Sampling method = matrix
                        } else if (dn_sampling_method_matrix){
                            sample_fraglib_matrix ( layer[k], j, scaf_link_sid, growing,
                                                    score, simplex, typer, last_layer );
                        // Sampling method = iso_ex
                        } else if (dn_sampling_method_isoswap){
                            sample_isofraglib_isoswap( layer[k], j, scaf_link_sid, growing,
                                                          score, simplex, typer, last_layer, iso_table,
                                                          counter );
                        } else {
                            cout <<"You chose...poorly" <<endl;
                            exit(0);
                        }
                    }
//TODO
// dummy_to_H(layer[k], layer[k].aps[j].heavy_atom, layer[k].aps[j].dummy_atom);  growing.push_back(layer[k]);

                    // clear layer and get ready to copy growing onto layer
                    layer.clear();

                    frag_sort(growing,fragment_sort);

                    if ( !dn_legacy_prune_hrmsd ) {
                        prune_h_rmsd_and_mw(growing);
                    }

                    root_mol_counter += growing.size();
                    num_all_frags_layer += growing.size();
//TODO these descriptors are now computed in sample_minimized torsions and filtering is done then
//     doing it again here may be redundant if that information is copied when you are creating
//     the torsions
                    // copy at most dn_max_layer_size back onto layer - at this point the loop will
                    // continue adding fragments to the remaining available attachment points to fill
                    // a layer CPC

                    // If soft cutoff is used, the variable below matters. I would have created another user-defined DN_Build attribute for this variable,
                    // I didn't want to do it right away due to time constraints. A grad student should do it to learn proper OOP. [GDRM]
                    bool evaluate_both{false};

                    //JDB JDB JDB
                    //Collect all complete molecules and remove them from the growing set.
                    if (j == 0) {
                        //vector< Fragment >::iterator pos = growing.begin();
                        //while (pos != growing.end()) {
                        //    if (dummy_in_mol((*pos).mol)){
                        //        pos++;
                        //    } else {
                        //        completed_molecules.push_back((*pos));
                        //        pos = growing.erase(pos);
                        //    }
                        //}
                        // TO-DO: Report complete
                        std::list<Fragment> tmp_linked_list(growing.begin(),growing.end());
                        for (std::list<Fragment>::iterator frag = tmp_linked_list.begin(); frag != tmp_linked_list.end();){
                            if ( dummy_in_mol((*frag).mol) ) {
                                ++frag;
                            } else {
                                completed_molecules.push_back(*frag);
                                num_completed +=1;
                                frag = tmp_linked_list.erase(frag);
                            }
                        }
                        growing.clear();
                        for (Fragment tmp_frag : tmp_linked_list){
                            growing.push_back(tmp_frag);
                        }
                        ////////PAK LIST/////

                        //JDB - Assess all complete molecules separately to incomplete ones.
                        //Keeps the # collected separate, leading to more potential roots.
                        vector<Fragment> tmp_completed;
                        for(int x=0; x < completed_molecules.size(); x++) {

                            calc_mol_wt(completed_molecules[x].mol);
                            calc_rot_bonds(completed_molecules[x].mol);
                            calc_formal_charge(completed_molecules[x].mol);
                            calc_num_HA_HD(completed_molecules[x].mol);

                            if (completed_molecules[x].mol.rot_bonds <= dn_constraint_rot_bon && 
                                fabs(completed_molecules[x].mol.formal_charge) <= (fabs(dn_constraint_formal_charge)+0.1) ){
                                // If hard cutoff is used, evaluate only upper limit
                                //cout << completed_molecules[x].mol.mol_wt << endl;

                                if (dn_MW_cutoff_type_hard && completed_molecules[x].mol.mol_wt <= dn_upper_constraint_mol_wt){
                                    tmp_completed.push_back( completed_molecules[x] );
                                
                                // If soft cutoff is used, evaluate. 
                                } else if (!dn_MW_cutoff_type_hard && mw_cutoff(completed_molecules[x],1)){
                                    tmp_completed.push_back( completed_molecules[x] );

                                // In any other situation, the molecule fails the cutoff requirements
                                } else {
                                    num_filtered_comp += 1;
                                    filtered_comp_counter_prune +=1;
                                    write_mol_to_file( dn_output_prefix, "filtered_comp",
                                                       completed_molecules[x], a, counter);
                                    //filtered_comp.push_back( completed_molecules[x] );
                                }
                                continue;
                            }

                            num_filtered_comp += 1;
                            filtered_comp_counter_prune +=1;
                            write_mol_to_file( dn_output_prefix, "filtered_comp",
                                               completed_molecules[x], a, counter);
                            
                        }
                        completed_molecules.clear();
                        for(int x=0; x < tmp_completed.size(); x++) {
                            completed_molecules.push_back( tmp_completed[x] );
                        }
                        tmp_completed.clear();

                        //JDB JDB JDB
                        //print out the entire contents of new_layer
                        //bickel_write_out(growing, "allmols", counter, false);
                    }
                    

                    int cur_att_layer_size =0;
                    //Assess constructed mols based on cutoffs.
                    for (int x=0; x < growing.size(); ++x){

                        // Calculate descriptors
                        calc_mol_wt(growing[x].mol);
                        calc_rot_bonds(growing[x].mol);
                        calc_formal_charge(growing[x].mol);
                        calc_num_HA_HD(growing[x].mol);
                        
                        //if (cur_att_layer_size < dn_max_layer_size)    {
                        if ( layer.size() < dn_max_layer_size)    {
                            // Only copy things that are within the user-defined property constraints
                            if (growing[x].mol.rot_bonds <= dn_constraint_rot_bon && 
                                fabs(growing[x].mol.formal_charge) <= (fabs(dn_constraint_formal_charge)+0.1) ){

                                // If hard cutoff is used, evaluate upper limit and pass if it is below
                                if (dn_MW_cutoff_type_hard && growing[x].mol.mol_wt <= dn_upper_constraint_mol_wt){
                                    layer.push_back( growing[x] );
                                    root_frag_counter+=1;

                                // If soft cutoff is used, evaluate. If molecule below limit, mw-cutoff returns true, otherwise, returns true with a probability
                                } else if (!dn_MW_cutoff_type_hard && mw_cutoff(growing[x],1)){
                                    layer.push_back( growing[x] );
                                    root_frag_counter+=1;

                                // In any other situation, the molecule fails the cutoff requirements
                                } else {
                                    num_pruned_roots+=1;
                                    prune_root_counter_prune +=1;

                                    write_mol_to_file( dn_output_prefix,
                                                       "prune_root",
                                                       growing[x],
                                                       a, counter );
                                }
                            } else {
                                num_pruned_roots+=1;
                                prune_root_counter_prune +=1;
                                write_mol_to_file( dn_output_prefix,
                                                   "prune_root",
                                                   growing[x],
                                                   a, counter );
                            }
                        } else {
                            num_pruned_roots+=1;
                            prune_root_counter_prune +=1;
                            write_mol_to_file( dn_output_prefix,
                                               "prune_root",
                                               growing[x],
                                               a, counter );
                        }
                    }

                    // clear growing to prepare for the next round of fragments
                    growing.clear();
    
                } // end for each attachment point in the fragment
 
                root_TIME_END = std::chrono::system_clock::now();
                std::chrono::duration<double> root_elapsed_seconds = root_TIME_END - root_TIME_START;

 
                // THIS IS A WELL FORMATTED COUT STATEMENT
                cout <<"    Root ["   << setw(3) << std::right << i <<"],"
                     <<" sccf atts "  << setw(6) << std::right << num_att_root 
                     <<" |acc mols "  << setw(4) << std::right << root_mol_counter
                     <<" |cmpd "      << setw(4) << std::right << num_completed
                     <<" |fltrd cmpd "<< setw(4) << std::right << num_filtered_comp 
                     <<" |Prnd "      << setw(4) << std::right << num_pruned_roots
                     <<" |Kpt "       << setw(4) << std::right << root_frag_counter  
                     <<" |S "         << std::fixed << std::setprecision(2) << setw(8) << std::right 
                                      << root_elapsed_seconds.count() 
                     <<" |A/S "       << std::fixed << std::setprecision(2) << setw(6) << std::right 
                                      << (float) num_att_root / root_elapsed_seconds.count() 
                     << endl;

                num_all_frags_layer -= root_frag_counter;

                // put all the contents of layer onto next_layer - this clears up layer for the next 
                // fragment in root, and later we will sort / prune all the contents of next_layer to
                // pick fragments to put back onto root
                
                #ifdef BUILD_DOCK_WITH_RDKIT
                if (dn_drive_clogp || dn_drive_esol || dn_drive_qed || dn_drive_sa || dn_drive_stereocenters || dn_drive_pains || dn_drive_tpsa) {
                    dn_normal = false;
                    // Start interface with RDKit
                    RDTYPER rdprops;
                    for (unsigned int i = 0; i < layer.size(); ++i){
                        rdprops.calculate_descriptors( layer[i].mol, fragMap, true, PAINSMap);
                    }

                    int num_rej_mol = 0;  
                    //if (rejected.size()!=0) { num_rej_mol = rejected.size();
                    //} else { num_rej_mol = 0;}
                    int num_rej_mol_in_attFrag = 0;
                    float rej_avg_mw   = 0.0;
                    float rej_avg_esol = 0.0;
                    float rej_avg_logp = 0.0;
                    float rej_avg_tpsa = 0.0;
                    float rej_avg_qed  = 0.0;
                    float rej_avg_sa   = 0.0;
                    float rej_avg_stereo = 0.0;
                    float rej_avg_pains = 0.0;
                    

                    for (unsigned int x = 0; x < layer.size(); ++x){   
                        // Figure out which molecules fail which test 
                        // Conditional that drives growth
                        if (dn_drive_verbose){ 
                            //stringstream tmp_string;
                            string out_string; 
                            
                            //tmp_string >> out_string; 
                            //cout << "      " <<setw(20) << std::right<< out_string; 
                        }
                        if ( (counter+1) >= dn_start_at_layer ){
                            bool drive;
                            if (dn_drive_verbose){
                                cout  <<"        attached_frag["<<to_string(x)<<"]"<<std::endl;
                            }
                            drive = drive_growth( layer[x] ); 
                            if (drive) { next_layer.push_back(layer[x]); }
                            else { 
                                rejected.push_back(layer[x]); 
                                num_rej_mol_in_attFrag++;
                            }
                        // Accept all molecules built before dn_start_at_layer
                        } else if ( (counter+1) < dn_start_at_layer ){ 
                            next_layer.push_back(layer[x]); 
                            if (dn_drive_verbose){
                                cout << "        attached_frag["<<to_string(x)<<"] molecule_accepted via not start on dn_start_at_layer" << endl;
                            }
                        }
                        // Reject all others
                        else {
                            rejected.push_back(layer[x]);
                            num_rej_mol_in_attFrag++;
                            if (dn_drive_verbose){
                                cout << "        attached_frag["<<to_string(x)<<"]" <<" rejected because notDrive "<< std::endl;;
                            }
                        }

                        if (num_rej_mol < rejected.size()) {
                            num_rej_mol = rejected.size();
                            total_rej_molecules = rejected.size();
                        }
                    }
                    if (dn_drive_verbose){  
                        for ( int u = 0; u < rejected.size(); u++){
                            
                            rej_avg_esol += rejected[u].mol.esol;
                            rej_avg_logp += rejected[u].mol.clogp;
                            rej_avg_tpsa += rejected[u].mol.tpsa;
                            rej_avg_qed  += rejected[u].mol.qed_score;
                            rej_avg_sa   += rejected[u].mol.sa_score;
                            rej_avg_pains += rejected[u].mol.pns;
                            rej_avg_stereo += rejected[u].mol.num_stereocenters;
                        }
                        rej_avg_esol = rej_avg_esol / rejected.size();
                        rej_avg_logp = rej_avg_logp / rejected.size(); 
                        rej_avg_tpsa = rej_avg_tpsa / rejected.size();
                        rej_avg_qed  = rej_avg_qed  / rejected.size();
                        rej_avg_sa   = rej_avg_sa   / rejected.size();
                        rej_avg_pains = rej_avg_pains / rejected.size(); 
                        rej_avg_stereo = rej_avg_stereo / rejected.size();


                        //cout <<"        From RD_layer[" << (counter+1)<<"], the number of rejected molecules is " << num_rej_mol_in_attFrag << endl;   
                        cout <<"        RDKIT: The number of Rejected Molecules is "<<setw(4)<< std::left<< num_rej_mol << endl; 
                        if (num_rej_mol != 0){
                            cout <<"        "<<setw(20)<<left<<"Rejected_avg_Stereocenter is " << rej_avg_stereo   << endl;
                            cout <<"        "<<setw(20)<<left<<"Rejected_avg_tpsa is "  << rej_avg_tpsa << endl;
                            cout <<"        "<<setw(20)<<left<<"Rejected_avg_clogp is "  << rej_avg_logp << endl;
                            cout <<"        "<<setw(20)<<left<<"Rejected_avg_logS is "   << rej_avg_esol << endl;
                            cout <<"        "<<setw(20)<<left<<"Rejected_avg_qed is "    << rej_avg_qed  << endl;
                            cout <<"        "<<setw(20)<<left<<"Rejected_avg_SynthA is " << rej_avg_sa   << endl;
                            cout <<"        "<<setw(20)<<left<<"Rejected_avg_numPains is " << rej_avg_pains   << endl;
                        }

                    } else {
                        //cout <<"        From RD_layer[" << (counter+1)<<"]"<<"anchor["<<a<<"]"<<"Root[" <<i <<"],the number of rejected molecules is " << num_rej_mol_in_attFrag << endl;
                        cout <<"        RDKIT: the number of rejected molecules is " << num_rej_mol << endl; 
                    }
                } else {
                // Regular de novo
                    dn_normal = true;
                    for (int x=0; x<layer.size(); x++){
                        next_layer.push_back(layer[x]);
                    }
                }
                #else
                for (int x=0; x<layer.size(); x++){
                    next_layer.push_back(layer[x]);
                }

                #endif
    
                // clear the contents of layer
                layer.clear();
    
            } // end for each fragment in root
    
            // at this point, we have built one complete layer - clear root and sort / prune the
            // contents of next_layer to figure out what frags will move onto the next layer of growth
            root.clear();

    
            // choose the sort function based on input parameters
            //frag_sort(next_layer,fragment_sort);
            cout <<"      From all attachments events, " << num_att_layer 
                 << " fragments were collected in layer " 
                 << counter+1 << " prior to all pruning/filtering." <<std::endl;
            cout <<"      From accepted frags ( N= " << num_all_frags_layer 
                 << " ) in layer "<< counter+1 << ", " <<std::endl;
            cout <<"          Collected a total of " << next_layer.size() 
                 <<" candidate root fragments,";
            // once next_layer is sorted, remove some redundancy based on the Hungarian RMSD heuristic
            prune_h_rmsd_and_mw(next_layer, prune_rmsd_mw, a, counter);
            prune_rmsd_mw.clear();

            cout <<" then that number goes down to " 
                 <<next_layer.size() << " due to pruning by the RMSD heuristic."<< endl;

    
            // prepare to calculate some ensemble properties for printing to stdout
            float avg_mol_wt = 0.0;
            float avg_rot_bonds = 0.0;
            float avg_formal_charge = 0.0;
            float avg_num_HA = 0.0;
            float avg_num_HD = 0.0;
    
            // for everything in next_layer...
            for (int i=0; i<next_layer.size(); i++){
    
                // calculate some ensemble properties
                calc_mol_wt(next_layer[i].mol);
                calc_rot_bonds(next_layer[i].mol);
                calc_formal_charge(next_layer[i].mol);
                calc_num_HA_HD(next_layer[i].mol);
    
                // add them to this running total
                avg_mol_wt += next_layer[i].mol.mol_wt;
                avg_rot_bonds += next_layer[i].mol.rot_bonds;
                avg_formal_charge += next_layer[i].mol.formal_charge;
                //avg_num_HA += next_layer[i].mol.hb_acceptors;
                //avg_num_HD += next_layer[i].mol.hb_donors;
            }
    
            // If for some reason the next_layer size was 0, you don't want to get an accidental
            // divide by 0 error here
            if (next_layer.size() != 0){
                avg_mol_wt /= next_layer.size();
                cout <<"          The average molecular weight is " <<avg_mol_wt <<endl; 
                avg_rot_bonds /= next_layer.size();
                cout <<"          The average number of rotatable bonds is " <<avg_rot_bonds <<endl; 
                avg_formal_charge /= next_layer.size();
                cout <<"          The average formal_charge is " <<avg_formal_charge <<endl; 
                // avg_num_HA /= next_layer.size();
                // cout <<"      From this ensemble, the average hbond acceptors is " <<avg_num_HA <<endl;
                // avg_num_HD /= next_layer.size();
                // cout <<"      From this ensemble, the average hbond donors is " <<avg_num_HA <<endl;
                #ifdef BUILD_DOCK_WITH_RDKIT
                    if (dn_drive_verbose){
                        cout << "          The total RD_D3N_rejected molecules is " <<total_rej_molecules<<endl;
                    }
                #endif 

            }

            //JDB BTB Prune step. 
            if (dn_make_unique) {
                make_unique(completed_molecules, duplicate_dump);
            }

            // fout_completed_name is the filestream for complete molecules
            ostringstream fout_completed_name;
            fout_completed_name << dn_output_prefix<<".completed.denovo_build.mol2";
            fstream fout_completed;
/* JDB JDB JDB 
            if (counter > 0){
                fout_completed.open(fout_completed_name.str().c_str(), fstream::out|fstream::app);
                
            } else {
                fout_completed.open(fout_completed_name.str().c_str(), fstream::out);
            }
*/
            fout_completed.open(fout_completed_name.str().c_str(), fstream::out|fstream::app);
            fout_completed_name.clear();

            // Keep track of rejected mol on lower bound mw
            int mw_low_prune_counter=0;
            int complete_write_to_file_counter=0;
            for (int i=0; i<completed_molecules.size(); i++){
                if ((dn_MW_cutoff_type_hard && completed_molecules[i].mol.mol_wt > dn_lower_constraint_mol_wt) || \
                    (!dn_MW_cutoff_type_hard && mw_cutoff(completed_molecules[i],0) )){
                    //LEP - 2018.06.20 
                    ostringstream new_comment;
                    new_comment << denovo_name<< "_layer"<<counter+1<<"_"<< complete_write_to_file_counter;
                    completed_molecules[i].mol.title = new_comment.str();
                    new_comment.clear();
                    // write mols to file
                    fout_completed << completed_molecules[i].mol.current_data;
                    //LEP added in unique denovo name and scientific notation to floating point for output
                    fout_completed << DELIMITER << setw(STRING_WIDTH)
                                   << "Name:" << setw(FLOAT_WIDTH) 
                                   << completed_molecules[i].mol.title <<endl;

                    fout_completed << DELIMITER << setw(STRING_WIDTH)
                                   << "Molecular_Weight:" << setw(FLOAT_WIDTH) 
                                   << completed_molecules[i].mol.mol_wt << endl;

                    fout_completed << DELIMITER << setw(STRING_WIDTH)
                                   << "DOCK_Rotatable_Bonds:" << setw(FLOAT_WIDTH)
                                   << completed_molecules[i].mol.rot_bonds << endl;

                    fout_completed << fixed << DELIMITER << setw(STRING_WIDTH)
                                   << "Formal_Charge:" << setw(FLOAT_WIDTH)
                                   << completed_molecules[i].mol.formal_charge << endl;

                    fout_completed << fixed << DELIMITER << setw(STRING_WIDTH)
                                   << "HBond_Acceptors:" << setw(FLOAT_WIDTH) 
                                   << completed_molecules[i].mol.hb_acceptors << endl;

                    fout_completed << fixed << DELIMITER << setw(STRING_WIDTH) 
                                   << "HBond_Donors:" << setw(FLOAT_WIDTH)
                                   << completed_molecules[i].mol.hb_donors << endl;
                    #ifdef BUILD_DOCK_WITH_RDKIT
                    // Descriptors (beginning)
                    if (!dn_normal) {
                        std::vector<std::string> vecpnstmp{};
                        vecpnstmp = completed_molecules[i].mol.pns_name;
                        std::string molpns_name = std::accumulate(vecpnstmp.begin(), vecpnstmp.end(), std::string(""));
                        if (completed_molecules[i].mol.pns == 0 && molpns_name.empty() == true){
                           molpns_name = "NO_PAINS"; 
                        } else if(completed_molecules[i].mol.pns != 0 && molpns_name.empty() == false){
                        } else {
                           molpns_name = "ERROR_IN_PAINS_MATCHING_PLEASE_INVESTIGATE";
                        }
                        fout_molecules << fixed << DELIMITER << setw(STRING_WIDTH) 
                                       << "RD_num_arom_rings:" << setw(FLOAT_WIDTH) 
                                       << completed_molecules[i].mol.num_arom_rings << endl;

                        fout_molecules << fixed << DELIMITER << setw(STRING_WIDTH) 
                                       << "RD_num_alip_rings:" << setw(FLOAT_WIDTH) 
                                       << completed_molecules[i].mol.num_alip_rings << endl;

                        fout_molecules << fixed << DELIMITER << setw(STRING_WIDTH) 
                                       << "RD_num_sat_rings:" << setw(FLOAT_WIDTH) 
                                       << completed_molecules[i].mol.num_sat_rings << endl;

                        fout_molecules << fixed << DELIMITER << setw(STRING_WIDTH) 
                                       << "RD_Stereocenters:" << setw(FLOAT_WIDTH) 
                                       << completed_molecules[i].mol.num_stereocenters << endl;

                        fout_molecules << fixed << DELIMITER << setw(STRING_WIDTH) 
                                       << "RD_Spiro_atoms:" << setw(FLOAT_WIDTH) 
                                       << completed_molecules[i].mol.num_spiro_atoms << endl;

                        fout_molecules << DELIMITER << setw(STRING_WIDTH) 
                                       << "RD_LogP:" << setw(FLOAT_WIDTH) 
                                       << completed_molecules[i].mol.clogp << endl;

                        fout_molecules << fixed << DELIMITER << setw(STRING_WIDTH) 
                                       << "RD_TPSA:" << setw(FLOAT_WIDTH) 
                                       << completed_molecules[i].mol.tpsa << endl;

                        fout_molecules << DELIMITER << setw(STRING_WIDTH) 
                                       << "RD_SYNTHA:" << setw(FLOAT_WIDTH) 
                                       << completed_molecules[i].mol.sa_score << endl;

                        fout_molecules << DELIMITER << setw(STRING_WIDTH) 
                                       << "RD_QED:" << setw(FLOAT_WIDTH) 
                                       << completed_molecules[i].mol.qed_score << endl;

                        fout_molecules << DELIMITER << setw(STRING_WIDTH)  
                                       << "RD_LogS:" << setw(FLOAT_WIDTH) 
                                       << completed_molecules[i].mol.esol << endl;

                        fout_molecules << DELIMITER << setw(STRING_WIDTH) 
                                       << "RD_num_of_PAINS:" << setw(FLOAT_WIDTH)  
                                       << completed_molecules[i].mol.pns << endl; 

                        fout_molecules << DELIMITER << setw(STRING_WIDTH) 
                                       << "RD_PAINS_names:" << setw(FLOAT_WIDTH) << molpns_name << endl; 

                        fout_molecules << DELIMITER << setw(STRING_WIDTH) 
                                       << "RD_SMILES:" << setw(FLOAT_WIDTH) 
                                       << completed_molecules[i].mol.smiles << endl;
                    }// Descriptors (end)
                    #endif   
                    fout_completed <<DELIMITER<<setw(STRING_WIDTH)
                                   <<"Layer_Completed:"<<setw(FLOAT_WIDTH)
                                   <<counter+1 <<endl;

                    fout_completed <<DELIMITER<<setw(STRING_WIDTH)
                                   <<"Frag_String:"<<setw(FLOAT_WIDTH)
                                   <<completed_molecules[i].mol.energy << "\n"<<endl;

                    Write_Mol2(completed_molecules[i].mol, fout_completed);
                    complete_write_to_file_counter++;
                } else {
                    mw_low_prune_counter++;
                    filtered_comp_counter_prune += 1;
                    write_mol_to_file( dn_output_prefix, "filtered_comp",
                                       completed_molecules[i], a, counter);
                    //filtered_comp.push_back( completed_molecules[i] );
                }
            }
            fout_completed.close();
            // prepare to write some molecules to file
            int write_to_file_counter=0;
            // Collect everything to use for the next layer before assessing writes
            for (int i=0; i<next_layer.size(); i++){
                ostringstream new_comment;
                new_comment << denovo_name<< "_layer"<<counter+1<<"_"<< write_to_file_counter;
                next_layer[i].mol.title = new_comment.str();
                new_comment.clear();

                //write growth trees if enabled
                if (dn_write_growth_trees){

                    // Increment the index so we have a unique filename
                    growth_tree_index++;

                    // Open up a file stream for growth trees
                    ostringstream fout_growth_tree_name;
                    fout_growth_tree_name <<dn_output_prefix <<".growth_tree_" <<growth_tree_index <<".mol2";
                    fstream fout_growth_tree;
                    fout_growth_tree.open (fout_growth_tree_name.str().c_str(), fstream::out|fstream::app);
                    fout_growth_tree_name.clear();

                    // Copy the final dockmol into the growth tree
                    next_layer[i].frag_growth_tree.push_back(next_layer[i].mol);

                    // Write the growth tree to file
                    for (int j=0; j<next_layer[i].frag_growth_tree.size(); j++){
                        activate_mol(next_layer[i].frag_growth_tree[j]);
                        fout_growth_tree <<next_layer[i].frag_growth_tree[j].current_data << endl;
                        Write_Mol2(next_layer[i].frag_growth_tree[j], fout_growth_tree);
                    }
            
                    // Close the growth tree file stream
                    fout_growth_tree.close();
                }
            }

            //If randomly selecting max_root_size of fragments for the next layer
            if (dn_random_root_selection){
                cout << "          Random root selection turned on - keeping up to " 
                     << dn_max_root_size << " random fragments for the next layer." << endl;

                //If there are less than max_root_size fragments, just push them all to root.
                //otherwise, perform random selection
                if (dn_max_root_size > next_layer.size()){
                    for (int num = 0; num < next_layer.size(); num++){
                        root.push_back(next_layer[num]);
                    }
                } else {
                    std::vector<int> next_layer_sel_integers;
                    int num_frags_to_select = dn_max_root_size;

                    //randomly generate values up to the max_root_size
                    while (next_layer_sel_integers.size() < num_frags_to_select) {
                        int sel_val = rand() % next_layer.size();
                        if ((std::find(next_layer_sel_integers.begin(), 
                                       next_layer_sel_integers.end(), 
                                       sel_val)) == next_layer_sel_integers.end()){

                            next_layer_sel_integers.push_back(sel_val);
                        }
                    }

                    //push all the selected values to root
                    for (int num = 0; num < next_layer_sel_integers.size(); num++) {
                        int select_value = next_layer_sel_integers[num];
                        if (root.size() < dn_max_root_size) {
                            root.push_back(next_layer[select_value]);
                        }
                    }

                    // push into cand_root_ign after collecting the randomly collected root cands
                    for (int num = 0; num < next_layer.size(); num++){
                        if ((std::find(next_layer_sel_integers.begin(), 
                                       next_layer_sel_integers.end(), 
                                       num)) == next_layer_sel_integers.end()){
   
                            cand_root_ign_counter_prune += 1;
                            write_mol_to_file ( dn_output_prefix,
                                                "cand_root_ign",
                                                next_layer[num],
                                                a, counter
                                              );
                            //cand_root_ign.push_back(next_layer[num]);

                        }
                    }
                }


            //otherwise, just scrape off the top (ie highest scoring fragments)
            } else {
                for (int i=0; i<next_layer.size(); i++){
                    if (root.size() < dn_max_root_size){
                         root.push_back(next_layer[i]);
                    } else {
                        cand_root_ign_counter_prune += 1;
                        write_mol_to_file ( dn_output_prefix,
                                            "cand_root_ign",
                                            next_layer[i],
                                            a, counter
                                          );
                    }
                }
            }

            // Keeping track of how many molecules failed to write to file on lower bound mw cutoff.
            // If dn_lower_constraint_mol_wt = 0 do not bother with this statement.
            if (dn_ga_flag == false && (dn_lower_constraint_mol_wt > 0)){
               cout << "          While writing molecules to file, " 
                    << mw_low_prune_counter << " molecules were rejected on the lower bound molecular weight." << endl;
            }
            
            #ifdef BUILD_DOCK_WITH_RDKIT
            if (dn_save_all_molecules) {


                // Print rejected molecules to file

                for (unsigned int i = 0; i < rejected.size(); ++i) {

                    bool rejected_flag = true;

                    for (int j=0; j<rejected[i].aps.size(); j++){

                        string bond_type;

                        // Loop over all bonds in frag1
                        for (int k=0; k<rejected[i].mol.num_bonds; k++){
                            if (rejected[i].mol.bonds_origin_atom[k] == rejected[i].aps[j].dummy_atom ||
                                rejected[i].mol.bonds_target_atom[k] == rejected[i].aps[j].dummy_atom   ){   
                                if (rejected[i].mol.bonds_origin_atom[k] == rejected[i].aps[j].heavy_atom ||
                                    rejected[i].mol.bonds_target_atom[k] == rejected[i].aps[j].heavy_atom   ){   

                                    // This is the bond you're looking for...
                                    bond_type = rejected[i].mol.bond_types[k];
                                    break;
                                }
                            }
                        }

                        // Only cap it if the dummy is on a single bond 
                        if ( bond_type == "1" ){

                            dummy_to_H( rejected[i], 
                                        rejected[i].aps[j].heavy_atom, 
                                        rejected[i].aps[j].dummy_atom );

                        // Otherwise, ignore this molecule
                        } else {
                            rejected_flag = false;
                            break;
                        }

                    }
                    //Run prune_dump mols through multi torenv check to not make garbage
                    if (valid_torenv_multi(rejected[i])){
                        // If there was no problem capping the molecules, then write to file
                        if (rejected_flag){
                            //reactivate all atoms so it can be recognized by  Write_Mol2() for mol2 writing 
                            for (int q =0; q<rejected[i].mol.num_atoms;q++){
                                rejected[i].mol.atom_active_flags[q] = true;
                            }

                            calc_mol_wt(rejected[i].mol);
                            calc_rot_bonds(rejected[i].mol);
                            calc_formal_charge(rejected[i].mol);
                            calc_num_HA_HD(rejected[i].mol);

                            // Prepare title 
                            ostringstream new_title;
                            new_title << denovo_name << "_layer" << counter + 1 << "_rejected_" << i;
                            rejected[i].mol.title = new_title.str();
                            new_title.clear();

                            std::vector<std::string> vecpnstmp{};
                            vecpnstmp = rejected[i].mol.pns_name;
                            std::string molpns_name = std::accumulate(vecpnstmp.begin(), vecpnstmp.end(), std::string(""));                
                            if (rejected[i].mol.pns == 0 && molpns_name.empty() == true){
                                molpns_name = "NO_PAINS"; 
                            }else if(rejected[i].mol.pns != 0 && molpns_name.empty() == false) {}
                            else{
                                molpns_name = "ERROR_IN_PAINS_MATCHING_PLEASE_INVESTIGATE";
                            }

                            // Write mols and attributes to file
                            fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Name:" << setw(FLOAT_WIDTH) << rejected[i].mol.title << endl;
                            fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Molecular_Weight:" << setw(FLOAT_WIDTH) << rejected[i].mol.mol_wt << endl;
                            fout_rejected << DELIMITER << setw(STRING_WIDTH) << "DOCK_Rotatable_Bonds:" << setw(FLOAT_WIDTH) << rejected[i].mol.rot_bonds << endl;
                            fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Formal_Charge:" << setw(FLOAT_WIDTH) << rejected[i].mol.formal_charge << endl;
                            fout_rejected << fixed << DELIMITER << setw(STRING_WIDTH) << "HBond_Acceptors:" << setw(FLOAT_WIDTH) << rejected[i].mol.hb_acceptors << endl;
                            fout_rejected << fixed << DELIMITER << setw(STRING_WIDTH) << "HBond_Donors:" << setw(FLOAT_WIDTH) << rejected[i].mol.hb_donors << endl;
                            fout_rejected << fixed << DELIMITER << setw(STRING_WIDTH) << "RD_num_arom_rings:" << setw(FLOAT_WIDTH) << rejected[i].mol.num_arom_rings << endl;
                            fout_rejected << fixed << DELIMITER << setw(STRING_WIDTH) << "RD_num_alip_rings:" << setw(FLOAT_WIDTH) << rejected[i].mol.num_alip_rings << endl;
                            fout_rejected << fixed << DELIMITER << setw(STRING_WIDTH) << "RD_num_sat_rings:" << setw(FLOAT_WIDTH) << rejected[i].mol.num_sat_rings << endl;
                            fout_rejected << fixed << DELIMITER << setw(STRING_WIDTH) << "RD_Stereocenters:" << setw(FLOAT_WIDTH) << rejected[i].mol.num_stereocenters << endl;
                            fout_rejected << fixed << DELIMITER << setw(STRING_WIDTH) << "RD_Spiro_atoms:" << setw(FLOAT_WIDTH) << rejected[i].mol.num_spiro_atoms << endl;
                            fout_rejected << DELIMITER << setw(STRING_WIDTH) << "RD_LogP:" << setw(FLOAT_WIDTH) << rejected[i].mol.clogp << endl;
                            fout_rejected << fixed << DELIMITER << setw(STRING_WIDTH) << "RD_TPSA:" << setw(FLOAT_WIDTH) << rejected[i].mol.tpsa << endl;
                            fout_rejected << DELIMITER << setw(STRING_WIDTH) << "RD_SYNTHA:" << setw(FLOAT_WIDTH) << rejected[i].mol.sa_score << endl;
                            fout_rejected << DELIMITER << setw(STRING_WIDTH) << "RD_QED:" << setw(FLOAT_WIDTH) << rejected[i].mol.qed_score << endl;
                            fout_rejected << DELIMITER << setw(STRING_WIDTH) << "RD_LogS:" << setw(FLOAT_WIDTH) << rejected[i].mol.esol << endl;
                            fout_rejected << DELIMITER << setw(STRING_WIDTH) << "RD_num_of_PAINS:" << setw(FLOAT_WIDTH) << rejected[i].mol.pns << endl; 
                            fout_rejected << DELIMITER << setw(STRING_WIDTH) << "RD_PAINS_names:" << setw(FLOAT_WIDTH) << molpns_name << endl; 
                            fout_rejected << DELIMITER << setw(STRING_WIDTH) << "RD_SMILES:" << setw(FLOAT_WIDTH) << rejected[i].mol.smiles << endl;

                            if (rejected[i].mol.fail_stereo) { fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Failed_stereo:" << setw(FLOAT_WIDTH) << "TRUE" << endl; }
                            else { fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Failed_stereo:" << setw(FLOAT_WIDTH) << "FALSE" << endl; }
                            if (rejected[i].mol.fail_clogp) { fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Failed_LogP:" << setw(FLOAT_WIDTH) << "TRUE" << endl; }
                            else { fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Failed_LogP:" << setw(FLOAT_WIDTH) << "FALSE" << endl; }                
                            if (rejected[i].mol.fail_tpsa) { fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Failed_TPSA:" << setw(FLOAT_WIDTH) << "TRUE" << endl; }
                            else { fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Failed_TPSA:" << setw(FLOAT_WIDTH) << "FALSE" << endl; }
                            if (rejected[i].mol.fail_sa) { fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Failed_SYNTHA:" << setw(FLOAT_WIDTH) << "TRUE" << endl; }
                            else { fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Failed_SYNTHA:" << setw(FLOAT_WIDTH) << "FALSE" << endl; }
                            if (rejected[i].mol.fail_qed) { fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Failed_QED:" << setw(FLOAT_WIDTH) << "TRUE" << endl; }
                            else { fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Failed_QED:" << setw(FLOAT_WIDTH) << "FALSE" << endl; }
                            if (rejected[i].mol.fail_esol) { fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Failed_LogS:" << setw(FLOAT_WIDTH) << "TRUE" << endl; }
                            else { fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Failed_LogS:" << setw(FLOAT_WIDTH) << "FALSE" << endl; }
                            if (rejected[i].mol.fail_pains) { fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Failed_numPains:" << setw(FLOAT_WIDTH) << "TRUE" << endl; }
                            else { fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Failed_numPains:" << setw(FLOAT_WIDTH) << "FALSE" << endl; }
                            fout_rejected << DELIMITER << setw(STRING_WIDTH) << "Frag_String:" << setw(FLOAT_WIDTH) << rejected[i].mol.energy << "\n"<<endl;

                            // Use dockmol Write_Mol2 function
                            Write_Mol2( rejected[i].mol, fout_rejected );   
                        }
                    }    
                }
            
            // Clear up rejected vector after writing to Mol2.
            rejected.clear();
            } else { rejected.clear(); }               
            #endif


            cout <<"          " << "-------------------------------------------" << std::endl;
            //where we start on writing out the results of the current layer
            if (dn_make_unique){
                cout <<"          " << setw(4) << std::left << duplicate_dump.size() <<" duplicate fragments were found."<< endl;
            }

            // if duplicate_
            if (dn_write_out_duplicates){
               for ( Fragment & dupfrag: duplicate_dump) {
	           write_mol_to_file(dn_output_prefix, "duplicate_dump", dupfrag, a, counter); 
               }
	    }
            duplicate_dump.clear();

            // write out the pruning stats here
            cout <<"          " << setw(4) << std::left 
                 << filtered_comp_counter_prune <<" complete molecules were filtered" << endl;
            cout <<"          " << setw(4) << std::left
                 << prune_root_counter_prune <<" root fragments were filtered" << endl;

            cout <<"          " << setw(4) << std::left << complete_write_to_file_counter 
                 << " completed molecules were written to file" << std::endl;
            cout <<"          " << "-------------------------------------------" << std::endl;
            
            cout <<"          " << setw(4) << std::left<< prune_rmsd_mw_counter_prune 
                 <<" root frags intended for the next layer were further " 
                 << "pruned via rmsd."<<endl;
            cout <<"          " << setw(4) << std::left<< cand_root_ign_counter_prune  
                 <<" candidate root frags intended for the next layer were ignored " 
                 << "due to limit on root size."<<endl;

            if ( failed_vtm_prune_rmsd_mw > 0 || failed_cap_H_prune_rmsd_mw > 0 )
                cout <<"          For " << setw(15) << std::right 
                     << "prune_rmsd_mw frags, " << setw(4) 
                     << std::right<< failed_cap_H_prune_rmsd_mw
                     << " frags couldn't be written out due to failed H-capping protocol\n"          
                     <<"                                     "<< setw(4) 
                     << std::right << failed_vtm_prune_rmsd_mw
                     << " frags couldn't be written out due to failed torsion checks on H-caps." 
                     << std::endl;
            if ( failed_vtm_prune_root > 0 || failed_cap_H_prune_root > 0) {
                cout <<"          For " << setw(15) << std::right 
                     << "prune_root frags, " << setw(4) 
                     << std::right<< failed_cap_H_prune_root
                     << " frags couldn't be written out due to failed H-capping protocol\n"          
                     <<"                                     "<< setw(4) 
                     << std::right << failed_vtm_prune_root
                     << " frags couldn't be written out due to failed torsion checks on H-caps." 
                     << std::endl;
            }    
            if ( failed_vtm_cand_root_ign > 0 || failed_cap_H_cand_root_ign > 0) {
                cout <<"          For " << setw(15) << std::right 
                     << "cand_root_ign frags, " << setw(4) 
                     << std::right<< failed_cap_H_cand_root_ign
                     << " frags couldn't be written out due to failed H-capping protocol\n"          
                     <<"                                     "<< setw(4) 
                     << std::right << failed_vtm_cand_root_ign
                     << " frags couldn't be written out due to failed torsion checks on H-caps." 
                     << std::endl;
            }
                
            if ( failed_vtm_filtered_comp > 0 || failed_cap_H_filtered_comp > 0 ) {
                cout <<"          For " << setw(15) << std::right 
                     << "filtered_comp frags, " << setw(4) 
                     << std::right<< failed_cap_H_filtered_comp
                     << " frags couldn't be written out due to failed H-capping protocol\n"          
                     <<"                                     "<< setw(4) 
                     << std::right << failed_vtm_filtered_comp
                     << " frags couldn't be written out due to failed torsion checks on H-caps." 
                     << std::endl;
            }

            if ( failed_vtm_duplicate_dump > 0 || failed_cap_H_duplicate_dump> 0 ) {
                cout <<"          For " << setw(15) << std::right 
                     << "duplicate_dump frags, " << setw(4) 
                     << std::right<< failed_cap_H_duplicate_dump
                     << " frags couldn't be written out due to failed H-capping protocol\n"          
                     <<"                                     "<< setw(4) 
                     << std::right << failed_vtm_duplicate_dump
                     << " frags couldn't be written out due to failed torsion checks on H-caps." 
                     << std::endl;
            }

            prune_rmsd_mw.clear();

            prune_root_counter_prune = 0;
            filtered_comp_counter_prune = 0;
            prune_rmsd_mw_counter_prune = 0;
            cand_root_ign_counter_prune = 0;

            failed_cap_H_prune_rmsd_mw = 0;
            failed_cap_H_prune_root = 0;
            failed_cap_H_cand_root_ign = 0;
            failed_cap_H_filtered_comp = 0;
            failed_cap_H_duplicate_dump = 0;

            failed_vtm_prune_rmsd_mw = 0;
            failed_vtm_prune_root = 0;
            failed_vtm_cand_root_ign = 0;
            failed_vtm_filtered_comp = 0;
            failed_vtm_duplicate_dump = 0;

            // clear the prune dump vector for the next round
            cout <<std::endl <<"    Layer " << counter+1 << " completed! Moving on the next layer with "<<root.size() 
                 <<" root fragments." << endl;
    
            // Print out all the contents of root, this will be useful for restarts CDS 
            if (dn_write_checkpoints || counter==(dn_max_grow_layers-1)){
    
                // Filename for filehandle
                ostringstream fout_root_name;
                fout_root_name <<dn_output_prefix <<".anchor_" <<(a+1) <<".root_layer_" <<counter+1 <<".mol2";
                fstream fout_root;
    
                // Open the filehandle
                fout_root.open (fout_root_name.str().c_str(), fstream::out|fstream::app);
                for (int i=0; i<root.size(); i++){
                    fout_root << root[i].mol.current_data <<endl;
                    Write_Mol2( root[i].mol, fout_root );
                }
    
                // Clear some memory
                fout_root_name.clear();
                fout_root.close();
            }
    
            // clear next_layer and increment the loop counter
            next_layer.clear();
            counter++;
    
            // Annealing schedule for sampling from the graph, scales to 0
            if (dn_sampling_method_graph){
    
                // dn_graph_temperature scales at [ Tn = Ti / 2^i ]
                dn_graph_temperature = dn_temp_begin / ( pow (2.0, counter) );

                // When the next layer is the final layer, set temp to 0
                if ( (dn_max_grow_layers-counter) == 1) { dn_graph_temperature = 0; }
            }

            // dn_pruning_conformer_score_cutoff scales at [ Nn = Ni / Input_param ^ i ]
            dn_pruning_conformer_score_cutoff = dn_pruning_conformer_score_cutoff_begin / ( pow (dn_pruning_conformer_score_scaling_factor, counter) );
  
            // END CLOCK NOW
            end_TIME = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end_TIME - start_TIME;
            
            std::cout << "    For Layer " << counter << ", time duration is: " 
                      << elapsed_seconds.count() << " seconds with " << (float) num_att_layer / elapsed_seconds.count() 
                      << " with attachments per second."<<std::endl <<std::endl;

         
 
        } // end main while loop
    
    
        // clear root for the next round of anchor orients
        root.clear();
    
    } // end for each anchor
    cout <<"There are no more anchors. Exiting." <<endl <<endl;

    ////////////////////////////////////////////////////////////////////////////
    // End Main Loop                                                          //
    ////////////////////////////////////////////////////////////////////////////


    // clean up any filehandles
    fout_molecules.close();

    return;

} // end DN_Build::build_molecules()
// ++++++++++++++++++++++++++++++++++++++++++++++++
// Select a fragment from the fragment library based on frequency
// by using a random number to sample from a series of bins.
// Only called when frag_frequency_selection is selected.
void
DN_Build::new_select_frag_by_frequency(vector <Fragment> & fraglib, vector<int> & available_fragment_indices,
                                       vector<int> & selected_fragments, int num_picks) {
    bool first_frag_picked = false;
    for (int pick = 0; pick < num_picks; pick++){
        //cout << pick << endl;
        double freq_rn = ((double)rand() / (RAND_MAX));
        int sel_val=0;
        for (int i=1; i<ordered_fragment_frequencies.size(); ++i) {
            //Assess what bin we're actually inside - right wall of bin is considered bin i,
            // left wall would be considered bin i-1.
            if (freq_rn <= ordered_fragment_frequencies[i] && freq_rn > ordered_fragment_frequencies[i-1]) {
                //cout << ordered_fragment_frequencies[i-1] << " < " << freq_rn 
                //    << " < " << ordered_fragment_frequencies[i] << endl;
                
                //cout << "TRACK: " << ordered_fragments[i] << endl; //JDB DEBUG

                vector < int >::iterator avail_pos = available_fragment_indices.begin();
                while (avail_pos != available_fragment_indices.end()){

                    //traverse available fragments until greater than or equal to index
                    if ((*avail_pos) < i-1) {
                        avail_pos++;
                        continue;
                    }
                    // If greater than i, and first bin selected, select again since there's no bin to go back to
                    if (((*avail_pos) > i-1) && (avail_pos == available_fragment_indices.begin())) {
                        pick -= 1;
                        break;
                    } 

                    //if greater than, go one selection back
                    if ((*avail_pos) > i-1) {
                        avail_pos -= 1;
                    }
                    //push and erase
                    selected_fragments.push_back((*avail_pos));
                    available_fragment_indices.erase(avail_pos);
                    break;
                }

                /*
                //traverse available fragments until you are greater than or equal to 
                while (((*avail_pos) < i-1) && (avail_pos != available_fragment_indices.end())) {
                    avail_pos++;
                }
                cout << i-1 << " < " << (*avail_pos) 
                    << " < " << i << endl;
                if ((*avail_pos) == i-1) {
                    sel_val = (*avail_pos);
                } else {
                    if (avail_pos == available_fragment_indices.begin()) {
                        pick -=1;
                        break;
                    }
                    sel_val = (*(avail_pos-1));
                }

                cout << (*avail_pos) << "\t" << i-1 << "\t" << fraglib[i-1].mol.energy << "\t" << freq_rn << endl;
                //cout << "The beginning" << endl;

                cout << "Pos before " << (*avail_pos) << "   " << sel_val << endl;
                selected_fragments.push_back((*avail_pos));
                //cout << "The middle " << endl;
                cout << "Pos after " << (*avail_pos) << "   " << sel_val << "   " << i << endl;
                //cout << (*(avail_pos-1)) << "    " << (*(avail_pos)) << "    " << (*(avail_pos+1)) << "   " << sel_val << endl;
                cout << ordered_fragment_frequencies.size() << endl;
                cout << (*(available_fragment_indices.end())) << endl;
                available_fragment_indices.erase(avail_pos);

                cout << "The end :( " << endl;
                break;*/
            }  
        }
    }
    return;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++
// Select a fragment from the fragment library based on frequency
// by using a random number to sample from a series of bins.
// Only called when frag_frequency_selection is selected.
int
DN_Build::select_frag_by_frequency(vector <Fragment> & fraglib ){
    double freq_rn = ((double)rand() / (RAND_MAX));
    int sel_val=0;
    for (int i=1; i<ordered_fragment_frequencies.size(); ++i) {
        //First case: We're in the first bin.
        if (i == 1) { 
            if (freq_rn < ordered_fragment_frequencies[i]){
                //    cout << "TRACK:" << ordered_fragments[i] << endl; //JDB DEBUG

                //compare selection against entire library to find frag obj
                for(int frag_ind=0; i<fraglib.size(); frag_ind++){
                    if (fraglib[frag_ind].mol.energy == ordered_fragments[i]){
                        sel_val=frag_ind;
                        break;
                    }
                }
            }
        }
        // Second case: We're in between two bins
        if (freq_rn < ordered_fragment_frequencies[i] && 
            freq_rn > ordered_fragment_frequencies[i-1]) {
           /*cout << ordered_fragment_frequencies[i-1] << " < " << freq_rn 
                << " < " << ordered_fragment_frequencies[i] << endl;
            */
            //cout << "TRACK: " << ordered_fragments[i] << endl; //JDB DEBUG

            //compare selection against entire library to find frag obj
            for(int frag_ind=0; i<fraglib.size(); frag_ind++){
                if (fraglib[frag_ind].mol.energy == ordered_fragments[i]){
                    sel_val=frag_ind;
                    break;
                }
            }
        }
    
    }
    return sel_val;
}

// +++++++++++++++++++++++++++++++++++++++++
// This function will orient a fragment and return filtered / clustered / pruned results on the
// root vector
void
DN_Build::orient_fragments( vector <Fragment> & root, Fragment & frag, Master_Score & score,
                            Simplex_Minimizer & simplex, AMBER_TYPER & typer, Orient & orient )
{

    // This vector will hold orients
    vector <Fragment> orients;

    // Generate the list of atom-sphere matches
    orient.match_ligand(frag.mol);

    // Transforms frag.mol to one of the matches in the list
    while (orient.new_next_orientation(frag.mol)){

        // Make sure this orient is scoreable, e.g., inside the grid
        if (score.use_score && score.compute_primary_score(frag.mol)){

            // Minimize the anchor position
            simplex.minimize_rigid_anchor(frag.mol, score);

            // If it is still a valid score, then push it on to orients
            if (score.compute_primary_score(frag.mol)){
 
                orients.push_back(frag);
            }
        }
    }

    // Sort the orients by score
    frag_sort(orients, fragment_sort);

    // First filter by the score
    // (the value of 1000.0 is also hard-coded in anchor and grow when clustering is on)
    for (int i=0; i<orients.size(); i++){
        if (orients[i].mol.current_score > 1000.0){ orients[i].used = true; }
    }

    // Cluster post-minimization orients by rank/RMSD heuristic
    float rmsd;

    // For every fragment in orients
    for (int i=0; i<orients.size(); i++){

        // If the fragment is not currently 'used'
        if (!orients[i].used){

            // Then check the next fragment
            for (int j=i+1; j<orients.size(); j++){

                // If it is also not 'used'
                if (!orients[j].used){

                    // Then calculate the rmsd between the two
                    rmsd = calc_fragment_rmsd(orients[i], orients[j]);
                    rmsd = MAX(rmsd, 0.001);

                    // If the rmsd is below the specified cut-off, they are in the same 'cluster'
                    if ((float) j / rmsd > dn_pruning_clustering_cutoff) {
                                orients[j].used = true;
                    }
                }
            }
        }
    }

    // Move up to N onto new vector. Here, N is max_orients, but it could also be max_root_size
    for (int i=0; i<orients.size(); i++){
        if (!orients[i].used && i<orient.max_orients){

            root.push_back(orients[i]);
        }
    }

    // Print some status updates to stdout
    //cout <<"Number of orients generated = " <<orients.size() <<endl;
    //cout <<"Number of orients that passed the filter and prune = " <<root.size() <<endl <<endl;

    // No longer need all orients, because the ones that will be used are now in root
    orients.clear();


    // Write orients to file if the user specified to do so
    if (dn_write_orients){

        // Open up a file stream for orients
        ostringstream fout_orients_name;
        fout_orients_name <<dn_output_prefix <<".orients.mol2";
        fstream fout_orients;
        fout_orients.open (fout_orients_name.str().c_str(), fstream::out|fstream::app);
        fout_orients_name.clear();

        // Write each orient to file
        for (int i=0; i<root.size(); i++){
            fout_orients <<root[i].mol.current_data << endl;
            Write_Mol2(root[i].mol, fout_orients);
        }

        // Close the orients file stream
        fout_orients.close();
    }

    return;

} // end DN_Build::orient_fragments()
void
DN_Build::accept_frags_by_frequency(vector <int> & selected_fragments){

    vector< int >::iterator pos = selected_fragments.begin();

    //because of the placeholder value, selected indices + 1 = number line index (i)
    //bin wall at i (which is the right wall) is considered bin selected index.
    // left wall is bin i-1.

    //Because each bin is DISCRETE, where one fragment keys to one bin, there's no
    // bin wall checking necessary like in torsion acceptance, just a simple less-than check.
    int number_removed = 0;
    //cout << endl << "Started with " << selected_fragments.size() << " fragments" << endl;
    while (pos != selected_fragments.end()){
        //cout << ordered_fragments[*pos] << endl;
        //cout << ordered_fragment_frequencies[*pos] << endl;

        //select random number between 1 and 0
        double freq_rn = ((double)rand() / (RAND_MAX));
        float right_wall_selected = ordered_fragment_frequencies[*pos+1];

        //check if random number is greater than or equal to right bin wall
        //If yes, keep fragment.
        if ((freq_rn >= right_wall_selected)) {
            //cout << "RANDOM NUMBER    " << freq_rn << " >= " << right_wall_selected << "   " << ordered_fragments[*pos+1] << endl;
            pos++;
        //if not, delete fragment and move to next.
        } else {
            //cout << "RANDOM NUMBER    " << freq_rn << " < " << right_wall_selected << endl;
            pos = selected_fragments.erase(pos);
            number_removed++;
        }
    }

    //cout << "Ended with " << selected_fragments.size() << " fragments" << endl;
    //cout << "Number removed: " << number_removed << endl << endl;
}

//generates initial list of available fragment indices for this set of attachment attempts
//allows us to get a solely unique set of random selections without doing a bunch of 
//iterating through vectors later
//If not first selection, will select num_picks more fragments to try from list of indices
void
DN_Build::select_frags_from_fraglib(vector <Fragment> & fraglib, vector<int> & available_fragment_indices,
                                    vector<int> & selected_fragments, int & num_picks, bool first_selection, bool last_layer ){
    
    //if this is a fresh set of selections, repopulate the fragment indices
    if (first_selection) {
        //ensure we're starting from fresh vectors
        available_fragment_indices.clear();
        for(int fill_val = 0; fill_val < fraglib.size(); fill_val++){
            //If we're at the last layer, we *only* want to attach sidechains (ie ap = 1)
            if (last_layer) {
                if (fraglib[fill_val].aps.size() > 1) {
                    continue;
                }
            }
            available_fragment_indices.push_back(fill_val);
        }
    }


    //selected_fragments is cleared anew each time, since we're pushing selections to it
    selected_fragments.clear();
    if (dn_sel_frag_by_freq_bool) {
        new_select_frag_by_frequency( fraglib,  available_fragment_indices, selected_fragments, num_picks);
    } else {
        //select num_picks worth of fragments and push them to the 
        for (int selection_loop_counter = 0; selection_loop_counter < num_picks;
            selection_loop_counter++) {
            //exit this loop if there are no more available fragments to pick from for some reason
            //can happen if fraglib.size() < num_picks
            if (available_fragment_indices.size() < 1){
                break;
            }

            int dnrn;
            dnrn = rand() % available_fragment_indices.size();

            //push the element from available to selected
            selected_fragments.push_back(available_fragment_indices[dnrn]);
            available_fragment_indices.erase(available_fragment_indices.begin()+dnrn);
        }       

    }

    if (dn_acc_frag_by_freq_bool) {
        accept_frags_by_frequency(selected_fragments);
        //cout << "Size outside " << selected_fragments.size() << endl;
    }
}

// +++++++++++++++++++++++++++++++++++++++++
// Sample the contents of a given isofragment library exhaustively
void
DN_Build::sample_isofraglib_isoswap( Fragment & layer_frag, int j, vector <Fragment> & fraglib,
                                        vector <Fragment> & growing_ref, Master_Score & score, 
                                        Simplex_Minimizer & simplex, AMBER_TYPER & typer, bool last_layer,
					Iso_Table::Iso_Tab & iso_table, int layer_num )
{
 
    // just for the attachments of the new isosteres
    double working_mat[3][3];

    // If the fragment library is empty, return without doing anything
    if ( fraglib.size() == 0 || num_att_root > dn_max_sccful_att_per_root ) { return; }


    // While N choices...
    int choice=0;
    int passes=0;
    int dnrn;
    int num_sid=0;
    int num_lnk=0;
    int num_scf=0;

    int num_sid_h=0;
    int num_lnk_h=0;
    int num_scf_h=0;
    int growing_ref_size_before_aps = growing_ref.size();     

    std::vector <int> picked {};

    Fragment layer_frag_iso = layer_frag;


    // number of attachments per root attachment point
    //int overall_built_counter = 0;
    //int overall_hypo_counter = 0;

    // needed if isoswap is turning on the verbose
    std::ofstream        fout_isofrag;
    std::ofstream        before_fout_isofrag;
    std::stringstream    base;
    std::stringstream    name_file;
    std::string          name_file_str;
    std::stringstream    before_name_file;
    std::string          before_name_file_str;
    std::string          boolean_write = "";
    int                  rank_num      = -1;
    float                rate_iso_att  = -1.0;
    if (dn_iso_print_out.compare("yes")==0) {
        base << dn_iso_output_verbose_path << "/verbose_iso_out";
        string base_str = base.str();
        mkdir(base_str.c_str(), 0744);
 
        //This where we make the file name
        name_file << base_str << "/after_ver_layer_"<< layer_num << ".mol2";
        name_file_str = name_file.str();
        before_name_file << base_str << "/before_ver_layer_"<< layer_num << ".mol2";
        before_name_file_str = before_name_file.str();

        fout_isofrag.open(name_file_str,std::ios_base::app);
        before_fout_isofrag.open(before_name_file_str,std::ios_base::app);
    }


    // For every fragment in the fragment library (i.e. sidechains, linkers, scaffolds)
    while (choice < dn_num_rand_head_picks){
  
        // The passes int counts every time we either keep or skip a fragment from the libary, so
        // once its size reaches the fraglib.size(), we have looked at everything
        // Also, if the num of att per root exceeds parameter, break
        if (passes >= fraglib.size() || num_att_root > dn_max_sccful_att_per_root ) { 
            break; 
        }

        // Choose a random starting position, fraglib[dnrn]
        // (rand has been seeded if and only if the minimizer is turned on)
        dnrn=rand() % fraglib.size();
        //LEP - verbose dn frag selection random seed output
        if (verbose) cout << "De novo fragment selection: " << dnrn << endl; 
        // Keep track of things we try so we don't try the same thing twice
        if (picked.size() > 0){
            if (picked.size() >= fraglib.size()){
               picked.clear();
               break; 
            }
            for (int p=0; p<picked.size(); p++){
                if (dnrn == picked[p]){ continue; }
            }
        }

        picked.push_back(dnrn);

        // (1) Would create too many scaffolds per layer
        if ( fraglib[dnrn].aps.size() > 2 &&
             layer_frag.scaffolds_this_layer >= dn_max_scaffolds_per_layer ){ 
           if (choice + 1 >= dn_num_rand_head_picks && 
                       (growing_ref.size() == growing_ref_size_before_aps)){
               passes++;
               continue;
           }
           passes++;
           choice++;
           continue; 
        }


        // (2) Would give too many current attachment points
        int temp_aps = fraglib[dnrn].aps.size() + layer_frag.aps.size() - 2;
        if (temp_aps > dn_max_current_aps) { 
           if (choice + 1 >= dn_num_rand_head_picks && 
                       (growing_ref.size() == growing_ref_size_before_aps)){
               passes++;
               continue;
           }
           passes++;
           choice++;
           continue; 
        }

        // (3) Rotatable bond cutoff
        // BCF assuming rot_bon does not take into account attachment points
        // will have at least one new rot bond per attachment point, on each fragment
        // except when combined two aps will become one rot bond.
        // Floor for remaining rot_bonds if all remaining aps are capped with sidechains.
        int temp_rbs = layer_frag.mol.rot_bonds + temp_aps;
        if (temp_rbs > dn_constraint_rot_bon) { 
           if (choice + 1 >= dn_num_rand_head_picks && 
                       (growing_ref.size() == growing_ref_size_before_aps)){
               passes++;
               continue;
           }
           passes++;
           choice++;
           continue; 
        }

        // (4) Molecular weight cutoff (MW computed for layer frag at last step, MW precomputed for all frags)
        //int temp_mw = fraglib[dnrn].mol.mol_wt + layer_frag.mol.mol_wt;
        //if (temp_mw > dn_constraint_mol_wt) { continue; }

        // (5) Formal Charge Cutoff
        int temp_fc = fraglib[dnrn].mol.formal_charge + layer_frag.mol.formal_charge;
        if (fabs(float(temp_fc)) > (fabs(float(dn_constraint_formal_charge))+0.1)) { 
           if (choice + 1 >= dn_num_rand_head_picks && 
                       (growing_ref.size() == growing_ref_size_before_aps)){
               passes++;
               continue;
           }
           passes++;
           choice++;
           continue; 
        }

        bool bond_compared = false;
        int dummy_atom_head = -1;
        int heavy_atom_tar = -1;
        float new_rad_head = 0.0;
        Fragment iso_head;

	// do the isoalignments PAK
	std::pair<std::vector<Iso_Acessory::Scored_Triangle>,
              std::vector<Fragment>> vec_iso_frags = {};

        // If frag doesn't exist in the iso_table
        // you must align and set them in the iso_table
        if (!iso_table.check_if(fraglib[dnrn]) && 
            dn_iso_fraglib.compare("onthefly")==0){
	    vec_iso_frags =  frag_iso_align(fraglib[dnrn]);
            iso_table.set(fraglib[dnrn],vec_iso_frags.second);
	} else {
            vec_iso_frags.second = iso_table.get(fraglib[dnrn], dn_iso_num_gets);
        } 
  
        bool iso_head_skip = false;
        iso_head = fraglib[dnrn];

        if (vec_iso_frags.second.size() <= 1){
            iso_head_skip = true;
        } else {
            if (iso_head.mol.title == "NOT_FOUND"){
                iso_head_skip = true;
            } else {
                std::stringstream newName;
                newName << "h_"<<iso_head.mol.energy;
                iso_head.mol.energy = newName.str();
            }
        }
    
        //Keep track on what gets picked 
        if ( iso_head.aps.size() == 1 ) {
            num_sid_h++;
        }else if ( iso_head.aps.size() == 2 ) {
            num_lnk_h++;
        }else{
            num_scf_h++;
        }

        // For each attachement point on that fragment
        for ( int m=0; m<iso_head.aps.size(); m++ ) {

            float degrees_turned_in_rad = 0.0; 
            Fragment copy_frag;

            // Check to see if the bond order is compatible
            // START ISOHEAD SAMPLING
            if ( !compare_dummy_bonds( layer_frag, j, iso_head, m ) ){
          
                continue;
            }
            // check growing_ref to see if it gets bigger
            int growing_ref_size_before = growing_ref.size();     

            // assign values for later use for the attachment of 
            // isosteres
            dummy_atom_head = iso_head.aps[m].dummy_atom;
            heavy_atom_tar = layer_frag.aps[j].heavy_atom;

            // If it is compatible, combine the two fragments to make a new, larger fragment
            std::pair<std::vector<std::vector<double>>,Fragment>
                combine_frag_info = combine_fragments_iso( layer_frag, j, iso_head, m );

            
            // check if the matrix is a 3x3
            if (combine_frag_info.first.size() == 3){

                if(combine_frag_info.first[0].size() != 3 ||
                   combine_frag_info.first[1].size() != 3 ||
                   combine_frag_info.first[2].size() != 3)
                {
                    std::cout <<"rot_mat in isofraglib is not a 3x3. exiting.." 
                            << std::endl;
                    exit(1);
                } else {

                    //if all checks out, assign the working_mat
                    working_mat[0][0] = combine_frag_info.first[0][0]; 
                    working_mat[0][1] = combine_frag_info.first[0][1]; 
                    working_mat[0][2] = combine_frag_info.first[0][2];

                    working_mat[1][0] = combine_frag_info.first[1][0]; 
                    working_mat[1][1] = combine_frag_info.first[1][1]; 
                    working_mat[1][2] = combine_frag_info.first[1][2];

                    working_mat[2][0] = combine_frag_info.first[2][0]; 
                    working_mat[2][1] = combine_frag_info.first[2][1]; 
                    working_mat[2][2] = combine_frag_info.first[2][2];
                }
            } else {
                std::cout <<"rot_mat in isofraglib is not a 3x3. exiting.." 
                        << std::endl;
                exit(1);
            }

            Fragment oriented_iso_head = combine_frag_info.second;
 
            combine_frag_info.first.clear();


            ////Grab the all torsions candidates
            std::pair<std::vector<float>,std::vector<Fragment> >
                   torsion_candidates = attach_and_get_torsions(layer_frag, j,
                                                                oriented_iso_head, m,
                                                                score,  simplex, typer);

            //Frgment final_iso_head; 
            float    final_angle = 0.0; 
            Fragment final_combined_head;
            bool     min_or_nah = false;
            bool     prep_or_nah = false;


            // If no torsion or dn_sample_torsion is false ( set in attach_and_get_torsions)
            if (torsion_candidates.second.size() == 0){
                final_angle      = 0.0;
                final_combined_head =
                    attach_isosteres(layer_frag, j, oriented_iso_head,m);
                min_or_nah = true;
                prep_or_nah = true;

            } else {
                std::pair<float,Fragment> best_iso_head_pair =
                    get_best_iso_head_tors(layer_frag, j ,torsion_candidates, m,
                                           score,  simplex, typer); 

                // If there are no viable candidate.
                if (best_iso_head_pair.first == -9999) {

                    final_angle      = 0.0;
                    final_combined_head =
                        attach_isosteres(layer_frag, j, oriented_iso_head,m);
                    min_or_nah = true;
                    prep_or_nah = true;
            
                } else {

                    final_angle         = best_iso_head_pair.first;
                    final_combined_head = best_iso_head_pair.second;
                    min_or_nah  = false;
                    prep_or_nah = false;
                }
            }
 
            torsion_candidates.first.clear();
            torsion_candidates.second.clear();

            if (dn_iso_print_out.compare("yes")==0 && vec_iso_frags.second.size() >1) {

                DOCKMol final_combined_head_mol_copy;
                copy_molecule(final_combined_head_mol_copy,final_combined_head.mol);

                before_fout_isofrag << "\n" <<
                    DELIMITER << std::setw(STRING_WIDTH) <<
                    "Name:" << std::setw(FLOAT_WIDTH) <<final_combined_head.mol.energy <<
                    endl <<

                    DELIMITER << std::setw(STRING_WIDTH) <<
                    "Head_name:" << std::setw(FLOAT_WIDTH) << iso_head.mol.energy <<
                    endl <<

                    DELIMITER << std::setw(STRING_WIDTH) <<
                    "att_point:" << std::setw(FLOAT_WIDTH) << j <<
                    endl <<

                    DELIMITER << std::setw(STRING_WIDTH) <<
                    "Iso_skip:" << std::setw(FLOAT_WIDTH) << dn_iso_skip <<
                    endl <<

                    DELIMITER << std::setw(STRING_WIDTH) <<
                    "Head_accepted:" << std::setw(FLOAT_WIDTH) << boolean_write <<
                    endl;

                //activate mol and write out mol2
                activate_mol(final_combined_head_mol_copy);
                Write_Mol2(final_combined_head_mol_copy,before_fout_isofrag);
                final_combined_head_mol_copy.clear_molecule();

            }
 
            // If the torenv is turned on...
            if (dn_use_torenv_table){
                bool valid_torenv_bool = false;
                // First check to see if the newest torsion environment is allowable
                // (note that there may be some dummy atoms included in the environment, and
                // even though it passes a check now, it may not pass the check later)
                if(dn_use_roulette){
                    valid_torenv_bool = roulette_valid_torenv( final_combined_head );
                } else {
                    valid_torenv_bool = valid_torenv( final_combined_head );
                }

                if (valid_torenv_bool){
                    // Second check to see whether all of the old torsion environments are
                    // still compatible. This is necessary because of the dummy as wild card
                    // thing
                    if (valid_torenv_multi( final_combined_head )){
                        molecule_counter++;
                        num_att_layer++;
                        num_att_root++;

    	        	just_minimize(final_combined_head, growing_ref, score, simplex, typer,
                        prep_or_nah, min_or_nah);
                        
                        boolean_write="TRUE";
                    } else {

                        // if the head torsions check is not compatible, skip the isosteres too 
                        if ( dn_iso_skip.compare("yes") == 0 ) {
                            iso_head_skip = true;
                        }
                        boolean_write="FALSE";
                    }
                } else {
                    if ( dn_iso_skip.compare("yes") == 0 ) {
                        iso_head_skip = true;
                    }
                    boolean_write="FALSE";
                }

            // If the torenv is not turned on, just skip straight to sampling torsions
            } else {
                molecule_counter++;
                num_att_layer++;
                num_att_root++;

    	        just_minimize( final_combined_head, growing_ref, score, simplex, typer, 
                    prep_or_nah, min_or_nah );
                boolean_write="TRUE";
            }

            if ( dn_iso_print_out.compare("yes")==0 && vec_iso_frags.second.size() > 1 ) {

                fout_isofrag << "\n" <<
                    DELIMITER << std::setw(STRING_WIDTH) <<
                    "Name:" << std::setw(FLOAT_WIDTH) <<final_combined_head.mol.energy <<
                    endl <<

                    DELIMITER << std::setw(STRING_WIDTH) <<
                    "Head_name:" << std::setw(FLOAT_WIDTH) << iso_head.mol.energy <<
                    endl <<

                    DELIMITER << std::setw(STRING_WIDTH) <<
                    "att_point:" << std::setw(FLOAT_WIDTH) << j <<
                    endl <<

                    DELIMITER << std::setw(STRING_WIDTH) <<
                    "Iso_skip:" << std::setw(FLOAT_WIDTH) << dn_iso_skip <<
                    endl <<

                    DELIMITER << std::setw(STRING_WIDTH) <<
                    "Head_accepted:" << std::setw(FLOAT_WIDTH) << boolean_write <<
                    endl;

                //activate mol and write out mol2
                activate_mol( final_combined_head.mol );
                Write_Mol2( final_combined_head.mol,fout_isofrag );
            }

            // If we got this far, and if the fragment we are looking at is a scaffold,
            // and if growing_ref increased in size, then we added a scaffold this layer
            if ( ( iso_head.aps.size() > 2 ) && 
                 ( growing_ref.size() > growing_ref_size_before ) ){ 

                layer_frag.scaffolds_this_layer += 1;

            }

            // If iso_head did not attach, please skip
            if ( iso_head_skip ) { continue; }
  
	    int c_for_topranked =    0;
            int iso_passes      =    0;
            int iso_it          =    0;
            rate_iso_att        = -1.0;
            int num_of_iso_att  =    0;
            std::vector <int> iso_picked {};

	    //This is the looping through of the all of the 
            //returned "isosteres" for alignment
            //START isoswapping sampling
            while (c_for_topranked < dn_num_iso_picks){

                // if the num of att per root exceeds parameter, break
                if ( num_att_root > dn_max_sccful_att_per_root ) { break; }

                Fragment isostere;

                //Get rand fragments
                if (dn_iso_pick_meth.compare("rand") == 0){
                    // The passes int counts every time we either keep or 
                    // skip a fragment from the libary, so
                    // once its size reaches the fraglib.size(), we have looked at everything
                    if ( iso_passes >= vec_iso_frags.second.size() ){ break; }

                    // Choose a random starting position, fraglib[iso_it]
                    // (rand has been seeded if and only if the minimizer is turned on)

                    rank_num = rand() % vec_iso_frags.second.size();
                    // verbose dn frag selection random seed output
                    if (verbose) cout << "De novo fragment selection: " << rank_num << endl; 

                    // Keep track of things we try so we don't try the same thing twice
                    if (iso_picked.size() > 0){
                        if (iso_picked.size() >= vec_iso_frags.second.size()){
                           iso_picked.clear();
                           break; 
                        }
                        for (int p=0; p<iso_picked.size(); p++){
                            if (rank_num == iso_picked[p]){ continue; }
                        }
                    }

                    iso_picked.push_back(rank_num);
                    isostere = vec_iso_frags.second[rank_num];

                    iso_it++;

                //Get the top rank then work your self up
                } else if (dn_iso_pick_meth.compare("top_rank") == 0){
                    if (iso_it == vec_iso_frags.second.size()){ 
                        break;
                    }
                    rank_num = iso_it;
                    // verbose dn frag selection random seed output
                    isostere = vec_iso_frags.second[iso_it];
                    //to just iterate right after you get the next isostere.
                    iso_it++;

                //Get the worst rank then work your self up
                } else if (dn_iso_pick_meth.compare("worst_rank") == 0){
                    if (iso_it == vec_iso_frags.second.size()){
                        break;
                    }
                    rank_num = vec_iso_frags.second.size() - 1 - iso_it;
                    // verbose dn frag selection random seed output
                    isostere = vec_iso_frags.second[vec_iso_frags.second.size() - 1 - iso_it];
                    //to just iterate right after you get the next isostere.
                    iso_it++;
                } else {
                    std::cout << " either incorrect iso pick method selected or none selected. exiting..."  <<
                        std::endl;
                    exit(1);
                }

                //FILTERING
                // (1) Would create too many scaffolds per layer
                if ( isostere.aps.size() > 2 &&
                     layer_frag.scaffolds_this_layer >= dn_max_scaffolds_per_layer ){ 
                    iso_passes++;
                    continue; 
                }

                // (2) Would give too many current attachment points
                int temp_aps = isostere.aps.size() + layer_frag.aps.size() - 2;
                if ( temp_aps > dn_max_current_aps ) { 
                    iso_passes++;
                    continue; 
                }

                // (3) Rotatable bond cutoff
                // BCF assuming rot_bon does not take into account attachment points
                // will have at least one new rot bond per attachment point, on each fragment
                // except when combined two aps will become one rot bond.
                // Floor for remaining rot_bonds if all remaining aps are capped with sidechains.
                int temp_rbs = layer_frag.mol.rot_bonds + temp_aps;
                if ( temp_rbs > dn_constraint_rot_bon ) { 
                    iso_passes++;
                    continue; 
                }

                // (5) Formal Charge Cutoff
                int temp_fc = isostere.mol.formal_charge + layer_frag.mol.formal_charge;
                if ( fabs ( float( temp_fc ) ) > ( fabs ( float (dn_constraint_formal_charge ) )+0.1) ) { 
                    iso_passes++;
                    continue; 
                }

	        bool compatible = false;

                if ( isostere.aps.size() == 1 ) {
                    num_sid++;
                } else if ( isostere.aps.size() == 2 ) {
                    num_lnk++;
                } else {
                    num_scf++;
                }

                //if dummy_atom_head and heavy_atom_tar is not messed up, keep going
                // and sample the isosteres
                if ( dummy_atom_head < 0 && heavy_atom_tar < 0 ) {
                    continue;
                }

                int du_tail_atom_num = -1;
                int hvy_tail_atom_num = -1;
                Fragment isofrag;
                
                //// get the du atom num that is associated with 
                //// the dummy atom of the iso frag
                for ( unsigned int z =0; z < isostere.best_tri.three_pairs.size(); z++ ) {

                    // return the atom num of the iso ref frag
                    // that is associated with the current working
                    // test frag that is trying to attach with the layer frag
                    int ref_atom_num = isostere.best_tri.three_pairs[z].get_ref();
                    // If the name of the du atom of the head frag is
                    // the same as the iso reference frag, 
                    // get the test frag atom num that is aligned with 
                    // the iso head du.
                    if ( dummy_atom_head == ref_atom_num ) {
                        du_tail_atom_num = isostere.best_tri.three_pairs[z].get_test();
                        
                        std::vector<int> tmp_nbrs = isostere.mol.get_atom_neighbors(du_tail_atom_num);
                        if (tmp_nbrs.size() != 1){
                            std::cout << "neighboring atom of a selected du atom is not size of 1. exiting...";
                            exit(0);
                        }
                        hvy_tail_atom_num = tmp_nbrs[0];
                        tmp_nbrs.clear();
                        break;
                    }
                }
                
                if ( du_tail_atom_num < 0 || hvy_tail_atom_num < 0 ) {
                    if ( verbose ) {
                        std::cout << "Isohead " << iso_head.mol.energy << " isostere " << isostere.mol.energy << std::endl;
                        std::cout << "du_tail_atom_num or hvy_tail_atom_num is -1" << std::endl;
                        std::cout << "It could be your mol2 files do no contain " <<
                                     "unique atom names" << std::endl;
                        std::cout << "Or, the du atom num of iso_heads that is chosen is not part of " << std::endl;
                        std::cout << " the three atom pairs for candidate isostere" << std::endl;
                    }
                    iso_passes++;
                    c_for_topranked++;
                    break;
                }

                // Calculate the desired bond length and remember as 'new_rad'
                float new_rad = calc_cov_radius(layer_frag.mol.atom_types[heavy_atom_tar]) +
                                calc_cov_radius(isostere.mol.atom_types[hvy_tail_atom_num]);

                // Calculate the x-y-z components of the current frag2 bond vector (bond_vec)
                DOCKVector bond_vec;
                bond_vec.x = isostere.mol.x[hvy_tail_atom_num] - 
                             isostere.mol.x[du_tail_atom_num];

                bond_vec.y = isostere.mol.y[hvy_tail_atom_num] - 
                             isostere.mol.y[du_tail_atom_num];

                bond_vec.z = isostere.mol.z[hvy_tail_atom_num] - 
                             isostere.mol.z[du_tail_atom_num];

                // Normalize the bond vector then multiply each component by new_rad so that it is the desired
                // length
                bond_vec = bond_vec.normalize_vector();
                bond_vec.x *= new_rad;
                bond_vec.y *= new_rad;
                bond_vec.z *= new_rad;

                // Change the coordinates of the frag2 dummy atom so that the bond length is correct
                isostere.mol.x[du_tail_atom_num] = 
                    isostere.mol.x[hvy_tail_atom_num] - bond_vec.x;

                isostere.mol.y[du_tail_atom_num] = 
                    isostere.mol.y[hvy_tail_atom_num] - bond_vec.y;

                isostere.mol.z[du_tail_atom_num] = 
                    isostere.mol.z[hvy_tail_atom_num] - bond_vec.z;

                // Figure out what translation is required to move the dummy atom to the origin
                // for the iso test frags
                DOCKVector trans1;
                trans1.x = -isostere.mol.x[du_tail_atom_num];
                trans1.y = -isostere.mol.y[du_tail_atom_num];
                trans1.z = -isostere.mol.z[du_tail_atom_num];

                // Use the dockmol function to translate the fragment so the dummy atom is at the origin
                isostere.mol.translate_mol(trans1);

                //apply the previous rotation matrix
                isostere.mol.rotate_mol( working_mat );

                ////// rotate the molecule at x axis and return back in place.
                rotate_on_x_axis(isostere.mol , du_tail_atom_num ,
                                 isostere.mol.get_atom_neighbors(du_tail_atom_num)[0],
                                 final_angle);

                // translate to the heavy atom neightboring du atom in 
                // the root layer fragment
                DOCKVector trans2;
                trans2.x = layer_frag.mol.x[heavy_atom_tar];
                trans2.y = layer_frag.mol.y[heavy_atom_tar];
                trans2.z = layer_frag.mol.z[heavy_atom_tar];

                isostere.mol.translate_mol(trans2);
    
                //get the heavy atom num that is neighboring to the du for the isostere
                // Then attach isofrag
                for (unsigned int y = 0; y < isostere.aps.size(); y++){
                    if (isostere.aps[y].dummy_atom == du_tail_atom_num){
                        //attach layer_frag with isosteres 
    	                isofrag = attach_isosteres(layer_frag, j,
    	                    		           isostere, y);
                        break;
                    } 

                    if ( y == isostere.aps.size()-1 ){
	    	        std::cout << "Exiting. can't get the heavy atom from du tail" << std::endl;
	    	        exit(1);
                    } 
                   
                }

                if (dn_iso_print_out.compare("yes")==0) {
                    DOCKMol isofrag_mol_copy;
                    copy_molecule(isofrag_mol_copy,isofrag.mol);
    
                    before_fout_isofrag << "\n" <<
                        DELIMITER << std::setw(STRING_WIDTH) <<
                        "Name:" << std::setw(FLOAT_WIDTH) <<isofrag.mol.energy <<
                        endl <<

                        DELIMITER << std::setw(STRING_WIDTH) <<
                        "Head_name:" << std::setw(FLOAT_WIDTH) << iso_head.mol.energy <<
                        endl <<

                        DELIMITER << std::setw(STRING_WIDTH) <<
                        "Tail_name:" << std::setw(FLOAT_WIDTH) << isostere.mol.energy <<
                        endl <<

                        DELIMITER << std::setw(STRING_WIDTH) <<
                        "Tail_size:" << std::setw(FLOAT_WIDTH) << vec_iso_frags.second.size()<<
                        endl <<

                        DELIMITER << std::setw(STRING_WIDTH) <<
                        "rank_num:" << std::setw(FLOAT_WIDTH) << rank_num <<
                        endl <<

                        DELIMITER << std::setw(STRING_WIDTH) <<
                        "att_point:" << std::setw(FLOAT_WIDTH) << j <<
                        endl <<

                        DELIMITER << std::setw(STRING_WIDTH) <<
                        "Method_of_sel:" << std::setw(FLOAT_WIDTH) << dn_iso_pick_meth <<
                        endl; 

                    activate_mol(isofrag_mol_copy);
                    Write_Mol2(isofrag_mol_copy,before_fout_isofrag);
                    isofrag_mol_copy.clear_molecule();

                }

                //Check torsion env or not...
    	        if (dn_use_torenv_table){
                        bool valid_torenv_bool = false;
                        // First check to see if the newest torsion environment is allowable
                        // (note that there may be some dummy atoms included in the environment, and
                        // even though it passes a check now, it may not pass the check later)

                        if(dn_use_roulette){
                            valid_torenv_bool=roulette_valid_torenv(isofrag);
                        } else {
                            valid_torenv_bool=valid_torenv(isofrag);
                        }

                        if (valid_torenv_bool){
                            // Second check to see whether all of the old torsion environments are
                            // still compatible. This is necessary because of the dummy as wild card
                            // thing
                            
                            if (valid_torenv_multi(isofrag)){
                                molecule_counter++;
                                num_att_layer++;
                                num_att_root++;

                                // If both checks pass, then this is a valid molecule. Sample torsions.
                                just_minimize(isofrag, growing_ref, score, simplex, typer, true, true );
                                boolean_write="TRUE";
                                num_of_iso_att++;
                            } else {
                                boolean_write="FALSE";
                            }
                        } else {
                            boolean_write="FALSE";
                        }
    	        } else {
                    molecule_counter++;
                    num_att_layer++;
                    num_att_root++;

    	            just_minimize(isofrag, growing_ref, score, simplex, typer, true, true );
                    boolean_write="TRUE";
                    num_of_iso_att++;
    	        }
//////////////////////
                //if (c_for_topranked == dn_num_iso_picks-1){
                //    if (c_for_topranked == vec_iso_frags.second.size()-1){ 
                //        rate_iso_att =  (float)num_of_iso_att / 
                //                              vec_iso_frags.second.size();
                //            std::cout << "non_filt_rates_successful_att: " << rate_iso_att << std::endl;

                //    } else if (c_for_topranked == (num_sid+num_lnk+num_scf)-1){
                //        rate_iso_att =  (float)num_of_iso_att / 
                //                              (num_sid+num_lnk+num_scf);
                //            std::cout << "rates_successful_att: " << rate_iso_att << std::endl;
                //    }
                //} else if (iso_it == vec_iso_frags.second.size()){
                //    rate_iso_att =  (float)num_of_iso_att / 
                //                          (vec_iso_frags.second.size());
                //        std::cout << "trunc_rates_successful_att: " << rate_iso_att << std::endl;
                //}

//////////////////////
                if (dn_iso_print_out.compare("yes")==0 || Parameter_Reader::verbosity_level() > 1 ) {
                    if (c_for_topranked == dn_num_iso_picks-1){
                        if (c_for_topranked == vec_iso_frags.second.size()-1){ 
                            rate_iso_att =  (float)num_of_iso_att / 
                                                  vec_iso_frags.second.size();
                            if (Parameter_Reader::verbosity_level() > 1) {
                               std::cout << "non_filt_rates_successful_att: " << rate_iso_att << std::endl;
                            }

                        } else if (c_for_topranked == (num_sid+num_lnk+num_scf)-1){
                            rate_iso_att =  (float)num_of_iso_att / 
                                                  (num_sid+num_lnk+num_scf);

                            if (Parameter_Reader::verbosity_level() > 1) {
                               std::cout << "rates_successful_att: " << rate_iso_att << std::endl;
                            }
                        }
                    } else if (iso_it == vec_iso_frags.second.size()){
                        rate_iso_att =  (float)num_of_iso_att / 
                                              (vec_iso_frags.second.size());
                        if (Parameter_Reader::verbosity_level() > 1) {
                            std::cout << "trunc_rates_successful_att: " << rate_iso_att << std::endl;
                        }
                    }

                    if (dn_iso_print_out.compare("yes")==0) {
                        fout_isofrag << "\n" <<
                            DELIMITER << std::setw(STRING_WIDTH) <<
                            "Name:" << std::setw(FLOAT_WIDTH) <<isofrag.mol.energy <<
                            endl <<

                            DELIMITER << std::setw(STRING_WIDTH) <<
                            "Head_name:" << std::setw(FLOAT_WIDTH) << iso_head.mol.energy <<
                            endl <<

                            DELIMITER << std::setw(STRING_WIDTH) <<
                            "Tail_name:" << std::setw(FLOAT_WIDTH) << isostere.mol.energy <<
                            endl <<

                            DELIMITER << std::setw(STRING_WIDTH) <<
                            "Tail_size:" << std::setw(FLOAT_WIDTH) << vec_iso_frags.second.size()<<
                            endl <<

                            DELIMITER << std::setw(STRING_WIDTH) <<
                            "rank_num:" << std::setw(FLOAT_WIDTH) << rank_num <<
                            endl <<

                            DELIMITER << std::setw(STRING_WIDTH) <<
                            "att_point:" << std::setw(FLOAT_WIDTH) << j <<
                            endl <<

                            DELIMITER << std::setw(STRING_WIDTH) <<
                            "Tail_accepted:" << std::setw(FLOAT_WIDTH) << boolean_write <<
                            endl <<

                            DELIMITER << std::setw(STRING_WIDTH) <<
                            "Tail_acc_rate:" << std::setw(FLOAT_WIDTH) << rate_iso_att <<
                            endl <<

                            DELIMITER << std::setw(STRING_WIDTH) <<
                            "Method_of_sel:" << std::setw(FLOAT_WIDTH) << dn_iso_pick_meth <<
                            endl; 

                        activate_mol(isofrag.mol);
                        Write_Mol2(isofrag.mol,fout_isofrag );
                    }
                }
    	        c_for_topranked++;
                iso_passes++;
            }//END isoswapping sampling

            num_sid =0;
            num_lnk =0;
            num_scf =0;
            rate_iso_att=-1;

        }//END ISOHEAD APS SAMPLING

        if ( choice + 1 >= dn_num_rand_head_picks && 
            ( growing_ref.size() == growing_ref_size_before_aps ) ){
            passes++;
            continue;
        }

        vec_iso_frags.first.clear();
        vec_iso_frags.second.clear();
        choice++;
        passes++;
    } // ... end While N choices

    dn_iso_write_it++;
    picked.clear();

    return;
}
// +++++++++++++++++++++++++++++++++++++++++
// Sample the contents of a given fragment library following the FragGraph
void
DN_Build::sample_fraglib_rand( Fragment & layer_frag, int j, vector <Fragment> & fraglib, vector <Fragment> & growing_ref,
                                Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & typer, 
                                vector<int> & available_fragment_indices, vector<int> & selected_fragments, 
                                bool last_layer )
{
    //cout << "-sample_fraglib_rand-" << endl;
    // If the fragment library is empty, return without doing anything
    if (fraglib.size() == 0 || num_att_root > dn_max_sccful_att_per_root ) { return; }
    //BCF Calc MW RB and FC for pruning later
    calc_mol_wt(layer_frag.mol);
    calc_rot_bonds(layer_frag.mol);
    calc_formal_charge(layer_frag.mol);

    vector< pair <Fragment, int> > initial_ensemble;
    //select a number of fragments from the fragment library...
    select_frags_from_fraglib(scaf_link_sid, available_fragment_indices, selected_fragments,
                              dn_num_random_picks, true, last_layer);

    //... and attach all of them to layer_frag. 
    attach_selected_frags( layer_frag, j, scaf_link_sid, selected_fragments, initial_ensemble);

    //If no molecules were successfully made, perform more rounds of selections
    //using a smaller selection pool (first pick == false in the select_frags function)
    //until something is made or all fragments exhausted.
    int number_of_picks = dn_num_random_picks;
    while (initial_ensemble.size() == 0) {
        if (available_fragment_indices.size() < 1 || num_att_root > dn_max_sccful_att_per_root ){
            break;
        }
        if (available_fragment_indices.size() < dn_num_random_picks) {
            number_of_picks = available_fragment_indices.size();
        } else {
            number_of_picks = dn_num_random_picks;
        }
        
        select_frags_from_fraglib(scaf_link_sid, available_fragment_indices, selected_fragments,
                                number_of_picks, false, last_layer);
        attach_selected_frags( layer_frag, j, scaf_link_sid, selected_fragments, initial_ensemble);
    }
    //From here, a vector < pair < Fragment, int > > that pairs all of the growing fragments with
    //their attempted fraglib index has been generated. These have already been pruned for overlaps
    //if enabled.
    molecule_counter += initial_ensemble.size();
    //num_att_layer += initial_ensemble.size();
    //Iterate through this entire vector, sampling torsions along the way
    vector< Fragment > ensemble_collector;
    vector< pair <Fragment, int> >::iterator pos = initial_ensemble.begin();
    while (pos != initial_ensemble.end()){
        vector< Fragment > torsions_collector;
        //1. INITIAL TORSION SAMPLING
        //cout << "Number of fragments before torsional sampling: " << torsions_collector.size() << endl;
        //sample_minimized_torsions((*pos).first, torsions_collector, score, simplex, typer);
        sample_minimized_torsions((*pos).first, torsions_collector, score, simplex, typer);
        //cout << "Number of fragments after torsional sampling: " << torsions_collector.size() << endl;
        if (torsions_collector.size() == 0) {
            pos++;
            continue;
        }

        //if we end up with no torsions, move to next fragment in the queue

        //2. PRUNE SAMPLED TORSIONS BY HRMSD
        //If enabled, sort and prune the torsions to account for any overlap
        //Activation is done in the sampling function - not necessary here
        if (dn_prune_individual_torsions && torsions_collector.size() > 1) {
            prune_h_rmsd_and_mw(torsions_collector);  
        }

        //cout << "Number of fragments after torsion pruning:  " << torsions_collector.size() << endl;
        //all of these torsions are valid for the moment - push them to the final ensemble for this round
        for (int mol = 0; mol < torsions_collector.size(); mol++){
            ensemble_collector.push_back(torsions_collector[mol]);
        }

        //clear torsions for next molecule, and increment iterator position
        //torsions_collector.clear();
        pos++;
    }

    //Prune the combined ensemble's output
    if (dn_prune_combined_torsions){
        prune_h_rmsd_and_mw(ensemble_collector);
    }

    //Collect the sampled molecules back onto growing
    for (int mol=0; mol<ensemble_collector.size(); mol++){
        growing_ref.push_back(ensemble_collector[mol]);
    }
    ensemble_collector.clear();
    return;
}

// +++++++++++++++++++++++++++++++++++++++++
// Sample the contents of a given fragment library following the FragGraph
void
DN_Build::sample_fraglib_graph( Fragment & layer_frag, int j, vector <Fragment> & fraglib,
                                vector <FragGraph> & fraggraph_vec, vector <Fragment> & growing_ref,
                                Master_Score & score, Simplex_Minimizer & simplex,
                                AMBER_TYPER & typer, vector<int> & available_fragment_indices, vector<int> & selected_fragments, 
                                bool last_layer )
{
    //cout << "-SAMPLE_FRAGLIB_GRAPH-" << endl;
    // If the fragment library is empty, return without doing anything
    if (fraglib.size() == 0 || num_att_root > dn_max_sccful_att_per_root) { return; }

    //BCF Calc MW RB and FC for pruning later
    calc_mol_wt(layer_frag.mol);
    calc_rot_bonds(layer_frag.mol);
    calc_formal_charge(layer_frag.mol);

    vector< pair <Fragment, int> > initial_ensemble;
    //select a number of fragments from the fragment library...
    select_frags_from_fraglib(scaf_link_sid, available_fragment_indices, selected_fragments,
                              dn_graph_max_picks, true, last_layer);

    //... and attach all of them to layer_frag. 
    attach_selected_frags( layer_frag, j, scaf_link_sid, selected_fragments, initial_ensemble);
    int number_of_picks = dn_graph_max_picks;
    //If no molecules were successfully made, perform more rounds of selections
    //using a smaller selection pool (first pick == false in the select_frags function)
    //until something is made or all fragments exhausted.
    while (initial_ensemble.size() == 0) {
        if (available_fragment_indices.size() < 1 || num_att_root > dn_max_sccful_att_per_root){
            break;
        }
        if (available_fragment_indices.size() < dn_graph_max_picks) {
            number_of_picks = available_fragment_indices.size();
        } else {
            number_of_picks = dn_graph_max_picks;
        }
        select_frags_from_fraglib(scaf_link_sid, available_fragment_indices, selected_fragments,
                              number_of_picks, false, last_layer);
        attach_selected_frags( layer_frag, j, scaf_link_sid, selected_fragments, initial_ensemble);
    }
    molecule_counter += initial_ensemble.size();
    //num_att_layer += initial_ensemble.size();
  
    //From here, a vector < pair < Fragment, int > > that pairs all of the growing fragments with
    //their attempted fraglib index has been generated. These have already been pruned for overlaps.

    //This is done as its own initial step to prune out initial attachments potentially having
    //geometric overlap due to fragment symmetry, which can lead to heavy duplicates

    //Iterate through this entire vector, with each position acting as a starting point
    //for the graph method
    vector< Fragment > torsions_collector;
    vector< Fragment > ensemble_collector;
    int max_tries = pow(float(dn_graph_breadth), dn_graph_depth);
    vector< pair <Fragment, int> >::iterator pos = initial_ensemble.begin();
    vector < pair <int, float> > sampling_queue;
    int test_counter=0;
    //cout << initial_ensemble.size() << endl;
    while (pos != initial_ensemble.end()){
        test_counter+=1;
        //assign the baseline as the initial root with no attachments
        float compare_score = layer_frag.mol.current_score;
        bool improvement_flag = false;
        bool first_selection = true;
        int tries = 0;

        sampling_queue.push_back(make_pair((*pos).second, layer_frag.mol.current_score));
        while (sampling_queue.size() > 0 && tries < max_tries){
            tries++;
            pair <int, float> this_dnrn;
            this_dnrn = sampling_queue.front();
            sampling_queue.erase(sampling_queue.begin());
            // Some things can be added to the queue more than once - e.g. they are already in the
            // queue, but not yet visited, so they are added again.
            if (fraggraph_vec[this_dnrn.first].visited == true){ continue; }

            // Mark this position as visited in the FragGraph, we are now "trying" a fragment

            fraggraph_vec[this_dnrn.first].visited = true;
            vector<Fragment> attachment_collector;
            if (!first_selection) {
                attach_single_frag(layer_frag, j, scaf_link_sid, this_dnrn.first, attachment_collector);
            } else {
                attachment_collector.push_back((*pos).first);
                first_selection=false;
            }
            //cout << "Number of attachments after initial selection " << attachment_collector.size() << endl;
            //1. INITIAL TORSION SAMPLING
            //cout << "Number of fragments before torsional sampling: " << torsions_collector.size() << endl;
            for(int mol = 0; mol < attachment_collector.size(); mol++) {
                sample_minimized_torsions(attachment_collector[mol], torsions_collector, score, simplex, typer);
            }
            //cout << "Number of fragments after torsional sampling: " << torsions_collector.size() << endl;
            //if we end up with no torsions, move to next fragment in the queue
            if (torsions_collector.size() == 0) {
                continue;
            }

            //2. PRUNE SAMPLED TORSIONS BY HRMSD
            //If enabled, sort and prune the torsions to account for any overlap
            //Activation is done in the sampling function - not necessary here
            if (dn_prune_individual_torsions && torsions_collector.size() > 1) {
                prune_h_rmsd_and_mw(torsions_collector);  
            }

            //cout << "Number of fragments after torsion pruning:  " << torsions_collector.size() << endl;
            //all of these torsions are valid for the moment - push them to the final ensemble for this round
            for (int mol = 0; mol < torsions_collector.size(); mol++){
                ensemble_collector.push_back(torsions_collector[mol]);
            }

            //3 GRAPH TRAVERSAL
            //get a score to compare against == the best score of initially sampled torsions
            float best_tors_score = dn_pruning_conformer_score_cutoff; // highest possible score

            //Because vector is sorted by score, the best score is the first torsion in the list
            //if something better was made, replace the score_cutoff value
            if (torsions_collector[0].mol.current_score < dn_pruning_conformer_score_cutoff){
                best_tors_score = torsions_collector[0].mol.current_score;
            }
            torsions_collector.clear();

            float temp1 = exp( (compare_score - best_tors_score) /
                                dn_graph_temperature );
            float temp2 = ((double) rand() / RAND_MAX);
            // Check the score of what we just kept, see if it is an improvement over
            // the previous best score according to {  exp(-(Ef-Ei)/T) > Rand }
            if (temp1 > temp2) {
                improvement_flag = true;
                    // Do this second check because you can be trying one new fragment in 
                    // multiple ways, and you just want to keep the best score
                    if (best_tors_score < compare_score){
                        compare_score = best_tors_score;
                    }

            }
            // If this fragment improved the score over the previous score, add its neighbors
            // to the sampling_queue (only the unvisited neighbors)
            if (improvement_flag){
                int temp_breadth = dn_graph_breadth;

                // If this is a very small frag graph or a very large breadth, reset the breadt
                // to be the size of the graph
                if (fraggraph_vec.size() <= temp_breadth){ temp_breadth = ( fraggraph_vec.size() - 1 ); }

                // For every child that will be sampled at this depth
                for (int n=1; n<temp_breadth+1; n++){

                    // BCF check to prevent seg fault
                    if (fraggraph_vec[this_dnrn.first].rankvec.size() == n) {break;}

                    // Make sure that dn_max_current_aps won't be exceeded by adding this fragment
                    int temp_aps = fraglib[ fraggraph_vec[this_dnrn.first].rankvec[n] ].aps.size() +
                               layer_frag.aps.size() - 2;
                    if (temp_aps > dn_max_current_aps){
                        fraggraph_vec[ fraggraph_vec[this_dnrn.first].rankvec[n] ].visited = true;
                    }

                    // Only check it if it has not already been visited
                    if ( !fraggraph_vec[ fraggraph_vec[this_dnrn.first].rankvec[n] ].visited ){

                        // Make a record of its position and the current score, and it the frag
                        // to the sampling queue
                        pair <int, float> temp_pair = make_pair(fraggraph_vec[this_dnrn.first].rankvec[n],
                                                                compare_score);
                        sampling_queue.push_back( temp_pair );
                    }
                }
            }

            improvement_flag = false;
        }
        sampling_queue.clear();
        pos++;

    }//end for each fragment in initial ensemble
    //Refresh the fraggraph for next time
    for (int i=0; i<fraggraph_vec.size(); i++){
        fraggraph_vec[i].visited = false;
    }

    //Prune the combined ensemble
    if (dn_prune_combined_torsions){
        prune_h_rmsd_and_mw(ensemble_collector);
    }

    //Collect the sampled molecules back onto growing
    for (int mol=0; mol<ensemble_collector.size(); mol++){
        growing_ref.push_back(ensemble_collector[mol]);
    }
    ensemble_collector.clear();
    //cout << "Number of mols after GRAPHing " << growing_ref.size() << endl; 
    return;
}

// +++++++++++++++++++++++++++++++++++++++++
// Sample the contents of a given fragment library exhaustively
void
DN_Build::sample_fraglib_exhaustive( Fragment & layer_frag, int j, vector <Fragment> & fraglib,
                                     vector <Fragment> & growing_ref, Master_Score & score, 
                                     Simplex_Minimizer & simplex, AMBER_TYPER & typer, bool last_layer )
{

    // If the fragment library is empty, return without doing anything
    if (fraglib.size() == 0 || num_att_root > dn_max_sccful_att_per_root) { return; }

    // For every fragment in the fragment library (i.e. sidechains, linkers, scaffolds)
    for (int l=0; l<fraglib.size(); l++){
  
        if (num_att_root > dn_max_sccful_att_per_root) { break; }

        // (1) Would create too many scaffolds per layer
        if ( fraglib[l].aps.size() > 2 &&
             layer_frag.scaffolds_this_layer >= dn_max_scaffolds_per_layer ){ continue; }


        // (2) Would give too many current attachment points
        int temp_aps = fraglib[l].aps.size() + layer_frag.aps.size() - 2;
        if (temp_aps > dn_max_current_aps) { continue; }

        // (3) Rotatable bond cutoff
        // BCF assuming rot_bon does not take into account attachment points
        // will have at least one new rot bond per attachment point, on each fragment
        // except when combined two aps will become one rot bond.
        // Floor for remaining rot_bonds if all remaining aps are capped with sidechains.
        int temp_rbs = layer_frag.mol.rot_bonds + temp_aps;
        if (temp_rbs > dn_constraint_rot_bon) { continue; }

        // (4) Molecular weight cutoff (MW computed for layer frag at last step, MW precomputed for all frags)
        //int temp_mw = fraglib[l].mol.mol_wt + layer_frag.mol.mol_wt;
        //if (temp_mw > dn_constraint_mol_wt) { continue; }

        // (5) Formal Charge Cutoff
        int temp_fc = fraglib[l].mol.formal_charge + layer_frag.mol.formal_charge;
        if (fabs(float(temp_fc)) > (fabs(float(dn_constraint_formal_charge))+0.1)) { continue; }

        // For each attachement point on that fragment
        for (int m=0; m<fraglib[l].aps.size(); m++){
            
            // Check to see if the bond order is compatible
            if (compare_dummy_bonds(layer_frag, j, fraglib[l], m)){

                // If it is compatible, combine the two fragments to make a new, larger fragment
                Fragment frag_combined = combine_fragments(layer_frag, j, fraglib[l], m);

                // check growing_ref to see if it gets bigger
                int growing_ref_size_before = growing_ref.size();        

                // If the torenv is turned on...
                if (dn_use_torenv_table){
                    bool valid_torenv_bool = false;
                    // First check to see if the newest torsion environment is allowable
                    // (note that there may be some dummy atoms included in the environment, and
                    // even though it passes a check now, it may not pass the check later)
                    if(dn_use_roulette){
                        valid_torenv_bool=roulette_valid_torenv(frag_combined);
                    } else {
                        valid_torenv_bool=valid_torenv(frag_combined);
                    }

                    if (valid_torenv_bool){
                        // Second check to see whether all of the old torsion environments are
                        // still compatible. This is necessary because of the dummy as wild card
                        // thing

                        if (valid_torenv_multi(frag_combined)){

                            molecule_counter++;
                            num_att_layer++;
                            num_att_root++;

                            // If both checks pass, then this is a valid molecule. Sample torsions.
                            sample_minimized_torsions(frag_combined, growing_ref, score, simplex, typer);
                        }
                    }

                // If the torenv is not turned on, just skip straight to sampling torsions
                } else {
                    molecule_counter++;
                    num_att_layer++;
                    num_att_root++;
                    sample_minimized_torsions(frag_combined, growing_ref, score, simplex, typer);
                }

                // If we got this far, and if the fragment we are looking at is a scaffold,
                // and if growing_ref increased in size, then we added a scaffold this layer
                if ( (fraglib[l].aps.size() > 2) && 
                     (growing_ref.size() > growing_ref_size_before)){ 

                    layer_frag.scaffolds_this_layer += 1;
                }
            }
        }
    }

    return;
 
} // end DN_Build::sample_fraglib_exhaustive()
// +++++++++++++++++++++++++++++++++++++++++
// Sample the contents of a given fragment library randomly
void
DN_Build::sample_fraglib_matrix( Fragment & layer_frag, int j, vector <Fragment> & fraglib,
                                 vector <Fragment> & growing_ref, Master_Score & score,
                                 Simplex_Minimizer & simplex, AMBER_TYPER & typer, bool last_layer )
{

    // If the fragment library is empty for some reason, return without doing anything
    if (fraglib.size() == 0 || num_att_root > dn_max_sccful_att_per_root) { return; }
    // While N choices...
    int choice=0;
    int passes=0;
    int dnrn;

    std::vector <int> picked;
    //cout << "START FRAG NAME: " << layer_frag.mol.energy << endl;

    vector<pair<string, float> > possible_attachments;
    //cout << "ap " << layer_frag.aps[j].frag_name << endl;
    for (int i=0; i<ordered_frag_list.size(); ++i){
        if (ordered_frag_list[i] == layer_frag.aps[j].frag_name){
            for (int x=0; x<fragment_frequency_matrix[i].size(); ++x){
                if (fragment_frequency_matrix[i][x] > 0){
                    //cout << ordered_frag_list[i] << " to " << ordered_frag_list[x] << ": " << fragment_frequency_matrix[i][x] << endl;

                    possible_attachments.push_back(make_pair(ordered_frag_list[x],rel_fragment_frequency_matrix[i][x]));
                }
            }
        }
    }
    /*
    for (int entry=0; entry<possible_attachments.size();++entry){
        cout << "This fragment can be attached to " << possible_attachments[entry].first << " with frequency " << possible_attachments[entry].second << endl;
    }
    */
    vector< pair <Fragment, int> > initial_ensemble;
    vector<int> selected_fragments;
    //cout << "Num att " << possible_attachments.size() << endl;
    int sel_cap = 25;
    for(int attempts=0; attempts < possible_attachments.size(); attempts++){
        if(attempts > sel_cap) {
            break;
        }
        for (int entry=0; entry<fraglib.size(); ++entry){
            if (fraglib[entry].mol.energy == possible_attachments[attempts].first){
                //cout << "Trying fragment " << fraglib[entry].mol.energy << endl;
                selected_fragments.push_back(entry);
            }
        }
    }
    attach_selected_frags(layer_frag, j, fraglib, selected_fragments, initial_ensemble);
    //From here, a vector < pair < Fragment, int > > that pairs all of the growing fragments with
    //their attempted fraglib index has been generated. These have already been pruned for overlaps
    //if enabled.
    molecule_counter += initial_ensemble.size();
    //num_att_layer += initial_ensemble.size();
    //Iterate through this entire vector, sampling torsions along the way
    vector< Fragment > ensemble_collector;
    vector< pair <Fragment, int> >::iterator pos = initial_ensemble.begin();
    while (pos != initial_ensemble.end()){
        vector< Fragment > torsions_collector;
        //1. INITIAL TORSION SAMPLING
        //cout << "Number of fragments before torsional sampling: " << torsions_collector.size() << endl;
        //sample_minimized_torsions((*pos).first, torsions_collector, score, simplex, typer);
        sample_minimized_torsions((*pos).first, torsions_collector, score, simplex, typer);
        //cout << "Number of fragments after torsional sampling: " << torsions_collector.size() << endl;
        if (torsions_collector.size() == 0) {
            pos++;
            continue;
        }

        //if we end up with no torsions, move to next fragment in the queue

        //2. PRUNE SAMPLED TORSIONS BY HRMSD
        //If enabled, sort and prune the torsions to account for any overlap
        //Activation is done in the sampling function - not necessary here
        if (dn_prune_individual_torsions && torsions_collector.size() > 1) {
            prune_h_rmsd_and_mw(torsions_collector);  
        }

        //cout << "Number of fragments after torsion pruning:  " << torsions_collector.size() << endl;
        //all of these torsions are valid for the moment - push them to the final ensemble for this round
        for (int mol = 0; mol < torsions_collector.size(); mol++){
            ensemble_collector.push_back(torsions_collector[mol]);
        }

        //clear torsions for next molecule, and increment iterator position
        //torsions_collector.clear();
        pos++;
    }

    //Prune the combined ensemble's output
    if (dn_prune_combined_torsions){
        prune_h_rmsd_and_mw(ensemble_collector);
    }

    //Collect the sampled molecules back onto growing
    for (int mol=0; mol<ensemble_collector.size(); mol++){
        growing_ref.push_back(ensemble_collector[mol]);
    }
    ensemble_collector.clear();
    return;

}
/*
// +++++++++++++++++++++++++++++++++++++++++
// Sample the contents of a given fragment library randomly
void
DN_Build::sample_fraglib_matrix( Fragment & layer_frag, int j, vector <Fragment> & fraglib,
                                 vector <Fragment> & growing_ref, Master_Score & score,
                                 Simplex_Minimizer & simplex, AMBER_TYPER & typer, bool last_layer )
{

    // If the fragment library is empty for some reason, return without doing anything
    if (fraglib.size() == 0) { return; }
    // While N choices...
    int choice=0;
    int passes=0;
    int dnrn;

    std::vector <int> picked;
    //cout << "START FRAG NAME: " << layer_frag.mol.energy << endl;

    vector<pair<string, float> > possible_attachments;

    for (int i=0; i<ordered_frag_list.size(); ++i){
        if (ordered_frag_list[i] == layer_frag.aps[j].frag_name){
            for (int x=0; x<fragment_frequency_matrix[i].size(); ++x){
                if (fragment_frequency_matrix[i][x] > 0){
                    cout << ordered_frag_list[i] << " to " << ordered_frag_list[x] << ": " << fragment_frequency_matrix[i][x] << endl;

                    possible_attachments.push_back(make_pair(ordered_frag_list[x],rel_fragment_frequency_matrix[i][x]));
                }
            }
        }
    }
    
    for (int entry=0; entry<possible_attachments.size();++entry){
        cout << "This fragment can be attached to " << possible_attachments[entry].first << " with frequency " << possible_attachments[entry].second << endl;
    }
    

    int attempts=0;
    //cout << "size " << possible_attachments.size() << endl;
    while (attempts < possible_attachments.size()){
        //cout << "Selection attempt " << attempts << endl;
        int index_selected;
        for (int entry=0; entry<fraglib.size(); ++entry){
            if (fraglib[entry].mol.energy == possible_attachments[attempts].first){
                //cout << "Trying fragment " << fraglib[entry].mol.energy << endl;
                index_selected = entry;
                break;
            }
        }
        for (int m=0; m<fraglib[index_selected].aps.size(); m++){
            //cout << "Attempts " << attempts << endl;
            if (verbose) cout << "Trying attachment point #" << m << endl; //LEP - 03.24.18
            //cout << "Frag_Name: " << fraglib[index_selected].mol.energy << endl;
            //cout << "AP FRAG: " << fraglib[index_selected].aps[m].frag_name << endl;
            //cout << "m " << m << endl;
            // Check to see if the bond order is compatible
            if (compare_dummy_bonds(layer_frag, j, fraglib[index_selected], m)) {
            
                // If it is compatible, combine the two fragments to make a new, larger fragment
                Fragment frag_combined = combine_fragments(layer_frag, j, fraglib[index_selected], m);
                //frag_combined.mol.cout_information();
                /*
                for (int i = 0; i<frag_combined.mol.num_atoms; i++) {
                    cout << frag_combined.mol.atom_number[i] << " " << frag_combined.mol.x[i] << " " 
                    << frag_combined.mol.y[i] << " " 
                    << frag_combined.mol.z[i] << endl;
                }
                
                // check growing_ref to see if it gets bigger
                int growing_ref_size_before = growing_ref.size(); 
                // If the torenv is turned on...
                if (dn_use_torenv_table){
                    bool valid_torenv_bool=false;
                    // First check to see if the newest torsion environment is allowable
                    // (note that there may be some dummy atoms included in the environment, and
                    // even though it passes a check now, it may not pass the check later)
                    if(dn_use_roulette){
                        valid_torenv_bool = roulette_valid_torenv(frag_combined);
                    } else {
                        valid_torenv_bool = valid_torenv(frag_combined);
                    }
                    if (valid_torenv_bool){

                        // Second check to see whether all of the old torsion environments are
                        // still compatible. This is necessary because of the dummy as wild card
                        // thing
                        if (valid_torenv_multi(frag_combined)){
                            // If both checks pass, then this is a valid molecule. Sample torsions.
                            sample_minimized_torsions(frag_combined, growing_ref, score, simplex, typer);
                        }
                    }
                // If the torenv is not turned on, just skip straight to sampling torsions
                } else {
                    sample_minimized_torsions(frag_combined, growing_ref, score, simplex, typer);
                }

                if (dn_prune_individual_torsions && torsions_collector.size() > 1) {
                    prune_h_rmsd_and_mw(torsions_collector);  
                }

                // If we got this far, and if the fragment we are looking at is a scaffold,
                // and if growing_ref increased in size, then we added a scaffold this layer
                if ( (fraglib[index_selected].aps.size() > 2) && 
                     (growing_ref.size() > growing_ref_size_before)){ 

                    layer_frag.scaffolds_this_layer += 1;
                }
            }
        }


        attempts++; 
    }   
    //Prune the combined ensemble's output
    if (dn_prune_combined_torsions){
        prune_h_rmsd_and_mw(ensemble_collector);
    }

    //Collect the sampled molecules back onto growing
    for (int mol=0; mol<ensemble_collector.size(); mol++){
        growing_ref.push_back(ensemble_collector[mol]);
    }

    ensemble_collector.clear();
    return;
} // end DN_Build::sample_fraglib_matrix()

*/
/*
                    // If the torenv is not turned on, just skip straight to sampling torsions                    
                    } else {
                        sample_minimized_torsions(frag_combined, growing_ref, score, simplex, typer);
                        if (!counted ) { //BCF TODO
                            i++; //increment random pick iterator
                            tries++; //increment graph pick iterator
                            counted = true; //don't want to count multiple aps on one frag
                        }
                    }

                    // If there were no torsions generated, go to the next iteration of this loop
                    if (growing_ref.size() == growing_ref_size_before){ continue; }

                    // If we got this far, and if the fragment we are looking at is a scaffold,
                    // and if growing_ref increased in size, then we added a scaffold this layer
                    if ( fraglib[this_dnrn.first].aps.size() > 2 && 
                         growing_ref.size() > growing_ref_size_before ){ 

                        layer_frag.scaffolds_this_layer += 1;
                    }

                    // 1/19/16 select energy of best conformer generated
                    float best_tors_score = dn_pruning_conformer_score_cutoff; //higher than all scores 
                    for (int grc = (growing_ref_size_before); grc < growing_ref.size(); grc++) {
                        if (growing_ref[grc].mol.current_score < best_tors_score) {
                            best_tors_score = growing_ref[grc].mol.current_score; }
                    }

//TODO
// Output statistics for this step - how many kept, how many rejected at each temperature

                    // For simulated annealing of the next score comparison. Old score minus new
                    // score will be a positive number if the change is good.
                    float temp1 = exp( (compare_score - best_tors_score) /
                                        dn_graph_temperature );
                    float temp2 = ((double) rand() / RAND_MAX);

                    // Check the score of what we just kept, see if it is an improvement over
                    // the previous best score according to {  exp(-(Ef-Ei)/T) > Rand }
                    if ( temp1 > temp2 ){
                        improvement_flag = true;

                        // Do this second check because you can be trying one new fragment in 
                        // multiple ways, and you just want to keep the best score
                        if (best_tors_score < compare_score){
                            compare_score = best_tors_score;
                        }
                    }
                }
            } // end for each attachment point on the candidate fragment

            // If this fragment improved the score over the previous score, add its neighbors
            // to the sampling_queue (only the unvisited neighbors)
            if (improvement_flag){
                int temp_breadth = dn_graph_breadth;

                // If this is a very small frag graph or a very large breadth, reset the breadt
                // to be the size of the graph
                if (fraggraph_vec.size() <= temp_breadth){ temp_breadth = ( fraggraph_vec.size() - 1 ); }

                // For every child that will be sampled at this depth
                for (int n=1; n<temp_breadth+1; n++){

                    // BCF check to prevent seg fault
                    if (fraggraph_vec[this_dnrn.first].rankvec.size() == n) {break;}

                    // Make sure that dn_max_current_aps won't be exceeded by adding this fragment
                    temp_aps = fraglib[ fraggraph_vec[this_dnrn.first].rankvec[n] ].aps.size() +
                               layer_frag.aps.size() - 2;
                    if (temp_aps > dn_max_current_aps){
                        fraggraph_vec[ fraggraph_vec[this_dnrn.first].rankvec[n] ].visited = true;
                    }

                    // Only check it if it has not already been visited
                    if ( !fraggraph_vec[ fraggraph_vec[this_dnrn.first].rankvec[n] ].visited ){

                        // Make a record of its position and the current score, and it the frag
                        // to the sampling queue
                        pair <int, float> temp_pair = make_pair(fraggraph_vec[this_dnrn.first].rankvec[n],
                                                                compare_score);
                        sampling_queue.push_back( temp_pair );
                    }
                }
            }

            improvement_flag = false;

        } // end while sampling_queue is not empty

        sampling_queue.clear();

    } // end for num starting points


    // reset the visited flags
    for (int i=0; i<fraggraph_vec.size(); i++){
        fraggraph_vec[i].visited = false;
    }

    //cout << "sampling complete" << endl;
    return;

} // end commented out Matrix code
*/

// +++++++++++++++++++++++++++++++++++++++++
// Given two fragments and connection point data, combine them into one and return it
std::pair<std::vector<std::vector<double>>,Fragment>
DN_Build::combine_fragments_iso( Fragment& frag1, int attachment_point_1,
                                 Fragment  frag2, int attachment_point_2)
{


    std::pair<std::vector<std::vector<double>>,Fragment> return_pair {};
    return_pair.first = {{0,0,0},
                         {0,0,0},
                         {0,0,0}};
    Fragment return_frag;
    return_pair.second = return_frag;

    //assign some variables from the input fragments
    //frag1 is the fragment we are attaching TO, and frag2
    //is the fragment we are trying to attach.
    int dummy1 = frag1.aps[attachment_point_1].dummy_atom;
    int dummy2 = frag2.aps[attachment_point_2].dummy_atom;
    int heavy1 = frag1.aps[attachment_point_1].heavy_atom;
    int heavy2 = frag2.aps[attachment_point_2].heavy_atom;


    // frag1 is in the correct position - the objective is to translate / rotate frag2
    // so that the bond from dummy2->heavy2 is overlapping the bond from heavy1->dummy1


    // Step 1. Adjust the bond length of frag2 to match covalent radii of new pair

    // Calculate the desired bond length and remember as 'new_rad'
    float new_rad = calc_cov_radius(frag1.mol.atom_types[heavy1]) +
                    calc_cov_radius(frag2.mol.atom_types[heavy2]);

    // Calculate the x-y-z components of the current frag2 bond vector (bond_vec)
    DOCKVector bond_vec;
    bond_vec.x = frag2.mol.x[heavy2] - frag2.mol.x[dummy2];
    bond_vec.y = frag2.mol.y[heavy2] - frag2.mol.y[dummy2];
    bond_vec.z = frag2.mol.z[heavy2] - frag2.mol.z[dummy2];

    // Normalize the bond vector then multiply each component by new_rad so that it is the desired
    // length
    bond_vec = bond_vec.normalize_vector();
    bond_vec.x *= new_rad;
    bond_vec.y *= new_rad;
    bond_vec.z *= new_rad;

    // Change the coordinates of the frag2 dummy atom so that the bond length is correct
    frag2.mol.x[dummy2] = frag2.mol.x[heavy2] - bond_vec.x;
    frag2.mol.y[dummy2] = frag2.mol.y[heavy2] - bond_vec.y;
    frag2.mol.z[dummy2] = frag2.mol.z[heavy2] - bond_vec.z;


    // Step 2. Translate dummy2 of frag2 to the origin

    // Figure out what translation is required to move the dummy atom to the origin
    DOCKVector trans1;
    trans1.x = -frag2.mol.x[dummy2];
    trans1.y = -frag2.mol.y[dummy2];
    trans1.z = -frag2.mol.z[dummy2];

    // Use the dockmol function to translate the fragment so the dummy atom is at the origin
    frag2.mol.translate_mol(trans1);


    // Step 3. Calculate dot product to determine theta (theta = angle between vec1 and vec2)

    // vec1 = vector pointing from heavy1 to dummy1 in frag1
    DOCKVector vec1;
    vec1.x = frag1.mol.x[dummy1] - frag1.mol.x[heavy1];
    vec1.y = frag1.mol.y[dummy1] - frag1.mol.y[heavy1];
    vec1.z = frag1.mol.z[dummy1] - frag1.mol.z[heavy1];

    // vec2 = vector pointing from dummy2 to heavy2 in frag2 (dummy2 is at the origin)
    DOCKVector vec2;
    vec2.x = frag2.mol.x[heavy2];
    vec2.y = frag2.mol.y[heavy2];
    vec2.z = frag2.mol.z[heavy2];

    // Declare some variables
    float dot = 0.0;          // dot product value of vec1 and vec2
    float vec1_magsq = 0.0;   // vec1 magnitude-squared
    float vec2_magsq = 0.0;   // vec2 magnitude-squared
    float cos_theta = 0.0;    // cosine of theta
    float sin_theta = 0.0;    // sine of theta

    //initialize final matrix
    double finalmat[3][3];

    // Compute the dot product using the function in utils.cpp
    dot = dot_prod(vec1, vec2);

    // Compute these magnitudes (squared)
    vec1_magsq = (vec1.x * vec1.x) + (vec1.y * vec1.y) + (vec1.z * vec1.z);
    vec2_magsq = (vec2.x * vec2.x) + (vec2.y * vec2.y) + (vec2.z * vec2.z);

    // Compute cosine and sine of theta (theta itself is not actually calculated)
    cos_theta = dot / (sqrt (vec1_magsq * vec2_magsq));
    //sin_theta = sqrt (1 - (cos_theta * cos_theta));

    // cos_theta should never be greater than one. 
    // if it is then it is a numeric issue. 
    if (cos_theta >= 1.0){
       cos_theta = 1.0;
       sin_theta = 0.0;
    }
    else if (cos_theta <= -1.0){ // likewise if it is less than negative one.  
       cos_theta = -1.0;
       sin_theta = 0.0;
    }
    else{
       sin_theta = sqrt (1 - (cos_theta * cos_theta));
    }

    //cout << "dot product info: " << dot << " " <<
    //        vec1_magsq << " " <<
    //        vec2_magsq << " " <<
    //        cos_theta << " " <<
    //        sin_theta <<endl;
    // dot product info: 3.60413 3.24744 4 1 -nan

    // Step 4. Rotate vec2 to be coincident with vec1

    // If cos_theta is -1, the vectors are parallel but in the opposite direction
    if (cos_theta == -1){

        //fill the rotation matrix and rotate frag2
        finalmat[0][0] = -1;  finalmat[0][1] =  0;  finalmat[0][2] =  0;
        finalmat[1][0] =  0;  finalmat[1][1] = -1;  finalmat[1][2] =  0;
        finalmat[2][0] =  0;  finalmat[2][1] =  0;  finalmat[2][2] = -1;

        frag2.mol.rotate_mol(finalmat);
    }

    // If cos_theta is 1, vec1 and vec2 are already parallel - only translation is needed.
    else if (cos_theta == 1){
        finalmat[0][0] =  1; finalmat[0][1] = 0; finalmat[0][2] = 0; 
        finalmat[1][0] =  0; finalmat[1][1] = 1; finalmat[1][2] = 0; 
        finalmat[2][0] =  0; finalmat[2][1] = 0; finalmat[2][2] = 1; 

        frag2.mol.rotate_mol(finalmat);

    // Otherwise, enter this loop and calculate out how to rotate frag2
    } else if (cos_theta != 1) {

        // Calculate cross product of vec1 and vec2 to get U (function from utils.cpp)
        DOCKVector normalU = cross_prod(vec1, vec2);

        // Calculate cross product of vec2 and U to get ~W
        DOCKVector normalW = cross_prod(vec2, normalU);

        // Normalize the vectors
        vec2 = vec2.normalize_vector();
        normalU = normalU.normalize_vector();
        normalW = normalW.normalize_vector();


        // (1) Make coordinate rotation matrix, which rotates {e1, e2, e3} coordinate to
        // {normalW, vec2, normalU} coordinate
        float coorRot[3][3];
        coorRot[0][0] = normalW.x;  coorRot[0][1] = vec2.x;  coorRot[0][2] = normalU.x;
        coorRot[1][0] = normalW.y;  coorRot[1][1] = vec2.y;  coorRot[1][2] = normalU.y;
        coorRot[2][0] = normalW.z;  coorRot[2][1] = vec2.z;  coorRot[2][2] = normalU.z;

        // (3) Make inverse  matrix of coorRot matrix - since coorRot is an orthogonal matrix,
        // the inverse is its transpose, (coorRot)^T
        float invcoorRot[3][3];
        invcoorRot[0][0] = coorRot[0][0];  
        invcoorRot[0][1] = coorRot[1][0];
        invcoorRot[0][2] = coorRot[2][0];
       
        invcoorRot[1][0] = coorRot[0][1];  
        invcoorRot[1][1] = coorRot[1][1];
        invcoorRot[1][2] = coorRot[2][1];

        invcoorRot[2][0] = coorRot[0][2];  
        invcoorRot[2][1] = coorRot[1][2];
        invcoorRot[2][2] = coorRot[2][2];

      /************************ TEB, 2020/04/03
       * Why the angle is never negative? 
       *
       * It is because the rotation matrices [ norm_W ,  norm_v2 , norm_U] transposed  
       * includes norm_U which is determined by taking the cross product of the vectors 
       * v1 and v2.  if the angle between v1 and v2 is positive (counterclockwise) then 
       * the vector is pointed out of the plane, if the angle is negative then the vector
       * (norm_U) is pointing into the plane.  This vector (norm_U) is placed onto the 
       * positive z axis.  Thus, even, if the the cross-poduct is pointing into the plane
       * when it is transformed it is point out. 
       *
       ************************  
       */

      /*
        // check the sign of the angle, and sign. TEB, 2020/04/02
        float v1[3], v2[3];
        v1[0] = vec1.x; v2[0] = vec2.x; 
        v1[1] = vec1.y; v2[1] = vec2.y; 
        v1[2] = vec1.z; v2[2] = vec2.z; 
        
        float sign = check_neg_angle(v1,v2,invcoorRot);
        //float sign = check_neg_angle(v2,v1,invcoorRot);
        // sin(-theata) = -sin(theata)
        sin_theta = sign * sin_theta;
       */

        // (2) Make rotation matrix, which rotates vec2 theta angle on a plane of vec2 and normalW
        // to the direction of normalW
        float planeRot[3][3];
        planeRot[0][0] =  cos_theta;  planeRot[0][1] = sin_theta;  planeRot[0][2] = 0;
        planeRot[1][0] =  -sin_theta;  planeRot[1][1] =  cos_theta;  planeRot[1][2] = 0;
        planeRot[2][0] =          0;  planeRot[2][1] =          0;  planeRot[2][2] = 1;


        // (4) Multiply three matrices together:  [coorRot * planeRot * invcoorRot]
        float temp[3][3];

        // First multiply coorRot * planeRot, save as temp
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                temp[i][j] = 0.0;
                for (int k=0; k<3; k++){
                    temp[i][j] += coorRot[i][k]*planeRot[k][j];
                }
            }
        }

        // Then multiply temp * invcoorRot, save as finalmat
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                finalmat[i][j] = 0.0;
                for (int k=0; k<3; k++){
                    finalmat[i][j] += temp[i][k]*invcoorRot[k][j];
                }
            }
        }

        // Rotate frag2 using finalmat[3][3]
        frag2.mol.rotate_mol(finalmat);
    }


    // Step 5. Translate frag2 to frag1

    // This is the translation vector to move dummy2 to heavy1 (dummy2 is at the origin)
    DOCKVector trans2;
    trans2.x = frag1.mol.x[heavy1];
    trans2.y = frag1.mol.y[heavy1];
    trans2.z = frag1.mol.z[heavy1];

    // Use the dockmol function to translate frag2
    frag2.mol.translate_mol(trans2);


    //// Step 6. Connect the two fragments and return one object
    //Fragment frag_combined = attach(frag1, dummy1, heavy1, frag2, dummy2, heavy2);

    //// Also, copy the growth tree over if those are turned on
    //if (dn_write_growth_trees){
    //    for (int i=0; i<frag1.frag_growth_tree.size(); i++){
    //        frag_combined.frag_growth_tree.push_back(frag1.frag_growth_tree[i]);
    //    }
    //}

    return_pair.first[0] = {finalmat[0][0], finalmat[0][1], finalmat[0][2]};
    return_pair.first[1] = {finalmat[1][0], finalmat[1][1], finalmat[1][2]};
    return_pair.first[2] = {finalmat[2][0], finalmat[2][1], finalmat[2][2]};

    return_pair.second = frag2;

    return return_pair;

} // end DN_Build::combine_fragments_iso()

// +++++++++++++++++++++++++++++++++++++++++
// Given two fragments and connection point data, combine them into one and return it
Fragment
DN_Build::combine_fragments( Fragment & frag1, int attachment_point_1,
                             Fragment frag2,   int attachment_point_2)
{

    //assign some variables from the input fragments
    //frag1 is the fragment we are attaching TO, and frag2
    //is the fragment we are trying to attach.
    int dummy1 = frag1.aps[attachment_point_1].dummy_atom;
    int dummy2 = frag2.aps[attachment_point_2].dummy_atom;
    int heavy1 = frag1.aps[attachment_point_1].heavy_atom;
    int heavy2 = frag2.aps[attachment_point_2].heavy_atom;
    // frag1 is in the correct position - the objective is to translate / rotate frag2
    // so that the bond from dummy2->heavy2 is overlapping the bond from heavy1->dummy1


    // Step 1. Adjust the bond length of frag2 to match covalent radii of new pair

    // Calculate the desired bond length and remember as 'new_rad'
    float new_rad = calc_cov_radius(frag1.mol.atom_types[heavy1]) +
                    calc_cov_radius(frag2.mol.atom_types[heavy2]);

    // Calculate the x-y-z components of the current frag2 bond vector (bond_vec)
    DOCKVector bond_vec;
    bond_vec.x = frag2.mol.x[heavy2] - frag2.mol.x[dummy2];
    bond_vec.y = frag2.mol.y[heavy2] - frag2.mol.y[dummy2];
    bond_vec.z = frag2.mol.z[heavy2] - frag2.mol.z[dummy2];

    // Normalize the bond vector then multiply each component by new_rad so that it is the desired
    // length
    bond_vec = bond_vec.normalize_vector();
    bond_vec.x *= new_rad;
    bond_vec.y *= new_rad;
    bond_vec.z *= new_rad;

    // Change the coordinates of the frag2 dummy atom so that the bond length is correct
    frag2.mol.x[dummy2] = frag2.mol.x[heavy2] - bond_vec.x;
    frag2.mol.y[dummy2] = frag2.mol.y[heavy2] - bond_vec.y;
    frag2.mol.z[dummy2] = frag2.mol.z[heavy2] - bond_vec.z;


    // Step 2. Translate dummy2 of frag2 to the origin

    // Figure out what translation is required to move the dummy atom to the origin
    DOCKVector trans1;
    trans1.x = -frag2.mol.x[dummy2];
    trans1.y = -frag2.mol.y[dummy2];
    trans1.z = -frag2.mol.z[dummy2];

    // Use the dockmol function to translate the fragment so the dummy atom is at the origin
    frag2.mol.translate_mol(trans1);


    // Step 3. Calculate dot product to determine theta (theta = angle between vec1 and vec2)

    // vec1 = vector pointing from heavy1 to dummy1 in frag1
    DOCKVector vec1;
    vec1.x = frag1.mol.x[dummy1] - frag1.mol.x[heavy1];
    vec1.y = frag1.mol.y[dummy1] - frag1.mol.y[heavy1];
    vec1.z = frag1.mol.z[dummy1] - frag1.mol.z[heavy1];

    // vec2 = vector pointing from dummy2 to heavy2 in frag2 (dummy2 is at the origin)
    DOCKVector vec2;
    vec2.x = frag2.mol.x[heavy2];
    vec2.y = frag2.mol.y[heavy2];
    vec2.z = frag2.mol.z[heavy2];

    // Declare some variables
    float dot;          // dot product value of vec1 and vec2
    float vec1_magsq;   // vec1 magnitude-squared
    float vec2_magsq;   // vec2 magnitude-squared
    float cos_theta;    // cosine of theta
    float sin_theta;    // sine of theta

    // Compute the dot product using the function in utils.cpp
    dot = dot_prod(vec1, vec2);

    // Compute these magnitudes (squared)
    vec1_magsq = (vec1.x * vec1.x) + (vec1.y * vec1.y) + (vec1.z * vec1.z);
    vec2_magsq = (vec2.x * vec2.x) + (vec2.y * vec2.y) + (vec2.z * vec2.z);

    // Compute cosine and sine of theta (theta itself is not actually calculated)
    cos_theta = dot / (sqrt (vec1_magsq * vec2_magsq));
    //sin_theta = sqrt (1 - (cos_theta * cos_theta));

    // cos_theta should never be greater than one. 
    // if it is then it is a numeric issue. 
    if (cos_theta >= 1.0){
       cos_theta = 1.0;
       sin_theta = 0.0;
    }
    else if (cos_theta <= -1.0){ // likewise if it is less than negative one.  
       cos_theta = -1.0;
       sin_theta = 0.0;
    }
    else{
       sin_theta = sqrt (1 - (cos_theta * cos_theta));
    }

    //cout << "dot product info: " << dot << " " <<
    //        vec1_magsq << " " <<
    //        vec2_magsq << " " <<
    //        cos_theta << " " <<
    //        sin_theta <<endl;
    // dot product info: 3.60413 3.24744 4 1 -nan

    // Step 4. Rotate vec2 to be coincident with vec1

    // If cos_theta is -1, the vectors are parallel but in the opposite direction
    if (cos_theta == -1){

        // Declare the rotation matrix and rotate frag2
        double finalmat[3][3] = { { -1, 0, 0}, {0, -1, 0}, {0, 0, -1} };
        frag2.mol.rotate_mol(finalmat);
    }

    // If cos_theta is 1, vec1 and vec2 are already parallel - only translation is needed.
    // Otherwise, enter this loop and calculate out how to rotate frag2
    else if (cos_theta != 1) {

        // Calculate cross product of vec1 and vec2 to get U (function from utils.cpp)
        DOCKVector normalU = cross_prod(vec1, vec2);

        // Calculate cross product of vec2 and U to get ~W
        DOCKVector normalW = cross_prod(vec2, normalU);

        // Normalize the vectors
        vec2 = vec2.normalize_vector();
        normalU = normalU.normalize_vector();
        normalW = normalW.normalize_vector();


        // (1) Make coordinate rotation matrix, which rotates {e1, e2, e3} coordinate to
        // {normalW, vec2, normalU} coordinate
        float coorRot[3][3];
        coorRot[0][0] = normalW.x;  coorRot[0][1] = vec2.x;  coorRot[0][2] = normalU.x;
        coorRot[1][0] = normalW.y;  coorRot[1][1] = vec2.y;  coorRot[1][2] = normalU.y;
        coorRot[2][0] = normalW.z;  coorRot[2][1] = vec2.z;  coorRot[2][2] = normalU.z;

        // (3) Make inverse  matrix of coorRot matrix - since coorRot is an orthogonal matrix,
        // the inverse is its transpose, (coorRot)^T
        float invcoorRot[3][3];
        invcoorRot[0][0] = coorRot[0][0];  invcoorRot[0][1] = coorRot[1][0];
        invcoorRot[0][2] = coorRot[2][0];
       
        invcoorRot[1][0] = coorRot[0][1];  invcoorRot[1][1] = coorRot[1][1];
        invcoorRot[1][2] = coorRot[2][1];

        invcoorRot[2][0] = coorRot[0][2];  invcoorRot[2][1] = coorRot[1][2];
        invcoorRot[2][2] = coorRot[2][2];

      /************************ TEB, 2020/04/03
       * Why the angle is never negative? 
       *
       * It is because the rotation matrices [ norm_W ,  norm_v2 , norm_U] transposed  
       * includes norm_U which is determined by taking the cross product of the vectors 
       * v1 and v2.  if the angle between v1 and v2 is positive (counterclockwise) then 
       * the vector is pointed out of the plane, if the angle is negative then the vector
       * (norm_U) is pointing into the plane.  This vector (norm_U) is placed onto the 
       * positive z axis.  Thus, even, if the the cross-poduct is pointing into the plane
       * when it is transformed it is point out. 
       *
       ************************  
       */

      /*
        // check the sign of the angle, and sign. TEB, 2020/04/02
        float v1[3], v2[3];
        v1[0] = vec1.x; v2[0] = vec2.x; 
        v1[1] = vec1.y; v2[1] = vec2.y; 
        v1[2] = vec1.z; v2[2] = vec2.z; 
        
        float sign = check_neg_angle(v1,v2,invcoorRot);
        //float sign = check_neg_angle(v2,v1,invcoorRot);
        // sin(-theata) = -sin(theata)
        sin_theta = sign * sin_theta;
       */

        // (2) Make rotation matrix, which rotates vec2 theta angle on a plane of vec2 and normalW
        // to the direction of normalW
        float planeRot[3][3];
        planeRot[0][0] =  cos_theta;  planeRot[0][1] = sin_theta;  planeRot[0][2] = 0;
        planeRot[1][0] = -sin_theta;  planeRot[1][1] = cos_theta;  planeRot[1][2] = 0;
        planeRot[2][0] =          0;  planeRot[2][1] =         0;  planeRot[2][2] = 1;


        // (4) Multiply three matrices together:  [coorRot * planeRot * invcoorRot]
        float temp[3][3];
        double finalmat[3][3];

        // First multiply coorRot * planeRot, save as temp
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                temp[i][j] = 0.0;
                for (int k=0; k<3; k++){
                    temp[i][j] += coorRot[i][k]*planeRot[k][j];
                }
            }
        }

        // Then multiply temp * invcoorRot, save as finalmat
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                finalmat[i][j] = 0.0;
                for (int k=0; k<3; k++){
                    finalmat[i][j] += temp[i][k]*invcoorRot[k][j];
                }
            }
        }

        // Rotate frag2 using finalmat[3][3]
        frag2.mol.rotate_mol(finalmat);
    }


    // Step 5. Translate frag2 to frag1

    // This is the translation vector to move dummy2 to heavy1 (dummy2 is at the origin)
    DOCKVector trans2;
    trans2.x = frag1.mol.x[heavy1];
    trans2.y = frag1.mol.y[heavy1];
    trans2.z = frag1.mol.z[heavy1];

    // Use the dockmol function to translate frag2
    frag2.mol.translate_mol(trans2);


    // Step 6. Connect the two fragments and return one object
    Fragment frag_combined = attach(frag1, dummy1, heavy1, frag2, dummy2, heavy2);

    // Also, copy the growth tree over if those are turned on
    if (dn_write_growth_trees){
        for (int i=0; i<frag1.frag_growth_tree.size(); i++){
            frag_combined.frag_growth_tree.push_back(frag1.frag_growth_tree[i]);
        }
    }

    return frag_combined;

} // end DN_Build::combine_fragments()



// +++++++++++++++++++++++++++++++++++++++++
// Return a covalent radius depending on atom type 
float
DN_Build::calc_cov_radius( string atom )
{
    // This function assumes that the atom type is a Sybyl atom type. These covalent
    // radii come from the CRC handbook, and are generalized by element (see header
    // file). Perhaps we could use some more rigorous numbers?

    if ( atom == "H")
        { return COV_RADII_H; }

    else if ( atom == "C.3" || atom == "C.2" || atom == "C.1" || atom == "C.ar" || atom  == "C.cat" )
        { return COV_RADII_C; }

    else if ( atom == "N.4" || atom == "N.3" || atom == "N.2" || atom == "N.1" || atom  == "N.ar" ||
              atom == "N.am" || atom == "N.pl3" )
        { return COV_RADII_N; }

    else if ( atom == "O.3" || atom == "O.2" || atom == "O.co2" )
        { return COV_RADII_O; }

    else if ( atom == "S.3" || atom == "S.2" || atom == "S.O" || atom == "S.o" || atom == "S.O2" ||
              atom == "S.o2" )
        { return COV_RADII_S; }

    else if ( atom == "P.3" )
        { return COV_RADII_P; }

    else if ( atom == "F" )
        { return COV_RADII_F; }

    else if ( atom == "Cl" )
        { return COV_RADII_CL; }

    else if ( atom == "Br" )
        { return COV_RADII_BR; }

    else
        { cout <<"Warning: Did not recognize the atom_type " <<atom 
               <<" in DN_Build::calc_cov_radius()"  <<endl; 
          return 0.71; }

} // end DN_Build::calc_cov_radius()



// +++++++++++++++++++++++++++++++++++++++++
// Attach two fragments that are already overlaid at an attachment point
Fragment
DN_Build::attach( Fragment & frag1, int dummy1, int heavy1, 
                  Fragment & frag2, int dummy2, int heavy2 )
{
    // Inactivate all atoms, then activate just the two heavy atoms that are about to be connected
    // together
    for (int i=0; i<frag1.mol.num_atoms; i++){
        frag1.mol.atom_active_flags[i] = 0;
    }
    for (int i=0; i<frag2.mol.num_atoms; i++){
        frag2.mol.atom_active_flags[i] = 0;
    }
    frag1.mol.atom_active_flags[heavy1] = 1;
    frag2.mol.atom_active_flags[heavy2] = 1;


    // Declare the dockmol that will be the combination of both fragments
    DOCKMol combined_mol;

    // Number of atoms and bonds that will be in newly combined fragment
    // (two atoms are deleted) (two bonds are deleted, and one new bond is formed)
    int natoms, nbonds;
    natoms = (frag1.mol.num_atoms + frag2.mol.num_atoms) - 2;
    nbonds = (frag1.mol.num_bonds + frag2.mol.num_bonds) - 1;


    // Allocate arrays and make assignments for <Tripos>Molecule section of new DOCKMol
    combined_mol.allocate_arrays(natoms, nbonds, 1);

    combined_mol.title         =  frag1.mol.title;
    combined_mol.comment1      =  frag1.mol.comment1;
    combined_mol.comment2      =  frag1.mol.comment2;
    combined_mol.comment3      =  frag1.mol.comment3;

    // Make a unique molecule name
    ostringstream new_name;
    new_name <<frag1.mol.energy <<"--" <<frag2.mol.energy;
    combined_mol.energy = new_name.str();
    if (verbose) cout << "Trying fragment: " << frag2.mol.energy << endl; //LEP - 03.24.18
    new_name.clear();

    // Fill in some more of the data
    combined_mol.num_atoms = natoms;
    combined_mol.num_bonds = nbonds;

    // Read in atom info of frag1, then frag2, ignoring dummy1 and dummy2
    int atom_index = 0;
    int dummy1_position;
    int largest_resid = 0;

    for (int i=0; i<frag1.mol.num_atoms; i++){

        // When the correct dummy is encountered, skip it, but remember its position
        if (i == dummy1){
            dummy1_position = i;
            continue;
        }

        combined_mol.x[atom_index] = frag1.mol.x[i];
        combined_mol.y[atom_index] = frag1.mol.y[i];
        combined_mol.z[atom_index] = frag1.mol.z[i];

        combined_mol.charges[atom_index]              = frag1.mol.charges[i];
        combined_mol.atom_names[atom_index]           = frag1.mol.atom_names[i];
        combined_mol.atom_types[atom_index]           = frag1.mol.atom_types[i];
        combined_mol.atom_number[atom_index]          = frag1.mol.atom_number[i];
        combined_mol.atom_residue_numbers[atom_index] = frag1.mol.atom_residue_numbers[i];
        if ( atoi(frag1.mol.atom_residue_numbers[i].c_str()) > largest_resid ){
            largest_resid = atoi(frag1.mol.atom_residue_numbers[i].c_str());
        }
        combined_mol.subst_names[atom_index]          = frag1.mol.subst_names[i];

        // Active atoms are the heavy atoms that were just attached together
        combined_mol.atom_active_flags[atom_index]    = frag1.mol.atom_active_flags[i];

        atom_index++;
    }    



    // Note that atom_index is not re-initialized, it is the running total between both fragments
    int dummy2_position;
    int next_resid = largest_resid + 1;
    combined_mol.num_residues = next_resid;
    stringstream temp_ss;
    temp_ss << next_resid;
    for (int i=0; i<frag2.mol.num_atoms; i++){

        // When the correct dummy is encountered, skip it, but remember its position
        if (i == dummy2){
            dummy2_position = i;
            continue;
        }

        combined_mol.x[atom_index] = frag2.mol.x[i];
        combined_mol.y[atom_index] = frag2.mol.y[i];
        combined_mol.z[atom_index] = frag2.mol.z[i];

        combined_mol.charges[atom_index]              = frag2.mol.charges[i];
        combined_mol.atom_names[atom_index]           = frag2.mol.atom_names[i];
        combined_mol.atom_types[atom_index]           = frag2.mol.atom_types[i];
        combined_mol.atom_number[atom_index]          = frag2.mol.atom_number[i];
        combined_mol.atom_residue_numbers[atom_index] = temp_ss.str();
        combined_mol.subst_names[atom_index]          = frag2.mol.subst_names[i];

        // Active atoms are the heavy atoms that were just attached together
        combined_mol.atom_active_flags[atom_index]    = frag2.mol.atom_active_flags[i];
        atom_index++;
    }


    // Now that we know dummy1 position, we can update frag1.torenv_recheck_indices to fix the
    // atom indices. Frag2 is the new thing, so there are no bonds to recheck within it. This 
    // will later be copied onto the final attached fragment.
    pair <int, int> temp_pair;
    vector < pair <int, int> > temp_torenv_recheck_indices;

    // For every pair of atoms that needs to be rechecked...
    for (int i=0; i<frag1.torenv_recheck_indices.size(); i++){

        // If the atoms appear before dummy1, no change required, otherwise subtract 1

        if (frag1.torenv_recheck_indices[i].first < dummy1_position){
            temp_pair.first = frag1.torenv_recheck_indices[i].first;
        } else {
            temp_pair.first = frag1.torenv_recheck_indices[i].first - 1;
        }
    
        if (frag1.torenv_recheck_indices[i].second < dummy1_position){
            temp_pair.second = frag1.torenv_recheck_indices[i].second;
        } else {
            temp_pair.second = frag1.torenv_recheck_indices[i].second - 1;
        }
    
        // Add the pair onto the vector
        temp_torenv_recheck_indices.push_back( temp_pair);
    }


    // Go through all of the bond information in the fragments and add it to the new dockmol object
    // (two bonds should be deleted, one should be added)
    int bond_index = 0;
    string new_bond_type = "";
    for (int i=0; i<frag1.mol.num_bonds; i++){

        // Here, determine the bond type for the new bond between fragments
        if (frag1.mol.bonds_origin_atom[i] == dummy1 || frag1.mol.bonds_target_atom[i] == dummy1){
            new_bond_type = frag1.mol.bond_types[i];
            continue;
        }

        int tmp_origin_atom;
        int tmp_target_atom;

        // The index for tmp_origin_atom will need to be corrected if it occurs after the position
        // of the previously deleted dummy atom
        if (frag1.mol.bonds_origin_atom[i] > dummy1_position){
            tmp_origin_atom = frag1.mol.bonds_origin_atom[i] - 1;
        } else {
            tmp_origin_atom = frag1.mol.bonds_origin_atom[i];
        }

        // The index for tmp_target_atom will need to be corrected if it occurs after the position
        // of the previously deleted dummy atom
        if (frag1.mol.bonds_target_atom[i] > dummy1_position){
            tmp_target_atom = frag1.mol.bonds_target_atom[i] - 1;
        } else {
            tmp_target_atom = frag1.mol.bonds_target_atom[i];
        }

        combined_mol.bonds_origin_atom[bond_index] = tmp_origin_atom;
        combined_mol.bonds_target_atom[bond_index] = tmp_target_atom;
        combined_mol.bond_types[bond_index] = frag1.mol.bond_types[i];

        bond_index++;
    }

    // Now look through all the bonds of frag2...
    for (int i=0; i<frag2.mol.num_bonds; i++){

        // We already determined what this bond type is going to be, we will form the bond later
        if (frag2.mol.bonds_origin_atom[i] == dummy2 || frag2.mol.bonds_target_atom[i] == dummy2){
            continue;
        }

        int tmp_origin_atom;
        int tmp_target_atom;

        // Correct the index of tmp_origin_atom according to whether 2 or 1 dummy atoms that were
        // deleted occurred before them in the file
        if (frag2.mol.bonds_origin_atom[i] > dummy2_position){
            tmp_origin_atom = frag2.mol.bonds_origin_atom[i] + frag1.mol.num_atoms - 2;
        } else {
            tmp_origin_atom = frag2.mol.bonds_origin_atom[i] + frag1.mol.num_atoms - 1;
        }

        // Correct the index of tmp_target_atom according to whether 2 or 1 dummy atoms that were
        // deleted occurred before them in the file
        if (frag2.mol.bonds_target_atom[i] > dummy2_position){
            tmp_target_atom = frag2.mol.bonds_target_atom[i] + frag1.mol.num_atoms - 2;
        } else {
            tmp_target_atom = frag2.mol.bonds_target_atom[i] + frag1.mol.num_atoms - 1;
        }

        combined_mol.bonds_origin_atom[bond_index] = tmp_origin_atom;
        combined_mol.bonds_target_atom[bond_index] = tmp_target_atom;
        combined_mol.bond_types[bond_index] = frag2.mol.bond_types[i];

        bond_index++;
    }

    // Create the last bond which will be the connection between the two fragments
    int last_ap_heavy;
    if (heavy1 > dummy1_position){
        combined_mol.bonds_origin_atom[bond_index] = heavy1 - 1;
        last_ap_heavy = heavy1 - 1;
    } else {
        combined_mol.bonds_origin_atom[bond_index] = heavy1;
        last_ap_heavy = heavy1;
    }

    if (heavy2 > dummy2_position){
        combined_mol.bonds_target_atom[bond_index] = heavy2 + frag1.mol.num_atoms - 2;
    } else {
        combined_mol.bonds_target_atom[bond_index] = heavy2 + frag1.mol.num_atoms - 1;
    }

    // Set the bond type of that new bond here
    combined_mol.bond_types[bond_index] = new_bond_type;


    // Make a new mol_info_line that contains total number of atoms and bonds
    stringstream mol_info;
    mol_info <<natoms <<"  " <<nbonds <<"  " <<largest_resid+1;
    combined_mol.mol_info_line = mol_info.str();
    mol_info.clear();


    // Create the Fragment object that will be returned by this function, and
    // copy the dockmol object into it
    Fragment final;
    copy_molecule(final.mol, combined_mol);
    final.last_ap_heavy = last_ap_heavy;
    final.mol.id_ring_atoms_bonds();

    // Copy the torenv_recheck_indices over
    for (int i=0; i<temp_torenv_recheck_indices.size(); i++){
        final.torenv_recheck_indices.push_back(temp_torenv_recheck_indices[i]);
    }

    // The old vector is no longer needed, clear from memory
    temp_torenv_recheck_indices.clear();

    // Go through that DOCKMol and identify remaining Dummy atoms, save in <aps>
    // Create temporary vectors of atom indices
    vector <int> tmp_dummy_atoms;
    vector <int> tmp_heavy_atoms;

    // Loop over every bond in the dockmol object
    for (int i=0; i<final.mol.num_bonds; i++){

        // Check if the bond origin is a dummy atom
        if (final.mol.atom_types[final.mol.bonds_origin_atom[i]].compare("Du") == 0){

            // if so, save origin as dummy_atom and target as heavy atom
            tmp_dummy_atoms.push_back(final.mol.bonds_origin_atom[i]);
            tmp_heavy_atoms.push_back(final.mol.bonds_target_atom[i]);
        }

        // Check if the bond target is a dummy atom
        if (final.mol.atom_types[final.mol.bonds_target_atom[i]].compare("Du") == 0){

            // if so, save target as dummy atom and origin as heavy atom
            tmp_dummy_atoms.push_back(final.mol.bonds_target_atom[i]);
            tmp_heavy_atoms.push_back(final.mol.bonds_origin_atom[i]);
        }
    }
    /*
    cout << endl;
    for (int i = 0; i < frag1.aps.size(); i++){
        cout << "Frag1 Heavy: " << frag1.aps[i].heavy_atom << "  |  "
         << "Frag1 Dummy " << frag1.aps[i].dummy_atom << endl;

        cout << "Heavy " << frag1.mol.x[frag1.aps[i].heavy_atom] << " " 
        << frag1.mol.y[frag1.aps[i].heavy_atom] << " " 
        << frag1.mol.z[frag1.aps[i].heavy_atom] << endl;

        cout << "Dummy " << frag1.mol.x[frag1.aps[i].dummy_atom] << " " 
        << frag1.mol.y[frag1.aps[i].dummy_atom] << " " 
        << frag1.mol.z[frag1.aps[i].dummy_atom] << endl;
    }
    
    cout << endl;
    for (int i = 0; i < frag2.aps.size(); i++){
        cout << "Frag2 Heavy: " << frag2.aps[i].heavy_atom << "  |  "
         << "Frag2 Dummy " << frag2.aps[i].dummy_atom << endl;

        cout << "Heavy " << frag2.mol.x[frag2.aps[i].heavy_atom] << " " 
        << frag2.mol.y[frag2.aps[i].heavy_atom] << " " 
        << frag2.mol.z[frag2.aps[i].heavy_atom] << endl;

        cout << "Dummy " << frag2.mol.x[frag2.aps[i].dummy_atom] << " " 
        << frag2.mol.y[frag2.aps[i].dummy_atom] << " " 
        << frag2.mol.z[frag2.aps[i].dummy_atom] << endl;
    }

    cout << endl;
    for (int i=0; i<tmp_dummy_atoms.size(); i++){
        cout << "Final dummy " << tmp_dummy_atoms[i] << endl;
        cout << final.mol.x[tmp_dummy_atoms[i]] << " " << final.mol.y[tmp_dummy_atoms[i]]
        << " " << final.mol.z[tmp_dummy_atoms[i]] << endl;
    }
    */
    // Create temporary attachment point object 
    AttPoint tmp_ap;

    // Loop over all of the dummy_atom / heavy_atom pairs
    for (int i=0; i<tmp_dummy_atoms.size(); i++){

        // Copy the heavy atom / dummy atom pairs onto the attachment point class
        tmp_ap.heavy_atom = tmp_heavy_atoms[i];
        tmp_ap.dummy_atom = tmp_dummy_atoms[i];
        // And add that information to the attachment points (aps) vector of the fragment
        final.aps.push_back(tmp_ap);
    }

    for (int i=0; i<final.aps.size(); i++){
        for (int x=0; x<frag1.aps.size(); x++){
            if (final.mol.x[final.aps[i].dummy_atom] == frag1.mol.x[frag1.aps[x].dummy_atom] &&
                final.mol.y[final.aps[i].dummy_atom] == frag1.mol.y[frag1.aps[x].dummy_atom] &&
                final.mol.z[final.aps[i].dummy_atom] == frag1.mol.z[frag1.aps[x].dummy_atom]) {
                
                //cout << "side1 " << frag1.aps[x].frag_name << endl;
                final.aps[i].frag_name = frag1.aps[x].frag_name;
            }
        }
        for (int x=0; x<frag2.aps.size(); x++){
            if (final.mol.x[final.aps[i].dummy_atom] == frag2.mol.x[frag2.aps[x].dummy_atom] &&
                final.mol.y[final.aps[i].dummy_atom] == frag2.mol.y[frag2.aps[x].dummy_atom] &&
                final.mol.z[final.aps[i].dummy_atom] == frag2.mol.z[frag2.aps[x].dummy_atom]) {

                //cout << "side2 " << frag2.aps[x].frag_name << endl;
                final.aps[i].frag_name = frag2.aps[x].frag_name;
            }
        }
    }
    // Clear some vectors
    tmp_dummy_atoms.clear();
    tmp_heavy_atoms.clear();

    //BCF new fragment needs to remember how many scaffolds were added this layer
    //it will be reset at the start of a new layer
    final.scaffolds_this_layer = frag1.scaffolds_this_layer;

    return final;

} // end DN_Build::attach()

void
DN_Build::attach_selected_frags( Fragment & layer_frag, int layer_frag_ap, vector <Fragment> & fraglib, 
                                vector<int> & selected_fragments, vector< pair <Fragment, int> > & successful_attachments)
{
    //for verbose, initialize these variables
    int num_succ_att = 0;
    bool success_att;


    //Check all fragments to see if they're valid before attempting an attachment
    vector< int >::iterator pos = selected_fragments.begin();
    while (pos != selected_fragments.end()){
        // (1) Would create too many scaffolds per layer
        if (fraglib[(*pos)].aps.size() > 2 &&
            layer_frag.scaffolds_this_layer >= dn_max_scaffolds_per_layer ) {
            //cout << *pos << " # of aps " << fraglib[(*pos)].aps.size() << endl;
            //cout << "Frag " << *pos << " rejected for too many scaffolds." << endl;
            pos = selected_fragments.erase(pos);
            continue;
        }

        // (2) Would give too many current attachment points
        int temp_aps = fraglib[(*pos)].aps.size() + layer_frag.aps.size() - 2;
        if (temp_aps > dn_max_current_aps) {
            //cout << "Frag " << *pos << " rejected for too many aps." << endl;
            pos = selected_fragments.erase(pos);
            continue;
        }

        // (3) Rotatable bond cutoff
        int temp_rbs = layer_frag.mol.rot_bonds + temp_aps;

        // (5) Formal Charge Cutoff
        int temp_fc = fraglib[(*pos)].mol.formal_charge + layer_frag.mol.formal_charge;
        if (fabs(float(temp_fc)) > (fabs(float(dn_constraint_formal_charge))+0.1)) { 
            //cout << "Frag " << *pos << " rejected for formal charge." << endl;
            pos = selected_fragments.erase(pos);
            continue; 
        }
        pos++;
    }

    //if (Parameter_Reader::verbosity_level() > 1 ) { num_succ_att = 0; }
     num_succ_att = 0; 
    //Anything left in the selected_fragments vector has passed initial checks - attempt combos
    pos = selected_fragments.begin(); 
    while (pos != selected_fragments.end()){

        if (num_att_root > dn_max_sccful_att_per_root) { break; }

        //if (Parameter_Reader::verbosity_level() > 1 ) { success_att = false; }
        success_att = false;
        for (int fraglib_ap = 0; fraglib_ap < fraglib[(*pos)].aps.size(); fraglib_ap++){
            //check to see if the bond order is compatible
            if (compare_dummy_bonds(layer_frag, layer_frag_ap, fraglib[*pos], fraglib_ap)) {
                //combine the fragments
                Fragment frag_combined = combine_fragments(layer_frag, layer_frag_ap,
                                                            fraglib[*pos], fraglib_ap);

                //check if the new torsion environment from the attachment is compatible
                if (dn_use_torenv_table){
                    bool valid_torenv_bool = false;
                    if(dn_use_roulette){ //JDB roulette acceptance
                        valid_torenv_bool = roulette_valid_torenv(frag_combined);
                    } else {
                        valid_torenv_bool = valid_torenv(frag_combined);
                    }

                    //check if *all* torsions are compatible from this pairing
                    if (valid_torenv_bool){
                        // Second check to see whether all of the old torsion environments are
                        // still compatible. This is necessary because of the dummy as wild card
                        // thing
                        if (valid_torenv_multi(frag_combined)){
                            successful_attachments.push_back(make_pair(frag_combined, *pos));
                            //if (Parameter_Reader::verbosity_level() > 1 ) { success_att = true; }
                            success_att = true; 
                        }
                    }
                } else {
                    successful_attachments.push_back(make_pair(frag_combined, *pos));
                    //if (Parameter_Reader::verbosity_level() > 1 ) { success_att = true; }
                    success_att = true; 
                }
            }
        }
  
        if (success_att == true) {
            num_succ_att++;
            num_att_root++;
            num_att_layer++;
        }

        //increment iterator
        pos++;
    }

        //float rate_succ_att = 0.0;
        //rate_succ_att = (float)num_succ_att/selected_fragments.size();
        //std::cout << "rates_successful_att: " << rate_succ_att << std::endl;

    if (Parameter_Reader::verbosity_level() > 1 && selected_fragments.size() > 0) { 
        float rate_succ_att = 0.0;
        rate_succ_att = (float)num_succ_att/selected_fragments.size();
        std::cout << "rates_successful_att: " << rate_succ_att << std::endl;
    }

    if (dn_prune_initial_sample) {
        //sort successful attachments by the scoring function, and then prune by hrmsd
        //done before torsional sampling to account for symmetric attachments
        //activation step done otherwise the pruning step will not work
        //cout << "Number of fragments before torsional sampling, before pruning: " << successful_attachments.size() << endl;
        for (int mol=0; mol<successful_attachments.size(); mol++){
            activate_mol(successful_attachments[mol].first.mol);
            calc_mol_wt(successful_attachments[mol].first.mol);
        }
        //frag_sort(successful_attachments, fragment_vec_pair_sort);
        prune_h_rmsd_vec_pair(successful_attachments);
    }
    //cout << "Number of fragments before torsional sampling, after pruning: " << successful_attachments.size() << endl;
}
void
DN_Build::attach_single_frag( Fragment & layer_frag, int layer_frag_ap, vector <Fragment> & fraglib, 
                                int & selected_fragment, vector<Fragment> & successful_attachments)
{
    //Check all fragments to see if they're valid before attempting an attachment
    // (1) Would create too many scaffolds per layer
    if (fraglib[selected_fragment].aps.size() > 2 &&
        layer_frag.scaffolds_this_layer >= dn_max_scaffolds_per_layer ) {
        //cout << selected_fragment << " # of aps " << fraglib[selected_fragment].aps.size() << endl;
        //cout << "Frag " << selected_fragment << " rejected for too many scaffolds." << endl;
        return;
    }

    // (2) Would give too many current attachment points
    int temp_aps = fraglib[selected_fragment].aps.size() + layer_frag.aps.size() - 2;
    if (temp_aps > dn_max_current_aps) {
        //cout << "Frag " << selected_fragment << " rejected for too many aps." << endl;
        return;
    }

    // (3) Rotatable bond cutoff
    int temp_rbs = layer_frag.mol.rot_bonds + temp_aps;

    // (5) Formal Charge Cutoff
    int temp_fc = fraglib[selected_fragment].mol.formal_charge + layer_frag.mol.formal_charge;
    if (fabs(float(temp_fc)) > (fabs(float(dn_constraint_formal_charge))+0.1)) { 
        //cout << "Frag " << selected_fragment << " rejected for formal charge." << endl;
        return; 
    }

    //If we're here, this is a valid potential attachment - try all attachments
    for (int fraglib_ap = 0; fraglib_ap < fraglib[selected_fragment].aps.size(); fraglib_ap++){
        //check to see if the bond order is compatible
        if (compare_dummy_bonds(layer_frag, layer_frag_ap, fraglib[selected_fragment], fraglib_ap)) {
            //combine the fragments
            Fragment frag_combined = combine_fragments(layer_frag, layer_frag_ap,
                                                        fraglib[selected_fragment], fraglib_ap);
            //check if the new torsion environment from the attachment is compatible
            if (dn_use_torenv_table){
                bool valid_torenv_bool = false;
                if(dn_use_roulette){ //JDB roulette acceptance
                    valid_torenv_bool = roulette_valid_torenv(frag_combined);
                } else {
                    valid_torenv_bool = valid_torenv(frag_combined);
                }
                //check if *all* torsions are compatible from this pairing
                if (valid_torenv_bool){
                    // Second check to see whether all of the old torsion environments are
                    // still compatible. This is necessary because of the dummy as wild card
                    // thing
                    if (valid_torenv_multi(frag_combined)){
                        successful_attachments.push_back(frag_combined);
                    }
                }
            } else {
                successful_attachments.push_back(frag_combined);
            }
        }
    }
    if (dn_prune_initial_sample) {
        //sort successful attachments by the scoring function, and then prune by hrmsd
        //done before torsional sampling to account for symmetric attachments
        //activation step done otherwise the pruning step will not work
        //cout << "Number of fragments before torsional sampling, before pruning: " << successful_attachments.size() << endl;
        for (int mol=0; mol<successful_attachments.size(); mol++){
            activate_mol(successful_attachments[mol].mol);
            calc_mol_wt(successful_attachments[mol].mol);
        }
        prune_h_rmsd_and_mw(successful_attachments);
    }
    //cout << "Number of fragments before torsional sampling, after pruning: " << successful_attachments.size() << endl;
}




// +++++++++++++++++++++++++++++++++++++++++
// Check new atom environments against table of allowable environments
bool
DN_Build::valid_torenv( Fragment & frag )
{

//TODO
// We should think about putting a monte carlo algorithm that checks the bond frequency in deciding whether to accept or reject the new connection

    // Create some temporary objects to compute torsion environments
    Fingerprint tmp_finger;
    pair<string, string> tmp_pair;

    // Will return atom environments only for active atoms. At this point, the only active atoms
    // should be the two that were just attached together. ALL BONDS SHOULD BE ACTIVE.
    tmp_finger.return_torsion_environments( frag.mol, tmp_pair );

    // Check to see if there is a dummy atom in either of the atom  environments in. If there is
    // one, then remember that the wild card will be invoked in compare_atom_environments(). This
    // means that this environment needs to be rechecked later in growth once more fragments are
    // added.
    bool wc_invoked = false;

    if (Parameter_Reader::verbosity_level() > 1) { 
        cout <<"VT: input torenv is " <<tmp_pair.first <<"-" <<tmp_pair.second <<endl;
        cout <<"VT: Starting a new check" <<endl;
        cout <<"VT:   tmp_pair.first.size() = " <<tmp_pair.first.size() <<endl;
    }
    // For every char in the string
    for (int i=0; i<tmp_pair.first.size(); i++){
        // If a dummy atom is encountered (Y in fingerprint.cpp)
        if (Parameter_Reader::verbosity_level() > 1) { 
            cout <<"VT:     tmp_pair.first[i] =" <<tmp_pair.first[i] <<endl; 
        }

        if (tmp_pair.first[i] == 'Y'){
            if (Parameter_Reader::verbosity_level() > 1) { 
                cout <<"VT:       found a match...breaking" <<endl; 
            }
            // Then make this flag true and break the loop
            wc_invoked = true;
            break;
        }
    }

    if (Parameter_Reader::verbosity_level() > 1) { 
        cout <<"VT:   tmp_pair.second.size() = " <<tmp_pair.second.size() <<endl;
    }

    // Iterate over the second string in the same way if needed
    if (!wc_invoked) {
        for (int i=0; i<tmp_pair.second.size(); i++){
            if (Parameter_Reader::verbosity_level() > 1) { 
                cout <<"VT:     tmp_pair.second[i] =" <<tmp_pair.second[i] <<endl; 
            }
            if (tmp_pair.second[i] == 'Y'){
                if (Parameter_Reader::verbosity_level() > 1) { 
                    cout <<"VT:       found a match...breaking" <<endl; 
                }
                wc_invoked = true;
                break;
            }
        }
    }

    if (Parameter_Reader::verbosity_level() > 1) { 
        cout <<"VT:         After both checks, wc_invoked = " <<wc_invoked <<endl <<endl;
    }


    // The returned tmp_pair is alphabetized, as is the torenv vector, so you will only ever have
    // to check tmp_pair.first against the origin_env
    for (int i=0; i<torenv_vector.size(); i++){

        if (Parameter_Reader::verbosity_level() > 1) {
            cout <<"VT: Checking to see if " <<tmp_pair.first <<" is in the fragment library " <<endl;
        }

        // If the first environment of the pair finds a match in the torenv_vector...
        if (compare_atom_environments(torenv_vector[i].origin_env, tmp_pair.first)){

            if (Parameter_Reader::verbosity_level() > 1) { 
                cout <<"VT: Yep it's in there, looking for a match... " <<endl;
            }

            // ...then start checking to see if the second environment in the pair has a 
            // corresponding match
            for (int j=0; j<torenv_vector[i].target_envs.size(); j++){

                if (Parameter_Reader::verbosity_level() > 1) { 
                    cout <<"VT:     Now checking if it is connected to a " <<tmp_pair.second <<endl;
                }

                // If the second env finds a corresponding match
                if (compare_atom_environments(torenv_vector[i].target_envs[j], tmp_pair.second)){

                    // If the wc was invoked, remember this pair for recheck later
                    if (wc_invoked){

                        // For every bond in the molecule
                        for (int k=0; k<frag.mol.num_bonds; k++){

                            // If the target and origin atoms are both active, then this is the torenv
                            // that we are currently checking
                            if (frag.mol.atom_active_flags[ frag.mol.bonds_origin_atom[k] ] &&
                                frag.mol.atom_active_flags[ frag.mol.bonds_target_atom[k] ]   ) {

                                // Record those atom indices for recheck later
                                pair <int, int> temp_pair;
                                temp_pair.first = frag.mol.bonds_origin_atom[k];
                                temp_pair.second = frag.mol.bonds_target_atom[k];
                                frag.torenv_recheck_indices.push_back(temp_pair);

                                // Once the correct bond is identified, don't need to check the rest
                                break;
                            }
                        }
                    }

                    // Return TRUE as this is a valid pair
                    return true;
                }
            }
        } 
    }

    // In the future, think about rewriting the torenv data structure to not be alphabitzed / order dependent.
    // Then we wouldn't need this second loop, but it might be kinda memory intensive.
    //
    //BCF with wildcards alphabetization might not be preserved, need to check other ordering of torenv
    for (int i=0; i<torenv_vector.size(); i++){

        if (Parameter_Reader::verbosity_level() > 1) { 
            cout <<"VT: Checking to see if " <<tmp_pair.second <<" is in the fragment library " <<endl;
        }

        if (compare_atom_environments(torenv_vector[i].origin_env, tmp_pair.second)){

            if (Parameter_Reader::verbosity_level() > 1) { 
                cout <<"VT: Yep it's in there, looking for a match... " <<endl;
            }

            for (int j=0; j<torenv_vector[i].target_envs.size(); j++){

                if (Parameter_Reader::verbosity_level() > 1) { 
                    cout <<"VT:     Now checking if it is connected to a " <<tmp_pair.first <<endl;
                }

                if (compare_atom_environments(torenv_vector[i].target_envs[j], tmp_pair.first)){

                    if (wc_invoked){

                        for (int k=0; k<frag.mol.num_bonds; k++){

                            if (frag.mol.atom_active_flags[ frag.mol.bonds_origin_atom[k] ] &&
                                frag.mol.atom_active_flags[ frag.mol.bonds_target_atom[k] ]   ) {

                                pair <int, int> temp_pair;
                                temp_pair.second = frag.mol.bonds_origin_atom[k];
                                temp_pair.first = frag.mol.bonds_target_atom[k];
                                frag.torenv_recheck_indices.push_back(temp_pair);

                                break;
                            }
                        }
                    }
                    return true;
                }
            }
        }
    }
    // end BCF part
    // If the end is reached, the newly created pair was not found in the torenv library, so return
    // false and do not keep the molecule
    return false;

} // end DN_Build::valid_torenv()




// +++++++++++++++++++++++++++++++++++++++++
// Check new atom environments against table of allowable environments
// and accept based on roulette probability
bool
DN_Build::roulette_valid_torenv( Fragment & frag )
{
    // Create some temporary objects to compute torsion environments
    Fingerprint tmp_finger;
    pair<string, string> tmp_pair;
    // Will return atom environments only for active atoms. At this point, the only active atoms
    // should be the two that were just attached together. ALL BONDS SHOULD BE ACTIVE.
    tmp_finger.return_torsion_environments( frag.mol, tmp_pair );
    
    // Check to see if there is a dummy atom in either of the atom environments in. If there is
    // one, then remember that the wild card will be invoked in compare_atom_environments(). This
    // means that this environment needs to be rechecked later in growth once more fragments are
    // added.
    bool wc_invoked = false;
    
    if (Parameter_Reader::verbosity_level() > 1) {
        cout <<"VT: input torenv is " <<tmp_pair.first <<"-" <<tmp_pair.second <<endl;
        cout <<"VT: Starting a new check" <<endl;
        cout <<"VT:   tmp_pair.first.size() = " <<tmp_pair.first.size() <<endl;
    }
    // For every char in the string
    for (int i=0; i<tmp_pair.first.size(); i++){
        // If a dummy atom is encountered (Y in fingerprint.cpp)
        if (Parameter_Reader::verbosity_level() > 1) {
            cout <<"VT:     tmp_pair.first[i] =" <<tmp_pair.first[i] <<endl;
        }
        
        if (tmp_pair.first[i] == 'Y'){
            if (Parameter_Reader::verbosity_level() > 1) { 
                cout <<"VT:       found a match...breaking" <<endl; 
            }
            // Then make this flag true and break the loop
            wc_invoked = true;
            break;
        }
    }
    
    if (Parameter_Reader::verbosity_level() > 1) {
        cout <<"VT:   tmp_pair.second.size() = " <<tmp_pair.second.size() <<endl;
    }
    
    // Iterate over the second string in the same way if needed
    if (!wc_invoked) {
        for (int i=0; i<tmp_pair.second.size(); i++){
            if (Parameter_Reader::verbosity_level() > 1) {
                cout <<"VT:     tmp_pair.second[i] =" <<tmp_pair.second[i] <<endl;
            }
            if (tmp_pair.second[i] == 'Y'){
                if (Parameter_Reader::verbosity_level() > 1) {
                    cout <<"VT:       found a match...breaking" <<endl; 
                }
                wc_invoked = true;
                break;
            }
        }
    }
    
    if (Parameter_Reader::verbosity_level() > 1) {
        cout <<"VT:         After both checks, wc_invoked = " <<wc_invoked <<endl <<endl;
    }
    
    
    // The returned tmp_pair is alphabetized, as is the torenv vector, so you will only ever have
    // to check tmp_pair.first against the origin_env
    bool first_tors_check=false, second_tors_check=false;
    tmp_first_check=false;
    tmp_second_check=false;
    for (int i=0; i<torenv_vector.size(); i++){
        
        if (Parameter_Reader::verbosity_level() > 1) {
            cout <<"VT: Checking to see if " <<tmp_pair.first <<" is in the fragment library " <<endl;
        }
        
        // If the first environment of the pair finds a match in the torenv_vector...
        tmp_first_check=true;
        if (compare_atom_environments(torenv_vector[i].origin_env, tmp_pair.first)){
            
            if (Parameter_Reader::verbosity_level() > 1) {
                cout <<"VT: Yep it's in there, looking for a match... " <<endl;
            }
            
            // ...then start checking to see if the second environment in the pair has a
            // corresponding match
            tmp_first_check=false;
            tmp_second_check=true;
            for (int j=0; j<torenv_vector[i].target_envs.size(); j++){
                
                if (Parameter_Reader::verbosity_level() > 1) {
                    cout <<"VT:     Now checking if it is connected to a " <<tmp_pair.second <<endl;
                }
                
                // If the second env finds a corresponding match
                if (compare_atom_environments(torenv_vector[i].target_envs[j], tmp_pair.second)){
                    //cout << "Corresponding Match Output forward: " << torenv_vector[i].target_envs[j] <<  endl << endl;
                    // If the wc was invoked, remember this pair for recheck later
                    if (wc_invoked){
                        
                        // For every bond in the molecule
                        for (int k=0; k<frag.mol.num_bonds; k++){
                            
                            // If the target and origin atoms are both active, then this is the torenv
                            // that we are currently checking
                            if (frag.mol.atom_active_flags[ frag.mol.bonds_origin_atom[k] ] &&
                                frag.mol.atom_active_flags[ frag.mol.bonds_target_atom[k] ]   ) {
                                
                                // Record those atom indices for recheck later
                                pair <int, int> temp_pair;
                                temp_pair.first = frag.mol.bonds_origin_atom[k];
                                temp_pair.second = frag.mol.bonds_target_atom[k];
                                frag.torenv_recheck_indices.push_back(temp_pair);
                                // Once the correct bond is identified, don't need to check the rest
                                break;
                            }
                        }
                    }
                    //if this point is reached, this is a valid torsion - set variable to true and break out
                    first_tors_check=true;
                    break;
                }
            }
        }
        //if the above variable is true, no further iterations of loop are necessary. break out.
        if(first_tors_check){
            break;
        }
    }
    // In the future, think about rewriting the torenv data structure to not be alphabitzed / order dependent.
    // Then we wouldn't need this second loop, but it might be kinda memory intensive.
    //
    //BCF with wildcards alphabetization might not be preserved, need to check other ordering of torenv
    if (!first_tors_check){
        //JDB reset variables for passing into compare_atom_environments
        tmp_first_check=false;
        tmp_second_check=false;
        for (int i=0; i<torenv_vector.size(); i++){
            
            if (Parameter_Reader::verbosity_level() > 1) {
                cout <<"VT: Checking to see if " <<tmp_pair.second <<" is in the fragment library " <<endl;
            }
            
            tmp_second_check=true;
            if (compare_atom_environments(torenv_vector[i].origin_env, tmp_pair.second)){
                
                if (Parameter_Reader::verbosity_level() > 1) {
                    cout <<"VT: Yep it's in there, looking for a match... " <<endl;
                }
                
                for (int j=0; j<torenv_vector[i].target_envs.size(); j++){
                    
                    if (Parameter_Reader::verbosity_level() > 1) {
                        cout <<"VT:     Now checking if it is connected to a " <<tmp_pair.first <<endl;
                    }
                    tmp_second_check=false;
                    tmp_first_check=true;
                    if (compare_atom_environments(torenv_vector[i].target_envs[j], tmp_pair.first)){
                        
                        if (wc_invoked){
                            
                            for (int k=0; k<frag.mol.num_bonds; k++){
                                
                                if (frag.mol.atom_active_flags[ frag.mol.bonds_origin_atom[k] ] &&
                                    frag.mol.atom_active_flags[ frag.mol.bonds_target_atom[k] ]   ) {
                                    
                                    pair <int, int> temp_pair;
                                    temp_pair.second = frag.mol.bonds_origin_atom[k];
                                    temp_pair.first = frag.mol.bonds_target_atom[k];
                                    frag.torenv_recheck_indices.push_back(temp_pair);
                                    
                                    break;
                                }
                            }
                        }
                        //if this point is reached, this is a valid torsion - set variable to true and break out
                        second_tors_check=true;
                        break;
                    }
                }
            }
            //if the above variable is true, no further iterations of loop are necessary. break out.
            if(second_tors_check){
                break;
            }
        }
    }
    std::string roulette_compare_environment;
    // end BCF part
    
    //if neither check returned true, then the torsion is invalid - return false
    if ((!first_tors_check) && (!second_tors_check)){
        return false;
        
        //if the first torsion ordering was found, enter this statement
    } else if (first_tors_check){
        //if a dummy atom was found, use the variable set from compare_atom_environments
        if (wc_invoked){
            //if statement for if both halves of the torsion contain a dummy atom
            if((tmp_pair.first.find("Y") != string::npos) && (tmp_pair.second.find("Y") != string::npos)){
                if(Parameter_Reader::verbosity_level() > 1)
                {
                    cout << "WC_Invoked, FirstTors; Replace string1&2    " << roulette_dummy_replace_1 << "   " << roulette_dummy_replace_2 << endl;
                }
                roulette_compare_environment=roulette_dummy_replace_1 + "-" + roulette_dummy_replace_2;
                //if statement for if only the first half contains a dummy atom
            } else if (tmp_pair.first.find("Y") != string::npos){
                if(Parameter_Reader::verbosity_level() > 1)
                {
                    cout << "WC_Invoked, FirstTors; Replace string1:     " << roulette_dummy_replace_1 << endl;
                    cout << "Old Torsion: " << tmp_pair.first + "-" + tmp_pair.second << endl;
                }
                roulette_compare_environment=roulette_dummy_replace_1 + "-" + tmp_pair.second;
                //if statement for if only the second half contains a dummy atom
            } else {
                if(Parameter_Reader::verbosity_level() > 1)
                {
                    cout << "WC_Invoked, FirstTors; Replace string2:     " << roulette_dummy_replace_2 << endl;
                    cout << "Old Torsion: " << tmp_pair.first + "-" + tmp_pair.second << endl;
                }
                roulette_compare_environment=tmp_pair.first + "-" + roulette_dummy_replace_2;
            }
            //if there is no dummy atom, just assign each half-torsion as is
        } else {
            roulette_compare_environment=tmp_pair.first + "-" + tmp_pair.second;
        }
        
        //if the second torsion ordering was found, enter this statement
    } else if (second_tors_check){
        //if a dummy atom was found, use the variables set from compare_atom_environments
        if (wc_invoked) {
            //if statement for if both halves of the torsion contain a dummy atom
            if((tmp_pair.first.find("Y") != string::npos) && (tmp_pair.second.find("Y") != string::npos)){
                roulette_compare_environment=roulette_dummy_replace_2 + "-" + roulette_dummy_replace_1;
                if(Parameter_Reader::verbosity_level() > 1)
                {
                    cout << "WC_Invoked, SecondTors; Replace string2&1    " << roulette_dummy_replace_2 << "    " << roulette_dummy_replace_1 << endl;
                }
                //if statement for if only the first half contains a dummy atom
            } else if(tmp_pair.first.find("Y") != string::npos){
                if(Parameter_Reader::verbosity_level() > 1)
                {
                    cout << "WC_Invoked, SecondTors; Replace string2:     " << roulette_dummy_replace_1 << endl;
                    cout << "Old Torsion: " << tmp_pair.second + "-" + tmp_pair.first;
                }
                roulette_compare_environment=tmp_pair.second + "-" + roulette_dummy_replace_1;
                //if statement for if only the second half contains a dummy atom
            } else {
                if(Parameter_Reader::verbosity_level() > 1)
                {
                    cout << "WC_Invoked, SecondTors; Replace string1:     " << roulette_dummy_replace_2 << endl;
                    cout << "Old Torsion: " << tmp_pair.second + "-" + tmp_pair.first << endl;
                }
                roulette_compare_environment=roulette_dummy_replace_2 + "-" + tmp_pair.first;
            }
            //if there is no dummy atom, just assign each half-torsion as is
        } else {
            roulette_compare_environment = tmp_pair.second + "-" + tmp_pair.first;
        }
    }
    
    //initialize variables - frequency is from the table, compare_frequency is a random number
    double roulette_frequency = 1.0; 
    double roulette_compare_frequency = ( (double) rand() / RAND_MAX) ;

    
    //checks for torsion in the torsion table, then pulls the frequency when it's found
    for(int i=0; i<torsion_vect.size();++i){
        if(torsion_vect[i] == roulette_compare_environment){
            roulette_frequency=roulette_vect[i];
            break;
        }
    }
    
    /*
     |----------|-----|---|-|
     0          1 ^   2   3
     ^ == random number selection point
     
     Torsion indexes are defined by the right wall of each bin on the number line.
     
     Case 1a: If torsion is bin 1, would accept as frequency < random number.
     Case 1b: If torsion is bin 0, accept automatically.
     If Case 1a and 1b are false, proceed to Operation 1.5
     
     Operation 1.5: Obtain the index of the torsion's right wall.
     
     Case 2: Perform check to see if the left of wall of selected torsion bin is less than the random number
     (ie, if the random number falls within the selected torsion's bin)
     Accept if the random number is in the bin, otherwise pass to Case 3.
     
     Case 3: The random number is not greater than the bin's right wall, nor falls within the bin.
     Reject torsion.
     */
    
    //CASE 1
    if((roulette_frequency <= roulette_compare_frequency) || (roulette_frequency == roulette_bin_vect[0])){
        
        //outputs the statistics for this selection.
        if(Parameter_Reader::verbosity_level() > 1)
        {
            cout << "FREQ < COMP   " << roulette_frequency << " < " << roulette_compare_frequency << endl;
            if(first_tors_check){
                cout << "Torsion Forward: " << roulette_compare_environment << "\t\tAccepted1." << endl;
            } else {
                cout << "Torsion Backward: " << roulette_compare_environment << "\t\tAccepted1." <<endl;
            }
            cout << "Roulette Frequency: " << roulette_frequency << endl;
            cout << "Roulette comparison: " << roulette_compare_frequency << endl;
            cout << "__________________" << endl << endl;
        }
        return true; //accept torsion
        
        //OPERATION 1.5
    } else { //if the frequency isn't directly less than the bin's higher end
        //initializes a variable that gets assigned the exact position the frequency is in the bin list
        int roulette_bin_index = roulette_bin_vect.size()-1;
        for(int i=0; i<roulette_bin_vect.size(); i++)
        {
            if(roulette_bin_vect[i] == roulette_frequency)
            {
                roulette_bin_index = i;
                break;
            }
        }
        
        //CASE 2
        if((roulette_bin_vect[roulette_bin_index-1] < roulette_compare_frequency) && (roulette_compare_frequency  <= roulette_frequency))
        {
            
            //outputs statistics for this selection.
            if(Parameter_Reader::verbosity_level() > 1)
            {
                cout << "FREQ < COMP   " << roulette_bin_vect[roulette_bin_index-1] << " < " << roulette_compare_frequency << " < " << roulette_frequency << endl;
                if(first_tors_check){
                    cout << "Torsion Forward: " << roulette_compare_environment << "\t\tAccepted3." << endl;
                } else {
                    cout << "Torsion Backward: " << roulette_compare_environment << "\t\tAccepted3." <<endl;
                }
                cout << "Roulette Frequency: " << roulette_frequency << endl;
                cout << "Roulette comparison: " << roulette_compare_frequency << endl;
                cout << "__________________" << endl << endl;
            }
            return true; //accept torsion
            
            //CASE 3
        } else {
            
            //outputs statistics for this selection
            if(Parameter_Reader::verbosity_level() > 1)
            {
                cout << "FREQ > COMP   " << roulette_frequency << " > " << roulette_compare_frequency << endl;
                cout << "Torsion: " << roulette_compare_environment << "\t\tRejected." << endl;
                cout << "Roulette Frequency: " << roulette_frequency << endl;
                cout << "Roulette comparison: " << roulette_compare_frequency << endl;
                cout << "___________________" << endl << endl;
            }
            
            return false; //reject torsion
        }
    }
    // If the end is reached, the newly created pair was not found in the torenv library, so return
    // false and do not keep the molecule
    return false;
    
} // end DN_Build::roulette_valid_torenv()



// +++++++++++++++++++++++++++++++++++++++++
// Check new atom environments against table of allowable environments
bool
DN_Build::valid_torenv_multi( Fragment & frag )
{

    // If there is nothing to recheck, return
    if (frag.torenv_recheck_indices.size() == 0){
        return true;
    }

    // Create some temporary objects to compute torsion environments
    Fingerprint tmp_finger;
    pair<string, string> tmp_pair;

    // Some variables to track how many of the rechecks pass the test
    int passed_check = 0;
    int expected = frag.torenv_recheck_indices.size();

    // For every atom environment pair that you want to recheck
    for (int i=0; i<frag.torenv_recheck_indices.size(); i++){
    
        // First inactivate all of the atoms
        for (int j=0; j<frag.mol.num_atoms; j++){
                frag.mol.atom_active_flags[j] = false;
        }
    
        // Then activate just the ones around this bond
        frag.mol.atom_active_flags[ frag.torenv_recheck_indices[i].first ] = true;
        frag.mol.atom_active_flags[ frag.torenv_recheck_indices[i].second ] = true;
    
        // Will return atom environments only for active atoms. At this point, the only active atoms
        // should be the two that were just attached together. ALL BONDS SHOULD BE ACTIVE.
        tmp_finger.return_torsion_environments( frag.mol, tmp_pair );
    
        bool double_break = false;
        // The returned tmp_pair is alphabetized, as is the torenv vector, so you will only ever have
        // to check tmp_pair.first against the origin_env
        for (int j=0; j<torenv_vector.size(); j++){
            if (double_break) {break;}
            if (Parameter_Reader::verbosity_level() > 1) { 
                cout <<"VTM: Checking to see if " <<tmp_pair.first <<" is in the fragment library " <<endl;
            }
    
            // If the first environment of the pair finds a match in the torenv_vector...
            if (compare_atom_environments(torenv_vector[j].origin_env, tmp_pair.first)){
                if (Parameter_Reader::verbosity_level() > 1) {
                    cout <<"VTM: Yep it's in there, looking for a match... " <<endl; 
                }
    
                // ...then start checking to see if the second environment in the pair has a 
                // corresponding match
                for (int k=0; k<torenv_vector[j].target_envs.size(); k++){
                    if (Parameter_Reader::verbosity_level() > 1) {
                        cout <<"VTM:     Now checking if it is connected to a " <<tmp_pair.second <<endl; 
                    }
    
                    // If the second env finds a corresponding match
                    if (compare_atom_environments(torenv_vector[j].target_envs[k], tmp_pair.second)){
    
                        // Then this one passed the check
                        passed_check++;
                        double_break = true;
                        break;
                    }
                }
            }
        }
//
// BCF Part
    // In the future, think about rewriting the torenv data structure to not be alphabitzed / order dependent.
    // Then we wouldn't need this second loop, but it might be kinda memory intensive.
    //
        if (!double_break) {
            for (int j=0; j<torenv_vector.size(); j++){
                if (double_break) {break;}
                if (Parameter_Reader::verbosity_level() > 1) {
                    cout <<"VTM: Checking to see if " <<tmp_pair.second <<" is in the fragment library " <<endl;
                }

                if (compare_atom_environments(torenv_vector[j].origin_env, tmp_pair.second)){
                    if (Parameter_Reader::verbosity_level() > 1) {
                        cout <<"VTM: Yep it's in there, looking for a match... " <<endl; 
                    }

                    for (int k=0; k<torenv_vector[j].target_envs.size(); k++){
                        if (Parameter_Reader::verbosity_level() > 1) {
                            cout <<"VTM:     Now checking if it is connected to a " <<tmp_pair.first <<endl; 
                        }

                        if (compare_atom_environments(torenv_vector[j].target_envs[k], tmp_pair.first)){
  
                            passed_check++;
                            double_break = true;
                            break;
                        }
                    }
                }
            }
        } //if !doublebreak we haven't found a match, switch the order and try again
// end BCF part
    }

    // If they all passed the check, return true, otherwise return false
    if (passed_check == expected){ return true; }
    else { return false; }

} // end DN_Build::valid_torenv_multi()



// +++++++++++++++++++++++++++++++++++++++++
// Function to compare two atom environments assuming that dummy atom Y is a wildcard
bool
DN_Build::compare_atom_environments( string s1, string s2 )
{
    // s1 is the reference string from the vector of torenvs
    // s2 is a candidate that we are checking

    if (Parameter_Reader::verbosity_level() > 1){
        cout <<"CAE: Comparing s1 = " <<s1 <<"\t and s2 = " <<s2 <<endl; 
    }

    // If they are exactly the same, then it's a match
    if (s1 == s2){
        if (Parameter_Reader::verbosity_level() > 1){
            cout <<"CAE: Matched s1 and s2: Strings are identical" <<endl <<endl; 
        }
        return true;
    }

    // The environments should not match if strings are not same length
    if (s1.size() != s2.size()) {
        return false;
    }

    // Remove carets (from s1 and s2) and dummies (from s2 only - s1 should never have dummies)
    int count = 0;
    int count_dummie = 0;
    string s1_new = "";
    string s2_new = "";

    // Remove carets from s1
    for (int i=0; i<s1.size(); i++){
        if (s1[i] != '^' && s1[i] != '#'){
            s1_new.push_back(s1[i]);
        } 
    }

    // Remove carets from s2 and count and remove dummies in s2
    for (int i=0; i<s2.size(); i++){
        if (s2[i] == 'Y'){
            count_dummie++;

        } else if (s2[i] != '^' && s2[i] != '#'){
            s2_new.push_back(s2[i]);
        } 
    }


    if (Parameter_Reader::verbosity_level() > 1){
        cout <<"CAE: Comparing s1_new = " <<s1_new <<"\t and s2_new = " <<s2_new <<endl; 
    }

    // Compare contents of s1 and s2
    for (int i=0; i<s1_new.size(); i+=2){
        if ( (s1_new.c_str()[i] == s2_new.c_str()[count]) && 
             (s1_new.c_str()[i+1] == s2_new.c_str()[count+1]) ){
            if (Parameter_Reader::verbosity_level() > 1) {
                cout <<"CAE: Matched  s1_new[" <<i <<"] = " <<s1_new.c_str()[i] 
                     <<" and s2_new[" <<count <<"] = " <<s2_new.c_str()[count] 
                     <<endl;
                cout <<"CAE: Matched  s1_new[" <<i+1 <<"] = " <<s1_new.c_str()[i+1] 
                     <<" and s2_new[" <<count+1 <<"] = " <<s2_new.c_str()[count+1] 
                     <<endl;
            }
            count+=2;
        } else if (count_dummie > 0){
            if (Parameter_Reader::verbosity_level() > 1){
                cout <<"CAE: Matched  s1_new[" <<i <<"] = " <<s1_new.c_str()[i] <<" and dummy" <<endl;
                cout <<"CAE: Matched  s1_new[" <<i+1 <<"] = " <<s1_new.c_str()[i+1] <<" and dummy" <<endl;
            }
            count_dummie--;
        } else {
            if (Parameter_Reader::verbosity_level() > 1){
                cout <<"CAE: Did not match s1 =\t" <<s1 <<"\tand s2 = " <<s2 <<endl <<endl;
            }
            return false;
        }
    }

    if (Parameter_Reader::verbosity_level() > 1) { 
        cout <<"CAE: Matched s1 =\t" <<s1 <<"\tand s2 = " <<s2 <<endl <<endl;
    }
    //JDB roul logic check for which variable is being passed, and appropriate variable assignment
    if (dn_use_roulette){
        if(tmp_first_check){
            roulette_dummy_replace_1=s1;
            if (Parameter_Reader::verbosity_level() > 1)
            {
                cout << "CAE: Roulette_Dummy_Replace_1 = " << roulette_dummy_replace_1 << endl << endl;
            }
        } else if(tmp_second_check) {
            roulette_dummy_replace_2=s1;
            if (Parameter_Reader::verbosity_level() > 1)
            {
                cout << "CAE: Roulette_Dummy_Replace_2 = " << roulette_dummy_replace_2 << endl << endl;
            }
        }
    }
    return true;

} // end DN_Build::compare_atom_environments()



// +++++++++++++++++++++++++++++++++++++++++
// Check to see if two dummy bond connections are the same bond order
bool
DN_Build::compare_dummy_bonds( Fragment & frag1, int attachment_point_1,
                             Fragment & frag2, int attachment_point_2 )
{

    //assign some variables from the input fragments
    //frag1 is the fragment we are attaching TO, and frag2
    //is the fragment we are trying to attach.
    int dummy1 = frag1.aps[attachment_point_1].dummy_atom;
    int dummy2 = frag2.aps[attachment_point_2].dummy_atom;
    int heavy1 = frag1.aps[attachment_point_1].heavy_atom;
    int heavy2 = frag2.aps[attachment_point_2].heavy_atom;


    // Check to see if bond order connecting heavy1 and dummy1 is same as bond order connecting
    // heavy2 and dummy2
    string bond_type1;
    string bond_type2;

    // Loop over all bonds in frag1
    for (int i=0; i<frag1.mol.num_bonds; i++){
        if (frag1.mol.bonds_origin_atom[i] == dummy1 ||
            frag1.mol.bonds_target_atom[i] == dummy1   ){
            if (frag1.mol.bonds_origin_atom[i] == heavy1 ||
                frag1.mol.bonds_target_atom[i] == heavy1   ){

                // This is the bond you're looking for...
                bond_type1 = frag1.mol.bond_types[i];
                break;
            }
        }
    }

    // Loop over all bonds in frag2
    for (int i=0; i<frag2.mol.num_bonds; i++){
        if (frag2.mol.bonds_origin_atom[i] == dummy2 ||
            frag2.mol.bonds_target_atom[i] == dummy2   ){
            if (frag2.mol.bonds_origin_atom[i] == heavy2 ||
                frag2.mol.bonds_target_atom[i] == heavy2   ){

                // This is the bond you're looking for...
                bond_type2 = frag2.mol.bond_types[i];
                break;
            }
        }
    }

    // If they match, return true
    return (bond_type1 == bond_type2);

} // end DN_Build::compare_dummy_bonds();

// +++++++++++++++++++++++++++++++++++++++++
// Given a fragment and an empty vector of fragments, sample torsions for the most recent addition
// and return the results in the vector.
void
DN_Build::sample_minimized_torsions( Fragment & frag1, vector <Fragment> & list_of_frags, 
                                     Master_Score & score, Simplex_Minimizer & simplex, 
                                     AMBER_TYPER & typer )
{
    // Activate all atoms and bonds prior to any scoring
    activate_mol(frag1.mol);

    // Create a temporary vector to populate with different torsions
    vector <Fragment> tmp_list;

    // This will be true if energy is calculated correctly
    bool valid_orient = false;

    // Declare a temporary fingerprint object and compute atom environments (necessary for Gasteiger)
    Fingerprint temp_finger;
    for (int i=0; i<frag1.mol.num_atoms; i++){
        frag1.mol.atom_envs[i] = temp_finger.return_environment(frag1.mol, i);
    }
                
    // Prepare the molecule using the amber_typer to assign bond types to each bond
    frag1.mol.prepare_molecule();
    typer.prepare_molecule(frag1.mol, true, score.use_chem, score.use_ph4, score.use_volume);

    // Compute the charges, saving them on the mol object
    float total_charges = 0;
    total_charges = compute_gast_charges(frag1.mol);

    //BCF  04/28/16 update values, pruning is now done in sample_fraglib_graph
    //after gasteiger charges were computed
    calc_mol_wt(frag1.mol);
    calc_rot_bonds(frag1.mol);
    calc_formal_charge(frag1.mol);

/*  BCF this check of mol_wt formal_charge and rot_bonds is now completed in sample_fraglib_graph 
    if (frag1.mol.mol_wt > dn_constraint_mol_wt ||
        frag1.mol.rot_bonds > dn_constraint_rot_bon ||
        fabs(frag1.mol.formal_charge) > (fabs(dn_constraint_formal_charge)+0.1) ){return;}
*/
    // If the total_charges exceeds user-defined cut off, forget this molecule
    // if (total_charges > dn_constraint_charge_groups){ return; }

    // Sample torsions for frag1, move torsions onto tmp_list
    if (dn_sample_torsions){
        frag_torsion_drive(frag1, tmp_list);
    } else {
        tmp_list.push_back(frag1);
    }
    // For all of the torsions that were just generated for the most recent fragment addition,
    for (int i=0; i<tmp_list.size(); i++){

        // Initially reset all of the values in bond_tors_vectors to -1 (bond direction does not
        // matter)
        bond_tors_vectors.clear();
        bond_tors_vectors.resize(tmp_list[i].mol.num_bonds, -1);

        // However, set bonds betweeen fragments to rotatable, adding degrees of freedom to the
        // minimizer. The value at rotatable bond index i in bond_tors_vectors should be the
        // atom_index of one atom involved in the bond, and it should be the atom that is *closer*
        // to the anchor fragment. amber_bt_id just says whether the bond is rotatable (positive
        // integer) or not (-1).
        for (int j=0; j<tmp_list[i].mol.num_bonds; j++){
            int res_origin = atoi(tmp_list[i].mol.atom_residue_numbers[tmp_list[i].mol.bonds_origin_atom[j]].c_str());
            int res_target = atoi(tmp_list[i].mol.atom_residue_numbers[tmp_list[i].mol.bonds_target_atom[j]].c_str());

            if (res_origin != res_target){

                if (res_origin < res_target){ bond_tors_vectors[j] = tmp_list[i].mol.bonds_origin_atom[j]; }
                else                        { bond_tors_vectors[j] = tmp_list[i].mol.bonds_target_atom[j]; }
    
                tmp_list[i].mol.amber_bt_id[j] = 1;
            } else {
               //BCF fix for dn refinement with rotatable bonds on user given anchors
                tmp_list[i].mol.amber_bt_id[j] = frag1.mol.amber_bt_id[j];
                //tmp_list[i].mol.amber_bt_id[j] = -1;
            }
        } 

        // Prepare for internal energy calculation for this fragment
        if (use_internal_energy){ prepare_internal_energy( tmp_list[i], score ); }
        // Copy the pre-min frag for the growth tree
        if (dn_write_growth_trees){
            tmp_list[i].frag_growth_tree.push_back(tmp_list[i].mol);
        }
        // Flexible minimization
        simplex.minimize_flexible_growth(tmp_list[i].mol, score, bond_tors_vectors);
        // Compute internal energy and primary score, store it in the dockmol
        if (score.use_primary_score) {
            valid_orient = score.compute_primary_score(tmp_list[i].mol);
        } else {
            valid_orient = true;
            tmp_list[i].mol.current_score = 0;
            tmp_list[i].mol.internal_energy = 0;
        }

        // If the energy or if the internal energy is greater than the cutoff, erase that conformer
        if ( tmp_list[i].mol.current_score > dn_pruning_conformer_score_cutoff  || 
             !valid_orient ){
            tmp_list[i].used = true;

// BCF part
            //tmp_list.erase(tmp_list.begin()+i);
            //i--;
        } else {
         
              if (use_internal_energy) {
                      if ( tmp_list[i].mol.internal_energy > ie_cutoff ) {
                          tmp_list[i].used = true;
                          //tmp_list.erase(tmp_list.begin()+i);
                          //i--;
                      }
              }
         }
// end BCF part
        bond_tors_vectors.clear();
    }
    // JDB 2022.11.16
    // Legacy pruning to keep only the best molecule if no dummy atoms are present
    bool dummy_found = true; // Defined as true so if legacy pruning is on, it'll use 6.9 behavior
    if (!dn_legacy_prune_complete){
        dummy_found = false;
        //Before RMSD pruning, check one mol to see if a dummy atom present
        for (int i=0; i<tmp_list[0].mol.num_atoms; i++) {
            if (tmp_list[0].mol.atom_types[i] == "Du"){
                dummy_found = true;
                break;
            }
        }
    }
    //If a dummy atom was found, this molecule is incomplete - cluster by
    //RMSD.
    if (dummy_found == true) {
        // Cluster post-minimization torsions by RMSD and remove those under a certain cutoff 
        float rmsd;
        // For every fragment in tmp_list
        for (int i=0; i<tmp_list.size(); i++){

            // If the fragment is not currently 'used'
            if (!tmp_list[i].used){

                // Then check the next fragment
                for (int j=i+1; j<tmp_list.size(); j++){

                    // If it is also not 'used'
                    if (!tmp_list[j].used){

                        // Then calculate the rmsd between the two
                        rmsd = calc_fragment_rmsd(tmp_list[i], tmp_list[j]);
                        rmsd = MAX(rmsd, 0.001);

                        // If the rmsd is below the specified cut-off, they are in the same 'cluster'
                        if ((float) j / rmsd > dn_pruning_clustering_cutoff) {
                            tmp_list[j].used = true;
                        }
                    }        
                }
            }
        }
    // TODO the above comparison is checking heavy atom RMSD for the WHOLE 
    // MOLECULE. In A&G, they only check heavy atom RMSD for the newly-added 
    // segment. We need to reconcile that

        int count_stuff = 0;
        // At this point, only cluster heads will not be flagged 'used', 
        //add those to list_of_frags
        for (int i=0; i<tmp_list.size(); i++){
            if (!tmp_list[i].used){
                list_of_frags.push_back(tmp_list[i]);
                count_stuff++;
            }
        }
    //If a dummy atom was not found, this molecule is complete. Only keep the
    // best one by score.
    } else if (dummy_found == false) {
        //if (Parameter_Reader::verbosity_level() > 1) {
        if (verbose) {
            cout << "Molecule complete! Keeping best in vector by score." << endl;
        }
        
        frag_sort(tmp_list, fragment_sort);
        list_of_frags.push_back(tmp_list[0]);
    }


    // Clear the memory of tmp_list
    tmp_list.clear();

    return;

} // end DN_Build::sample_minimized_torsions()

// +++++++++++++++++++++++++++++++++++++++++
// Given a fragment and an empty vector of fragments, sample torsions for the most recent addition
// and return the results in the vector.
float
DN_Build::sample_minimized_torsions_iso( Fragment & frag1, vector <Fragment> & list_of_frags, 
                                     Master_Score & score, Simplex_Minimizer & simplex, 
                                     AMBER_TYPER & typer )
{
    // Activate all atoms and bonds prior to any scoring
    activate_mol(frag1.mol);

    // Create a temporary vector to populate with different torsions
    vector <Fragment> tmp_list;
    std::vector<float> tmp_rot_list;

    // This will be true if energy is calculated correctly
    bool valid_orient = false;

    // Declare a temporary fingerprint object and compute atom environments (necessary for Gasteiger)
    Fingerprint temp_finger;
    for (int i=0; i<frag1.mol.num_atoms; i++){
        frag1.mol.atom_envs[i] = temp_finger.return_environment(frag1.mol, i);
    }
                
    // Prepare the molecule using the amber_typer to assign bond types to each bond
    frag1.mol.prepare_molecule();
    typer.prepare_molecule(frag1.mol, true, score.use_chem, score.use_ph4, score.use_volume);

    // Compute the charges, saving them on the mol object
    float total_charges = 0;
    total_charges = compute_gast_charges(frag1.mol);

    //BCF  04/28/16 update values, pruning is now done in sample_fraglib_graph
    //after gasteiger charges were computed
    calc_mol_wt(frag1.mol);
    calc_rot_bonds(frag1.mol);
    calc_formal_charge(frag1.mol);

/*  BCF this check of mol_wt formal_charge and rot_bonds is now completed in sample_fraglib_graph 
    if (frag1.mol.mol_wt > dn_constraint_mol_wt ||
        frag1.mol.rot_bonds > dn_constraint_rot_bon ||
        fabs(frag1.mol.formal_charge) > (fabs(dn_constraint_formal_charge)+0.1) ){return;}
*/
    // If the total_charges exceeds user-defined cut off, forget this molecule
    // if (total_charges > dn_constraint_charge_groups){ return; }

    // Sample torsions for frag1, move torsions onto tmp_list
    if (dn_sample_torsions){
        // return list of degrees
        tmp_rot_list = frag_torsion_drive_iso(frag1, tmp_list);
    } else {

        // if not just fill tmp_list with a frags
        // if you don't want to sample torsions 
        tmp_list.push_back(frag1);
    }
    // For all of the torsions that were just generated for the most recent fragment addition,
    for (int i=0; i<tmp_list.size(); i++){

        // Initially reset all of the values in bond_tors_vectors to -1 (bond direction does not
        // matter)
        bond_tors_vectors.clear();
        bond_tors_vectors.resize(tmp_list[i].mol.num_bonds, -1);

        // However, set bonds betweeen fragments to rotatable, adding degrees of freedom to the
        // minimizer. The value at rotatable bond index i in bond_tors_vectors should be the
        // atom_index of one atom involved in the bond, and it should be the atom that is *closer*
        // to the anchor fragment. amber_bt_id just says whether the bond is rotatable (positive
        // integer) or not (-1).
        for (int j=0; j<tmp_list[i].mol.num_bonds; j++){
            int res_origin = atoi(tmp_list[i].mol.atom_residue_numbers[tmp_list[i].mol.bonds_origin_atom[j]].c_str());
            int res_target = atoi(tmp_list[i].mol.atom_residue_numbers[tmp_list[i].mol.bonds_target_atom[j]].c_str());

            if (res_origin != res_target){

                if (res_origin < res_target){ bond_tors_vectors[j] = tmp_list[i].mol.bonds_origin_atom[j]; }
                else                        { bond_tors_vectors[j] = tmp_list[i].mol.bonds_target_atom[j]; }
    
                tmp_list[i].mol.amber_bt_id[j] = 1;
            } else {
               //BCF fix for dn refinement with rotatable bonds on user given anchors
                tmp_list[i].mol.amber_bt_id[j] = frag1.mol.amber_bt_id[j];
                //tmp_list[i].mol.amber_bt_id[j] = -1;
            }
        } 

        // Prepare for internal energy calculation for this fragment
        if (use_internal_energy){ prepare_internal_energy( tmp_list[i], score ); }
        // Copy the pre-min frag for the growth tree
        if (dn_write_growth_trees){
            tmp_list[i].frag_growth_tree.push_back(tmp_list[i].mol);
        }
        // Flexible minimization
        simplex.minimize_flexible_growth(tmp_list[i].mol, score, bond_tors_vectors);
        // Compute internal energy and primary score, store it in the dockmol
        if (score.use_primary_score) {
            valid_orient = score.compute_primary_score(tmp_list[i].mol);
        } else {
            valid_orient = true;
            tmp_list[i].mol.current_score = 0;
            tmp_list[i].mol.internal_energy = 0;
        }

        // If the energy or if the internal energy is greater than the cutoff, erase that conformer
        if ( tmp_list[i].mol.current_score > dn_pruning_conformer_score_cutoff  || 
             !valid_orient ){
            tmp_list[i].used = true;

// BCF part
            //tmp_list.erase(tmp_list.begin()+i);
            //i--;
        } else {
         
              if (use_internal_energy) {
                      if ( tmp_list[i].mol.internal_energy > ie_cutoff ) {
                          tmp_list[i].used = true;
                          //tmp_list.erase(tmp_list.begin()+i);
                          //i--;
                      }
              }
         }
// end BCF part
        bond_tors_vectors.clear();
    }

    // JDB 2022.11.16
    // Legacy pruning to keep only the best molecule if no dummy atoms are present
    bool dummy_found = true; // Defined as true so if legacy pruning is on, it'll use 6.9 behavior
    if (!dn_legacy_prune_complete){
        dummy_found = false;
        //Before RMSD pruning, check one mol to see if a dummy atom present
        for (int i=0; i<tmp_list[0].mol.num_atoms; i++) {
            if (tmp_list[0].mol.atom_types[i] == "Du"){
                dummy_found = true;
                break;
            }
        }
    }
    //If a dummy atom was found, this molecule is incomplete - cluster by
    //RMSD.
    if (dummy_found == true) {
        // Cluster post-minimization torsions by RMSD and remove those under a certain cutoff 
        float rmsd;
        // For every fragment in tmp_list
        for (int i=0; i<tmp_list.size(); i++){

            // If the fragment is not currently 'used'
            if (!tmp_list[i].used){

                // Then check the next fragment
                for (int j=i+1; j<tmp_list.size(); j++){

                    // If it is also not 'used'
                    if (!tmp_list[j].used){

                        // Then calculate the rmsd between the two
                        rmsd = calc_fragment_rmsd(tmp_list[i], tmp_list[j]);
                        rmsd = MAX(rmsd, 0.001);

                        // If the rmsd is below the specified cut-off, they are in the same 'cluster'
                        if ((float) j / rmsd > dn_pruning_clustering_cutoff) {
                            tmp_list[j].used = true;
                        }
                    }        
                }
            }
        }
    // TODO the above comparison is checking heavy atom RMSD for the WHOLE 
    // MOLECULE. In A&G, they only check heavy atom RMSD for the newly-added 
    // segment. We need to reconcile that
        int count_stuff = 0;
        // At this point, only cluster heads will not be flagged 'used', 
        //add those to list_of_frags
        if (dn_sample_torsions){
            std::vector<Fragment> refined_list_of_frags {};
            std::vector<float> tmp_degree_list_for_sort {};
            for (int i=0; i<tmp_list.size(); i++){
                if (!tmp_list[i].used){
                    refined_list_of_frags.push_back(tmp_list[i]);
                    tmp_degree_list_for_sort.push_back(tmp_rot_list[i]);
                    count_stuff++;
                }
            }
        
            if (tmp_degree_list_for_sort.size() !=
                refined_list_of_frags.size()){
                std::cout << "In DN_Build::sample_minimized_torsions_iso, " <<
                             "the degrees list and list_of_frags do not " <<
                             "have the same size" << std::endl;
                exit(0);
            }
        
            // If the lists have no candidate torsions, just return
            // 0.0 radians turned. 
            if (tmp_degree_list_for_sort.size() == 0 && 
                refined_list_of_frags.size() == 0) {
               return 0.0;
            }
      
            // Sort and access by index number
            std::vector<int> indices(refined_list_of_frags.size());
            std::iota(indices.begin(),indices.end(),0);
            std::sort(indices.begin(), indices.end(),
                      [&](int A, int B) -> bool {
                           return refined_list_of_frags[A].mol.current_score
                               < refined_list_of_frags[B].mol.current_score;
                       });

            // This is the final list you wanna push your best
            // scoring torsion turn
            list_of_frags.push_back(refined_list_of_frags[indices[0]]);

            // Clear the memory of tmp_list
            tmp_list.clear();
            refined_list_of_frags.clear();

            // Return the degrees turned that is associated
            // with the best scoring torsion turned
            return tmp_degree_list_for_sort[indices[0]];
        } else {

            frag_sort(tmp_list, fragment_sort);
            list_of_frags.push_back(tmp_list[0]);
            return 0.0;
        }

    //If a dummy atom was not found, this molecule is complete. Only keep the
    // best one by score.
    } else if (dummy_found == false) {
        //if (Parameter_Reader::verbosity_level() > 1) {
        if (verbose) {
            cout << "Molecule complete! Keeping best in vector by score." << endl;
        }

        if (dn_sample_torsions){

            if ((tmp_list.size() !=
                tmp_rot_list.size())){
                std::cout << "In DN_Build::sample_minimized_torsions_iso, " <<
                             "the degrees list and list_of_frags do not " <<
                             "have the same size" << std::endl;
                exit(0);
            }
      
            std::vector<int> indices(tmp_list.size());
            std::iota(indices.begin(),indices.end(),0);
            std::sort(indices.begin(), indices.end(),
                      [&](int A, int B) -> bool {
                           return tmp_list[A].mol.current_score
                               < tmp_list[B].mol.current_score;
                       });

            list_of_frags.push_back(tmp_list[indices[0]]);

            // Clear the memory of tmp_list
            tmp_list.clear();

            // Sample torsions for frag1, move torsions onto tmp_list
            return tmp_rot_list[indices[0]];
        } else {

            frag_sort(tmp_list, fragment_sort);
            list_of_frags.push_back(tmp_list[0]);
            return 0.0;
        }
    }

    return 0.0;

} // end DN_Build::sample_minimized_torsions_iso()

void
DN_Build::just_minimize(Fragment & frag1, vector <Fragment> & list_of_frags, 
                                     Master_Score & score, Simplex_Minimizer & simplex, 
                                     AMBER_TYPER & typer, bool prepare, bool minimize )
{
   

    // Activate all atoms and bonds prior to any scoring
    activate_mol(frag1.mol);

    // This will be true if energy is calculated correctly
    bool valid_orient = false;

    if ( prepare ) {
        // Declare a temporary fingerprint object and compute atom environments (necessary for Gasteiger)
        Fingerprint temp_finger;
        for (int i=0; i<frag1.mol.num_atoms; i++){
            frag1.mol.atom_envs[i] = temp_finger.return_environment(frag1.mol, i);
        }
                    
        // Prepare the molecule using the amber_typer to assign bond types to each bond
        frag1.mol.prepare_molecule();
        typer.prepare_molecule(frag1.mol, true, score.use_chem, score.use_ph4, score.use_volume);

        // Compute the charges, saving them on the mol object
        float total_charges = 0;
        total_charges = compute_gast_charges(frag1.mol);

        //BCF  04/28/16 update values, pruning is now done in sample_fraglib_graph
        //after gasteiger charges were computed
        calc_mol_wt(frag1.mol);
        calc_rot_bonds(frag1.mol);
        calc_formal_charge(frag1.mol);
    }


    // Initially reset all of the values in bond_tors_vectors to -1 (bond direction does not
    // matter)
    bond_tors_vectors.clear();
    bond_tors_vectors.resize(frag1.mol.num_bonds, -1);

    // However, set bonds betweeen fragments to rotatable, adding degrees of freedom to the
    // minimizer. The value at rotatable bond index i in bond_tors_vectors should be the
    // atom_index of one atom involved in the bond, and it should be the atom that is *closer*
    // to the anchor fragment. amber_bt_id just says whether the bond is rotatable (positive
    // integer) or not (-1).
    for (int j=0; j<frag1.mol.num_bonds; j++){
        int res_origin = atoi(frag1.mol.atom_residue_numbers[frag1.mol.bonds_origin_atom[j]].c_str());
        int res_target = atoi(frag1.mol.atom_residue_numbers[frag1.mol.bonds_target_atom[j]].c_str());

        if (res_origin != res_target){

            if (res_origin < res_target){ bond_tors_vectors[j] = frag1.mol.bonds_origin_atom[j]; }
            else                        { bond_tors_vectors[j] = frag1.mol.bonds_target_atom[j]; }
    
            frag1.mol.amber_bt_id[j] = 1;
        } else {
           //BCF fix for dn refinement with rotatable bonds on user given anchors
            frag1.mol.amber_bt_id[j] = frag1.mol.amber_bt_id[j];
            //frag1.mol.amber_bt_id[j] = -1;
        }
    } 
    

    // Prepare for internal energy calculation for this fragment
    if (use_internal_energy){ prepare_internal_energy( frag1, score ); }

    // Copy the pre-min frag for the growth tree
    if (dn_write_growth_trees){
        frag1.frag_growth_tree.push_back(frag1.mol);
    }
 
    if ( minimize ) { 
        // Flexible minimization
        simplex.minimize_flexible_growth(frag1.mol, score, bond_tors_vectors);
    }

    // Compute internal energy and primary score, store it in the dockmol
    if (score.use_primary_score) {
        valid_orient = score.compute_primary_score(frag1.mol);
    } else {
        valid_orient = true;
        frag1.mol.current_score = 0;
        frag1.mol.internal_energy = 0;
    }

    // If the energy or if the internal energy is greater than the cutoff, erase that conformer
    if ( frag1.mol.current_score > dn_pruning_conformer_score_cutoff  || 
         !valid_orient ){
        frag1.used = true;

// Bart
        //tmp_list.erase(tmp_list.begin()+i);
        //i--;
    } else {
     
          if (use_internal_energy) {
                  if ( frag1.mol.internal_energy > ie_cutoff ) {
                      frag1.used = true;
                  }
          }
     }
// eCF part
    bond_tors_vectors.clear();

    // At this point, only cluster heads will not be flagged 'used', add those to list_of_frags
    if (!frag1.used){
        list_of_frags.push_back(frag1);
    }

} // end DN_build::just_minimize

// +++++++++++++++++++++++++++++++++++++++++
// This function is used for a soft cutoff of MW and will accept molecules at a rate inversely proportional to how much they exceed the cutoff
bool DN_Build::mw_cutoff( Fragment &temp_molecule, int upper_lower_swtch ){ 
    // Define variables
    float  rand_num{};
    float  rand_num_dec{};
    double excessMW{};
    double Z_scoreExcess{};
    double acceptRate{};

    // if exceeds upper boundary accept or reject with some probability
    if (upper_lower_swtch == 1 && temp_molecule.mol.mol_wt > dn_upper_constraint_mol_wt){
        rand_num = (rand() % 100 +1) ; // This generates a random number from 1 to 100
        rand_num_dec = rand_num / 100;
        excessMW = temp_molecule.mol.mol_wt - dn_upper_constraint_mol_wt; //How much the growing molecule exceeds cutoff by
        Z_scoreExcess = excessMW / dn_MW_std_dev; // Similar to a Z score (how many std dev the excess mw is from the cutoff)
        acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess); //e raised to the expression inside the parentheses (similar to Metropolis)
        if (rand_num_dec < acceptRate) {
            return true;
        } else {
            return false;
        }
    }

    // if exceeds lower boundary accept or reject with some probability
    else if (upper_lower_swtch == 0 && temp_molecule.mol.mol_wt < dn_lower_constraint_mol_wt){
        rand_num = (rand() % 100 +1) ; // This generates a random number from 1 to 100
        rand_num_dec = rand_num / 100;
        excessMW = dn_lower_constraint_mol_wt - temp_molecule.mol.mol_wt; //How much the growing molecule exceeds cutoff by
        Z_scoreExcess = excessMW / dn_MW_std_dev; // Similar to a Z score (how many std dev the excess mw is from the cutoff)
        acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess); //e raised to the expression inside the parentheses (similar to Metropolis)
        if (rand_num_dec < acceptRate) {
            return true;
        } else {
            return false;
        }
    }

    // This last else should be entered if a soft cut off is used and the mol wt is within accepted MW bounds
    else {
        return true;
    }
} // End DN_Build::mw_cutoff()

// +++++++++++++++++++++++++++++++++++++++++
#ifdef BUILD_DOCK_WITH_RDKIT
bool DN_Build::clogp_cutoff( Fragment &temp_molecule, ostringstream &verb_str){//, Fragment &temp_molecule ){
    // Define variables
    float  rand_num{};
    float  rand_num_dec{};
    double excessCLOGP{};
    double Z_scoreExcess{};
    double acceptRate{};
    
    // directing clogp

    // soft cutoff
    if (temp_molecule.mol.clogp > dn_lower_clogp) {

        // Check if value is higher than upper boundary
        if (temp_molecule.mol.clogp > dn_upper_clogp) {

            // Do soft cutoff procedure
            rand_num = (rand() % 100 + 1);
            rand_num_dec = rand_num / 100;
            excessCLOGP = temp_molecule.mol.clogp - dn_upper_clogp;
            Z_scoreExcess = excessCLOGP / dn_clogp_std_dev;
            acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
            if (rand_num_dec < acceptRate) {
                if (dn_drive_verbose) {

                    cout<<setw(28)<<std::right <<" clogp_pseudoacc upperbound "<<verb_str.str();                                                    

                }
                return true; // molecule accepted and given to push_back
            } else {
                if (dn_drive_verbose) {

                    cout <<setw(28)<<std::right <<" clogp_REJECTED upperbound "<<verb_str.str();                                                    

                }
                return false;
            }
        } else { // mol.clogp > dn_lower_clogp && mol.clogp < dn_upper_clogp
            if (dn_drive_verbose) {

                cout << setw(28)<<std::right <<" clogp_acc "<<verb_str.str();
                                                                    
            } 
            return true;
        }
    } else { // mol.clogp < dn_lower_clogp. Needs soft cutoff

        // Do soft cutoff procedure
        rand_num = (rand() % 100 + 1);
        rand_num_dec = rand_num / 100;
        excessCLOGP = dn_lower_clogp - temp_molecule.mol.clogp;
        Z_scoreExcess = excessCLOGP / dn_clogp_std_dev;
        acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
        if (rand_num_dec < acceptRate) {
            if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" clogp_pseudoacc lowerbound "<<verb_str.str();                                                    
            }
            return true; // molecule accepted and given to push_back
        } else {
            if (dn_drive_verbose) {

                cout<<setw(28)<<std::right <<" clogp_REJECTED lowerbound " <<verb_str.str();
            }
            return false;
        }
    }
                                                           
} // End DN_Build::clogp_cutoff()


bool DN_Build::esol_cutoff( Fragment &temp_molecule, ostringstream & verb_str){//, Fragment &temp_molecule ){
    // Initialize variables
    float  rand_num{};
    float  rand_num_dec{};
    double excessESOL{};
    double Z_scoreExcess{};
    double acceptRate{};

    // soft cutoff
    if (temp_molecule.mol.esol > dn_lower_esol) {

        // Check if value is higher than upper boundary
        if (temp_molecule.mol.esol > dn_upper_esol) {

            // Do soft cutoff procedure
            rand_num = (rand() % 100 + 1);
            rand_num_dec = rand_num / 100;
            excessESOL = temp_molecule.mol.esol - dn_upper_esol;
            Z_scoreExcess = excessESOL / dn_esol_std_dev;
            acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
            if (rand_num_dec < acceptRate) {
                if (dn_drive_verbose) { 
                    cout <<setw(28)<<std::right <<" logs_pseudoacc upperbound "<< verb_str.str(); 

                }
                return true; // molecule accepted and given to push_back
            } else {
                if (dn_drive_verbose) { 

                   cout <<setw(28)<<std::right <<" logs_REJECTED upperbound "<<verb_str.str();                                                    

                } 
                return false;
            }
        } else { // mol.esol > dn_lower_esol && mol.esol < dn_upper_esol
            if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" logs_acc "<<verb_str.str();
            }
            return true;
        }
    } else { // mol.esol < dn_lower_esol. Needs soft cutoff

        // Do soft cutoff procedure
        rand_num = (rand() % 100 + 1);
        rand_num_dec = rand_num / 100;
        excessESOL = dn_lower_esol - temp_molecule.mol.esol;
        Z_scoreExcess = excessESOL / dn_esol_std_dev;
        acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
        if (rand_num_dec < acceptRate) {
            if (dn_drive_verbose) { 

                cout <<setw(28)<<std::right <<" logs_pseudoacc lowerbound "<< verb_str.str();
            }
            return true; // molecule accepted and given to push_back
        } else {
            if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" logs_REJECTED lowerbound "<< verb_str.str();
            }
            return false;                                           
        }
    }

} // End DN_Build::esol_cutoff()

bool DN_Build::tpsa_cutoff( Fragment &temp_molecule,ostringstream& verb_str ){//, Fragment &temp_molecule ){
    // Initialize variables
    float  rand_num{};
    float  rand_num_dec{};
    double excessTPSA{};
    double Z_scoreExcess{};
    double acceptRate{};

    // soft cutoff
    if (temp_molecule.mol.tpsa > dn_lower_tpsa) {

        // Check if value is higher than upper boundary
        if (temp_molecule.mol.tpsa > dn_upper_tpsa) {

            // Do soft cutoff procedure
            rand_num = (rand() % 100 + 1);
            rand_num_dec = rand_num / 100;
            excessTPSA = temp_molecule.mol.tpsa - dn_upper_tpsa;
            Z_scoreExcess = excessTPSA / dn_tpsa_std_dev;
            acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
            if (rand_num_dec < acceptRate) {
                if (dn_drive_verbose) {
                    cout <<setw(28)<<std::right <<" tpsa_pseudoacc upperbound "<<verb_str.str();

                }
                return true; // molecule accepted and given to push_back
            } else {
                if (dn_drive_verbose) {

                    cout <<setw(28)<<std::right <<" tpsa_REJECTED upperbound "<<verb_str.str();                                                    

                }
                return false;
            }
        } else { // mol.tpsa > dn_lower_tpsa && mol.tpsa < dn_upper_esol
            if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" tpsa_acc "<<verb_str.str();
                                                                    
                                                                    
            }
            return true;
        }
    } else { // mol.tpsa < dn_lower_tpsa. Needs soft cutoff

        // Do soft cutoff procedure
        rand_num = (rand() % 100 + 1);
        rand_num_dec = rand_num / 100;
        excessTPSA = dn_lower_tpsa - temp_molecule.mol.tpsa;
        Z_scoreExcess = excessTPSA / dn_tpsa_std_dev;
        acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
        if (rand_num_dec < acceptRate) {
            if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" tpsa_pseudoacc lowerbound "<< verb_str.str();
                                                                    
            }
            return true; // molecule accepted and given to push_back
        } else {
            if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" tpsa_REJECTED lowerbound "<<verb_str.str();
                                                                    
            }
            return false;
        }
    }

} // End DN_Build::tpsa_cutoff()

bool DN_Build::qed_cutoff( Fragment &temp_molecule,ostringstream& verb_str ){//, Fragment &temp_molecule ){
    // Initialize variables
    float  rand_num{};
    float  rand_num_dec{};
    double excessQED{};
    double Z_scoreExcess{};
    double acceptRate{};

    // directing qed
    if (temp_molecule.mol.qed_score < dn_lower_qed) {
        rand_num = (rand() % 100 + 1);
        rand_num_dec = rand_num / 100;
        excessQED = dn_lower_qed - temp_molecule.mol.qed_score;
        Z_scoreExcess = excessQED  / dn_qed_std_dev;
        acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
        if (rand_num_dec < acceptRate) {
            if (dn_drive_verbose) { 

                cout <<setw(28)<<std::right <<" qed_pseudoacc "<<verb_str.str();                                                    
                                                                    
            }
            return true; // molecule accepted and given to push_back
        } else {
            if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" qed_REJECTED "<<verb_str.str();
            }
            return false;
        }
    } else {
        if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" qed_acc "<<verb_str.str();
        }
        return true; // molecule accepted and given to push_back
    }
} // End DN_Build::qed_cutoff()


bool DN_Build::sa_cutoff( Fragment &temp_molecule, ostringstream & verb_str ){//, Fragment &temp_molecule ){
    // Initialize variables
    float  rand_num{};
    float  rand_num_dec{};
    double excessSA{};
    double Z_scoreExcess{};
    double acceptRate{};

    // directing sa
    if (temp_molecule.mol.sa_score > dn_upper_sa) {
        
        // proper method
        rand_num = (rand() % 100 + 1);
        rand_num_dec = rand_num / 100;
        excessSA = temp_molecule.mol.sa_score - dn_upper_sa;
        Z_scoreExcess = excessSA / dn_sa_std_dev;
        acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
        if (rand_num_dec < acceptRate) {
            if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" SynthA_pseudoacc " << verb_str.str();
                                                                    
            }
            return true;
        } else {
            if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" SynthA_REJECTED " << verb_str.str();
            }
            return false;
        }
    } else {
        if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" SynthA_acc " << verb_str.str();
        }
        return true;
    }
} // End DN_Build::sa_cutoff() 


bool DN_Build::stereocenter_cutoff( Fragment &temp_molecule,ostringstream & verb_str ){//, Fragment &temp_molecule ){
    // Initialize variables
    float  rand_num1{};
    float  rand_num2{};
    float  rand_num1_dec{};
    float  rand_num2_dec{};

    // smoothing #stereocenters
    if (temp_molecule.mol.num_stereocenters > dn_upper_stereocenter) {
        rand_num1 = (rand() % 100 + 1);
        rand_num1_dec = rand_num1 / 100;
        rand_num2 = (rand() % 100 + 1) / (temp_molecule.mol.num_stereocenters - 1);
        rand_num2_dec = rand_num2 / 100;
        if (rand_num1_dec < rand_num2_dec) {
            if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" stereocenter_pseudoacc "<<verb_str.str();
            }
            return true;
        } else {
            if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" stereocenter_REJECTED " <<verb_str.str();
            }
            return false;
        }
    } else {
        if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" stereocenter_acc " << verb_str.str();
        }
        return true;
    }
} // End DN_Build::stereocenter_cutoff()

bool DN_Build::pains_cutoff( Fragment &temp_molecule,ostringstream & verb_str ){//, Fragment &temp_molecule ){
    // Initialize variables
    float  rand_num1{};
    float  rand_num2{};
    float  rand_num1_dec{};
    float  rand_num2_dec{};

    // smoothing #pains
    if (temp_molecule.mol.pns > dn_upper_pains) {
        rand_num1 = (rand() % 100 + 1);
        rand_num1_dec = rand_num1 / 100; 
        rand_num2 = (rand() % 100 + 1) / (temp_molecule.mol.pns - 1);
        rand_num2_dec = rand_num2 / 100; 
        if (rand_num1_dec < rand_num2_dec) {
            if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" PAINS_pseudoacc " <<verb_str.str();
            }
            return true;
        } else {
            if (dn_drive_verbose) {

               cout <<setw(28)<<std::right <<" PAINS_REJECTED " <<verb_str.str();
            }
            return false;
        }
    } else {
        if (dn_drive_verbose) {

                cout <<setw(28)<<std::right <<" PAINS_acc " <<verb_str.str();
                
        }
        return true;
    }    
} // End DN_Build::pains_cutoff()



// bool dn_drive_clogp, bool dn_drive_esol, bool dn_drive_qed, bool dn_drive_sa, 
// and bool dn_drive_stereocenters are attributes of the DN_Build object
bool DN_Build::drive_growth( Fragment &temp_molecule ){
    //Initialize LOCAL variables    
    bool drive_clogp;
    bool drive_esol;
    bool drive_qed;
    bool drive_sa;
    bool drive_stereocenters;
    bool drive_pains;
    bool drive_tpsa;
    ostringstream verbose_comment; 
    if (dn_drive_verbose) {  
        verbose_comment<< "#stereo " << setw(3)<<fixed<< temp_molecule.mol.num_stereocenters
                       << " tpsa "  <<setw(10)<<fixed<< setprecision(6)<< temp_molecule.mol.tpsa
                       << " clogp "  <<setw(10)<<fixed<< setprecision(6)<< temp_molecule.mol.clogp
                       << " logs "    <<setw(10)<<fixed<< setprecision(6)<<temp_molecule.mol.esol
                       << " qed "    <<setw(10)<<fixed<< setprecision(6)<< temp_molecule.mol.qed_score
                       << " SynthA " <<setw(10)<<fixed<< setprecision(6)<< temp_molecule.mol.sa_score
                       << " #Pains " <<setw(3)<<fixed<< setprecision(6)<< temp_molecule.mol.pns <<endl;
    }
    // Does temp_molecule passes filters?
    //
    // Remember that only one false is enough to make this 
    // method return `false`. That's why the `else` clause
    // should be always true. Remember of TRUTH TABLES.
    //
    // STEREOCENTERS
    if (dn_drive_stereocenters) {
        drive_stereocenters = stereocenter_cutoff( temp_molecule,verbose_comment);
        temp_molecule.mol.fail_stereo = !drive_stereocenters;
    } else { drive_stereocenters = true; }
    // TPSA
    if (dn_drive_tpsa) {
        drive_tpsa = tpsa_cutoff( temp_molecule,verbose_comment);
        temp_molecule.mol.fail_tpsa = !drive_tpsa;
    } else { drive_tpsa = true; } 
    // CLOGP
    if (dn_drive_clogp) {
        drive_clogp = clogp_cutoff( temp_molecule,verbose_comment);
        temp_molecule.mol.fail_clogp = !drive_clogp;
    } else { drive_clogp = true; }
    // ESOL
    if (dn_drive_esol) {
        drive_esol = esol_cutoff( temp_molecule,verbose_comment);
        temp_molecule.mol.fail_esol = !drive_esol;
    } else { drive_esol = true; }
    // QED
    if (dn_drive_qed) {
        drive_qed = qed_cutoff( temp_molecule,verbose_comment);
        temp_molecule.mol.fail_qed = !drive_qed;
    } else { drive_qed = true; }
    // SA
    if (dn_drive_sa) {
        drive_sa = sa_cutoff( temp_molecule,verbose_comment);
        temp_molecule.mol.fail_sa = !drive_sa;
    } else { drive_sa = true; }
    // PAINS
    if (dn_drive_pains) {
        drive_pains = pains_cutoff( temp_molecule,verbose_comment);
        temp_molecule.mol.fail_pains = !drive_pains;
    } else { drive_pains = true; } 

    // Analyze truth value and return result
    if (drive_clogp && drive_esol && drive_tpsa && drive_qed &&
        drive_sa && drive_stereocenters && drive_pains) {
        // These are AND operations. If a single one of these booleans is
        // `false,` this method returns `false.`
        return true;
    } else { return false; }
} // End DN_Build::drive_growth()

#endif 

                                                              
// +++++++++++++++++++++++++++++++++++++++++
// Compute the RMSD between two fragments considering all atoms except H
// (equivalent to heavy atoms + dummy atoms)
float
DN_Build::calc_fragment_rmsd( Fragment & a, Fragment & b )
{
    // Declare rmsd and total number of heavy atoms
    float rmsd = 0.0;
    int heavy_total = 0;

    // Iterate through every atom in the fragment (pose a and pose b have the same atoms)
    for (int i=0; i<a.mol.num_atoms; i++){

        // If the atom is not hydrogen (it is important to include dummy positions in this)
        if (a.mol.atom_active_flags[i] && a.mol.atom_types[i].compare("H") != 0){

            // Then compute rmsd and add to running total
            rmsd +=  (((a.mol.x[i] - b.mol.x[i])*(a.mol.x[i] - b.mol.x[i])) +
                      ((a.mol.y[i] - b.mol.y[i])*(a.mol.y[i] - b.mol.y[i])) + 
                      ((a.mol.z[i] - b.mol.z[i])*(a.mol.z[i] - b.mol.z[i])));

            // Increment the count of atoms
            heavy_total++;
        }
    }

    // Make sure not to divide by 0, and compute final rmsd
    if (heavy_total > 0){
        rmsd = sqrt(rmsd / (float) heavy_total);
    } else {
        rmsd = 0.0;
    }

    return rmsd;

} // end DN_Build::calc_fragment_rmsd()



// +++++++++++++++++++++++++++++++++++++++++
// Populate a list of frags with different torsions. In A&G, this is called segment_torsion_drive
void
DN_Build::frag_torsion_drive( Fragment & new_frag, vector <Fragment> & frag_list )
{
    // Note: Every time two frags are combined into one, the last bond in the list is the new bond
    // between the two fragments. This is the bond about which torsions should be generated
    int num_torsions;
    float new_angle;    

    // Identify four atoms that make up dihedral. Atom1 / atom2 should be in original frag, atom3 /
    // atom4 should be in the newly-added frag. This way, the newly-added atoms will be rotated 
    // instead of the original anchor.
    int atom1 = -1;
    int atom2 =  new_frag.mol.bonds_origin_atom[ (new_frag.mol.num_bonds-1) ];
    int atom3 =  new_frag.mol.bonds_target_atom[ (new_frag.mol.num_bonds-1) ];
    int atom4 = -1;

    // Atom1 is any atom bound to atom2, except for atom3 - break as soon as one is identified
    for (int i=0; i<new_frag.mol.num_bonds; i++){
        if (new_frag.mol.bonds_origin_atom[i] == atom2){ 
            atom1 = new_frag.mol.bonds_target_atom[i];
            break;
        }
        if (new_frag.mol.bonds_target_atom[i] == atom2){ 
            atom1 = new_frag.mol.bonds_origin_atom[i];
            break;
        }
    }

    // Atom4 is any atom bound to atom3, except for atom2 - break as soon as one is identified
    for (int i=0; i<new_frag.mol.num_bonds; i++){
        if (new_frag.mol.bonds_origin_atom[i] == atom3){
            atom4 = new_frag.mol.bonds_target_atom[i];
            break;
        }
        if (new_frag.mol.bonds_target_atom[i] == atom3){
            atom4 = new_frag.mol.bonds_origin_atom[i];
            break;
        }
    }

    // Make sure all atoms have been assigned
    if (atom1 == -1 || atom4 == -1) {
        cout << "Warning: could not find all the atoms needed for torsions in "
             << "DN_Build::frag_torsion_drive" <<endl;
        return;
    }

    // Make sure there is no overlap in atom assignment - This might happen for adding new 
    // fragments with only two atoms (Du-X, for example). In that case, no torsions need to be
    // sampled.
    if (atom1 == atom3 || atom4 == atom2) {
        cout <<"Warning: atoms were assigned in a weird way in DN_Build::frag_torsion_drive" <<endl;
        cout <<"Perhaps a fragment was added that does not need to sample torsions." <<endl;
        return;
    }


    // Number of torsions comes from amber_bt_torsion_total (returns a different number for a given
    // bond type)
    num_torsions = new_frag.mol.amber_bt_torsion_total[ new_frag.mol.num_bonds-1 ];

    // Create torsions at appropriate intervals
    for (int i=0; i<num_torsions; i++){

        // Figure out the angle for this torsion and convert it to radians
        new_angle = new_frag.mol.amber_bt_torsions[ new_frag.mol.num_bonds-1 ][i];
        new_angle = (PI / 180) * new_angle;

        if (Parameter_Reader::verbosity_level() > 1){
            cout <<"FTD: amber_bt_torsions[" <<i <<"] = " 
                 <<new_frag.mol.amber_bt_torsions[new_frag.mol.num_bonds-1][i] <<endl;
            cout <<"FTD: new_angle = " <<new_angle <<endl;
        }

        // Push the fragment onto the end of the frag_list vector
        frag_list.push_back(new_frag);

        // Then adjust the torsions of that fragment according to the angle at this iteration
        frag_list[frag_list.size()-1].mol.set_torsion(atom1, atom2, atom3, atom4, new_angle);
    }

    return;

} // end DN_Build::frag_torsion_drive()

// +++++++++++++++++++++++++++++++++++++++++
void
DN_Build::just_turn_one_tors( Fragment & new_frag, float rad_turned)
{
    // Note: Every time two frags are combined into one, the last bond in the list is the new bond
    // between the two fragments. This is the bond about which torsions should be generated
    int num_torsions;
    float new_angle;    

    // Identify four atoms that make up dihedral. Atom1 / atom2 should be in original frag, atom3 /
    // atom4 should be in the newly-added frag. This way, the newly-added atoms will be rotated 
    // instead of the original anchor.
    int atom1 = -1;
    int atom2 =  new_frag.mol.bonds_origin_atom[ (new_frag.mol.num_bonds-1) ];
    int atom3 =  new_frag.mol.bonds_target_atom[ (new_frag.mol.num_bonds-1) ];
    int atom4 = -1;

    // Atom1 is any atom bound to atom2, except for atom3 - break as soon as one is identified
    for (int i=0; i<new_frag.mol.num_bonds; i++){
        if (new_frag.mol.bonds_origin_atom[i] == atom2){ 
            atom1 = new_frag.mol.bonds_target_atom[i];
            break;
        }
        if (new_frag.mol.bonds_target_atom[i] == atom2){ 
            atom1 = new_frag.mol.bonds_origin_atom[i];
            break;
        }
    }

    // Atom4 is any atom bound to atom3, except for atom2 - break as soon as one is identified
    for (int i=0; i<new_frag.mol.num_bonds; i++){
        if (new_frag.mol.bonds_origin_atom[i] == atom3){
            atom4 = new_frag.mol.bonds_target_atom[i];
            break;
        }
        if (new_frag.mol.bonds_target_atom[i] == atom3){
            atom4 = new_frag.mol.bonds_origin_atom[i];
            break;
        }
    }

    // Make sure all atoms have been assigned
    if (atom1 == -1 || atom4 == -1) {
        cout << "Warning: could not find all the atoms needed for torsions in "
             << "DN_Build::just_turn_one_tors" 
             << "Exiting..."<<endl;
        exit(1);
    }

    // Make sure there is no overlap in atom assignment - This might happen for adding new 
    // fragments with only two atoms (Du-X, for example). In that case, no torsions need to be
    // sampled.
    if (atom1 == atom3 || atom4 == atom2) {
        cout <<"Warning: atoms were assigned in a weird way in DN_Build::just_turn_one_tors" <<endl;
        cout <<"Perhaps a fragment was added that does not need to sample torsions."
             << "Exiting..."<<endl;
        exit(1);
    }

    if (Parameter_Reader::verbosity_level() > 1){
        cout <<"JFTD: new_angle = " <<rad_turned<<endl;
        cout <<"JFTD: For the fragment" <<new_frag.mol.title<<endl;
    }

    // Then adjust the torsions of that fragment according to the angle at this iteration
    new_frag.mol.set_torsion(atom1, atom2, atom3, atom4, rad_turned);

    return;

} // end DN_Build::just_turn_one_tors()


// +++++++++++++++++++++++++++++++++++++++++
// Populate a list of frags with different torsions. In A&G, this is called segment_torsion_drive
std::vector<float>
DN_Build::frag_torsion_drive_iso( Fragment & new_frag, vector <Fragment> & frag_list )
{
    // Note: Every time two frags are combined into one, the last bond in the list is the new bon
    // between the two fragments. This is the bond about which torsions should be generated
    int num_torsions;
    float new_angle;    
    std::vector<float> degree_list {};

    // Identify four atoms that make up dihedral. Atom1 / atom2 should be in original frag, atom3 /
    // atom4 should be in the newly-added frag. This way, the newly-added atoms will be rotated 
    // instead of the original anchor.
    int atom1 = -1;
    int atom2 =  new_frag.mol.bonds_origin_atom[ (new_frag.mol.num_bonds-1) ];
    int atom3 =  new_frag.mol.bonds_target_atom[ (new_frag.mol.num_bonds-1) ];
    int atom4 = -1;

    // Atom1 is any atom bound to atom2, except for atom3 - break as soon as one is identified
    for (int i=0; i<new_frag.mol.num_bonds; i++){
        if (new_frag.mol.bonds_origin_atom[i] == atom2){ 
            atom1 = new_frag.mol.bonds_target_atom[i];
            break;
        }
        if (new_frag.mol.bonds_target_atom[i] == atom2){ 
            atom1 = new_frag.mol.bonds_origin_atom[i];
            break;
        }
    }

    // Atom4 is any atom bound to atom3, except for atom2 - break as soon as one is identified
    for (int i=0; i<new_frag.mol.num_bonds; i++){
        if (new_frag.mol.bonds_origin_atom[i] == atom3){
            atom4 = new_frag.mol.bonds_target_atom[i];
            break;
        }
        if (new_frag.mol.bonds_target_atom[i] == atom3){
            atom4 = new_frag.mol.bonds_origin_atom[i];
            break;
        }
    }

    // Make sure all atoms have been assigned
    if (atom1 == -1 || atom4 == -1) {
        cout << "Warning: could not find all the atoms needed for torsions in "
             << "DN_Build::frag_torsion_drive_iso" <<endl;
        return {};
    }

    // Make sure there is no overlap in atom assignment - This might happen for adding new 
    // fragments with only two atoms (Du-X, for example). In that case, no torsions need to be
    // sampled.
    if (atom1 == atom3 || atom4 == atom2) {
        cout <<"Warning: atoms were assigned in a weird way in DN_Build::frag_torsion_drive_iso" <<endl;
        cout <<"Perhaps a fragment was added that does not need to sample torsions." <<endl;
        return {};
    }


    // Number of torsions comes from amber_bt_torsion_total (returns a different number for a given
    // bond type)
    num_torsions = new_frag.mol.amber_bt_torsion_total[ new_frag.mol.num_bonds-1 ];

    // Create torsions at appropriate intervals
    for (int i=0; i<num_torsions; i++){

        // Figure out the angle for this torsion and convert it to radians
        new_angle = new_frag.mol.amber_bt_torsions[ new_frag.mol.num_bonds-1 ][i];
        new_angle = (PI / 180) * new_angle;

        if (Parameter_Reader::verbosity_level() > 1){
            cout <<"FTD: amber_bt_torsions[" <<i <<"] = " 
                 <<new_frag.mol.amber_bt_torsions[new_frag.mol.num_bonds-1][i] <<endl;
            cout <<"FTD: new_angle = " <<new_angle <<endl;
        }

        // Push the fragment onto the end of the frag_list vector
        frag_list.push_back(new_frag);

        // Then adjust the torsions of that fragment according to the angle at this iteration
        frag_list[frag_list.size()-1].mol.set_torsion(atom1, atom2, atom3, atom4, new_angle);
        degree_list.push_back(new_angle);
    }

    return degree_list;

} // end DN_Build::frag_torsion_drive_iso()

// +++++++++++++++++++++++++++++++++++++++++
// Calculate the internal energy of a fragment
void
DN_Build::prepare_internal_energy( Fragment & frag1, Master_Score & score )
{
    // Pass the method ('2') to base_score
    score.primary_score->method = 2;

    // If using internal energy, initialize it here
    if (score.use_primary_score) {

        // Initialize internal energy
        score.primary_score->use_internal_energy = use_internal_energy;

        if (use_internal_energy) {
            // Pass vdw parameters to base_score
            score.primary_score->ie_att_exp = ie_att_exp;
            score.primary_score->ie_rep_exp = ie_rep_exp;
            score.primary_score->ie_diel = ie_diel;

            // Need to use DOCKMol with radii and segments assigned
            // it does not matter which atoms are labeled active
            score.primary_score->initialize_internal_energy(frag1.mol);
        }
    } 

    return;

} // end DN_Build::prepare_internal_energy()


// +++++++++++++++++++++++++++++++++++++++++
// Function to prune fragments by Hungarian RMSD using a generic best-first clustering
// algorithm
void
DN_Build::prune_h_rmsd_vec_pair( vector < pair <Fragment, int> > & list_of_frags )
{ 
    // Declare the Hungarian_RMSD object
    Hungarian_RMSD h;
    pair <double, int> result;
    // For every fragment in layer
    for (int i=0; i<list_of_frags.size(); i++){

        // If the fragment is not currently 'used'
        if (!list_of_frags[i].first.used){

            // Then check the next fragment
            for (int j=i+1; j<list_of_frags.size(); j++){

                // If it is also not 'used'
                if (!list_of_frags[j].first.used){

                    // Then calculate the RMSD between the two
                    result = h.calc_Hungarian_RMSD_dissimilar( list_of_frags[i].first.mol, 
                                                               list_of_frags[j].first.mol);

                    // If the RMSD is above the specified cut-off, they are in the same 'cluster'
                    if ( result.first < dn_heur_matched_rmsd && 
                         result.second < dn_heur_unmatched_num ) {
                        list_of_frags[j].first.used = true;
                    //    cout << result.first << "    " << result.second << endl; // JDB JDB JDB debug
                    }
                }
            }
        }
    }

    // Now that we know what the cluster heads are, clear the original vector, and copy all cluster
    // heads back onto original vector

    //iterate through the entire vector of fragments

    std::list<std::pair<Fragment, int> > tmp_linked_list(list_of_frags.begin(), list_of_frags.end());

    for(std::list< std::pair <Fragment, int> >::iterator pos = tmp_linked_list.begin(); 
                                                         pos != tmp_linked_list.end();) {
        //if they're listed as 'used' from above, remove them
        if ((*pos).first.used){
            pos = tmp_linked_list.erase(pos);
        //otherwise, move to next fragment
        } else {
            pos++;
        }
    }

    list_of_frags.clear();
    for (const auto& tmp : tmp_linked_list) {
        list_of_frags.push_back(tmp);
    }
    return;

} // end DN_Build::prune_h_rmsd_vec_pair()

// +++++++++++++++++++++++++++++++++++++++++
// Function to prune fragments by Hungarian RMSD using a generic best-first clustering
// algorithm, but splits the vector into MW chunks before doing so.
void
DN_Build::prune_h_rmsd_and_mw( vector <Fragment> & list_of_frags )
{ 
    //sort this list by mw, then by score 
    frag_sort(list_of_frags, mw_then_score_sort_fragments);
    vector<Fragment> pruned_collector;
    int mol = 0;
    float tolerance = 0.001;
    while (mol < list_of_frags.size()) {
        int pos1 = mol;
        //collect all molecules of the same MW - things are already score sorted internally
        vector<Fragment> same_mw_frags_collector;
        if (list_of_frags[pos1].mol.mol_wt == 0) {
            cout << "prune_h_rmsd: ~Assessing vector of fragments with MW 0~" << endl << endl;
        }

        //Check each set of values based on the current position - 
        //looking for all mols where compare_val <= mol.mol_wt <= compare_val + tolerance
        float mol_wt_compare_val = list_of_frags[pos1].mol.mol_wt;
        while (mol_wt_compare_val <= list_of_frags[mol].mol.mol_wt  && 
            list_of_frags[mol].mol.mol_wt <= mol_wt_compare_val + tolerance) {
            same_mw_frags_collector.push_back(list_of_frags[mol]);
            mol++;
            if (mol == list_of_frags.size()){
                break;
            }
        }
        //prune everything with the same MW by h_rmsd, and collect the remainder
        prune_h_rmsd(same_mw_frags_collector);

        if (same_mw_frags_collector.size() == 1){
            pruned_collector.push_back(same_mw_frags_collector[0]);
            same_mw_frags_collector.clear();
        } else { 

            for (unsigned int i =0;i<same_mw_frags_collector.size();i++){
                pruned_collector.push_back(same_mw_frags_collector[i]);
            }

            //vector< Fragment >::iterator pos = same_mw_frags_collector.begin();
            //while (pos != same_mw_frags_collector.end()) {
            //    pruned_collector.push_back((*pos));
            //    pos = same_mw_frags_collector.erase(pos);
            //}
        }
    }

    //clear the input vector and push the remaining molecules back onto it
    list_of_frags.clear();
    for (unsigned int i =0;i<pruned_collector.size();i++){
        list_of_frags.push_back(pruned_collector[i]);
    }

    //vector< Fragment >::iterator pos = pruned_collector.begin();
    //while (pos != pruned_collector.end()) {
    //    list_of_frags.push_back((*pos));
    //    pos = pruned_collector.erase(pos);
    //}

    //Sort the output so it's fully sorted by score
    frag_sort(list_of_frags,fragment_sort);
    return;

} // end DN_Build::prune_h_rmsd_and_mw()

// +++++++++++++++++++++++++++++++++++++++++
// Function to prune fragments by Hungarian RMSD using a generic best-first clustering
// algorithm, but splits the vector into MW chunks before doing so.
void
DN_Build::prune_h_rmsd_and_mw( vector <Fragment> & list_of_frags , 
                               vector <Fragment> & prune_rmsd_mw
                             )
{ 
    //sort this list by mw, then by score 
    frag_sort(list_of_frags, mw_then_score_sort_fragments);
    vector<Fragment> pruned_collector;
    int mol = 0;
    float tolerance = 0.001;
    while (mol < list_of_frags.size()) {
        int pos1 = mol;
        //collect all molecules of the same MW - things are already score sorted internally
        vector<Fragment> same_mw_frags_collector;
        if (list_of_frags[pos1].mol.mol_wt == 0) {
            cout << "prune_h_rmsd: ~Assessing vector of fragments with MW 0~" << endl << endl;
        }

        //Check each set of values based on the current position - 
        //looking for all mols where compare_val <= mol.mol_wt <= compare_val + tolerance
        float mol_wt_compare_val = list_of_frags[pos1].mol.mol_wt;
        while (mol_wt_compare_val <= list_of_frags[mol].mol.mol_wt  && 
            list_of_frags[mol].mol.mol_wt <= mol_wt_compare_val + tolerance) {
            same_mw_frags_collector.push_back(list_of_frags[mol]);
            mol++;
            if (mol == list_of_frags.size()){
                break;
            }
        }
        //prune everything with the same MW by h_rmsd, and collect the remainder
        prune_h_rmsd ( same_mw_frags_collector , prune_rmsd_mw );

        if (same_mw_frags_collector.size() == 1){
            pruned_collector.push_back(same_mw_frags_collector[0]);
            same_mw_frags_collector.clear();
        } else { 

            for (unsigned int i =0;i<same_mw_frags_collector.size();i++){
                pruned_collector.push_back(same_mw_frags_collector[i]);
            }

            //vector< Fragment >::iterator pos = same_mw_frags_collector.begin();
            //while (pos != same_mw_frags_collector.end()) {
            //    pruned_collector.push_back((*pos));
            //    pos = same_mw_frags_collector.erase(pos);
            //}
        }
    }

    //clear the input vector and push the remaining molecules back onto it
    list_of_frags.clear();
    for (unsigned int i =0;i<pruned_collector.size();i++){
        list_of_frags.push_back(pruned_collector[i]);
    }

    //vector< Fragment >::iterator pos = pruned_collector.begin();
    //while (pos != pruned_collector.end()) {
    //    list_of_frags.push_back((*pos));
    //    pos = pruned_collector.erase(pos);
    //}

    //Sort the output so it's fully sorted by score
    frag_sort(list_of_frags,fragment_sort);
    return;

} // end DN_Build::prune_h_rmsd_and_mw()

//Overloaded prunign function that
//writes out to mol2 on the fly
void
DN_Build::prune_h_rmsd_and_mw( vector <Fragment> & list_of_frags , 
                               vector <Fragment> & prune_rmsd_mw ,
                               int anchor_num,
                               int lay_num
                             )
{ 
    //sort this list by mw, then by score 
    frag_sort(list_of_frags, mw_then_score_sort_fragments);
    vector<Fragment> pruned_collector;
    int mol = 0;
    float tolerance = 0.001;
    while (mol < list_of_frags.size()) {
        int pos1 = mol;
        //collect all molecules of the same MW - things are already score sorted internally
        vector<Fragment> same_mw_frags_collector;
        if (list_of_frags[pos1].mol.mol_wt == 0) {
            cout << "prune_h_rmsd: ~Assessing vector of fragments with MW 0~" << endl << endl;
        }

        //Check each set of values based on the current position - 
        //looking for all mols where compare_val <= mol.mol_wt <= compare_val + tolerance
        float mol_wt_compare_val = list_of_frags[pos1].mol.mol_wt;
        while (mol_wt_compare_val <= list_of_frags[mol].mol.mol_wt  && 
            list_of_frags[mol].mol.mol_wt <= mol_wt_compare_val + tolerance) {
            same_mw_frags_collector.push_back(list_of_frags[mol]);
            mol++;
            if (mol == list_of_frags.size()){
                break;
            }
        }
        //prune everything with the same MW by h_rmsd, and collect the remainder
        prune_h_rmsd ( same_mw_frags_collector , 
                       prune_rmsd_mw , 
                       anchor_num , 
                       lay_num );

        if (same_mw_frags_collector.size() == 1){
            pruned_collector.push_back(same_mw_frags_collector[0]);
            same_mw_frags_collector.clear();
        } else { 

            for (unsigned int i =0;i<same_mw_frags_collector.size();i++){
                pruned_collector.push_back(same_mw_frags_collector[i]);
            }

            //vector< Fragment >::iterator pos = same_mw_frags_collector.begin();
            //while (pos != same_mw_frags_collector.end()) {
            //    pruned_collector.push_back((*pos));
            //    pos = same_mw_frags_collector.erase(pos);
            //}
        }
    }

    //clear the input vector and push the remaining molecules back onto it
    list_of_frags.clear();
    for (unsigned int i =0;i<pruned_collector.size();i++){
        list_of_frags.push_back(pruned_collector[i]);
    }

    //vector< Fragment >::iterator pos = pruned_collector.begin();
    //while (pos != pruned_collector.end()) {
    //    list_of_frags.push_back((*pos));
    //    pos = pruned_collector.erase(pos);
    //}

    //Sort the output so it's fully sorted by score
    frag_sort(list_of_frags,fragment_sort);
    return;

} // end DN_Build::prune_h_rmsd_and_mw()

// +++++++++++++++++++++++++++++++++++++++++
// Function to prune fragments by Hungarian RMSD using a generic best-first clustering
// algorithm
void
DN_Build::prune_h_rmsd( vector <Fragment> & list_of_frags )
{ 
    // Declare the Hungarian_RMSD object
    Hungarian_RMSD h;
    pair <double, int> result;

    // For every fragment in layer
    for (int i=0; i<list_of_frags.size(); i++){

        // If the fragment is not currently 'used'
        if (!list_of_frags[i].used){

            // Then check the next fragment
            for (int j=i+1; j<list_of_frags.size(); j++){

                // If it is also not 'used'
                if (!list_of_frags[j].used){

                    // Then calculate the RMSD between the two
                    result = h.calc_Hungarian_RMSD_dissimilar( list_of_frags[i].mol, 
                                                               list_of_frags[j].mol);

                    // If the RMSD is above the specified cut-off, they are in the same 'cluster'
                    if ( result.first < dn_heur_matched_rmsd && 
                         result.second < dn_heur_unmatched_num ) {
                        list_of_frags[j].used = true;
                        //cout << result.first << "    " << result.second << endl; // JDB JDB JDB debug
                    }
                }
            }
        }
    }

    // Now that we know what the cluster heads are, clear the original vector, and copy all cluster
    // heads back onto original vector
    vector <Fragment> temp_vec;
    temp_vec = list_of_frags;
    list_of_frags.clear();

    // For every fragment on temp_vec
    for (int i=0; i<temp_vec.size(); i++){
        // If it is not used, push it back onto list_of_frags
        if (!temp_vec[i].used){
            list_of_frags.push_back(temp_vec[i]);
        }
    }

    temp_vec.clear();

    return;

} // end DN_Build::prune_h_rmsd()

// +++++++++++++++++++++++++++++++++++++++++
// Function to prune fragments by Hungarian RMSD using a generic best-first clustering
// algorithm
void
DN_Build::prune_h_rmsd( vector <Fragment> & list_of_frags ,
                        vector <Fragment> & prune_rmsd_mw
                      )
{ 
    // Declare the Hungarian_RMSD object
    Hungarian_RMSD h;
    pair <double, int> result;

    // For every fragment in layer
    for (int i=0; i<list_of_frags.size(); i++){

        // If the fragment is not currently 'used'
        if (!list_of_frags[i].used){

            // Then check the next fragment
            for (int j=i+1; j<list_of_frags.size(); j++){

                // If it is also not 'used'
                if (!list_of_frags[j].used){

                    // Then calculate the RMSD between the two
                    result = h.calc_Hungarian_RMSD_dissimilar( list_of_frags[i].mol, 
                                                               list_of_frags[j].mol);

                    // If the RMSD is above the specified cut-off, they are in the same 'cluster'
                    if ( result.first < dn_heur_matched_rmsd && 
                         result.second < dn_heur_unmatched_num ) {
                        list_of_frags[j].used = true;
                        //cout << result.first << "    " << result.second << endl; // JDB JDB JDB debug
                    }
                }
            }
        }
    }

    // Now that we know what the cluster heads are, clear the original vector, and copy all cluster
    // heads back onto original vector
    vector <Fragment> temp_vec;
    temp_vec = list_of_frags;
    list_of_frags.clear();

    // For every fragment on temp_vec
    for (int i=0; i<temp_vec.size(); i++){
        // If it is not used, push it back onto list_of_frags
        if (!temp_vec[i].used){
            list_of_frags.push_back(temp_vec[i]);
        } else {
            prune_rmsd_mw.push_back(temp_vec[i]);
        }
    }

    temp_vec.clear();

    return;

} // end DN_Build::prune_h_rmsd()

// +++++++++++++++++++++++++++++++++++++++++
// Function to prune fragments by Hungarian RMSD using a generic best-first clustering
// algorithm
// Also this writes molecules on the fly
void
DN_Build::prune_h_rmsd( vector <Fragment> & list_of_frags ,
                        vector <Fragment> & prune_rmsd_mw ,
                        int anchor_num ,
                        int lay_num 
                      )
{ 
    // Declare the Hungarian_RMSD object
    Hungarian_RMSD h;
    pair <double, int> result;

    // For every fragment in layer
    for (int i=0; i<list_of_frags.size(); i++){

        // If the fragment is not currently 'used'
        if (!list_of_frags[i].used){

            // Then check the next fragment
            for (int j=i+1; j<list_of_frags.size(); j++){

                // If it is also not 'used'
                if (!list_of_frags[j].used){

                    // Then calculate the RMSD between the two
                    result = h.calc_Hungarian_RMSD_dissimilar( list_of_frags[i].mol, 
                                                               list_of_frags[j].mol);

                    // If the RMSD is above the specified cut-off, they are in the same 'cluster'
                    if ( result.first < dn_heur_matched_rmsd && 
                         result.second < dn_heur_unmatched_num ) {
                        list_of_frags[j].used = true;
                        //cout << result.first << "    " << result.second << endl; // JDB JDB JDB debug
                    }
                }
            }
        }
    }

    // Now that we know what the cluster heads are, clear the original vector, and copy all cluster
    // heads back onto original vector
    vector <Fragment> temp_vec;
    temp_vec = list_of_frags;
    list_of_frags.clear();

    // For every fragment on temp_vec
    for (int i=0; i<temp_vec.size(); i++){
        // If it is not used, push it back onto list_of_frags
        if (!temp_vec[i].used){
            list_of_frags.push_back(temp_vec[i]);
        } else {
            prune_rmsd_mw_counter_prune += 1;
            write_mol_to_file( dn_output_prefix,  
                               "prune_rmsd_mw", 
                               temp_vec[i], 
                               anchor_num,
                               lay_num
                             ); 
            //prune_rmsd_mw.push_back(temp_vec[i]);
        }
    }

    temp_vec.clear();

    return;

} // end DN_Build::prune_h_rmsd()


// +++++++++++++++++++++++++++++++++++++++++
// Calculate the molecular weight of the dockmol object of a fragment
// (modified from amber_typer.cpp)
void
DN_Build::calc_mol_wt( DOCKMol & mol )
{
    // atomic weights from General Chemistry 3rd Edition Darrell D. Ebbing
    float mw = 0.0;

    // For every atom in the dockmol object
    for (int i=0; i<mol.num_atoms; i++) {

        string atom = mol.atom_types[i];

        if ( atom == "H")
            { mw += ATOMIC_WEIGHT_H; }

        else if ( atom == "C.3" || atom == "C.2" || atom == "C.1" || atom == "C.ar" ||
                  atom  == "C.cat" )
            { mw += ATOMIC_WEIGHT_C; }

        else if ( atom == "N.4" || atom == "N.3" || atom == "N.2" || atom == "N.1" ||
                  atom  == "N.ar" || atom == "N.am" || atom == "N.pl3" )
            { mw += ATOMIC_WEIGHT_N; }

        else if ( atom == "O.3" || atom == "O.2" || atom == "O.co2" )
            { mw += ATOMIC_WEIGHT_O; }

        else if ( atom == "S.3" || atom == "S.2" || atom == "S.O" || atom == "S.o" ||
                  atom == "S.O2" || atom == "S.o2" )
            { mw += ATOMIC_WEIGHT_S; }

        else if ( atom == "P.3" )
            { mw += ATOMIC_WEIGHT_P; }

        else if ( atom == "F" )
            { mw += ATOMIC_WEIGHT_F; }

        else if ( atom == "Cl" )
            { mw += ATOMIC_WEIGHT_Cl; }

        else if ( atom == "Br" )
            { mw += ATOMIC_WEIGHT_Br; }

        else if ( atom == "I" )
            { mw += ATOMIC_WEIGHT_I; }

        else if ( atom == "Du" )
            { mw += 0.0; }

        else
            { cout <<"Warning: Did not recognize the atom_type " <<atom 
                   <<" in DN_Build::calc_mol_wt()" <<endl; }

    }

    // Assign the molecular weight to the mol object
    mol.mol_wt = mw;

    return;

} // end DN_Build::calc_mol_wt()



// +++++++++++++++++++++++++++++++++++++++++
// Given a fragment, populate the rot_bonds field
void
DN_Build::calc_rot_bonds( DOCKMol & mol )
{
    // The number of rotatable bonds
    int counter = 0;

    // Iterate over all bonds, check which ones are rotatable
    for (int i=0; i<mol.num_bonds; i++){
        if (mol.bond_is_rotor(i)){
            counter++;
        }
    }

    // Assign it directly to the referenced mol object
    mol.rot_bonds = counter;

    return;

} // end DN_Build::calc_rot_bonds();



// +++++++++++++++++++++++++++++++++++++++++
// Given a fragment, populate the formal_charge field
void
DN_Build::calc_formal_charge( DOCKMol & mol )
{
    float charge = 0.0;

    // Iterate over all atoms, find the partial charge
    for (int i=0; i<mol.num_atoms; i++){
        charge += mol.charges[i];
    }

    // Assign it directly to the referenced mol object
    mol.formal_charge = charge;

    return;

} // end DN_Build::calc_formal_charge();

// +++++++++++++++++++++++++++++++++++++++++
// LEP: Given a mol, calculate the hydrogen acceptors and donor
void
DN_Build::calc_num_HA_HD( DOCKMol & mol )
{
    // Populate HD fields
    int counter = 0;
    for (int i=0; i<mol.num_atoms; i++){
        if (mol.flag_acceptor[i] == true){
           counter++;
        }
    }
    mol.hb_acceptors = counter;

    // Populate HA fields
    counter = 0;
    for (int i=0; i<mol.num_atoms; i++){
        if (mol.flag_donator[i] == true){
           counter++;
        }
    }
    mol.hb_donors = counter;

   return;

} //end DN_Build::calc_num_HA_HD();


// +++++++++++++++++++++++++++++++++++++++++
// Activate all of the atoms AND bonds in a given DOCKMol
void
DN_Build::activate_mol( DOCKMol & mol )
{
     // Iterate through all atoms and set atom_active_flag to true
     for (int i=0; i<mol.num_atoms; i++){
               mol.atom_active_flags[i] = true;
     }

     // Iterate through all bonds and set bond_active_flag to true
     for (int i=0; i<mol.num_bonds; i++){
               mol.bond_active_flags[i] = true;
     }

    return;

} // end DN_Build::activate_mol()

// +++++++++++++++++++++++++++++++++++++++++
// deactivate all of the atoms AND bonds in a given DOCKMol
void
DN_Build::deactivate_mol( DOCKMol & mol )
{
     // Iterate through all atoms and set atom_active_flag to true
     for (int i=0; i<mol.num_atoms; i++){
               mol.atom_active_flags[i] = false;
     }

     // Iterate through all bonds and set bond_active_flag to true
     for (int i=0; i<mol.num_bonds; i++){
               mol.bond_active_flags[i] = false;
     }

    return;

} // end DN_Build::activate_mol()


// +++++++++++++++++++++++++++++++++++++++++
// Check to see if there are still unsatisfied attachment points in the fragment
bool
DN_Build::dummy_in_mol( DOCKMol & mol )
{
    // For every atom, if one of them is a 'Du', return true
    for (int i=0; i<mol.num_atoms; i++){
        if (mol.atom_types[i].compare("Du") == 0){
            return true; 
        }
    }

    return false;

} // end DN_Build::dummy_in_mol()



// +++++++++++++++++++++++++++++++++++++++++
// Given a specific dummy atom, change it to a hydrogen
void
DN_Build::dummy_to_H( Fragment & frag, int heavy, int dummy )
{
    // Change the atom type to H
    frag.mol.atom_types[dummy] = "H";

    // Re-index the names of all Hydrogen atoms
    // (this can be done later)

    // Calculate the desired bond length and remember as 'new_rad'
    float new_rad = calc_cov_radius(frag.mol.atom_types[heavy]) +
                    calc_cov_radius(frag.mol.atom_types[dummy]);

    // Calculate the x-y-z components of the bond vector (bond_vec)
    DOCKVector bond_vec;
    bond_vec.x = frag.mol.x[heavy] - frag.mol.x[dummy];
    bond_vec.y = frag.mol.y[heavy] - frag.mol.y[dummy];
    bond_vec.z = frag.mol.z[heavy] - frag.mol.z[dummy];

    // Normalize the bond vector then multiply each component by new_rad so that it is the desired 
    // length
    bond_vec = bond_vec.normalize_vector();
    bond_vec.x *= new_rad;
    bond_vec.y *= new_rad;
    bond_vec.z *= new_rad;

    // Change the coordinates of the Hydrogen atom so that the bond length is correct
    frag.mol.x[dummy] = frag.mol.x[heavy] - bond_vec.x;
    frag.mol.y[dummy] = frag.mol.y[heavy] - bond_vec.y;
    frag.mol.z[dummy] = frag.mol.z[heavy] - bond_vec.z;

    return;

} // end DN_Build::dummy_to_H()



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Accessory, debugging, or unused functions                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by energy score
bool 
fragment_sort( const Fragment & a, const Fragment & b )
{
    return (a.mol.current_score < b.mol.current_score);
}

// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by energy score
bool 
fragment_vec_pair_sort( const pair<Fragment, int> & a, const pair<Fragment, int> & b )
{
    return (a.first.mol.current_score < b.first.mol.current_score);
}


// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting FragGraph Tanimotos
bool 
fgpair_sort( const pair<float,int> & a, const pair<float,int> & b )
{
    return (a.first > b.first);
}

// +++++++++++++++++++++++++++++++++++++++++
// // LEP - Converts H's to Du's
void
DN_Build::convert_H_to_Du ( DOCKMol & mol, int a )
{
    Trace trace( "DN_Build::convert_H_to_Du( DOCKMol & mol, int a )" );
    // Replace with Du
    cout << "before atom type: " << mol.atom_types[a] << endl;
    mol.atom_types[a] = "Du";
    cout << "after atom type: " << mol.atom_types[a] << endl;
    //Add new charge
    mol.charges[a] = 0.0;
    //set new atom name
    stringstream dn_rf_rand;
    dn_rf_rand << "Du1";
    mol.atom_names[a] = dn_rf_rand.str();
    cout << "dummy name: "<< mol.atom_names[a] << endl;
    dn_rf_rand.clear();

    return;
}// end DN_Build::convert_H_to_Du()


// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by number of heavy atoms
// (more heavy atoms comes first), then consider number of att points
bool 
size_sort( const Fragment & frag1, const Fragment & frag2 )
{
    if (frag1.mol.heavy_atoms != frag2.mol.heavy_atoms)
        return (frag1.mol.heavy_atoms > frag2.mol.heavy_atoms);
    else
        return (frag1.aps.size() > frag2.aps.size());
}

// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors first by score, then by
// molecular weight JDB JDB JDB
bool 
mw_then_score_sort_fragments( const Fragment & frag1, const Fragment & frag2 )
{
    //First assessment: sort by molecular weight
    if (frag1.mol.mol_wt < frag2.mol.mol_wt) { return true; }
    if (frag2.mol.mol_wt < frag1.mol.mol_wt) { return false; }

    //If the MW cases are equal, then sort by score
    if (frag1.mol.current_score < frag2.mol.current_score) { return true; }
    if (frag2.mol.current_score < frag1.mol.current_score) { return false; }

    //a catch for if both mw and score are equivalent
    return false;
}

// +++++++++++++++++++++++++++++++++++++++++
// Debugging function to print the TorEnv data structure
void
DN_Build::print_torenv( vector <TorEnv> tmp_torenv_vector )
{
    int total = 0;
    cout <<endl;

    for (int i=0; i<tmp_torenv_vector.size(); i++){
        cout <<tmp_torenv_vector[i].origin_env <<":" <<endl;

        for (int j=0; j<tmp_torenv_vector[i].target_envs.size(); j++){

            cout <<"\t\t" <<tmp_torenv_vector[i].target_envs[j] <<"\t" 
                 <<tmp_torenv_vector[i].target_freqs[j] <<endl;
            //cout <<tmp_torenv_vector[i].origin_env <<"-" 
            //     <<tmp_torenv_vector[i].target_envs[j] <<endl;

        }

        cout <<"\t\tTotal: " <<tmp_torenv_vector[i].target_envs.size() <<endl;
        total = total+tmp_torenv_vector[i].target_envs.size();
    }

    cout <<"Final Total = " <<total <<endl;

    return;

} // end DN_Build::print_torenv()

void DN_Build::bickel_write_out(std::vector <Fragment> write_vector, std::string filename, int out_counter, bool overwrite) {
    ostringstream fout_name;
    fout_name << filename << out_counter << ".mol2";
    fstream fout_file;

    if (overwrite) { 
        fout_file.open(fout_name.str().c_str(), fstream::out);
    } else {
        fout_file.open(fout_name.str().c_str(), fstream::out|fstream::app);
    }

    int write_to_file_counter=0;
    for (int i=0; i<write_vector.size(); i++){
        //LEP - 2018.06.20 
        ostringstream new_comment;
        new_comment << denovo_name<< "_layer"<<out_counter+1<<"_"<< write_to_file_counter;
        write_vector[i].mol.title = new_comment.str();
        new_comment.clear();
        // write mols to file
        fout_file << write_vector[i].mol.current_data;
        //LEP added in unique denovo name and scientific notation to floating point for output
        fout_file << DELIMITER << setw(STRING_WIDTH) << "Name:" << setw(FLOAT_WIDTH) << write_vector[i].mol.title <<endl;
        fout_file <<DELIMITER<<setw(STRING_WIDTH)<<"Molecular_Weight:"<<setw(FLOAT_WIDTH)<<write_vector[i].mol.mol_wt <<endl;
        fout_file <<DELIMITER<<setw(STRING_WIDTH)<<"DOCK_Rotatable_Bonds:"<<setw(FLOAT_WIDTH)<<write_vector[i].mol.rot_bonds <<endl;
        fout_file << fixed <<DELIMITER<<setw(STRING_WIDTH)<<"Formal_Charge:"<<setw(FLOAT_WIDTH)<<write_vector[i].mol.formal_charge <<endl;
        fout_file << fixed <<DELIMITER<<setw(STRING_WIDTH)<<"HBond_Acceptors:"<<setw(FLOAT_WIDTH)<<write_vector[i].mol.hb_acceptors <<endl;
        fout_file << fixed <<DELIMITER<<setw(STRING_WIDTH)<<"HBond_Donors:"<<setw(FLOAT_WIDTH)<<write_vector[i].mol.hb_donors <<endl;
        fout_file <<DELIMITER<<setw(STRING_WIDTH)<<"Layer_Completed:"<<setw(FLOAT_WIDTH)<<out_counter+1 <<endl;
        fout_file <<DELIMITER<<setw(STRING_WIDTH)<<"Frag_String:"<<setw(FLOAT_WIDTH)<<write_vector[i].mol.energy << "\n"<<endl;
        Write_Mol2(write_vector[i].mol, fout_file);
        write_to_file_counter++;
    }
    fout_file.close();
    return;
}

// +++++++++++++++++++++++++++++++++++++++++
// Debugging function to print the FragGraph data structure
void
DN_Build::print_fraggraph( )
{

    // First print some data about what is connected to stdout:

    // Scaf_link_sid
    cout <<"scaf_link_sid_graph summary:" <<endl;
    for (int i=0; i<scaf_link_sid_graph.size(); i++){
        cout <<"i="  <<i <<"  j="; 
        for (int j=0; j<scaf_link_sid_graph[i].rankvec.size(); j++){ 
            cout <<scaf_link_sid_graph[i].rankvec[j] << ",";
        } cout <<endl;
    } cout <<endl;

    // Then write multi-mol2 files, one for each fragment, containing the fragment itself as the
    // parent, then the top N similar things along with Tanimoto. N right now is pulled from the
    // input file as dn_graph_breadth

    // Scaf_link_sid
    for (int i=0; i<scaf_link_sid_graph.size(); i++){
        ostringstream fout_name;
        if      (i<10)  { fout_name <<"scaf_link_sid_00" <<i <<".mol2"; }
        else if (i<100) { fout_name <<"scaf_link_sid_0" <<i <<".mol2"; }
        else if (i<1000){ fout_name <<"scaf_link_sid_" <<i <<".mol2"; }
        else { cout <<"this is probably too many to write to file" <<endl; exit(0); }
        fstream fout;
        fout.open(fout_name.str().c_str(), fstream::out|fstream::app);
        fout_name.clear();

        // First write the parent
        fout <<DELIMITER<<setw(STRING_WIDTH)<<"Tanimoto_to_parent:"<<setw(FLOAT_WIDTH)<<"1.0" <<endl;
        Write_Mol2(scaf_link_sid[i].mol, fout);

        // Then write the top N children
        for (int j=0; j<dn_graph_breadth; j++){
            // Remember: the tanvec and rankvec are offset by one
            fout <<DELIMITER<<setw(STRING_WIDTH)<<"Tanimoto_to_parent:"<<setw(FLOAT_WIDTH)<<scaf_link_sid_graph[i].tanvec[j+1].first <<endl;
            Write_Mol2(scaf_link_sid[ scaf_link_sid_graph[i].rankvec[j] ].mol, fout);
        }
        fout.close();
    }

    // We just wrote a ton of stuff to disk, so exit:
    exit(0);
} // end DN_Build::print_fraggraph()

//JDB - calculate pairwise distances between all atoms in a molecule
void DN_Build::calc_pairwise_distance(DOCKMol & a)
{

    //Declare RMSD and total number of heavy atoms
    float distance = 0.0;
    int heavy_total = 0;
    bool short_bond_found = false;
    //iterate through every atom
    for (int i=0; i<a.num_atoms; i++){
        for (int j=i+1; j<a.num_atoms; j++){
            distance = 0.0;
            //calculates the distance between two atoms           
            distance += sqrt((a.x[i]-a.x[j]) * (a.x[i]-a.x[j]) +
                        (a.y[i]-a.y[j]) * (a.y[i]-a.y[j]) +
                        (a.z[i]-a.z[j]) * (a.z[i]-a.z[j]));
            //if the distance is below a hardcoded cutoff, output an error
            if (distance < 0.1){

                //cout << a.x[i] << " " << a.y[i] << " " << a.z[i] << endl;
                short_bond_found = true;
                //cout << "SMALL BOND SMALL BOND SMALL BOND" << endl;
                cout << "SMALL BOND SMALL BOND: " << distance << endl;

            }
        }
    }


    if (short_bond_found){
        ostringstream bad_mols;
        bad_mols << "bad_mols.mol2";
        fstream fout_file;
        fout_file.open(bad_mols.str().c_str(), fstream::out|fstream::app);
        //fout_file << a.title << endl;
        Write_Mol2(a, fout_file);
        fout_file.close();
    }
}

void DN_Build::print_out_coordinates(DOCKMol & a, string append_val)
{
    for(int atom=0; atom < a.num_atoms; atom++){
        cout << a.atom_names[atom] << " |  " <<
        a.x[atom] << " | " << a.y[atom] << " | " << a.z[atom] << " | "<< append_val << endl;
    }
    cout << endl;

}
void DN_Build::frag_sort(std::vector<Fragment> & vec_frag, std::function<bool(const Fragment&, const Fragment&)> func){

    std::list<Fragment> list_growing(vec_frag.begin(), vec_frag.end());
    list_growing.sort(func);

    vec_frag.clear();
    for (Fragment tmp_frag : list_growing){
        vec_frag.push_back(tmp_frag); 
    }

}

void DN_Build::frag_sort(std::vector< std::pair <Fragment, int> > & vec_frag, 
                         std::function<bool(const std::pair <Fragment, int> &, const std::pair <Fragment, int> &)> func){

    std::list< pair < Fragment, int >> list_growing(vec_frag.begin(), vec_frag.end());
    list_growing.sort(func);
    for (pair<Fragment,int> tmp_frag : list_growing){
        vec_frag.push_back(tmp_frag); 
    }

}
// end DN_Build::frag_sort

/*
// +++++++++++++++++++++++++++++++++++++++++
// Given a dockmol, return true if it is within the desired molecular weight range, else false
bool
DN_Build::prune_molecular_weight( DOCKMol & mol )
{
    // This needs to be tested. Actually it probably needs to be rewritten.
    return (mol.mol_wt >= dn_target_mol_wt_lower && mol.mol_wt <= dn_target_mol_wt_upper);

} // end DN_Build::prune_molecular_weight();
*/

//By: Brock Boysan - 2022.07.26
//Sorts the list by molecular weight, then ensures there is only one occurance of each molecular weight
//This function will overwrite the first vector that is passed to it with the uniqufied vector, it will also return that vector
//This function will store pruned molecules in the second vector passed to it
void DN_Build::make_unique(vector <Fragment> &list_of_frags, vector <Fragment> &dn_prune_dump_list)
{
    //hehehehe, now I don't have to fu- mess around with all the other print outs
    int local_verbosity = 0;
    if (Parameter_Reader::verbosity_level() > 1 ){
        local_verbosity = 1;
    }
    vector <Fragment> final_list;
    vector <vector <Fragment>> tmp_list;
    int i = list_of_frags.size()-1;
    float tanimoto;
    int initial_size = list_of_frags.size();
    Fingerprint finger;
    //the leeway on MW matching
    float tolerance = 0.001; 
    //MW sort
    frag_sort(list_of_frags, mw_then_score_sort_fragments);
    if (local_verbosity > 0){
        for (int j=0; j<list_of_frags.size(); j++){
            cout << endl << "Index " << j << " INITIAL MW: " << list_of_frags[j].mol.mol_wt << endl;
        }
    }

    while(i > 0){
        if (local_verbosity > 0 ){
            cout << endl << "Processing index " << i << " MW: " << list_of_frags[i].mol.mol_wt << endl;
        }
        //case 1: no mw match & no duplicates stored
        if ((list_of_frags[i].mol.mol_wt < list_of_frags[i-1].mol.mol_wt - tolerance || list_of_frags[i].mol.mol_wt > list_of_frags[i-1].mol.mol_wt + tolerance) && tmp_list.size() == 0){
            //push the mol to the final list (vector)
            if (local_verbosity > 0){
                cout << endl << "Index " << i << " is a non duplicate frag" << endl;
            }
            final_list.push_back(list_of_frags[i]);
            list_of_frags.pop_back();
            i--;
            continue;
        }

        //case 2: mw weights match 
        if (list_of_frags[i].mol.mol_wt >= list_of_frags[i-1].mol.mol_wt - tolerance && list_of_frags[i].mol.mol_wt <= list_of_frags[i-1].mol.mol_wt + tolerance){
            if (local_verbosity > 0){
                cout << endl << "Index " << i << " possible duplicate frag w/ MW: " << list_of_frags[i].mol.mol_wt << endl;
            }
            //cout << endl << "Index " << i << "Processed frags tanimoto: " << tanimoto << endl;

            //if this list of mw ensembles is empty, initialize it
            if (tmp_list.size() == 0){
                vector <Fragment> tmp_vec;
                //store the mol in its own list as the head node
                tmp_list.push_back(tmp_vec);
                tmp_list[0].push_back(list_of_frags[i]);
                list_of_frags.pop_back();
            }
            //check to see if this is a duplicate against already stored mols, otherwise make it a head node of its own ensemble
            else{
                bool matched = false;
                for (int j=0; j<tmp_list.size(); j++){
                    tanimoto = finger.compute_tanimoto(list_of_frags[i].mol, tmp_list[j][0].mol);
                    if (tanimoto == 1) {
                        if (local_verbosity > 0) {
                            cout << endl << "Duplicate Frag Detected" << endl;
                        }
                        //store the mol in the list with its other duplicate(s)
                        tmp_list[j].push_back(list_of_frags[i]);
                        list_of_frags.pop_back();
                        matched = true;
                        break;
                    }
                }
                //otherwise make it a head node of its own ensemble
                if (!matched){
                    if (local_verbosity > 0){
                                cout << endl << "New head node found! It's tanimoto didn't match any of the previous but its MW did!" << endl;
                    }
                    vector <Fragment> tmp_vec;
                    //store the mol in itws own list as the head node
                    tmp_list.push_back(tmp_vec);
                    tmp_list[tmp_list.size()-1].push_back(list_of_frags[i]);
                    list_of_frags.pop_back();
                }
                
            }
            i--;
            continue;
        }
         
        //case 3: no mw match & duplicates stored
        if ((list_of_frags[i].mol.mol_wt < list_of_frags[i-1].mol.mol_wt - tolerance || list_of_frags[i].mol.mol_wt > list_of_frags[i-1].mol.mol_wt + tolerance) && tmp_list.size() > 0){
            if (local_verbosity > 0){
                cout << endl << i << " index frag last of matching mw: " << list_of_frags[i].mol.mol_wt << endl;
            }
            bool matched = false;
            //loop through and see if our last mol is a duplicate
            for (int j=0; j<tmp_list.size(); j++){
                //cout << endl << " loop at " << j << endl;
                tanimoto = finger.compute_tanimoto(list_of_frags[i].mol, tmp_list[j][0].mol);
                if (tanimoto == 1){
                    if (local_verbosity > 0){
                         cout << endl << "matching tanimoto, appending last molecule of matching MW to ensemble of duplicates " << endl;
                    }
                    tmp_list[j].push_back(list_of_frags[i]);
                    list_of_frags.pop_back();
                    matched = true;
                    break;
                }
            }
            //if our mol is not a duplicate, push it to the final list of mols
            if (!matched){
                final_list.push_back(list_of_frags[i]);
                list_of_frags.pop_back();
            }
            //go through the stored duplicates, for each ensemble only save the top scoring one
            for (int j=0; j<tmp_list.size(); j++){
                frag_sort(tmp_list[j], fragment_sort);
                if (local_verbosity > 0){
                     cout << endl << " sorting bin " << j << " and pushing to output" << endl;
                }
                for (int k=0; k<dn_max_duplicates_per_molecule + 1; k++){
                    final_list.push_back(tmp_list[j][tmp_list[j].size()-1]);
                    tmp_list[j].pop_back();
                }

                dn_prune_dump_list.insert(dn_prune_dump_list.end(), tmp_list[j].begin(), tmp_list[j].end() );

                tmp_list[j].clear();
            }

            tmp_list.clear();
            i--;
            continue;
        }
        
    }

    if (local_verbosity > 0){
        cout << endl << " processing last mol" << endl;
    }
    //end while loop
    //deal with the last molecule in the original list
    if (tmp_list.size() > 0 && i == 0){
        bool matched = false;
        for (int j=0; j<tmp_list.size(); j++){
            tanimoto = finger.compute_tanimoto(list_of_frags[i].mol, tmp_list[j][0].mol);
            if (tanimoto == 1){
                tmp_list[j].push_back(list_of_frags[i]);
                list_of_frags.pop_back();
                matched = true;
                break;
            }
        }
        if (!matched){
            final_list.push_back(list_of_frags[i]);
            list_of_frags.pop_back();
        }
        for (int j=0; j<tmp_list.size(); j++){
            frag_sort(tmp_list[j], fragment_sort);
            for (int k=0; k<dn_max_duplicates_per_molecule + 1; k++){
                final_list.push_back(tmp_list[j][tmp_list[j].size()-1]);
                tmp_list[j].pop_back();
            }

            dn_prune_dump_list.insert(dn_prune_dump_list.end(), tmp_list[j].begin(), tmp_list[j].end() );

            tmp_list[j].clear();
        }
        tmp_list.clear();

    }
    else if(i==0){
        if (local_verbosity > 0){
             cout << endl << " last mol was not a mw match" << endl;
        }
        final_list.push_back(list_of_frags[0]);
        list_of_frags.pop_back();
    }
      
    if (local_verbosity > 0){
        cout << endl << " finished main body of make_unique function " << endl;
    }
 
    //purge the original list for debugging purposes 
    list_of_frags.clear();
    frag_sort(final_list, fragment_sort);
    if (local_verbosity > 0 ){
         cout <<"Starting vector length after clearing: " << list_of_frags.size() <<endl;
         cout <<"Starting vector length: " << initial_size <<endl;
         cout <<"Uniqueified vector length: " << final_list.size() <<endl;
         cout <<"dn duplicates length: " << dn_prune_dump_list.size() << endl;
         ////print out the survivors' molecular weights
         //for (int i = 0; i < final_list.size(); i++){
         //    cout <<endl << " ### " << i << " MW: " << final_list[i].mol.mol_wt <<endl;
         //}
    }

    for (Fragment tmp_frag : final_list){
        list_of_frags.push_back(tmp_frag);
    }
    
    //list_of_frags = final_list; 
    final_list.clear();
    return;

}  //end make_unique()

//################ WRITE OUT MOLS FUNCTION

//write ONE molecule to file to save on memory
void DN_Build::write_mol_to_file ( string dn_output_prefix, 
                                   string frag_mol_name,  
                                   Fragment & frag_mol, 
                                   int anchor_num, int counter ) {

  if ( !dn_write_prune_dump || !dn_write_out_duplicates) { return; } 

  // Filename for filehandle
  ostringstream fout_frag_mol_name;
  fout_frag_mol_name <<dn_output_prefix <<".anchor_" <<(anchor_num+1) <<"." << frag_mol_name << "_layer_" <<counter+1 <<".mol2";
  fstream fout_frag_mol;

  bool overwrite = false;

  // Open the filehandle
  if (overwrite) {
      fout_frag_mol.open(fout_frag_mol_name.str().c_str(), fstream::out);
  } else {
      fout_frag_mol.open(fout_frag_mol_name.str().c_str(), fstream::out|fstream::app);
  }

  int write_to_file_counter = 0;
  int failed_cap_H  = 0;
  int failed_vtm  = 0;

  bool frag_mol_flag = true;
  // Iterate over each attachment point
  for (int j=0; j<frag_mol.aps.size(); j++){

      string bond_type;

      // Loop over all bonds in frag1
      for (int k=0; k<frag_mol.mol.num_bonds; k++){
          if (frag_mol.mol.bonds_origin_atom[k] == frag_mol.aps[j].dummy_atom ||
              frag_mol.mol.bonds_target_atom[k] == frag_mol.aps[j].dummy_atom   ){
              if (frag_mol.mol.bonds_origin_atom[k] == frag_mol.aps[j].heavy_atom ||
                  frag_mol.mol.bonds_target_atom[k] == frag_mol.aps[j].heavy_atom   ){

                  // This is the bond you're looking for...
                  bond_type = frag_mol.mol.bond_types[k];
                  break;
              }
          }
      }

      // Only cap it if the dummy is on a single bond 
      if ( bond_type == "1" ){

          dummy_to_H(frag_mol, frag_mol.aps[j].heavy_atom, frag_mol.aps[j].dummy_atom);

          frag_mol.mol.mol_wt += ATOMIC_WEIGHT_H;
          continue;

      // Otherwise, ignore this molecule
      } 

      frag_mol_flag = false;
      failed_cap_H +=1;
      break;
  }

  // LEP - run frag mol through multi torenv check to not make garbage
  if (valid_torenv_multi(frag_mol)){

      // If there was no problem capping the molecules, then write to file
      if (frag_mol_flag){
          //PAK-reactivate all atoms so it can be recognized by  Write_Mol2() for mol2 writing
          for (int q =0; q<frag_mol.mol.num_atoms;q++){
              frag_mol.mol.atom_active_flags[q] = true;
          }

          calc_mol_wt(frag_mol.mol);
          calc_rot_bonds(frag_mol.mol);
          calc_formal_charge(frag_mol.mol);

          //TODO there is an extra space being printed after current_data
          fout_frag_mol << frag_mol.mol.current_data <<endl;
          fout_frag_mol << DELIMITER << setw(STRING_WIDTH)
                        << "Molecular_Weight:" << setw(FLOAT_WIDTH)<<frag_mol.mol.mol_wt <<endl;
          fout_frag_mol << DELIMITER<<setw(STRING_WIDTH)
                        << "DOCK_Rotatable_Bonds:" << setw(FLOAT_WIDTH)<<frag_mol.mol.rot_bonds <<endl;
          fout_frag_mol << DELIMITER<<setw(STRING_WIDTH) 
                        << "Formal_Charge:"<<setw(FLOAT_WIDTH)<<frag_mol.mol.formal_charge <<endl;
          fout_frag_mol << DELIMITER<<setw(STRING_WIDTH)
                        << "Layer_Completed:"<<setw(FLOAT_WIDTH)<<counter+1 <<endl;
          fout_frag_mol << DELIMITER<<setw(STRING_WIDTH)
                        << "Frag_String:"<<setw(FLOAT_WIDTH)<<frag_mol.mol.energy <<endl<<endl;
          Write_Mol2( frag_mol.mol, fout_frag_mol );
          write_to_file_counter++;
      }
  } else { failed_vtm +=1; }

  if (frag_mol_name ==  "prune_rmsd_mw") {
      failed_cap_H_prune_rmsd_mw += failed_cap_H;
      failed_vtm_prune_rmsd_mw   += failed_vtm;
  }
  if (frag_mol_name ==  "prune_root") {
      failed_cap_H_prune_root += failed_cap_H;
      failed_vtm_prune_root   += failed_vtm;
  }
  if (frag_mol_name == "cand_root_ign") {
      failed_cap_H_cand_root_ign += failed_cap_H;
      failed_vtm_cand_root_ign   += failed_vtm;
  }
  if (frag_mol_name ==  "filtered_comp") {
      failed_cap_H_filtered_comp += failed_cap_H;
      failed_vtm_filtered_comp   += failed_vtm;
  }
  if (frag_mol_name ==  "duplicate_dump"){
    failed_cap_H_duplicate_dump += failed_cap_H;
    failed_vtm_duplicate_dump   += failed_vtm;
  }
  
}

void DN_Build::write_mols_to_file(string dn_output_prefix, 
                                  string frag_vec_name, 
                                  vector<Fragment>  &frag_vec, 
                                  int anchor_num, int counter){
  // Filename for filehandle
  ostringstream fout_frag_vec_name;
  fout_frag_vec_name <<dn_output_prefix <<".anchor_" <<(anchor_num+1) <<"." << frag_vec_name << "_layer_" <<counter+1 <<".mol2";
  fstream fout_frag_vec;

  bool overwrite = false;

  // Open the filehandle
  if (overwrite) {
        fout_frag_vec.open(fout_frag_vec_name.str().c_str(), fstream::out);
    } else {
        fout_frag_vec.open(fout_frag_vec_name.str().c_str(), fstream::out|fstream::app);
    }

  int write_to_file_counter = 0;
  int failed_cap_H  = 0;
  int failed_vtm  = 0;
  // Check the molecules and cap as appropriate
  for (int i=0; i<frag_vec.size(); i++){

      bool frag_vec_flag = true;

      // Iterate over each attachment point
      for (int j=0; j<frag_vec[i].aps.size(); j++){

          string bond_type;

          // Loop over all bonds in frag1
          for (int k=0; k<frag_vec[i].mol.num_bonds; k++){
              if (frag_vec[i].mol.bonds_origin_atom[k] == frag_vec[i].aps[j].dummy_atom ||
                  frag_vec[i].mol.bonds_target_atom[k] == frag_vec[i].aps[j].dummy_atom   ){
                  if (frag_vec[i].mol.bonds_origin_atom[k] == frag_vec[i].aps[j].heavy_atom ||
                      frag_vec[i].mol.bonds_target_atom[k] == frag_vec[i].aps[j].heavy_atom   ){

                      // This is the bond you're looking for...
                      bond_type = frag_vec[i].mol.bond_types[k];
                      break;
                  }
              }
          }

          // Only cap it if the dummy is on a single bond 
          if ( bond_type == "1" ){

              dummy_to_H(frag_vec[i], frag_vec[i].aps[j].heavy_atom, frag_vec[i].aps[j].dummy_atom);

              frag_vec[i].mol.mol_wt += ATOMIC_WEIGHT_H;
              continue;
          // Otherwise, ignore this molecule
          } 

          frag_vec_flag = false;
          failed_cap_H +=1;
          break;
      }

      // LEP - run frag_vec mols through multi torenv check to not make garbage
      if (valid_torenv_multi(frag_vec[i])){

          // If there was no problem capping the molecules, then write to file
          if (frag_vec_flag){
              //PAK-reactivate all atoms so it can be recognized by  Write_Mol2() for mol2 writing
              for (int q =0; q<frag_vec[i].mol.num_atoms;q++){
                  frag_vec[i].mol.atom_active_flags[q] = true;
              }

              calc_mol_wt(frag_vec[i].mol);
              calc_rot_bonds(frag_vec[i].mol);
              calc_formal_charge(frag_vec[i].mol);

              //TODO there is an extra space being printed after current_data
              fout_frag_vec << frag_vec[i].mol.current_data <<endl;
              fout_frag_vec <<DELIMITER<<setw(STRING_WIDTH)<<"Molecular_Weight:"<<setw(FLOAT_WIDTH)<<frag_vec[i].mol.mol_wt <<endl;
              fout_frag_vec <<DELIMITER<<setw(STRING_WIDTH)<<"DOCK_Rotatable_Bonds:"<<setw(FLOAT_WIDTH)<<frag_vec[i].mol.rot_bonds <<endl;
              fout_frag_vec <<DELIMITER<<setw(STRING_WIDTH)<<"Formal_Charge:"<<setw(FLOAT_WIDTH)<<frag_vec[i].mol.formal_charge <<endl;
              fout_frag_vec <<DELIMITER<<setw(STRING_WIDTH)<<"Layer_Completed:"<<setw(FLOAT_WIDTH)<<counter+1 <<endl;
              fout_frag_vec <<DELIMITER<<setw(STRING_WIDTH)<<"Frag_String:"<<setw(FLOAT_WIDTH)<<frag_vec[i].mol.energy <<endl<<endl;
              Write_Mol2( frag_vec[i].mol, fout_frag_vec );
              write_to_file_counter++;
          }
      }
  }

  failed_vtm = frag_vec.size() - write_to_file_counter - failed_cap_H;

  if (failed_vtm > 0 || failed_cap_H > 0){ 
      cout <<"          For " << setw(15) << std::right << frag_vec_name << " frags, " << setw(4) << std::right<< failed_cap_H << " frags couldn't be written out due to failed H-capping protocol\n" 
           <<"                                     "<< setw(4) << std::right << failed_vtm << " frags couldn't be written out due to failed torsion checks." << std::endl;

  }
  // Clear some memory
  fout_frag_vec_name.clear();
  fout_frag_vec.close();


}

std::pair<std::vector<float>,std::vector<Fragment> >
DN_Build::attach_and_get_torsions( Fragment layer_frag, int att_1, Fragment att_frag, int att_2,
                                   Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER &  typer){

    std::vector<Fragment>  return_vec {};
    std::vector<Fragment>  tmp_list {};
    std::vector<float> tmp_rot_list {};
 
    //frag1 is the fragment we are attaching TO, and frag2
    //is the fragment we are trying to attach.
    int dummy1 = layer_frag.aps[att_1].dummy_atom;
    int dummy2 = att_frag.aps[att_2].dummy_atom;
    int heavy1 = layer_frag.aps[att_1].heavy_atom;
    int heavy2 = att_frag.aps[att_2].heavy_atom;

    Fragment att_frag_copy = att_frag; 
    Fragment layer_frag_copy = layer_frag; 
    Fragment frag_combined = attach(layer_frag_copy, dummy1, heavy1, att_frag_copy, dummy2, heavy2);


    ////// Activate all atoms and bonds prior to any scoring
    //activate_mol(frag_combined.mol);

    //// Declare a temporary fingerprint object and compute atom environments (necessary for Gasteiger)
    //Fingerprint temp_finger;
    //for (int i=0; i<frag_combined.mol.num_atoms; i++){
    //    frag_combined.mol.atom_envs[i] = temp_finger.return_environment(frag_combined.mol, i);
    //}
    //            
    ////// Prepare the molecule using the amber_typer to assign bond types to each bond
    frag_combined.mol.prepare_molecule();
    //typer.prepare_molecule(frag_combined.mol, true, score.use_chem, score.use_ph4, score.use_volume);
    typer.prepare_for_torsions(frag_combined.mol);

    ////// Compute the charges, saving them on the mol object
    //float total_charges = 0;
    //total_charges = compute_gast_charges(frag_combined.mol);

    ////BCF  04/28/16 update values, pruning is now done in sample_fraglib_graph
    ////after gasteiger charges were computed
    //calc_mol_wt(frag_combined.mol);
    //calc_rot_bonds(frag_combined.mol);
    //calc_formal_charge(frag_combined.mol);

    // Sample torsions for frag_combined,
    // move torsions onto tmp_list
    if (dn_sample_torsions){
        // return list of degrees
        tmp_rot_list = frag_torsion_drive_iso(frag_combined, tmp_list);
    } else {

        // if not just fill tmp_list with a frags
        // if you don't want to sample torsions 
        return std::make_pair(tmp_rot_list, tmp_list);
    }    

    for (float angle: tmp_rot_list){
        Fragment att_frag_copy_tor = att_frag; 
 
        // Figure out what translation is required to move the dummy atom to the origin
        // for the iso test frags
        DOCKVector trans1;
        trans1.x = -att_frag_copy_tor.mol.x[dummy2];
        trans1.y = -att_frag_copy_tor.mol.y[dummy2];
        trans1.z = -att_frag_copy_tor.mol.z[dummy2];

        att_frag_copy_tor.mol.translate_mol(trans1);

        rotate_on_x_axis(att_frag_copy_tor.mol,dummy2,heavy2,angle);

        DOCKVector trans2;
        trans2.x = layer_frag.mol.x[heavy1];
        trans2.y = layer_frag.mol.y[heavy1];
        trans2.z = layer_frag.mol.z[heavy1];

        att_frag_copy_tor.mol.translate_mol(trans2);
    
        return_vec.push_back(att_frag_copy_tor);
    }
    

    
    return std::make_pair(tmp_rot_list, return_vec);
}


std::pair<float,Fragment>
DN_Build::get_best_iso_head_tors( Fragment layer_frag, int att_1, std::pair<std::vector<float>,std::vector<Fragment> > torsion_cand , int att_2 ,
                                  Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER &  typer){

    //left best isohead , right combined frag
    std::pair<float,Fragment> best_orient_pair;
    best_orient_pair.first  = -9999;


    std::vector <Fragment> tmp_list {};
    std::vector <float >   angle_list = {};
    std::vector <Fragment> original_att_frag_list {};

    // Create a temporary vector to populate with different torsions
    std::vector<float> tmp_rot_list;

    // Compute the the finger print once then
    // then save that information for the next coming for loop
    Fragment tors_cand_first     = torsion_cand.second[0];
    Fragment frag_combined_first = attach_isosteres(layer_frag ,att_1, tors_cand_first, att_2);


    Fingerprint temp_finger_first;
    for (int i=0; i<frag_combined_first.mol.num_atoms; i++){
        frag_combined_first.mol.atom_envs[i] = temp_finger_first.return_environment(frag_combined_first.mol, i);
    }

    frag_combined_first.mol.prepare_molecule();
    typer.prepare_molecule(frag_combined_first.mol, true, score.use_chem, score.use_ph4, score.use_volume);
    float first_total_charges = compute_gast_charges(frag_combined_first.mol);

    calc_mol_wt(frag_combined_first.mol);
    calc_rot_bonds(frag_combined_first.mol);
    calc_formal_charge(frag_combined_first.mol);

    //for (Fragment tors_cand: torsion_cand){
    for (int tors_ind = 0; tors_ind< torsion_cand.second.size(); tors_ind++){
        Fragment tors_cand = torsion_cand.second[tors_ind];
        Fragment frag_combined = attach_isosteres(layer_frag ,att_1, tors_cand, att_2);

        //activate_mol(frag_combined.mol);

        // Declare a temporary fingerprint object and compute atom environments (necessary for Gasteiger)
        for (int i=0; i<frag_combined.mol.num_atoms; i++){
            frag_combined.mol.atom_envs[i] = frag_combined_first.mol.atom_envs[i];
        }
                    
        //// Prepare the molecule using the amber_typer to assign bond types to each bond
        frag_combined.mol.prepare_molecule();
        typer.prepare_molecule(frag_combined.mol, true, score.use_chem, score.use_ph4, score.use_volume);

        float total_charges = first_total_charges;
   
        for (int j=0;j<frag_combined.mol.num_atoms; j++){
            frag_combined.mol.charges[j] = frag_combined_first.mol.charges[j];
        }

        frag_combined.mol.mol_wt = frag_combined_first.mol.mol_wt;
        frag_combined.mol.rot_bonds = frag_combined_first.mol.rot_bonds;
        frag_combined.mol.formal_charge = frag_combined_first.mol.formal_charge;
        
        original_att_frag_list.push_back( frag_combined );
        tmp_list.push_back( frag_combined );

    }

 
    // For all of the torsions that were just generated for the most recent fragment addition,
    for (int i=0; i<tmp_list.size(); i++){
 
        bool valid_orient = false;

        // Initially reset all of the values in bond_tors_vectors to -1 (bond direction does not
        // matter)
        bond_tors_vectors.clear();
        bond_tors_vectors.resize(tmp_list[i].mol.num_bonds, -1);

        // However, set bonds betweeen fragments to rotatable, adding degrees of freedom to the
        // minimizer. The value at rotatable bond index i in bond_tors_vectors should be the
        // atom_index of one atom involved in the bond, and it should be the atom that is *closer*
        // to the anchor fragment. amber_bt_id just says whether the bond is rotatable (positive
        // integer) or not (-1).
        for (int j=0; j<tmp_list[i].mol.num_bonds; j++){
            int res_origin = atoi(tmp_list[i].mol.atom_residue_numbers[tmp_list[i].mol.bonds_origin_atom[j]].c_str());
            int res_target = atoi(tmp_list[i].mol.atom_residue_numbers[tmp_list[i].mol.bonds_target_atom[j]].c_str());

            if (res_origin != res_target){

                if (res_origin < res_target){ bond_tors_vectors[j] = tmp_list[i].mol.bonds_origin_atom[j]; }
                else                        { bond_tors_vectors[j] = tmp_list[i].mol.bonds_target_atom[j]; }

                tmp_list[i].mol.amber_bt_id[j] = 1;
            } else {
               //BCF fix for dn refinement with rotatable bonds on user given anchors
               tmp_list[i].mol.amber_bt_id[j] = original_att_frag_list[i].mol.amber_bt_id[j];
               //tmp_list[i].mol.amber_bt_id[j] = -1;
            }
        }

        // Prepare for internal energy calculation for this fragment
        if (use_internal_energy){ prepare_internal_energy( tmp_list[i], score ); }

        // Flexible minimization
        simplex.minimize_flexible_growth(tmp_list[i].mol, score, bond_tors_vectors);
        // Compute internal energy and primary score, store it in the dockmol
        if (score.use_primary_score) {
            valid_orient = score.compute_primary_score(tmp_list[i].mol);
        } else {
            valid_orient = true;
            tmp_list[i].mol.current_score = 0;
            tmp_list[i].mol.internal_energy = 0;
        }

        // If the energy or if the internal energy is greater than the cutoff, erase that conformer
        if ( tmp_list[i].mol.current_score > dn_pruning_conformer_score_cutoff  ||
             !valid_orient ){
            tmp_list[i].used = true;

        } else {

              if (use_internal_energy) {
                      if ( tmp_list[i].mol.internal_energy > ie_cutoff ) {
                          tmp_list[i].used = true;
                          //tmp_list.erase(tmp_list.begin()+i);
                          //i--;
                      }
              }
         }
// end BCF part
        bond_tors_vectors.clear();

    }
    
    // TODO the above comparison is checking heavy atom RMSD for the WHOLE 
    // MOLECULE. In A&G, they only check heavy atom RMSD for the newly-added 
    // segment. We need to reconcile that
    int count_stuff = 0;
    // At this point, only cluster heads will not be flagged 'used', 
    //add those to list_of_frags

    if (verbose) {
        cout << "Molecule complete! Keeping best in vector by score." << endl;
    }

    if ((tmp_list.size() !=
        torsion_cand.second.size())){
        std::cout << "in dn_build::sample_minimized_torsions_iso, " <<
                     "the tmp_list and torsion_cand do not " <<
                     "have the same size" << std::endl;
        exit(0);
    }

    std::vector<int> indices(tmp_list.size());
    std::iota(indices.begin(),indices.end(),0);
    std::sort(indices.begin(), indices.end(),
              [&](int a, int b) -> bool {
                   return tmp_list[a].mol.current_score
                       < tmp_list[b].mol.current_score;
               });
    for (int index: indices){
        if ( !tmp_list[index].used ){ 
            best_orient_pair.first  = torsion_cand.first[indices[index]];
            //best_orient_pair.second = torsion_cand.second[indices[index]];
            best_orient_pair.second = tmp_list[index];
            break; 
        }

    }

    // sample torsions for frag1, move torsions onto tmp_list
    return best_orient_pair;

}

// +++++++++++++++++++++++++++++++++++++++++
// Given two fragments and connection point data, combine them into one and return it
Fragment
DN_Build::attach_isosteres( Fragment& frag1, int attachment_point_1,
                             Fragment  frag2, int attachment_point_2)
{


    std::pair<std::vector<std::vector<double>>,Fragment> return_pair;
    return_pair.first = {{0,0,0},
                         {0,0,0},
                         {0,0,0}};

    //assign some variables from the input fragments
    //frag1 is the fragment we are attaching TO, and frag2
    //is the fragment we are trying to attach.
    int dummy1 = frag1.aps[attachment_point_1].dummy_atom;
    int dummy2 = frag2.aps[attachment_point_2].dummy_atom;
    int heavy1 = frag1.aps[attachment_point_1].heavy_atom;
    int heavy2 = frag2.aps[attachment_point_2].heavy_atom;


    // Step 6. Connect the two fragments and return one object
    Fragment frag_combined = attach(frag1, dummy1, heavy1, frag2, dummy2, heavy2);

    // Also, copy the growth tree over if those are turned on
    if (dn_write_growth_trees){
        for (int i=0; i<frag1.frag_growth_tree.size(); i++){
            frag_combined.frag_growth_tree.push_back(frag1.frag_growth_tree[i]);
        }
    }

    return frag_combined;

} // end DN_Build::combine_fragments_iso()
