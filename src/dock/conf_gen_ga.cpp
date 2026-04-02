#include "conf_gen_ga.h"
#include "conf_gen_ag.h"
//#include "conf_gen_dn.h"
#include "conf_gen_dn_ga.h"
#include "fingerprint.h"
#include "hungarian.h"
#include "gasteiger.h"
#include "trace.h"
//#include "master_score.h"


#ifdef BUILD_DOCK_WITH_RDKIT
#include "rdtyper.h"
class RDTYPER;
#endif

using namespace std;

class Bump_Filter;
class Fingerprint;
class AMBER_TYPER;
class Master_Score;

// +++++++++++++++++++++++++++++++++++++++++
// static member initializers

// These are the same as those in Base_Score.
const string GA_Recomb::DELIMITER    = "########## ";
const int    GA_Recomb::FLOAT_WIDTH  = 20;
const int    GA_Recomb::STRING_WIDTH = 17 + 19;



// +++++++++++++++++++++++++++++++++++++++++
// Some constructors and destructors

GA_Recomb::GA_Recomb(){
    used = false;
    ga_selection_method_elitism = false;
    ga_selection_method_tournament = false;
    ga_selection_method_roulette = false;
    ga_selection_method_sus = false;
    ga_selection_method_metropolis = false;
    ga_niching = false;
    ga_niche_sharing = false;
    ga_selection_extinction = false;
    ga_extinction_on = false;
    ga_extinction_switch_selection = false;
    ga_secondary_selection_method_elitism = false;
    ga_secondary_selection_method_tournament = false;
    ga_secondary_selection_method_roulette = false;
    ga_secondary_selection_method_sus = false;
    ga_secondary_selection_method_metropolis = false;
    num_segs_removed = 0;
}
GA_Recomb::~GA_Recomb(){
    parents.clear();
    children.clear();
    pruned_children.clear();
    scored_generation.clear();
    tmp_parents.clear();
    mutants.clear();
    mutated_parents.clear();
}


Tor_Env::Tor_Env(){
    origin_env = "";
}
Tor_Env::~Tor_Env(){
    target_envs.clear();
    target_freqs.clear();
}



// +++++++++++++++++++++++++++++++++++++++++
// Read parameters from the dock.in file
void
GA_Recomb::input_parameters( Parameter_Reader & parm )
{
    cout << "\nGenetic Algorithm Recombination Parameters\n";
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    // BLOCK 1: INPUT FILES
    // Parent ensemble input file
    ga_molecule_file = parm.query_param( "ga_molecule_file",  "ga_molecule_file.mol2" );

    // LEP enter ga utilities
    ga_utilities = (parm.query_param( "ga_utilities", "no", "yes no" ) == "yes");
    if (ga_utilities){
        ga_calc_parent_pairwise_tan = (parm.query_param( "ga_calc_parent_pairwise_tan", "no", "yes no" ) == "yes");
        if (ga_calc_parent_pairwise_tan){
            ga_tan_similarity_prune = (parm.query_param( "ga_tan_similarity_prune", "no", "yes no" ) == "yes");
            ga_tan_best_first_clustering = (parm.query_param( "ga_tan_best_first_clustering", "no", "yes no" ) == "yes");
        }
        ga_calc_parent_pairwise_hms = (parm.query_param( "ga_calc_parent_pairwise_hms", "no", "yes no" ) ==     "yes");
        ga_charge_parent_gasteiger = (parm.query_param( "ga_charge_parent_gasteiger", "no", "yes no" ) == "yes");
    }else{
        ga_calc_parent_pairwise_tan = false;
        ga_calc_parent_pairwise_hms = false;
        ga_charge_parent_gasteiger = false;
    }

    // De Novo fragment library files
    ga_fraglib_scaffold_file  = parm.query_param("ga_fraglib_scaffold_file",  "fraglib_scaffold.mol2");
    ga_fraglib_linker_file    = parm.query_param("ga_fraglib_linker_file",    "fraglib_linker.mol2");
    ga_fraglib_sidechain_file = parm.query_param("ga_fraglib_sidechain_file", "fraglib_sidechain.mol2");

    #ifdef BUILD_DOCK_WITH_RDKIT
    ga_sa_fraglib_path = parm.query_param("sa_fraglib_path","sa_fraglib.dat");
    if (ga_fragMap.empty() == true){
        std::ifstream fin(ga_sa_fraglib_path);
        double key; 
        double val; 
        while (fin >> key >> val) {
            ga_fragMap[key] = val; 
        }
        fin.close();
    }
    ga_PAINS_path = parm.query_param("PAINS_path","pains_table.dat");
    std::vector<std::string> PAINStmp;
    if (PAINStmp.empty() == true){
        std::ifstream fin(ga_PAINS_path);
        std::string tmp_string;
        while (fin >> tmp_string) {
            PAINStmp.push_back(tmp_string);
        }
        fin.close();
    }
    int numofvectors = PAINStmp.size() - 1; 
    for (int i{0}; i < numofvectors; ) {
        ga_PAINSMap[PAINStmp[i]] = " " + PAINStmp[i + 1];
        i += 2; 
    }
    PAINStmp.clear();     
    #endif 
      
    // Use allowable bond environment table. Use is recommended
    ga_use_torenv_table = true;
    ga_torenv_table = parm.query_param( "ga_torenv_table", "fraglib_torenv.dat" );
    //}
    
    // BLOCK 2: TERMINATION CHECK
    // Max number of generations to be completed
    max_generations = atof(parm.query_param( "ga_max_generations", "100" ).c_str());


    // BLOCK 3: CROSSOVER 
    // Pick GA sampling method
    // Turn crossover on/off CS - LEP
    ga_xover_on = (parm.query_param( "ga_xover_on", "yes", "yes no" ) == "yes");
    
    // only show xover options when xover turned on
    ga_check_only = false; //preset as false JDB
    if (ga_xover_on){
        ga_xover_sampling_method_rand = 
		(parm.query_param( "ga_xover_sampling_method_rand", "yes", "yes no" ) == "yes");

        // Max number of offspring generated 
        ga_xover_max = atof(parm.query_param( "ga_xover_max", "150" ).c_str());

        // Sets allowable distance between atoms on overlapping parent bonds  
        ga_bond_tolerance = atof(parm.query_param( "ga_bond_tolerance", "0.5" ).c_str());

        // Set bond angle 
        ga_angle_cutoff = atof(parm.query_param( "ga_angle_cutoff", "0.14" ).c_str());

        // Check for molecules with overlap
        ga_check_overlap = (parm.query_param( "ga_check_overlap", "no", "yes no" ) == "yes");
        // LEP - leave hardcoded for developers only
        if ( ga_check_overlap ){
            ga_check_only = false;
            // Reset max gen to 1
            if ( ga_check_only ){ max_generations = 1; }
        }
    } else {
        ga_xover_sampling_method_rand = true;
        ga_xover_max = 150;
        ga_bond_tolerance = 0.5;
        ga_angle_cutoff = 0.14;
        ga_check_overlap = false;
        ga_check_only = false;
    }

    // BLOCK 4: MUTATIONS
    //ga_mutations = (parm.query_param( "ga_mutations", "yes", "yes no" ) == "yes");
    ga_mutations = true;
    if ( ga_mutations ){
        ga_mutate_addition = (parm.query_param( "ga_mutate_addition", "yes", "yes no" ) == "yes");
        ga_mutate_deletion = (parm.query_param( "ga_mutate_deletion", "yes", "yes no" ) == "yes");
        ga_mutate_substitution = (parm.query_param( "ga_mutate_substitution", "yes", "yes no" ) == "yes");
        ga_mutate_replacement = (parm.query_param( "ga_mutate_replacement", "yes", "yes no" ) == "yes");
        // use dn roulette - LEP
        //if (ga_mutate_addition || ga_mutate_deletion || ga_mutate_substitution || ga_mutate_replacement){
            //ga_use_dn_roulette = (parm.query_param("ga_use_dn_roulette", "no", "yes no") == "yes");
            ga_use_dn_roulette = false;
            // If mutating parents, set mutation rate
            ga_mutate_parents = (parm.query_param( "ga_mutate_parents", "no", "yes no" ) == "yes");
            if (ga_mutate_parents == true){
                ga_pmut_rate = atof(parm.query_param( "ga_pmut_rate", "0.3" ).c_str());
                //LEP 09.12.18
                if (ga_pmut_rate >= 1){
                    cout << "ERROR: Parent Mutation Rate must be a less than 1" << endl;
                    exit(0);
                }
            }

            // Set offspring mutation rate
            if ( ga_xover_sampling_method_rand == true ){
                ga_omut_rate = atof(parm.query_param( "ga_omut_rate", "0.7" ).c_str());
            }

            // Set number of mutation cycle 
            ga_max_mut_cycles = atoi(parm.query_param( "ga_max_mut_cycles", "5" ).c_str());
    
            // De novo sampling method parameters
            ga_mut_sampling_method = parm.query_param("ga_mut_sampling_method", "rand", "rand | graph");
            ga_mut_sampling_method_rand = false;
            ga_mut_sampling_method_graph = false;

            if (ga_mut_sampling_method.compare("rand") == 0){
                ga_mut_sampling_method_rand = true;
            } else if (ga_mut_sampling_method.compare("graph") == 0){
                ga_mut_sampling_method_graph = true;
            } else {
                cout <<"You chose...poorly." <<endl;
                exit(0);
            }

            // If you are using DN random method, limit number of fragment picks
            if (ga_mut_sampling_method_rand){
                // Default in De novo is 20, so the GA default will be lower
                ga_num_random_picks = atoi(parm.query_param("ga_num_random_picks", "15").c_str());
            }

            // If using DN graph sampling method
            if (ga_mut_sampling_method_graph){
                // The number of random starting points in the graph - DN default = 10
                ga_graph_max_picks = atoi(parm.query_param("ga_graph_max_picks", "10").c_str());
                // The number of children explored relative to a starting point - DN default = 5
                ga_graph_breadth = atoi(parm.query_param("ga_graph_breadth", "3").c_str());
                // The number of generations, e.g. 2=explore children, then children's children
                ga_graph_depth = atoi(parm.query_param("ga_graph_depth", "2").c_str());
                // A starting value for the graph temp [ P=exp(-(Ef-Ei)/T) ]
                ga_graph_temperature = atof(parm.query_param("ga_graph_temperature", "100").c_str());
            }
    
            // Max root size for DN sampling
            ga_max_root_size = atoi(parm.query_param("ga_max_root_size", "5").c_str());
        
       // }else{

         //   ga_use_dn_roulette = false;
         //   ga_mutate_parents = false;
         //   ga_mut_sampling_method_rand = false;
         //   ga_mut_sampling_method_graph = false;
         //   ga_max_mut_cycles = 5;
       // }
    }//end of ga mutations on parm loop


    // BLOCK 5: PRUNING
    // Energy cutoff (in grid or continuous energy) that will be used for conformer pruning
    ga_energy_cutoff = atof(parm.query_param( "ga_energy_cutoff", "100" ).c_str());

    // RMSD horizontal pruning heuristic using Hungarian algorithm
    ga_heur_unmatched_num = atoi(parm.query_param("ga_heur_unmatched_num", "1").c_str());
    ga_heur_matched_rmsd = atof(parm.query_param("ga_heur_matched_rmsd", "0.5").c_str());

    // Based on Leeson et al, Nature., 2012, 481, pp 455-456
    // Upper bounds for molecular weight
    ga_constraint_mol_wt = atoi(parm.query_param("ga_constraint_mol_wt", "500").c_str());

    // Based on Veber et al, J. Med. Chem., 2002, 45(12), pp 2615-2623
    // Upper bounds for rotatable bonds
    ga_constraint_rot_bon = atoi(parm.query_param("ga_constraint_rot_bon", "10").c_str());

    // Upper bounds for hydrogen acceptors
    ga_constraint_H_accept = atoi(parm.query_param("ga_constraint_H_accept", "10").c_str());

    // Upper bounds for hydrogen donors
    ga_constraint_H_donor = atoi(parm.query_param("ga_constraint_H_don", "5").c_str());

    // Upper bounds for formal charge. The negative will used as the lower bounds
    ga_user_constraint_formal_charge = atoi(parm.query_param("ga_constraint_formal_charge", "2").c_str());
    ga_constraint_formal_charge = (ga_user_constraint_formal_charge +0.1);


    // BLOCK 6: SELECTION METHODS
    // Max parent ensemble size
    ga_ensemble_size = atoi(parm.query_param("ga_ensemble_size", "200").c_str());

    // Choice selection method
    // metropolis/sus hidden untested/broken LEP
    //ga_selection_method = parm.query_param("ga_selection_method", "elitism", "elitism | tournament | roulette | sus | metropolis");
    ga_selection_method = parm.query_param("ga_selection_method", "elitism", "elitism | tournament | roulette");
   
    // Initialize options
    ga_selection_method_elitism = false;
    ga_selection_method_tournament = false;
    ga_selection_method_roulette = false;
    ga_selection_method_sus = false;
    ga_selection_method_metropolis = false;

    if (ga_selection_method.compare("elitism") == 0){
        ga_selection_method_elitism = true;
    } else if (ga_selection_method.compare("tournament") == 0){
        ga_selection_method_tournament = true;
    } else if (ga_selection_method.compare("roulette") == 0){
        ga_selection_method_roulette = true;
    } else if (ga_selection_method.compare("sus") == 0){
        ga_selection_method_sus = true;
    } else if (ga_selection_method.compare("metropolis") == 0){
        ga_selection_method_metropolis = true;
    } else {
        cout <<"You chose...poorly." <<endl;
        exit(0);
    }

    input_parameters_selection(parm);

    // Diversity
    //ga_niching = (parm.query_param( "ga_niching", "no", "yes no" ) == "yes");
    ga_niching = false; //JDB disabling for 6.10
    if (ga_niching){
       // If Elitism or Tournament, there are two options
       if (ga_selection_method_elitism || ga_selection_method_tournament){
          // Choose fitness sharing or crowding method
          ga_niching_method = parm.query_param("ga_niching_method", "sharing", "sharing | crowding");
          if (ga_niching_method.compare("sharing") == 0){
             ga_niche_sharing = true;
          }else{ ga_niche_sharing = false; }
       }
       else{ ga_niche_sharing = true; }
    }
    //ga_selection_extinction = (parm.query_param( "ga_selection_extinction", "no", "yes no" ) == "yes");
    ga_selection_extinction = false; //JDB disabling for 6.10
    if (ga_selection_extinction){
       // Number of generation to turn on - counter is not interrupted by generneations of extinction
       ga_extinction_frequency = atoi(parm.query_param("ga_extinction_frequency", "200").c_str());
       // Number of generations to remain on
       ga_extinction_duration = atoi(parm.query_param("ga_extinction_duration", "50").c_str());              
       // Number of molecules to keep - recommend 10% of cutoff
       ga_extinction_keep = atoi(parm.query_param("ga_extinction_keep", "20").c_str());              

       // If you want to change mutation type 
       // NOTE: You have to change the entire type, cannot be a change in the subtype of the save primary method
       ga_extinction_switch_selection = (parm.query_param( "ga_extinction_switch_selection", "no", "yes no" ) == "yes");

       if (ga_extinction_switch_selection){
          //ga_secondary_selection_method = parm.query_param("ga_secondary_selection_method", "elitism", "elitism | tournament | roulette | sus");
          ga_secondary_selection_method = parm.query_param("ga_secondary_selection_method", "elitism", "elitism | tournament | roulette");
          // Initialize options
          ga_secondary_selection_method_elitism = false;
          ga_secondary_selection_method_tournament = false;
          ga_secondary_selection_method_roulette = false;
          ga_secondary_selection_method_sus = false;
          ga_secondary_selection_method_metropolis = false;

          if (ga_secondary_selection_method.compare("elitism") == 0){
             ga_secondary_selection_method_elitism = true;
          } else if (ga_secondary_selection_method.compare("tournament") == 0){
             ga_secondary_selection_method_tournament = true;
          } else if (ga_secondary_selection_method.compare("roulette") == 0){
             ga_secondary_selection_method_roulette = true;
          } else if (ga_secondary_selection_method.compare("sus") == 0){
             ga_secondary_selection_method_sus = true;
          } else if (ga_secondary_selection_method.compare("metropolis") == 0){
             ga_secondary_selection_method_metropolis = true;
          } else {
	         cout <<"You chose...poorly." <<endl;
             exit(0);
          }
          input_parameters_selection(parm);
       }
    }

    ga_no_xover_exit_gen = atoi(parm.query_param("ga_max_num_gen_with_no_crossover", "25").c_str()); //2018.02.25 - YZ
    ga_name_identifier = parm.query_param( "ga_name_identifier", "ga" );

    ga_output_prefix = parm.query_param("ga_output_prefix", "ga_output"); //2018.01.23 - LEP 
    return;

} // end GA_Recomb::input_parameters()




// +++++++++++++++++++++++++++++++++++++++++
// Read parameters from the dock.in file
void
GA_Recomb::input_parameters_selection( Parameter_Reader & parm )
{
    // If selecting the top scored molecules as survivors,
    if ( (ga_selection_method_elitism && !ga_extinction_switch_selection) || ga_secondary_selection_method_elitism){
        ga_elitism_combined = (parm.query_param( "ga_elitism_combined", "yes", "yes no" ) == "yes");
        
        // Choose an option and initialize
        ga_elitism_option = parm.query_param("ga_elitism_option", "max", "percent | number | max");
        
        if (ga_elitism_option.compare("percent") == 0){
            ga_elitism_percent = atof(parm.query_param("ga_elitism_percent", "0.20").c_str());
            ga_elitism_number = 0;
        } else if (ga_elitism_option.compare("number") == 0){
            // Number of ancestrual parents maintained - This does not limit ensemble size
            ga_elitism_number = atoi(parm.query_param("ga_elitism_number", "20").c_str());
            ga_elitism_percent = 0.0;
        } else if (ga_elitism_option.compare("max") == 0){
            ga_elitism_number = ga_ensemble_size;
            ga_elitism_percent = 0.0;
        } else{
            cout << "No elitism selection method option was specified. Try again." << endl;
            exit(0);
        }
    }

    // If completing tournament, there are two options: (1) Parents Vs Children to compete & (2) Parent Vs Parent + Child Vs Child
    if ( (ga_selection_method_tournament && !ga_extinction_switch_selection) || ga_secondary_selection_method_tournament){
       ga_tournament_p_vs_c = (parm.query_param( "ga_tournament_p_vs_c", "yes", "yes no" ) == "yes");
    }

    // If compling roulette, there are two options: (1) Parents Vs Children to compete & (2) Parent Vs Parent + Child Vs Child
    if ( (ga_selection_method_roulette && !ga_extinction_switch_selection) || ga_secondary_selection_method_roulette){
       ga_roulette_separate = (parm.query_param( "ga_roulette_separate", "yes", " yes no" ) == "yes");
    }

    // If compling roulette, there are two options: (1) Parents Vs Children to compete & (2) Parent Vs Parent + Child Vs Child
    if ( (ga_selection_method_sus && !ga_extinction_switch_selection) || ga_secondary_selection_method_sus){
       ga_sus_separate = (parm.query_param( "ga_sus_separate", "yes", " yes no" ) == "yes");
    }

    // If using metropolis,
    if ( (ga_selection_method_metropolis && !ga_extinction_switch_selection) || ga_secondary_selection_method_metropolis){
       // A starting value for the graph temp [ P=exp(-(Ff-Fi)/(k*T)) ]
       // 373 kelvin is the boiling point of sea water
       ga_mc_initial_temp = atof(parm.query_param("ga_mc_intial_temp", "500").c_str());
       // Temperature is calculated via an expotential decay function: Tf = Ti(1-R)^(generation) 
       // Where R = c/max_generations
       ga_mc_rate_constant = atof(parm.query_param("ga_mc_rate_constant", "2").c_str());
       ga_mc_k = atof(parm.query_param("ga_mc_k", "1").c_str());
       ga_metropolis_separate = (parm.query_param( "ga_metropolis_separate", "yes", " yes no" ) == "yes");
    }
    // For verbose statistics - LEP 02.27.2018
    verbose = 0 != parm.verbosity_level();

    return;

} // end GA_Recomb::input_parameters_selection()




// +++++++++++++++++++++++++++++++++++++++++
// The initialize function tells the user the filenames it is reading and populates that parents library vector
void
GA_Recomb::initialize()
{
    cout << "\nInitializing Genetic Algorithms Routines...\n";

    cout << " Reading the molecule library from " << ga_molecule_file <<"...";
    read_library( parents, ga_molecule_file );
    cout << "Done (#=" << parents.size() <<")\n";

    // Set DNM booleans, determine if DNM present for this run
    set_DNM_bools_start( parents );
    if (dnm_enabled){
        cout << "DNM atoms were specified in the parent ensemble." << endl;
    }
    DN_GA_Build c_dn;
    cout <<" Reading the scaffold library from " <<ga_fraglib_scaffold_file <<"...";
    c_dn.read_library( tmp_scaffolds, ga_fraglib_scaffold_file );
    cout <<"Done (#=" <<tmp_scaffolds.size() <<")" <<endl;

    cout <<" Reading the linker library from " <<ga_fraglib_linker_file <<"...";
    c_dn.read_library( tmp_linkers, ga_fraglib_linker_file );
    cout <<"Done (#=" <<tmp_linkers.size() <<")" <<endl;

    cout <<" Reading the sidechain library from " <<ga_fraglib_sidechain_file <<"...";
    c_dn.read_library( tmp_sidechains, ga_fraglib_sidechain_file );
    cout <<"Done (#=" <<tmp_sidechains.size() <<")" <<endl;

    // Additional Preparation steps for fragment replacement
    initialize_fraglib(tmp_scaffolds);
    initialize_fraglib(tmp_linkers);

    /*if ( ga_use_torenv_table ){
        cout << " Reading the torenv table from " << ga_torenv_table <<"...";
        read_torenv_table( ga_torenv_table );
    }*/
    return;

} // end GA_Recomb::initialize()




// +++++++++++++++++++++++++++++++++++++++++
// Prepare scaffolds for hopping in mutation - complete this once
// Reference ReCore - NOTE: Not completed
void
GA_Recomb::initialize_fraglib( std::vector<Fragment> & scaffolds )
{

   Trace trace( "GA_Recomb::initialize_fraglib()" );
   // Initialize the size of the scaffold

   // STEP 1: Define scaffold size 
   for (int i = 0; i<scaffolds.size(); i++){

       // Number of atoms is already defined. We need to know how many are in a ring and how many are dummy
       // DELETE
       //cout <<endl <<"####MOLECULE # " << i << endl << "num atoms " << scaffolds[i].mol.num_atoms <<endl;

       // Ensure that atom_ring_flag is being properly populated
       scaffolds[i].mol.id_ring_atoms_bonds();

       // Define size of heavy atoms only (non-ring) 
       for (int j=0; j<scaffolds[i].mol.num_atoms;j++){
           //cout << " scaffold flag " << i << " " <<tmp_scaffolds[j].mol.atom_ring_flags[i] <<endl;
           if ((scaffolds[i].mol.atom_types[j] != "H") && (scaffolds[i].mol.atom_types[j] != "Du")){
              // If the atoms are not in a ring (for heavy atoms only)
              // Check for rings
              if (scaffolds[i].mol.atom_ring_flags[j] == true){
                 scaffolds[i].ring_size++;
              }
              else{ scaffolds[i].size++; }
           }
       }//cout << "size " << scaffolds[i].size << endl; //DELETE
       //cout << "ring_size " << scaffolds[i].ring_size << endl; 


       //------------------------------------------------------------------------
       // STEP 2: Define distances between each aps

       // Detemine the size of the aps dist matrix
       std::vector < std::vector <float> > tmp_aps_dist(scaffolds[i].aps.size(), std::vector <float> (scaffolds[i].aps.size()));

       // For each aps
       for (int j=0; j<scaffolds[i].aps.size(); j++){ 
           for (int k=j+1; k<scaffolds[i].aps.size(); k++){ 
               DOCKVector vec1;
               vec1.x = scaffolds[i].mol.x[scaffolds[i].aps[j].heavy_atom] - scaffolds[i].mol.x[scaffolds[i].aps[k].heavy_atom];
               vec1.y = scaffolds[i].mol.y[scaffolds[i].aps[j].heavy_atom] - scaffolds[i].mol.y[scaffolds[i].aps[k].heavy_atom];
               vec1.z = scaffolds[i].mol.z[scaffolds[i].aps[j].heavy_atom] - scaffolds[i].mol.z[scaffolds[i].aps[k].heavy_atom];

               // Measure the distance between the heavy atoms of the attachment points
               float dist1 = (vec1.x*vec1.x) +  ( vec1.y*vec1.y) + ( vec1.z*vec1.z);
               tmp_aps_dist[j][k] = tmp_aps_dist[k][j] = dist1;
               // DELETE
           /*cout <<"vec1.x = " << vec1.x << " vec1.y " << vec1.y << " vec1.z " << vec1.z << " atom type j " << scaffolds[i].mol.atom_types[scaffolds[i].aps[j].heavy_atom] << " k " << scaffolds[i].mol.atom_types[scaffolds[i].aps[k].heavy_atom] << endl;
           cout << " distances k " << scaffolds[i].mol.x[scaffolds[i].aps[k].heavy_atom] << " " << scaffolds[i].mol.y[scaffolds[i].aps[k].heavy_atom] << " " << scaffolds[i].mol.z[scaffolds[i].aps[k].heavy_atom] << endl;
           cout << " distances j " << scaffolds[i].mol.x[scaffolds[i].aps[j].heavy_atom] << " " << scaffolds[i].mol.y[scaffolds[i].aps[j].heavy_atom] << " " << scaffolds[i].mol.z[scaffolds[i].aps[j].heavy_atom] << endl;*/
         }
     }
     // Save date, change to sum
     scaffolds[i].aps_dist = tmp_aps_dist;

     // DELETE
     /*for (int j=0; j<scaffolds[i].aps.size(); j++){ 
         for (int k=j+1; k<scaffolds[i].aps.size(); k++){ 
             cout << "IN APS DIST " << scaffolds[i].aps_dist[j][k] << " " << scaffolds[i].aps_dist[k][j] <<endl;
         }
     }*/
   }
   return;
}//end GA_Recomb::inititalize_scaffolds()




// +++++++++++++++++++++++++++++++++++++++++
// Get the internal energy parameters from Master_Conformer_Search
void
GA_Recomb::initialize_internal_energy_parms( bool uie, int rep_exp, int att_exp, float diel, float iec )
{
   Trace trace( "GA_Recomb::initialize_internal_energy_parms()" );
   // These parameters will be defined before the minimzer is called
   std::string              ga_molecule_file;     //filename of molecules taken from the dock.in file    

   use_internal_energy = uie;         // boolean, use internal energy in ligand minimization?
   ie_rep_exp          = rep_exp;     // int. energy vdw repulsive exponent (default 12)
   ie_att_exp          = att_exp;     // int. energy vdw attractive exponent (6)
   ie_diel             = diel;        // int. energy dielectric (4.0)
   ie_cutoff           = iec;         // BCF internal energy cutoff

   return;
}//end GA_Recomb::inititalize_internal_energy_parms()




// +++++++++++++++++++++++++++++++++++++++++
// Reads parents library file and saves DOCKMOL objects in parents vector
void
GA_Recomb::read_library( vector <DOCKMol> & parents, string filename )
{
    Trace trace( "GA_Recomb::read_library()" );
    // Create temporary DOCKMol object
    DOCKMol temp_mol;

    // Create filehandle for reading mol2 library
    ifstream fin_mols;

    // Open the filehandle for reading parents library
    fin_mols.open( filename.c_str() );
    if ( fin_mols.fail() ){
        cout <<"Could not open " <<filename <<". Program will terminate.\n";
        exit(0);
    }

    // Read the mol2 file and copy DOCKMol objects into a vector DOCKMol objects
    while ( Read_Mol2(temp_mol, fin_mols, false, false, false) ){

        // Add the complete DOCKMol object to the appropriate vector
        parents.push_back( temp_mol ); 
    }


    // Close the filehandle
    fin_mols.close();



    return;

} // end GA_Recomb::read_library()




// +++++++++++++++++++++++++++++++++++++++++
// Main function that will call the breeding function  based on the user-specified number of generations
// and fitness parameters - called in dock.cpp
void
GA_Recomb::max_breeding( Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & c_typer, Orient & orient)
{
    cout <<endl<< "----- ENTERING DOCK GENETIC ALGORITHM ------" << endl;
    //cout << ""
    Trace trace( "GA_Recomb::max_breeding()" );
    double start_time;
    double stop_time;

    // The following if statement will seg fault if prepare_molecule is not called first to initialize the values.
    c_typer.prepare_molecule( parents[0], true, score.use_chem, score.use_ph4, score.use_volume );
    // This will check if the dummy atoms have vdw param = 0 and will exit if so
    if ( (c_typer.atom_typer.types[c_typer.atom_typer.atom_types[0]].radius == 0.0) || (c_typer.atom_typer.types[c_typer.atom_typer.atom_types[0]].well_depth == 0.0) ){
         cout <<"Dummy atoms must have well depth > 0 and radius > 0 . Use vdw_de_novo.defn in parameters folder. Program will terminate." <<endl;
         exit(0);
    }

    srand(simplex.random_seed);
    //int lauren_rand = (simplex.random_seed + i)*10;
    //cout << "GA Random Seed:" << lauren_rand << endl;
    double start_time_selection;
    double stop_time_selection;
    //int lauren_rand = 0;
    //cout << "GA random seed: " << lauren_rand << endl;
    // STEP 0: Initial MC temperature, if being used
    // Use generation to vary the temperature
    int gen_counter = 0;
    if (ga_selection_method_metropolis){
       // Temperature is calculated via an expotential decay function: Tf = Ti(1-R)^(generation) 
       // Where R = c/max_generations
       float rate = ga_mc_rate_constant/max_generations;
       new_temp = ga_mc_initial_temp*pow(float(1 - rate), gen_counter);
       cout << "#### Metropolis Round 1: new temp " << new_temp<< endl;
    }

    int extinction_cycle = 1;
    int extinction_duration = 0;
    // Define a preset vector for on off extinction time
    // NOTE: This should start based on user in put or be based on the 1st cycle of selection
    std::vector  <bool>  extinction_on;
    if (ga_selection_extinction){
       for (int i=0; i<max_generations; i++){
          // If extinction starts 
          if ( i+1 == ga_extinction_frequency * extinction_cycle){
             extinction_on.push_back(1);
             extinction_cycle++;
             extinction_duration = 1; //reset for each iteration
          }
          // Else turn off
          else{ extinction_on.push_back(0); }
          // Update with duration
          if (( i > 0 ) && (extinction_on[i-1] == 1) && (extinction_duration < ga_extinction_duration)){
             extinction_on[i] = 1;
             extinction_duration++;
          }
       }       
    }
    for (int i=0; i<extinction_on.size(); i++){
       if (verbose) cout << "extinction " << i << " " << extinction_on[i] <<endl;
    }    

    // Initialize the selection methods, if using extinction & changing type
    std::vector  <int>  extinction_selection_type;
    if (ga_extinction_switch_selection){
       // Save selection type based on number
       //1: Elitism, 2: Tournament, 3: Roulette, 4:sus, 5:metropolis 
       for (int i=0; i<extinction_on.size(); i++){
          if ( extinction_on[i] == 0 ){
             if ( ga_selection_method_elitism ){
                extinction_selection_type.push_back(1);
             }
             if ( ga_selection_method_tournament ){
                extinction_selection_type.push_back(2);
             }
             if ( ga_selection_method_roulette ){
                extinction_selection_type.push_back(3);
             }
             if ( ga_selection_method_sus ){
                extinction_selection_type.push_back(4);
             }
             if ( ga_selection_method_metropolis ){
                extinction_selection_type.push_back(5);
             }
          }
          else{
             if ( ga_secondary_selection_method_elitism ){
                extinction_selection_type.push_back(1);
             }
             if ( ga_secondary_selection_method_tournament ){
                extinction_selection_type.push_back(2);
             }
             if ( ga_secondary_selection_method_roulette ){
                extinction_selection_type.push_back(3);
             }
             if ( ga_secondary_selection_method_sus ){
                extinction_selection_type.push_back(4);
             }
             if ( ga_secondary_selection_method_metropolis ){
                extinction_selection_type.push_back(5);
             }
          }
       }       
    }
    //-------------------------------------------------------------------------
    // STEP 1: Iterate over user specified number of generations
    for (int i=0; i<max_generations; i++){  
    //cout << "GA random seed: " << lauren_rand << endl;	
	//Pass some variables to simplex for debugging LEP
	simplex.ga_gen = i;
    simplex.simplex_ga_flag = true;
    //srand((simplex.random_seed + i)*10);
    //int lauren_rand = (simplex.random_seed + i)*10;
    //cout << "GA Random Seed:" << lauren_rand << endl;

        // If the first generation
        if ( i == 0 ){
           // GENERATION 1:
           cout << endl << "##### Initializing Initial Parent Generation #####" <<endl;

           // Initialize time counter
           double          start_time = time_seconds();

           // STEP 1: Initialize parents
           // STEP 1A: Before GA runs, the parents names are initialized using the substructure names
     
           // This variable will hold the subst label string
           ostringstream subst_name;
           ostringstream new_title;
           ostringstream res_num;

           // Reset Atom residue numbers
           res_num << "1";

           // For each parent
           for (int j=0; j<parents.size(); j++){
               // The new subst label will be "p" and the number of the parent in the list starting with p1     
               subst_name << "p" << j+1;
       
	       if (!ga_utilities){ 
                   // Change title
                   new_title << "p" << j+1;
                   parents[j].title = new_title.str();
	       }
               // Add the label to each atom in the parent 
               for (int k=0; k<parents[j].num_atoms; k++){
                   parents[j].subst_names[k] = subst_name.str();
                   // Update residue numbers to prevent issues with DN minimization for parent mutations
                   parents[j].atom_residue_numbers[k] = res_num.str();
               }
               // Clear the label
               subst_name.str("");
               new_title.str("");
           }
           
	   //
 
           for ( int m=0; m<parents.size(); m++){
               if ( ga_use_torenv_table ){
                   prepare_mol_torenv( parents[m] );
                   c_typer.skip_verbose_flag = true;
                   c_typer.prepare_molecule( parents[m], true, score.use_chem, score.use_ph4, score.use_volume );
            
                   calc_descriptors(parents[m]);

                   if (parents[m].mol_wt > ga_constraint_mol_wt ){
                       cout << "Parent " << m<< " did not pass MW check." << endl;
                   }
                   if ( (parents[m].rot_bonds > ga_constraint_rot_bon) ){
                       cout << "Parent "<< m<< " did not pass RB check." << endl;
                   }
                   if ((parents[m].hb_acceptors > ga_constraint_H_accept) ){
                       cout << "Parent " << m << " did not pass HA check." << endl;
                   }
                   if ((parents[m].hb_donors > ga_constraint_H_donor) ){
                       cout << "Parent "<< m << " did not pass HD check." << endl;
                   }
                   if (((parents[m].formal_charge > ga_constraint_formal_charge) || 
                        (parents[m].formal_charge < -ga_constraint_formal_charge)) ){
                        cout << "Parent " << m << " did not pass FC check." << endl;
                   }

                   DN_GA_Build c_dn;
                   Fragment tmp_parent = mol_to_frag(parents[m]);
                   c_dn.dn_ga_flag = true;
                   c_dn.dn_use_roulette = ga_use_dn_roulette;
                   c_dn.dn_use_torenv_table = ga_use_torenv_table;
                   c_dn.dn_torenv_table = ga_torenv_table;
                   c_dn.read_torenv_table(c_dn.dn_torenv_table);
                   c_dn.verbose = 0;
                   c_dn.dn_MW_cutoff_type_hard = true;
                   c_dn.dn_MW_cutoff_type_soft = false;

                   prepare_torenv_indices(tmp_parent);
                   if ( c_dn.valid_torenv_multi( tmp_parent ) ){
                       cout << "Parent " << m<< " passed initial torenv check." << endl;
                   }else{
                       cout << endl<< "WARNING: Parent " << m << " did not pass initial torenv check." << endl;
                       cout << "Consider adding initial parent torsion environments to the torsion table."<< endl<< endl;
                   }
                   
               }
           }

           // STEP 1B: Sort + score Parents        
           score_parents(parents, score, c_typer);

	   //GA utils
           if (ga_calc_parent_pairwise_tan){
           	pairwise_tanimoto_calc(parents);
           }
           if (ga_calc_parent_pairwise_hms){
           	pairwise_hms_calc(parents);
           }



           sort(parents.begin(), parents.end(), mol_sort);
   
           // If only two molecules are parents, DO Not Prune because you want to make analogs
           if (parents.size() > 2){
              unique_parents(parents, score, simplex, c_typer);
           }
     
           // Fitness score if being used
           if ( ga_niching ){
               if ( ga_niche_sharing ){
                  niche_sharing( parents, score );
               }
               else{
                  niche_crowding( parents, score );
               }
           }

           // Print unique molecules to file with appropriate descriptors
           ostringstream fout_prefix_name;
           fout_prefix_name << ga_output_prefix << ".restart0000.mol2";
           fstream fout_molecules;
           fout_molecules.open (fout_prefix_name.str().c_str(), fstream::out|fstream::app );
           for (int j=0; j<parents.size(); j++){
               // Only Update name if a new molecule
               naming_function(parents[j], i, j);
               // Print necessary descriptors
               calc_descriptors(parents[j]);
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Name:" 
                              << setw(FLOAT_WIDTH) << parents[j].title <<endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Generation:" 
                              << setw(FLOAT_WIDTH) << i <<endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Molecular_Weight:" 
                              << setw(FLOAT_WIDTH) << fixed << setprecision(3) << parents[j].mol_wt <<endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "DOCK_Rotatable_Bonds:" 
                              << setw(FLOAT_WIDTH) << parents[j].rot_bonds <<endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Hydrogen_Acceptors:" 
                              << setw(FLOAT_WIDTH) << parents[j].hb_acceptors <<endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Hydrogen_Donors:" 
                              << setw(FLOAT_WIDTH) << parents[j].hb_donors <<endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Formal_Charge:" 
                              << setw(FLOAT_WIDTH) << fixed << setprecision(3) << round(parents[j].formal_charge) <<endl;

               #ifdef BUILD_DOCK_WITH_RDKIT
               RDTYPER rdprops;
               rdprops.calculate_descriptors(parents[j], ga_fragMap, true, ga_PAINSMap);

               std::vector<std::string> vecpnstmp{};
               vecpnstmp = parents[j].pns_name;
               std::string molpns_name = std::accumulate(vecpnstmp.begin(), vecpnstmp.end(), std::string(""));

               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "num_arom_rings:"
                              << setw(FLOAT_WIDTH) << parents[j].num_arom_rings << endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "num_alip_rings:"
                              << setw(FLOAT_WIDTH) << parents[j].num_alip_rings << endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "num_sat_rings:"
                              << setw(FLOAT_WIDTH) << parents[j].num_sat_rings << endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Stereocenters:"
                              << setw(FLOAT_WIDTH) << parents[j].num_stereocenters << endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Spiro_atoms:"
                              << setw(FLOAT_WIDTH) << parents[j].num_spiro_atoms << endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "TPSA:"
                              << setw(FLOAT_WIDTH) << parents[j].tpsa << endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "cLogP:"
                              << setw(FLOAT_WIDTH) << parents[j].clogp << endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "ESOL:"
                              << setw(FLOAT_WIDTH) << parents[j].esol << endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "QED:"
                              << setw(FLOAT_WIDTH) << parents[j].qed_score << endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "SA_Score:"
                              << setw(FLOAT_WIDTH) << parents[j].sa_score << endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "SMILES: "
                              << setw(FLOAT_WIDTH) << parents[j].smiles << endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "num_of_PAINS:"
                              << setw(FLOAT_WIDTH) << parents[j].pns << endl;
               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "PAINS_names:"
                              << setw(FLOAT_WIDTH) << molpns_name << endl;
               #endif
	       // Calculate the fitness score for the initial ensemble if niching is being used
               if ( ga_niching ){
                  if ( ga_niche_sharing ){
                      fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Fitness_Score:" 
                                     << setw(FLOAT_WIDTH) << fixed << setprecision(3) << parents[j].fitness <<endl;
                  }
                  else{
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Rank:"
                                    << setw(FLOAT_WIDTH) << parents[j].rank <<endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Crowding_Distance:" 
                                    << setw(FLOAT_WIDTH) << fixed << setprecision(3) << parents[j].crowding_dist <<endl;
                 }
               }
               //fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Parents:" 
               //               << setw(FLOAT_WIDTH) << parents[j].energy <<endl;
               fout_molecules <<parents[j].current_data << endl;
               Write_Mol2(parents[j], fout_molecules);
           }
           fout_molecules.close();
           fout_prefix_name.clear();

           // Average score of parents
           float score_parent = 0.00;
           for (int j=0; j<parents.size(); j++){
               score_parent = score_parent + parents[j].current_score;
           }

           // Print average parents score
           float ave_parent_score = score_parent/parents.size();
           cout << "Average score of parents for generation " << i+1 << ": " << ave_parent_score << endl;
            cout << "For generation " << i+1 << ", new size of parents: " << parents.size() << endl;

        }// end if Initial generation

        // If not the first gen, check to make sure there are more than 2 parents
        if ( (( i > 0 ) && ( parents.size() >= 1)) || ( i == 0 ) ){
 
            cout <<endl << "##### Entering Generation " << i+1 << " #####"<<endl;
            //cout << "GA random seed: " << lauren_rand << endl;
            start_time = time_seconds();
            if (parents.size() == 1) cout << "WARNING: parent ensemble size is 1" <<endl;
            /// Print size of new parent vector
            cout << "For generation " << i+1 << ", new size of parents: " << parents.size() <<endl;

            // STEP 2: Run the GA
            // Clear DOCKMol vectors that are filled from breeding
            pruned_children.clear();
            scored_generation.clear();
            children.clear();
            mutants.clear();

            if ( ga_check_overlap && ( ga_check_only || ga_xover_sampling_method_rand ) ){
                check_exhaustive(score, simplex, c_typer, orient);
                // Save the final unique vector to file
                ostringstream fout_molecules_name;
                fout_molecules_name << ga_output_prefix << "_unique_crossover" << i << ".mol2"; 
                //^2018.01.23 - LEP - user defined output 
                fstream fout_molecules;
                fout_molecules.open (fout_molecules_name.str().c_str(), fstream::out|fstream::app);

                for (int j=0; j<xover_parents.size(); j++){
                    activate_mol(xover_parents[j]);
                    calc_descriptors(xover_parents[j]);
                    fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Name:" 
                                   << setw(FLOAT_WIDTH) << xover_parents[j].title <<endl;
                    fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Generation:" 
                                   << setw(FLOAT_WIDTH) << i << endl;
                    fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Molecular_Weight:" 
                                   << setw(FLOAT_WIDTH) << fixed << setprecision(3) << xover_parents[j].mol_wt <<endl;
                    fout_molecules << DELIMITER << setw(STRING_WIDTH) << "DOCK_Rotatable_Bonds:" 
                                   << setw(FLOAT_WIDTH) << xover_parents[j].rot_bonds <<endl;
                    fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Hydrogen_Acceptors:" 
                                   << setw(FLOAT_WIDTH) << xover_parents[j].hb_acceptors <<endl;
                    fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Hydrogen_Donors:" 
                                   << setw(FLOAT_WIDTH) << xover_parents[j].hb_donors <<endl;
                    fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Formal_Charge:" 
                                   << setw(FLOAT_WIDTH) << fixed << setprecision(3) << round(xover_parents[j].formal_charge) <<endl;

	            // Calculate the fitness score for the initial ensemble if niching is being used
                    if ( ga_niching ){
                       if ( ga_niche_sharing ){
                          fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Fitness_Score:" 
                                         << setw(FLOAT_WIDTH) << fixed << setprecision(3) << xover_parents[j].fitness <<endl;
                       }
                       else{
                          fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Rank:" 
                                         << setw(FLOAT_WIDTH) << xover_parents[j].rank <<endl;
                          fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Crowding Distance:" 
                                         << setw(FLOAT_WIDTH) << fixed << setprecision(3) << xover_parents[j].crowding_dist <<endl;
                       }
                   }
                   fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Status:" 
                                  << setw(FLOAT_WIDTH) << xover_parents[j].energy <<endl;
                   fout_molecules <<xover_parents[j].current_data << endl;
                   Write_Mol2(xover_parents[j], fout_molecules);
               }
               xover_parents.clear();
            }

            // Call breeding function
            if ( !ga_check_only ){
               if (ga_xover_sampling_method_rand){
                  cout << endl<< "#### Entering Random Breeding Routine ####" << endl;
                  breeding_rand(score, simplex, c_typer, orient, i);
               }
               else{ 
                  if ( ga_check_overlap ){
                     cout <<endl <<  "#### Entering Exhaustive Breeding Routine ####" << endl;
                     breeding_exhaustive(score, simplex, c_typer, orient, i);
                     // Save the final unique vector to file
                     ostringstream fout_molecules_name;
                     //2018.01.23  - LEP - user defined output files
                     fout_molecules_name << ga_output_prefix << "_unique_crossover" << i << ".mol2";
                     fstream fout_molecules;
                     fout_molecules.open (fout_molecules_name.str().c_str(), fstream::out|fstream::app);

                     for (int j=0; j<xover_parents.size(); j++){
                         calc_descriptors(xover_parents[j]);
                         fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Name:" 
                                        << setw(FLOAT_WIDTH) << xover_parents[j].title <<endl;
                         fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Generation:"
                                        << setw(FLOAT_WIDTH) << i << endl;
                         fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Molecular_Weight:"
                                        << setw(FLOAT_WIDTH) << fixed << setprecision(3) << xover_parents[j].mol_wt <<endl;
                         fout_molecules << DELIMITER << setw(STRING_WIDTH) << "DOCK_Rotatable_Bonds:"
                                        << setw(FLOAT_WIDTH) << xover_parents[j].rot_bonds <<endl;
                         fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Hydrogen_Acceptors:"
                                        << setw(FLOAT_WIDTH) << xover_parents[j].hb_acceptors <<endl;
                         fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Hydrogen_Donors:"
                                        << setw(FLOAT_WIDTH) << xover_parents[j].hb_donors <<endl;
                         fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Formal_Charge:"
                                        << setw(FLOAT_WIDTH) << fixed << setprecision(3) << round(xover_parents[j].formal_charge) <<endl;

              	         // Calculate the fitness score for the initial ensemble if niching is being used
                         if ( ga_niching ){
                            if ( ga_niche_sharing ){
                               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Fitness_Score:"
                                              << setw(FLOAT_WIDTH) << fixed << setprecision(3) << xover_parents[j].fitness <<endl;
                            }
                            else{
                               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Rank:"
                                              << setw(FLOAT_WIDTH) << xover_parents[j].rank <<endl;
                               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Crowding Distance:"
                                              << setw(FLOAT_WIDTH) << fixed << setprecision(3) << xover_parents[j].crowding_dist <<endl;
                            }
                         }
                         fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Status:"
                                        << setw(FLOAT_WIDTH) << xover_parents[j].energy <<endl;
                         fout_molecules <<xover_parents[j].current_data << endl;
                         Write_Mol2(xover_parents[j], fout_molecules);
                     }
                     xover_parents.clear();
                  }
               }
           }
        

           // STEP 3: Print to screen
           // STEP 3A: Save offspring information
           if ( scored_generation.size() > 0 ){
              float score_counter = 0.00;
              score_counter = 0.00;
              for (int j=0; j<scored_generation.size(); j++){
                  score_counter = score_counter + scored_generation[j].current_score;
              }
              float ave_score = score_counter/scored_generation.size();
              cout << "Average Score offspring gen " << i+1 << ": " << ave_score << endl;
           }
           else{
              if ( i >= ga_no_xover_exit_gen){
                  cout << "WARNING: THERE ARE NO OFFSPRING (from pruning and crossover) AND THE PROGRAM WILL EXIT)"<<endl;
                  exit(0);
              }
           }
              // STEP 8B: Print molecules to file
              //fstream fout_molecules;
              /*fout_molecules.open ( "zzz.scored_gen.mol2", fstream::out|fstream::app );
              for (int j=0; j<scored_generation.size(); j++){
                  calc_descriptors(scored_generation[j]);
                  fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Name:"
                                 << setw(FLOAT_WIDTH) << scored_generation[i].energy <<endl;
                  fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Generation:"
                                 << setw(FLOAT_WIDTH) << i << endl;
                  fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Molecular Weight:"
                                 << setw(FLOAT_WIDTH) << scored_generation[j].mol_wt <<endl;
                  fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Rotatable Bonds:"
                                 << setw(FLOAT_WIDTH) << scored_generation[j].rot_bonds <<endl;
                  fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Hydrogen Acceptors:"
                                 << setw(FLOAT_WIDTH) << scored_generation[j].hb_acceptors <<endl;
                  fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Hydrogen Donors:"
                                 << setw(FLOAT_WIDTH) << scored_generation[j].hb_donors <<endl;
                  fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Formal Charge:"
                                 << setw(FLOAT_WIDTH) << round(scored_generation[j].formal_charge) <<endl;
                  if ( ga_niching ){
                     if ( ga_niche_sharing ){
                        fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Fitness Score:"
                                       << setw(FLOAT_WIDTH) << scored_generation[j].fitness <<endl;
                     }
                     else{
                        fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Rank:"
                                       << setw(FLOAT_WIDTH) << scored_generation[j].rank <<endl;
                        fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Crowding Distance:"
                                       << setw(FLOAT_WIDTH) << scored_generation[j].crowding_dist <<endl;
                    }
                 }
                 fout_molecules << scored_generation[j].current_data << endl;
                 Write_Mol2( scored_generation[j], fout_molecules );
              }
              fout_molecules.close();*/

              // STEP 4: Perform selection
              // STEP 4A: Update the temperature
              gen_counter++;
              /*if (ga_selection_method_metropolis){
                 // Temperature is calculated via an expotential decay function: Tf = Ti(1-R)^(generation) 
                 // where R = c/max_generations
                 float rate = ga_mc_rate_constant/max_generations;
                 new_temp = ga_mc_initial_temp*pow(float(1 - rate), gen_counter);               
                 cout << "#### Metropolis Round " << i << " new temp: " << new_temp<< endl;
              }*/


              // STEP 4B: Perform selection
              if ( !ga_check_only ){
                 start_time_selection = time_seconds();
                 // Need to rescore parents to correct scoring components for niching only
                 if ( ga_niching ){
                    score_parents(parents, score, c_typer);
                 }

                 // If using extinction
                 if ( ga_selection_extinction ){
                    // Update on off
                    ga_extinction_on = extinction_on[i];
                    // Check to see if the selection methods should be switched
                    if (ga_extinction_switch_selection){
                       // Reset selection method
                       ga_selection_method_elitism = false;
                       ga_selection_method_tournament = false;
                       ga_selection_method_roulette = false;
                       ga_selection_method_sus = false;
                       ga_selection_method_metropolis = false;
                       
                       // Reset selection based on type
                       // 1:elitism, 2:tournament, 3: roulette, 4: sus, 5:metropolis
                       if (extinction_selection_type[i] == 1){
                          ga_selection_method_elitism = true;
                       }
                       if (extinction_selection_type[i] == 2){
                          ga_selection_method_tournament = true;
                       }
                       if (extinction_selection_type[i] == 3){
                          ga_selection_method_roulette = true;
                       }
                       if (extinction_selection_type[i] == 4){
                          ga_selection_method_sus = true;
                       }
                       if (extinction_selection_type[i] == 5){
                          ga_selection_method_metropolis = true;
                       }
                    }
                 }
                 // Currently starting based on first frequency
                 //score_parents(parents, score, c_typer);
                 selection_method(score, c_typer);

                 stop_time_selection = time_seconds();
                 if (verbose) cout << "\n" "-----------------------------------" "\n";
                 if (verbose) cout << " Elapsed time selection:\t" << stop_time_selection - start_time_selection
                      << " seconds for GA Generation "<< i+1<< "\n\n";
                 if (verbose )cout <<endl << " #### " << parents.size() << " molecules in new ensemble " << i+1 << " of GA" <<endl;


                 // STEP 4C: Print average survivor score to screen
                 float score_ensemble = 0.00;
                 for (int j=0; j<parents.size(); j++){
                     score_ensemble = score_ensemble + parents[j].current_score;
                 }

                 // Print average parents score
                 float ave_parent_score = score_ensemble/parents.size();
                 cout << "Average score survivors gen " << i+1 << ": " << ave_parent_score << endl;
                 cout << "\n" "-----------------------------------" "\n";

                 // STEP 5: Print molecules to file 
                 // Survivors
                 /*
                 ostringstream fout_molecules_name69;
                 fstream fout_molecules69;
                 fout_molecules69.open (fout_molecules_name69.str().c_str(), fstream::out|fstream::app);
                 for (int j=0; j<some_molecules_go_here.size(); j++){
                     calc_descriptors(some_molecules_go_here[j]);
                     fout_molecules69 << DELIMITER << setw(STRING_WIDTH) << "Name:"
                                    << setw(FLOAT_WIDTH) << some_molecules_go_here[j].title <<endl;
                     fout_molecules69 << DELIMITER << setw(STRING_WIDTH) << "Generation:"
                                    << setw(FLOAT_WIDTH) << i+1 << endl;
                     fout_molecules69 << DELIMITER << setw(STRING_WIDTH) << "Molecular_Weight:"
                                    << setw(FLOAT_WIDTH) << fixed << setprecision(3) << some_molecules_go_here[j].mol_wt <<endl;
                     fout_molecules69 << DELIMITER << setw(STRING_WIDTH) << "DOCK_Rotatable_Bonds:"
                                    << setw(FLOAT_WIDTH) << some_molecules_go_here[j].rot_bonds <<endl;
                     fout_molecules69 << DELIMITER << setw(STRING_WIDTH) << "Hydrogen_Acceptors:"
                                    << setw(FLOAT_WIDTH) << some_molecules_go_here[j].hb_acceptors <<endl;
                     fout_molecules69 << DELIMITER << setw(STRING_WIDTH) << "Hydrogen_Donors:"
                                    << setw(FLOAT_WIDTH) << some_molecules_go_here[j].hb_donors <<endl;
                     fout_molecules69 << DELIMITER << setw(STRING_WIDTH) << "Formal_Charge:"
                                    << setw(FLOAT_WIDTH) << fixed << setprecision(3) << round(some_molecules_go_here[j].formal_charge) <<endl;
                     fout_molecules69 << DELIMITER << setw(STRING_WIDTH) << "Type:"
                                   << setw(FLOAT_WIDTH) << some_molecules_go_here[j].energy <<endl;
                     fout_molecules69 << some_molecules_go_here[j].current_data << endl;
                     Write_Mol2(some_molecules_go_here[j], fout_molecules69);
                 }
                 fout_molecules69.close();
                 fout_molecules_name69.clear();*/
                 //some_molecules_go_here.clear();
                 //if (lauren_var) { some_molecules_go_here.clear();}
                   
             /*    ostringstream fout_molecules_name699;
                 fout_molecules_name699 << ga_output_prefix << ".good_molecule" <<i+1 <<".mol2";
                 fstream fout_molecules699;
                 fout_molecules699.open (fout_molecules_name699.str().c_str(), fstream::out|fstream::app);
                 for (int j=0; j<good_molecules_go_here.size(); j++){
                     calc_descriptors(good_molecules_go_here[j]);
                     fout_molecules699 << DELIMITER << setw(STRING_WIDTH) << "Name:"
                                    << setw(FLOAT_WIDTH) << good_molecules_go_here[j].title <<endl;
                     fout_molecules699 << DELIMITER << setw(STRING_WIDTH) << "Generation:"
                                    << setw(FLOAT_WIDTH) << i+1 << endl;
                     fout_molecules699 << DELIMITER << setw(STRING_WIDTH) << "Molecular_Weight:"
                                    << setw(FLOAT_WIDTH) << fixed << setprecision(3) << good_molecules_go_here[j].mol_wt <<endl;
                     fout_molecules699 << DELIMITER << setw(STRING_WIDTH) << "DOCK_Rotatable_Bonds:"
                                    << setw(FLOAT_WIDTH) << good_molecules_go_here[j].rot_bonds <<endl;
                     fout_molecules699 << DELIMITER << setw(STRING_WIDTH) << "Hydrogen_Acceptors:"
                                    << setw(FLOAT_WIDTH) << good_molecules_go_here[j].hb_acceptors <<endl;
                     fout_molecules699 << DELIMITER << setw(STRING_WIDTH) << "Hydrogen_Donors:"
                                    << setw(FLOAT_WIDTH) << good_molecules_go_here[j].hb_donors <<endl;
                     fout_molecules699 << DELIMITER << setw(STRING_WIDTH) << "Formal_Charge:"
                                    << setw(FLOAT_WIDTH) << fixed << setprecision(3) << round(good_molecules_go_here[j].formal_charge) <<endl;
                     fout_molecules699 << DELIMITER << setw(STRING_WIDTH) << "Parents:"
                                   << setw(FLOAT_WIDTH) << good_molecules_go_here[j].energy <<endl;
                     fout_molecules699 << good_molecules_go_here[j].current_data << endl;
                     Write_Mol2(good_molecules_go_here[j], fout_molecules699);
                 }
                 fout_molecules699.close();
                 fout_molecules_name699.clear();
                 good_molecules_go_here.clear();
               */  

                 ostringstream fout_molecules_name;
                 // "I want the files to be listed in numerical order on my laptop" - Rob //LEP
                 if ( i < 9 ){
                    fout_molecules_name << ga_output_prefix << ".restart000" <<i+1 <<".mol2";
                 }
                 else if ( i < 99 ){
                    fout_molecules_name << ga_output_prefix << ".restart00" <<i+1 <<".mol2";
                 } 
                 else if ( i < 999 ){
                    fout_molecules_name << ga_output_prefix << ".restart0" <<i+1 <<".mol2";
                 }
                 else {
                    fout_molecules_name << ga_output_prefix << ".restart" <<i+1 <<".mol2";
                 }
                 fstream fout_molecules;
                 fout_molecules.open (fout_molecules_name.str().c_str(), fstream::out|fstream::app);

                 for (int j=0; j<parents.size(); j++){
                     calc_descriptors(parents[j]);
                     // Only Update name if a new molecule
                     if ( parents[j].parent == 0 ){
                        naming_function(parents[j], i+1, j);
                     }
                     //calc_pairwise_distance(parents[j]); 
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Name:"
                                    << setw(FLOAT_WIDTH) << parents[j].title <<endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Generation:"
                                    << setw(FLOAT_WIDTH) << i+1 << endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Molecular_Weight:"
                                    << setw(FLOAT_WIDTH) << fixed << setprecision(3) << parents[j].mol_wt <<endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "DOCK_Rotatable_Bonds:"
                                    << setw(FLOAT_WIDTH) << parents[j].rot_bonds <<endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Hydrogen_Acceptors:"
                                    << setw(FLOAT_WIDTH) << parents[j].hb_acceptors <<endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Hydrogen_Donors:"
                                    << setw(FLOAT_WIDTH) << parents[j].hb_donors <<endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Formal_Charge:"
                                    << setw(FLOAT_WIDTH) << fixed << setprecision(3) << round(parents[j].formal_charge) <<endl;
                     #ifdef BUILD_DOCK_WITH_RDKIT
                     RDTYPER rdprops;
                     rdprops.calculate_descriptors(parents[j], ga_fragMap, true, ga_PAINSMap);

                     std::vector<std::string> vecpnstmp{};
                     vecpnstmp = parents[j].pns_name;
                     std::string molpns_name = std::accumulate(vecpnstmp.begin(), vecpnstmp.end(), std::string(""));

                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "num_arom_rings:" 
                                    << setw(FLOAT_WIDTH) << parents[j].num_arom_rings << endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "num_alip_rings:"
                                    << setw(FLOAT_WIDTH) << parents[j].num_alip_rings << endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "num_sat_rings:" 
                                    << setw(FLOAT_WIDTH) << parents[j].num_sat_rings << endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Stereocenters:"
                                    << setw(FLOAT_WIDTH) << parents[j].num_stereocenters << endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Spiro_atoms:"
                                    << setw(FLOAT_WIDTH) << parents[j].num_spiro_atoms << endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "TPSA:"
                                    << setw(FLOAT_WIDTH) << parents[j].tpsa << endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "cLogP:"
                                    << setw(FLOAT_WIDTH) << parents[j].clogp << endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "ESOL:"
                                    << setw(FLOAT_WIDTH) << parents[j].esol << endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "QED:"
                                    << setw(FLOAT_WIDTH) << parents[j].qed_score << endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "SA_Score:"
                                    << setw(FLOAT_WIDTH) << parents[j].sa_score << endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "SMILES: "
                                    << setw(FLOAT_WIDTH) << parents[j].smiles << endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "num_of_PAINS:"
                                    << setw(FLOAT_WIDTH) << parents[j].pns << endl;
                     fout_molecules << DELIMITER << setw(STRING_WIDTH) << "PAINS_names:"
                                    << setw(FLOAT_WIDTH) << molpns_name << endl;
                     #endif

                     if ( ga_niching ){
                        if ( ga_niche_sharing ){
                           fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Fitness_Score:"
                                          << setw(FLOAT_WIDTH) << fixed << setprecision(3) << parents[j].fitness <<endl;
                        }
                        else{
                           fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Rank:"
                                          << setw(FLOAT_WIDTH) << parents[j].rank <<endl;
                           fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Crowding_Distance:"
                                          << setw(FLOAT_WIDTH) << fixed << setprecision(3) << parents[j].crowding_dist <<endl;
                       }
                    }
                    fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Type:"
                                   << setw(FLOAT_WIDTH) << parents[j].energy <<endl;
                    fout_molecules << parents[j].current_data << endl;
                    Write_Mol2(parents[j], fout_molecules);
                 }
                 fout_molecules.close();
                 fout_molecules_name.clear();
              }//if check only

              // STEP 6: Print Generation timing to screen
              stop_time = time_seconds();
              if (verbose) cout << "\n" "-----------------------------------" "\n";
              if (verbose) cout << " Elapsed time:\t" << stop_time - start_time 
                   << " seconds for GA Generation " << i+1 << "\n\n";
        }// if
        // If less than two parents, exit
        else { exit(0); }

        // Reset the random seed so that it is the same for every anchor
        // This behavior is quasi-consistent with A&G
        //simplex.initialize();
        //srand((simplex.random_seed + i)*10);
        //lauren_rand = (simplex.random_seed + i)*10;
	    //cout << "GA Random Seed:" << lauren_rand << endl; 
        
    }// end for max generation
    return;
} // end GA_Recomb::max_generations()



// +++++++++++++++++++++++++++++++++++++++++
// Determine the molecules with overlapping bonds and returns the set as reference 
void
GA_Recomb::check_exhaustive( Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & c_typer, Orient & orient )
{

    Trace trace( "GA_Recomb::check_exhaustive()" );
    // Step 1: Prepare parents for breeding
    // Activate parents before amber_typer
    for ( int i=0; i<parents.size(); i++ ){
        for ( int j=0; j<parents[i].num_bonds; j++ ){
            parents[i].bond_active_flags[j]=true;
        }
    }

    for ( int i=0; i<parents.size(); i++ ){
        for ( int j=0; j<parents[i].num_atoms; j++ ){
            parents[i].atom_active_flags[j]=true;
        }
    }

    // Prepare amber typer for parents
    for ( int i=0; i<parents.size(); i++ ){
        //cout << "Entering Amber typer" << endl;
        parents[i].prepare_molecule();
        c_typer.skip_verbose_flag = true;
        c_typer.prepare_molecule( parents[i], true, score.use_chem, score.use_ph4, score.use_volume );
        parents[i].used = false;
    } 

    // Deactivate all atoms and then bonds in each parent molecule so that later 
    // only the substructures destined to make the children can be activated
    for ( int i=0; i<parents.size(); i++ ){
        for ( int j=0; j<parents[i].num_bonds; j++ ){
            parents[i].bond_active_flags[j]=false;
            parents[i].bond_keep_flags[j]=false;
        }
    }

    for ( int i=0; i<parents.size(); i++ ){
        for ( int j=0; j<parents[i].num_atoms; j++ ){
            parents[i].atom_active_flags[j]=false;
        }
    }

    // Reset the parents comments (energy)
    for ( int i=0; i<parents.size(); i++ ){
        ostringstream new_energy;
        new_energy << "parent";
        parents[i].energy = new_energy.str();
        parents[i].parent = 1;
    }

    // Step 2: Determine the bonds and coorsponding atoms that overlap
    // Step 2A: Determine if the new bond to be formed during breeding is of the appropriate length

    // Set variables that will hold distances between the origin and target atoms
    // of the bonds on both parents that will be involved in breeding
    float dist13,
          dist14,
          dist23,
          dist24;

    // Calculate the distance between bonds on two different parents, 
    // then determine if the bonds are close enough for breeding

    // Check number of overlapping bonds per molecule
    int num_overlapping_bonds = 0;
    int num_bond_pairs_sampled = 0;

    int pair_counter = 0;
    // For each parent
    for ( int i=0; i<parents.size(); i++ ){
        // For each bond j
        for ( int j=0; j<parents[i].num_bonds; j++ ){
             // Identify if the bond is rotatable
             if ( parents[i].amber_bt_id[j] != -1 ){
                // Compare to the rotatable bond in parent i+1
                for ( int k=i+1; k<parents.size(); k++ ){
                    for ( int m=0; m<parents[k].num_bonds; m++ ){
                        if ( parents[k].amber_bt_id[m] != -1 ){
                           if ( parents[i].bond_types[j] == parents[k].bond_types[m] ){
    // Increment the number of overlapping bonds
    num_bond_pairs_sampled++; 

    // Define distance vectors
    DOCKVector vec13;
    vec13.x = parents[i].x[parents[i].bonds_origin_atom[j]] - parents[k].x[parents[k].bonds_origin_atom[m]];
    vec13.y = parents[i].y[parents[i].bonds_origin_atom[j]] - parents[k].y[parents[k].bonds_origin_atom[m]];
    vec13.z = parents[i].z[parents[i].bonds_origin_atom[j]] - parents[k].z[parents[k].bonds_origin_atom[m]];

    DOCKVector vec14;
    vec14.x = parents[i].x[parents[i].bonds_origin_atom[j]] - parents[k].x[parents[k].bonds_target_atom[m]];
    vec14.y = parents[i].y[parents[i].bonds_origin_atom[j]] - parents[k].y[parents[k].bonds_target_atom[m]];
    vec14.z = parents[i].z[parents[i].bonds_origin_atom[j]] - parents[k].z[parents[k].bonds_target_atom[m]];

    DOCKVector vec23;
    vec23.x = parents[i].x[parents[i].bonds_target_atom[j]] - parents[k].x[parents[k].bonds_origin_atom[m]];
    vec23.y = parents[i].y[parents[i].bonds_target_atom[j]] - parents[k].y[parents[k].bonds_origin_atom[m]];
    vec23.z = parents[i].z[parents[i].bonds_target_atom[j]] - parents[k].z[parents[k].bonds_origin_atom[m]];

    DOCKVector vec24;
    vec24.x = parents[i].x[parents[i].bonds_target_atom[j]] - parents[k].x[parents[k].bonds_target_atom[m]];
    vec24.y = parents[i].y[parents[i].bonds_target_atom[j]] - parents[k].y[parents[k].bonds_target_atom[m]];
    vec24.z = parents[i].z[parents[i].bonds_target_atom[j]] - parents[k].z[parents[k].bonds_target_atom[m]];

    // Measure the distance between the origin and target atoms on the rotatable single bonds in the parents
    // ie the magnitude of the vectors
    dist13 = ( vec13.x*vec13.x) +  ( vec13.y*vec13.y) + ( vec13.z*vec13.z);

    dist14 = ( vec14.x*vec14.x) +  ( vec14.y*vec14.y) + ( vec14.z*vec14.z);

    dist23 = ( vec23.x*vec23.x) +  ( vec23.y*vec23.y) + ( vec23.z* vec23.z);
  
    dist24 = ( vec24.x*vec24.x) +  ( vec24.y*vec24.y) + ( vec24.z*vec24.z);

    // Define vectors that coorespond to the idenfied bond in parent
    DOCKVector vec1;    // vector containing the *direction* from origin to 
                        // target atom in parent 1 around identified bond
    DOCKVector vec2;    // vector containing parent 2 bond information
    float vec1_mag;     // vector magnitude 
    float vec2_mag;

    // Variables to determine if the identified bonds will be used for breeding
    int if_bool = 0;         // boolean that identifies the *direction* of the bonds in parents
    float dot;               // dot product value of vectors vec1 & vec2 
    float cos_theta1;        // *angle* between parents' bonds (used to determine the *direction*) 

    // If all approx distances measured for the new children bonds are within the bond tolerance
    // ie if the two parent bonds are overlapping or nearly overlapped as defined by user defined bond_tolerance
    if ( dist13<ga_bond_tolerance && dist23<ga_bond_tolerance && dist14<ga_bond_tolerance && dist24<ga_bond_tolerance ){

        // Determine which distances provide the best overlap
        if ( (dist13 + dist24) < (dist14 + dist23) ){

            // If the if_bool is true, then the origin-to-origin & target-to-target atom dist is the shortest;
            // Therefore, the *direction* of the origin and target atoms are in the same direction
            if_bool = 1;
        }

        // All other conditions (if (dist14 + dist23) is greater or if the combination of distances are equal) will enter the else condition
        else { if_bool = 2; }
    }
    
    // Else if the origin-to-origin and target-to-target distances are within the tolerance
    else if ( (dist13<ga_bond_tolerance && dist24<ga_bond_tolerance) ){
          if_bool = 1;
    }
    else if ( dist14<ga_bond_tolerance && dist23<ga_bond_tolerance ){
          if_bool = 2;
    }



    //-------------------------------------------------------------------------
    // STEP 2B: Determine if the overlap angle is close to zero and run breeding functions

    // If the if_bool is true, ga_bond_tolerance has been met but the origin & target atoms must be reversed for parent 2
    if ( if_bool == 1 ){
   
         // Populate vectors with the bond length of each parent molecule
         // Note: the vector directions are the same for vec1 and vec2           
         vec1.x = parents[i].x[parents[i].bonds_origin_atom[j]] - parents[i].x[parents[i].bonds_target_atom[j]];  
         vec1.y = parents[i].y[parents[i].bonds_origin_atom[j]] - parents[i].y[parents[i].bonds_target_atom[j]];  
         vec1.z = parents[i].z[parents[i].bonds_origin_atom[j]] - parents[i].z[parents[i].bonds_target_atom[j]];

         vec2.x = parents[k].x[parents[k].bonds_origin_atom[m]] - parents[k].x[parents[k].bonds_target_atom[m]];   
         vec2.y = parents[k].y[parents[k].bonds_origin_atom[m]] - parents[k].y[parents[k].bonds_target_atom[m]];   
         vec2.z = parents[k].z[parents[k].bonds_origin_atom[m]] - parents[k].z[parents[k].bonds_target_atom[m]];
   
         // Magnitude of above vectors
         vec1_mag = ( vec1.x*vec1.x) + ( vec1.y*vec1.y) + ( vec1.z*vec1.z);
         vec2_mag = ( vec2.x*vec2.x) + ( vec2.y*vec2.y) + ( vec2.z*vec2.z);

         // Compute the dot product using the function in utils.cpp
         dot = dot_prod(vec1, vec2);

         // Compute cosine theta
         cos_theta1 = dot / (sqrt (vec1_mag * vec2_mag));

         // If cos_theta1 is close to +1, then the parent bonds are going in the same direction
         if ( cos_theta1 > cos( ga_angle_cutoff ) ){

             // MARK - DELETE: check the cos theta calculated
             //cout << "This bonds passed the overlap criteria with cos theta = " << cos_theta1 <<endl; 
             
             // Increment overlapping bond variable
             num_overlapping_bonds++;

             // Save to unique vector
             if ( parents[i].used == false ){ 
                xover_parents.push_back(parents[i]); 
                parents[i].used = true;
             }
             if ( parents[k].used == false ){ 
                xover_parents.push_back(parents[k]); 
                parents[k].used = true;
             }
         }
    }

    // Repeat if origin-to-target & target-to-origin atom distances is within ga_bond_tolerance
    else if ( if_bool == 2 ){

        // Populate bond vectors
        vec1.x = parents[i].x[parents[i].bonds_origin_atom[j]] - parents[i].x[parents[i].bonds_target_atom[j]];  
        vec1.y = parents[i].y[parents[i].bonds_origin_atom[j]] - parents[i].y[parents[i].bonds_target_atom[j]];  
        vec1.z = parents[i].z[parents[i].bonds_origin_atom[j]] - parents[i].z[parents[i].bonds_target_atom[j]];

        vec2.x = parents[k].x[parents[k].bonds_origin_atom[m]] - parents[k].x[parents[k].bonds_target_atom[m]];   
        vec2.y = parents[k].y[parents[k].bonds_origin_atom[m]] - parents[k].y[parents[k].bonds_target_atom[m]];   
        vec2.z = parents[k].z[parents[k].bonds_origin_atom[m]] - parents[k].z[parents[k].bonds_target_atom[m]];
   
        // Calculate magnitude of above vectors
        vec1_mag = ( vec1.x*vec1.x) + ( vec1.y*vec1.y) + ( vec1.z*vec1.z);
        vec2_mag = ( vec2.x*vec2.x) + ( vec2.y*vec2.y) + ( vec2.z*vec2.z);

        // Compute the dot product using the function in utils.cpp
        dot = dot_prod( vec1, vec2 );

        // Compute cosine theta of the angle
        cos_theta1 = dot / ( sqrt ( vec1_mag * vec2_mag ) );

        // If cos_theta1 is close to -1, then create the children molecules
        if ( cos_theta1 < -cos( ga_angle_cutoff ) ){

            // Increment overlapping bond variable
            num_overlapping_bonds++;
            // Save to unique vector
            if ( parents[i].used == false ){ 
               xover_parents.push_back(parents[i]); 
               parents[i].used = true;
            }
            if ( parents[k].used == false ){ 
               xover_parents.push_back(parents[k]); 
               parents[k].used = true;
            }
        }
    
    // Clear tmp vectors with manipulated parents from breeding
    tmp_parent1.clear_molecule();
    tmp_parent2.clear_molecule();
    }    
                            }
                        }
                    }
                }
            }
        }
    }
    cout << endl << "#### Crossover Preparation ####" << endl;
    cout << "Unique molecules that overlap: " << xover_parents.size() << endl;
    cout << "Number of overlapping bonds this generation: " << num_overlapping_bonds <<" of " << num_bond_pairs_sampled<< " bonds sampled" <<endl;

     return;
} //end GA_Recomb:check_exhaustive()



// +++++++++++++++++++++++++++++++++++++++++
// Breeding is the unoffical main genetic algorithm function for breeding unique child molecules 
// from two parents - called in max_generations
// NOTE: Needs to be updated
void
GA_Recomb::breeding_exhaustive( Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & c_typer, Orient & orient, int gen_num )
{
    Trace trace( "GA_Recomb::breeding_exhaustive()" );
    // Step 1: Prepare parents for breeding
    // Activate parents before amber_typer
    for ( int i=0; i<parents.size(); i++ ){
        for ( int j=0; j<parents[i].num_bonds; j++ ){
            parents[i].bond_active_flags[j]=true;
        }
    }

    for ( int i=0; i<parents.size(); i++ ){
        for ( int j=0; j<parents[i].num_atoms; j++ ){
            parents[i].atom_active_flags[j]=true;
        }
    }

    // Prepare amber typer for parents
    for ( int i=0; i<parents.size(); i++ ){
        //cout << "Entering Amber typer" << endl;
        parents[i].prepare_molecule();
        c_typer.skip_verbose_flag = true;
        c_typer.prepare_molecule( parents[i], true, score.use_chem, score.use_ph4, score.use_volume );
        parents[i].used = false;
    } 

    // Deactivate all atoms and then bonds in each parent molecule so that later 
    // only the substructures destined to make the children can be activated
    for ( int i=0; i<parents.size(); i++ ){
        for ( int j=0; j<parents[i].num_bonds; j++ ){
            parents[i].bond_active_flags[j]=false;
            parents[i].bond_keep_flags[j]=false;
        }
    }

    for ( int i=0; i<parents.size(); i++ ){
        for ( int j=0; j<parents[i].num_atoms; j++ ){
            parents[i].atom_active_flags[j]=false;
        }
    }

    // Reset the parents comments (energy)
    for ( int i=0; i<parents.size(); i++ ){
        ostringstream new_energy;
        new_energy << "parent";
        parents[i].energy = new_energy.str();
        parents[i].parent = 1;
    }

    // Step 2: Determine the bonds and coorsponding atoms that overlap
    // Step 2A: Determine if the new bond to be formed during breeding is of the appropriate length

    // Set variables that will hold distances between the origin and target atoms
    // of the bonds on both parents that will be involved in breeding
    float dist13,
          dist14,
          dist23,
          dist24;

    // Calculate the distance between bonds on two different parents, 
    // then determine if the bonds are close enough for breeding

    // Check number of overlapping bonds per molecule
    int num_overlapping_bonds = 0;
    int num_bond_pairs_sampled = 0;

    int pair_counter = 0;
    // For each parent
    for ( int i=0; i<parents.size(); i++ ){
        // For each bond j
        for ( int j=0; j<parents[i].num_bonds; j++ ){
             // Identify if the bond is rotatable
             if ( parents[i].amber_bt_id[j] != -1 ){
                // Compare to the rotatable bond in parent i+1
                for ( int k=i+1; k<parents.size(); k++ ){
                    for ( int m=0; m<parents[k].num_bonds; m++ ){
                        if ( parents[k].amber_bt_id[m] != -1 ){
                           if ( parents[i].bond_types[j] == parents[k].bond_types[m] ){
    // Increment the number of overlapping bonds
    num_bond_pairs_sampled++; 

    // Define distance vectors
    DOCKVector vec13;
    vec13.x = parents[i].x[parents[i].bonds_origin_atom[j]] - parents[k].x[parents[k].bonds_origin_atom[m]];
    vec13.y = parents[i].y[parents[i].bonds_origin_atom[j]] - parents[k].y[parents[k].bonds_origin_atom[m]];
    vec13.z = parents[i].z[parents[i].bonds_origin_atom[j]] - parents[k].z[parents[k].bonds_origin_atom[m]];

    DOCKVector vec14;
    vec14.x = parents[i].x[parents[i].bonds_origin_atom[j]] - parents[k].x[parents[k].bonds_target_atom[m]];
    vec14.y = parents[i].y[parents[i].bonds_origin_atom[j]] - parents[k].y[parents[k].bonds_target_atom[m]];
    vec14.z = parents[i].z[parents[i].bonds_origin_atom[j]] - parents[k].z[parents[k].bonds_target_atom[m]];

    DOCKVector vec23;
    vec23.x = parents[i].x[parents[i].bonds_target_atom[j]] - parents[k].x[parents[k].bonds_origin_atom[m]];
    vec23.y = parents[i].y[parents[i].bonds_target_atom[j]] - parents[k].y[parents[k].bonds_origin_atom[m]];
    vec23.z = parents[i].z[parents[i].bonds_target_atom[j]] - parents[k].z[parents[k].bonds_origin_atom[m]];

    DOCKVector vec24;
    vec24.x = parents[i].x[parents[i].bonds_target_atom[j]] - parents[k].x[parents[k].bonds_target_atom[m]];
    vec24.y = parents[i].y[parents[i].bonds_target_atom[j]] - parents[k].y[parents[k].bonds_target_atom[m]];
    vec24.z = parents[i].z[parents[i].bonds_target_atom[j]] - parents[k].z[parents[k].bonds_target_atom[m]];

    // Measure the distance between the origin and target atoms on the rotatable single bonds in the parents
    // ie the magnitude of the vectors
    dist13 = ( vec13.x*vec13.x) +  ( vec13.y*vec13.y) + ( vec13.z*vec13.z);

    dist14 = ( vec14.x*vec14.x) +  ( vec14.y*vec14.y) + ( vec14.z*vec14.z);

    dist23 = ( vec23.x*vec23.x) +  ( vec23.y*vec23.y) + ( vec23.z*vec23.z);
  
    dist24 = ( vec24.x*vec24.x) +  ( vec24.y*vec24.y) + ( vec24.z*vec24.z);

    // Define vectors that coorespond to the idenfied bond in parent
    DOCKVector vec1;    // vector containing the *direction* from origin to 
                        // target atom in parent 1 around identified bond
    DOCKVector vec2;    // vector containing parent 2 bond information
    float vec1_mag;     // vector magnitude 
    float vec2_mag;

    // Variables to determine if the identified bonds will be used for breeding
    int if_bool = 0;         // boolean that identifies the *direction* of the bonds in parents
    float dot;               // dot product value of vectors vec1 & vec2 
    float cos_theta1;        // *angle* between parents' bonds (used to determine the *direction*) 

    // MARK - DELETE
    // Check the distances between the bonds 
    /*cout << "distances " << dist13 << " " << dist14 << " " << dist23 << " " << dist24 << endl;
    cout << "vec13 x " << vec13.x << " y " << vec13.y << " z " << vec13.z << endl;
    cout << "vec14 x " << vec14.x << " y " << vec14.y << " z " << vec14.z << endl;
    cout << "vec23 x " << vec23.x << " y " << vec23.y << " z " << vec23.z << endl;
    cout << "vec24 x " << vec24.x << " y " << vec24.y << " z " << vec24.z << endl;

    cout << " parent " << i << " parent " << k <<endl; 
    cout << " coordinates x o " << parents[i].x[parents[i].bonds_origin_atom[j]] << " " << parents[k].x[parents[k].bonds_origin_atom[m]] << endl 
         << " corrdinates y o " << parents[i].y[parents[i].bonds_origin_atom[j]] << " " << parents[k].y[parents[k].bonds_origin_atom[m]] << endl
         << " coordinates z o " << parents[i].z[parents[i].bonds_origin_atom[j]] << " " << parents[k].z[parents[k].bonds_origin_atom[m]] << endl
         << " coordinates x t " << parents[i].x[parents[i].bonds_target_atom[j]] << " " << parents[k].x[parents[k].bonds_target_atom[m]] << endl
         << " coordinates y t " << parents[i].y[parents[i].bonds_target_atom[j]] << " " << parents[k].y[parents[k].bonds_target_atom[m]] << endl
         << " coordinates z y " << parents[i].z[parents[i].bonds_target_atom[j]] << " " << parents[k].z[parents[k].bonds_target_atom[m]] << endl;

    cout << " parent " << i << " bond " << j << " parent " << k << " bond "<< m <<endl;
    cout <<" dist13: " << dist13 << " dist24: " << dist24 << " dist24: " << dist24 << " dist23: " << dist23 <<endl;*/
    // If all approx distances measured for the new children bonds are within the bond tolerance
    // ie if the two parent bonds are overlapping or nearly overlapped as defined by user defined bond_tolerance
    if ( dist13<ga_bond_tolerance && dist23<ga_bond_tolerance && dist14<ga_bond_tolerance && dist24<ga_bond_tolerance ){

        // Determine which distances provide the best overlap
        if ( (dist13 + dist24) < (dist14 + dist23) ){

            // If the if_bool is true, then the origin-to-origin & target-to-target atom dist is the shortest;
            // Therefore, the *direction* of the origin and target atoms are in the same direction
            if_bool = 1;
        }

        // All other conditions (if (dist14 + dist23) is greater or if the combination of distances are equal) will enter the else condition
        else { if_bool = 2; }
    }
    
    // Else if the origin-to-origin and target-to-target distances are within the tolerance
    else if ( (dist13<ga_bond_tolerance && dist24<ga_bond_tolerance) ){
          if_bool = 1;
    }
    else if ( dist14<ga_bond_tolerance && dist23<ga_bond_tolerance ){
          if_bool = 2;
    }



    //-------------------------------------------------------------------------
    // STEP 2B: Determine if the overlap angle is close to zero and run breeding functions

    // If the if_bool is true, ga_bond_tolerance has been met but the origin & target atoms must be reversed for parent 2
    if ( if_bool == 1 ){
   
        // Print distance cutoff to screen
        /*cout <<endl <<"#### GA distance tolerance was met at: " << dist13 << " and " << dist24<<endl;

        // Check the properties of the bonds that have met necessary criteria
        cout << " parent " << i <<  " " << parents[i].energy << " parent " << k << " " <<  parents[k].energy <<endl; 
        cout << "amber bond id " << parents[i].amber_bt_id[j] <<" " << parents[k].amber_bt_id[m] <<endl;
        cout << "o1 " <<parents[i].bonds_origin_atom[j] <<" " <<parents[i].atom_types[parents[i].bonds_origin_atom[j]] 
             << " t1 " <<parents[i].bonds_target_atom[j] << " " <<parents[i].atom_types[parents[i].bonds_target_atom[j]];
         cout<< " parent " << i << " bond " << j << " o2 "<<parents[k].bonds_origin_atom[m] << " " << parents[k].atom_types[parents[k].bonds_origin_atom[m]] 
             << " t2 " <<parents[k].bonds_target_atom[m] << " " << parents[k].atom_types[parents[k].bonds_target_atom[m]] << " parent " << k << " bond "<< m <<endl;*/
      
         // Populate vectors with the bond length of each parent molecule
         // Note: the vector directions are the same for vec1 and vec2           
         vec1.x = parents[i].x[parents[i].bonds_origin_atom[j]] - parents[i].x[parents[i].bonds_target_atom[j]];  
         vec1.y = parents[i].y[parents[i].bonds_origin_atom[j]] - parents[i].y[parents[i].bonds_target_atom[j]];  
         vec1.z = parents[i].z[parents[i].bonds_origin_atom[j]] - parents[i].z[parents[i].bonds_target_atom[j]];

         vec2.x = parents[k].x[parents[k].bonds_origin_atom[m]] - parents[k].x[parents[k].bonds_target_atom[m]];   
         vec2.y = parents[k].y[parents[k].bonds_origin_atom[m]] - parents[k].y[parents[k].bonds_target_atom[m]];   
         vec2.z = parents[k].z[parents[k].bonds_origin_atom[m]] - parents[k].z[parents[k].bonds_target_atom[m]];
   
         // Magnitude of above vectors
         vec1_mag = ( vec1.x*vec1.x) + ( vec1.y*vec1.y) + ( vec1.z*vec1.z);
         vec2_mag = ( vec2.x*vec2.x) + ( vec2.y*vec2.y) + ( vec2.z*vec2.z);

         // Compute the dot product using the function in utils.cpp
         dot = dot_prod(vec1, vec2);

         // Compute cosine theta
         cos_theta1 = dot / (sqrt (vec1_mag * vec2_mag));

         // If cos_theta1 is close to +1, then the parent bonds are going in the same direction
         if ( cos_theta1 > cos( ga_angle_cutoff ) ){

             // MARK - DELETE: check the cos theta calculated
             //cout << "This bonds passed the overlap criteria with cos theta = " << cos_theta1 <<endl; 
             
             // Increment overlapping bond variable
             num_overlapping_bonds++;

             // Save to unique vector
             if ( ga_check_overlap ){
                if ( parents[i].used == false ){ 
                   xover_parents.push_back(parents[i]); 
                   parents[i].used = true;
                }
                if ( parents[k].used == false ){ 
                   xover_parents.push_back(parents[k]); 
                   parents[k].used = true;
                }
             }

             // Make temp parent molecules for manipulation
             copy_molecule_shallow( tmp_parent1, parents[i] );
             copy_molecule_shallow( tmp_parent2, parents[k] );
             //copy_molecule( new_tmp_parent1, parents[i] );
             copy_molecule_shallow( new_tmp_parent2, parents[k] );

          
             // Determine whether the smallest side of the parent molecules are the same
             bool similar = false;
             similar = similarity_compare( tmp_parent1, tmp_parent1.bonds_origin_atom[i], tmp_parent1.bonds_target_atom[i], 
                                           tmp_parent2, tmp_parent2.bonds_target_atom[j], tmp_parent2.bonds_origin_atom[j],
                                           false, score, simplex, c_typer );
             
             // If the smaller side is not the same, then compare that larger halves                    
             if (similar == false ){ 
                similar = similarity_compare( tmp_parent1, tmp_parent1.bonds_origin_atom[i], tmp_parent1.bonds_target_atom[i], 
                                              tmp_parent2, tmp_parent2.bonds_target_atom[j], tmp_parent2.bonds_origin_atom[j],
                                              true, score, simplex, c_typer );
             }
    
             // If neither side of the molecule is the same
             // TODO: REDUCE REDUNDANT MAPPING
             if (similar == false ){ 
                // For parent mutation, update used flag
                parents[i].used = true; 
                parents[j].used = true; 

                // Activates one half of each parent (the halves opposite one another)
                activate_half_mol( tmp_parent1, tmp_parent1.bonds_origin_atom[j], tmp_parent1.bonds_target_atom[j],
                                   tmp_parent2, tmp_parent2.bonds_target_atom[m], tmp_parent2.bonds_origin_atom[m], false );

                // Creates a new bond between the halves
                recomb_mols_xover( tmp_parent1, tmp_parent1.bonds_origin_atom[j], tmp_parent1.bonds_target_atom[j],
                             tmp_parent2, tmp_parent2.bonds_target_atom[m], tmp_parent2.bonds_origin_atom[m],
                             children );

                // Activates the other half of the parents and called recomb_mols to create child
                switch_active_halves ( tmp_parent1, tmp_parent1.bonds_origin_atom[j], tmp_parent1.bonds_target_atom[j],
                                       new_tmp_parent2, new_tmp_parent2.bonds_target_atom[m], new_tmp_parent2.bonds_origin_atom[m],
                                        children );
             }
        }

       // Clear tmp vectors with parents that have been manipulated during breeding
       tmp_parent1.clear_molecule();
       tmp_parent2.clear_molecule();
       new_tmp_parent2.clear_molecule();
    }

    // Repeat if origin-to-target & target-to-origin atom distances is within ga_bond_tolerance
    else  if ( if_bool == 2 ){

        // Print distance cutoff to screen
        /*cout <<endl <<"#### GA distance tolerance was met at: " << dist14 << " and " << dist23<<endl;

        cout << " parent " << i << " " << parents[i].energy << " parent " << k << " " << parents[k].energy <<endl; 
        cout << "o1 " <<parents[i].bonds_origin_atom[j] <<" " <<parents[i].atom_types[parents[i].bonds_origin_atom[j]] 
             << " t1 " <<parents[i].bonds_target_atom[j] << " " <<parents[i].atom_types[parents[i].bonds_target_atom[j]] 
             << " parent " << i << " bond " << j <<" o2 "<<parents[k].bonds_origin_atom[m] << " " << parents[k].atom_types[parents[k].bonds_origin_atom[m]] 
             << " t2 " <<parents[k].bonds_target_atom[m] << " " << parents[k].atom_types[parents[k].bonds_target_atom[m]] << " parent " << k << " bond "<< m <<endl;

        cout <<"amber bond id " << parents[i].amber_bt_id[j] <<" " << parents[k].amber_bt_id[m] <<endl;*/

        // Populate bond vectors
        vec1.x = parents[i].x[parents[i].bonds_origin_atom[j]] - parents[i].x[parents[i].bonds_target_atom[j]];  
        vec1.y = parents[i].y[parents[i].bonds_origin_atom[j]] - parents[i].y[parents[i].bonds_target_atom[j]];  
        vec1.z = parents[i].z[parents[i].bonds_origin_atom[j]] - parents[i].z[parents[i].bonds_target_atom[j]];

        vec2.x = parents[k].x[parents[k].bonds_origin_atom[m]] - parents[k].x[parents[k].bonds_target_atom[m]];   
        vec2.y = parents[k].y[parents[k].bonds_origin_atom[m]] - parents[k].y[parents[k].bonds_target_atom[m]];   
        vec2.z = parents[k].z[parents[k].bonds_origin_atom[m]] - parents[k].z[parents[k].bonds_target_atom[m]];
   
        // Calculate magnitude of above vectors
        vec1_mag = ( vec1.x*vec1.x) + ( vec1.y*vec1.y) + ( vec1.z*vec1.z);
        vec2_mag = ( vec2.x*vec2.x) + ( vec2.y*vec2.y) + ( vec2.z*vec2.z);

        // Compute the dot product using the function in utils.cpp
        dot = dot_prod( vec1, vec2 );

        // Compute cosine theta of the angle
        cos_theta1 = dot / ( sqrt ( vec1_mag * vec2_mag ) );

        // If cos_theta1 is close to -1, then create the children molecules
        if ( cos_theta1 < -cos( ga_angle_cutoff ) ){

            // Increment overlapping bond variable
            num_overlapping_bonds++;

            // Save to unique vector
            if ( ga_check_overlap ){
               if ( parents[i].used == false ){ 
                  xover_parents.push_back(parents[i]); 
                  parents[i].used = true;
               }
               if ( parents[k].used == false ){ 
                  xover_parents.push_back(parents[k]); 
                  parents[k].used = true;
               }
            }

            // Make temp parent molecules for manipulation
            copy_molecule_shallow( tmp_parent1, parents[i] );
            copy_molecule_shallow( tmp_parent2, parents[k] );
            copy_molecule_shallow( new_tmp_parent2, parents[k] );
          
             // Detemine whether the smallest side of the parent molecules are the same
             bool similar = false;
             similar = similarity_compare( tmp_parent1, tmp_parent1.bonds_origin_atom[i], tmp_parent1.bonds_target_atom[i], 
                                           tmp_parent2, tmp_parent2.bonds_origin_atom[j], tmp_parent2.bonds_target_atom[j],
                                           false, score, simplex, c_typer );

             // If the smaller side is not the same, then compare that larger halves                    
             if (similar == false ){ 
                similar = similarity_compare( tmp_parent1, tmp_parent1.bonds_origin_atom[i], tmp_parent1.bonds_target_atom[i], 
                                              tmp_parent2, tmp_parent2.bonds_origin_atom[j], tmp_parent2.bonds_target_atom[j],
                                              true, score, simplex, c_typer );
             }

             // If neither side of the molecule is the same
             // TODO: REDUCE REDUNDANT MAPPING
             if (similar == false){
                // For parent mutation, update used flag
                parents[i].used = true; 
                parents[j].used = true; 

                activate_half_mol( tmp_parent1, tmp_parent1.bonds_origin_atom[j], tmp_parent1.bonds_target_atom[j],
                                   tmp_parent2, tmp_parent2.bonds_origin_atom[m], tmp_parent2.bonds_target_atom[m], false );

                recomb_mols_xover( tmp_parent1, tmp_parent1.bonds_origin_atom[j], tmp_parent1.bonds_target_atom[j],
                             tmp_parent2, tmp_parent2.bonds_origin_atom[m], tmp_parent2.bonds_target_atom[m],
                             children );

                switch_active_halves ( tmp_parent1, tmp_parent1.bonds_origin_atom[j], tmp_parent1.bonds_target_atom[j],
                                       new_tmp_parent2, new_tmp_parent2.bonds_origin_atom[m], new_tmp_parent2.bonds_target_atom[m],
                                       children );
            }
        }
    
    // Clear tmp vectors with manipulated parents from breeding
    tmp_parent1.clear_molecule();
    tmp_parent2.clear_molecule();
    //new_tmp_parent1.clear_molecule();
    new_tmp_parent2.clear_molecule();
    }    
                            }
                        }
                    }
                }
            }
        }
    }

    cout <<endl << "#### In breeding, size of children " << children.size() << endl;
    cout << "Number of overlapping bonds this generation " << num_overlapping_bonds <<" of " << num_bond_pairs_sampled<< " bonds sampled" <<endl;
    cout << "Number of parents that can crossover " << xover_parents.size() <<endl;

    //-------------------------------------------------------------------------
    // STEP 3: Minimize and Pruning

    // Activate parents and children for pruning
    activate_vector( parents );
    activate_vector( children );

    // Minimize children with valid torenv
    int num_valid = 0;
    int num_invalid = children.size();
    
    for ( int i=0; i<children.size(); i++){
        // If using torenv table
        if ( ga_use_torenv_table ){
            // Prepare molecule for fingerprint/bond environ generation
            prepare_mol_torenv( children[i] );
            //In case torenv depends on atom types?
            //cout << "Entering Amber typer" << endl;
            c_typer.skip_verbose_flag = true;
            c_typer.prepare_molecule( children[i], true, score.use_chem, score.use_ph4, score.use_volume );

            DN_GA_Build c_dn;
            Fragment tmp_child = mol_to_frag(children[i]);
            c_dn.dn_ga_flag = true;
            c_dn.dn_use_roulette = ga_use_dn_roulette;
            c_dn.dn_use_torenv_table = ga_use_torenv_table;
            c_dn.dn_torenv_table = ga_torenv_table; 
            c_dn.read_torenv_table(c_dn.dn_torenv_table);
            c_dn.verbose = 0;
            c_dn.dn_MW_cutoff_type_hard = true;
            c_dn.dn_MW_cutoff_type_soft = false;

            if ( c_dn.valid_torenv( tmp_child ) ){
	       // Check all rot bonds 
	       prepare_torenv_indices(tmp_child);
               if ( c_dn.valid_torenv_multi( tmp_child ) ){
                  // Minimize and increment counter
                  num_valid++;
                  minimize_children( children[i], pruned_children, score, simplex, c_typer );
 
                  //Clear vectors for torenv
                  //c_dn.torenv_vector.target_envs.clear();
                  //c_dn.torenv_vector.target_freqs.clear();
                  c_dn.torenv_vector.clear();
                  tmp_child.torenv_recheck_indices.clear();
              }
           }
        }
        else{
            minimize_children( children[i], pruned_children, score, simplex, c_typer );
                     

        }
    }

    // REMOVE    
    /*fstream fout_molecules;
    fout_molecules.open ( "zzz.minimized_children.mol2", fstream::out|fstream::app );
    for( int i=0; i<pruned_children.size(); i++ ){
        fout_molecules << pruned_children[i].current_data << endl;
        Write_Mol2( pruned_children[i], fout_molecules );
    }
    fout_molecules.close();*/

    //Clear children      
    children.clear();
    cout << "The number of children with valid tornev are " << num_valid <<endl;
    cout << "The number of children deleted: " << num_invalid - num_valid  <<endl; 

    // This function will select the best scoring conformer and use in the next generation
    fitness_pruning( pruned_children, scored_generation, score, simplex, c_typer );
    cout << "The number of scored children are " << scored_generation.size() <<endl;

    //Check if size is 0 then exit
    if (scored_generation.size() == 0){
       cout << "WARNING: THERE ARE NO OFFSPRING (from pruning and crossover) AND THE PROGRAM WILL EXIT)"<<endl;
       if (ga_mutate_parents == false){
          cout << "ERROR: Parent mutation needed" << endl;
          exit(0);
       }
    }

    // MARK - Modify function, check that it is working and implement
    if (scored_generation.size()>0){
        hard_filter(scored_generation); 
    }
    // Clear vectors to free up memory
    pruned_children.clear(); 
 

    //-------------------------------------------------------------------------
    // STEP 4: Mutations
    // Determine which parents did not contribute to the final ensemble

    // Reset used
    for (int i=0; i<parents.size(); i++){
        parents[i].used = false;
    }

    // Search through the comments of each molecule
    for ( int i=0; i<scored_generation.size(); i++ ){
        int p1 = 0;
        int p2 = 0;
        istringstream tmp(scored_generation[i].energy);
        tmp>> p1 >> p2;
 
        // For parent mutation, update used flag
        parents[p1].used = true; 
        parents[p2].used = true; 

        // Save parent index information in energy
        ostringstream new_comment;
        new_comment << p1 << "_" << p2;
        scored_generation[i].energy = new_comment.str();
        new_comment.clear();
    }

    // Print the number of unique parents involved in offspring crossover  
    int xover_parents = 0;
    for (int i=0; i<parents.size(); i++){
        // If the molecule was not used for crossover
        if ( parents[i].used == true ){
           xover_parents++;
        }
    }
    //cout <<"Total number of parents involved in xover: " << xover_parents <<endl;
     
    // STEP 4A: If mutating parents
    if ( ga_mutations ){
        if (ga_mutate_parents == true ){
            std::vector <DOCKMol> pmuts;
            std::vector <DOCKMol> ppruned; 
            std::vector <DOCKMol> pfinal; 
      
            // Save used status 
            for (int i=0; i<parents.size(); i++){
                // If the molecule was not used for crossover
                if ( parents[i].used == false ){
                pmuts.push_back( parents[i] );
                }
            }
            master_mut_exhaustive( pmuts, saved_parents, ppruned, true, score, simplex, c_typer, orient, gen_num );

            pmuts.clear();

            // If mutations were successful
            if (ppruned.size() !=0){
                // After scoring, prune the lot
                fitness_pruning(ppruned, pfinal, score, simplex, c_typer);
                //cout << "pfinal, parents" << endl;
                uniqueness_prune_mut(pfinal, parents, score, c_typer);

                // Combine final mutated parents with original 
                for (int i=0; i<pfinal.size(); i++ ){
                    mutated_parents.push_back(pfinal[i]);
                }
 
                ppruned.clear();
                pfinal.clear();
            }
            saved_parents.clear();

        }
    
    
        // STEP 4B: Mutate offspring
        if (scored_generation.size() != 0){
            master_mut_exhaustive( scored_generation, saved_parents, pruned_children, false, score, simplex, c_typer, orient, gen_num );
            scored_generation.clear();

            // For Mutants: After scoring in mutants function, prune the lot
            fitness_pruning(pruned_children, scored_generation, score, simplex, c_typer);

            // Combine with offspring from crossover
            // Run hard Filters on Crossover offspring
            if (saved_parents.size()>0){
                hard_filter(saved_parents); 
            }
            for (int i=0; i<saved_parents.size(); i++ ){
                scored_generation.push_back(saved_parents[i]);
            }

            saved_parents.clear();
            pruned_children.clear();

            // Prune again - with all offspring
            for (int i=0; i<scored_generation.size(); i++){
                   pruned_children.push_back(scored_generation[i]);
            } 
            scored_generation.clear();

            fitness_pruning(pruned_children, scored_generation, score, simplex, c_typer);
            pruned_children.clear();


            // If mutating the parents
            if ((ga_mutate_parents == true) && (mutated_parents.size() !=0)){
                  //cout << "scoredgen, mutatedparents" << endl;   
                  uniqueness_prune_mut(scored_generation, mutated_parents, score, c_typer);
            }

           // Prune offspring that are similar to the parents
            //cout << "scoredgen, parents" << endl;
            uniqueness_prune_mut(scored_generation, parents, score, c_typer);

        }
        // Save mutated parents to offspring vector
        if ((ga_mutate_parents == true) && (mutated_parents.size() != 0)){
            for ( int i=0; i<mutated_parents.size(); i++){
                scored_generation.push_back(mutated_parents[i]);
            }
        }
        // Clear necessary vectors
        mutated_parents.clear();

        return;
    }// end if ga mutations on loop
} // end GA_Recomb::breeding_exhaustive()


 
// +++++++++++++++++++++++++++++++++++++++++
// Breeding is the unoffical main genetic algorithm function for breeding unique child molecules 
// from two parents - called in max_generations
void
GA_Recomb::breeding_rand( Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & c_typer, Orient & orient, int gen_num )
{

    Trace trace( "GA_Recomb::breeding_rand()" );
    // Step 1: Prepare parents for breeding
    // Activate parents before amber_typer
    for ( int i=0; i<parents.size(); i++ ){
        for ( int j=0; j<parents[i].num_bonds; j++ ){
            parents[i].bond_active_flags[j]=true;
            parents[i].bond_keep_flags[j]=false;
        }
    }

    for ( int i=0; i<parents.size(); i++ ){
        for ( int j=0; j<parents[i].num_atoms; j++ ){
            parents[i].atom_active_flags[j]=true;
        }
    }

    // Prepare amber typer for parents
    for ( int i=0; i<parents.size(); i++ ){
        //cout << "Entering Amber typer" << endl;
        parents[i].prepare_molecule();
        c_typer.skip_verbose_flag = true;
        c_typer.prepare_molecule( parents[i], true, score.use_chem, score.use_ph4, score.use_volume );
    } 

    // Deactivate all atoms and then bonds in each parent molecule so that later 
    // only the substructures destined to make the children can be activated
    for ( int i=0; i<parents.size(); i++ ){
        for ( int j=0; j<parents[i].num_bonds; j++ ){
            parents[i].bond_active_flags[j]=false;
            parents[i].bond_keep_flags[j]=false;
        }
    }

    for ( int i=0; i<parents.size(); i++ ){
        for ( int j=0; j<parents[i].num_atoms; j++ ){
            parents[i].atom_active_flags[j]=false;
        }
    }

    // Reset the parents comments (energy)
    for ( int i=0; i<parents.size(); i++ ){
        ostringstream new_energy;
        new_energy << "parent";
        parents[i].energy    = new_energy.str();
    }

    // Step 2: Select parents for breeding
    if ( ga_xover_on ){

        // Step 2A: Set search parameters
        // Determine the max number of pairs if exhaustive sampling, but without redundant pairs
        int pair_max =  (parents.size() * (parents.size() - 1))/2;
    
        // Set counter for number of pairs
        int pair = 0;

        // Define a matrix to determine whether a pair has already been used
        int ppairs[parents.size()][parents.size()];
        // Intialize the matrix
        for (int i=0; i<parents.size(); i++){
            for (int j=0; j<parents.size(); j++){
                ppairs[i][j] = 0;
            }
        }

        //Should add verbose flag
        double start_time_step2 = time_seconds();


        // Step 2B: Infastructure for 2 random parents
        // If you have not exceeded max number of pairs or offspring
        // Total number of mol sampled
        int total_num_mols_made = 0;
        int valid_torenv = 0;
        int valid_descriptors = 0;
        int hungarian_similarity_prune = 0;

    	while ((pair < pair_max) & (children.size() < ga_xover_max)){
    
       	    // Select two parents for breeding randomly
       	    int index_p1 = (rand() % parents.size()); 
       	    int index_p2 = (rand() % parents.size());

            // Check to make sure that two parents are not the same
            if (index_p1 != index_p2){
              // If the parent pair has not been used
              if (ppairs[index_p1][index_p2] == 0){

                 // Increment pair counter
                 pair++;
                 // Update ppairs so that p1/p2 and p2/p1 cannot be pairs
                 ppairs[index_p1][index_p2] = 1;
                 ppairs[index_p2][index_p1] = 1;

                 // Exit if we have exceeded pair_max
                 if (pair > pair_max) {break;}

                 // Complete exaustive crossover
                 // Set variables that will hold distances between the origin and target atoms
                 // of the bonds on both parents that will be involved in breeding
                 float dist13,
                       dist14,
                       dist23,
                       dist24;

                 // Calculate the distance between bonds on two different parents, 
                 // then determine if the bonds are close enough for breeding

                 // Check number of overlapping bonds per molecule
                 int num_overlapping_bonds = 0;
                 int num_bond_pairs_sampled = 0;
                 // For each bond in parent1 
                 for ( int i=0; i<parents[index_p1].num_bonds; i++ ){
                     // Identify if the bond is rotatable
                     if ( parents[index_p1].amber_bt_id[i] != -1 ){
                        // Compare to the rotatable bond in parent 2
                        for ( int j=0; j<parents[index_p2].num_bonds; j++ ){
                            if ( parents[index_p2].amber_bt_id[j] != -1 ){
                               if ( parents[index_p1].bond_types[i] == parents[index_p2].bond_types[j] ){
                               
                                  // Exit if the children size exceeds cutoff
                                  if ( (children.size()+ tmp_children.size()) >= ga_xover_max) {break;}
       
                                  // Increment the number of overlapping bonds
                                  num_bond_pairs_sampled++;

                                  // Define distance vectors
                                  DOCKVector vec13;
                                  vec13.x = parents[index_p1].x[parents[index_p1].bonds_origin_atom[i]] 
                                          - parents[index_p2].x[parents[index_p2].bonds_origin_atom[j]];
                                  vec13.y = parents[index_p1].y[parents[index_p1].bonds_origin_atom[i]] 
                                          - parents[index_p2].y[parents[index_p2].bonds_origin_atom[j]];
                                  vec13.z = parents[index_p1].z[parents[index_p1].bonds_origin_atom[i]] 
                                          - parents[index_p2].z[parents[index_p2].bonds_origin_atom[j]];

                                  DOCKVector vec14;
                                  vec14.x = parents[index_p1].x[parents[index_p1].bonds_origin_atom[i]] 
                                          - parents[index_p2].x[parents[index_p2].bonds_target_atom[j]];
                                  vec14.y = parents[index_p1].y[parents[index_p1].bonds_origin_atom[i]] 
                                          - parents[index_p2].y[parents[index_p2].bonds_target_atom[j]];
                                  vec14.z = parents[index_p1].z[parents[index_p1].bonds_origin_atom[i]] 
                                          - parents[index_p2].z[parents[index_p2].bonds_target_atom[j]];

                                  DOCKVector vec23;
                                  vec23.x = parents[index_p1].x[parents[index_p1].bonds_target_atom[i]] 
                                          - parents[index_p2].x[parents[index_p2].bonds_origin_atom[j]];
                                  vec23.y = parents[index_p1].y[parents[index_p1].bonds_target_atom[i]] 
                                          - parents[index_p2].y[parents[index_p2].bonds_origin_atom[j]];
                                  vec23.z = parents[index_p1].z[parents[index_p1].bonds_target_atom[i]] 
                                          - parents[index_p2].z[parents[index_p2].bonds_origin_atom[j]];


                                  DOCKVector vec24;
                                  vec24.x = parents[index_p1].x[parents[index_p1].bonds_target_atom[i]] 
                                          - parents[index_p2].x[parents[index_p2].bonds_target_atom[j]];
                                  vec24.y = parents[index_p1].y[parents[index_p1].bonds_target_atom[i]] 
                                          - parents[index_p2].y[parents[index_p2].bonds_target_atom[j]];
                                  vec24.z = parents[index_p1].z[parents[index_p1].bonds_target_atom[i]] 
                                          - parents[index_p2].z[parents[index_p2].bonds_target_atom[j]];

                                  // Measure the distance between the origin and target atoms on the rotatable 
                                  // single bonds in the parents ie the magnitude of the vectors
                                  dist13 = ( vec13.x*vec13.x) +  ( vec13.y*vec13.y) + ( vec13.z*vec13.z);
                                  dist14 = ( vec14.x*vec14.x) +  ( vec14.y*vec14.y) + ( vec14.z*vec14.z);
                                  dist23 = ( vec23.x*vec23.x) +  ( vec23.y*vec23.y) + ( vec23.z*vec23.z);
                                  dist24 = ( vec24.x*vec24.x) +  ( vec24.y*vec24.y) + ( vec24.z*vec24.z);

                                  // Define vectors that coorespond to the idenfied bond in parent
                                  DOCKVector vec1;    // vector containing the *direction* from origin to 
                                                      // target atom in parent 1 around identified bond
                                  DOCKVector vec2;    // vector containing parent 2 bond information
                                  float vec1_mag;     // vector magnitude 
                                  float vec2_mag;

                                  // Variables to determine if the identified bonds will be used for breeding
                                  int if_bool = 0;         // boolean that identifies the *direction* of the bonds
                                  float dot;               // dot product value of vectors vec1 & vec2 
                                  float cos_theta1;        // *angle* between parents' bonds (used to determine 
                                                         //the *direction*) 


                              //cout << " parent " << index_p1 << " bond " << i << " parent " << index_p2 << " bond "<< j <<endl;
                              //cout <<" dist13: " << dist13 << " dist24: " << dist24 << " dist24: " << dist24 << " dist23: " << dist23 <<endl;
                              // If all approx distances measured for the new children bonds are within the 
                              // bond tolerance ie if the two parent bonds are overlapping or nearly overlapped
                              // as defined by user defined bond_tolerance
                                  if ( dist13<ga_bond_tolerance && dist23<ga_bond_tolerance 
                                    && dist14<ga_bond_tolerance && dist24<ga_bond_tolerance ){

                                     // Determine which distances provide the best overlap
                                     if ( (dist13 + dist24) < (dist14 + dist23) ){

                                     // If the if_bool is true, then the origin-to-origin & target-to-target 
                                     // atom dist is the shortest; Therefore, the *direction* of the origin 
                                     // and target atoms are in the same direction
                                         if_bool = 1;
                                     }

                                     // All other conditions (if (dist14 + dist23) is greater or if the 
				     // combination of distances are equal) will enter the else condition
                                     else { if_bool = 2; }
                                 }
    
                                 // Else if the origin-to-origin and target-to-target distances are within the tolerance
                                 else if ( (dist13<ga_bond_tolerance && dist24<ga_bond_tolerance) ){
                                       if_bool = 1;
                                 }
                                 else if ( dist14<ga_bond_tolerance && dist23<ga_bond_tolerance ){
                                       if_bool = 2;
                                 }



                             // STEP 2B: Determine if the overlap angle is close to zero and run breeding 

                             // If the if_bool is true, ga_bond_tolerance has been met but the origin & target atoms must be reversed for parent 2
                                 if ( if_bool == 1 ){
                                //cout << "Entering if_bool=1" << endl;
                                // Print distance cutoff to screen
                                // DELETE - REMOVE
                                if (verbose) cout <<endl <<"#### if GA distance tolerance was met at: " << dist13 << " and " << dist24<<endl;

                                // Check the properties of the bonds that have met necessary criteria
                                //cout << " LAUREN: parent " << index_p1 <<  " " << parents[index_p1].title << " parent " << index_p2 << " " <<  parents[index_p2].title <<endl; 
      
                                    // Populate vectors with the bond length of each parent molecule
                                    // Note: the vector directions are the same for vec1 and vec2           
                                    vec1.x = parents[index_p1].x[parents[index_p1].bonds_origin_atom[i]] 
                                           - parents[index_p1].x[parents[index_p1].bonds_target_atom[i]];  
                                    vec1.y = parents[index_p1].y[parents[index_p1].bonds_origin_atom[i]] 
                                           - parents[index_p1].y[parents[index_p1].bonds_target_atom[i]];  
                                    vec1.z = parents[index_p1].z[parents[index_p1].bonds_origin_atom[i]] 
                                           - parents[index_p1].z[parents[index_p1].bonds_target_atom[i]];
 
                                    vec2.x = parents[index_p2].x[parents[index_p2].bonds_origin_atom[j]] 
                                           - parents[index_p2].x[parents[index_p2].bonds_target_atom[j]];  
                                    vec2.y = parents[index_p2].y[parents[index_p2].bonds_origin_atom[j]] 
                                           - parents[index_p2].y[parents[index_p2].bonds_target_atom[j]];  
                                    vec2.z = parents[index_p2].z[parents[index_p2].bonds_origin_atom[j]] 
                                           - parents[index_p2].z[parents[index_p2].bonds_target_atom[j]];
   
                                    // Magnitude of above vectors
                                    vec1_mag = ( vec1.x*vec1.x) + ( vec1.y*vec1.y) + ( vec1.z*vec1.z);
                                    vec2_mag = ( vec2.x*vec2.x) + ( vec2.y*vec2.y) + ( vec2.z*vec2.z);
 
                                    // Compute the dot product using the function in utils.cpp
                                    dot = dot_prod(vec1, vec2);

                                    // Compute cosine theta
                                    cos_theta1 = dot / (sqrt (vec1_mag * vec2_mag));

                                    // If cos_theta1 is close to +1, then the parent 
				                    // bonds are going in the same direction
                                    if ( cos_theta1 > cos( ga_angle_cutoff ) ){
 
                                       // Increment overlapping bond variable
                                       num_overlapping_bonds++;
                                       // Make temp parent molecules for manipulation
                                       copy_molecule_shallow( tmp_parent1, parents[index_p1] );
                                       copy_molecule_shallow( tmp_parent2, parents[index_p2] );
                                       // LEP - double copy to remove rotation error
                                       //copy_molecule( new_tmp_parent1, parents[index_p1] );
                                       copy_molecule_shallow( new_tmp_parent2, parents[index_p2] );
                                       // Detemine whether the smallest side of the parent molecules are the same
                                       bool similar = false;
                                       similar = similarity_compare( tmp_parent1, tmp_parent1.bonds_origin_atom[i], tmp_parent1.bonds_target_atom[i], 
                                                                     tmp_parent2, tmp_parent2.bonds_target_atom[j], tmp_parent2.bonds_origin_atom[j],
                                                                     false, score, simplex, c_typer );
                
                                       if (similar == true){ hungarian_similarity_prune++; }

                                       // If the smaller side is not the same, then compare that larger halves
                                       else{ 
                                          similar = similarity_compare( tmp_parent1, tmp_parent1.bonds_origin_atom[i], tmp_parent1.bonds_target_atom[i], 
                                                                        tmp_parent2, tmp_parent2.bonds_target_atom[j], tmp_parent2.bonds_origin_atom[j],
                                                                        true, score, simplex, c_typer );
               
                                       }

                                       // If neither side of the molecule is the same
                                       // TODO: REDUCE REDUNDANT MAPPING
                                       if (similar == true){ hungarian_similarity_prune++; }
                                       else{ 
                                        //cout << "Entering activate" << endl;
                                          // Activates one half of each parent (the halves opposite one another)
                                          activate_half_mol( tmp_parent1, tmp_parent1.bonds_origin_atom[i], tmp_parent1.bonds_target_atom[i],
                                                             tmp_parent2, tmp_parent2.bonds_target_atom[j], tmp_parent2.bonds_origin_atom[j], false );
                                       //cout << "Entering recomb_mols" << endl;         
                                            //cout << "children size before recomb" << tmp_children.size() << endl;
                                          // Creates a new bond between the halves
                                          recomb_mols_xover( tmp_parent1, tmp_parent1.bonds_origin_atom[i], tmp_parent1.bonds_target_atom[i],
                                                       tmp_parent2, tmp_parent2.bonds_target_atom[j], tmp_parent2.bonds_origin_atom[j],
                                                       tmp_children );
                                            //cout <<"Entering switch" << endl;
                                            //cout << "children size after recomb" << tmp_children.size() << endl;
                                          // Activates the other half of the parents and called recomb_mols to create child
                                          switch_active_halves ( tmp_parent1, tmp_parent1.bonds_origin_atom[i], tmp_parent1.bonds_target_atom[i],
                                                                 new_tmp_parent2, new_tmp_parent2.bonds_target_atom[j], new_tmp_parent2.bonds_origin_atom[j],
                                                                 tmp_children );
                                            //cout << "children size after switch" << tmp_children.size()<< endl;
                                       }
                                            
                                    }

                                    // Clear tmp vectors with parents that have been manipulated during breeding
                                    tmp_parent1.clear_molecule();
                                    tmp_parent2.clear_molecule();
                                    //LEP - remove double copies
                                    //new_tmp_parent1.clear_molecule();
                                    new_tmp_parent2.clear_molecule();

                                 }// end if bool
   

                                 // Repeat if origin-to-target & target-to-origin atom distances is within ga_bond_tolerance
                                 else if ( if_bool == 2 ){
                                //cout << "Entering if_bool=2" << endl;
            // Print distance cutoff to screen
            if (verbose) cout <<endl <<"#### else if GA distance tolerance was met at: " << dist14 << " and " << dist23<<endl;

                                    // Populate bond vectors
                                    vec1.x = parents[index_p1].x[parents[index_p1].bonds_origin_atom[i]] 
                                           - parents[index_p1].x[parents[index_p1].bonds_target_atom[i]];  
                                    vec1.y = parents[index_p1].y[parents[index_p1].bonds_origin_atom[i]] 
                                           - parents[index_p1].y[parents[index_p1].bonds_target_atom[i]];  
                                    vec1.z = parents[index_p1].z[parents[index_p1].bonds_origin_atom[i]] 
                                           - parents[index_p1].z[parents[index_p1].bonds_target_atom[i]];

                                    vec2.x = parents[index_p2].x[parents[index_p2].bonds_origin_atom[j]] 
                                           - parents[index_p2].x[parents[index_p2].bonds_target_atom[j]];   
                                    vec2.y = parents[index_p2].y[parents[index_p2].bonds_origin_atom[j]] 
                                           - parents[index_p2].y[parents[index_p2].bonds_target_atom[j]];   
                                    vec2.z = parents[index_p2].z[parents[index_p2].bonds_origin_atom[j]] 
                                           - parents[index_p2].z[parents[index_p2].bonds_target_atom[j]];
   
                                // Calculate magnitude of above vectors
                                    vec1_mag = ( vec1.x*vec1.x) + ( vec1.y*vec1.y) + ( vec1.z*vec1.z);
                                    vec2_mag = ( vec2.x*vec2.x) + ( vec2.y*vec2.y) + ( vec2.z*vec2.z);

                                    // Compute the dot product using the function in utils.cpp
                                    dot = dot_prod( vec1, vec2 );
 
                                    // Compute cosine theta of the angle
                                    cos_theta1 = dot / ( sqrt ( vec1_mag * vec2_mag ) );

                                    // If cos_theta1 is close to -1, then create the children molecules
                                    if ( cos_theta1 < -cos( ga_angle_cutoff ) ){

                                       // Increment overlapping bond variable
                                       num_overlapping_bonds++;
 
                                       // Make temp parent molecules for manipulation
                                       copy_molecule_shallow( tmp_parent1, parents[index_p1] );
                                       copy_molecule_shallow( tmp_parent2, parents[index_p2] );
                                       //LEP - make double copies 
                                       //copy_molecule( new_tmp_parent1, parents[index_p1] );
                                       copy_molecule_shallow( new_tmp_parent2, parents[index_p2] );
                                       // Detemine whether the smallest side of the parent molecules are the same
                                       bool similar = false;
                                       similar = similarity_compare( tmp_parent1, tmp_parent1.bonds_origin_atom[i], tmp_parent1.bonds_target_atom[i], 
                                                                     tmp_parent2, tmp_parent2.bonds_origin_atom[j], tmp_parent2.bonds_target_atom[j],
                                                                     false, score, simplex, c_typer );
  
                                       if (similar == true){ hungarian_similarity_prune++; }

                                       // If the smaller side is not the same, then compare that larger halves
                                       if (similar == false ){ 
                                          similar = similarity_compare( tmp_parent1, tmp_parent1.bonds_origin_atom[i], tmp_parent1.bonds_target_atom[i], 
                                                                        tmp_parent2, tmp_parent2.bonds_origin_atom[j], tmp_parent2.bonds_target_atom[j],
                                                                        true, score, simplex, c_typer );
             
                                       }

                                       // If neither side of the molecule is the same
                                       // TODO: REDUCE REDUNDANT MAPPING
                                       if (similar == true){ hungarian_similarity_prune++; }
                                       else{
                                          //cout << "Entering activate" << endl;
                                          activate_half_mol( tmp_parent1, tmp_parent1.bonds_origin_atom[i], tmp_parent1.bonds_target_atom[i],
                                                              tmp_parent2, tmp_parent2.bonds_origin_atom[j], tmp_parent2.bonds_target_atom[j], false );
                                        //cout << "Entering recob_mols" << endl;
                                        //cout << "children size before recomb" << tmp_children.size() << endl;
                                          recomb_mols_xover( tmp_parent1, tmp_parent1.bonds_origin_atom[i], tmp_parent1.bonds_target_atom[i],
                                                       tmp_parent2, tmp_parent2.bonds_origin_atom[j], tmp_parent2.bonds_target_atom[j],
                                                       tmp_children );
                                        //cout << "Entering switch" << endl;
                                        //cout << "children size after recomb" << tmp_children.size() << endl;
                                          switch_active_halves ( tmp_parent1, tmp_parent1.bonds_origin_atom[i], tmp_parent1.bonds_target_atom[i],
                                                                 new_tmp_parent2, new_tmp_parent2.bonds_origin_atom[j], new_tmp_parent2.bonds_target_atom[j],
                                                                 tmp_children );
                                        //cout << "children size after switch" << tmp_children.size() << endl;
                                       }
                                    }
                                    
                                    // Clear tmp vectors with manipulated parents from breeding
                                    tmp_parent1.clear_molecule();
                                    tmp_parent2.clear_molecule();
                                    //LEP - clear double copie
                                    //new_tmp_parent1.clear_molecule();
                                    new_tmp_parent2.clear_molecule();

                                 }// end else if
                              }//end if same bond order
                           }// end if rot
                        }// end for p1 bonds
                     }// end if rot
                  }// end for p1 bonds
                 //Check torenv
                 //cout << "TOTAL KIDS FROM XOVER PRE-TORENV:"<< tmp_children.size() << endl;
                 int num_valid = 0;
                 //cout << "checking torenvs..." << endl;
                 for ( int i=0; i<tmp_children.size(); i++){
                     // If using torenv table
                     if ( ga_use_torenv_table ){
                        // Prepare molecule for fingerprint/bond environ generation
                        prepare_mol_torenv( tmp_children[i] );
                        //In case torenv depends on atom types?
                        //cout << "Entering Amber typer" << endl;
                        c_typer.skip_verbose_flag = true;
                        c_typer.prepare_molecule( tmp_children[i], true, score.use_chem, score.use_ph4, score.use_volume );

                        DN_GA_Build c_dn;
                        Fragment tmp_child = mol_to_frag(tmp_children[i]);
                        c_dn.dn_ga_flag = true;
                        c_dn.dn_use_roulette = ga_use_dn_roulette;
                        c_dn.dn_use_torenv_table = ga_use_torenv_table;
                        c_dn.dn_torenv_table = ga_torenv_table; 
                        c_dn.read_torenv_table(c_dn.dn_torenv_table);
                        c_dn.verbose = 0;
                        c_dn.dn_MW_cutoff_type_hard = true;
                        c_dn.dn_MW_cutoff_type_soft = false;

                        // Increment for number of molecules made            
                        total_num_mols_made++;
    
                        if ( c_dn.valid_torenv( tmp_child ) ){
    	                   // Check all rot bonds 
	                       prepare_torenv_indices(tmp_child);
                           if ( c_dn.valid_torenv_multi( tmp_child ) ){
                              //Push back to real children vector
                              num_valid++;
                              valid_torenv++;

                              // Activate all atoms and bonds prior to any scoring
                              activate_mol( tmp_children[i] );

                              // Declare a temporary fingerprint object and compute atom environments 
                              // (necessary for Gasteiger)
                              Fingerprint temp_finger;
                              for ( int j=0; j<tmp_children[i].num_atoms; j++ ){
                                  tmp_children[i].atom_envs[j] = temp_finger.return_environment( tmp_children[i], j );
                              } 

                              // Prepare the molecule using the amber_typer to assign bond types to each bond
                              //cout << "Entering Amber typer" << endl;
                              tmp_children[i].prepare_molecule();
                              c_typer.skip_verbose_flag = true;
                              c_typer.prepare_molecule( tmp_children[i], true, score.use_chem, score.use_ph4, score.use_volume );

                              // Check to see if the molecule exceeds the cutoffs (if not true, true)
                              if (hard_filter_mol(tmp_children[i]) == true){

                                 // Make energy the index of the parent molecules
                                 stringstream new_energy;
                                 new_energy << index_p1 <<" "<< index_p2;
                                 tmp_children[i].energy = new_energy.str();
                                 new_energy.clear();

                                 children.push_back(tmp_children[i]);
                                 valid_descriptors++;
                              }
                           }
                        }
                        //Clear vectors for torenv
                        c_dn.torenv_vector.clear();
                        tmp_child.torenv_recheck_indices.clear();
                     }
                     else{ children.push_back(tmp_children[i]);}
                           //good_molecules_go_here.push_back(tmp_children[i]); }
                  } 
                  // Delete molecules from tmp_children
                  tmp_children.clear();
                  //cout << "TOTAL KIDS FROM XOVER POST-TORENV:"<< children.size() << endl;
               } //end if ppairs
            } //end if index     
        }// End while crossover pair selection

    // Clear ppairs
    // Print ppairs
    /*for (int i=0; i<parents.size(); i++){
        cout << i << " ";
    }cout <<endl;
    for (int i=0; i<parents.size(); i++){
        for (int j=0; j<parents.size(); j++){
            cout << ppairs[i][j] << " ";
        }cout <<endl;
    }*/

        for (int i=0; i<parents.size(); i++){
           for (int j=0; j<parents.size(); j++){
               ppairs[i][j] = 0;
           }    
        }

        //for ( int i=0; i<children.size(); i++){
        //    calc_pairwise_distance(children[i]);
        //}

        // Print Breeding Stats to file
        //cout <<endl << "--------------------------------------------------" << endl;
        //cout <<"### End of Random Breeding Statistics" << endl; 
        cout <<"Size of offspring ensemble retained: " << children.size() << endl;
        cout <<"Total number of molecules made: " << total_num_mols_made << endl
             <<"Total number with valid torenv: " << valid_torenv << endl
             <<"Total number with valid descriptors: " << valid_descriptors <<endl;
        cout <<"Number of unique pairs sampled " << pair<< " of " << pair_max <<endl; 
        // Print Hungarian similarity data
        cout <<"Total number of xovers prevented by Hungarian Similarity: " << hungarian_similarity_prune <<endl; 

        double stop_time_step2 = time_seconds();
        //cout << "\n" "-----------------------------------" "\n";
        if (verbose) cout << " Elapsed time Breeding:\t" << stop_time_step2 - start_time_step2
             << " seconds for GA breeding\n\n";



        //-------------------------------------------------------------------------
        // STEP 3: Minimize and Pruning
        cout << endl<< "#### Entering Minimization Routine for Offspring ####" << endl;
        double start_time_step3 = time_seconds();
        // Minimize children with valid torenv
        for ( int i=0; i<children.size(); i++){
            // Minimize and increment counter
            minimize_children( children[i], pruned_children, score, simplex, c_typer );
        }
        for ( int i=0; i<pruned_children.size(); i++){
            calc_pairwise_distance(pruned_children[i]);
        }

        //Clear children      
        children.clear();
        // Print number of molecules after minimization
        cout <<"Total number of molecules after minimization: " << pruned_children.size() <<endl;

        cout << endl<<"#### Entering Fitness Pruning Routine for Offspring ####" << endl;
        // This function will select the best scoring conformer and use in the next generation
        fitness_pruning( pruned_children, scored_generation, score, simplex, c_typer );
        //cout << "Total number of molecules after fitness pruning: " << scored_generation.size() <<endl;

        // Clear vectors to free up memory
        pruned_children.clear(); 
        //Check if size is 0 then exit
        if (scored_generation.size() == 0){
           cout << endl << "WARNING: THERE ARE NO OFFSPRING (after crossover and pruning)" <<endl;
           if (ga_mutate_parents == false){
              cout << "ERROR: THE PROGRAM WILL EXIT. Consider using parent mutations."<<endl;
              exit(0);
           }
        }

	// Time for pruning from crossover CS LEP
        double stop_time_step3 = time_seconds();
        if (verbose) cout << "\n" "-----------------------------------" "\n";
        if (verbose) cout << " Elapsed time Minimize + Prune:\t" << stop_time_step3 - start_time_step3
            << " seconds for GA breeding\n\n";
        
    }// end if xover on CS LEP
    // If only completing mutations CS LEP
    if ( !ga_xover_on ){
        // Copy parents to scored_generation
        for (int i=0; i<parents.size(); i++){
            int false_counter=0;
            scored_generation.push_back(parents[i]);
            //cout << "Lauren: pushing parents to scored_gen" << endl;
        }
        activate_vector( scored_generation );
    } // end if !xover CS LEP
 


    //-------------------------------------------------------------------------
    // STEP 4: Mutations
    // Determine which parents did not contribute to the final ensemble

    // Reset used
    for (int i=0; i<parents.size(); i++){
        parents[i].used = false;
    }
    // Search through the comments of each molecule
    for ( int i=0; i<scored_generation.size(); i++ ){
        int p1 = 0;
        int p2 = 0;
        istringstream tmp(scored_generation[i].energy);
        tmp>> p1 >> p2;
        string p11 = parents[p1].title.substr(ga_name_identifier.size()+1);
        string p22 = parents[p2].title.substr(ga_name_identifier.size()+1);
        // For parent mutation, update used flag
        parents[p1].used = true; 
        parents[p2].used = true; 

        // Save parent index information 
        //ostringstream new_comment;
        //new_comment << p1 << "_" << p2;
        scored_generation[i].energy = p11 + "_X_" + p22;
        //new_comment.clear();
    }

    // Print the number of unique parents involved in offspring crossover  
    if ( ga_xover_on ){
        int xover_parents = 0;
        for (int i=0; i<parents.size(); i++){
            // If the molecule was not used for crossover
            if ( parents[i].used == true ){ xover_parents++; }
        }
        //cout <<"#### Total number of parents involved in xover: " << xover_parents <<endl;
    }

    //double stop_time_step3 = time_seconds();
    //if (verbose) cout << "\n" "-----------------------------------" "\n";
    //if (verbose) cout << " Elapsed time Minimize + Prune:\t" << stop_time_step3 - start_time_step3
      //   << " seconds for GA breeding\n\n";

    // STEP 4A: If mutating parents
    if (ga_mutate_parents == true ){
       double start_time_step4 = time_seconds();
       std::vector <DOCKMol> pmuts;
       std::vector <DOCKMol> ppruned; 
       std::vector <DOCKMol> pfinal; 
      
       // Save used status 
       for (int i=0; i<parents.size(); i++){
           // If the molecule was not used for crossover
           if ( parents[i].used == false ){
              pmuts.push_back( parents[i] );
              //cout << "Lauren: pushing parents to pmuts in breeding_rand" << endl;
           }
       }

       // Use random parent sampling:
       if ( (pmuts.size() > 10) ){ 
          cout<< endl << "#### Entering Random Parent Mutations ####" << endl; 
          cout<< "#### Parent Ensemble Size: " << pmuts.size() <<endl;
          master_mut_rand( pmuts, saved_parents, ppruned, true, score, simplex, c_typer, orient, gen_num );
       }
       else{ 
          cout<< endl << "#### Entering Exhaustive Parent Mutations ####" << endl; 
          cout << "Parent Ensemble Size: " << pmuts.size() <<endl;
          master_mut_exhaustive( pmuts, saved_parents, ppruned, true, score, simplex, c_typer, orient, gen_num ); 
       }
       pmuts.clear();

       double stop_time_step4 = time_seconds();
       //cout << "\n" "-----------------------------------" "\n";
       if (verbose) cout << " Elapsed time Parents Mutations:\t" << stop_time_step4 - start_time_step4
            << " seconds for GA breeding\n\n";

       // If there are successful mutations
       double start_time_step5 = time_seconds();
       if (ppruned.size() != 0){
          // After scoring, prune the lot
          fitness_pruning(ppruned, pfinal, score, simplex, c_typer);
          //cout << "pfinal, parents" << endl;
          uniqueness_prune_mut(pfinal, parents, score, c_typer);

          // Combine final mutated parents with original 
          for (int i=0; i<pfinal.size(); i++ ){
              mutated_parents.push_back(pfinal[i]);
          }
          ppruned.clear();
          pfinal.clear();
       }
       saved_parents.clear();

       double stop_time_step5 = time_seconds();
       if (verbose) cout << "\n" "-----------------------------------" "\n";
       if (verbose) cout << " Elapsed time Pruning-Mutate Parents:\t" << stop_time_step5 - start_time_step5
            << " seconds for GA breeding\n\n";
    }

    // STEP 4B: Mutate offspring
    if (scored_generation.size() > 0){
       double start_time_step6 = time_seconds();
       cout<< endl <<  "#### Entering Random Offspring mutations ####" << endl;
       cout << "Offspring Ensemble Size: " << scored_generation.size() <<endl;
       master_mut_rand( scored_generation, saved_parents, pruned_children, false, score, simplex, c_typer, orient, gen_num );

       // Clear the possibility manipulated molecules
       scored_generation.clear();

       double stop_time_step6 = time_seconds();
       if (verbose) cout << "\n" "-----------------------------------" "\n";
       if (verbose) cout << " Elapsed time Offspring Mutations:\t" << stop_time_step6 - start_time_step6
            << " seconds for GA breeding\n\n";

       // STEP 4C: Pruning
       // After scoring, prune the lot
       double start_time_step7 = time_seconds();
       fitness_pruning(pruned_children, scored_generation, score, simplex, c_typer);
  
       // Combine unmanipulated offspring to scored_gen
       // Run All hard filters on crossover offspring
       hard_filter(saved_parents); 
       for (int i=0; i<saved_parents.size(); i++ ){
           scored_generation.push_back(saved_parents[i]);
       }

       // Clear save vector
       saved_parents.clear();
       pruned_children.clear();

       // Combine all molecules onto pruned_children
       for (int i=0; i<scored_generation.size(); i++){
           pruned_children.push_back(scored_generation[i]);
       } 
       scored_generation.clear();

       // Prune again with all offspring and mutated offspring
       fitness_pruning(pruned_children, scored_generation, score, simplex, c_typer);
       pruned_children.clear();
       for (int i = 0; i<scored_generation.size(); i++){
            int falsecounter = 0;
            for (int atom = 0; atom<scored_generation[i].num_atoms; atom++){
                if (scored_generation[i].atom_active_flags[atom] == false) {
                    falsecounter++;
                    //scored_generation[i].atom_active_flags[atom] = true;
                }
            }
        }

       // If mutating the parents
       if ((ga_mutate_parents == true) && (mutated_parents.size() != 0)){
          //cout << "scoredgen, mutated parents"<< endl;    
          uniqueness_prune_mut(scored_generation, mutated_parents, score, c_typer);
       }
       // Prune offspring that are similar to the parents
       cout << "scoredgen, parents" << endl;
       uniqueness_prune_mut(scored_generation, parents, score, c_typer);
       double stop_time_step7 = time_seconds();
       if (verbose) cout << "\n" "-----------------------------------" "\n";
       if (verbose) cout << " Elapsed time Pruning-Mutate Offspring:\t" << stop_time_step7 - start_time_step7
            << " seconds for GA breeding\n\n";
    }
    // Save mutated parents to offspring vector
    if ((ga_mutate_parents == true) && (mutated_parents.size() != 0)){
       for ( int i=0; i<mutated_parents.size(); i++){
           scored_generation.push_back(mutated_parents[i]);
       }
    }
    // Clear necessary vectors
    mutated_parents.clear();

    return;

} // end GA_Recomb::breeding_rand()




            //////////////////////////////////////////
            //*********CROSSOVER FUNCTIONS**********//
            /////////////////////////////////////////



// ++++++++++++++++++++++++++++++++++++++++
// Compare the molecules based on hungarian rmsd
bool
GA_Recomb::similarity_compare( DOCKMol &  mol1, int mol1o, int mol1t, DOCKMol & mol2, int mol2o, int mol2t, bool second, Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & typer )
{
    Trace trace( "GA_Recomb::similarity_compare()" );
    // Set all atom_active_flags to false
    for ( int i=0; i<mol1.num_atoms; i++ ){
        mol1.atom_active_flags[i]=false;
    }
    for ( int i=0; i<mol2.num_atoms; i++ ){
        mol2.atom_active_flags[i]=false;
    }

    // STEP 1: Activate the smaller side of each parent (on the same side)
    // Define vectors to hold the atom and bond indices in activated halves of each parent
    vector <int> atom_vec1;
    vector <int> atom_vec2;

    // STEP 1A: Find the smaller side of the molecule
    // Save all of the atoms connected to the bond origin atom of each parent in their respective vectors
    // (Note: get_atom_children(int, int) obtains the offspring of each bond involving mol1o
    // and excludes mol1o and mol1t)
    atom_vec1 = mol1.get_atom_children( mol1t, mol1o );
    
    // If atom_vec1 is the smaller side then continue
    if ( ( (second != true) && (atom_vec1.size() <= (mol1.num_atoms/2)) ) 
      || ( (second == true) && (atom_vec1.size() >= (mol1.num_atoms/2)) )){
       // Define atom_vec2
       atom_vec2 = mol2.get_atom_children( mol2o, mol2t );
 
       // Activate all atoms
       // Activate the origin bond atoms around the bonds involved in breeding for each parent
       mol1.atom_active_flags[mol1o] = true;

       // Activate the atoms in atom_vec1
       for ( int i=0; i<atom_vec1.size(); i++ ){
           mol1.atom_active_flags[atom_vec1[i]]=true;
       } 

       // Activate the origin bond atoms around the bonds involved in breeding for each parent
       mol2.atom_active_flags[mol2t] = true;
       // Activate the atoms in atom_vec1
       for ( int i=0; i<atom_vec2.size(); i++ ){
           mol2.atom_active_flags[atom_vec2[i]]=true;
       } 
    }
    else{
       // Switch sides
       atom_vec1 = mol1.get_atom_children( mol1o, mol1t );

       // Define atom_vec2
       atom_vec2 = mol2.get_atom_children( mol2t, mol2o );
 
       // Activate all atoms
       // Activate the origin bond atoms around the bonds involved in breeding for each parent
       mol1.atom_active_flags[mol1t] = true;

       // Activate the atoms in atom_vec1
       for ( int i=0; i<atom_vec1.size(); i++ ){
           mol1.atom_active_flags[atom_vec1[i]]=true;
       } 

       // Activate the origin bond atoms around the bonds involved in breeding for each parent
       mol2.atom_active_flags[mol2o] = true;

       // Activate the atoms in atom_vec1
       for ( int i=0; i<atom_vec2.size(); i++ ){
           mol2.atom_active_flags[atom_vec2[i]]=true;
       } 
    }

    // Print molecules to file
    // DELETE - REMOVE
    /*fstream fout_molecules;
    fout_molecules.open ( "zzz.sameside.mol2", fstream::out|fstream::app);
    fout_molecules << mol1.current_data << endl;
    Write_Mol2( mol1, fout_molecules );
    fout_molecules << mol2.current_data <<endl;
    Write_Mol2(  mol2, fout_molecules );
    fout_molecules.close();*/
 

    // STEP 1B: Idenfity whether the active sections are the same
    Hungarian_RMSD h;
    pair <double, int> result;

    result = h.calc_Hungarian_RMSD_dissimilar( mol1, mol2 );
        
    if ((result.first < ga_heur_matched_rmsd) && (result.second < ga_heur_unmatched_num)) { 

       /*fstream fout_molecules;
       activate_mol(mol1);
       activate_mol(mol2);

       fout_molecules.open ( "zzz.hung.mol2", fstream::out|fstream::app);
       fout_molecules << mol1.current_data << endl;
       Write_Mol2( mol1, fout_molecules );
       fout_molecules << mol2.current_data <<endl;
       Write_Mol2(  mol2, fout_molecules );
       fout_molecules.close();*/

       return true;
   }

   return false;
} // end GA_Recomb::similarity_compare()



// +++++++++++++++++++++++++++++++++++++++++
// This function will activate the half of each parent that will be combined in Recomb_mols to form each child
void
GA_Recomb::activate_half_mol( DOCKMol &  mol1, int mol1o, int mol1t, DOCKMol &  mol2, int mol2o, int mol2t,
                              bool replace  )
{
    Trace trace( "GA_Recomb::activate_half_mol()" );
    // Set all atom_active_flags to false
    for ( int i=0; i<mol1.num_atoms; i++ ){
        mol1.atom_active_flags[i]=false;
    }
   
    for ( int i=0; i<mol1.num_bonds; i++ ){
        mol1.bond_keep_flags[i]=false;
    }

    // Activate the origin bond atoms around the bonds involved in breeding for each parent
    if ( replace == false ){ mol1.atom_active_flags[mol1o] = true; }
    else{ mol1.atom_active_flags[mol1t] = true; }
    
    // Define vectors to hold the atom and bond indices in activated halves of each parent
    vector <int> atom_vec1;
    vector <int> atom_vec2;
    vector <int> bond_vec1;
    vector <int> bond_vec2;

    // Save all of the atoms connected to the bond origin atom of each parent in their respective vectors
    // (Note: get_atom_children(int, int) obtains the offspring of each bond involving mol1o
    // and excludes mol1o and mol1t)
    // If completing a replacement, then the origin atom will be the dummy atom
    if (replace == false ){ atom_vec1 = mol1.get_atom_children( mol1t, mol1o ); }
    else{ atom_vec1 = mol1.get_atom_children( mol1o, mol1t ); }

    // Activate the atoms in atom_vec1
    for ( int i=0; i<atom_vec1.size(); i++ ){
        mol1.atom_active_flags[atom_vec1[i]]=true;
    } 
        
    // Identify the bonds between active atoms and add them to bond_vec1
    for ( int i=0; i<atom_vec1.size(); i++ ){
        vector <int> temp_bond_vec1;
        temp_bond_vec1 = mol1.get_bond_neighbors(atom_vec1[i]);


        for ( int j=0; j<temp_bond_vec1.size(); j++ ){       
            bond_vec1.push_back( temp_bond_vec1[j] );
        }
        
    temp_bond_vec1.clear();
    }

    // MARK: Later make bond vec1 unique to not active the same bonds
    for ( int i=0; i<bond_vec1.size(); i++){
        mol1.bond_keep_flags[bond_vec1[i]]=true;
    } 


    // In order to use the same function for crossovers and mutations (where there is only one molecule),
    // an additional step is present to ensure that there is a second molecule
    if ( mol2o != -1 ){
        ///cout << "activate_half_mol for nonmutation" << endl;
        // Set all atom_active_flags to false
        for ( int i=0; i<mol2.num_atoms; i++ ){
            mol2.atom_active_flags[i]=false;
        }
   
        for ( int i=0; i<mol2.num_bonds; i++ ){
            mol2.bond_keep_flags[i]=false;
        }

        // Activate the origin bond atoms around the bonds involved in breeding for each parent
        mol2.atom_active_flags[mol2o] = true;

        // Save all of the atoms connected to the bond origin atom of each parent in their respective vectors
        // (Note: get_atom_children(int, int) obtains the offspring of each bond involving mol1o
        // and excludes mol1o and mol1t)
        atom_vec2 = mol2.get_atom_children( mol2t, mol2o ); 

    
        // Repeat actions above for atom_vec1 and bonc_vec2 for atom_vec2 and bond_vec2
        for ( int i=0; i<atom_vec2.size(); i++ ){
            mol2.atom_active_flags[atom_vec2[i]]=true;
        } 

        for ( int i=0; i<atom_vec2.size(); i++ ){
            vector <int> temp_bond_vec2;
            temp_bond_vec2 = mol2.get_bond_neighbors(atom_vec2[i]);
    
            for ( int j=0; j<temp_bond_vec2.size(); j++ ){       
                bond_vec2.push_back( temp_bond_vec2[j] );
            }

            temp_bond_vec2.clear();
       }

       for ( int i=0; i<bond_vec2.size(); i++ ){
            mol2.bond_keep_flags[bond_vec2[i]]=true;
       } 
    }

    // Clear atom and bond vectors used to active half of each parent mol
    atom_vec1.clear();
    atom_vec2.clear();
    bond_vec1.clear();
    bond_vec2.clear();
 
    return;

} // end GA_Recomb::activate_half_mol() 




// +++++++++++++++++++++++++++++++++++++++++
// Return a covalent radius depending on atom type - called in recomb_mols
float
GA_Recomb::calc_cov_radius( string atom )
{
    Trace trace( "GA_Recomb::calc_cov_radius()" );
    // This function assumes that the atom type is a Sybyl atom type. These covalent
    // radii come from the CRC handbook, and are generalized by element (see header
    // file). MARK - Perhaps we could use some more rigorous numbers?

    if ( atom == "H" )
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
        { cout <<"WARNING: Did not recognize the atom_type " <<atom <<" in DN_GA_Build::calc_cov_radius()\n"; return 0.71; }

} // end GA_Recomb::calc_cov_radius()

    


// +++++++++++++++++++++++++++++++++++++++++
// The objective of this function is to is to translate / rotate activated halve of parents[mol2]
// so that the bonds overlap with the activated half of parents[mol1]
void
GA_Recomb::recomb_mols( DOCKMol & mol1, int mol1o, int mol1t, DOCKMol &  mol2, int mol2o, int mol2t, 
                        std::vector <DOCKMol> & keep)
{
    Trace trace( "GA_Recomb::recomb_mols()" );
    // STEP 1: Adjust the bond length of mol2 to match covalent radii of new pair

    // Calculate the desired bond length and remember as 'new_rad'
    float new_rad = calc_cov_radius(mol1.atom_types[mol1o]) +
                    calc_cov_radius(mol2.atom_types[mol2o]);
    //cout << "new_rad: " << new_rad << endl; 
    // Calculate the x-y-z components of the current mol2 bond vector (bond_vec)
    DOCKVector bond_vec;
    //cout << "mol2.x[mol2o]:" << mol2.x[mol2o] << endl;
    //cout << "mol2.x[mol2t]:" << mol2.x[mol2t] << endl;
    //cout << "mol2.y[mol2o]:" << mol2.y[mol2o] << endl;
    //cout << "mol2.y[mol2t]:" << mol2.y[mol2t] << endl;
    //cout << "mol2.z[mol2o]:" << mol2.z[mol2o] << endl;
    //cout << "mol2.z[mol2t]:" << mol2.z[mol2t] << endl;

    bond_vec.x = mol2.x[mol2o] - mol2.x[mol2t];
    //cout << "bond_vec.x" << bond_vec.x << endl;
    bond_vec.y = mol2.y[mol2o] - mol2.y[mol2t];
    //cout << "bond_vec.y" << bond_vec.y << endl;
    bond_vec.z = mol2.z[mol2o] - mol2.z[mol2t];
    //cout << "bond_vec.z" << bond_vec.z << endl;

    // Normalize the bond vector then multiply each component by new_rad so that it is the desired length
    bond_vec = bond_vec.normalize_vector();
    //cout << "Normalized bond_vec xyz:" << bond_vec.x << "," << bond_vec.y << "," << bond_vec.z << endl;
    
    bond_vec.x *= new_rad;
    bond_vec.y *= new_rad;
    bond_vec.z *= new_rad;

    //cout << "New bond_vec size xyz:" << bond_vec.x << "," << bond_vec.y << "," << bond_vec.z << endl;

    // Define the amount that mol2o should move to form desired bond length
    DOCKVector move_vec;
    move_vec.x = ( mol2.x[mol2o] - mol2.x[mol2t] ) - bond_vec.x;
    move_vec.y = ( mol2.y[mol2o] - mol2.y[mol2t] ) - bond_vec.y;
    move_vec.z = ( mol2.z[mol2o] - mol2.z[mol2t] ) - bond_vec.z;

    //cout << "move_vec xyz:" << move_vec.x << "," << move_vec.y << "," << move_vec.z << endl;
    // Change the coordinates of all the active atoms connected to mol2o in order to have the correct
    // bond length without disturbing the other bond orientations
    for ( int i=0; i<mol2.num_atoms; i++ ){
        if ( mol2.atom_active_flags[i] == true ){       
            mol2.x[i] = mol2.x[i] - move_vec.x;
            mol2.y[i] = mol2.y[i] - move_vec.y;
            mol2.z[i] = mol2.z[i] - move_vec.z;
        }
    }



    // STEP 2: Translate dummy2 of frag2 to the origin

    // Figure out what translation is required to move the dummy atom to the origin
    DOCKVector trans1;
    trans1.x = -mol2.x[mol2t];
    trans1.y = -mol2.y[mol2t];
    trans1.z = -mol2.z[mol2t];
    //cout <<"trans required to move Du to origin xyz: " << trans1.x << ","<< trans1.y << "," << trans1.z << endl;
    // Use the dockmol function to translate the fragment so the dummy atom is at the origin
    mol2.translate_mol(trans1);
    //cout <<"fragment should be at origin xyz: " << trans1.x << ","<< trans1.y << "," << trans1.z << endl;


    // STEP 3: Calculate dot product to determine theta (theta = angle between vec1 and vec2)

    // vec1 = vector pointing from origin to target atom in mol1
    DOCKVector vec1;
    vec1.x = mol1.x[mol1t] - mol1.x[mol1o];
    vec1.y = mol1.y[mol1t] - mol1.y[mol1o];
    vec1.z = mol1.z[mol1t] - mol1.z[mol1o];

    // vec2 = vector pointing from origin to target in mol2
    DOCKVector vec2;
    vec2.x = mol2.x[mol2o] - mol2.x[mol2t];
    vec2.y = mol2.y[mol2o] - mol2.y[mol2t];
    vec2.z = mol2.z[mol2o] - mol2.z[mol2t];

    // Declare some variables
    float dot;          // dot product value of vec1 and vec2
    float vec1_magsq;   // vec1 magnitude-squared
    float vec2_magsq;   // vec2 magnitude-squared
    float cos_theta;    // cosine of theta
    float sin_theta;    // sine of theta

    // Compute the dot product using the function in utils.cpp
    dot = dot_prod( vec1, vec2 );

    // Compute these magnitudes (squared)
    vec1_magsq = ( vec1.x * vec1.x ) + ( vec1.y * vec1.y ) + ( vec1.z * vec1.z );
    vec2_magsq = ( vec2.x * vec2.x ) + ( vec2.y * vec2.y ) + ( vec2.z * vec2.z );

    // Compute cosine and sine of theta (theta itself is not actually calculated)
     cos_theta = dot / ( sqrt ( vec1_magsq * vec2_magsq ) );
     //cout << "cos_theta:" << cos_theta << endl;
     if (cos_theta >= 1.0){
        cos_theta = 1.0;
        sin_theta = 0.0;
     }
     if (cos_theta <= -1.0){
        cos_theta = -1.0;
        sin_theta = 0.0;
     }
     sin_theta = sqrt ( 1 - ( cos_theta * cos_theta ) );
     //cout << "sin_theta:" << sin_theta << endl;



    // STEP 4: Rotate vec2 to be coincident with vec1
    // If cos_theta is less than the cut off, enter the rotation function
    float cos_cut_off = 0.9999999404;
    if ( cos_theta < cos_cut_off ){
        rotate( cos_theta, sin_theta, vec1, mol2, vec2 );
    }


    
    // STEP 5: Translate frag2 to frag1

    // This is the translation vector to move mol2o to mol1o (mol2o is at the origin)
    DOCKVector trans2;
    trans2.x = mol1.x[mol1o];
    trans2.y = mol1.y[mol1o];
    trans2.z = mol1.z[mol1o];

    // Use the dockmol function to translate frag2
    mol2.translate_mol( trans2 );

    int bond_num1 = mol1.get_bond( mol1o, mol1t );
    int bond_num2 = mol2.get_bond( mol2o, mol2t );


    //cout << "Entering recomb preattach" << endl;
    //cout << "...for mol1" <<endl;
    //calc_pairwise_distance(mol1);
    //cout << "...for mol2" << endl;
    //calc_pairwise_distance(mol2);      
    // STEP 6: Connect the two fragments and return one object
    DOCKMol tmp_child = attach( mol1, bond_num1, mol2, bond_num2 );
    int counter_active = 0;
    for ( int i=0; i<tmp_child.num_atoms; i++){
        if (tmp_child.atom_active_flags[i] == true){
           counter_active++;
        }
    }
    keep.push_back( tmp_child );
    tmp_child.clear_molecule();
    
    //cout << "Entering recomb postattach" << endl;
    //for ( int i=0; i<keep.size(); i++){
    //        calc_pairwise_distance(keep[i]);
    //}

    
    return;

}  // end GA_Recomb::recomb_mols()
// +++++++++++++++++++++++++++++++++++++++++
// The objective of this function is to is to translate / rotate activated halve of parents[mol2]
// so that the bonds overlap with the activated half of parents[mol1]
// JDB This version of the function is specifically so that crossover can take advantage
// of the new Do Not Modify (DNM) mechanism JDB
void
GA_Recomb::recomb_mols_xover( DOCKMol & mol1, int mol1o, int mol1t, DOCKMol &  mol2, int mol2o, int mol2t, 
                        std::vector <DOCKMol> & keep)
{
    Trace trace( "GA_Recomb::recomb_mols_xover()" );
    // JDB STEP 0: Checks for DNM logic on active and inactive halves of both parents
    if (dnm_enabled){
        //Double check to make sure this matters for this set of parents
        if ( mol1.mol_dnm_flag || mol2.mol_dnm_flag ) {
            int child_dnm_counter=0;
            int p1_dnm_counter=0;
            int p2_dnm_counter=0;
            //Checks all inactive of mol1 for DNM - doesn't accept if so.
            //Keeps a counter of all active DNM atoms along the way for a final
            //check to make sure the molecule is valid (ie doesn't GAIN any DNM atoms)
            if ( mol1.mol_dnm_flag ) {
                for ( int i=0; i<mol1.num_atoms; i++){
                //    //If the atom is inactive AND marked as DNM
                //    if ( mol1.atom_dnm_flag[i] && !mol1.atom_active_flags[i] ){
                //        cout << "Crossover would result in loss of DNM segment. This crossover is invalid." << endl;
                //        return;
//
                //    //If the atom is active AND marked as DNM
                //    } else if ( mol1.atom_dnm_flag[i] && mol1.atom_active_flags[i] ){
                //        child_dnm_counter++;
                //        p1_dnm_counter++;
                //    }
                //}
                    if ( mol1.atom_dnm_flag[i] && mol1.atom_active_flags[i] ){
                        child_dnm_counter++;
                        p1_dnm_counter++;
                    }
                }
            }
            //Checks all inactive of mol2 for DNM - doesn't accept if so.
            if ( mol2.mol_dnm_flag ) {
                for ( int i=0; i<mol2.num_atoms; i++ ){
                   ////If the atom is inactive AND marked as DNM
                   //if ( mol2.atom_dnm_flag[i] && !mol2.atom_active_flags[i] ){
                   //    cout << "Crossover would result in loss of DNM segment. This crossover is invalid." << endl;
                   //    return;   
                   ////If the atom is active AND marked as DNM
                   //} else if ( mol2.atom_dnm_flag[i] && mol2.atom_active_flags[i] ){
                   //    child_dnm_counter++;
                   //    p2_dnm_counter++;
                   //}

                    if ( mol2.atom_dnm_flag[i] && mol2.atom_active_flags[i] ){
                            child_dnm_counter++;
                            p2_dnm_counter++;
                    }
                }
            }
            if (child_dnm_counter != p1_dnm_counter || child_dnm_counter != p2_dnm_counter){
                cout << "Crossover would result in a differing number of DNM atoms. This crossover is invalid." << endl;
                return;
            }
        }
    } // End Step 0
    
    // STEP 1: Adjust the bond length of mol2 to match covalent radii of new pair

    // Calculate the desired bond length and remember as 'new_rad'
    float new_rad = calc_cov_radius(mol1.atom_types[mol1o]) +
                    calc_cov_radius(mol2.atom_types[mol2o]);
    //cout << "new_rad: " << new_rad << endl; 
    // Calculate the x-y-z components of the current mol2 bond vector (bond_vec)
    DOCKVector bond_vec;
    //cout << "mol2.x[mol2o]:" << mol2.x[mol2o] << endl;
    //cout << "mol2.x[mol2t]:" << mol2.x[mol2t] << endl;
    //cout << "mol2.y[mol2o]:" << mol2.y[mol2o] << endl;
    //cout << "mol2.y[mol2t]:" << mol2.y[mol2t] << endl;
    //cout << "mol2.z[mol2o]:" << mol2.z[mol2o] << endl;
    //cout << "mol2.z[mol2t]:" << mol2.z[mol2t] << endl;

    bond_vec.x = mol2.x[mol2o] - mol2.x[mol2t];
    //cout << "bond_vec.x" << bond_vec.x << endl;
    bond_vec.y = mol2.y[mol2o] - mol2.y[mol2t];
    //cout << "bond_vec.y" << bond_vec.y << endl;
    bond_vec.z = mol2.z[mol2o] - mol2.z[mol2t];
    //cout << "bond_vec.z" << bond_vec.z << endl;

    // Normalize the bond vector then multiply each component by new_rad so that it is the desired length
    bond_vec = bond_vec.normalize_vector();
    //cout << "Normalized bond_vec xyz:" << bond_vec.x << "," << bond_vec.y << "," << bond_vec.z << endl;
    
    bond_vec.x *= new_rad;
    bond_vec.y *= new_rad;
    bond_vec.z *= new_rad;

    //cout << "New bond_vec size xyz:" << bond_vec.x << "," << bond_vec.y << "," << bond_vec.z << endl;

    // Define the amount that mol2o should move to form desired bond length
    DOCKVector move_vec;
    move_vec.x = ( mol2.x[mol2o] - mol2.x[mol2t] ) - bond_vec.x;
    move_vec.y = ( mol2.y[mol2o] - mol2.y[mol2t] ) - bond_vec.y;
    move_vec.z = ( mol2.z[mol2o] - mol2.z[mol2t] ) - bond_vec.z;

    //cout << "move_vec xyz:" << move_vec.x << "," << move_vec.y << "," << move_vec.z << endl;
    // Change the coordinates of all the active atoms connected to mol2o in order to have the correct
    // bond length without disturbing the other bond orientations
    for ( int i=0; i<mol2.num_atoms; i++ ){
        if ( mol2.atom_active_flags[i] == true ){       
            mol2.x[i] = mol2.x[i] - move_vec.x;
            mol2.y[i] = mol2.y[i] - move_vec.y;
            mol2.z[i] = mol2.z[i] - move_vec.z;
        }
    }

    // STEP 2: Translate dummy2 of frag2 to the origin

    // Figure out what translation is required to move the dummy atom to the origin
    DOCKVector trans1;
    trans1.x = -mol2.x[mol2t];
    trans1.y = -mol2.y[mol2t];
    trans1.z = -mol2.z[mol2t];
    //cout <<"trans required to move Du to origin xyz: " << trans1.x << ","<< trans1.y << "," << trans1.z << endl;
    // Use the dockmol function to translate the fragment so the dummy atom is at the origin
    mol2.translate_mol(trans1);
    //cout <<"fragment should be at origin xyz: " << trans1.x << ","<< trans1.y << "," << trans1.z << endl;


    // STEP 3: Calculate dot product to determine theta (theta = angle between vec1 and vec2)

    // vec1 = vector pointing from origin to target atom in mol1
    DOCKVector vec1;
    vec1.x = mol1.x[mol1t] - mol1.x[mol1o];
    vec1.y = mol1.y[mol1t] - mol1.y[mol1o];
    vec1.z = mol1.z[mol1t] - mol1.z[mol1o];

    // vec2 = vector pointing from origin to target in mol2
    DOCKVector vec2;
    vec2.x = mol2.x[mol2o] - mol2.x[mol2t];
    vec2.y = mol2.y[mol2o] - mol2.y[mol2t];
    vec2.z = mol2.z[mol2o] - mol2.z[mol2t];

    // Declare some variables
    float dot;          // dot product value of vec1 and vec2
    float vec1_magsq;   // vec1 magnitude-squared
    float vec2_magsq;   // vec2 magnitude-squared
    float cos_theta;    // cosine of theta
    float sin_theta;    // sine of theta

    // Compute the dot product using the function in utils.cpp
    dot = dot_prod( vec1, vec2 );

    // Compute these magnitudes (squared)
    vec1_magsq = ( vec1.x * vec1.x ) + ( vec1.y * vec1.y ) + ( vec1.z * vec1.z );
    vec2_magsq = ( vec2.x * vec2.x ) + ( vec2.y * vec2.y ) + ( vec2.z * vec2.z );

    // Compute cosine and sine of theta (theta itself is not actually calculated)
     cos_theta = dot / ( sqrt ( vec1_magsq * vec2_magsq ) );
     //cout << "cos_theta:" << cos_theta << endl;
     if (cos_theta >= 1.0){
        cos_theta = 1.0;
        sin_theta = 0.0;
     }
     if (cos_theta <= -1.0){
        cos_theta = -1.0;
        sin_theta = 0.0;
     }
     sin_theta = sqrt ( 1 - ( cos_theta * cos_theta ) );
     //cout << "sin_theta:" << sin_theta << endl;



    // STEP 4: Rotate vec2 to be coincident with vec1
    // If cos_theta is less than the cut off, enter the rotation function
    float cos_cut_off = 0.9999999404;
    if ( cos_theta < cos_cut_off ){
        rotate( cos_theta, sin_theta, vec1, mol2, vec2 );
    }

    // STEP 5: Translate frag2 to frag1

    // This is the translation vector to move mol2o to mol1o (mol2o is at the origin)
    DOCKVector trans2;
    trans2.x = mol1.x[mol1o];
    trans2.y = mol1.y[mol1o];
    trans2.z = mol1.z[mol1o];

    // Use the dockmol function to translate frag2
    mol2.translate_mol( trans2 );

    int bond_num1 = mol1.get_bond( mol1o, mol1t );
    int bond_num2 = mol2.get_bond( mol2o, mol2t );


    //cout << "Entering recomb preattach" << endl;
    //cout << "...for mol1" <<endl;
    //calc_pairwise_distance(mol1);
    //cout << "...for mol2" << endl;
    //calc_pairwise_distance(mol2);  


    // STEP 6: Connect the two fragments and return one object

    DOCKMol tmp_child = attach( mol1, bond_num1, mol2, bond_num2 );
    int counter_active = 0;
    for ( int i=0; i<tmp_child.num_atoms; i++){
        if (tmp_child.atom_active_flags[i] == true){
           counter_active++;
        }
    }

    
    if (dnm_enabled) {
        set_DNM_bools ( tmp_child );
    }
    keep.push_back( tmp_child );
    tmp_child.clear_molecule();
    
    //cout << "Entering recomb postattach" << endl;
    //for ( int i=0; i<keep.size(); i++){
    //        calc_pairwise_distance(keep[i]);
    //}

    
    return;

}  // end GA_Recomb::recomb_mols_xover()



// +++++++++++++++++++++++++++++++++++++++++
// Rotate function that completes Step 4 of the recomb_mols function - obtained function from conf_gen_den.cpp
void
GA_Recomb::rotate( float cos_theta, float sin_theta, DOCKVector vec1, DOCKMol &  mol2, DOCKVector vec2 )
{
    Trace trace( "GA_Recomb::rotate()" );
    // If cos_theta is -1, the vectors are parallel but in the same direction
    // therefore, rotate mol2 180 degrees
    if ( cos_theta == -1 ) {
       // Declare the rotation matrix and rotate  mol2
       double finalmat[3][3] = { { -1, 0, 0}, {0, -1, 0}, {0, 0, -1} };
       mol2.rotate_mol( finalmat );
    }

    // If cos_theta is 1, vec1 and vec2 are already parallel and in opposite directions - only translation is needed.
    // Otherwise, enter this loop and calculate out how to rotate frag2
    else if ( cos_theta != -1 ) {

        //cout << "entered not 1, cos_theta= " << cos_theta << endl; /* MARK - DELETE */

        // Calculate cross product of vec1 and vec2 to get U (function from utils.cpp)
        DOCKVector normalU = cross_prod( vec1, vec2 );

        // Calculate cross product of vec2 and U to get ~W
        DOCKVector normalW = cross_prod( vec2, normalU );

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


        // (2) Make rotation matrix, which rotates vec2 theta angle on a plane of vec2 and normalW
        // to the direction of normalW
        float planeRot[3][3];
        planeRot[0][0] =  cos_theta;  planeRot[0][1] = sin_theta;  planeRot[0][2] = 0;
        planeRot[1][0] = -sin_theta;  planeRot[1][1] = cos_theta;  planeRot[1][2] = 0;
        planeRot[2][0] =          0;  planeRot[2][1] =         0;  planeRot[2][2] = 1;


        // (3) Make inverse  matrix of coorRot matrix - since coorRot is an orthogonal matrix,
        // the inverse is its transpose, (coorRot)^T
        float invcoorRot[3][3];
        invcoorRot[0][0] = coorRot[0][0];  invcoorRot[0][1] = coorRot[1][0];  invcoorRot[0][2] = coorRot[2][0];
        invcoorRot[1][0] = coorRot[0][1];  invcoorRot[1][1] = coorRot[1][1];  invcoorRot[1][2] = coorRot[2][1];
        invcoorRot[2][0] = coorRot[0][2];  invcoorRot[2][1] = coorRot[1][2];  invcoorRot[2][2] = coorRot[2][2];

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

        // Rotate mol2 using finalmat[3][3]
        mol2.rotate_mol(finalmat);
    }
    return;

} //end GA_recomb::rotate()




// +++++++++++++++++++++++++++++++++++++++++
// The attach function will make a new bond between the active halves of each parent
// and create a new DOCKMOL for the child molcule - called in recomb_mols
DOCKMol
GA_Recomb::attach( DOCKMol & mol1, int bond_num1, DOCKMol &  mol2, int bond_num2 )
{
    Trace trace( "GA_Recomb::attach()" );
    // Vectors to hold active atoms and bonds from each parent
    vector <int> atom_vec1;
    vector <int> atom_vec2;
    vector <int> bond_vec1;
    vector <int> bond_vec2;

    // Populate each vector with the active atoms and then bonds for each parent
    for ( int i=0; i<mol1.num_atoms; i++ ){ 
        if ( mol1.atom_active_flags[i] == true ) {
            atom_vec1.push_back( mol1.atom_active_flags[i] );
        }
    }

    for (  int i=0; i<mol1.num_bonds; i++ ){ 
        if ( mol1.bond_keep_flags[i] == true ) {
            bond_vec1.push_back( mol1.bond_keep_flags[i] );
        }
    }

    for ( int i=0; i<mol2.num_atoms; i++ ){
        if ( mol2.atom_active_flags[i] == true ) {
            atom_vec2.push_back( mol2.atom_active_flags[i] );
        }
    }
 
    for ( int i=0; i<mol2.num_bonds; i++ ){ 
        if ( mol2.bond_keep_flags[i] == true ) {
            bond_vec2.push_back( mol2.bond_keep_flags[i] );
        }

    }
  
    // Declare the dockmol that will be the combination of the two active halves
    DOCKMol combined_mol;

    // Number of atoms and bonds that will be in newly combined mol
    // Same number of atoms as tmp_parents and one more bond
    int natoms, nbonds;
    natoms = ( atom_vec1.size() + atom_vec2.size() );
    nbonds = ( bond_vec1.size() + bond_vec2.size() + 1 );

    // Allocate arrays and make assignments for <Tripos>Molecule section of new DOCKMol
    combined_mol.allocate_arrays( natoms, nbonds, 1 );

    // MARK: make a title for the new molecules
    combined_mol.title         =  mol1.title;
    combined_mol.comment1      =  mol1.comment1;
    combined_mol.comment2      =  mol1.comment2;
    combined_mol.comment3      =  mol1.comment3;
    combined_mol.energy        =  mol1.energy;
    
    // Fill in some more of the data
    combined_mol.num_atoms = natoms;
    combined_mol.num_bonds = nbonds;

    if (mol1.mol_dnm_flag || mol2.mol_dnm_flag){
        combined_mol.mol_dnm_flag = true; // JDB
    }
    
    // Initialize renumber vectors
    int atom_index = 0;
    vector <int> renumber1;
    vector <int> renumber2;

    // Renumber the active atoms of both tmp_parents so that the offspring has continous atom numbering
    renumber1.resize( mol1.num_atoms + 1 );
    renumber2.resize( mol2.num_atoms + 1 );
   
    ostringstream res_num;
    // Reset Atom residue numbers to prevent issues with DN minimization - rot bond initialization
    res_num << "1";
    for ( int i=0; i<mol1.num_atoms; i++ ){
        if ( mol1.atom_active_flags[i] == true ){
            combined_mol.x[atom_index]                    = mol1.x[i];
            combined_mol.y[atom_index]                    = mol1.y[i];
            combined_mol.z[atom_index]                    = mol1.z[i];

            combined_mol.charges[atom_index]              = mol1.charges[i];
            combined_mol.atom_names[atom_index]           = mol1.atom_names[i];
            combined_mol.atom_types[atom_index]           = mol1.atom_types[i];
            combined_mol.atom_number[atom_index]          = mol1.atom_number[i];
            combined_mol.atom_residue_numbers[atom_index] = res_num.str();
            combined_mol.subst_names[atom_index]          = mol1.subst_names[i];
            combined_mol.flag_acceptor[atom_index]        = mol1.flag_acceptor[i]; 
            combined_mol.flag_donator[atom_index]         = mol1.flag_donator[i]; 
            combined_mol.acc_heavy_atomid[atom_index]     = mol1.acc_heavy_atomid[i]; 

            if (dnm_enabled) {
                combined_mol.atom_dnm_flag[atom_index]        = mol1.atom_dnm_flag[i]; // JDB
            }
            
            renumber1[i] = atom_index;
            atom_index++;
        }
    }

    //Read atom info of active atom information for mol2
    for ( int i=0; i<mol2.num_atoms; i++ ){
        if ( mol2.atom_active_flags[i] == true ){
            combined_mol.x[atom_index]                    = mol2.x[i];
            combined_mol.y[atom_index]                    = mol2.y[i];
            combined_mol.z[atom_index]                    = mol2.z[i];
            
            combined_mol.charges[atom_index]              = mol2.charges[i];
            combined_mol.atom_names[atom_index]           = mol2.atom_names[i];
            combined_mol.atom_types[atom_index]           = mol2.atom_types[i];
            combined_mol.atom_number[atom_index]          = mol2.atom_number[i];
            combined_mol.atom_residue_numbers[atom_index] = res_num.str();
            combined_mol.subst_names[atom_index]          = mol2.subst_names[i];
            combined_mol.flag_acceptor[atom_index]        = mol2.flag_acceptor[i]; 
            combined_mol.flag_donator[atom_index]         = mol2.flag_donator[i]; 
            combined_mol.acc_heavy_atomid[atom_index]     = mol2.acc_heavy_atomid[i]; 

            if (dnm_enabled) {
                combined_mol.atom_dnm_flag[atom_index]        = mol2.atom_dnm_flag[i]; // JDB
            }
            

            renumber2[i] = atom_index;
            atom_index++;
        }
    }
 
    // Amend bond information from mol1 to combined_mol
    // bond_index serves as a place holder to determine the bond number of the new bond
    int bond_index = 0;
    string new_bond_type = "";
    for ( int i=0; i<mol1.num_bonds; i++ ){
        if ( mol1.bond_keep_flags[i] == true ){

            combined_mol.bonds_origin_atom[bond_index] = renumber1[mol1.bonds_origin_atom[i]];
            combined_mol.bonds_target_atom[bond_index] = renumber1[mol1.bonds_target_atom[i]];
            combined_mol.bond_types[bond_index]        = mol1.bond_types[i];
            combined_mol.bond_ring_flags[bond_index]   = mol1.bond_ring_flags[i];
            combined_mol.bond_keep_flags[bond_index]   = mol1.bond_keep_flags[i];

            bond_index++;

        }
    }

    // Now look through all the bonds of the seocnd parent
    for ( int i=0; i<mol2.num_bonds; i++ ){
        if ( mol2.bond_keep_flags[i] == true ){

            combined_mol.bonds_origin_atom[bond_index] = renumber2[mol2.bonds_origin_atom[i]];
            combined_mol.bonds_target_atom[bond_index] = renumber2[mol2.bonds_target_atom[i]];
            combined_mol.bond_types[bond_index]        = mol2.bond_types[i];
            combined_mol.bond_ring_flags[bond_index]   = mol2.bond_ring_flags[i];
            combined_mol.bond_keep_flags[bond_index]   = mol2.bond_keep_flags[i];

            bond_index++;

        }
    }

    // Determine the bond type for  the last bond which will be the connection between the two active halves
    new_bond_type = mol1.bond_types[bond_num1];  

    // This set of if-else statements will determine which atoms will be the origin and target atoms for the new bond
    if ( renumber1[mol1.bonds_origin_atom[bond_num1]] != 0 ){
        combined_mol.bonds_origin_atom[bond_index] = renumber1[mol1.bonds_origin_atom[bond_num1]];
    }
    else {
        combined_mol.bonds_origin_atom[bond_index] = renumber1[mol1.bonds_target_atom[bond_num1]];}

    
    if ( renumber2[mol2.bonds_origin_atom[bond_num2]] != 0 ){
        combined_mol.bonds_target_atom[bond_index] = renumber2[mol2.bonds_origin_atom[bond_num2]];
     }
     else {
         combined_mol.bonds_target_atom[bond_index] = renumber2[mol2.bonds_target_atom[bond_num2]];}
    

    // Set the bond type of that new bond here
    combined_mol.bond_types[bond_index] = new_bond_type;

    // Make a new mol_info_line that contains total number of atoms and bonds
    stringstream mol_info;
    mol_info <<" mol info " <<natoms <<"  " <<nbonds;
    combined_mol.mol_info_line = mol_info.str();
    mol_info.clear();
    // Update ring flags
    combined_mol.id_ring_atoms_bonds();
   

    atom_vec1.clear();
    atom_vec2.clear();
    bond_vec1.clear();
    bond_vec2.clear();
    renumber1.clear();
    renumber2.clear();
 
    //MARK - DELETE
    /*fstream fout_molecules;
    fout_molecules.open ( "zzz.xover.mol2", fstream::out|fstream::app );
    fout_molecules <<combined_mol.current_data << endl;
    Write_Mol2(combined_mol, fout_molecules);
    fout_molecules.close();*/

    return combined_mol;

} // end GA_Recomb::attach()




// +++++++++++++++++++++++++++++++++++++++++
// The objective of this function is to is to translate / rotate activated halve of parents[mol2]
// so that the bonds overlap with the activated half of parents[mol1]
// Will return the new atom numbers for dummy and heavy atoms provided
// If none, provided, will return a null vector
std::vector <std::pair <int, int> >
GA_Recomb::replace_combine_mols( DOCKMol & mol1, int mol1o, int mol1t, DOCKMol &  mol2, int mol2o, int mol2t, 
                         std::vector < std::pair <int,int> > & du_heavy, std::vector <DOCKMol> & keep)
{
    Trace trace( "GA_Recomb::replace_combine_mols()" );
    // ** ROTATION AND TRANSLATION ** //
    // STEP 1: Adjust the bond length of mol2 to match covalent radii of new pair

    // Calculate the desired bond length and remember as 'new_rad'
    float new_rad = calc_cov_radius(mol1.atom_types[mol1o]) +
                    calc_cov_radius(mol2.atom_types[mol2o]);
 
    // Calculate the x-y-z components of the current mol2 bond vector (bond_vec)
    DOCKVector bond_vec;
    bond_vec.x = mol2.x[mol2o] - mol2.x[mol2t];
    bond_vec.y = mol2.y[mol2o] - mol2.y[mol2t];
    bond_vec.z = mol2.z[mol2o] - mol2.z[mol2t];

    // Normalize the bond vector then multiply each component by new_rad so that it is the desired length
    bond_vec = bond_vec.normalize_vector();
    bond_vec.x *= new_rad;
    bond_vec.y *= new_rad;
    bond_vec.z *= new_rad;

    // Define the amount that mol2o should move to form desired bond length
    DOCKVector move_vec;
    move_vec.x = ( mol2.x[mol2o] - mol2.x[mol2t] ) - bond_vec.x;
    move_vec.y = ( mol2.y[mol2o] - mol2.y[mol2t] ) - bond_vec.y;
    move_vec.z = ( mol2.z[mol2o] - mol2.z[mol2t] ) - bond_vec.z;

    // Change the coordinates of all the active atoms connected to mol2o in order to have the correct
    // bond length without disturbing the other bond orientations
    for ( int i=0; i<mol2.num_atoms; i++ ){
        if ( mol2.atom_active_flags[i] == true ){       
            mol2.x[i] = mol2.x[i] - move_vec.x;
            mol2.y[i] = mol2.y[i] - move_vec.y;
            mol2.z[i] = mol2.z[i] - move_vec.z;
        }
    }



    // STEP 2: Translate dummy2 of frag2 to the origin

    // Figure out what translation is required to move the dummy atom to the origin
    DOCKVector trans1;
    trans1.x = -mol2.x[mol2t];
    trans1.y = -mol2.y[mol2t];
    trans1.z = -mol2.z[mol2t];


    // Use the dockmol function to translate the fragment so the dummy atom is at the origin
    mol2.translate_mol(trans1);


    // STEP 3: Calculate dot product to determine theta (theta = angle between vec1 and vec2)

    // vec1 = vector pointing from origin to target atom in mol1
    DOCKVector vec1;
    vec1.x = mol1.x[mol1t] - mol1.x[mol1o];
    vec1.y = mol1.y[mol1t] - mol1.y[mol1o];
    vec1.z = mol1.z[mol1t] - mol1.z[mol1o];

    // vec2 = vector pointing from origin to target in mol2
    DOCKVector vec2;
    vec2.x = mol2.x[mol2o] - mol2.x[mol2t];
    vec2.y = mol2.y[mol2o] - mol2.y[mol2t];
    vec2.z = mol2.z[mol2o] - mol2.z[mol2t];

    // Declare some variables
    float dot;          // dot product value of vec1 and vec2
    float vec1_magsq;   // vec1 magnitude-squared
    float vec2_magsq;   // vec2 magnitude-squared
    float cos_theta;    // cosine of theta
    float sin_theta;    // sine of theta

    // Compute the dot product using the function in utils.cpp
    dot = dot_prod( vec1, vec2 );

    // Compute these magnitudes (squared)
    vec1_magsq = ( vec1.x * vec1.x ) + ( vec1.y * vec1.y ) + ( vec1.z * vec1.z );
    vec2_magsq = ( vec2.x * vec2.x ) + ( vec2.y * vec2.y ) + ( vec2.z * vec2.z );

    // Compute cosine and sine of theta (theta itself is not actually calculated)
    cos_theta = dot / ( sqrt ( vec1_magsq * vec2_magsq ) );
    sin_theta = sqrt ( 1 - ( cos_theta * cos_theta ) );



    // STEP 4: Rotate vec2 to be coincident with vec1
    // If cos_theta is less than the cut off, enter the rotation function
    float cos_cut_off = 0.9999999404;
    if ( cos_theta < cos_cut_off ){
        rotate( cos_theta, sin_theta, vec1, mol2, vec2 );
    }

    
    // STEP 5: Translate frag2 to frag1

    // This is the translation vector to move mol2o to mol1o (mol2o is at the origin)
    DOCKVector trans2;
    trans2.x = mol1.x[mol1o];
    trans2.y = mol1.y[mol1o];
    trans2.z = mol1.z[mol1o];

    // Use the dockmol function to translate frag2
    mol2.translate_mol( trans2 );

    int bond_num1 = mol1.get_bond( mol1o, mol1t );
    int bond_num2 = mol2.get_bond( mol2o, mol2t );



    // ** RECOMBINATION ** //
    // STEP 6: Attachment
    // Vectors to hold active atoms and bonds from each parent
    vector <int> atom_vec1;
    vector <int> atom_vec2;
    vector <int> bond_vec1;
    vector <int> bond_vec2;

    // Populate each vector with the active atoms and then bonds for each parent
    for ( int i=0; i<mol1.num_atoms; i++ ){ 
        if ( mol1.atom_active_flags[i] == true ) {
            atom_vec1.push_back( mol1.atom_active_flags[i] );
        }
    }

    for (  int i=0; i<mol1.num_bonds; i++ ){ 
        if ( mol1.bond_keep_flags[i] == true ) {
            bond_vec1.push_back( mol1.bond_keep_flags[i] );
        }
    }

    for ( int i=0; i<mol2.num_atoms; i++ ){
        if ( mol2.atom_active_flags[i] == true ) {
            atom_vec2.push_back( mol2.atom_active_flags[i] );
        }
    }
 
    for ( int i=0; i<mol2.num_bonds; i++ ){ 
        if ( mol2.bond_keep_flags[i] == true ) {
            bond_vec2.push_back( mol2.bond_keep_flags[i] );
        }

    }
  
    // Declare the dockmol that will be the combination of the two active halves
    DOCKMol combined_mol;

    // Number of atoms and bonds that will be in newly combined mol
    // Same number of atoms as tmp_parents and one more bond
    int natoms, nbonds;
    natoms = ( atom_vec1.size() + atom_vec2.size() );
    nbonds = ( bond_vec1.size() + bond_vec2.size() + 1 );

    // Allocate arrays and make assignments for <Tripos>Molecule section of new DOCKMol
    combined_mol.allocate_arrays( natoms, nbonds, 1 );

    // MARK: make a title for the new molecules
    combined_mol.title         =  mol1.title;
    combined_mol.comment1      =  mol1.comment1;
    combined_mol.comment2      =  mol1.comment2;
    combined_mol.comment3      =  mol1.comment3;
    combined_mol.energy        =  mol1.energy;
    
    // Fill in some more of the data
    combined_mol.num_atoms = natoms;
    combined_mol.num_bonds = nbonds;

    if (mol1.mol_dnm_flag || mol2.mol_dnm_flag){ // JDB
       combined_mol.mol_dnm_flag = true;
    }
    
    // Initialize renumber vectors
    int atom_index = 0;
    vector <int> renumber1;
    vector <int> renumber2;

    // Renumber the active atoms of both tmp_parents so that the offspring has continous atom numbering
    renumber1.resize( mol1.num_atoms + 1 );
    renumber2.resize( mol2.num_atoms + 1 );
   
    ostringstream res_num;
    // Reset Atom residue numbers to prevent issues with DN minimization - rot bond initialization
    res_num << "1";

    for ( int i=0; i<mol1.num_atoms; i++ ){
        if ( mol1.atom_active_flags[i] == true ){
            combined_mol.x[atom_index]                    = mol1.x[i];
            combined_mol.y[atom_index]                    = mol1.y[i];
            combined_mol.z[atom_index]                    = mol1.z[i];

            combined_mol.charges[atom_index]              = mol1.charges[i];
            combined_mol.atom_names[atom_index]           = mol1.atom_names[i];
            combined_mol.atom_types[atom_index]           = mol1.atom_types[i];
            combined_mol.atom_number[atom_index]          = mol1.atom_number[i];
            combined_mol.atom_residue_numbers[atom_index] = res_num.str();
            combined_mol.subst_names[atom_index]          = mol1.subst_names[i];
            combined_mol.flag_acceptor[atom_index]        = mol1.flag_acceptor[i]; 
            combined_mol.flag_donator[atom_index]         = mol1.flag_donator[i]; 
            combined_mol.acc_heavy_atomid[atom_index]     = mol1.acc_heavy_atomid[i]; 

            if (dnm_enabled){
                combined_mol.atom_dnm_flag[atom_index]        = mol1.atom_dnm_flag[i]; // JDB
            }

            renumber1[i] = atom_index;
            atom_index++;
        }
    }

    //Read atom info of active atom information for mol2
    for ( int i=0; i<mol2.num_atoms; i++ ){
        if ( mol2.atom_active_flags[i] == true ){
            combined_mol.x[atom_index]                    = mol2.x[i];
            combined_mol.y[atom_index]                    = mol2.y[i];
            combined_mol.z[atom_index]                    = mol2.z[i];
            
            combined_mol.charges[atom_index]              = mol2.charges[i];
            combined_mol.atom_names[atom_index]           = mol2.atom_names[i];
            combined_mol.atom_types[atom_index]           = mol2.atom_types[i];
            combined_mol.atom_number[atom_index]          = mol2.atom_number[i];
            combined_mol.atom_residue_numbers[atom_index] = res_num.str();
            combined_mol.subst_names[atom_index]          = mol2.subst_names[i];
            combined_mol.flag_acceptor[atom_index]        = mol2.flag_acceptor[i]; 
            combined_mol.flag_donator[atom_index]         = mol2.flag_donator[i]; 
            combined_mol.acc_heavy_atomid[atom_index]     = mol2.acc_heavy_atomid[i]; 

            if (dnm_enabled){
                combined_mol.atom_dnm_flag[atom_index]        = mol2.atom_dnm_flag[i]; // JDB
            }

            renumber2[i] = atom_index;
            atom_index++;
        }
    }
 
    // Amend bond information from mol1 to combined_mol
    // bond_index serves as a place holder to determine the bond number of the new bond
    int bond_index = 0;
    string new_bond_type = "";
    for ( int i=0; i<mol1.num_bonds; i++ ){
        if ( mol1.bond_keep_flags[i] == true ){

            combined_mol.bonds_origin_atom[bond_index] = renumber1[mol1.bonds_origin_atom[i]];
            combined_mol.bonds_target_atom[bond_index] = renumber1[mol1.bonds_target_atom[i]];
            combined_mol.bond_types[bond_index]        = mol1.bond_types[i];
            combined_mol.bond_ring_flags[bond_index]   = mol1.bond_ring_flags[i];
            combined_mol.bond_keep_flags[bond_index]   = mol1.bond_keep_flags[i];
            bond_index++;
        }
    }

    // Now look through all the bonds of the seocnd parent
    for ( int i=0; i<mol2.num_bonds; i++ ){
        if ( mol2.bond_keep_flags[i] == true ){

            combined_mol.bonds_origin_atom[bond_index] = renumber2[mol2.bonds_origin_atom[i]];
            combined_mol.bonds_target_atom[bond_index] = renumber2[mol2.bonds_target_atom[i]];
            combined_mol.bond_types[bond_index]        = mol2.bond_types[i];
            combined_mol.bond_ring_flags[bond_index]   = mol2.bond_ring_flags[i];
            combined_mol.bond_keep_flags[bond_index]   = mol2.bond_keep_flags[i];
            bond_index++;

        }
    }

    // Determine the bond type for  the last bond which will be the connection between the two active halves
    new_bond_type = mol1.bond_types[bond_num1];  

    // This set of if-else statements will determine which atoms will be the origin and target atoms for the new bond
    if ( renumber1[mol1.bonds_origin_atom[bond_num1]] != 0 ){
        combined_mol.bonds_origin_atom[bond_index] = renumber1[mol1.bonds_origin_atom[bond_num1]];
    }
    else {
        combined_mol.bonds_origin_atom[bond_index] = renumber1[mol1.bonds_target_atom[bond_num1]];}

    
    if ( renumber2[mol2.bonds_origin_atom[bond_num2]] != 0 ){
        combined_mol.bonds_target_atom[bond_index] = renumber2[mol2.bonds_origin_atom[bond_num2]];
     }
     else {
         combined_mol.bonds_target_atom[bond_index] = renumber2[mol2.bonds_target_atom[bond_num2]];}
    

    // Set the bond type of that new bond here
    combined_mol.bond_types[bond_index] = new_bond_type;

    // Renumber the dummy and heavy atoms
    std::vector <std::pair <int, int> > renumbered_du_heavy;
    std::pair <int, int> tmp_du_heavy;

    // If there are entrie sin the pair vector
    if ( du_heavy.size() > 0 ){
       // For each pair
       for ( int i=0; i<du_heavy.size(); i++ ){
          // All of the dummy/heavy atoms should be on mol1
          // In renumbered, the old position is the index and the new position is the value
          tmp_du_heavy.first = renumber1[du_heavy[i].first];
          tmp_du_heavy.second = renumber1[du_heavy[i].second];
          renumbered_du_heavy.push_back(tmp_du_heavy);
       }
    }
 
    // Make a new mol_info_line that contains total number of atoms and bonds
    stringstream mol_info;
    mol_info <<" mol info " <<natoms <<"  " <<nbonds;
    combined_mol.mol_info_line = mol_info.str();
    mol_info.clear();
    // Update ring flags
    combined_mol.id_ring_atoms_bonds();
   

    atom_vec1.clear();
    atom_vec2.clear();
    bond_vec1.clear();
    bond_vec2.clear();
    renumber1.clear();
    renumber2.clear();
 
    //MARK - DELETE
    /*fstream fout_molecules;
    fout_molecules.open ( "zzz.replace_combine.mol2", fstream::out|fstream::app );
    fout_molecules <<combined_mol.current_data << endl;
    Write_Mol2(combined_mol, fout_molecules);
    fout_molecules.close();*/

    keep.push_back( combined_mol );

    return renumbered_du_heavy;
} // end GA_Recomb::replace_combine_mols()




// +++++++++++++++++++++++++++++++++++++++++
// This function will create the second child from the two parents around the same overlapping bond - called in breeding
void
GA_Recomb::switch_active_halves( DOCKMol & mol1, int mol1o, int mol1t, DOCKMol &  mol2, int mol2o, int mol2t, 
                                 std::vector <DOCKMol> & keep ) 
{
    Trace trace( "GA_Recomb::switch_active_halves()" );
    // Set all atom_active_flags to false
    for ( int i=0; i<mol1.num_atoms; i++ ){
        mol1.atom_active_flags[i]=false;
    }

    for ( int i=0; i<mol2.num_atoms; i++ ){
        mol2.atom_active_flags[i]=false;
    }

    // Set all bond_keep_flags to false
    for ( int i=0; i<mol1.num_bonds; i++ ){
        mol1.bond_keep_flags[i]=false;
    }

    for ( int i=0; i<mol2.num_bonds; i++ ){
        mol2.bond_keep_flags[i]=false;
    }
    //cout << "Entering switch activate" << endl;
    // Activate the opposite portion of the molecule than was previously activated
    activate_half_mol( mol1, mol1t, mol1o, mol2, mol2t, mol2o, false );
    //cout << "Entering switch recomb" << endl;
    // Return a DOCKMol of the new child 
    recomb_mols_xover( mol1, mol1t, mol1o, mol2, mol2t, mol2o, keep );

    return;

} // end GA_Recomb::switch_active_halves()




// +++++++++++++++++++++++++++++++++++++++++
// Comverts a mol to fragment class and initializes parameter for torenv check
Fragment
GA_Recomb::mol_to_frag( DOCKMol & mol )
{
    Trace trace( "GA_Recomb::mol_to_frag()" );
    Fragment frag;
    frag.mol = mol;
    frag.used = false;
    frag.tanimoto = 0.0;
    frag.last_ap_heavy = -1;
    frag.scaffolds_this_layer = 0;    
    frag.torenv_recheck_indices.clear();
    frag.size = 0;
    frag.ring_size = 0;
    frag.frag_growth_tree.clear();
    frag.aps_bonds_type.clear();
    frag.aps_cos.clear();
    frag.aps_dist.clear();

    return frag;
} // end GA_Rcomb::mol_to_frag()




// +++++++++++++++++++++++++++++++++++++++++
// Activate all of the atoms AND bonds in a given DOCKMol - called in sample_minimized_torsions
void
GA_Recomb::prepare_mol_torenv( DOCKMol & mol )
{
    Trace trace( "GA_Recomb::prepare_mol_torenv()" );
    // For torenv comparison, only the children origin and target atoms of the bond 
    // created in breeding should be active
    // Iterate through all atoms and set atom_active_flag to false
    for ( int i=0; i<mol.num_atoms; i++ ){
        mol.atom_active_flags[i] = false;
    }

    // Iterate through all bonds and set bond_active_flag to false
    for ( int i=0; i<mol.num_bonds; i++ ){
        mol.bond_active_flags[i] = true;
    }

    //Activate the origin and target atoms around the bond that was just involved in breeding
    mol.atom_active_flags[mol.bonds_origin_atom[mol.num_bonds-1]] =  mol.atom_active_flags[mol.bonds_target_atom[mol.num_bonds-1]] = true;
    
    return;

} // end GA_Recomb::prepare_mol_torenv()




// +++++++++++++++++++++++++++++++++++++++++
// Activate all of the atoms AND bonds in a given DOCKMol - called in sample_minimized_torsions
void
GA_Recomb::prepare_mut_torenv( DOCKMol & mol, int bond )
{
    Trace trace( "GA_Recomb::prepare_mut_torenv()" );
    // For torenv comparison, only the children origin and target atoms of the bond 
    // created in breeding should be active
    // Iterate through all atoms and set atom_active_flag to false
    for ( int i=0; i<mol.num_atoms; i++ ){
        mol.atom_active_flags[i] = false;
    }

    // Iterate through all bonds and set bond_active_flag to false
    for ( int i=0; i<mol.num_bonds; i++ ){
        mol.bond_active_flags[i] = true;
    }

    //Activate the origin and target atoms around the bond that was just involved in breeding
    mol.atom_active_flags[mol.bonds_origin_atom[bond]] =  mol.atom_active_flags[mol.bonds_target_atom[bond]] = true;
    
    return;

} // end GA_Recomb::prepare_mut_torenv()





// +++++++++++++++++++++++++++++++++++++++++
// Activate all of the atoms AND bonds in a given DOCKMol - called in sample_minimized_torsions
void
GA_Recomb::prepare_torenv_indices( Fragment & frag)
{
    Trace trace( "GA_Recomb::prepare_torenv_indices()" );
    // For bonds that are rotatable
    for (int i=0; i<frag.mol.num_bonds; i++){
        // if bond is rotor
        if (frag.mol.amber_bt_id[i] != -1){
           pair <int, int> temp_pair;
           temp_pair.second = frag.mol.bonds_origin_atom[i];
           temp_pair.first = frag.mol.bonds_target_atom[i];
           frag.torenv_recheck_indices.push_back(temp_pair);
        }
    }
    return;

} // end GA_Recomb::prepare_torenv_indices()




            //////////////////////////////////////////
            //**********MUTATION FUNCTIONS**********//
            /////////////////////////////////////////
            



// +++++++++++++++++++++++++++++++++++++++++
// Exhaustive mutation will be attempted until max_mut_cycles is reached  
// Only one mutation type will be perfomed on each molecule and only one mutant will be maintained
void
GA_Recomb::master_mut_exhaustive( std::vector <DOCKMol> & start, std::vector <DOCKMol> & save, std::vector <DOCKMol> & output, 
                                  bool parents, Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & c_typer, Orient & orient, int gen_num)
{
    Trace trace( "GA_Recomb::master_mut_exhaustive()" );
    // STEP 4: Mutations
    // Vector of unmanipulated molecules
    for (int i=0; i<start.size(); i++ ){
        save.push_back(start[i]);
    }

    //Initialize the vector to keep track of molecules that have been mutations
    // Clear vectors that hold mutation success rate information
    success_del.clear();
    success_add.clear();
    success_sub.clear();
    success_replace.clear();

    // In formation on number of cycles and random picks
    int mut_cycles = 0;
    int num_picks = 0;

    //int max_parents_mut = start.size();
    // Tmp object to hold unmanipulated molecule
    DOCKMol copy_mol;

    // While final size is less than expected size
    while (output.size() < start.size()){

       // For each molecule
       for ( int i=0; i<start.size(); i++ ){
          
          num_picks++;
          // Copy the molecule for safe keeping
          copy_molecule( copy_mol, start[i] );

          cout << endl << "Gen "<< gen_num+1 << ": Mutation Attempt "<< num_picks <<": Mutating Molecule #" << i << endl;

          // Save number of molecules that have been generated
          int prev_size = output.size();
    
          // Initialize segment size
          int num_segments = 0;

          // Prepare molecule for segment id
          segment_id( start[i], c_typer );
  
          // Perform mutation
          mutation_selection( start[i], output, score, simplex, c_typer, orient, gen_num, num_picks);
          //cout << "in master_mut_exhaustive output size: " << output.size() << endl; 
          //cout << "in master_mut_exhaustive start size: " << start.size() << endl;
            // Delete the molecule from the vector
          start.erase(start.begin()+i);
            
          // Only one offspring is kept from each parent regardless
          // of how many mutations or mutation rounds - LEP

          // If there was not a successful mutation
          if (output.size() == prev_size){
             // Add the molecule to the end of the vector
             start.push_back(copy_mol);
          }

          // Delete the copy
          copy_mol.clear_molecule();

          // If we have reached the expected size, then break
          if (output.size() == start.size()){break;
            cout << "we have reached the expected size - break" << endl;
            }
            //LEP wtf is this
       }
       // Increment num cycles
       mut_cycles++;

       if (mut_cycles > ga_max_mut_cycles){break;
        cout << "mut_cycles is greater than max mut cycles - break" << endl;
        }
        cout << "output size: " << output.size() << endl;
    }

    // If mutated parents
    if (parents == true){
       // Replace start with original saved molecules
       cout << "For this generation, parents' mutations: " << output.size() << endl;
       cout <<"Expected number of parent mutations: " << save.size() << ", mut cycles " << mut_cycles<<endl;

       //Print mutation success rates
       int dels = 0;
       for (int i=0; i<success_del.size(); i++){
           if (success_del[i] == 1){
              dels++;
           }
       }
       cout << setw(STRING_WIDTH) << "Successful parent deletions: " << dels << " / " << success_del.size() <<endl;

       int add = 0;
       for (int i=0; i<success_add.size(); i++){
           if (success_add[i] == 1){
              add++;
           }
       }
       //string add_make_pretty;
       //string add_make_pretty  = add.str() /*<< " / " << string(success_add.size()) << endl*/; 
       cout << setw(STRING_WIDTH) << "Successful parent additions: " << add<< " / " << success_add.size()  <<endl;
       //add_make_pretty.clear();
       int subsit = 0;
       for (int i=0; i<success_sub.size(); i++){
           if (success_sub[i] == 1){
              subsit++;
           }
       }
       cout << setw(STRING_WIDTH) << "Successful parent substitutions: " << subsit << " / " << success_sub.size() <<endl;

       int replace = 0;
       for (int i=0; i<success_replace.size(); i++){
           if (success_replace[i] == 1){
              replace++;
           }
       }
       cout << setw(STRING_WIDTH) << "Successful parent replacements: " << replace << " / " << success_replace.size() <<endl;
    }
    else {
       // Replace start with original saved molecules
       cout << "For this generation, offspring mutations: " << output.size() << endl;
       cout <<"Expected number of offspring mutations: " << save.size() << " mut cycles " << mut_cycles<<endl;

       //Print mutation success rates
       int dels = 0;
       for (int i=0; i<success_del.size(); i++){
           if (success_del[i] == 1){
              dels++;
           }
       }
       cout << setw(STRING_WIDTH)<< "Successful offspring deletions: " << dels << " / " << success_del.size() <<endl;

       int add = 0;
       for (int i=0; i<success_add.size(); i++){
           if (success_add[i] == 1){
              add++;
           }
       }
       cout <<setw(STRING_WIDTH)<< "Successful offspring additions: " << add << " / " << success_add.size() <<endl;

       int subsit = 0;
       for (int i=0; i<success_sub.size(); i++){
           if (success_sub[i] == 1){
              subsit++;
           }
       }
       cout <<setw(STRING_WIDTH)<< "Successful offspring substitutions: " << subsit << " / " << success_sub.size() <<endl;

       int replace = 0;
       for (int i=0; i<success_replace.size(); i++){
           if (success_replace[i] == 1){
              replace++;
           }
       }
       cout <<setw(STRING_WIDTH)<< "Successful offspring replacements: " << replace << " / " << success_replace.size() <<endl;
    }
    // Clear vectors that hold mutation success rate information
    success_del.clear();
    success_add.clear();
    success_sub.clear();
    success_replace.clear();
    return;
} // end GA_Recomb::master_mut_exhaustive()



// +++++++++++++++++++++++++++++++++++++++++
// Mutate molecules by randomly selecting a molecule until the expected number of mutants is generationed
// or until max_mut_cycles is reacted 
// Only one mutation type will be perfomed on each molecule and only one mutant will be maintained
void
GA_Recomb::master_mut_rand( std::vector <DOCKMol> & start, std::vector <DOCKMol> & save, std::vector <DOCKMol> & output, 
                       bool parents, Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & c_typer, Orient & orient, int gen_num)
{
    Trace trace( "GA_Recomb::master_mut_rand()" );
    // STEP 4: Mutations
    // User-defined a mutation rate so need the number of mutatants expected
    // If mutation rate define onlyd
    int expected_muts = 0;
    if (parents == true){
        if (start.size() == 1){
            expected_muts = 1;
        }else{
            expected_muts = ga_pmut_rate * start.size();
        }
       // If there are no expected mutations, exit function
       if (expected_muts == 0){
          cout <<"#### Your mutation rate (" << ga_pmut_rate << ") resulted in no mutations. Perhaps consider a different rate." <<endl;
          return;
       }
    } 
    else{ 
        if (start.size() == 1){
            expected_muts = 1;
        }else {
            expected_muts = ga_omut_rate * start.size();
        }
       // If there are no expected mutations, exit function
       if (expected_muts == 0){
          cout <<"#### Your mutation rate (" << ga_omut_rate << ") resulted in no mutations. Perhaps consider a different rate." <<endl;
          return;
       }
    }

    for (int i=0; i<start.size(); i++ ){
        save.push_back(start[i]);
    }

    // Complete the appropriate number of mutations

    // Initialize the vector to keep track of molecules that have been mutations
    // Clear vectors that hold mutation success rate information
    success_del.clear();
    success_add.clear();
    success_sub.clear();
    success_replace.clear();

    // In formation on number of cycles and random picks
    int num_picks = 0;
    int total_num_picks = 0;
    int mut_cycles = 0;
    DOCKMol copy_mol;

    while (output.size() < expected_muts){

       // Select a random molecule to mutation
       int index_mut = (rand() % start.size()); 

       // Increment the number of picks
       num_picks++;
       total_num_picks++; 
       // Copy the molecule for safe keeping
       copy_molecule( copy_mol, start[index_mut] );
       cout << endl << "Gen "<< gen_num+1 << ": Mutation Attempt "<< total_num_picks <<": Mutating Molecule #" << index_mut << endl;
       //cout << " MOLECULE #" << index_mut << endl;

       // Save number of molecules that have been generated
       int prev_size = output.size();
    
       // Initialize segment size
       int num_segments = 0;

       // Prepare molecule for segment id
       segment_id( start[index_mut], c_typer );
  
       // Preform mutation
       mutation_selection( start[index_mut], output, score, simplex, c_typer, orient, gen_num, total_num_picks );
   
       // Delete the molecule from the vector
       start.erase(start.begin()+index_mut);

       // If there was not a successful mutation
       if (output.size() == prev_size){
          // Add the molecule to the end of the vector
          start.push_back(copy_mol);
       }

       // Delete the copy
       copy_mol.clear_molecule();

       // If we have reached the expected size, then break
       if (output.size() == expected_muts){break;}
       if (verbose) cout <<"Number of expected " << expected_muts << " size " <<output.size() << " mut cycle " << mut_cycles<<endl;

       // When the number of random picks == the original size of the offspring, that is 1 cycle
       // Only a certain number of cycles can be completed 
       if ( num_picks == save.size() ){
          // Increment num cycles
          mut_cycles++;

          // Clear the num picks
          num_picks = 0;
       }
       if (mut_cycles > ga_max_mut_cycles){break;}
    }

    if (parents == true){
       cout <<"For this generation, parent mutations: " << output.size() << endl;
       cout <<"Expected number of parent mutations: " << expected_muts <<  " size " <<output.size() << " mut cycle " << mut_cycles<<endl;

       //Print mutation success rates
       int dels = 0;
       for (int i=0; i<success_del.size(); i++){
           if (success_del[i] == 1){
              dels++;
           }
       }
       cout << setw(STRING_WIDTH)<< "Successful parent deletions: " << dels << " / " << success_del.size() <<endl;

       int add = 0;
       for (int i=0; i<success_add.size(); i++){
           if (success_add[i] == 1){
              add++;
           }
       }
       cout << setw(STRING_WIDTH)<< "Successful parent additions: " << add << " / " << success_add.size() <<endl;

       int subsit = 0;
       for (int i=0; i<success_sub.size(); i++){
           if (success_sub[i] == 1){
              subsit++;
           }
       }
       cout << setw(STRING_WIDTH)<< "Successful parent substitutions: " << subsit << " / " << success_sub.size() <<endl;

       int replace = 0;
       for (int i=0; i<success_replace.size(); i++){
           if (success_replace[i] == 1){
              replace++;
           }
       }
       cout << setw(STRING_WIDTH)<< "Successful parent replacements: " << replace << " / " << success_replace.size() <<endl;
    }
    else{
       cout <<"For this generation, offspring mutations: " << output.size() << endl;
       cout <<"Expected number of offspring mutations: " << expected_muts <<  " size " <<output.size() << " mut cycle " << mut_cycles<<endl;

       //Print mutation success rates
       int dels = 0;
       for (int i=0; i<success_del.size(); i++){
           if (success_del[i] == 1){
              dels++;
           }
       }
       cout << setw(STRING_WIDTH)<< "Successful offspring deletions: " << dels << " / " << success_del.size() <<endl;

       int add = 0;
       for (int i=0; i<success_add.size(); i++){
           if (success_add[i] == 1){
              add++;
           }
       }
       cout << setw(STRING_WIDTH)<< "Successful offspring additions: " << add << " / " << success_add.size() <<endl;

       int subsit = 0;
       for (int i=0; i<success_sub.size(); i++){
           if (success_sub[i] == 1){
              subsit++;
           }
       }
       cout << setw(STRING_WIDTH)<< "Successful offspring substitutions: " << subsit << " / " << success_sub.size() <<endl;

       int replace = 0;
       for (int i=0; i<success_replace.size(); i++){
           if (success_replace[i] == 1){
              replace++;
           }
       }
       cout << setw(STRING_WIDTH)<< "Successful offspring replacements: " << replace << " / " << success_replace.size() <<endl;
    }

    // Clear vectors that hold mutation success rate information
    success_del.clear();
    success_add.clear();
    success_sub.clear();
    success_replace.clear();
    return;
} // end GA_Recomb::master_mut_rand()



// +++++++++++++++++++++++++++++++++++++++++
// Prepares molecule so that segments and intersegment bonds can be easily called for mutations
void
GA_Recomb::segment_id( DOCKMol & gen, AMBER_TYPER & typer )
{
    Trace trace( "GA_Recomb::segment_id()" );
    double start_time = time_seconds();
    // Prepare molecules with amber typer, etc
    activate_mol(gen);
    //cout << "Entering Amber typer" << endl;
    gen.prepare_molecule();
    
    // Identify the molecule segments using a function from AG
    AG_Conformer_Search c_ag;

    // Initialize necessary "user" input 
    c_ag.user_specified_anchor = false; 
    c_ag.limit_max_anchors = false;
    c_ag.anchor_size = 40;
    c_ag.max_anchor_num = 10;

    // Id the atoms and bonds in segments
    //cout << "Entering Amber typer" << endl;
    c_ag.prepare_molecule(gen);

    // Re-initialize num_segments 
    num_segments = 0;

    // Save number of segments
    for ( int i=0; i<gen.num_atoms; i++ ){
        if ( gen.atom_segment_ids[i] > num_segments ){
           num_segments = gen.atom_segment_ids[i];
        }
    }

    // Add 1 for proper C++ indexing
    num_segments = num_segments +1;
   
    // Save segment info for easy retrieval
    // Add empty tmp_segment for each segment in mol
    SEGMENTS tmp_segment;
    for ( int i=0; i<num_segments; i++ ){
        orig_segments.push_back(tmp_segment);
    }

    // Save the atoms associated with a segment
    for ( int i=0; i<gen.num_atoms; i++ ){
        for ( int j=0; j< num_segments; j++ ){
            if ( gen.atom_segment_ids[i] == j ){
                orig_segments[j].atoms.push_back(i);
                orig_segments[j].num_atoms++;
            }
        }
    }
  
    // List of all bonds in all segments
    INTVec  in_seg_bonds;
   
    // Save the bonds associated with a segment
    // For each bond in mol
    for ( int i=0; i<gen.num_bonds; i++ ){
        // For each segment
        for ( int j=0; j< num_segments; j++ ){
            // Associate the bond with its segment via segment id
            if ( gen.bond_segment_ids[i] == j ){
                // Save bond to segment bond list
                orig_segments[j].bonds.push_back(i);
                // Save num bonds
                orig_segments[j].num_bonds++;
                // Add bond to list of all bonds in segments
                in_seg_bonds.push_back(i);
            }
        }
    }

    // Identify intersegment bonds
    // for each bond in gen
    for ( int i=0; i<gen.num_bonds; i++ ){
        // If the bond is present in in_seg_bonds
        if ( find ( in_seg_bonds.begin(), in_seg_bonds.end(), i ) != in_seg_bonds.end() ){
        }
        // If no present in a segment, add to no_seg_bonds
        else{
            no_seg_bonds.push_back(i);
        }
    }

    // Save origin/target segment ids for intersegment bonds and update num attachment points per segment
    INTPair tmp;
    for( int i=0; i<no_seg_bonds.size(); i++){
       tmp.first = gen.atom_segment_ids[gen.bonds_origin_atom[no_seg_bonds[i]]];
       tmp.second = gen.atom_segment_ids[gen.bonds_target_atom[no_seg_bonds[i]]];
      
       orig_segments[tmp.first].num_aps++;
       orig_segments[tmp.second].num_aps++;
       bond_btwn_seg.push_back(tmp);
    }

    // Print to file
    // DELETE - REMOVE
    /*fstream fout_molecules;
    fout_molecules.open ( "zzz.gen.mol2", fstream::out|fstream::app );
    fout_molecules <<gen.current_data << endl;
    Write_Mol2(gen, fout_molecules);
    fout_molecules.close();*/

    // Clear vectors
    in_seg_bonds.clear();

    double stop_time = time_seconds();
    if (verbose) cout << "\n" "-----------------------------------" "\n";
    if (verbose) cout << " Elapsed time segment_ids:\t" << stop_time - start_time
         << " seconds for GA Mutations\n\n";

    return;
}//end GA_Recomb:segment_id()





// +++++++++++++++++++++++++++++++++++++++++
// This function will randomly select the mutation type, identify bond (& atoms) for mutation,
// and call the appropriate mutation functions (delete, replace, insert, or substitute)
void
GA_Recomb::mutation_selection( DOCKMol & gen, std::vector <DOCKMol> & mfinal, Master_Score & score, Simplex_Minimizer & simplex, 
                               AMBER_TYPER & typer, Orient & orient, int gen_num, int mut_attempt )
{
    Trace trace( "GA_Recomb::mutation_selection()" );
    double start_time = time_seconds();
   // cout << "Picking a mutation type" << endl;
    // STEP 1: Randomly select a mutation type (from 4 options)
    // MARK : need other mutation selection types and to make a user parameter to choose one
    int mutation_type = (rand() % (int)(4));
    cout << "Random number - mutation_type (mutation_selection): " << mutation_type << endl; 
    //cout << "Mutation type " <<mutation_type <<endl;
    if (mutation_type == 0){
        cout << "Trying mutation type: deletion" << endl;
    }
    if (mutation_type == 1){
        cout << "Trying mutation type: addition" << endl;
    }
    
    if (mutation_type == 2){
        cout << "Trying mutation type: substitution" << endl;
    }
    if (mutation_type == 3){
        cout << "Trying mutation type: replacement" << endl;
    }
    //mutation_type = 1;


    //-------------------------------------------------------------------------
    // STEP 2: Pick a scaffolds, linkers, and sidechains of the molecule
    // (1) Identify the number of scaffolds, linkers and sidechains
    INTVec sidechain_ids;
    INTVec linker_ids;
    INTVec scaffold_ids;
    INTVec rigid_ids;
    //for (int i=0; i<num_segments; i++){
    //    for (int atom_num=0; atom_num<orig_segments[i].num_atoms; ++atom_num){
    //        //checks all atoms to see if they are mutatable
    //        cout << "START_MUT: "<< gen.atom_dnm_flag[orig_segments[i].atoms[atom_num]] << endl;
    //    }
    //}
//JDB start for DNM
    if (dnm_enabled){
        seg_exclude_indices.clear(); //clear from last molecule
        bool dnm_encountered=false;
        
        // For all segments in mol
        for (int i=0; i<num_segments; i++){
          bool seg_exclude=false; //for exclusion for all but addition
          bool seg_exclude_add=false; //for exclusion for addition
          int dna_counter = 0;


          //runs through all the atoms in the segment
          for (int atom_num=0; atom_num<orig_segments[i].num_atoms; ++atom_num){
            //checks all atoms to see if they are mutatable
            if (gen.atom_dnm_flag[orig_segments[i].atoms[atom_num]]) {
              //if true, add to counter and exclude for all mutations but addition
              seg_exclude=true;
              dna_counter+=1;
            }
          }
            //if all atoms are marked as not mutatable, exclude for addition
          if (orig_segments[i].num_atoms == dna_counter){
            seg_exclude_add=true;
          } 

          //if excluded for addition and addition selected
          if ((seg_exclude_add == true) && (mutation_type == 1)){
            cout << "Will not select segment " << i << " as all atoms are labeled for no mutation."<< endl;
            seg_exclude_indices.push_back(i);
            seg_exclude_add = false;
            seg_exclude = false;

          //if excluded for all but addition and addition is selected
          } else if ((seg_exclude == true) && (mutation_type == 1)) { 
            cout << "Segment " << i << " labeled for no mutation, but addition still possible."<< endl;

          //if excluded for all but addition and something other than addition chosen
          } else if ((seg_exclude == true) && (mutation_type != 1)) {
            cout << "Segment " << i << " labeled for no mutation." << endl;
            seg_exclude_indices.push_back(i);
            seg_exclude=false;
          }
        }
    } //JDB end for DNM


    // For all segments in mol
    for (int i=0; i<num_segments; i++){
        // If the segment has 1 attachment point
        if (orig_segments[i].num_aps == 1){
           // Save to vector
          sidechain_ids.push_back(i);
          terminal_seg.push_back(i); // This is the global version needed to speed up deletions
        }
        else if ( orig_segments[i].num_aps == 2 ){
           linker_ids.push_back(i);   
           /*  
           //JDB - Atom Names per segment debugging
            for (int x=0; x<orig_segments[i].num_atoms; x++){
                cout << gen.atom_names[orig_segments[i].atoms[x]] << ", ";
                cout << gen.x[orig_segments[i].atoms[x]] << ", " 
                     << gen.y[orig_segments[i].atoms[x]] << ", " 
                     << gen.z[orig_segments[i].atoms[x]] << endl;
            }
            cout << orig_segments[i].num_aps;
            cout << endl << endl;      */ 
        }
        else if ( orig_segments[i].num_aps >= 3 ){
           scaffold_ids.push_back(i);
        }
        // If there are no attachment points
        else {
           rigid_ids.push_back(i);   
        }
    }
    //JDB
    //if all segment ids are empty, then the molecule is invalid for any mutations. Program will exit with error.
    if ( (sidechain_ids.size() == 0) && (linker_ids.size() == 0) && (scaffold_ids.size() == 0) && \
      (terminal_seg.size()) == 0 && (rigid_ids.size() == 0) ){
        cout << "No segments were kept as valid. Possible all segments were selected as DNM." << endl;
        exit(1);
    }
    // (2) Pick a random scaffold, sidechain and linker for mutation 
    int segment_types = 0;
    int final_segment_type = 0;
    // If the molecule is not rigid
    if ( (rigid_ids.size() == 0) ){
       // If there are multiple segment types
       if ( (sidechain_ids.size() > 0) && (linker_ids.size() > 0) && (scaffold_ids.size() > 0) ){
          segment_types = 3;
          // Randomly select a fragment type
          final_segment_type = (rand() % (int)(segment_types));
          //cout << "Random number - segment_type (mutation_selection): " << final_segment_type << endl;
       }
       // If only the sidechains and linkers are present
       else if ( (sidechain_ids.size() > 0) && (linker_ids.size() > 0) ){
          segment_types = 2; 
          // Randomly select a fragment type
          final_segment_type = (rand() % (int)(segment_types));
          //cout << "Random number - segment_type (mutation_selection): " << final_segment_type << endl;  
       }
       // If only if sidechains and scaffolds
       else if ( (sidechain_ids.size() > 0) && (scaffold_ids.size() > 0) ){
           segment_types = 2;
           final_segment_type = (rand() % (int)(segment_types));
           //cout << "Random number - segment_type (mutation_selection): " << final_segment_type << endl;
           // Here, option 0 = sidechains 7 option 1 = scaffolds
           if ( final_segment_type == 1 ){ final_segment_type = 2; 
                //cout << "changing segment type 1 to 2 because no linkers but there are scaffolds"<< endl; 
           }
       }
       // If there are only sidechains
       else { final_segment_type = 0; 
        cout << "there are only sidechains" << endl;
        }
       if (final_segment_type == 0){
          cout << "Segment Type for Mutation: sidechains" << endl;
       }
       if (final_segment_type == 1){
          cout << "Segment Type for Mutation: linkers" << endl;
       }
       if (final_segment_type == 2){
          cout << "Segment Type for Mutation: scaffolds" << endl;
       }
       if (final_segment_type == 3){
          cout << "Segment Type for Mutation: rigid" << endl;
       }

       //cout << "Segment Type for Mutation: " << final_segment_type << endl;
       /*cout <<"Molecular framework:  scaffolds " << scaffold_ids.size() 
            << ", linkers " << linker_ids.size() << ", sidechain " << sidechain_ids.size() 
            << ", rigid " << rigid_ids.size() <<endl;*/
    }    
    // If rigid, do an addition
    else{
       mutation_type = 1;
       cout << "Mutation type changed to addition because molecule is rigid" << endl;
       cout << " scaffolds " << scaffold_ids.size() << " linkers " << linker_ids.size() << " sidechain " << sidechain_ids.size() << " rigid " << rigid_ids.size() <<endl;
    }


    // (3) If not rigid, pick a segment from the available type
    int index_seg = -1;
    int final_seg = -1;
    if ( rigid_ids.size() == 0 ){
      //JDB DNM functionality
      bool dna_selected = false;
      int number_reselect = 10;
      int select_counter = 0;

      // If type = 0, select a sidechain
      if ( final_segment_type == 0 ){
        //cout << "type = 0, select sidechain" << endl;
        index_seg = (rand() % (int)(sidechain_ids.size()));
        //cout << "Random number - index_seg (mutation_selection): " << index_seg << endl;
        //JDB check chosen segment against all DNM excluded segments
        if (dnm_enabled) {
          for (int x=0; x<seg_exclude_indices.size(); ++x){
            //if the the included segment matches any of the excluded, break loop
            if (seg_exclude_indices[x] == sidechain_ids[index_seg]){
              dna_selected = true;
              //cout << "DNA TRUE 1 " << endl;//DEBUG 
              break;
            }
          }
        }
        // JDB if segment is DNM, and we haven't reached hard coded reselection threshhold
        while ((dna_selected == true) && (select_counter <= number_reselect)){
          dna_selected = false;
          //select another segment
          index_seg = (rand() % (int)(sidechain_ids.size()));

          //and check again
          if (dnm_enabled) {
            for (int x=0; x<seg_exclude_indices.size(); ++x){
              if (seg_exclude_indices[x] == sidechain_ids[index_seg]){
                dna_selected = true;
                //cout << "DNA TRUE 2 " << endl; //DEBUG
                break;
            
              }
            }
          }
          select_counter++;
        } // DNM
        final_seg = sidechain_ids[index_seg];
        //cout << "num sidechains " << sidechain_ids.size() << " index " << index_seg 
        //     << " working segment " << final_seg <<endl;
 
        // If the mutation type was 3 (replacement), change the mutation type to 2 (subsit)
        if ( mutation_type == 3 ){
           mutation_type = 2;
           cout << " Changed mutation type from replacement to substitution" <<endl;
        }
      }
      // If type = 1, select a linker
      if ( final_segment_type == 1 ){
        //cout << "final segment type = 1, select a linker" << endl;
        index_seg = (rand() % (int)(linker_ids.size()));
        cout << "Random number - index_seg (mutation_selection): " << index_seg << endl;
        //JDB check chosen segment against all DNM excluded segments
        if (dnm_enabled) {
          for (int x=0; x<seg_exclude_indices.size(); ++x){
            if (seg_exclude_indices[x] == linker_ids[index_seg]){
              dna_selected = true;
              //cout << "DNA TRUE 3 " << endl; //DEBUG
              break; //break if DNM
            }
          }
        // JDB if segment is DNM, and we haven't reached hard coded reselection threshhold
          while ((dna_selected == true) && (select_counter <= number_reselect)){
            dna_selected = false;
            //select another segment
            index_seg = (rand() % (int)(linker_ids.size()));

            //and try again
            for (int x=0; x<seg_exclude_indices.size(); ++x){
              if (seg_exclude_indices[x] == linker_ids[index_seg]){
                dna_selected = true;
                //cout << "DNA TRUE 4 " << endl; //DEBUG
                break; //break if DNM
              }
            }
            select_counter++;
          }  //DNM
        }

        final_seg = linker_ids[index_seg];
        //cout << "num linkers " << linker_ids.size() << " index " << index_seg 
        //     << " working segment " << final_seg <<endl;
      }
      // If type = 2, select a scaffolds
      if ( final_segment_type == 2 ){
        //cout << "final segment type = 2, select a scaffold" << endl;
        index_seg = (rand() % (int)(scaffold_ids.size()));
        cout << "Random number - index_seg (mutation_selection): " << index_seg << endl;

        if (dnm_enabled) {
          for (int x=0; x<seg_exclude_indices.size(); ++x){
            if (seg_exclude_indices[x] == scaffold_ids[index_seg]){
              dna_selected = true;
              //cout << "DNA TRUE 5 " << endl; //DEBUG
              break; //break if DNM
            }
          }
          // JDB if segment is DNA, and we haven't reached hard coded reselection threshhold
          while ((dna_selected == true) && (select_counter <= number_reselect)){
            dna_selected = false;
            //select another segment
            index_seg = (rand() % (int)(scaffold_ids.size()));

            //and try again
            for (int x=0; x<seg_exclude_indices.size(); ++x){
              if (seg_exclude_indices[x] == scaffold_ids[index_seg]){
                dna_selected = true;
                //cout << "DNA TRUE 6 " << endl; //DEBUG
                break; //break if DNA
              }
            }
            select_counter++;
          }
        }

        final_seg = scaffold_ids[index_seg];
        //cout << "num scaffold " << scaffold_ids.size() << " index " << index_seg 
        //     << " working segment " << final_seg <<endl;
      }

       // JDB
       // if mass reselection still yields a DNM fragment, perform cleanup and return function
       // before mutation performed
      if (dnm_enabled) {
        if (dna_selected == true) {
          cout << "Final segment for this attempt is not allowed for mutation type by user input. Moving to next attempt." << endl << endl;


          orig_segments.clear();
          bond_btwn_seg.clear();
          no_seg_bonds.clear();
          mutants.clear();

          sidechain_ids.clear();
          linker_ids.clear();
          scaffold_ids.clear();
          rigid_ids.clear();

          return;
         }
      }
    } else { // If rigid update index
      final_seg = rigid_ids[0];
    }

    //-------------------------------------------------------------------------
    // STEP 3: Map the aps bond and atom vec for mutations
    INTVec working_bonds; // all intersegment bonds associated with selected segment
    INTVec atom_start;    // origin and target atoms, respectively, associated with bond
    INTVec atom_term;

    // While there are attachment points for the working segment
    int i = 0; 
    while ( i < orig_segments[final_seg].num_aps ){
       // For each intersegment bond (and pair of segments it is between)
       for ( int j=0; j<bond_btwn_seg.size(); j++){
           // If the working segment is associated with the intersegment bond
           if ( final_seg == bond_btwn_seg[j].first 
             || final_seg == bond_btwn_seg[j].second ){
              // Add bond number to working list
              working_bonds.push_back(no_seg_bonds[j]);
              i++;
            }
       }
    }
    //JDB outputs the atom names of segments that will be mutated this round
    cout << "Atoms of segment for mutation: ";
    for (int x = 0; x < orig_segments[final_seg].num_atoms; x++) {
      cout << gen.atom_names[orig_segments[final_seg].atoms[x]];

      if (x < orig_segments[final_seg].num_atoms - 1) {
        cout << ",";
      } else {
        cout << endl;
        continue;
      }
    }
    // Find the atoms associated with working segment and its adjacent segment
    for ( int i=0; i<working_bonds.size(); i++ ){
        if ( gen.atom_segment_ids[gen.bonds_origin_atom[working_bonds[i]]] == final_seg ){
           atom_start.push_back(gen.bonds_origin_atom[working_bonds[i]]);
           atom_term.push_back(gen.bonds_target_atom[working_bonds[i]]);
        }
        else if ( gen.atom_segment_ids[gen.bonds_target_atom[working_bonds[i]]] == final_seg ){
            atom_start.push_back(gen.bonds_target_atom[working_bonds[i]]);
            atom_term.push_back(gen.bonds_origin_atom[working_bonds[i]]);
        }
    }

    //-------------------------------------------------------------------------
    // STEP 4: Call the appropriate mutation functions
    // (1) To delete a portion of the molecule and replace with a H - yields valid molecule, when applicable
    if ( (mutation_type == 0) && (num_segments > 2) && (ga_mutate_deletion)){ 
        cout << "Calling Mutation Function: Deletion" << endl;
        deletion_mutation(gen, working_bonds, atom_start, atom_term, false, score, simplex, typer );
        if (dnm_enabled) {
          if (dnm_encountered == true){
              cout << "DNM Encountered in Deletion step. Exiting substitution." << endl;
              atom_start.clear();
              atom_term.clear();
              working_bonds.clear();
              orig_segments.clear();
              bond_btwn_seg.clear();
              no_seg_bonds.clear();
              mutants.clear();

              sidechain_ids.clear();
              linker_ids.clear();
              scaffold_ids.clear();
              rigid_ids.clear();
              return;
          }
        }
        //if (mutants.size() == 0){
        //    cout << "no molecules made from deletion" << endl;
        //}
        if (mutants.size()>0){
            hard_filter(mutants); 
        }
        if ( mutants.size() > 0 ){ 
            if (mutants[0].energy == "parent"){
                ostringstream parent_field;
                parent_field << gen.title << "-d";
                mutants[0].energy = parent_field.str();
            }else{
                ostringstream xover_parent_field;
                xover_parent_field << mutants[0].energy << "-d";
                mutants[0].energy = xover_parent_field.str();
            }
            success_del.push_back(1); 
            //cout << "Lauren:" << gen.title << endl;
        } 
        else { success_del.push_back(0); }
    }//else{
     //   cout << "pushed to mutants" << endl;
     //   mutants.push_back(gen);
   // }
  
    // (2) To add to the terminal of a molecule - yields a valid, minimized, and scored molecule every time
    else if (( mutation_type == 1) && (ga_mutate_addition == true) ){
       cout << "Calling Mutation Function: Addition" << endl;
       ostringstream parent_field;
        if (gen.energy == "parent"){
           parent_field << gen.title << "-a";
        } else {
            parent_field << gen.energy << "-a";
        }
       additions(gen, final_seg, score, simplex, typer, orient, false, gen_num, mut_attempt);
       if (mutants.size()>0){
            hard_filter(mutants); 
       }
       if ( mutants.size() > 0 ){
           //if (mutants[0].energy == "parent"){
             //   ostringstream parent_field;
               // parent_field << gen.title << "-a";
                mutants[0].energy = parent_field.str();
           // }else{
             //   ostringstream xover_parent_field;
              //  xover_parent_field << mutants[0].energy << "-a";
              //  mutants[0].energy = xover_parent_field.str();
           // }
            //cout << "Lauren:" << gen.title << endl;
          for (int i=0; i<mutants[0].num_atoms; i++){
              // Reset Atom residue numbers to prevent issues with DN minimization - rot bond initialization
              ostringstream res_num;
              res_num << "1";
              mutants[0].atom_residue_numbers[i] = res_num.str(); 
          }
          success_add.push_back(1); 
       } else{ success_add.push_back(0); }
    }//else{
     //   cout << "pushed to mutants" << endl;
     //   mutants.push_back(gen);
   // }      

    // (3) To replace any portion of the molecule - yields a valid, minimized, and scored molecule (only tries twice)
    else if (( mutation_type == 2) && (ga_mutate_substitution == true) ){
      cout << "Calling Mutation Function: Substitution" << endl;
      ostringstream parent_field;

      if (gen.energy == "parent"){
        parent_field << gen.title << "-s";
      } else {
        parent_field << gen.energy << "-s";
      }

      deletion_mutation(gen, working_bonds, atom_start, atom_term, true, score, simplex, typer );
      if (dnm_enabled) {
        if (dnm_encountered == true){
          cout << "DNM Encountered in Deletion step. Exiting substitution." << endl;
          atom_start.clear();
          atom_term.clear();
          working_bonds.clear();
          orig_segments.clear();
          bond_btwn_seg.clear();
          no_seg_bonds.clear();
          mutants.clear();

          sidechain_ids.clear();
          linker_ids.clear();
          scaffold_ids.clear();
          rigid_ids.clear();
          return;
        }
      }
       // Note this number is not entirely accurate should get number of deletions
       //ostringstream parent_field;
       //parent_field << gen.title << "-s";

      // JDB - necessary change for DNM to check for mutant size before addition - otherwise it can
      // pass an empty vector and seg fault
      if (mutants.size() > 0) {
        additions(mutants[0], final_seg, score, simplex, typer, orient, true, gen_num, mut_attempt);
      }
       
       if (mutants.size()>0){
            hard_filter(mutants); 
       }
       if ( mutants.size() > 0 ){
                //parent_field << gen.title << "-s";
                mutants[0].energy = parent_field.str();
                //cout << "Lauren:" << gen.title << endl;
                //ostringstream xover_parent_field;
                //xover_parent_field << mutants[0].energy << "-s";
                //mutants[0].energy = xover_parent_field.str();
            for (int i=0; i<mutants[0].num_atoms; i++){
                // Reset Atom residue numbers to prevent issues with DN minimization - rot bond initialization
                ostringstream res_num;
                res_num << "1";
                mutants[0].atom_residue_numbers[i] = res_num.str(); 
            }
            success_sub.push_back(1); 
       } else { success_sub.push_back(0); }
    }//else{
      //  cout << "pushed to mutants" << endl;
      //  mutants.push_back(gen);
   // }

    // (4) Replaces one fragment from a molecule, when applicable
    else if (( mutation_type == 3 ) && (ga_mutate_replacement == true) ){
       cout << "Calling Mutation Function: Replacement" << endl; 
       replacement( gen, final_seg, working_bonds, atom_start, atom_term, score, simplex, typer, orient ); 
        
       if ( mutants.size() > 0 ){
           if (mutants[0].energy == "parent"){
                ostringstream parent_field;
                parent_field << gen.title << "-r";
                mutants[0].energy = parent_field.str();
            } else {
                ostringstream xover_parent_field;
                xover_parent_field << mutants[0].energy << "-r";
                mutants[0].energy = xover_parent_field.str();
            }
          success_replace.push_back(1);
       }
       else{ success_replace.push_back(0); }
    }else{
        //cout << "pushed to mutants" << endl;
        mutants.push_back(gen);
    }
    cout << "Offspring size: " << mutants.size() << ", only single best scoring offspring retained" << endl;
    
    if ( mutants.size() > 0 ){
        // If the molecules have not been minimized, then pass to minimization function
        //cout << "mutants size " << mutants.size() << " mutants score " << mutants[0].current_score << endl;
        if (dnm_enabled) {
            for (int mutant_ind=0; mutant_ind<mutants.size(); mutant_ind++) {
                set_DNM_bools( mutants[mutant_ind] );
            }
        }  
        //DNM Debugging JDB JDB JDB
        //for(int i=0; i<mutants[0].num_atoms; i++){
        //    cout << "POST MUTATION BOOL: " << mutants[0].atom_dnm_flag[i] << endl;
        //}
        mfinal.push_back(mutants[0]);     
    }

    //-------------------------------------------------------------------------
    // STEP 4: Print to file
    // DELETE - REMOVE
    //if (mutants.size() == 1){
    /*for (int i=0; i<mutants.size(); i++){
       fstream fout_molecules;
       fout_molecules.open ( "zzz.gen.mol2", fstream::out|fstream::app );
       fout_molecules <<mutants[i].current_data << endl;
       Write_Mol2(mutants[i], fout_molecules);
       fout_molecules.close();
    }*/

    // Clear all vectors
    atom_start.clear();
    atom_term.clear();
    working_bonds.clear();
    orig_segments.clear();
    bond_btwn_seg.clear();
    no_seg_bonds.clear();
    mutants.clear();

    sidechain_ids.clear();
    linker_ids.clear();
    scaffold_ids.clear();
    rigid_ids.clear();

    double stop_time = time_seconds();
    if (verbose) cout << "\n" "-----------------------------------" "\n";
    if (verbose) cout << " Elapsed time mutation_selection:\t" << stop_time - start_time
         << " seconds for GA Mutations\n\n";
    return;
}//end GA_Recomb::mutation_selection()




// +++++++++++++++++++++++++++++++++++++++++
// This function inactivates a portion of the molecule and replaces the inactive section with a H 
void
GA_Recomb::deletion_mutation( DOCKMol & gen, INTVec bonds, INTVec origV, INTVec termV, bool subst, Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & typer )
{
    Trace trace( "GA_Recomb::deleteion_mutation()" );
    double start_time = time_seconds();

    // STEP 1: Select 1 bond for deletion 
    // Make an empty DOCKMol object to serve as a place holder for activation
    DOCKMol empty;

    // Bond index of bond for deletion
    int sel_bond = 0;  // Bond number 
    int orig = 0;      // origin atom
    int term = 0;      // target atom

    // Id the bond that will be used
    // MARK : If the segment has two bonds, select a bond a random?
    if ( bonds.size() >> 1 ){
 
        // Randomly select an index location in bonds
        int index = rand() % (int)(bonds.size());

        // Populate the appropriate variables
        sel_bond = bonds[index];
        orig = origV[index];
        term = termV[index];       
    }
    // If there is only 1 bond 
    else{
        sel_bond = bonds[0];
        orig = origV[0];
        term = termV[0];
    }
  

    // STEP 2: Inactivate the appropriate portion of the molecule
    // False if at end of mol (terminal seg)
    bool next_segment = true; 

    // (1)  If the bond is associated with a terminal segment 
    for ( int i=0; i<terminal_seg.size(); i++ ){
        if ( gen.atom_segment_ids[orig] == terminal_seg[i] ){
          // Count number of atoms from term to end of molecule
          float term_term = orig_segments[ terminal_seg[i] ].num_atoms;

          // Save the remaining number of atoms
          activate_half_mol( gen, term, orig, empty, -1, -1, false);
          if (dnm_enabled) {
            for ( int i=0; i<gen.num_atoms; i++ ){
              if ( gen.atom_active_flags[i] == false ){
                if ( gen.atom_dnm_flag[i] ){ //DNM Check
                  //cout << "false_Name: " << gen.atom_names[i] << endl;
                  //cout << "DNA10 Found in Deletion Segment." << endl;
                  dnm_encountered=true;
                  return;
                }
              }
            }
          }

          add_H(gen,term, orig, subst, score, simplex, typer);
          next_segment = false;
          num_segs_removed = 1;
          if (!subst){
            ostringstream new_energy;
            new_energy << gen.energy << "-del";
            gen.energy    = new_energy.str();
          }
           double stop_time = time_seconds();
           if (verbose) cout << "\n" "-----------------------------------" "\n";
           if (verbose) cout << " Elapsed time deletions:\t" << stop_time - start_time
                << " seconds for GA Mutations\n\n";
           return;
        } 
    }

    //-------------------------------------------------------------------------
    // (2) If non-terminal, inactivate the side with fewest atoms
    // Initialize variables
    INTVec term_segs;            // Segment numbers on term atom side of molecule
    int start = term;            // Atom on "target" side of bond
    int prev = orig;             // Previous target atom
    float term_atoms = 0.000;    // Number of atoms on term atom side of molecule
    int orig_atoms = 0;          // Number of atoms on orig atom side
   
    // Add the term atom (target of bond) to list 
    term_segs.push_back( gen.atom_segment_ids[start] );

    // While at a non-terminal segment
    while ( next_segment ){
    
       // Variable to break out of loop once ided appropriate segment pair
       bool found = false;

       //------------------------------------------------------------------------
       // (A) If the next segment has 2 attachment points
       if ( orig_segments[gen.atom_segment_ids[start]].num_aps == 2 ) {

           // (a) If start segment is first in bond_btwn_seg list
           // Find the other segment it is attached to
           for ( int i=0; i<bond_btwn_seg.size(); i++ ){
               // ID bond segment pairs associated with start
               if ( bond_btwn_seg[i].first == gen.atom_segment_ids[start] ){
                   // If you have not identified the previous segment
                   if ( bond_btwn_seg[i].second != gen.atom_segment_ids[prev] ){
                       // Save segment id to list
                       term_segs.push_back(bond_btwn_seg[i].second);
                       // Save the appropriate atom associated with segment 
                       if ( gen.atom_segment_ids[gen.bonds_origin_atom[no_seg_bonds[i]]] == bond_btwn_seg[i].second ){
                           // Update variables
                           prev = start;
                           start = gen.bonds_origin_atom[ no_seg_bonds[i] ];
                           found = true;
                       }
                       else{ 
                           prev = start;
                           start = gen.bonds_target_atom[ no_seg_bonds[i] ];
                           found = true;
                      } 
                   }
               } if (found) { break;} // breaks out of for loop

               // (b) If start segment is second
               if ( bond_btwn_seg[i].second == gen.atom_segment_ids[start]){
                    // If you have not identified the previous segment
                    if (bond_btwn_seg[i].first != gen.atom_segment_ids[prev] ){
                       // Save segment id to list
                       term_segs.push_back(bond_btwn_seg[i].first);
                       // Save the appropriate atom associated with segment 
                       if ( gen.atom_segment_ids[gen.bonds_origin_atom[ no_seg_bonds[i]] ] == bond_btwn_seg[i].first){
                           prev = start;
                           start = gen.bonds_origin_atom[ no_seg_bonds[i] ];
                           found = true;
                       }
                       else{ 
                           prev = start;
                           start = gen.bonds_target_atom[ no_seg_bonds[i] ];
                           found = true;
                       }
                    }
               }  if (found) { break;} // breaks out of for 
           } 
       } // end if two attachment points

       
       //------------------------------------------------------------------------
       // (B) If 1 aps, ie terminal segment 
       // This is the last segment on the molecule in this direction, 
       // so it will only be attached to the previous segment
       if ( orig_segments[ gen.atom_segment_ids[start] ].num_aps == 1 ) {
           num_segs_removed = term_segs.size();

           next_segment = false;

          // Count number of atoms from term to end of molecule
          for ( int i=0; i<term_segs.size(); i++ ){
              term_atoms = orig_segments[ term_segs[i] ].num_atoms + term_atoms;
          }       

          // (a) Remove term side of molecule if fits parameters
          if ( double(term_atoms/gen.num_atoms) <= 0.5) {

             // Active/ inactivate the appropriate side ofthe molecule
             activate_half_mol(gen, orig, term, empty, -1, -1, false); 
             for ( int i=0; i<gen.num_atoms; i++ ){
                if ( gen.atom_active_flags[i] == false ){
                  if ( gen.atom_dnm_flag[i] ){
                    //cout << "false_Name: " << gen.atom_names[i] << endl;
                    //cout << "DNA1 Found in Deletion Segment." << endl;
                    dnm_encountered=true;
                      return;
                    }
                }
             }
             // Add H to complete molecule
             add_H(gen, orig, term, subst, score, simplex, typer);
          }
          // (b) Switch sides if necessary
          else{
            num_segs_removed = num_segments - term_segs.size();

            activate_half_mol(gen, term, orig, empty, -1, -1, false);
            if (dnm_enabled) {
              for ( int i=0; i<gen.num_atoms; i++ ){
                if ( gen.atom_active_flags[i] == false ){
                  if ( gen.atom_dnm_flag[i] ){
                    //cout << "false_Name: " << gen.atom_names[i] << endl;
                    //cout << "DNA2 Found in Deletion Segment." << endl;
                    dnm_encountered=true;
                    return;
                  }
                }
              } 
            }
            add_H(gen, term, orig, subst, score, simplex, typer);
          }
          double stop_time = time_seconds();
          if (verbose) cout << "\n" "-----------------------------------" "\n";
          if (verbose) cout << " Elapsed time deletions:\t" << stop_time - start_time
               << " seconds for GA Mutations\n\n";
          return;

       }// end if one attachment point

     
       //------------------------------------------------------------------------
       // (C) If at a scaffold, stop moving and delete that arm of the scaffold 
       if ( orig_segments[gen.atom_segment_ids[start]].num_aps > 2 ) {
           num_segs_removed = term_segs.size();

           // Update bool
           next_segment = false;

           // Count number of atoms on term arm of scaffold and remainder of molecule
           for ( int i=0; i<term_segs.size(); i++ ){
               term_atoms = orig_segments[term_segs[i]].num_atoms + term_atoms;
           }       
           orig_atoms = gen.num_atoms - term_atoms;

           // Save initial size before deactivation for check
           double orig_size = gen.num_atoms;

           // (a) Remove term side of molecule if fits parameters
           if ( double(term_atoms/gen.num_atoms) <= 0.5) {
 
              // MARK : NEED to add a bool here so that this will not be called for insertion( if else)
              // Active/ inactivate the appropriate side ofthe molecule
              activate_half_mol(gen, term, orig, empty, -1, -1, false); 
              if (dnm_enabled) {
              for ( int i=0; i<gen.num_atoms; i++ ){
                  if ( gen.atom_active_flags[i] == false ){
                    if ( gen.atom_dnm_flag[i] ){
                      //cout << "false_Name: " << gen.atom_names[i] << endl;
                      //cout << "DNA3 Found in Deletion Segment." << endl;
                      dnm_encountered=true;
                      return;
                    }
                  }
                }
              }
              // DELETE - REMOVE
              //fstream fout_molecules;
              //fout_molecules.open ( "zzz.scaffold_del.mol2", fstream::out|fstream::app );
              //fout_molecules <<gen.current_data << endl;
              //Write_Mol2(gen, fout_molecules);
              //fout_molecules.close();
 
              double mol_size =0.0;
              for ( int i=0; i<gen.num_atoms; i++){
                  if (gen.atom_active_flags[i] == true){
                      mol_size++;         
                  }
              }

              // If the size is greater than or equal to 50% or substitution, keep
              if ((subst == true) || ((mol_size/orig_size) >= 0.5)){
                 // Add H to complete molecule
                 add_H(gen, term, orig, subst, score, simplex, typer);
              }

              // Else switch sides
              else if ( (subst == false) && ((mol_size/orig_size) <= 0.5)) {
                //Switch side
                activate_half_mol(gen, orig, term, empty, -1, -1, false); 
                if (dnm_enabled) {
                  for ( int i=0; i<gen.num_atoms; i++ ){
                    if ( gen.atom_active_flags[i] == false ){
                      if ( gen.atom_dnm_flag ){
                        //cout << "false_Name: " << gen.atom_names[i] << endl;
                        //cout << "DNA4 Found in Deletion Segment." << endl;
                        dnm_encountered=true;
                        return;
                      }
                    }
                  }
                }
                add_H(gen, orig, term, subst, score, simplex, typer);
              }
            }

          
           // (b) Switch sides if necessary
           else{
              num_segs_removed = num_segments - term_segs.size();

              activate_half_mol(gen, orig, term, empty, -1, -1, false); 
              for ( int i=0; i<gen.num_atoms; i++ ){
                if ( gen.atom_active_flags[i] == false ){
                  if ( gen.atom_dnm_flag[i] ){
                    //cout << "false_Name: " << gen.atom_names[i] << endl;
                    //cout << "DNA5 Found in Deletion Segment." << endl;
                    dnm_encountered=true;
                    return;
                  }
                }
              }    
              double mol_size =0.0;
              for ( int i=0; i<gen.num_atoms; i++){
                if (gen.atom_active_flags[i] == true){
                  mol_size++;         
                }
              }

              // You could remove all of the molecule because you will be adding something new
              if ((subst == true) || ((mol_size/orig_size) >= 0.5)){
                if (dnm_enabled) {
                  for ( int i=0; i<gen.num_atoms; i++ ){
                    if ( gen.atom_active_flags[i] == false ){
                      if ( gen.atom_dnm_flag[i] ){
                        //cout << "false_Name: " << gen.atom_names[i] << endl;
                        //cout << "DNA7 Found in Deletion Segment." << endl;
                        dnm_encountered=true;
                        return;
                      }
                    }
                  }
                }
                // Add H to complete molecule
                add_H(gen, orig, term, subst, score, simplex, typer);
              }
              // Else switch sides
              else if ( subst == false && ((mol_size/orig_size) <= 0.5)) {
                // Switch side
                activate_half_mol(gen, term, orig, empty, -1, -1, false);
                if (dnm_enabled) {
                  for ( int i=0; i<gen.num_atoms; i++ ){
                    if ( gen.atom_active_flags[i] == false ){
                      if ( gen.atom_dnm_flag[i] ){
                        //cout << "false_Name: " << gen.atom_names[i] << endl;
                        //cout << "DNA4 Found in Deletion Segment." << endl;
                        dnm_encountered=true;
                        return;
                      } // end if for DNM check
                    } // end if for active_atom check
                  }  // end overall atom check for loop
                }
                add_H(gen, term, orig, subst, score, simplex, typer);
              }
           }
           double stop_time = time_seconds();
           cout << "\n" "-----------------------------------" "\n";
           if (verbose) cout << " Elapsed time deletions:\t" << stop_time - start_time
                << " seconds for GA Mutations\n\n";
           return;
        } // end if 3+ attachment points
        if (!subst){
                ostringstream new_energy;
                new_energy << gen.energy << "-del";
                gen.energy    = new_energy.str();
        }

    }// close while 

    double stop_time = time_seconds();
    if (verbose) cout << "\n" "-----------------------------------" "\n";
    if (verbose) cout << " Elapsed time deletions:\t" << stop_time - start_time
         << " seconds for GA Mutations\n\n";

   return;
}//end GA_Recomb::deletion_mutation()




// +++++++++++++++++++++++++++++++++++++++++
// This function will prepare the molecule for the addition of a fragment by the DN code
void
GA_Recomb::additions ( DOCKMol & anchor, int segment, Master_Score & score, Simplex_Minimizer & simplex, 
                       AMBER_TYPER & typer, Orient & orient, bool subst, int gen_num, int mut_attempt )
{
    Trace trace( "GA_Recomb::additions()" );
    double start_time = time_seconds();

    DN_GA_Build c_dn;
    Fragment frag;
    // STEP 1: Generate a fragment with Du, if necessary
    if (subst == true){
       frag.mol = anchor;
    }else{
         bool found = rand_H_to_Du (anchor, segment );
         // If a H was present
         if ( found ){
            frag.mol = anchor;
         }
         // Else exit function
         else{ return; }
    }
    // Clear mutants - has the anchor if we are doing a subt
    mutants.clear();

    // Populate fragment class information
    Fragment tmp_frag = frag.read_mol(frag);
    frag.used = false;
    frag.tanimoto = 0.0;
    frag.last_ap_heavy = -1;
    frag.scaffolds_this_layer = 0;    
    frag.torenv_recheck_indices.clear();
    frag.size = 0;
    frag.ring_size = 0;

    // If there are no aps, exit function
    // DELETE REMOVE
    if (tmp_frag.aps.size() == 0){
        cout <<"There are no attachment points. NO addition will be performed." <<endl;
        // DELETE - REMOVE
        /*fstream fout_molecules;
        fout_molecules.open ( "zzz.faulty-du.mol2", fstream::out|fstream::app );
        fout_molecules <<tmp_frag.mol.current_data << endl;
        Write_Mol2(tmp_frag.mol, fout_molecules);
        fout_molecules.close();*/
        return;
    } 
    c_dn.dn_user_specified_anchor = true;
    c_dn.anchors.push_back(tmp_frag);
    c_dn.dn_MW_cutoff_type_hard = true;
    c_dn.dn_MW_cutoff_type_soft = false;
    c_dn.dn_ga_gen = gen_num;
    c_dn.dn_ga_mut_attempt = mut_attempt;

    //------------------------------------------------------------------------
    // STEP 2: Pass necessary files and flags to DN
    // (a) Pass fragments to DN
    for (int i=0; i<tmp_sidechains.size(); i++){
        c_dn.scaf_link_sid.push_back( tmp_sidechains[i] );
    }

    if (subst == true){
        for (int i=0; i<tmp_linkers.size(); i++){
            c_dn.scaf_link_sid.push_back( tmp_linkers[i] );
        }
    }

    if (num_segs_removed > 3){
       for (int i=0; i<tmp_scaffolds.size(); i++){
           c_dn.scaf_link_sid.push_back( tmp_scaffolds[i] );
       }
    }

    // (b) Compute fragment descriptors 
    for (int i=0; i<c_dn.scaf_link_sid.size(); i++) {
        // Declare a temporary fingerprint object and compute atom environments (necessary for Gasteiger)
        activate_mol(c_dn.scaf_link_sid[i].mol);
        Fingerprint temp_finger;
        for (int j=0; j<c_dn.scaf_link_sid[i].mol.num_atoms; j++){
            c_dn.scaf_link_sid[i].mol.atom_envs[j] = temp_finger.return_environment(c_dn.scaf_link_sid[i].mol, j);
        }

        //cout << "Entering Amber typer" << endl;
        c_dn.scaf_link_sid[i].mol.prepare_molecule();
        // Compute the charges, saving them on the mol object
        float total_charges = 0;
        total_charges = compute_gast_charges(c_dn.scaf_link_sid[i].mol);
        c_dn.calc_mol_wt(c_dn.scaf_link_sid[i].mol);
        c_dn.calc_formal_charge(c_dn.scaf_link_sid[i].mol);
        c_dn.calc_rot_bonds(c_dn.scaf_link_sid[i].mol);
        c_dn.scaf_link_sid[i].mol.rot_bonds = 0;
    }

    // (c) Populate parameters for DN
    c_dn.dn_ga_flag = true;
    c_dn.dn_use_torenv_table = true;
    c_dn.dn_use_roulette = ga_use_dn_roulette;
    c_dn.dn_torenv_table = ga_torenv_table; 
    c_dn.read_torenv_table(c_dn.dn_torenv_table);
    c_dn.dn_sampling_method_ex = false;
    c_dn.dn_sampling_method_rand = ga_mut_sampling_method_rand;
    c_dn.dn_sampling_method_graph = ga_mut_sampling_method_graph;
    c_dn.dn_MW_cutoff_type_hard = true;
    c_dn.dn_MW_cutoff_type_soft = false;
    c_dn.dn_ga_mut_attempt = mut_attempt;

    if (c_dn.dn_sampling_method_rand){
        c_dn.dn_num_random_picks = ga_num_random_picks;
    }
    if (c_dn.dn_sampling_method_graph){
       c_dn.dn_graph_max_picks = ga_graph_max_picks;
       c_dn.dn_graph_breadth = ga_graph_breadth;
       c_dn.dn_graph_depth = ga_graph_depth;
       c_dn.dn_graph_temperature = ga_graph_temperature;
       c_dn.dn_temp_begin = c_dn.dn_graph_temperature;

       c_dn.prepare_fragment_graph( c_dn.scaf_link_sid, c_dn.scaf_link_sid_graph );
   }

    c_dn.dn_pruning_conformer_score_cutoff = 100.0;
    c_dn.dn_pruning_conformer_score_cutoff_begin = c_dn.dn_pruning_conformer_score_cutoff;
    c_dn.dn_pruning_conformer_score_scaling_factor = 1.0;

    c_dn.dn_pruning_clustering_cutoff = 100.0;
    c_dn.dn_constraint_mol_wt = ga_constraint_mol_wt;
    c_dn.dn_constraint_rot_bon = ga_constraint_rot_bon;
    c_dn.dn_constraint_formal_charge = ga_constraint_formal_charge;
    //c_dn.dn_constraint_charge_groups = 4;

    c_dn.dn_heur_unmatched_num = ga_heur_unmatched_num; 
    c_dn.dn_heur_matched_rmsd =  ga_heur_matched_rmsd;

    c_dn.dn_unique_anchors = ga_max_root_size;
    if ( subst ){
       c_dn.dn_max_grow_layers = num_segs_removed;
       if (verbose) cout << "#### Number of attempted layers of DN growth " << num_segs_removed << " " << c_dn.dn_max_grow_layers <<endl;
       c_dn.dn_max_root_size = ga_max_root_size;
    }
    else { 
       c_dn.dn_max_grow_layers = 1;
       c_dn.dn_max_root_size = ga_max_root_size;
       } 
    if (c_dn.dn_sampling_method_rand){
       c_dn.dn_max_layer_size = ga_num_random_picks;
    }
    else{
       c_dn.dn_max_layer_size = ga_graph_max_picks;
    }
    c_dn.dn_max_current_aps = 5;
    c_dn.dn_max_scaffolds_per_layer = 1;
    c_dn.dn_write_checkpoints = 0;
    c_dn.dn_write_prune_dump = 0;
    c_dn.dn_write_orients = 0;
    c_dn.dn_write_growth_trees = 0;
    c_dn.dn_output_prefix = "output";
    c_dn.verbose = 0;
    

    c_dn.use_internal_energy = use_internal_energy;
    c_dn.ie_att_exp = ie_att_exp;
    c_dn.ie_rep_exp = ie_rep_exp;
    c_dn.ie_diel = ie_diel;
    c_dn.ie_cutoff = ie_cutoff;

    // DELETE - REMOVE
    /*fstream fout_molecules;
    fout_molecules.open ( "zzz.du_add.mol2", fstream::out|fstream::app );
    for ( int i=0; i<c_dn.anchors.size(); i++ ){
        fout_molecules << c_dn.anchors[i].mol.current_data << endl;
        Write_Mol2( c_dn.anchors[i].mol, fout_molecules );
    }
    fout_molecules.close();*/


    //------------------------------------------------------------------------
    // STEP 3: Complete addition (Enter DN function)
    double start_time2 = time_seconds();
    c_dn.build_molecules( score, simplex, typer, orient );
    double stop_time2 = time_seconds();
    if (verbose) cout << "\n" "-----------------------------------" "\n";
    if (verbose) cout << " Elapsed time in DN only:\t" << stop_time2 - start_time2
         << " seconds for GA Mutations\n\n";

    // Clear vectors
    c_dn.anchors.clear();
    c_dn.scaf_link_sid.clear();
    
    // DELETE cout <<"number of tmp molecules " << c_dn.tmp_mutants.size()<<endl;

    //------------------------------------------------------------------------
    // STEP 4: Maintain the best scored molecules 
    // In DN, we save the results from each layer onto tmp_mutants - 
    // so the information from the last layer will be the final input
    bool found = false;
    int mol_counter = 0;

    // Sort molecules by score
    if (c_dn.tmp_mutants.size() > 0){
       if (c_dn.tmp_mutants.size() > 1){
          sort(c_dn.tmp_mutants.begin(), c_dn.tmp_mutants.end(), ga_fragment_sort);
       }
       // While you have yet to find a valid molecule
        while (mutants.size() == 0){
           // (a) Create a complete molecule 
           for ( int i=0; i<c_dn.tmp_mutants.size(); i++){
               mol_counter++;
               vector <int> heavy;

               // For each atom of the molecule - change all DU
               for ( int l=0; l<c_dn.tmp_mutants[i].mol.num_atoms; l++){

                   // Convert DU to H
                   if (c_dn.tmp_mutants[i].mol.atom_types[l] == "Du"){
                      heavy = c_dn.tmp_mutants[i].mol.get_atom_neighbors(l);

                      // Make a complete molecule
                      c_dn.dummy_to_H(c_dn.tmp_mutants[i], heavy[0], l);
                   }
               }

               // (b) There are no DU atoms so check the torenv of the entire molecule, if applicable
               // NOTE: DN only checked the torenvs of the bonds that were made in the design
               if (ga_use_torenv_table){
                  //cout << "Entering Amber typer" << endl;
                  c_dn.tmp_mutants[i].mol.prepare_molecule();
                  typer.skip_verbose_flag = true;
                  typer.prepare_molecule( c_dn.tmp_mutants[i].mol, true, score.use_chem, score.use_ph4, score.use_volume );

                  // Checking all rot bonds
	          prepare_torenv_indices(c_dn.tmp_mutants[i]);
                  if ( c_dn.valid_torenv_multi( (c_dn.tmp_mutants[i]) )){
                     // Activate and retain
                     activate_mol(c_dn.tmp_mutants[i].mol);
                     mutants.push_back(c_dn.tmp_mutants[i].mol);
                     found = true;
                  }

                  // Clear torenv 
                  c_dn.torenv_vector.clear();
                  c_dn.tmp_mutants[i].torenv_recheck_indices.clear();
               }
               // If not using torenv, save molecule
               else{
                  mutants.push_back(c_dn.tmp_mutants[i].mol);
                  found = true;
               }

               // Exit loop if there are no more molecules 
               if (i == c_dn.tmp_mutants.size()-1){found = true;}
               if (found) {break;}
           }
           // Exit loop if there are no more molecules 
           if (mol_counter == c_dn.tmp_mutants.size()-1){found = true;}
           if(found){break;}
        }
     }
     //cout << "mutants size " << mutants.size() << endl;// DELETE REMOVE
            
    // If there are mutants
    //if (mutants.size() > 0 ){
       // DELETE - REMOVE
       //fstream fout_molecules;
       /*fout_molecules.open ( "zzz.successful-addition.mol2", fstream::out|fstream::app );
       fout_molecules << mutants[0].current_data << endl;
       Write_Mol2( mutants[0], fout_molecules );
       fout_molecules.close();*/

       // Update energy to include addition
       //ostringstream new_energy;
       //new_energy << mutants[0].energy << "-add";
       //mutants[0].energy    = new_energy.str();
    //}

    // Clear molecules
    c_dn.torenv_vector.clear();
    c_dn.anchors.clear();
    c_dn.tmp_mutants.clear();
    frag.mol.clear_molecule();
    tmp_frag.mol.clear_molecule();

    double stop_time = time_seconds();
    if (verbose) cout << "\n" "-----------------------------------" "\n";
    if (verbose) cout << " Elapsed time additions:\t" << stop_time - start_time
         << " seconds for GA Mutations\n\n";
    return;
} // end GA_recomb:additions()



// +++++++++++++++++++++++++++++++++++++++++
// Replaces a scaffold or linker with a fragment from the library - maintains to best scored option
void
GA_Recomb::replacement( DOCKMol & gen, int final_seg, INTVec bonds, INTVec origV, INTVec termV, Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & typer, Orient & orient )
{
    Trace trace( "GA_Recomb::replacement()" );
    double start_time = time_seconds();

    // Convert mol to fragment 
    Fragment frag = mol_to_frag(gen);

    // Populate additional fragment information (size, etc)
    prepare_replacement_segment(frag, final_seg, bonds, origV, termV, typer);
    
    // Call replacement function
    if ( frag.aps.size() == 2 ){
       rand_replacement (frag, final_seg, tmp_linkers, score, simplex, typer, orient);
    }else{
       rand_replacement (frag, final_seg, tmp_scaffolds, score, simplex, typer, orient);
    }

    // Sort results and keep best score
    sort(mutants.begin(), mutants.end(), mol_sort);

    /*if ( mutants.size() > 0 ){
       // Update energy to include replacement
       ostringstream new_energy;
       new_energy << mutants[0].energy << "-repl";
       mutants[0].energy    = new_energy.str();
 
       // DELETE - REMOVE
       fstream fout_molecules;
       fout_molecules.open ( "zzz.successful_replace.mol2", fstream::out|fstream::app );
       fout_molecules << gen.current_data << endl;
       Write_Mol2( gen, fout_molecules );
       for ( int i=0; i<mutants.size(); i++ ){
           fout_molecules << mutants[i].current_data << endl;
           Write_Mol2( mutants[i], fout_molecules );
       }
       fout_molecules.close();
    }
    // DELETE
    else{
       // DELETE - REMOVE
       fstream fout_molecules;
       fout_molecules.open ( "zzz.unsuccessful_replace.mol2", fstream::out|fstream::app );
       fout_molecules << gen.current_data << endl;
       Write_Mol2( gen, fout_molecules );
       deactivate_mol(gen);
       for (int i=0; i<orig_segments[final_seg].num_atoms; i++){
           gen.atom_active_flags[orig_segments[final_seg].atoms[i]] = true;
       }
       for (int i=0; i<orig_segments[final_seg].num_bonds; i++){
           gen.bond_active_flags[orig_segments[final_seg].bonds[i]] = true;
       }
       fout_molecules << gen.current_data << endl;
       Write_Mol2( gen, fout_molecules );
       fout_molecules.close();

    }*/
    double stop_time = time_seconds();
    if (verbose) cout << "\n" "-----------------------------------" "\n";
    if (verbose) cout << " Elapsed time replacement:\t" << stop_time - start_time
         << " seconds for GA Mutations\n\n";
    return;
} //end GA_Recomb::replace()




// +++++++++++++++++++++++++++++++++++++++++
// Prepare the segment for replacement
// Reference ReCore
void
GA_Recomb::prepare_replacement_segment( Fragment & frag, int segment, INTVec aps_bonds, INTVec origV, INTVec targV, AMBER_TYPER & typer )
{
    Trace trace( "GA_Recomb::prepare_replacement_segment()" );
    // STEP 1: Define segment size 
    // Ensure that atom_ring_flag is being properly populated
    frag.mol.id_ring_atoms_bonds();

    // For the segment's atoms
    // Define size of heavy atoms only (non-ring) 
    for (int i=0; i<orig_segments[segment].num_atoms; i++){
        if (( frag.mol.atom_types[orig_segments[segment].atoms[i]] != "H") && (frag.mol.atom_types[orig_segments[segment].atoms[i]] != "Du")){
           // If the atoms are not in a ring (for heavy atoms only)
           // Check for rings
           if ( frag.mol.atom_ring_flags[orig_segments[segment].atoms[i]] == true){
              frag.ring_size++;
           }
           else{ frag.size++; }
       }
    }//cout << "size " << frag.size << " ring_size " << frag.ring_size << " "; //DELETE - REMOVE 


    //------------------------------------------------------------------------
    // STEP 2: Characterize the "aps" for the molecule
    // For each bond that should be replaced
    for (int i=0; i<aps_bonds.size(); i++){ 
        // Do not know if the orig or target atoms are on the current segment
        for ( int j=0; j<orig_segments[segment].num_atoms; j++ ){
            if ( (orig_segments[segment].atoms[j] == origV[i] ) ){
               // Update the aps information
               // Create temporary attachment point object
               AttPoint tmp_ap;
               
               // Copy the heavy atom / dummy atom pairs onto the attachment point class
               tmp_ap.heavy_atom = origV[i];
               tmp_ap.dummy_atom = targV[i];

               // And add that information to the attachment points (aps) vector of the scaffold
               frag.aps.push_back(tmp_ap);

               // Create temp pair for bond information - CS: 09/19/16
               std::pair <int, std::string> bond_info;
               // Update the bond number and type - CS: 09/19/16
               bond_info.first = aps_bonds[i];
               bond_info.second = frag.mol.bond_types[aps_bonds[i]];
               frag.aps_bonds_type.push_back(bond_info);                  
            }
            else if ( (orig_segments[segment].atoms[j] == targV[i]) ){
               // Update the aps information
               // Create temporary attachment point object
               AttPoint tmp_ap;
               
               // Copy the heavy atom / dummy atom pairs onto the attachment point class
               tmp_ap.heavy_atom = targV[i];
               tmp_ap.dummy_atom = origV[i];

               // And add that information to the attachment points (aps) vector of the scaffold
               frag.aps.push_back(tmp_ap);

               // Create temp pair for bond information - CS: 09/19/16
               std::pair <int, std::string> bond_info;
               // Update the bond number and type - CS: 09/19/16
               bond_info.first = aps_bonds[i];
               bond_info.second = frag.mol.bond_types[aps_bonds[i]];
               frag.aps_bonds_type.push_back(bond_info);                  
            }
        }
    }
    //cout << "num aps " << frag.aps.size() <<endl;

    //------------------------------------------------------------------------
    // STEP 3: Define the distances between each aps
    // Detemine the size of the aps dist matrix
    std::vector < std::vector <float> > tmp_aps_dist(frag.aps.size(), std::vector <float> (frag.aps.size()));

    DOCKVector vec1;
    // For each aps
    for (int i=0; i<aps_bonds.size(); i++){
        for (int j=i+1; j<aps_bonds.size(); j++){
            vec1.x = frag.mol.x[frag.aps[j].heavy_atom] - frag.mol.x[frag.aps[i].heavy_atom];
            vec1.y = frag.mol.y[frag.aps[j].heavy_atom] - frag.mol.y[frag.aps[i].heavy_atom];
            vec1.z = frag.mol.z[frag.aps[j].heavy_atom] - frag.mol.z[frag.aps[i].heavy_atom];

            // Measure the distance between the heavy atoms of the attachment points
            float dist1 = ( vec1.x*vec1.x) +  ( vec1.y*vec1.y) + ( vec1.z*vec1.z);
            tmp_aps_dist[i][j] = tmp_aps_dist[j][i] = dist1;
            /*cout <<"vec1.x = " << vec1.x << " vec1.y " << vec1.y << " vec1.z " << vec1.z
                 << " atom type j " << frag.mol.atom_types[frag.aps[j].heavy_atom] << " i " << frag.mol.atom_types[frag.aps[i].heavy_atom] << endl;
         cout << " distances j " << frag.mol.x[frag.aps[j].heavy_atom] << " " << frag.mol.y[frag.aps[j].heavy_atom] << " " << frag.mol.z[frag.aps[j].heavy_atom] << endl;
          cout << " distances i " << frag.mol.x[frag.aps[i].heavy_atom] << " " << frag.mol.y[frag.aps[i].heavy_atom] << " " << frag.mol.z[frag.aps[i].heavy_atom] << endl;*/
       }
    }
    frag.aps_dist = tmp_aps_dist;

    // REMOVE 
    // Deactive dummy atoms and bonds of the frag from fraglib
    /*for ( int i=0; i<frag.aps.size(); i++ ){
        frag.mol.atom_active_flags[frag.aps[i].dummy_atom] = false;
    }
    fstream fout_molecules;
    fout_molecules.open ( "zzz.prepare_frag.mol2", fstream::out|fstream::app );
    fout_molecules << frag.mol.current_data << endl;
    Write_Mol2( frag.mol, fout_molecules );
    fout_molecules.close();*/
  
    return;
} // end GA_Recomb::prepare_replacement_segment()



// +++++++++++++++++++++++++++++++++++++++++
// Complete a single scaffold/linker replacement by randomly selecting a user-specified number of fragments
// The same fragment cannot be replaced with itself
// Only rings can replace rings and non-rings can replace rings
// Reference ReCore
void
GA_Recomb::rand_replacement( Fragment & ref, int final_seg, std::vector <Fragment> & fraglib, Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & typer, Orient & orient )
{
    Trace trace( "GA_Recomb::rand_replacement()" );
    double          start_time = time_seconds();
    int invalid_torenv = 0;
    // If the fragment library is empty for some reason, return without doing anything
    if (fraglib.size() == 0) { return; }

    // Create a temp Frag vector to hold the frag library
    std::vector <Fragment> tmp_fraglib;
    for ( int i=0; i<fraglib.size(); i++ ){
        tmp_fraglib.push_back(fraglib[i]);
    }
    
    // While N choices...
    int current = 0;
    int passes = 0;

    // This loop will iterate ga_num_random_picks times where every possible 
    // allowed conformation of the selected fragment will be generated. 
    // Current is incremeted whether the fragment is allowed or results in a molecule.
    while (current < ga_num_random_picks){

        // The passes int counts every time we either keep or skip a fragment from the libary, so
        // once its size reaches the fraglib.size(), we have looked at everything
        if (passes >= fraglib.size()){ return; }

        // STEP 1: ID a compatible fragment for replacement

        // (1) Choose a random starting position,
        // (rand has been seeded if and only if the minimizer is turned on)
        int frag_index = rand() % tmp_fraglib.size();
        cout << "Selecting frag: " << tmp_fraglib[frag_index].mol.energy << endl;
        // (2) Only replace rings with other rings
        if ( ref.ring_size > 0 ){
           if ( tmp_fraglib[frag_index].ring_size == 0 ){ 
           // Delete the molecule from the vector
           tmp_fraglib.erase(tmp_fraglib.begin()+frag_index);
               current++;
               passes++;
               cout << "deleted frag from vector because no ring" << endl;
               continue; 
            }
        }
        // (3) Replace non-rings with non-rings
        // NOTE - techically a ring where all of the aps are on 1 heavy atom would work
        else if ( ref.ring_size == 0 ){ 
           if ( tmp_fraglib[frag_index].ring_size != 0 ){ 
           // Delete the molecule from the vector
           tmp_fraglib.erase(tmp_fraglib.begin()+frag_index);
              current++;
              passes++;
              cout << "deleted frag from vector because it has a ring" << endl;
              continue; 
           }
        }

        
        // (4) Determine if the aps have the same bond types
        // Will compare all combinations of bond types and save to as pair with aps number
        std::vector <std::vector <int> > compare_aps_bond_types = compare_bond_types( ref, tmp_fraglib[frag_index] );
        // If the bond types to not match, remove the fragment and continue
        if ( compare_aps_bond_types.size() == 0 ){
           compare_aps_bond_types.clear();
           // Delete the molecule from the vector
           tmp_fraglib.erase(tmp_fraglib.begin()+frag_index);
           cout << "Deleted fragment because bond types did not match" << endl;
           current++;
           passes++;
           continue; 
        }

        // (5) Make sure the fragment is different from the segment
        // Extract, via activation, the segment of the current mol
        // NOTE: Hungarian only ignores H atoms, not Du 
        deactivate_mol(ref.mol);
        // The dummy atoms of the ref are not activated
        for (int i=0; i<orig_segments[final_seg].num_atoms; i++){
            ref.mol.atom_active_flags[orig_segments[final_seg].atoms[i]] = true;
        }
        for (int i=0; i<orig_segments[final_seg].num_bonds; i++){
            ref.mol.bond_active_flags[orig_segments[final_seg].bonds[i]] = true;
        }

        // Deactive dummy atoms and bonds of the frag from fraglib
        for ( int i=0; i<tmp_fraglib[frag_index].aps.size(); i++ ){
            tmp_fraglib[frag_index].mol.atom_active_flags[tmp_fraglib[frag_index].aps[i].dummy_atom] = false;
        }
        
        // Idenfity whether the active sections are the same
        Hungarian_RMSD h;
        pair <double, int> result;

        result = h.calc_Hungarian_RMSD_dissimilar( ref.mol, tmp_fraglib[frag_index].mol );
        //DELETE cout << "calc H rmsd = " << result.first << " unmatched = " << result.second <<endl;
        
        // Hard coded - so that we do not want to pick the same segment
        if ((result.second < 1)) { 
           compare_aps_bond_types.clear();
           // Delete the molecule from the vector
           tmp_fraglib.erase(tmp_fraglib.begin()+frag_index);
           cout << "Deleted fragment because its too similar to chosen segment" << endl;
           current++;
           passes++;
           continue;
        }

        // (6) Make sure that the allowed bond combinations can result in a molecule
        std::vector <int> allowed_results = allowed_bond_combos ( ref, tmp_fraglib[frag_index], compare_aps_bond_types );

        // If the allowed bool (third term) is false, remove the fragment and continue
        if ( allowed_results[2] == 0 ){
           allowed_results.clear();
           // Delete the molecule from the vector
           tmp_fraglib.erase(tmp_fraglib.begin()+frag_index);
           cout << "Deleted fragment because not matching bondtypes" << endl;
           current++; 
           passes++;
           continue; 
        }
        
        //------------------------------------------------------------------------
        // STEP 2: For non-ring linkers and scaffolds, complete recombination
        // Allowed fragment
        bool allowed = true;
        if ( ref.ring_size > 0 ){
           // (1) Determine if the distances between the same bt aps in each molecule are similar
           // NOTE: This step does not need to be completed if there is only one heavy atom in the segment
           // Since all distances will be zero

           // Matrix of possible matches for a given ref aps
           std::vector < std::vector <int> > final_compare_aps;

           // NOTE: If there are only two aps, then 1 check is necessary
           // For each allowed pair
           for ( int i=0; i<compare_aps_bond_types.size(); i++ ){
               float ref_total = 0;
               for ( int k=0; k<ref.aps.size()-1; k++ ){
                     ref_total += ref.aps_dist[i][k];
               } 
               float delta_dist = 0;
               // Save old values
               float min_delta = 1000;

               for ( int j=0; j<compare_aps_bond_types[i].size(); j++ ){

                   // Compute the change in distance
                   // If the aps pair is allowed
                   if (compare_aps_bond_types[i][j] == 1){
                      // Compute the difference between all of the distances the selected aps
                      // and the rest of the aps in the respective molecule
                      //Try adding the together first!
                      // -1 since the distance does not include itself
                      float frag_total = 0;
                      for ( int k=0; k<tmp_fraglib[frag_index].aps.size()-1; k++ ){
                          frag_total += tmp_fraglib[frag_index].aps_dist[j][k];
                      } 

                      delta_dist = abs( ref_total - frag_total );
                      // Update previous data entry with smaller distance if applicable
                      if ( delta_dist < min_delta ){ min_delta = delta_dist; }
                   }
              }
              
              // Check to make sure that each aps has an allowable attachment
              // Check the smallest distance only
              if ( min_delta > 1.0 ){ 
                 allowed = false;
                 cout << "Deleted frag because of distance cutoff" << endl; 
                 break;
              }
           }
        }

        //------------------------------------------------------------------------
        // STEP 3: Find the largest side of the molecule, which will not move during attachement
        if ( allowed == true ){
           // Save the fractional size of the molecule on each aps
           DOCKMol tmp;
           float mol_size = 0;
           int aps_num = -1;
           
           // (1) Determine which aps is on the larger side of the molecule
           for ( int i=0; i<ref.aps.size(); i++){
               activate_half_mol(ref.mol, ref.aps[i].heavy_atom, ref.aps[i].dummy_atom, tmp, -1, -1, true);
               int active_atoms = 0;
               // Proportion of active atoms
               for ( int j=0; j<ref.mol.num_atoms; j++){
                   if ( ref.mol.atom_active_flags[j] == true ){
                      active_atoms++;
                   }
               }
               if ( (active_atoms > mol_size) ){
                  // Update mol_size and the aps number
                  mol_size = active_atoms;
                  aps_num = i;
               }
           }

           // Define DOCK vec to hold the frowign molecule
           std::vector <DOCKMol> recomb;

           // (2) Orient the fraglib to the reference
           std::vector <DOCKMol> final_orients = orient_frag_to_ref ( ref, 
                                    tmp_fraglib[frag_index], final_seg, score, typer, orient );
           //cout << "Trying fragment: " << tmp_fraglib[frag_index].mol.energy << endl; //LEPLEPLEP
           // If there are no orients, continue to next fragment
           if ( final_orients.size() == 0 ){
              // Delete the molecule from the vector
              tmp_fraglib.erase(tmp_fraglib.begin()+frag_index);
              passes++;
              cout << "unsuccessful orienting" << endl;
              // DELETE
              /*fstream fout_molecules;
              fout_molecules.open ( "zzz.could_not_orient.mol2", fstream::out|fstream::app );
              activate_mol(ref.mol);
              fout_molecules << ref.mol.current_data << endl;
              Write_Mol2( ref.mol, fout_molecules );
              activate_mol(tmp_fraglib[frag_index].mol);
              fout_molecules << tmp_fraglib[frag_index].mol.current_data << endl;
              Write_Mol2( tmp_fraglib[frag_index].mol, fout_molecules );
              fout_molecules.close();*/
              continue;
           }

           //------------------------------------------------------------------------
           // STEP 4: Complete recombination, keep mols with allowed bonds, and minimize
           // (3) Generate 
           //cout << "### Entering Replacement for aps: " << ref.aps.size() <<endl; 
           //cout << "Number of molecules prior to replace_scaffold: " << final_orients.size() << endl;
           recomb = replace_scaffold ( ref, aps_num, final_seg, allowed_results[0], allowed_results[1], 
                                       tmp_fraglib[frag_index], final_orients, typer, score );
           //cout << "Number of molecules prior to multi torenv check: "<< recomb.size() << endl; 
           // (4) Check each rotatable bond's torenv
           for ( int k=0; k<recomb.size(); k++ ){
               if ( ga_use_torenv_table ){
                  //In case torenv depends on atom types?
                  Fragment new_recomb = mol_to_frag (recomb[k]);
                  //cout << "Entering Amber typer" << endl;
                  new_recomb.mol.prepare_molecule();
                  typer.skip_verbose_flag = true;
                  typer.prepare_molecule( new_recomb.mol, true, score.use_chem, score.use_ph4, score.use_volume );
	              prepare_torenv_indices(new_recomb);
                  DN_GA_Build c_dn;
                  c_dn.dn_ga_flag = true;
                  c_dn.dn_use_roulette = ga_use_dn_roulette;
                  c_dn.dn_use_torenv_table = ga_use_torenv_table;
                  c_dn.dn_torenv_table = ga_torenv_table; 
                  c_dn.read_torenv_table(c_dn.dn_torenv_table);
                  c_dn.verbose = 0;
                  c_dn.dn_MW_cutoff_type_hard = true;
                  c_dn.dn_MW_cutoff_type_soft = false;                  
 
                  if ( c_dn.valid_torenv_multi( new_recomb ) ){
                     // Minimize and save to mutants vector
                     //cout << "entering minimize children1" << endl;
                     activate_mol(recomb[k]);
                     //cout << "score before: " << recomb[k].current_score << endl;
                     minimize_children( recomb[k], mutants, score, simplex, typer );
                     //cout << "score after: " << recomb[k].current_score << endl;
                     //cout << "Number of molecules after minimization: " << mutants.size() << endl;
                 }
                 else{ invalid_torenv++; 
                cout << "Invalid torenv" << endl;}
                //cout << "Invalid torenv" << endl;
                 new_recomb.mol.clear_molecule();
                 c_dn.torenv_vector.clear();
                 new_recomb.torenv_recheck_indices.clear();
              }
              else {
                 // Minimize and save to mutants vector
                 //cout << "entering minimize children2" << endl;
                 //cout << "score before: " << recomb[k].current_score << endl;
                 minimize_children( recomb[k], mutants, score, simplex, typer );
                 //cout << "Number of molecules after minimization: " << mutants.size() << endl;
              }
           }
           cout << "Number of molecules after minimization: " << mutants.size() << endl;
        }
        //cout << "Number of molecules after minimization: " << mutants.size() << endl;
        // Check ligand constraints (Mw, RB, etc)
        //cout << "Number of molecules in mutants: " << mutants.size() << endl;
        /*if (mutants.size() > 0){
            hard_filter(mutants); 
        }*/
        // Number of picks increments regardless of success
        current++;
        //cout << "counter: " << current << endl;
        passes++;
        // Delete the molecule from the vector
        compare_aps_bond_types.clear();
        allowed_results.clear();
        tmp_fraglib.erase(tmp_fraglib.begin()+frag_index);
    } 
    if (mutants.size() > 0){
        hard_filter(mutants);
    }

    //cout << "#### It took " << passes << " number of passes to reach " << current << " num molecules. Num picks: " << ga_num_random_picks <<endl;
    cout << "Number with invalid torenv " << invalid_torenv <<endl;
    double stop_time = time_seconds();
    if (verbose) cout << "\n" "-----------------------------------" "\n";
    if (verbose) cout << " Elapsed time rand_replacement:\t" << stop_time - start_time
         << " seconds for GA Mutations\n\n";
   return;
}//end GA_Recomb::rand_replacement()







// +++++++++++++++++++++++++++++++++++++++++
// This function will populate a bool matrix for aps indices with the same bond type 
// between the reference segment and new fragment. True = 1.
// It will return a matrix of size zero if one bond type cannot be found.
std::vector < std::vector <int> >
GA_Recomb::compare_bond_types ( Fragment & ref, Fragment & frag )
{
    Trace trace( "GA_Recomb::compare_bond_types()" );
    // Define a matrix to store the bond number with the same type as ref
    std::vector < std::vector <int> > aps_bond_types(ref.aps.size(), std::vector <int> (frag.aps.size()));
 
    // For each attachment point in the reference,
    for ( int i=0; i<ref.aps.size(); i++ ){
        // Temp int to hold information for each aps
        int same_order = 0;

        // Compare with each aps of the frag, until you find a match
        for ( int j=0; j<frag.aps.size(); j++ ){
            // If they are the same, save to the temp vector
            if ( ref.aps_bonds_type[i].second == frag.aps_bonds_type[j].second ){
               aps_bond_types[i][j] = 1;
               same_order++;
            }
        }

        // After comparing, determine if there are any aps that are then same
        if ( same_order == 0 ){
           // Erase all previous data, until the row before the current
           aps_bond_types.clear();
           return aps_bond_types;
        }
    }  
    return aps_bond_types;
} //end GA_Recomb::compare_bond_types()




// +++++++++++++++++++++++++++++++++++++++++
// Function will deteremine if all of the reference aps have the same or unique bond types,
// will compare the fragments bond types with the reference, and return three bool values
// (1) all_same_bt, (2) all_diff_bt, (3) true is ref matches the frag bt
std::vector <int>
GA_Recomb::allowed_bond_combos ( Fragment & ref, Fragment & frag, std::vector < std::vector <int> > & compare_aps_bond_types )
{
    Trace trace( "GA_Recomb::allowed_bond_combos()" );
    // Initialize a vector of ints to be returned (Positions represent 1. all_same, 2. all_diff, 3. allowed)
    std::vector <int> true_comparison_vec;

    // Save a matrix of similar bond types
    std::vector < std::vector <int> > same_bond_types(ref.aps.size(), std::vector <int>(ref.aps.size()));
    // Initialize diagonal
    for ( int i=0; i<ref.aps.size(); i++ ){
        same_bond_types[i][i] = 0;
    }


    // STEP 1: Determine whether all of the aps are the same, all different, or a combination
    bool all_same_bt = true;
    bool all_diff_bt = true;

    // Determine if all same or different bt
    // For each bond type
    for ( int i=0; i<ref.aps.size(); i++ ){
        // Compare to each bond type 
        for ( int j=i+1; j<ref.aps.size(); j++ ){
            // If a pair is the same, then it is possible that all are the same
            // but not that all are different
            if ( ref.aps_bonds_type[i].second == ref.aps_bonds_type[j].second ){

               // Update index
               same_bond_types[i][j] = same_bond_types[j][i] = 1;   
 
               // However, if the same was previously false, then all cannot be the same
               if ( all_same_bt == false ){
                  all_same_bt = false;
               }
               all_diff_bt = false;
            }

            // If a pair are different, then it is possible that all are different
            // but not all the same
            else{
               // Update index
               same_bond_types[i][j] = same_bond_types[j][i] = 0;   
               // However, if all diff was previously false, then all cannot be different
               if ( all_diff_bt == false ){
                  all_diff_bt = false;
               }
               all_same_bt = false;
            }
        }
    }

    // Update the comparison vector
    true_comparison_vec.push_back(all_same_bt);
    true_comparison_vec.push_back(all_diff_bt);

    //------------------------------------------------------------------------
    // STEP 2: Compare the bond types between the ref and frag
    // (a) If all the same, then the frag should have all the same bt of the same order
    // therefore, the first aps of the ref should match all of the aps of the frag

    // If all the same, only check the first aps
    if ( all_same_bt == true ){
       // If there is a match, the sum of the bools for the compare_aps_bond_types vec
       // should match with each aps of the ref
       int bool_sum = 0;
       for ( int i=0; i<compare_aps_bond_types[0].size(); i++ ){
           if ( compare_aps_bond_types[0][i] == 1 ){
              bool_sum++;
           }
       }
      
       // If the frag does not have the same bt as the ref,
       if ( bool_sum != ref.aps.size() ){  
          // Update comparison vector
          true_comparison_vec.push_back(0);
          // Exit function
          return true_comparison_vec;
       } 
    }

    // (b) If all different, then each aps should only match with 1 aps of the frag
    else if ( all_diff_bt == true ){
       // The sum of the bools for the compare_aps_bond_types vec
       // should be 1

       // Hold frag aps that match ref
       INTVec used;
       for ( int i=0; i<frag.aps.size(); i++ ){
           used.push_back(0);
       }
      
       // Check all aps, until criteria not matched
       for ( int i=0; i<ref.aps.size(); i++){
           int bool_sum = 0;
           for ( int j=0; j<compare_aps_bond_types[i].size(); j++ ){
               if ( compare_aps_bond_types[i][j] == 1 ){
                  bool_sum++;
                  used[j] = 1;
               }
           }
           // If the ref aps matches with more than one frag aps, exit
           if ( bool_sum > 1 ){  
              // Update comparison vector
              true_comparison_vec.push_back(0);
              return true_comparison_vec;
           }
           // If the ref aps only matches with 1 frag aps, but it is redundant exit
           int used_sum=0;
           for ( int j=0; j<used.size(); j++ ){
               if ( used[j] == 1 ){ used_sum++; }
           }
           if ( used_sum < i+1 ){
              // Update comparison vector
              true_comparison_vec.push_back(0);
              return true_comparison_vec;
           }
       } 
    }

    // (c) If there are a combination, then check all of the bonds
    else{
       INTVec used;
       for ( int i=0; i<ref.aps.size(); i++ ){
           used.push_back(0);
       }
 
       for ( int i=0; i<ref.aps.size(); i++ ){
           // If it is not used
           if ( used[i] == 0 ){
              // Determine whether the total matched bonds are the same as the total similar
              int bool_sum = 0;
              for ( int j=0; j<compare_aps_bond_types[i].size(); j++ ){
                  if ( compare_aps_bond_types[i][j] == 1 ){ bool_sum++; }
              }
 
               // Find sum of number of aps (ref) with the same bt
              int same_aps = 0;
              for ( int j=0; j<same_bond_types[i].size(); j++ ){
                  if ( same_bond_types[i][j] == 1 ){ same_aps++; }
              }
              if ( bool_sum != same_aps ){
                 // Update comparison vector
                 true_comparison_vec.push_back(0);
                 return true_comparison_vec;
              }
              // Else it works so update used vector
              else{
                 used[i] = 1;
                 // Also include bonds of the same atom types
                 for (int j=0; j<same_bond_types[i].size(); j++ ){
                    if ( same_bond_types[i][j] == 1 ){ used[i] = 1; }
                 }
             }
          } 
       }
    }

    // If you get here, then the pair is allowed
    true_comparison_vec.push_back(1);
    return true_comparison_vec;
}// end GA_Recomb::allowed_bond_combos()



// +++++++++++++++++++++++++++++++++++++++++
// Generates orients of the frag to match the reference and prunes orients that do not
// have a good overlap with all DU of the ref
std::vector <DOCKMol>
GA_Recomb::orient_frag_to_ref ( Fragment & ref, Fragment & frag, int final_seg, 
                                Master_Score & score, AMBER_TYPER & typer, Orient & orient )
{
    Trace trace( "GA_Recomb::orient_frag_to_ref()" );

    // STEP 1: Initialize default input parameters
    orient.use_chemical_matching = false;
    orient.orient_ligand = true;
    orient.automated_matching = true;
    orient.tolerance = 0.25;
    orient.dist_min = 2.0;
    orient.min_nodes = 3;
    orient.max_nodes = 10;
    orient.max_orients = 50;
    orient.critical_points = true; // Requires a DU atom be be in the sphere matches
    orient.use_chemical_matching = false;
    // This refers othe the frag input sphere file
    orient.use_ligand_spheres = false;
    orient.verbose = 0;   

    //------------------------------------------------------------------------
    // STEP 2: Prepare and conduct matching
    // (a) Prepare and activate fragment
    //cout << "Entering Amber typer" << endl;
    frag.mol.prepare_molecule();
    typer.skip_verbose_flag = true;
    typer.prepare_molecule( frag.mol, true, score.use_chem, score.use_ph4, score.use_volume );
    activate_mol(frag.mol);
 
    // (b) Only activate the ref mol, including dummy atoms 
    deactivate_mol(ref.mol);
    // Activate the segment
    for (int i=0; i<orig_segments[final_seg].num_atoms; i++){
         ref.mol.atom_active_flags[orig_segments[final_seg].atoms[i]] = true;
    }
    for (int i=0; i<orig_segments[final_seg].num_bonds; i++){
        ref.mol.bond_active_flags[orig_segments[final_seg].bonds[i]] = true;
    }
    // Also activate the aps and its dummy atoms
    for (int i=0; i<ref.aps.size(); i++){
        // Activate dummy
        ref.mol.atom_active_flags[ref.aps[i].dummy_atom] = true;
        // Activate bond
        ref.mol.bond_active_flags[ref.aps_bonds_type[i].first] = true;
    } 
   
    // (c) Generate the spheres from the ref mol
    orient.get_lig_reference_spheres(ref.mol);

    // (d) Generate the orients
    // Generate the list of atom-sphere matches
    std::vector <DOCKMol> tmp;
    orient.match_ligand(frag.mol);
                         
    // Transforms frag.mol to one of the matches in the list
    // DELETE - REMOVE
    /*fstream fout_molecules;
    fout_molecules.open ( "zzz.all_orients.mol2", fstream::out|fstream::app );
    fout_molecules <<ref.mol.current_data << endl;
    Write_Mol2(ref.mol, fout_molecules);
    fout_molecules.close();*/

    // Generate and save orients
    while (orient.new_next_orientation(frag.mol)){
       // DELETE - REMOVE
       /*fstream fout_molecules;
       fout_molecules.open ( "zzz.all_orients.mol2", fstream::out|fstream::app );
       fout_molecules <<frag.mol.current_data << endl;
       Write_Mol2(frag.mol, fout_molecules);
       fout_molecules.close();*/

       tmp.push_back(frag.mol);             
    }
    // Clear orient vectors
    orient.clean_up();

    // DELETE
    //cout << "Number of orients (pre-pruning): " << tmp.size() << endl;

    //------------------------------------------------------------------------
    // STEP 3: Remove orients that do not have the dummys overlapped
    std::vector <DOCKMol> pruned_orients;
    for ( int i=0; i<tmp.size(); i++ ){

        // Initialize bool vec to detemine which aps overlap
        INTVec used_aps;
        for ( int j=0; j<ref.aps.size(); j++ ){
            used_aps.push_back(0);
        }

        // Compare with each aps of the reference
        for ( int j=0; j<ref.aps.size(); j++ ){
            std::vector <std::pair <int, float> > results =  overlapping_aps ( ref, j, frag, tmp[i], typer );
            // If there are no matches
            if ( results.size() == 0 ){ break; }
            // If there was a match, check to see if it was used
            else {
               // Update used set for all matches - usually only 1 match
               for ( int k=0; k<results.size(); k++ ){
                   // Hard-coded rmsd 
                   if ( results[k].second < 0.25 ) {
                      used_aps[results[k].first] = 1;
                   }
               }
            }
        }
        // If all of the aps are used ( ie close to a ref du), then keep the orient for further analysis
        int sum_aps = 0;
        for ( int j=0; j<used_aps.size(); j++ ){
            if ( used_aps[j] == 1 ){ sum_aps++; }
        }
   
        if ( sum_aps == ref.aps.size() ){ pruned_orients.push_back(tmp[i]);  }
    }

    //cout << "Number of orients (post-pruning): " << pruned_orients.size() <<endl;
    // DELETE - REMOVE
    /*fout_molecules.open ( "zzz.pruned_orients.mol2", fstream::out|fstream::app );
    activate_mol(ref.mol);
    fout_molecules <<ref.mol.current_data << endl;
    Write_Mol2(ref.mol, fout_molecules);
    for ( int i=0; i<pruned_orients.size(); i++ ){
        activate_mol(pruned_orients[i]);
        fout_molecules <<pruned_orients[i].current_data << endl;
        Write_Mol2(pruned_orients[i], fout_molecules);
    }
    fout_molecules.close();*/

    //------------------------------------------------------------------------
    // STEP 4: Clear orient parameters
    orient.use_chemical_matching = false;
    orient.orient_ligand = false;
    orient.automated_matching = false;
    orient.critical_points = false;
    orient.use_chemical_matching = false;
    orient.use_ligand_spheres = false;
    orient.verbose = 0;   

    // Clear vectors in orient
    orient.spheres.clear();
    orient.centers.clear();

    delete[]orient.sph_dist_mat;//Release memory
    orient.sph_dist_mat = NULL; //Release memory from pointer

    return pruned_orients;
}// end GA_Recomb::orient_frag_to_ref() 



// +++++++++++++++++++++++++++++++++++++++++
// Measure the Hungarian RMSD of every DU of the fragment tot eh given ref aps
std::vector <std::pair <int, float> >
GA_Recomb::overlapping_aps ( Fragment & ref, int ref_aps, Fragment & fraglib, DOCKMol & oriented_mol, AMBER_TYPER & typer )
{ 
    Trace trace( "GA_Recomb::overlapping_aps()" );
    std::vector <std::pair <int, float> > overlapping;
    std::pair <int, float> hungarian_results;
 
    // For each aps of the reference
    for ( int i=0; i<fraglib.aps.size(); i++){

       // Declare rmsd and total number of heavy atoms
       float rmsd = 0.0;
       int heavy_total = 0;

       // Iterate through every atom in the ref
       // Then compute rmsd and add to running total
       // DELETE - REMOVE
       //cout << " a " << ref.mol.x[i] << " " << ref.mol.y[i] << " " << ref.mol.z[i] << " b " <<  oriented_mol.x[fraglib.aps[i].dummy_atom] << " " << oriented_mol.y[fraglib.aps[i].dummy_atom] << " " <<  oriented_mol.z[fraglib.aps[i].dummy_atom] <<endl;
       rmsd += (((ref.mol.x[ref.aps[ref_aps].dummy_atom] - oriented_mol.x[fraglib.aps[i].dummy_atom])*(ref.mol.x[ref.aps[ref_aps].dummy_atom] - oriented_mol.x[fraglib.aps[i].dummy_atom])) 
              + ((ref.mol.y[ref.aps[ref_aps].dummy_atom] - oriented_mol.y[fraglib.aps[i].dummy_atom])*(ref.mol.y[ref.aps[ref_aps].dummy_atom] - oriented_mol.y[fraglib.aps[i].dummy_atom])) 
              + ((ref.mol.z[ref.aps[ref_aps].dummy_atom] - oriented_mol.z[fraglib.aps[i].dummy_atom])*(ref.mol.z[ref.aps[ref_aps].dummy_atom] - oriented_mol.z[fraglib.aps[i].dummy_atom]))); 
        heavy_total++;

       // Make sure not to divide by 0, and compute final rmsd
       rmsd = sqrt(rmsd / (float) heavy_total);
       // DELETE - REMOVE
       //cout << "overlap rmsd " << rmsd << " ref " << ref.mol.atom_names[ref.aps[ref_aps].dummy_atom] << " frag " << oriented_mol.atom_names[fraglib.aps[i].dummy_atom]<<endl;

       // Hard coded - only need to check distance and bond types
       if ((rmsd < 0.5)) { 
          // Save pair information
          hungarian_results.first = i;
          hungarian_results.second = rmsd;
          overlapping.push_back(hungarian_results);
       }   
    }

    return overlapping;
}// end GA_Recomb::overlapping_aps()






// +++++++++++++++++++++++++++++++++++++++++
// Main replacement function that exchanges unique orients with the ref segment
// After each attachment, the torenv environment is computed and orients that do not pass are removed
std::vector <DOCKMol>
GA_Recomb::replace_scaffold ( Fragment & ref, int aps_num, int ref_seg, bool all_same_bt, bool all_diff_bt, 
                              Fragment & frag, std::vector <DOCKMol> & orients, 
                              AMBER_TYPER & typer, Master_Score & score )
{
    Trace trace( "GA_Recomb::replace_scaffold()" );
    std::vector <DOCKMol> recomb;
    std::vector <DOCKMol> mol_results;

    // STEP 1: Detemine how many unique orients are expected
    int expected = 0;  
    // (a) If no rings
    if ( frag.ring_size == 0 ){
       // If the bonds are all the same or all different, only one answer is important
       if ( (all_same_bt == true || all_diff_bt == true) ){ expected = 1; }
       else { expected = 2; }
    }
    // (b) If rings
    else{
       if ( all_diff_bt == true ){ expected = 1; }
       // Else try as many orients as possible
       else if ( ref.aps.size() == 2 ){ expected = 4; }
       else { expected = ref.aps.size(); }
    }

    // If two aps and a ring, calc norm
    bool check_norm = false; // If the ref aps are across from one another, then check frag
    std::vector <float> ref_norm; 
    // Calculate norm only if there are two aps and a ring
    // MARK: There may be an issue if there are two fused rings
    if ( (ref.aps.size() == 2) && (ref.ring_size > 0) ){
       // Calc cos between heavy atoms and norm is applicable
       ref_norm = calc_norm_ref ( ref, ref_seg ); 
       if ( ref_norm.size() > 0 ){ check_norm = true; }
    }

    //------------------------------------------------------------------------
    // STEP 2: Generate unique, allowed molecules through the attachment of 1 aps at a time
    int current = 0;  // Count number of orients that are valid 
    int invalid_torenv = 0; // Count the number of orients with invalid torenv

    // For each orient
    for ( int i=0; i<orients.size(); i++ ){
        // Skip move with bad hungarian
        bool skip = false;

        // (a) Pick a non-redundant fragment
        //cout << "NEW ORIENT " << i << " total num orients " << orients.size() << endl;//DELETE
        // Exit if at end of orients file
        if ( i == orients.size() ){ break; }
        // Exit if the orient is not valid
        if ( invalid_torenv == expected ){ break; }

        // Do not pick a redundant frag
        if ( i != 0 ){
           Hungarian_RMSD h;
           pair <double, int> result;

           // DELETE - REMOVE
           //fstream fout_molecules;
           //fout_molecules.open ( "zzz.hungarian_compare.mol2", fstream::out|fstream::app );
           //activate_mol(orients[i]);
           //fout_molecules << orients[i].current_data << endl;
           //Write_Mol2( orients[i], fout_molecules );
           //fout_molecules.close();

           // Compare with all previous orients used
           for ( int j=0; j<i; j++ ){
               // DELETE - REMOVE
               /*fout_molecules.open ( "zzz.hungarian_compare.mol2", fstream::out|fstream::app );
               activate_mol(orients[j]);
               fout_molecules << orients[j].current_data << endl;
               Write_Mol2( orients[j], fout_molecules );
               fout_molecules.close();*/
               //cout << "Entering orient clustering" << endl;
               
               //cout << "Comparing hungarian of orient " << i << " with previous " << j <<endl;
               result = h.calc_Hungarian_RMSD_dissimilar( orients[i], orients[j] );
               //cout << "H rmsd = " << result.first << " unmatched = " << result.second <<endl;
         
               // Hard-coded rmsd cutoff so that input parameters will not effect results
               if (result.first < 1.0) { 
                  skip = true;
                  //cout << "line: 6288 rmsd cutoff" << endl;
                  break; 
               } 
           }

        }

        // Skip orient, if similar
        if ( skip ) {continue; 
            cout << "Skipping orient because too similar" << endl;
        }

        // Skip the molecule if it does not lay in the plane with ref
        // MARK: There may be an issue if there are two fused rings
        float cos_ref = exit_vector (ref.mol, ref.aps[0].heavy_atom, ref.aps[0].dummy_atom, ref.aps[1].heavy_atom);
        float cos_frag = exit_vector (orients[i], frag.aps[0].heavy_atom, frag.aps[0].dummy_atom, frag.aps[1].heavy_atom);
        // DELETE cout << "The cos of the fragment " << cos_frag << " of the ref " << cos_ref <<endl;
        if ( (cos_ref < 0.99) && (cos_frag > 0.99) ){ cout << "WARNING: would not have checked the frag aps" <<endl;}
        if ( check_norm ){
           bool allowed  = compare_norm(ref, ref_seg, ref_norm, frag, orients[i]);
           if ( allowed == false ) { continue; }
        }
        //cout << "Number of orients after hungarian: " << orients
        // (b) Attach the fragment to the largest side
        // (1) Need to re-prepare the fragment to update the atom neighbors list for get_atom_children
        //cout << "Entering Amber typer" << endl;
        orients[i].prepare_molecule();

        // Work with a ref copy so that it is not translated or rotated
        DOCKMol tmp_ref;
        copy_molecule( tmp_ref, ref.mol);
        Fragment new_ref = mol_to_frag (tmp_ref);
        tmp_ref.clear_molecule();
  
        // (2) Identify which fragment aps is close to the current ref aps
        // For the pair, first is the aps index and second is the rmsd
        std::vector <std::pair <int, float> > overlapping_results = overlapping_aps ( ref, aps_num, 
                                                                         frag, orients[i], typer );
        // If there are no close matches 
        if ( overlapping_results.size() == 0 ){ continue; }

        // ID the fragment aps for attachment
        int frag_aps = compare_rmsd_bt ( ref, frag, all_same_bt, aps_num, overlapping_results );
        // Clear overlap vector
        overlapping_results.clear();

        // If there is no close aps, continue to next orient
        if ( frag_aps == -1 ){continue; }

        // (3) Activate the correct portions of the molecule
        activate_half_mol( new_ref.mol, ref.aps[aps_num].heavy_atom, ref.aps[aps_num].dummy_atom,
                           orients[i], frag.aps[frag_aps].heavy_atom, frag.aps[frag_aps].dummy_atom, true);

        // (4) Creates a new bond between the halves
        recomb_mols( new_ref.mol, ref.aps[aps_num].dummy_atom, ref.aps[aps_num].heavy_atom,
                     orients[i], frag.aps[frag_aps].heavy_atom, frag.aps[frag_aps].dummy_atom,recomb );

        // (5) Determine whether the new attachment is allowed
        // Generate the fragment
        Fragment f;
        Fragment tmp_growing = mol_to_frag (recomb[0]);
        Fragment growing = f.read_mol(tmp_growing);
        //cout << "Entering Amber typer" << endl;
        growing.mol.prepare_molecule();
        tmp_growing.mol.clear_molecule();

        // If using the torenv
        if ( ga_use_torenv_table ){
           // Prepare molecule for fingerprint/bond environ generation
           prepare_mol_torenv( growing.mol );
           //cout << "Entering Amber typer" << endl;
           typer.skip_verbose_flag = true;
           typer.prepare_molecule( growing.mol, true, score.use_chem, score.use_ph4, score.use_volume );

           DN_GA_Build c_dn;
           c_dn.dn_ga_flag = true;
           c_dn.dn_use_roulette = ga_use_dn_roulette;
           c_dn.dn_use_torenv_table = ga_use_torenv_table;
           c_dn.dn_torenv_table = ga_torenv_table; 
           c_dn.read_torenv_table(c_dn.dn_torenv_table);
           c_dn.verbose = 0;
           c_dn.dn_MW_cutoff_type_hard = true;
           c_dn.dn_MW_cutoff_type_soft = false;           

           // If not valid continue to next orient
           if ( c_dn.valid_torenv( growing ) == false ){ 
              growing.mol.clear_molecule();
              c_dn.torenv_vector.clear();
              growing.torenv_recheck_indices.clear();
              invalid_torenv++;
              cout << "invalid_torenv counter: " << invalid_torenv << endl;
              continue;
           }
        }

        //------------------------------------------------------------------------
        // STEP 3: Since the new molecule was renumbered, find the dummy atom & aps bond
        // (a) For a linker of 1 heavy atom only:
        if ( (ref.aps.size() == 2) ){
           // Need to use the aps of the orig ref that was not used
           for ( int j=0; j<ref.aps.size(); j++){
               // If the aps was not previously used
               if ( j != aps_num ){
                  // (1) Work with a ref copy so that it is not translated or rotated
                  DOCKMol tmp_ref;
                  copy_molecule( tmp_ref, ref.mol);
                  new_ref.mol.clear_molecule();
                  Fragment new_ref = mol_to_frag (tmp_ref);
                  tmp_ref.clear_molecule();

                  // (2) Activate the appropriate side of the growing ref & orig ref
                  activate_half_mol( growing.mol, growing.aps[0].dummy_atom, growing.aps[0].heavy_atom,
                                     new_ref.mol, ref.aps[j].dummy_atom, ref.aps[j].heavy_atom, true);

                  // (3) Creates a new bond between the halves
                  // In this case, we are making a bond between the DU of growing and heavy of tmp
                  // But, the dummy of growing has atom type Du, which has no radius
                  // Adjust the atom type to that atom type of the heavy atom of orig ref
                  growing.mol.atom_types[growing.aps[0].dummy_atom] = ref.mol.atom_types[ref.aps[j].heavy_atom];
                  recomb_mols( growing.mol, growing.aps[0].dummy_atom, growing.aps[0].heavy_atom,
                               new_ref.mol, ref.aps[j].heavy_atom, ref.aps[j].dummy_atom, recomb);


                  // (4) Check torenv if used
                  if ( ga_use_torenv_table ){
                     // Prepare the ref and convert to frag
                     //cout << "Entering Amber typer" << endl;
                     recomb[0].prepare_molecule();
                     typer.skip_verbose_flag = true;
                     typer.prepare_molecule( recomb[0], true, score.use_chem, score.use_ph4, score.use_volume );
                     prepare_mol_torenv( recomb[0] );
                     Fragment new_result = mol_to_frag (recomb[0]);

                     DN_GA_Build c_dn;
                     c_dn.dn_ga_flag = true;
                     c_dn.dn_use_roulette = ga_use_dn_roulette;
                     c_dn.dn_use_torenv_table = ga_use_torenv_table;
                     c_dn.dn_torenv_table = ga_torenv_table; 
                     c_dn.read_torenv_table(c_dn.dn_torenv_table);
                     c_dn.verbose = 0;
                     c_dn.dn_MW_cutoff_type_hard = true;
                     c_dn.dn_MW_cutoff_type_soft = false;

                     // If not valid, continue to next orient
                     if ( c_dn.valid_torenv( new_result ) == false ){ 
                        new_result.mol.clear_molecule();
                        c_dn.torenv_vector.clear();
                        new_result.torenv_recheck_indices.clear();
                        
                        invalid_torenv++;
                        cout << "invalid_torenv2: " << invalid_torenv << endl;
                        continue;
                     }
                  }
               }
           }
        }

        //------------------------------------------------------------------------
        // (b) For a scaffold:
        else {
           // (1) Recalc which aps of the ref matches with the growing molecule
           // because of the atom renumbering
           std::vector <std::pair <int, int> > final_matches;
           std::pair <int, int> tmp_pair;

           // For the remainder of the ref's aps
           for ( int j=0; j<ref.aps.size(); j++ ){
               // If the aps was not previously used
               if ( j !=aps_num ){
                  // Determine which frag aps is close
                  // For the pair, first is the aps index and second is the rmsd
                  std::vector <std::pair <int, float> > overlapping_results = overlapping_aps ( ref, j, 
                                                                           growing, recomb[0], typer );

                  // If there are no close matches 
                  if ( overlapping_results.size() == 0 ){ 
                     // Clear new molecule
                     recomb.clear();
                     growing.mol.clear_molecule();
                     break; 
                  }

                  // ID the closest atom to the reference
                  frag_aps = compare_rmsd_bt ( ref, growing, all_same_bt, j, overlapping_results );
                  overlapping_results.clear();
                  //cout << "overlapping " << overlapping_results.size() << " frag aps " << frag_aps <<endl;// DELETE
                  
                  // If the closest did not work, continue
                  if ( frag_aps == -1 ){ 
                     recomb.clear();
                     growing.mol.clear_molecule();
                     break; 
                  } 
                 
                  // Determine if the aps has already been used
                  bool found = false;
                  if ( final_matches.size() > 0 ){
                     for ( int k=0; k<final_matches.size(); k++ ){
                         if ( final_matches[k].second == frag_aps ){ 
                            recomb.clear();
                            growing.mol.clear_molecule();
                            found = true;
                            break; 
                         }
                     }
                  }
                 
                  // If already used break
                  if ( found == true ){ break; }

                  // Else, add pair to final list
                  tmp_pair.first = j;
                  tmp_pair.second = frag_aps;
                  final_matches.push_back(tmp_pair);
               } 
           }

           //------------------------------------------------------------------------
           // (2) Attach the remaining portion of the ref to the growing ref using final_matches 

           // If the aps matches are correct
           if ( final_matches.size() == ref.aps.size() ){
              // (2a) Make a vector of pairs for the dummy and heavy atoms that will 
              // be renumbered during attachment 
              std::vector < std::pair <int, int> > du_heavy;
              std::pair <int, int> tmp_du_heavy;

              // For the remaining aps in final_matches
              for ( int j=0; j<final_matches.size(); j++){
                  // Save index of DU and heavy atom
                  tmp_du_heavy.first = growing.aps[final_matches[j].second].dummy_atom;
                  tmp_du_heavy.second = growing.aps[final_matches[j].second].heavy_atom;
                  du_heavy.push_back(tmp_du_heavy);
              }

              // (2b) For each aps of the growing ref, attach a portion of ref
              for ( int j=1; j<growing.aps.size()+1; j++ ){
                  // Work with a ref copy so that it is not translated or rotated
                  DOCKMol tmp_ref;
                  copy_molecule( tmp_ref, ref.mol);
                  new_ref.mol.clear_molecule();
                  Fragment new_ref = mol_to_frag (tmp_ref);
                  tmp_ref.clear_molecule();

                  // First, prepare recomb mol each time for activation
                  //cout << "Entering Amber typer" << endl;
                  recomb[j-1].prepare_molecule();

                  // DELETE - REMOVE
                  /*cout << "CHECK DU HEAVY index " << du_heavy[j].first << " type" << recomb[j-1].atom_types[du_heavy[j].first] << " names " << recomb[j-1].atom_names[du_heavy[j].first] <<endl; 
                  cout << "CHECK DU HEAVY index " << du_heavy[j].second << " type" << recomb[j-1].atom_types[du_heavy[j].second] << " names " << recomb[j-1].atom_names[du_heavy[j].second] <<endl; 
                  for ( int k=0; k<recomb[j-1].num_bonds; k++ ){
                      if ( recomb[j-1].bonds_origin_atom[k] ==  du_heavy[j].first ){
                          cout << "found dummy heavy is " << recomb[j-1].bonds_target_atom[k] <<endl;
                      }
                      if ( recomb[j-1].bonds_target_atom[k] ==  du_heavy[j].first ){
                          cout << "found dummy heavy is " << recomb[j-1].bonds_origin_atom[k] <<endl;
                      }
                  }*/

                  // Second, activate the appropriate side of the growing ref & orig ref
                  activate_half_mol( recomb[j-1], du_heavy[j].first, du_heavy[j].second, new_ref.mol, 
                                     ref.aps[final_matches[j].first].dummy_atom, 
                                     ref.aps[final_matches[j].first].heavy_atom, true);

                   // Third, creates a new bond between the halves
                   // In this case, we are making a bond between the DU of growing and heavy of tmp
                   // But, the dummy of growing has atom type Du, which has no radius
                   // Adjust the atom type to that atom type of the heavy atom of orig ref
                   recomb[j-1].atom_types[du_heavy[j].first] = 
                                 ref.mol.atom_types[ref.aps[final_matches[j].first].heavy_atom];

                   // Fourth, save the DU and heavy atoms for the remaining aps
                   std::vector <std::pair <int, int> > remaining; 
                   if ( j+1 < growing.aps.size()+1 ){
                      std::pair <int, int> tmp_remain;
                      tmp_remain.first = du_heavy[j+1].first;
                      tmp_remain.second = du_heavy[j+1].second;
                      remaining.push_back(tmp_remain);
                   } 

                   // Fifth, create a new molecule and return the renumbered DU and heavy atoms
                   std::vector <std::pair <int, int> > new_du_heavy = 
                   replace_combine_mols( recomb[j-1], du_heavy[j].first, du_heavy[j].second, 
                                         new_ref.mol, ref.aps[final_matches[j].first].heavy_atom, 
                                         ref.aps[final_matches[j].first].dummy_atom, remaining, recomb);
                   // Clear previous atom index vector
                   remaining.clear();

                   // Sixth, if using torenv, check bond environment
                   if ( ga_use_torenv_table ){
                      //cout << "Entering Amber typer" << endl;
                      recomb[j].prepare_molecule();
                      typer.skip_verbose_flag = true;
                      typer.prepare_molecule( recomb[j], true, score.use_chem, score.use_ph4, score.use_volume );
                      prepare_mol_torenv( recomb[j] );
                      DN_GA_Build c_dn;
                      Fragment new_result = mol_to_frag (recomb[j]);
                      c_dn.dn_ga_flag = true;
                      c_dn.dn_use_roulette = ga_use_dn_roulette;
                      c_dn.dn_use_torenv_table = ga_use_torenv_table;
                      c_dn.dn_torenv_table = ga_torenv_table; 
                      c_dn.read_torenv_table(c_dn.dn_torenv_table);
                      c_dn.verbose = 0;
                      c_dn.dn_MW_cutoff_type_hard = true;
                      c_dn.dn_MW_cutoff_type_soft = false;

                      if ( c_dn.valid_torenv( new_result ) == false ){ 
                         new_result.mol.clear_molecule();
                         c_dn.torenv_vector.clear();
                        new_result.torenv_recheck_indices.clear();
                         invalid_torenv++;
                        continue;
                      }
                   }

                   // Update dummy/heavy atoms for growing ref, if applicable
                   // The list will include the current aps
                   //cout << "new du heavy size " << new_du_heavy.size() <<endl; //DELETE
                   if ( new_du_heavy.size() > 0 ){
                      // For the remaining aps
                      for ( int k=j+1; k<growing.aps.size()+1; k++ ){
                          du_heavy[k].first = new_du_heavy[k-(j+1)].first;
                          du_heavy[k].second = new_du_heavy[k-(j+1)].second;
                      }
                      new_du_heavy.clear();
                   }
                   // DELETE
                   //cout << "end for growing aps number " << j << endl;
              }
           }
           else{ continue; }
        }
        // DELETE
        /*fstream fout_molecules;
        fout_molecules.open ( "zzz.final_recomb.mol2", fstream::out|fstream::app );
        activate_mol(ref.mol);
        fout_molecules << ref.mol.current_data << endl;
        Write_Mol2( ref.mol, fout_molecules );
        
        activate_mol(orients[i]);
        fout_molecules << orients[i].current_data << endl;
        Write_Mol2( orients[i], fout_molecules );
      
        for ( int j=0; j<recomb.size(); j++){
        activate_mol(recomb[j]);
        fout_molecules << recomb[j].current_data << endl;
        Write_Mol2( recomb[j], fout_molecules );
        }
        fout_molecules.close();*/

        // Save the complete molecule
        mol_results.push_back(recomb[recomb.size()-1]);
        // Clear vectors/mols
        recomb.clear();
        growing.mol.clear_molecule();
        new_ref.mol.clear_molecule();

        // Update counter
        current++;
        // If current == expected, exit
        if ( current == expected ){ break; }
        // Exit if at end of orients file
        if ( i == orients.size() ){ break; }// DELETE
    }
    if (verbose) cout << "#### Found " << current << " of " << expected << " expected orients for a segment of size " << frag.mol.num_atoms << " with " << frag.aps.size() << "aps. Total orients " << orients.size()<<endl;
    if (verbose) cout << "#### Number with invalid torenv " << invalid_torenv <<endl;
    return mol_results;
}// end GA_Recomb::replace_scaffold





// +++++++++++++++++++++++++++++++++++++++++
// Calculates the angle between both heavy atoms of aps and one dummy
// Computes and returns the norm vector for the ref
// NOTE: this is only called when there are two aps
std::vector <float>
GA_Recomb::calc_norm_ref ( Fragment & ref, int seg ) 
{
   Trace trace( "GA_Recomb::calc_norm_ref()" );
   // Calculate the angle of the ref aps ( heavy atoms and one dummy atom )
   float cos_ref = exit_vector (ref.mol, ref.aps[0].heavy_atom, ref.aps[0].dummy_atom, ref.aps[1].heavy_atom);

   if ( abs(cos_ref) > 0.99 ){ 
      // Save ring atoms to vector
      std::vector <int> heavy_atoms;
      for ( int i=0; i<orig_segments[seg].num_atoms; i++ ){
          heavy_atoms.push_back(orig_segments[seg].atoms[i]);
      }
      std::vector <float> norm_crd = calc_norm( ref.mol, heavy_atoms ); 
      heavy_atoms.clear();
       return norm_crd;
   }
   // If not in a straight line, return empty vector
   else{ 
      std::vector <float> empty;
      return empty;
   }
} //end calc_norm_ref()




// +++++++++++++++++++++++++++++++++++++++++
// Computes and compares the norm of frag to ref only if the frag has has aps that are in a straight line
bool
GA_Recomb::compare_norm ( Fragment & ref, int seg, std::vector <float> ref_norm, Fragment & frag, DOCKMol & orient ) 
{
   Trace trace( "GA_Recomb::compare_norm()" );
   bool allowed = false;
   // Calculate the angle of the ref aps ( heavy atoms and one dummy atom )
   float cos_ref = exit_vector (ref.mol, ref.aps[0].heavy_atom, ref.aps[0].dummy_atom, ref.aps[1].heavy_atom);
   float cos_frag = exit_vector (orient, frag.aps[0].heavy_atom, frag.aps[0].dummy_atom, frag.aps[1].heavy_atom);

   double v_proj = 0;
   if ( abs(cos_frag) > 0.99 ){ 
      if ( abs(cos_ref) < 0.99 ) { cout << "WARNING>>> THIS IS A MISMATCH IN COS" <<endl; }
      // Include on the ring atoms for norm generations
      std::vector <int> heavy_atoms;
      for ( int i=0; i<orient.num_atoms; i++ ){
          if ( ( i != frag.aps[0].dummy_atom ) || ( i != frag.aps[1].dummy_atom ) ){
             heavy_atoms.push_back(i);
          }
      }
      std::vector <float> frag_norm = calc_norm(orient, heavy_atoms);
      heavy_atoms.clear();

      // Compare norm angles
      // v_proj cutoff, 60 deg:0.5, 45 deg:0.7071, 30 deg:0.8660.
      //double v_proj = ref_norm[0] * frag_norm[0]
      v_proj = ref_norm[0] * frag_norm[0]
                    + ref_norm[1] * frag_norm[1]
                    + ref_norm[2] * frag_norm[2];
      //DELETE cout << "v proj " << v_proj <<endl;
      if ( abs(v_proj) >= 0.7071 ){ // if greater than 45 deg (hard coded)
         allowed = true;
      }
   }
   else { allowed = true; }
   /*if ( allowed ){
        fstream fout_molecules;
        fout_molecules.open ( "zzz.cos_allowed.mol2", fstream::out|fstream::app );
        activate_mol(ref.mol);
        fout_molecules << DELIMITER << setw(STRING_WIDTH) << "V_proj:"
                       << setw(FLOAT_WIDTH) << abs(v_proj) <<endl;
        fout_molecules << ref.mol.current_data << endl;
        Write_Mol2( ref.mol, fout_molecules );
        fout_molecules << DELIMITER << setw(STRING_WIDTH) << "V_proj:"
                       << setw(FLOAT_WIDTH) << abs(v_proj) <<endl;
        fout_molecules << orient.current_data << endl;
        Write_Mol2( orient, fout_molecules );
        fout_molecules.close();
  }
   if ( !allowed ){
        fstream fout_molecules;
        fout_molecules.open ( "zzz.cos_notallowed.mol2", fstream::out|fstream::app );
        activate_mol(ref.mol);
        fout_molecules << DELIMITER << setw(STRING_WIDTH) << "V_proj:"
                       << setw(FLOAT_WIDTH) << abs(v_proj) <<endl;
        fout_molecules << ref.mol.current_data << endl;
        Write_Mol2( ref.mol, fout_molecules );
        fout_molecules << DELIMITER << setw(STRING_WIDTH) << "V_proj:"
                       << setw(FLOAT_WIDTH) << abs(v_proj) <<endl;
        fout_molecules << orient.current_data << endl;
        Write_Mol2( orient, fout_molecules );
        fout_molecules.close();
   }*/
   return allowed;
}// end calc_norm () 



// +++++++++++++++++++++++++++++++++++++++++
// Calculates the norm vector of a mol for the given number of atoms
std::vector <float>
GA_Recomb::calc_norm ( DOCKMol & mol, std::vector <int> heavy_atoms ) 
{
   Trace trace( "GA_Recomb::calc_norm()" );
   // Save the combined x, y, z coordinates for each atom of the ring
   std::vector <float> ring_center;
   for ( int i=0; i<3; i++){
       ring_center.push_back(0);
   }

   for ( int i=0; i<heavy_atoms.size(); i++ ){
       ring_center[0] += mol.x[heavy_atoms[i]];
       ring_center[1] += mol.y[heavy_atoms[i]];
       ring_center[2] += mol.z[heavy_atoms[i]];
   }


   // Calculate the center of the ring
   ring_center[0] /= float(heavy_atoms.size());
   ring_center[1] /= float(heavy_atoms.size());
   ring_center[2] /= float(heavy_atoms.size());

   // Compute the dot product (norm vector) of adjacent center-vertex vectors
   std::vector <float> norm_crd; // Will hold final normal vector
   std::vector <float> current_vertex; // vector from center of ring to current atom
   std::vector <float> next_vertex; // vector from ring center to next atom
   std::vector <float> current_norm; // norm vector of adjacent center-vertex vectors

   for ( int i=0; i<3; i++){
       norm_crd.push_back(0);
       current_vertex.push_back(0);
       next_vertex.push_back(0);
       current_norm.push_back(0);
   }
      
   for (int i=0; i<heavy_atoms.size()-1; i++) {
       // Save coordinates for current atom
       current_vertex[0] = mol.x[heavy_atoms[i]]-ring_center[0];
       current_vertex[1] = mol.y[heavy_atoms[i]]-ring_center[1];
       current_vertex[2] = mol.z[heavy_atoms[i]]-ring_center[2];

       // Save coordinates for adjacent atom
       next_vertex[0] = mol.x[heavy_atoms[i+1]]-ring_center[0];
       next_vertex[1] = mol.y[heavy_atoms[i+1]]-ring_center[1];
       next_vertex[2] = mol.z[heavy_atoms[i+1]]-ring_center[2];

       // Compute vector distance for pair (norm vector)
       current_norm[0] =(current_vertex[1]*next_vertex[2]) - (current_vertex[2]*next_vertex[1]);
       current_norm[1] =(current_vertex[2]*next_vertex[0]) - (current_vertex[0]*next_vertex[2]);
       current_norm[2] =(current_vertex[0]*next_vertex[1]) - (current_vertex[1]*next_vertex[0]);
 
       // Save to final norm vector
       norm_crd[0] += current_norm[0];
       norm_crd[1] += current_norm[1];
       norm_crd[2] += current_norm[2];
    }

    // Compute norm for last and first atom
    // Save coordinates for the last atom
    current_vertex[0] = mol.x[heavy_atoms[heavy_atoms.size()-1]]-ring_center[0];
    current_vertex[1] = mol.y[heavy_atoms[heavy_atoms.size()-1]]-ring_center[1];
    current_vertex[2] = mol.z[heavy_atoms[heavy_atoms.size()-1]]-ring_center[2];

    // Save coordinates for the first atom
    next_vertex[0] = mol.x[heavy_atoms[0]]-ring_center[0];
    next_vertex[1] = mol.y[heavy_atoms[0]]-ring_center[1];
    next_vertex[2] = mol.z[heavy_atoms[0]]-ring_center[2];

    // Compute vector distance for pair (norm vector)
    current_norm[0] =(current_vertex[1]*next_vertex[2]) - (current_vertex[2]*next_vertex[1]);
    current_norm[1] =(current_vertex[2]*next_vertex[0]) - (current_vertex[0]*next_vertex[2]);
    current_norm[2] =(current_vertex[0]*next_vertex[1]) - (current_vertex[1]*next_vertex[0]);

    // Save to final norm vector
    norm_crd[0] += current_norm[0];
    norm_crd[1] += current_norm[1];
    norm_crd[2] += current_norm[2];


    // Compute the average norm vector
    norm_crd[0] /= float (heavy_atoms.size());
    norm_crd[1] /= float (heavy_atoms.size());
    norm_crd[2] /= float (heavy_atoms.size());

    // Calculate distance of average norm vector
    float norm_radius = (norm_crd[0]*norm_crd[0]) + (norm_crd[1]*norm_crd[1]) + (norm_crd[2]*norm_crd[2]);
    norm_radius = sqrt(norm_radius);

    // Normalize the average norm vector by its length
    norm_crd[0] /= norm_radius;
    norm_crd[1] /= norm_radius;
    norm_crd[2] /= norm_radius;
    return norm_crd;
} //end GA_Recomb::calc_norm()



// +++++++++++++++++++++++++++++++++++++++++
// Extract the length and cosine of a bond between any two vectors
// In order to run this function, a mol2 and 3 atom numbers must be provided
float
GA_Recomb::exit_vector( DOCKMol & tmp_mol, int origin, int a, int b )
{
    Trace trace( "GA_Recomb::exit_vector()" );
    // Variables for exit_vector
    DOCKVector vecA;
    DOCKVector vecB;
    float vecA_mag;
    float vecB_mag;
    float dot;
    float cos_exit;
   
    // The length that will be used as a reference is between origin and a
    vecA.x = tmp_mol.x[origin] - tmp_mol.x[a];
    vecA.y = tmp_mol.y[origin] - tmp_mol.y[a];
    vecA.z = tmp_mol.z[origin] - tmp_mol.z[a];

    // Define distance vectors
    vecB.x = tmp_mol.x[origin] - tmp_mol.x[b];
    vecB.y = tmp_mol.y[origin] - tmp_mol.y[b];
    vecB.z = tmp_mol.z[origin] - tmp_mol.z[b];

    vecA_mag = ( vecA.x*vecA.x) + ( vecA.y*vecA.y) + ( vecA.z*vecA.z);
    vecB_mag = ( vecB.x*vecB.x) + ( vecB.y*vecB.y) + ( vecB.z*vecB.z);

    // Compute the dot product using the function in utils.cpp
    dot = dot_prod(vecA, vecB);

    // Compute cosine theta
    // If dot product is zero ( ie 90degrees) return 0
    if ( dot == 0 ){
       cos_exit = 0;
    } else{ cos_exit = dot / (sqrt (vecA_mag * vecB_mag)); }

    return cos_exit;
}//end GA_Recomb::exit_vector()



// +++++++++++++++++++++++++++++++++++++++++
// Return the index of the frag aps that overlaps with the given ref aps
// Also, takes into consideration the same bond types
int
GA_Recomb::compare_rmsd_bt ( Fragment & ref, Fragment & frag, bool all_same_bt, int aps_num,
                             std::vector <std::pair <int, float> > & overlapping_results )
{
    Trace trace( "GA_Recomb::compare_rmsd_bt()" );
    int frag_aps = -1;

    // If more than 1 aps is close to current dummy
    if ( overlapping_results.size() > 1 ){
       // Use the closer aps
       for ( int j=0; j<overlapping_results.size()-1; j++ ){

           if ( overlapping_results[j].second > overlapping_results[j+1].second ){
              frag_aps = overlapping_results[j+1].first;

              // If the bond types are different & the ref/frag bond types do not match
              if ( all_same_bt == false ){
                 if (ref.aps_bonds_type[aps_num].second != frag.aps_bonds_type[frag_aps].second ){
                    frag_aps = -1;
                    continue;  
                 }
              }
           }
           else { 
              frag_aps = overlapping_results[j].first; 
              // If the bond types are different & the ref/frag bond types do not match
              if ( all_same_bt == false ){
                 //DELETE cout << " ref bt " << ref.aps_bonds_type[aps_num].second << " frag " << frag.aps_bonds_type[frag_aps].second <<endl;
                 if ( ref.aps_bonds_type[aps_num].second != frag.aps_bonds_type[frag_aps].second ){
                    frag_aps = -1;
                    continue;    
                 }
              }
           }
        }
     }

     // Else, if only 1 orientation, use the one aps present
     else {
        frag_aps = overlapping_results[0].first;
        if ( all_same_bt == false ){
           if ( ref.aps_bonds_type[aps_num].second != frag.aps_bonds_type[frag_aps].second ){
               frag_aps = -1;
           }
        } 
     }
     return frag_aps;
} // end GA_Recomb::compare_rmsd_bt()




// +++++++++++++++++++++++++++++++++++++++++
// Convert a H to a Du to generate an anchor for DN design
// and returns true if a H was found
// NOTE: Picks a H at random from all present in the segment
bool
GA_Recomb::rand_H_to_Du ( DOCKMol & mol, int segment )
{
    Trace trace( "GA_Recomb::rand_H_to_Du()" );
    // Bool used to exit loop once the appropriate H is found
    bool found = false;
    // Save the H atom index
    INTVec h_index; 

    // For each atom on the molecular segment
    for (int i=0; i<orig_segments[segment].num_atoms; i++){
        
        // make sure H isn't a DNM moiety LEP
        if ( mol.atom_names[orig_segments[segment].atoms[i]] != "DNM"){
        // Find the hydrogens and save the index  
            if ( mol.atom_types[orig_segments[segment].atoms[i]] == "H" ){
                h_index.push_back(orig_segments[segment].atoms[i]);
            }
        }
    }
   
    // If there are hydrogen atoms
    if ( h_index.size() > 0 ){
       // If there is more than 1
       if ( h_index.size() > 1 ){
          // Randomly pick one of the hydrogens
          int sel_index=rand() % h_index.size();

          // Change atom type and charge 
          found = true;
          cout << "H is being turned into a Du: " << mol.atom_names[h_index[sel_index]] << endl;
          H_to_Du( mol, h_index[sel_index]);
       }

       // If only one H
       else{
          // Change atom type and charge 
          found = true;
          cout << "H is being turned into a Du: " << mol.atom_names[h_index[0]] << endl;
          H_to_Du( mol, h_index[0]);
       }
    }
    // Else, there are no hydrogens
    else {
      // Print warning - DELETE
      if (verbose) cout << "NOTE: At present the GA code is unable to add a fragment to a segment "
           << "that does not contain a hydrogen atom." <<endl;
      found = false;
    }

    return found;
}//end GA_Recomb::rand_H_to_Du()




// +++++++++++++++++++++++++++++++++++++++++
// Converts a specific H atom to Du
void
GA_Recomb::H_to_Du ( DOCKMol & mol, int a )
{
    Trace trace( "GA_Recomb::H_to_Du()" );
    /*fstream fout_molecules;
    fout_molecules.open ( "zzz.du.mol2", fstream::out|fstream::app );
    fout_molecules <<mol.current_data << endl;
    Write_Mol2(mol, fout_molecules);
    fout_molecules.close();*/
    // Replace wih Du - regardless of bond order
    mol.atom_types[a] = "Du";

    // Also change the paritial charge of Du to 0.0
    mol.charges[a] = 0.0;

    stringstream ss;
    ss << "Du1";
    mol.atom_names[a] = ss.str();
    ss.clear();

    // DELETE
    //fstream fout_molecules;
    /*fout_molecules.open ( "zzz.du.mol2", fstream::out|fstream::app );
    fout_molecules <<mol.current_data << endl;
    Write_Mol2(mol, fout_molecules);
    fout_molecules.close();*/

    return;
}// end GA_Recomb::H_to_Du()





// +++++++++++++++++++++++++++++++++++++++++
// Generates a new DOCKMol where a molecule is capped with a hydrogen
// In order to run this function, a mol2 and 3 atom numbers must be provided
void
GA_Recomb::add_H( DOCKMol & mol, int origin, int a, bool subst,  Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & typer)
{
    Trace trace( "GA_Recomb::add_H()" );
 
    // STEP 1: Change the atom from delete to H and correct bond length
    // Activate this atom
    mol.atom_active_flags[a] = true;
    mol.atom_types[a] = "H";
    mol.atom_names[a] = "Hd";

    // DELETE
    //fstream fout_molecules;
    //fout_molecules.open ( "zzz.predel.mol2", fstream::out|fstream::app );
    //fout_molecules << mol.current_data << endl;
    //Write_Mol2( mol, fout_molecules );
    //fout_molecules.close();

    // Re-index the names of all Hydrogen atoms
    // (this can be done later)
    
    // Calculate the desired bond length and remember as 'new_rad'
    float new_rad = calc_cov_radius(mol.atom_types[origin]) +
    calc_cov_radius(mol.atom_types[a]);
    
    // Calculate the x-y-z components of the bond vector (bond_vec)
    DOCKVector bond_vec;
    bond_vec.x = mol.x[origin] - mol.x[a];
    bond_vec.y = mol.y[origin] - mol.y[a];
    bond_vec.z = mol.z[origin] - mol.z[a];
    
    // Normalize the bond vector then multiply each component by new_rad so that it is the desired length
    bond_vec = bond_vec.normalize_vector();
    bond_vec.x *= new_rad;
    bond_vec.y *= new_rad;
    bond_vec.z *= new_rad;
    
    // Change the coordinates of the Hydrogen atom so that the bond length is correct
    mol.x[a] = mol.x[origin] - bond_vec.x;
    mol.y[a] = mol.y[origin] - bond_vec.y;
    mol.z[a] = mol.z[origin] - bond_vec.z;
  

    //------------------------------------------------------------------------
    // STEP 2: Update active bond flags and keep track of the bond associated 
    // with the H to ensure addition is allowed

    // Save bond number for activation later
    int to_activate;

    for ( int i=0; i<mol.num_bonds; i++){

        // (a) Update the bond active flags using the keep bonds flags
        mol.bond_active_flags[i] = false;
        if (mol.bond_keep_flags[i] == 1){
           mol.bond_active_flags[i] = true;
        }

        // (b) Activate the bond associated with the H
        if ( mol.bonds_origin_atom[i] == a ){
            if ( mol.atom_active_flags[mol.bonds_origin_atom[i]] == 1 &&  mol.atom_active_flags[mol.bonds_target_atom[i]] == 1){
               // If the bond is not a single bond & not subst, exit
               if ( subst == false) {
                  if (mol.bond_types[i] != "1"){
                     activate_mol( mol);
                     // Deactivate specific bond
                     mol.bond_active_flags[i] = false;

                     // DELETE
                     /*cout << "WARNING: deletion was not conducted - bond type " << mol.bond_types[i] <<endl;
                     fstream fout_molecules;
                     fout_molecules.open ( "zzz.faulty-deletion.mol2", fstream::out|fstream::app );
                     fout_molecules << "########## Deletion failed at bond number " << i << " bond type " << mol.bond_types[i] <<endl;
                     fout_molecules << mol.current_data << endl;
                     Write_Mol2( mol, fout_molecules );
                     fout_molecules.close();*/
                     return;
                  }
                  // else if single bond save bond index
                  else{
                     mol.bond_active_flags[i] = true;

                     //Save bond index
                     to_activate = i;
    
                     // DELETE
                     //cout << "orig bond " << i << " orig " << mol.bonds_origin_atom[i] << " target " << mol.bonds_target_atom[i]  <<endl;
                  }
               }
               // If subs
               else if ( subst ){
                  mol.bond_active_flags[i] = true;

                  //Save bond index
                  to_activate = i;
    
                  // DELETE
                  //cout << "orig bond " << i << " orig " << mol.bonds_origin_atom[i] << " target " << mol.bonds_target_atom[i]  <<endl;
               }
            }
        }    

        // (c) Activate the bond associated with the H if target
        else if ( mol.bonds_target_atom[i] == a ){
            if ( mol.atom_active_flags[mol.bonds_origin_atom[i]] == 1 &&  mol.atom_active_flags[mol.bonds_target_atom[i]] == 1){
               // If the bond is not a single bond, exit
               if (subst == false){
                  if (mol.bond_types[i] != "1"){
                     activate_mol( mol);
                     // Deactivate specific bond
                     mol.bond_active_flags[i] = false;

                     // DELETE
                     cout << "WARNING: deletion was not conducted - bond type " << mol.bond_types[i] <<endl;
                     /*
		     fstream fout_molecules;
                     fout_molecules.open ( "zzz.faulty-deletion.mol2", fstream::out|fstream::app );
                     fout_molecules << "########## Deletion failed at bond number " << i << " bond type " << mol.bond_types[i] <<endl;
                     fout_molecules << mol.current_data << endl;
                     Write_Mol2( mol, fout_molecules );
                     fout_molecules.close();*/
                     return;
                  }
                  // else if single bond save bond index
                  else{
                     mol.bond_active_flags[i] = true;

                     //Save bond index
                     to_activate = i;
    
                     // DELETE
                     //cout << "orig bond " << i << " orig " << mol.bonds_origin_atom[i] << " target " << mol.bonds_target_atom[i]  <<endl;
                  }
               } 
               else if ( subst ){
                  mol.bond_active_flags[i] = true;
                  //Save information to intvec
                  to_activate = i;
             
                  // DELETE
                  //cout << "orig bond " << i << " orig " << mol.bonds_origin_atom[i] << " target " << mol.bonds_target_atom[i]  <<endl;
               }
            }  
        }    
    }
 
    //------------------------------------------------------------------------
    // STEP 3: Create a new DOCKMol
    DOCKMol update;

    int natoms = 0;
    int nbonds = 0;

    for( int i=0; i< mol.num_atoms; i++){
       if( mol.atom_active_flags[i] == true){
          natoms++;
       } 
    }

    for( int i=0; i< mol.num_bonds; i++){
       if( mol.bond_active_flags[i] == true){
          nbonds++;
       } 
    }

    update.allocate_arrays(natoms, nbonds, 1);

    // Update name to include deletion
    //ostringstream new_energy;

    //new_energy << mol.energy << "-del";
    //update.energy    = new_energy.str();
    update.comment1 = mol.comment1;
    update.comment2 = mol.comment2;
    update.comment3 = mol.comment3;
    update.energy   = mol.energy;

    update.num_atoms = natoms;
    update.num_bonds = nbonds;

    if (dnm_enabled){
        update.mol_dnm_flag = mol.mol_dnm_flag;
    }
    


    int atom_index = 0;
    INTVec renumber;
    renumber.resize( mol.num_atoms );

    for( int i=0; i< mol.num_atoms; i++){
       if ( mol.atom_active_flags[i] == 1 ){
          update.x[atom_index] = mol.x[i];
          update.y[atom_index] = mol.y[i];
          update.z[atom_index] = mol.z[i];

          update.charges[atom_index]                = mol.charges[i];
          update.atom_names[atom_index]             = mol.atom_names[i];
          update.atom_types[atom_index]             = mol.atom_types[i];
          update.atom_number[atom_index]            = mol.atom_number[i];
          update.atom_residue_numbers[atom_index]   = mol.atom_residue_numbers[i]; 
          update.subst_names[atom_index]            = mol.subst_names[i];
          update.flag_acceptor[atom_index]          = mol.flag_acceptor[i];
          update.flag_donator[atom_index]           = mol.flag_donator[i];
          update.acc_heavy_atomid[atom_index]       = mol.acc_heavy_atomid[i];

          if (dnm_enabled){
            update.atom_dnm_flag[atom_index]          = mol.atom_dnm_flag[i];
          }
          

          renumber[i] = atom_index;
          atom_index++;
       }
    }        
   
    int bond_index = 0;
    for (int i=0; i<mol.num_bonds; i++){
         if ( mol.bond_active_flags[i] == true ){
             update.bonds_origin_atom[bond_index] = renumber[mol.bonds_origin_atom[i]];
             update.bonds_target_atom[bond_index] = renumber[mol.bonds_target_atom[i]];
             update.bond_types[bond_index]        = mol.bond_types[i];
             update.bond_ring_flags[bond_index]   = mol.bond_ring_flags[i];

             // Update the index of the bond associated with the H atom
             if ( i == to_activate ){
                   to_activate = bond_index;
             }
             bond_index++;
         }
    }

    
    // Prepare molecule to repopulate rot bond flags
    activate_mol(update);
    //cout << "Entering Amber typer" << endl;
    update.prepare_molecule();
    typer.skip_verbose_flag = true;
    typer.prepare_molecule( update, true, score.use_chem, score.use_ph4, score.use_volume );

    // DELETE - REMOVE
    //fstream fout_molecules;
    //fout_molecules.open ( "zzz.deletion.mol2", fstream::out|fstream::app );
    //fout_molecules << update.current_data << endl;
    //Write_Mol2( update, fout_molecules );
    //fout_molecules.close();


    //------------------------------------------------------------------------
    // STEP 4: Remove molecules that are no longer valid, unless completing subst
    // (a) If doing substitution, make anchor and exit
    if ( subst == true ){
       H_to_Du (update, renumber[a] );
       mol.clear_molecule();
       mutants.clear();
       mutants.push_back(update);
       activate_mol(mutants[0]);
       return;
    }

    // (b) If using the torenv table, id adjacent rotatable bonds near new H bond
    // NOTE: The bond saved above, in variable to_activate, is a bond between a heavy atom and H
    // Since the torenv does not include hydrogens, the adjacent heavy - heavy atom bonds will be checked 
    
    if ( (ga_use_torenv_table == true) ){
       // Save adjacent rot bonds
       INTVec heavy_nbrs;

       // Find adjacent bonds but save only rot bonds 
       if ( update.bonds_origin_atom[to_activate] != renumber[a] ){
          vector <int> nbrs = update.get_bond_neighbors(update.bonds_origin_atom[to_activate]);
          for (int i=0; i<nbrs.size(); i++){
              // If rot 
              if (update.amber_bt_id[nbrs[i]] != -1){
                 heavy_nbrs.push_back(nbrs[i]);
              }
          }
       }

       // Find adjacent bonds but save only rot bonds 
       else if ( update.bonds_target_atom[to_activate] != renumber[a] ){
          vector <int> nbrs = update.get_bond_neighbors(update.bonds_target_atom[to_activate]);
          for (int i=0; i<nbrs.size(); i++){
              if (update.amber_bt_id[nbrs[i]] != -1 ){
                 heavy_nbrs.push_back(nbrs[i]);
              }
          }
       }

       // (c) If there are rotatable adjacent bonds, check the torenv
       //cout << "heavy atoms size " << heavy_nbrs.size() <<endl;// DELETE
       if ( heavy_nbrs.size() != 0 ){
          //cout << "heavy_nbrs.size" << endl;
          DN_GA_Build c_dn; 
          // Prepare molecule for fingerprint/bond environ generation
          int valid = 0;
          Fragment tmp_frag = mol_to_frag(update);
          c_dn.dn_ga_flag = true;
          c_dn.dn_use_roulette = ga_use_dn_roulette;
          c_dn.dn_use_torenv_table = ga_use_torenv_table;
          c_dn.dn_torenv_table = ga_torenv_table; 
          c_dn.read_torenv_table(c_dn.dn_torenv_table);
          c_dn.verbose = 0;
          c_dn.dn_MW_cutoff_type_hard = true;
          c_dn.dn_MW_cutoff_type_soft = false;


          // Check each adjacent rotbonds torenv 
          for (int i=0; i<heavy_nbrs.size(); i++){
              prepare_mut_torenv( update, heavy_nbrs[i] );

              // If the new mol has valid bonds, increment counter
              if ( c_dn.valid_torenv( tmp_frag )){ valid++; }
              // DELETE
              else{
                 if (verbose) cout << "#### The bond that did not pass was: " << heavy_nbrs[i] <<endl;
                 if (verbose) cout << "#### The orig atom: " << update.bonds_origin_atom[heavy_nbrs[i]] << " target " << update.bonds_target_atom[heavy_nbrs[i]] <<endl;
              }
          }

          // If each of the adjacent bonds is valid, check all of the bonds
          if ( valid == heavy_nbrs.size()){
             //cout << "valid == heavy_nbrs.size" << endl;
             prepare_torenv_indices(tmp_frag);
             if ( c_dn.valid_torenv_multi( tmp_frag ) ){
                cout << "Passed multi torenv check" << endl;
                // If valid, minimize & save to mutants
                minimize_children(update, mutants, score, simplex, typer);
                activate_mol(mutants[0]);

                // DELETE - REMOVE
                //fstream fout_molecules;
                //fout_molecules.open ( "zzz.valid-deletion.mol2", fstream::out|fstream::app );
                //fout_molecules << mutants[0].current_data << endl;
                //Write_Mol2( mutants[0], fout_molecules );
                //fout_molecules.close();
             }
          }
          //Clear vectors for torenv
          c_dn.torenv_vector.clear();
          tmp_frag.torenv_recheck_indices.clear();
          heavy_nbrs.size();
       }//end if there are heavy_nbrs()

       // (d) If no heavy rot atoms, check all bonds
       else{
          DN_GA_Build c_dn;
          Fragment tmp_frag = mol_to_frag(update);
          c_dn.dn_ga_flag = true;
          c_dn.dn_use_roulette = ga_use_dn_roulette;
          c_dn.dn_use_torenv_table = ga_use_torenv_table;
          c_dn.dn_torenv_table = ga_torenv_table; 
          c_dn.read_torenv_table(c_dn.dn_torenv_table);
          c_dn.verbose = 0;
          c_dn.dn_MW_cutoff_type_hard = true;
          c_dn.dn_MW_cutoff_type_soft = false;

          prepare_torenv_indices(tmp_frag);
          if ( c_dn.valid_torenv_multi( tmp_frag ) ){
             // If valid, minimize & save to mutants
             minimize_children(update, mutants, score, simplex, typer);
 
             // DELETE - REMOVE
             //fout_molecules.open ( "zzz.noheavy.mol2", fstream::out|fstream::app );
             //fout_molecules << update.current_data << endl;
             //Write_Mol2( update, fout_molecules );
             //fout_molecules.close();
          }
          //Clear vectors for torenv
          c_dn.torenv_vector.clear();
          tmp_frag.torenv_recheck_indices.clear();
       }
    }
    // (e) If not using torenv table
    else{
       minimize_children(update, mutants, score, simplex, typer);
    }

    // Clear update/ vectors
    update.clear_molecule();

    // If there are mutants DELETE
    /*if (mutants.size() > 0 ){
       fstream fout_molecules;
       fout_molecules.open ( "zzz.successful-deletion.mol2", fstream::out|fstream::app );
       fout_molecules << mutants[0].current_data << endl;
       Write_Mol2( mutants[0], fout_molecules );
       fout_molecules.close();
    }*/
    return;
    
} //end GA_Recomb::add_H()




            //////////////////////////////////////////
            //****ENERGY MINIMIZATION FUNCTIONS*****//
            /////////////////////////////////////////
            



// +++++++++++++++++++++++++++++++++++++++++
// Calculate the internal energy of molecule
void
GA_Recomb::prepare_internal_energy( DOCKMol & tmp_child, Master_Score & score )
{
    Trace trace( "GA_Recomb::prepare_internal_energy()" );

    // Pass the method ('2') to base_score
    score.primary_score->method = 2;

    // If using internal energy, initialize it here
    if ( score.use_primary_score ) {

        // Initialize internal energy
        score.primary_score->use_internal_energy = use_internal_energy;

        if ( use_internal_energy ) {
            // Pass vdw parameters to base_score
            score.primary_score->ie_att_exp = ie_att_exp;
            score.primary_score->ie_rep_exp = ie_rep_exp;
            score.primary_score->ie_diel = ie_diel;

            // Need to use DOCKMol with radii and segments assigned
            // it does not matter which atoms are labeled active
            score.primary_score->initialize_internal_energy(tmp_child);
        }
    }

    return;
} // end GA_Recomb::prepare_internal_energy()





// +++++++++++++++++++++++++++++++++++++++++
// Given a fragment and an empty vector of fragments, sample torsions for the most recent addition
// and return the results in the vector. NOTE: DOES NOT MINIMIZE TORSIONS
void
GA_Recomb::minimize_children( DOCKMol & child, std::vector <DOCKMol> & molvec, Master_Score & score,Simplex_Minimizer & simplex, AMBER_TYPER & typer )
{
    Trace trace( "GA_Recomb::minimize_children()" );
    // Activate all atoms and bonds prior to any scoring
    activate_mol( child );

    // This will be true if energy is calculated correctly
    bool valid_orient = false;

    // Declare a temporary fingerprint object and compute atom environments (necessary for Gasteiger)
    Fingerprint temp_finger;
    for ( int i=0; i<child.num_atoms; i++ ){
        child.atom_envs[i] = temp_finger.return_environment( child, i );
    }

    // Prepare the molecule using the amber_typer to assign bond types to each bond
    //cout << "Entering Amber typer" << endl;
    child.prepare_molecule();
    typer.skip_verbose_flag = true;
    typer.prepare_molecule( child, true, score.use_chem, score.use_ph4, score.use_volume );

    // Compute the charges, saving them on the mol object
    float total_charges = 0;
    total_charges = compute_gast_charges(child);

    // Initially reset all of the values in bond_tors_vectors to -1 (bond direction does not matter)
    bond_tors_vectors.clear();
    bond_tors_vectors.resize( child.num_bonds, -1 );

    // Prepare for internal energy calculation for all conformers 
    if (use_internal_energy){ prepare_internal_energy( child, score ); }

    // Flexible minimization
    simplex.minimize_flexible_growth( child, score, bond_tors_vectors );
       
    // Compute internal energy and primary score, store it in the dockmol
    if ( score.use_primary_score ) {
       valid_orient = score.compute_primary_score( child );
         
    } else {
        valid_orient = true;
        child.current_score = 0;
        child.internal_energy = 0;
    }

    // If the energy or if the internal energy is greater than the energy_cutoff (100 is default)
    if ( valid_orient ){
        molvec.push_back( child );
    }
    else { cout <<"WARNING: Could not minimize structure" << endl;}

    bond_tors_vectors.clear();
    return;

} // end GA_Recomb::minimize_children()





            //////////////////////////////////////////
            //****FITNESS FUNCTIONS (PRUNING)*****//
            /////////////////////////////////////////



// +++++++++++++++++++++++++++++++++++++++++
// // Prune children by similarity with parents DOCKMol vec using Hungarian Score
// // If there are redundant molecules, keep the one with the best score between the
// // children and parents
void GA_Recomb::uniqueness_prune_mut( std::vector <DOCKMol> & children, std::vector <DOCKMol> & parents, Master_Score & score, AMBER_TYPER & typer)
{
    Trace trace( "GA_Recomb::uniqueness_prune_mut()" );
    double start_time = time_seconds();
    std::vector <DOCKMol> saved;

    // These children need to be prepared with amber typer before rmsd calculation
    for ( int i=0; i<children.size(); i++ ){
         //cout << "Entering Amber typer" << endl;
         typer.skip_verbose_flag = true;
         typer.prepare_molecule(children[i], true, score.use_chem, score.use_ph4, score.use_volume);
         children[i].used = false;
    }

    // Activate parents, there were deactivated at the top of breeding()
    activate_vector( parents);
    for ( int i=0; i<parents.size(); i++ ){
         parents[i].used = false;
    }
    
    // Integers to count the number of compounds pruned at each step
    int initial_parents = 0; // number the same as parents
    int hungarian_prune_cc = 0; // number of children that are the same
    
    Hungarian_RMSD h;
    pair <double, int> result;

    // If the two mols do not have the same number of atoms, rmsd is -1000
    double torsion_rmsd_cutoff = 2.0; //MARK: user-specified?

    // Compare offspring to parents
    for ( int i=0; i<parents.size(); i++ ){
        // For each molecule listed below it
        for ( int j=0; j<children.size(); j++ ){
            // If the below molecule has not been marked used
            if ( !children[j].used ){
               // Compare score of parents and offspring
               if (parents[i].current_score < children[j].current_score){
                  // Calculate and save hunarian rmsd
                  result = h.calc_Hungarian_RMSD_dissimilar( parents[i], children[j] );
                  //cout << "calc H rmsd = " << result.first <<endl;
                  //cout << "num unmatched = " << result.second << endl;
                  // If the two confomers are very similar, mark the one with higher energy as used
                  // If the two molecules do not have the same number of atoms than rmsd is -1000
                  if (result.first < ga_heur_matched_rmsd && result.second < ga_heur_unmatched_num) {
                     children[j].used = true;
                  }
               }
               // If child score is better
               else{
                   result = h.calc_Hungarian_RMSD_dissimilar(children[j], parents[i] );
                   //cout << "calc H rmsd = " << result.first <<endl;
                   //cout << "num unmatched = " << result.second << endl;

                  if (result.first < ga_heur_matched_rmsd && result.second < ga_heur_unmatched_num) {
                     // Check to see if parent has already been used!
                     if (!parents[i].used){
                        parents[i].used = true;
                     }
                     // if Parent used already make child as used
                     else{
                        children[j].used = true;
                     }
                  }

               }
            }
        }
    }
    // At this point, only cluster heads will not be flagged 'used', add those to scored_generation
    for (int i=0; i<children.size(); i++){
        if (!children[i].used){
           saved.push_back(children[i]);
        }
    }

    // DELETE
    cout <<"### Offspring_size_before_parents_pruning " << children.size();
    children.clear();
    for (int i=0; i<saved.size(); i++){
       children.push_back(saved[i]);
    }
    cout << " offspring_size_after_pruning " << children.size() << endl;
    saved.clear();

    for (int i=0; i<parents.size();i++){
       if (!parents[i].used){
          saved.push_back(parents[i]);
       }
    }
    cout << "### Parents size before P-O pruning " << parents.size();
    parents.clear();

    for (int i=0; i<saved.size(); i++){
       parents.push_back(saved[i]);
    }
    saved.clear();
    cout << " Parents size after pruning " << parents.size() << endl;

    double stop_time = time_seconds();
    cout << "\n" "-----------------------------------" "\n";
    if (verbose) cout << " Elapsed time P-O pruning:\t" << stop_time - start_time
         << " seconds for GA\n\n";

return;
}//end GA_Recomb::uniqueness_prune_mut()

// +++++++++++++++++++++++++++++++++++++++++
// Prune children by similarity with parents DOCKMol vec using Hungarian Score 
// If there are redundant molecules, keep the one with the best score between the
// children and parents
void
GA_Recomb::uniqueness_prune( std::vector <DOCKMol> & children, std::vector <DOCKMol> & parents, Master_Score & score, AMBER_TYPER & typer)
{
    Trace trace( "GA_Recomb::uniqueness_prune()" );

    double start_time = time_seconds();
    std::vector <DOCKMol> saved;

    // These children need to be prepared with amber typer before rmsd calculation
    for ( int i=0; i<children.size(); i++ ){
         //cout << "Entering Amber typer" << endl;
         typer.skip_verbose_flag = true;
         typer.prepare_molecule(children[i], true, score.use_chem, score.use_ph4, score.use_volume);
         children[i].used = false;
         /*fstream fout_molecules;
         fout_molecules.open ( "zzz.kids_for_rob.mol2", fstream::out|fstream::app );
         fout_molecules <<children[i].current_data << endl;
         Write_Mol2(children[i], fout_molecules);
         fout_molecules.close();*/
    
    }
  
    // Activate parents, there were deactivated at the top of breeding()
    activate_vector( parents);
    for ( int i=0; i<parents.size(); i++ ){
         parents[i].used = false;
    }

    // Integers to count the number of compounds pruned at each step
    int initial_parents = 0; // number the same as parents
    int hungarian_prune_cc = 0; // number of children that are the same
     
    Hungarian_RMSD h; 
    pair <double, int> result;

    // If the two mols do not have the same number of atoms, rmsd is -1000
    double torsion_rmsd_cutoff = 2.0; //MARK: user-specified?
    
    // Compare offspring to parents
    for ( int i=0; i<parents.size(); i++ ){
        // For each molecule listed below it
        for ( int j=0; j<children.size(); j++ ){
            // If the below molecule has not been marked used
            if ( !children[j].used ){
               // Compare score of parents and offspring
               if (parents[i].current_score < children[j].current_score){
                  // Calculate and save hunarian rmsd
                  result = h.calc_Hungarian_RMSD_dissimilar( parents[i], children[j] );
                  //cout << "calc H rmsd = " << result.first <<endl;
                  //cout << "num unmatched = " << result.second << endl;

                  // If the two confomers are very similar, mark the one with higher energy as used
                  // If the two molecules do not have the same number of atoms than rmsd is -1000
                  if (result.first < ga_heur_matched_rmsd && result.second < ga_heur_unmatched_num) {
                     children[j].used = true;
                  }
               }
               // If child score is better
               else{
                   result = h.calc_Hungarian_RMSD_dissimilar(children[j], parents[i] );
                   //cout << "calc H rmsd = " << result.first <<endl;
                   //cout << "num unmatched = " << result.second << endl;

                  if (result.first < ga_heur_matched_rmsd && result.second < ga_heur_unmatched_num) {
                     // Check to see if parent has already been used!
                     if (!parents[i].used){
                        parents[i].used = true;
                     }
                     // if Parent used already make child as used
                     else{
                        children[j].used = true;
                     }
                  }

               }
            }
        }
    }
    // At this point, only cluster heads will not be flagged 'used', add those to scored_generation
    for (int i=0; i<children.size(); i++){
        if (!children[i].used){
           saved.push_back(children[i]);
        }
    }

    // DELETE
    cout <<"### Offspring_size_before_parents_pruning " << children.size();
    children.clear();
    for (int i=0; i<saved.size(); i++){
       children.push_back(saved[i]);
    }
    cout << " offspring_size_after_pruning " << children.size() << endl;
    saved.clear();

    for (int i=0; i<parents.size();i++){
       if (!parents[i].used){
          saved.push_back(parents[i]);
       }
    }
    cout << "### Parents size before P-O pruning " << parents.size();
    parents.clear();

    for (int i=0; i<saved.size(); i++){
       parents.push_back(saved[i]);
    }
    saved.clear();
    cout << " Parents size after pruning " << parents.size() << endl;

    double stop_time = time_seconds();
    cout << "\n" "-----------------------------------" "\n";
    if (verbose) cout << " Elapsed time P-O pruning:\t" << stop_time - start_time
         << " seconds for GA\n\n";

return;
}//end GA_Recomb::uniqueness_prune()




// +++++++++++++++++++++++++++++++++++++++++
// Fitness function that includes pruning by internal energy, scoring functions (as determined by dock input file)
// and by compounds that could not be properly scored
// MARK: This function requires the user to use a scoring function
void
GA_Recomb::fitness_pruning( std::vector <DOCKMol> & pruned, std::vector <DOCKMol> & mfinal, Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & typer )
{
    Trace trace( "GA_Recomb::fitness_pruning()" );
    double start_time = time_seconds();

    // Counter for number removed in fitness pruning
    int poor_score = 0;
    int poor_ie = 0;
    int rmsd_overlap = 0;

    // Sort compounds by score
    for (int i=0; i<pruned.size(); i++){
         pruned[i].used = false;
         sort(pruned.begin(), pruned.end(), mol_sort);
    }

    // Step 1: Prune by energy score
    // NOTE: the energy function will be defined as the user specified in the input
    // MARK - Need to define a currention to the energy_cutoff with different combs of energies
    cout << "Initial Ensemble Size: " << pruned.size() << endl;
    for ( int i=0; i<pruned.size(); i++ ){
        if ( pruned[i].current_score > ga_energy_cutoff ){
             poor_score++;
             pruned.erase(pruned.begin()+i);
             i--;
        }
    }
    for ( int i=0; i<pruned.size(); i++ ){
        if (use_internal_energy){
            if (pruned[i].internal_energy > ie_cutoff ){
                poor_ie++;
                pruned.erase(pruned.begin()+i);
                i--;
            }
        }
    }
    //cout << "### Pruning by Score/Internal Energy ###" << endl; 
    cout << "Number pruned for poor score: " << poor_score << endl;
    cout << "Number pruned for poor internal energy: " << poor_ie << endl;
    cout <<"Size after score/internal energy pruning: " << pruned.size() << endl;

    //cout << "### Pruning by Similarity ###" << endl;
    //cout << "Ensemble size before similarity pruning: " << pruned.size() << endl;
    // Step 2: Prune DOCKMols using Hungarian Score
    Hungarian_RMSD h;
    pair <double, int> result;

    // These children need to be prepared with amber typer before rmsd calculation
    for ( int i=0; i<pruned.size(); i++ ){
         //cout << "Entering Amber typer" << endl;
         typer.skip_verbose_flag = true;
         typer.prepare_molecule(pruned[i], true, score.use_chem, score.use_ph4, score.use_volume);
         pruned[i].used = false;
    }

    // For each molecule in pruned
    for ( int i=0; i<pruned.size(); i++ ){
        // If it has not been marked used
        if ( !pruned[i].used ){
           // For each molecule listed below it
            for ( int j=i+1; j<pruned.size(); j++ ){
                // If the below molecule has not been marked used
                if ( !pruned[j].used ){
                   // Calculate and save hunarian rmsd
                   result = h.calc_Hungarian_RMSD_dissimilar( pruned[i], pruned [j] );

                   // If the two confomers are very similar, mark the one with higher energy as used
                   // If the two molecules do not have the same number of atoms than rmsd is -1000
                   //if ( rmsd < torsion_rmsd_cutoff & rmsd != -1000 ){
                   if (result.first < ga_heur_matched_rmsd && result.second < ga_heur_unmatched_num) { 
                       pruned[j].used = true;
                       rmsd_overlap++;
                   //cout << "Entering fitness pruning" << endl; 
                   //cout << "results.first: " << result.first << endl;
                   //cout << "results.second: " << result.second << endl;
                   }
                }
            }
        }
    }

    // At this point, only cluster heads will not be flagged 'used', add those to scored_generation
    for (int i=0; i<pruned.size(); i++){
        if (!pruned[i].used){
           mfinal.push_back(pruned[i]);
        }
    }
    
    // Print the number removed from fitness pruning
   // cout << endl << "#### Entering Similarity Pruning Routine ####" << endl;
    cout << "Number pruned for Hungarian similarity: " << (pruned.size() - mfinal.size()) << endl;
    cout << "Final Offspring Ensemble Size: " << mfinal.size() << endl;
    double stop_time = time_seconds();
    //cout << "\n" "-----------------------------------" "\n";
    if (verbose) cout << " Elapsed time fitness prune:\t" << stop_time - start_time
         << " seconds for GA\n\n";

    return;
}// end GA_Recomb::fitness_pruning()





// +++++++++++++++++++++++++++++++++++++++++
// Prune a vector of DOCKMols by ligand properties - used in breeding_exhaustive
void
GA_Recomb::hard_filter( vector <DOCKMol> & temp_vec )
{
    Trace trace( "GA_Recomb::hard_filter()" );
    // Make binary vector to keep status of molecules
    INTVec binary;
    for ( int i=0; i< temp_vec.size(); i++){
        // Zero means that the molecule is ok!
        binary.push_back(0);
    }

    //counter for number that fit the qualifications
    int invalid_mw = 0;
    int invalid_rot = 0;
    int invalid_HA = 0;
    int invalid_HD = 0;
    int invalid_formal = 0;

    if (temp_vec.size() == 0){
        cout << "No molecules made it to hard filter" << endl;
        return;
    }
    for (int i=0; i<temp_vec.size(); i++){
        // All below functions will assign the appropriate values to the DOCKMol object
        calc_descriptors(temp_vec[i]);

        // Prune if MW greater than constraint
        if (temp_vec[i].mol_wt > ga_constraint_mol_wt ){
           cout << "Did not pass MW check." << endl;
           invalid_mw++;
           binary[i] = 1;
        }

        // Prune if Rot bonds are greater
        if ( (temp_vec[i].rot_bonds > ga_constraint_rot_bon) ){
           cout << "Did not pass RB check." << endl;
           invalid_rot++;
           binary[i] = 1;
        }

        // Prune if hydrogen acceptors are greater than constraint
        if ((temp_vec[i].hb_acceptors > ga_constraint_H_accept) ){
           cout << "Did not pass HA check." << endl;
           invalid_HA++;
           binary[i] = 1;
        }
     
        // Prune if Hdonors are greater than constraint
        if ((temp_vec[i].hb_donors > ga_constraint_H_donor) ){
           cout << "Did not pass HD check." << endl;
           invalid_HD++;
           binary[i] = 1;
        }

        // Prune by formal charge range
        if (((temp_vec[i].formal_charge > ga_constraint_formal_charge) || (temp_vec[i].formal_charge < -ga_constraint_formal_charge)) ){
           cout << "Did not pass FC check." << endl;
           invalid_formal++;
           binary[i] = 1;
        }
        //cout << "HARD_FILTER: MOLECULE # " << i <<endl;
        //cout << "MW " << temp_vec[i].mol_wt << " Rot bonds " << temp_vec[i].rot_bonds << 
        //       " HA " << temp_vec[i].hb_acceptors << " HD " << temp_vec[i].hb_donors << 
        //       " FC " << temp_vec[i].formal_charge <<endl;
    } 

    // Copy temp_vec to new storage vector
    std::vector <DOCKMol> tmp;
    for (int i=0; i<temp_vec.size(); i++){
        tmp.push_back(temp_vec[i]);
    }
    temp_vec.clear();

    int erased = 0;
    for (int i=0; i<binary.size(); i++){
        if ( binary[i] == 0){
           temp_vec.push_back(tmp[i]);
        }
        else { erased++; }
    }
    // Print flags to show how many compounds would be removed if implementd
    cout << endl << invalid_mw << " compounds exceeded MW cutoff of " << ga_constraint_mol_wt<<endl;
    cout << invalid_rot << " compounds exceeded the rot bond cutoff of " << ga_constraint_rot_bon <<endl;
    cout << invalid_HA << " compounds exceeded the HA acceptor # of " << ga_constraint_H_accept <<endl;
    cout << invalid_HD << " compounds exceeded the HB donor # of " << ga_constraint_H_donor <<endl;
    cout << invalid_formal << " compounds exceeded the formal charge of " << ga_constraint_formal_charge <<endl;
    cout << erased << " compounds were removed " << endl;
    return;

} // end GA_Recomb::hard_filter()




// +++++++++++++++++++++++++++++++++++++++++
// Prune a single molecule by ligand properties (MARK-HA/HD is commented out)
bool
GA_Recomb::hard_filter_mol( DOCKMol & temp_vec )
{
    Trace trace( "GA_Recomb::hard_filter_mol()" );
    // Make binary vector to keep status of molecules
    int binary = 0;

    //counter for number that fit the qualifications
    int invalid_mw = 0;
    int invalid_rot = 0;
    int invalid_HA = 0;
    int invalid_HD = 0;
    int invalid_formal = 0;

    //All below functions will assign the appropriate values to the DOCKMol object
    calc_descriptors(temp_vec);
    
    // Prune by MW 
    if (temp_vec.mol_wt > ga_constraint_mol_wt ){
       invalid_mw++;
       binary = 1;
    }

    // Prune by rot bonds
    if ( (temp_vec.rot_bonds > ga_constraint_rot_bon) ){
       invalid_rot++;
       binary = 1;
    }

    // Prune by hydrogen donors
    /*if ( (temp_vec.hb_acceptors > ga_constraint_H_accept) ){
        invalid_HA++;
        binary = 1;
    }
     
    // Prune by hydrogen acceptors
    if ( (temp_vec.hb_donors > ga_constraint_H_donor) ){
       invalid_HD++;
       binary = 1;
    }*/

    // Prune by formal charge range   
    if (((temp_vec.formal_charge > ga_constraint_formal_charge) || (temp_vec.formal_charge < -ga_constraint_formal_charge)) ){
       invalid_formal++;
       binary = 1;
    }

    // DELETE REMOVE
    /*cout << "MW " << temp_vec.mol_wt << " Rot bonds " << temp_vec.rot_bonds << 
           " HA " << temp_vec.hb_acceptors << " HD " << temp_vec.hb_donors << 
           " FC " << temp_vec.formal_charge <<endl;*/

    // Print flags to show how many compounds would be removed if implementd
    if (verbose) cout << endl << "#### " << invalid_mw << " compounds exceeded MW cutoff of " << ga_constraint_mol_wt<<endl;
    if (verbose) cout << invalid_rot << " compounds exceeded the rot bond cutoff of " << ga_constraint_rot_bon <<endl;
    if (verbose) cout << invalid_HA << " compounds exceeded the HA acceptor # of " << ga_constraint_H_accept <<endl;
    if (verbose) cout << invalid_HD << " compounds exceeded the HB donor # of " << ga_constraint_H_donor <<endl;
    if (verbose) cout << invalid_formal << " compounds exceeded the formal charge of " << ga_constraint_formal_charge <<endl;
    //cout << erased << " compounds were removed " << endl;

    // False means that the molecule exceeds the cutoffs
    if ( binary == 1){ return false; }
    return true;
} // end GA_Recomb::hard_filter_mol()




// +++++++++++++++++++++++++++++++++++++++++
// Calculate the molecular weight of the dockmol object of a fragment
// (modified from amber_typer.cpp)
void
GA_Recomb::calc_descriptors( DOCKMol & mol )
{
   Trace trace( "GA_Recomb::calc_descriptors()" );
   calc_mol_wt (mol);   
   calc_rot_bonds (mol);
   num_HA_HD (mol);
   calc_formal_charge(mol);

  return;
} // end GA_Recomb::calc_descriptors()




// +++++++++++++++++++++++++++++++++++++++++
// Calculate the molecular weight of the dockmol object of a fragment
// (modified from amber_typer.cpp)
void
GA_Recomb::calc_mol_wt( DOCKMol & mol )
{
    Trace trace( "GA_Recomb::calc_mol_wt()" );
    //atomic weights from General Chemistry 3rd Edition Darrell D. Ebbing
    float mw = 0.0;

    for (int i=0; i<mol.num_atoms; i++) {

        string atom = mol.atom_types[i];

        if ( atom == "H")
            { mw += 1.00794; }

        else if ( atom == "C.3" || atom == "C.2" || atom == "C.1" || atom == "C.ar" || atom  == "C.cat" )
            { mw += 12.011; }

        else if ( atom == "N.4" || atom == "N.3" || atom == "N.2" || atom == "N.1" || atom  == "N.ar" ||
                  atom == "N.am" || atom == "N.pl3" )
            { mw += 14.00674; }

        else if ( atom == "O.3" || atom == "O.2" || atom == "O.co2" )
            { mw += 15.9994; }

        else if ( atom == "S.3" || atom == "S.2" || atom == "S.O" || atom == "S.o" || atom == "S.O2" ||
                  atom == "S.o2" )
            { mw += 32.066; }

        else if ( atom == "P.3" )
            { mw += 30.973762; }

        else if ( atom == "F" )
            { mw += 18.9984032; }

        else if ( atom == "Cl" )
            { mw += 35.4527; }

        else if ( atom == "Br" )
            { mw += 79.904; }

        else if ( atom == "I" )
            { mw += 126.90447; }

        else if ( atom == "Du" )
            { mw += 0; }

        else
            { cout <<"WARNING: Did not recognize the atom_type " <<atom <<" in DN_GA_Build::calc_mol_wt()\n"; }

    }

    mol.mol_wt = mw;
    return;

} // end GA_Recomb::calc_mol_wt()




// +++++++++++++++++++++++++++++++++++++++++
// Given a mol, populate the rot_bonds field
void
GA_Recomb::calc_rot_bonds( DOCKMol & mol )
{
    Trace trace( "GA_Recomb::calc_rot_bonds()" );
    // The number of rotatable bonds
    int counter = 0;

    for (int i=0; i<mol.num_bonds; i++){
        if (mol.amber_bt_id[i] != -1 ){
            counter++;
        }
    }

    // Assign it directly to the referenced mol object
    mol.rot_bonds = counter;

    return;

} // end GA_Recomb::calc_rot_bonds();




// +++++++++++++++++++++++++++++++++++++++++
// Given a mol, calculate the hydrogen acceptors and donor 
void
GA_Recomb::num_HA_HD( DOCKMol & mol )
{
    Trace trace( "GA_Recomb::num_HA_HD()" );
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

} //end GA_Recomb::num_HA_HD();





// +++++++++++++++++++++++++++++++++++++++++
// Given a mol, populate the formal_charge field
// NOTE: Molecules must be precharged with gastegier
void
GA_Recomb::calc_formal_charge( DOCKMol & mol )
{
    Trace trace( "GA_Recomb::calc_formal_charge()" );
    float charge = 0.0;

    // Iterate over all atoms, find the partial charge
    for (int i=0; i<mol.num_atoms; i++){
        charge += mol.charges[i];
    }

    // Assign it directly to the referenced mol object
    mol.formal_charge = charge;

    return;
} // end DN_GA_Build::calc_formal_charge();





            //////////////////////////////////////////
            //****SELECTION FUNCTIONS *****//
            /////////////////////////////////////////



// +++++++++++++++++++++++++++++++++++++++++
// Selection Methods - this function will call the appropriate functions
void
GA_Recomb::selection_method ( Master_Score & score, AMBER_TYPER & typer )
{
    Trace trace( "GA_Recomb::selection_method()" );
    // Make a tmp parent holder
    std::vector< DOCKMol > tmp_parents;

    // Update if parent status
    for (int i=0; i<parents.size(); i++){
        parents[i].parent = 1;
        tmp_parents.push_back(parents[i]);
    }
    parents.clear();

    for (int i=0; i<scored_generation.size(); i++){
        scored_generation[i].parent = 0;
    }

    // Print ensemble size of offspring
    cout <<"#### Final parent size before selection: " << tmp_parents.size() <<endl;
    cout <<"#### Final offspring size before selection: " << scored_generation.size() <<endl;
    // Call the appropriate function
    if (ga_selection_method_elitism ){
       cout << "Selection Method: Elitism" <<endl;
       selection_elite(tmp_parents, score, typer);
    }
    else if (ga_selection_method_tournament ){
       cout << "Selection Method: Tournament" <<endl;
        selection_tournament(tmp_parents, score, typer);
    }
    else if (ga_selection_method_roulette ){
       cout << "Selection Method: Roulette" <<endl;
        selection_roulette(tmp_parents, score, typer);
    }
    else if (ga_selection_method_sus ){
       cout << "Selection Method: SUS" <<endl;
        selection_sus(tmp_parents, score, typer);
    }
    else if (ga_selection_method_metropolis ){
       cout << "Selection Method: Metropolis" <<endl;
        selection_metropolis(tmp_parents, score, typer);
    }


    // Activate all molecules
    activate_vector(parents);
    return;
} //end GA_Recomb::selection_method();




// +++++++++++++++++++++++++++++++++++++++++
// The final ensemble will be limited by two methods - the ga_ensemble_size and the elitism method
void
GA_Recomb::selection_elite( std::vector< DOCKMol > & tmp_parents,  Master_Score & score, AMBER_TYPER & typer )
{
    Trace trace( "GA_Recomb::selection_elite()" );
    // If Extinction then update the size
    int final_ensemble_size = ga_ensemble_size;
    if (ga_extinction_on ){ final_ensemble_size = ga_extinction_keep; }


   // STEP 1: PERCENTAGE
   // If limiting the percent of parents for the next generation
   if ( ga_elitism_percent > 0 ){
       //cout <<" INSIDE ELITE PERCENT " <<endl;

       // Determine the number of parents that will be maintained
       float carry = 0.00;
       carry = float(tmp_parents.size()) *  ga_elitism_percent;
       // DELETE cout << "ga percent " <<   ga_elitism_percent << " carry " << carry << " parents " << tmp_parents.size() << " ensmble size " << ga_ensemble_size <<endl;
       
       // Only run if you are over the number of max parents for GA run
       //DELETE cout << " combined " << scored_generation.size() + int(carry) << endl;
       
       if ( (scored_generation.size() + tmp_parents.size()) > final_ensemble_size){
           
           // If using niching, rank molecules by descriptors
           if ( ga_niching ){
              // Crowding?
              if ( !ga_niche_sharing ){
                 // Compute rank and crowdingdistances
                 niche_crowding( tmp_parents, score );
                 niche_crowding( scored_generation, score );
                 // Check to make sure the correct scoring functions are being used
                 if ( score.c_desc.use_score == true || score.c_mg_nrg.use_score == true ){
                    sort( tmp_parents.begin(), tmp_parents.end(), compare_rank );
                    sort( scored_generation.begin(), scored_generation.end(), compare_rank );
                 }
              }
              else{
                 // Rescore parents to reset their score parameters for rescaling
                 score_parents(tmp_parents, score, typer);
                 niche_sharing( tmp_parents, score );
                 niche_sharing( scored_generation, score );
                 // Check to make sure the correct scoring functions are being used
                 if ( score.c_desc.use_score == true || score.c_mg_nrg.use_score == true ){
                    sort( tmp_parents.begin(), tmp_parents.end(), fitness_sort );
                    sort( scored_generation.begin(), scored_generation.end(), fitness_sort );
                 }
              }
           }

           if (! ga_extinction_on){
              // Add the appropriate number of parents
              for (int i=0; i<int(carry); i++){
                  parents.push_back(tmp_parents[i]);
              }
              tmp_parents.clear();

              // Add the number of offpsring to complete ensemble size
              int remainder_parents = final_ensemble_size - parents.size();

              // If you have more offspring than the desired size
              if ( remainder_parents < scored_generation.size()){
                  for ( int i=0; i<remainder_parents; i++ ){
                      parents.push_back( scored_generation[i] );
                  }
              }
              // Otherwise, add all the offspring ???
              else{
                   for ( int i=0; i<scored_generation.size(); i++){
                       parents.push_back( scored_generation[i] );
                   }
              }
           }
           else{
              // If extinction is on, keep number of parents if less than half of current new size
              if (int(carry) < ga_extinction_keep/2){
                 for (int i=0; i<int(carry); i++){
                    parents.push_back(tmp_parents[i]);
                 }
              }
	      else{
                  for (int i=0; i<ga_extinction_keep/2; i++){
                     parents.push_back(tmp_parents[i]);
                  }
              }
              // Add the number of offpsring to complete ensemble size
              int remainder_parents = ga_extinction_keep - parents.size();
              // If you have more offspring than the desired size
              if ( remainder_parents < scored_generation.size()){
                 for ( int i=0; i<remainder_parents; i++ ){
                     parents.push_back( scored_generation[i] );
                 }
              }
              tmp_parents.clear();
           }
       } 
       // Else, combine all parents with offspring 
       else{
           cout <<"No selection performed" <<endl;
           for ( int i=0; i<tmp_parents.size(); i++ ){
               parents.push_back( tmp_parents[i] );
           }
           tmp_parents.clear();
           for ( int i=0; i<scored_generation.size(); i++ ){
               parents.push_back( scored_generation[i] );
           }
           // Fitness score if being used
           if ( ga_niching ){
               if ( ga_niche_sharing ){
                  niche_sharing( parents, score );
               }
               else{
                  niche_crowding( parents, score );
               }
           }
       }
   }


    //------------------------------------------------------------------------
   // NUMBER
   // If using the top number of parents for next generation
   else if ( ga_elitism_number > 0 ){
       // DELETE cout <<" INSIDE ELITE NUM " <<endl;
       // If greater than ensemble size
       if ( (scored_generation.size() + tmp_parents.size()) > final_ensemble_size ) {
          // DELETE cout << " combined " << scored_generation.size() + tmp_parents.size() << endl;

          // If using niching, rank molecules by descriptors
          if ( ga_niching && !ga_elitism_combined ){
             if ( ga_niche_sharing == false ){
                // Compute rank and crowdingdistances
                niche_crowding( tmp_parents, score );
                niche_crowding( scored_generation, score );
                // Check to make sure the correct scoring functions are being used
                if ( score.c_desc.use_score == true || score.c_mg_nrg.use_score == true ){
                   sort( tmp_parents.begin(), tmp_parents.end(), compare_rank );
                   sort( scored_generation.begin(), scored_generation.end(), compare_rank );
                }
             }
             else{
                // Rescore parents to reset their score parameters for rescaling
                score_parents(tmp_parents, score, typer);
                niche_sharing( tmp_parents, score );
                niche_sharing( scored_generation, score );
                // Check to make sure the correct scoring functions are being used
                if ( score.c_desc.use_score == true || score.c_mg_nrg.use_score == true ){
                   sort( tmp_parents.begin(), tmp_parents.end(), fitness_sort );
                   sort( scored_generation.begin(), scored_generation.end(), fitness_sort );
                }
             }
          }

          // If using Max number of parents (ga_elitism_number is ga_ensmble size)
          // DELETE cout << "GA ELITE NUM " <<ga_elitism_number << " " << ga_ensemble_size << endl;
          if ( ga_elitism_number == ga_ensemble_size){
                
             // If combined
             if (ga_elitism_combined == true){
                // new DOCKMol vector
                vector <DOCKMol> tmp_gen;

                // Add parents (original and mutated)
                for ( int i=0; i<tmp_parents.size(); i++ ){
                    tmp_gen.push_back( tmp_parents[i] );
                }

                // Add offspring (from xover and mutaton)
                for ( int i=0; i<scored_generation.size(); i++ ){
                    tmp_gen.push_back( scored_generation[i] );
                }

                // If using niching, rank molecules by descriptors
                if ( ga_niching ){
                   if ( ga_niche_sharing == false ){
                      // Compute rank and crowdingdistances
                      niche_crowding( tmp_gen, score );
                      // Check to make sure the correct scoring functions are being used
                      if ( score.c_desc.use_score == true || score.c_mg_nrg.use_score == true ){
                         sort( tmp_gen.begin(), tmp_gen.end(), compare_rank );
                      }
                   }
                   else{
                      // Rescore parents to reset their score parameters for rescaling
                      score_parents(tmp_parents, score, typer);
                      // Compute rank and crowdingdistances
                      niche_sharing( tmp_gen, score );
                      // Check to make sure the correct scoring functions are being used
                      if ( score.c_desc.use_score == true || score.c_mg_nrg.use_score == true ){
                         sort( tmp_gen.begin(), tmp_gen.end(), fitness_sort );
                      }
                   }
                }
                else{
                   sort( tmp_gen.begin(), tmp_gen.end(), mol_sort );
                }

                if (! ga_extinction_on){
                   // Add molecules to parents vector
                   for ( int i=0; i<ga_ensemble_size; i++ ){
                       parents.push_back( tmp_gen[i] );
                   }
                }
                else{
                   // If using extinction
                   for ( int i=0; i<ga_extinction_keep; i++ ){
                       parents.push_back( tmp_gen[i] );
                   }
                }
                // Clear tmp_gen
                tmp_gen.clear();
             }
             else{
                //DELETE cout <<"inside elitism not combined" <<endl;
 
                if (! ga_extinction_on){
                   // Do not know if we fit the criteria to half each pop
                   int p_diff = tmp_parents.size() - ga_elitism_number/2;
                   int c_diff = scored_generation.size() - ga_elitism_number/2;
                   //DELETE cout << " p diff " << p_diff << " c diff " << c_diff << endl;
                   //OPTION 1:
                   if ( (p_diff >= 0) && (c_diff < 0) ){ 
                      // Add all of offspring 
                      for ( int i=0; i<scored_generation.size(); i++ ){
                          parents.push_back( scored_generation[i] );
                      }
                      // Add the appropriate number of parents 
                      int remainder_parents = ga_ensemble_size - parents.size();
                      //DELETE cout << "remainder " << remainder_parents << endl; 
                      for ( int i=0; i<remainder_parents; i++){
                          parents.push_back(tmp_parents[i]);
                      }
                   }
 
                   //OPTION 2:
                   else if ( (c_diff >= 0) && (p_diff < 0) ){ 
                      // Add all parents 
                      for ( int i=0; i<tmp_parents.size(); i++ ){
                          parents.push_back( tmp_parents[i] );
                      }
                      // Add the appropriate number of offspring
                      int remainder_parents = ga_ensemble_size - parents.size();
                      //DELETE cout << "remainder " << remainder_parents << endl; 
                      for ( int i=0; i<remainder_parents; i++){
                          parents.push_back(scored_generation[i]);
                      }
                   }

                   //OPTION 3:
                   else { 
                      //DELETE cout <<"Equal numbers of parents + offspring" <<endl;
                      // Add equal number of offspring and parents
                      for ( int i=0; i<ga_ensemble_size/2; i++ ){
                          parents.push_back( scored_generation[i] );
                      }
                      for ( int i=0; i<ga_ensemble_size/2; i++){
                          parents.push_back(tmp_parents[i]);
                      }
                   }
                }
                else{
                   // If using extinction
                   if (tmp_parents.size() > ga_extinction_keep/2){
                      for (int i=0; i<ga_extinction_keep/2; i++){
                         parents.push_back(tmp_parents[i]);
                      }
                   }
                   else{
                      for (int i=0; i<tmp_parents.size(); i++){
                         parents.push_back(tmp_parents[i]);
                      }
                   }
                   tmp_parents.clear();

                   // Equivalent for Offspring
                   if (scored_generation.size() > ga_extinction_keep/2){
                      for (int i=0; i<ga_extinction_keep/2; i++){
                         parents.push_back(scored_generation[i]);
                      }
                   }
                   else{
                      for (int i=0; i<scored_generation.size(); i++){
                         parents.push_back(scored_generation[i]);
                      }
                   }
                }
             }//else
          }// end ELITISM MAX

          // If the user set a number of parents that want
          else{
             if (! ga_extinction_on){
                //DELETE cout <<"if user specificed number "<<endl;
                if ( ga_elitism_number < tmp_parents.size()){
                   for ( int i=0; i<ga_elitism_number; i++){
                       parents.push_back( tmp_parents[i] );
                   }
                   //DELETE cout << " ga elitism number " << ga_elitism_number << " parents size " << parents.size() << endl;
                   // Add the number of offpsring to complete ensemble size
                   int remainder_parents = ga_ensemble_size - ga_elitism_number;
                   //DELETE cout << "num parents " << remainder_parents <<endl;
                   if (remainder_parents < scored_generation.size()){
                      for ( int i=0; i<remainder_parents; i++){
                          parents.push_back(scored_generation[i]);
                      } 
                   }
                   else{
                      for ( int i=0; i<scored_generation.size(); i++){
                          parents.push_back(scored_generation[i]);
                      } 
                   }
                }
                else{
                   for ( int i=0; i<tmp_parents.size(); i++){
                       parents.push_back( tmp_parents[i] );
                   }
                   // Add the number of offpsring to complete ensemble size
                   int remainder_parents = ga_ensemble_size - parents.size();
                   //DELETE cout << "num parents " << remainder_parents <<endl;
                   if (remainder_parents < scored_generation.size()){
                      for ( int i=0; i<remainder_parents; i++){
                          parents.push_back(scored_generation[i]);
                      }
                   }
                   else{
                      for ( int i=0; i<scored_generation.size(); i++){
                          parents.push_back(scored_generation[i]);
                      } 
                   }
                }
             }
             // If using extinction
             else{
                int carry = 0;
                if (tmp_parents.size() > ga_extinction_keep/2){
                   if (ga_elitism_number < ga_extinction_keep/2){
                      carry = ga_elitism_number;
                   }
                   else{
                      carry = ga_extinction_keep/2;
                   }
                   for ( int i=0; i<carry; i++){
                       parents.push_back( tmp_parents[i] );
                   }
                }
                else{
                   for ( int i=0; i<tmp_parents.size(); i++){
                       parents.push_back( tmp_parents[i] );
                   }
                }
                tmp_parents.clear();
       
                // For the offspring
                // Add the number of offpsring to complete ensemble size
                int remainder_parents = ga_extinction_keep - parents.size();
                //DELETE cout << "num parents " << remainder_parents <<endl;
                if (remainder_parents < scored_generation.size()){
                   for ( int i=0; i<remainder_parents; i++){
                       parents.push_back(scored_generation[i]);
                   }
                }
             }
          }// end ELITISM NUM
       }// end IF greater than CUTOFF     
           
       // If not over the cutoff  
       else{
          cout <<"No selection performed - not over initial cutoff" <<endl;
          for ( int i=0; i<tmp_parents.size(); i++ ){
              parents.push_back( tmp_parents[i] );
          }
          tmp_parents.clear();
          for ( int i=0; i<scored_generation.size(); i++ ){
              parents.push_back( scored_generation[i] );
          }
          // Fitness score if being used
          if ( ga_niching ){
              if ( ga_niche_sharing ){
                 niche_sharing( parents, score );
              }
              else{
                 niche_crowding( parents, score );
              }
          }
       }
   }
   // Sort parents
   sort(parents.begin(), parents.end(), mol_sort);

   // Clear ancestors and offspring
   tmp_parents.clear();
   //scored_generation.clear();


  return;

} //end GA_Recomb::selection_elite();




// +++++++++++++++++++++++++++++++++++++++++
// Two molecules will compete to move to the next generation - survivor has the best score
void
GA_Recomb::selection_tournament( std::vector< DOCKMol > & tmp_parents, Master_Score & score, AMBER_TYPER & typer )
{
   Trace trace( "GA_Recomb::selection_tournament()" );
   // There are two competition methods - parents vs children or parent vs parent. If parent vs children
   // is selected, you will not guarantee a certain number of any group
  
   // Determine which is the limiting vector
   int p_diff = 0;
   int c_diff = 0;
   if (verbose) cout << "TOURNAMENT ga_extinction_on: "<< ga_extinction_on <<endl;
   p_diff = tmp_parents.size() - ga_ensemble_size;
   c_diff = scored_generation.size() - ga_ensemble_size;

   // If Extinction then update the size
   int final_ensemble_size = ga_ensemble_size;
   if (ga_extinction_on ){ final_ensemble_size = ga_extinction_keep; }

   // If using niching
   if ( ga_niching ){
      if ( ga_niche_sharing ){
         // Rescore parents to reset their score parameters for rescaling
         score_parents(tmp_parents, score, typer);
         // Compute rank and crowding distances
         if (verbose) cout << "Tournament niching: extinction:" << ga_extinction_on << ", keep: " << final_ensemble_size <<endl;
         niche_sharing( tmp_parents, score );
         niche_sharing( scored_generation, score );

         // Determine which score is being used
         if ( score.c_desc.use_score == true || score.c_mg_nrg.use_score == true ){
            sort( tmp_parents.begin(), tmp_parents.end(), fitness_sort );
            sort( scored_generation.begin(), scored_generation.end(), fitness_sort );
         }
      }
      else{
         // Compute rank and crowding distances
         niche_crowding( tmp_parents, score );
         niche_crowding( scored_generation, score );

         // Determine which score is being used
         if ( score.c_desc.use_score == true || score.c_mg_nrg.use_score == true ){
            sort( tmp_parents.begin(), tmp_parents.end(), compare_rank );
            sort( scored_generation.begin(), scored_generation.end(), compare_rank );
         }
      }
   }

   //------------------------------------------------------------------------
   // If Parents vs Children
   // The number of tournaments is limited to the smallest ensemble
   if (ga_tournament_p_vs_c){
       if ( (tmp_parents.size() + scored_generation.size()) > final_ensemble_size ){

          // OPTION 1:
          // If one is less than max number of picks
          if ( (p_diff < 0) || (c_diff < 0) ){

             // OPTION 1A:
             //Determine which is smallest ( most negative)
             if ( p_diff < c_diff ){
                if (verbose) cout << "1A: " << p_diff << " " << c_diff << endl;
      if (verbose) cout << "Tournament: ga_ensemble: " << ga_ensemble_size << 
      "tmp_parents " << tmp_parents.size() << " scored " << scored_generation.size() <<endl;
                while ( (parents.size() < tmp_parents.size())){
                   //cout << "size" << parents.size() <<endl;
                   // Define index for random picks
                   int index_p = (rand() % tmp_parents.size()); 
                   int index_c = (rand() % scored_generation.size()); 
          
                   // If using niching
                   if ( ga_niching ){
                      if ( ga_niche_sharing ){
                         tournament_trad( tmp_parents, index_p, scored_generation, index_c, parents );
                      }
                      else{ tournament_niche( tmp_parents, index_p, scored_generation, index_c, parents ); }
                   }
            
                   // If not using niching, Pick the molecule with the best score
                   else {
                      tournament_trad( tmp_parents, index_p, scored_generation, index_c, parents );
                   }
                }
             }

             // OPTION 1B:
             //Determine which is smallest
             else{
                //cout << "other " << p_diff << " " << c_diff << endl;
                if (verbose) cout << "1B: " << p_diff << " " << c_diff << endl;
      if (verbose) cout << "Tournament: ga_ensemble: " << ga_ensemble_size << 
      "tmp_parents " << tmp_parents.size() << " scored " << scored_generation.size() <<endl;
                while ( (parents.size() < scored_generation.size())){
                   //cout << "size" << parents.size() <<endl;
                   // Define index for random picks
                   int index_p = (rand() % tmp_parents.size()); 
                   int index_c = (rand() % scored_generation.size()); 
          
                   // If using niching
                   if ( ga_niching ){
                      if ( ga_niche_sharing ){
                         tournament_trad( tmp_parents, index_p, scored_generation, index_c, parents );
                      }
                      else{ tournament_niche( tmp_parents, index_p, scored_generation, index_c, parents ); }
                   }
            
                   // If not completing niching
                   else {
                      tournament_trad( tmp_parents, index_p, scored_generation, index_c, parents );
                   }
                }
             }
          }// end OPTION 1

          // OPTION 2:
          // They are both over the cut off
          else{
              while (parents.size() < final_ensemble_size ){
                 // Define index for random picks
                 int index_p = (rand() % tmp_parents.size()); 
                 int index_c = (rand() % scored_generation.size()); 
          
                 // If using niching
                 if ( ga_niching ){
                    if ( ga_niche_sharing ){
                       tournament_trad( tmp_parents, index_p, scored_generation, index_c, parents );
                    }
                    else{ tournament_niche( tmp_parents, index_p, scored_generation, index_c, parents ); }
                 }
            
                 // If not using niching
                 else{
                    tournament_trad( tmp_parents, index_p, scored_generation, index_c, parents );
                 }
             }
          }//end OPTION 2
      }
      else{
         // Combine parents and offspring
         for (int i=0; i<tmp_parents.size(); i++){
             parents.push_back(tmp_parents[i]);
         }
         for (int i=0; i<scored_generation.size(); i++){
             parents.push_back(scored_generation[i]);
         }
         // Fitness score if being used
         if ( ga_niching ){
            if ( ga_niche_sharing ){
               niche_sharing( parents, score );
            }
            else{
               niche_crowding( parents, score );
            }
         }
      }    
   }

   //------------------------------------------------------------------------
   // If P vs P and C vs C
   else{
      if ( (tmp_parents.size() + scored_generation.size()) > final_ensemble_size ){
         // For the parents
         if ( tmp_parents.size() > final_ensemble_size/2){
            while ( parents.size() < final_ensemble_size/2 ){
               // Define index for random picks
               int index_p1 = (rand() % tmp_parents.size()); 
               int index_p2 = (rand() % tmp_parents.size()); 

               if (index_p1 != index_p2){
                  // If using niching
                  if ( ga_niching ){
                     if ( ga_niche_sharing ){
                        tournament_trad( tmp_parents, index_p1, tmp_parents, index_p2, parents );
                     }
                     else{ tournament_niche( tmp_parents, index_p1, tmp_parents, index_p2, parents ); }
                  }
            
                  // If not completing niching
                  else{
                     tournament_trad( tmp_parents, index_p1, tmp_parents, index_p2, parents );
                  }
               }
               if ( tmp_parents.size() < 2 ){ break; }
            }
         }
         // else add all parents to new ensemble
         else{
            for ( int i=0; i<tmp_parents.size(); i++){
               parents.push_back(tmp_parents[i]);
            }
         }

         // For the children
         if ( scored_generation.size() > final_ensemble_size/2 ){
            while ( parents.size() < final_ensemble_size ){
               // Define idex for random picks
               int index_c1 = (rand() % scored_generation.size()); 
               int index_c2 = (rand() % scored_generation.size()); 

               if (index_c1 != index_c2){
                  // If using niching
                  if ( ga_niching ){
                     if ( ga_niche_sharing ){
                        tournament_trad( scored_generation, index_c1, scored_generation, index_c2, parents );
                      }
                      else{ tournament_niche( scored_generation, index_c1, scored_generation, index_c2, parents ); }
                  }
            
                  // If not completing niching
                  else{
                     tournament_trad( scored_generation, index_c1, scored_generation, index_c2, parents );
                  }
               }
               if ( scored_generation.size() < 2 ){ break; }
            } 
         }
         // else add all offspring to new ensemble
         else{
            for ( int i=0; i<scored_generation.size(); i++){
               parents.push_back(scored_generation[i]);
            }
         }
      }

      else{
         // Combine parents and offspring
         for (int i=0; i<tmp_parents.size(); i++){
             parents.push_back(tmp_parents[i]);
         }
         for (int i=0; i<scored_generation.size(); i++){
             parents.push_back(scored_generation[i]);
         }
         // Fitness score if being used
         if ( ga_niching ){
            if ( ga_niche_sharing ){
               niche_sharing( parents, score );
            }
            else{
               niche_crowding( parents, score );
            }
         }
      }    
   }

   // Sort parents
   sort(parents.begin(), parents.end(), mol_sort);

   // Clear necessary vector
   tmp_parents.clear();
   //scored_generation.clear();


   return;
} //end GA_Recomb::selection_tournament()




// +++++++++++++++++++++++++++++++++++++++++
// Tournament function
void
GA_Recomb::tournament_trad( std::vector <DOCKMol> & mol1, int index1, std::vector <DOCKMol> & mol2, int index2, std::vector <DOCKMol> & save )
{
   Trace trace( "GA_Recomb::tournament_trad()" );
   // If mol1 score is better than mol2,  
   if (mol1[index1].current_score < mol2[index2].current_score ){
      save.push_back(mol1[index1]);
      mol1.erase(mol1.begin()+index1);
   }
   else{
      save.push_back(mol2[index2]);
      mol2.erase(mol2.begin()+index2);
   }

   return;
} // end tournament_trad()



// +++++++++++++++++++++++++++++++++++++++++
// Tournmanet niching
void
GA_Recomb::tournament_niche( std::vector <DOCKMol> & mol1, int index1, std::vector <DOCKMol> & mol2, int index2, std::vector <DOCKMol> & save )
{
   Trace trace( "GA_Recomb::tournament_niche()" );
   // Pick the molecule with the best rank
   if ( mol1[index1].rank < mol2[index2].rank ){
      save.push_back(mol1[index1]);
      mol1.erase(mol1.begin()+index1);
   }
   else if ( mol1[index1].rank > mol2[index2].rank ){
      save.push_back(mol2[index2]);
      mol2.erase(mol2.begin()+index2);
   }
   else if ( mol1[index1].rank == mol2[index2].rank ){
      // Check the crowding distance
      if ( mol1[index1].crowding_dist > mol2[index2].crowding_dist ){
         save.push_back(mol1[index1]);
         mol1.erase(mol1.begin()+index1);
      }
      else{
         save.push_back(mol2[index2]);
         mol2.erase(mol2.begin()+index2);
      }
  }
            
   return;
} // emd GA_Recomb::tournament_niche()




// +++++++++++++++++++++++++++++++++++++++++
// Two molecules will compete to move to the next generation - survivor has the best score
void
GA_Recomb::selection_roulette( std::vector< DOCKMol > & tmp_parents, Master_Score & score, AMBER_TYPER & typer )
{
   Trace trace( "GA_Recomb::selection_roulette()" );

   // If Extinction then update the size
   int final_ensemble_size = ga_ensemble_size;
   if ( ga_extinction_on ){ final_ensemble_size = ga_extinction_keep; }

   if (ga_roulette_separate){
      cout << "Running Roulette Separate" <<endl;
      // If using niching
      if ( ga_niching ){
         // Rescore parents to reset their score parameters for rescaling
         score_parents(tmp_parents, score, typer);
         if (verbose) cout << "BEFORE NICHING" <<endl;
         // Compute rank and crowding distances
         niche_sharing( tmp_parents, score );
         niche_sharing( scored_generation, score );
         if (verbose) cout << "AFTER NICHING" <<endl;

         // Determine which score is being used
         if ( score.c_desc.use_score == true || score.c_mg_nrg.use_score == true ){
            sort( tmp_parents.begin(), tmp_parents.end(), fitness_sort );
            sort( scored_generation.begin(), scored_generation.end(), fitness_sort );
         }
      }

      if (tmp_parents.size() > final_ensemble_size/2){
         int cycle = 0;
         int size = 0;
         while ( parents.size() < final_ensemble_size/2 ){ 
            roulette(tmp_parents,final_ensemble_size/2, false);
         
            if ( parents.size() == size ){cycle++;}
            if ( cycle > ga_max_mut_cycles ){break;}
	    if ( tmp_parents.size() < 2 ) {break;}
            size = parents.size();
         }
      }
      else{
         for (int i=0; i<tmp_parents.size(); i++){
             parents.push_back(tmp_parents[i]);
           }
      }

      if (scored_generation.size() > final_ensemble_size/2){
         int cycle = 0;
         int size = 0;
         while ( parents.size() < final_ensemble_size ){ 
            roulette(scored_generation,final_ensemble_size, false);
            // If parents size does not increase X number of times, quit
            if ( parents.size() == size ){cycle++;}
            if ( cycle > ga_max_mut_cycles ){break;}
	    if ( scored_generation.size() < 2 ) {break;}
            size = parents.size();
         }
      }
      else{
          for (int i=0; i<scored_generation.size(); i++){
              parents.push_back(scored_generation[i]);
	  }
      }
      if (verbose) cout << "FINISHED ROULETTE" <<endl;
   }
   else{
      if ( tmp_parents.size() + scored_generation.size() > final_ensemble_size ){
	 cout << "Running Roulette Combined" <<endl;

         //Combine parents and offspring
         std::vector <DOCKMol> temp;
         for (int i=0; i<tmp_parents.size(); i++){
            temp.push_back(tmp_parents[i]);
         }
         for (int i=0; i<scored_generation.size(); i++){
            temp.push_back(scored_generation[i]);
         }
         // Rescore to reset their score parameters for rescaling
         score_parents(temp, score, typer);

         // If using niching
         if ( ga_niching ){
            // Compute fitness and crowding distances
            niche_sharing( temp, score );

            // Determine which score is being used
            if ( score.c_desc.use_score == true || score.c_mg_nrg.use_score == true ){
               sort( temp.begin(), temp.end(), fitness_sort );
            }
         }
         else{
            // Sort parents
            sort(temp.begin(), temp.end(), mol_sort);
         }

         int cycle = 0;
         int size = 0;
         while ( parents.size() < final_ensemble_size ){ 
            roulette(temp, final_ensemble_size, false);
            // If parents size does not increase X number of times, quit
            if ( parents.size() == size ){cycle++;}
            if ( cycle > ga_max_mut_cycles ){break;}
	    if ( temp.size() < 2 ) {break;}
            size = parents.size();
         }

         //Clear
         temp.clear();
      }
      else{
        for (int i=0; i<tmp_parents.size(); i++){
           parents.push_back(tmp_parents[i]);
        }
         
         for (int i=0; i<scored_generation.size(); i++){
            parents.push_back(scored_generation[i]);
         }
         // Fitness score if being used
         if ( ga_niching ){
            if ( ga_niche_sharing ){
               niche_sharing( parents, score );
            }
            else{
               niche_crowding( parents, score );
            }
         }
      }
   }

   // Sort parents
   sort(parents.begin(), parents.end(), mol_sort);

   // Clear necessary vector
   tmp_parents.clear();
   //scored_generation.clear();
  
   return;
}//end GA_Recomb::selection_roulette()




// +++++++++++++++++++++++++++++++++++++++++
// Two molecules will compete to move to the next generation - survivor has the best score
void
GA_Recomb::selection_sus( std::vector< DOCKMol > & tmp_parents, Master_Score & score, AMBER_TYPER & typer )
{
   Trace trace( "GA_Recomb::selection_sus()" );
   if (ga_sus_separate){
      // If using niching
      if ( ga_niching ){
         // Rescore parents to reset their score parameters for rescaling
         //score_parents(tmp_parents, score, typer);
         // Compute rank and crowding distances
         niche_sharing( tmp_parents, score );
         niche_sharing( scored_generation, score );

         // Determine which score is being used
         if ( score.c_desc.use_score == true || score.c_mg_nrg.use_score == true ){
            sort( tmp_parents.begin(), tmp_parents.end(), fitness_sort );
            sort( scored_generation.begin(), scored_generation.end(), fitness_sort );
         }
      }

      if (tmp_parents.size() >= ga_ensemble_size/2){
         roulette(tmp_parents,ga_ensemble_size/2, true);
      }
      else{
        for (int i=0; i<tmp_parents.size(); i++){
           parents.push_back(tmp_parents[i]);
        }
      }

      if (scored_generation.size() >= ga_ensemble_size/2){
         roulette(scored_generation,ga_ensemble_size/2, true);
      }
      else{
         for (int i=0; i<scored_generation.size(); i++){
            parents.push_back(scored_generation[i]);
         }
      }
   }
   else{
      if ( tmp_parents.size() + scored_generation.size() > ga_ensemble_size ){
         /*if ( ga_niching ){
            // Rescore parents to reset their score parameters for rescaling
            score_parents(tmp_parents, score, typer);
         }*/
         //Combine parents and offspring
         std::vector <DOCKMol> temp;
         for (int i=0; i<tmp_parents.size(); i++){
            temp.push_back(tmp_parents[i]);
         }
         for (int i=0; i<scored_generation.size(); i++){
            temp.push_back(scored_generation[i]);
         }

         // If using niching
         if ( ga_niching ){
            // Compute rank and crowding distances
            niche_sharing( temp, score );

            // Determine which score is being used
            if ( score.c_desc.use_score == true || score.c_mg_nrg.use_score == true ){
               sort( temp.begin(), temp.end(), fitness_sort );
            }
         }
         else{
            // Sort parents
            sort(temp.begin(), temp.end(), mol_sort);
         }

         roulette(temp, ga_ensemble_size, true);

         //Clear
         temp.clear();
      }
      else{
        for (int i=0; i<tmp_parents.size(); i++){
           parents.push_back(tmp_parents[i]);
        }
         
         for (int i=0; i<scored_generation.size(); i++){
            parents.push_back(scored_generation[i]);
         }
         // Fitness score if being used
         if ( ga_niching ){
            if ( ga_niche_sharing ){
               niche_sharing( parents, score );
            }
            else{
               niche_crowding( parents, score );
            }
         }
      }
   }

   // Sort parents
   sort(parents.begin(), parents.end(), mol_sort);

   // Clear necessary vector
   tmp_parents.clear();
   //scored_generation.clear();
  
   return;
}//end GA_Recomb::selection_roulette()





// +++++++++++++++++++++++++++++++++++++++++
// The roulette function that will select molecules based on necessary input
void
GA_Recomb::roulette( std::vector< DOCKMol > & tmp, int num_mols, bool sus )
{
   Trace trace( "GA_Recomb::roulette()" );
   // Initialize vector to hold probabilitis
   std::vector <float> probs;
   std::vector <float> roulette;
   INTVec used_mol;
   float fitness_sum = 0.000;

   if ( ga_niching ){
      for (int i=0; i<tmp.size(); i++){
          fitness_sum += tmp[i].fitness;
      }
      for (int i=0; i<tmp.size(); i++){
          probs.push_back(tmp[i].fitness/fitness_sum);
          //DELETE cout << " individial prob " << tmp[i].fitness/fitness_sum << endl;
      }
   } 
   else{
      for (int i=0; i<tmp.size(); i++){
          fitness_sum += tmp[i].current_score;
      }
      //DELETE cout << "fitness sum " << fitness_sum << endl;

      for (int i=0; i<tmp.size(); i++){
          probs.push_back(tmp[i].current_score/fitness_sum);
          //cout << " individial prob " << tmp[i].current_score/fitness_sum << endl;
      }
   }
   // Create Roulette wheel for parents
   for (int i=0; i<probs.size(); i++){
       // Add the probability to the previous probability
       if (i == 0){
          roulette.push_back(probs[0]);
          //cout << " Roulette wheel: " << roulette[i] << endl;
       }
       else {
          roulette.push_back(probs[i] + roulette[i-1]);
          //cout << " Roulette wheel: " << roulette[i] << endl;
       }
   }

   //------------------------------------------------------------------------
   // Pick appropriate number of mols based on ga_ensemble_size
   if (sus == false){
      // Pick a random float between 0 and 1
      //float rand_prob = float(rand()) % float(roulette[roulette.size()]); 
      float rand_prob = ((float) rand() / RAND_MAX);
      //cout << "rand prob " << rand_prob << endl;

      bool less_than = false;
      int picks = 0;
 
      while (less_than == false){
           for ( int i=0; i<roulette.size(); i++){
               // If the rand number is less than current
               if ((rand_prob < roulette[i])){
                  less_than = true;
                     parents.push_back(tmp[i]);
                     tmp.erase(tmp.begin()+i);
                     //cout << " i roulette " << i << endl;
                     picks++;
               }
               if (less_than){break;}
            }   
            //if (extinction){
            //   if (picks == max_extinction_num){break;}//EXTINCTION Not Implemented
            //}
      }
   }
   //------------------------------------------------------------------------
   // THIS HAS TO BE REDONE
   else{
     // Complete SUS again
     bool end = false;
     bool start_over = false;

     // Keep track of cycles
     int cycle = 0;
     // While not at max num_mols and not at end of list
     while ( (parents.size() < num_mols) || (cycle < ga_max_mut_cycles) ){
        // DELETE cout <<"here f parents size " << parents.size() <<endl;
        // Pick a random float between 0 and 1
        float rand_prob = ((float) rand() / RAND_MAX);
        //DELETE cout << "rand prob " << rand_prob <<endl;
        // Define the distance used to select molecules
        float dist = roulette[roulette.size()-1]/(float(num_mols));
        //DELETE cout << "distance " << dist <<endl;

        // Pick molecules from list starting at rand_prob
        for ( int i=0; i<roulette.size(); i++){
            // If the random number is less than current location
            if ( rand_prob < roulette[i] ){
               // Update used and add to parent
               parents.push_back(tmp[i]);
               used_mol[i] = 1;

               //DELETE cout << "found starting molecule" <<endl;
               // Update rand number
               rand_prob += dist;           
               // While not at max num_mols and not at end of list
               while ( (parents.size() < num_mols) ){
                  // DELETE cout << "while less than num " << parents.size() <<endl;
                  // For each molecule past the starting location
                  for ( int j=i+1; j<roulette.size(); j++){
                      if ( rand_prob < roulette[j] ) {
                         // If not used
                         if ( used_mol[j] == 0 ){
                            // Save molecule
                            parents.push_back(tmp[j]);
                            used_mol[j] = 1;
                            // Update rand number
                            rand_prob += dist;           
                            //DELETE cout << "new prob " << rand_prob <<endl;
                            // If max number, break
                            if ( parents.size() == num_mols ){end = true;}
                         }
                      }
                      //DELETE cout << "j " << j << " rand prob " << rand_prob <<endl;
                      // If past last prob or the lenght of the wheel and not fulfilled all molecules
                      // redundant
                      if ( ((j == roulette.size() -1) || ( rand_prob >= roulette[roulette.size()-1] )) 
                         && ( end == false ) ){
                         // DELETE cout << " over list j " << j << " rand prob " << rand_prob << " final " << roulette[roulette.size()-1]<<endl;
                         // Start at the beginning
                         rand_prob = rand_prob - roulette[roulette.size()-1];
                         // DELETE cout << "new prob " << rand_prob <<endl;
                         for ( int k=0; k<roulette.size(); k++){
                             if ( rand_prob < roulette[j] ) {
                                // If not used
                                if ( used_mol[j] == 0 ){
                                   // Save molecule
                                   parents.push_back(tmp[j]);
                                   used_mol[j] = 1;
                                   // Update rand number
                                   rand_prob += dist;           
                                   //DELETE cout << "new prob " << rand_prob <<endl;
                                   // If max number, break
                                   if ( parents.size() == num_mols ){end = true;}
                                }
                                // If used, then picking the same molecules
                                else{
                                  // Change random number and start over
                                  start_over = true;  
                                  end = true;
                                  // DELETE cout << "had to start from the beginning" <<endl;
                               }
                             }// end if
                             if (end || start_over){break;}
                         }// end second for
                      }// end if need to repeat  
                      if (end || start_over){break;}
                  }// end for
                  if (end || start_over){break;}
               }// end while
            }// end first if     
            // If i = max size of roulette then exit
            if (i == roulette.size()-1){end = true;}
            if (end || start_over){break;}
        }// end first for
        if ( cycle < ga_max_mut_cycles){
           cycle += 1;
        } else { break; }
        if (end || start_over){break;}
     }// end main while
   }// if using sus
   //Clear vectors
   probs.clear();
   roulette.clear();
   used_mol.clear();

   return;
}//end GA_Recomb::roulette()








// +++++++++++++++++++++++++++++++++++++++++
// Activate all of the atoms AND bonds in a given DOCKMol - called in sample_minimized_torsions 
void
GA_Recomb::selection_metropolis( std::vector <DOCKMol> & tmp_parents, Master_Score & score, AMBER_TYPER & typer )
{
   Trace trace( "GA_Recomb::selection_metropolis()" );
   // STEP 1: Parents and Offspring selected separately
   if (ga_metropolis_separate){
      // If using niching
      if ( ga_niching ){
         // Rescore parents to reset their score parameters for rescaling
         score_parents(tmp_parents, score, typer);
         // Compute rank and crowding distances
         niche_sharing( tmp_parents, score );
         niche_sharing( scored_generation, score );

         // Determine which score is being used
         if ( score.c_desc.use_score == true || score.c_mg_nrg.use_score == true ){
            sort( tmp_parents.begin(), tmp_parents.end(), fitness_sort );
            sort( scored_generation.begin(), scored_generation.end(), fitness_sort );
         }
      }

      if (tmp_parents.size() >= ga_ensemble_size/2){
         //DELETE cout << "metropolis for parents " << tmp_parents.size() <<endl;
         //mc_metropolis(tmp_parents,ga_ensemble_size/2);
         true_mc_metropolis(tmp_parents,ga_ensemble_size/2);
      }
      else{
        for (int i=0; i<tmp_parents.size(); i++){
           parents.push_back(tmp_parents[i]);
        }
        //DELETE cout << "no metropolis for parents " << parents.size() <<endl;
      }

      if (scored_generation.size() >= ga_ensemble_size/2){
         //DLETE cout << "metropolis for offspring " << scored_generation.size() <<endl;
         //mc_metropolis(scored_generation,ga_ensemble_size/2);
         true_mc_metropolis(scored_generation,ga_ensemble_size/2);
      }
      else{
         for (int i=0; i<scored_generation.size(); i++){
            parents.push_back(scored_generation[i]);
         }
         //DELETE cout <<"no metropolis for offspring " << parents.size() <<endl;
      }
   }

   //------------------------------------------------------------------------
   // STEP 2: Combined
   else{
      if (tmp_parents.size() + scored_generation.size() > ga_ensemble_size){
         // If using niching
         if ( ga_niching ){
            // Rescore parents to reset their score parameters for rescaling
            score_parents(tmp_parents, score, typer);
         }

         //Combine parents and offspring
         std::vector <DOCKMol> temp;
         for (int i=0; i<tmp_parents.size(); i++){
             temp.push_back(tmp_parents[i]);
         }
         for (int i=0; i<scored_generation.size(); i++){
             temp.push_back(scored_generation[i]);
         }

         // If using niching
         if ( ga_niching ){
            // Compute rank and crowding distances
            niche_sharing( temp, score );
         
            // Determine which score is being used
            if ( score.c_desc.use_score == true || score.c_mg_nrg.use_score == true ){
               sort( temp.begin(), temp.end(), fitness_sort );
            }
         }

         if ( !ga_niching ){
            // Sort parents
            sort(temp.begin(), temp.end(), mol_sort);
         }

         //mc_metropolis(temp, ga_ensemble_size);
         true_mc_metropolis(tmp_parents,ga_ensemble_size/2);

         //Clear
         temp.clear();
      }
      else{
        for (int i=0; i<tmp_parents.size(); i++){
           parents.push_back(tmp_parents[i]);
        }
         
         for (int i=0; i<scored_generation.size(); i++){
            parents.push_back(scored_generation[i]);
         }
         // Fitness score if being used
         if ( ga_niching ){
            if ( ga_niche_sharing ){
               niche_sharing( parents, score );
            }
            else{
               niche_crowding( parents, score );
            }
         }
      }
   }

   // Sort parents
   sort(parents.begin(), parents.end(), mol_sort);


   // Clear necessary vector
   tmp_parents.clear();
   //scored_generation.clear();
  
   return;
} // end GA_Recomb:: selection_metropolis()





// +++++++++++++++++++++++++++++++++++++++++
// Activate all of the atoms AND bonds in a given DOCKMol - called in sample_minimized_torsions 
void
GA_Recomb::mc_metropolis( std::vector <DOCKMol> & tmp, int num_mols )
{
   Trace trace( "GA_Recomb::mc_metropolis()" );
   // First, select a reference fitness (the mean of all ensemble energies)
   float ave = 0.00;

   // If using niche_sharing
   if ( ga_niching ){
      for (int i=0; i<tmp.size(); i++){
          ave = ave + tmp[i].fitness;
      }
   }
   else{
      for (int i=0; i<tmp.size(); i++){
          ave = ave + tmp[i].current_score;
      }
   }

   float fitness_ref = ave/tmp.size();

   // Keep all molecules with fitness above mean
   // Keep track of all mols that have been selected
   INTVec selected;
 
   for (int i=0; i<tmp.size(); i++){
       if ( ga_niching ){
          if ( tmp[i].fitness < fitness_ref ){
             parents.push_back(tmp[i]);
             selected.push_back(1);
          }
          else{
            selected.push_back(0);
         }
       }
       // If not using niching
       else{
          if ( tmp[i].current_score < fitness_ref ){
             parents.push_back(tmp[i]);
             selected.push_back(1);
          }
          else{
            selected.push_back(0);
         }
      }
   }
   //DELETE cout << "number of mols above the mean " << parents.size() << endl;
   int initial_num = parents.size();


   // For the remainder of the molecules, maintain a few based on probability
   // P = e^(-(Fmol - Fref)/(k*T)) > rand (random number between 0 and 1)

   //DELETE cout << " num mols " << num_mols << " parents " << parents.size() << endl;

   int cycle_counter = 0;
   while (parents.size() < num_mols){
      //cout << " num mols " << num_mols << " parents " << parents.size() << endl;
      for (int i=initial_num-1; i<tmp.size(); i++){
         if (selected[i] == 0){

            float prob = 0;
            if ( ga_niching ){
               prob = exp(-(tmp[i].fitness - fitness_ref)/(ga_mc_k*new_temp));
               //DELTE cout << "prob " << prob << endl;
            }
            else{
               prob = exp(-(tmp[i].current_score - fitness_ref)/(ga_mc_k*new_temp));
               //DELTE cout << "prob " << prob << endl;
            }

            float rand_prob = ((float) rand() / RAND_MAX);
            //DELETE cout << " rand prob " << rand_prob << endl;

            if ( prob > rand_prob ){
               parents.push_back(tmp[i]);
            }
            if (parents.size() == num_mols){ break;}
         }
      }
      cycle_counter++;
      //DELETE cout << " cycle counter " << cycle_counter << endl;
      if ( cycle_counter == 3 ) { break;}
   }
   return;
} // end GA_Recomb::mc_metropolis()





// +++++++++++++++++++++++++++++++++++++++++
// Activate all of the atoms AND bonds in a given DOCKMol - called in sample_minimized_torsions 
void
GA_Recomb::true_mc_metropolis( std::vector <DOCKMol> & tmp, int num_mols )
{
   Trace trace( "GA_Recomb::true_mc_metropolis()" );
   // First, select a reference fitness (the mean of all ensemble energies)
   float ave = 0.00;

   // If using niche_sharing
   if ( ga_niching ){
      for (int i=0; i<tmp.size(); i++){
          ave = ave + tmp[i].fitness;
      }
   }
   else{
      for (int i=0; i<tmp.size(); i++){
          ave = ave + tmp[i].current_score;
      }
   }

   float fitness_ref = ave/tmp.size();

   // Keep all molecules with fitness above mean
   // Keep track of all mols that have been selected
   INTVec selected;
 
   for (int i=0; i<tmp.size(); i++){
      selected.push_back(0);
   }


   // For the remainder of the molecules, maintain a few based on probability
   // P = e^(-(Fmol - Fref)/(k*T)) > rand (random number between 0 and 1)

   int cycle_counter = 0;
   int starting_size = parents.size();
   while (parents.size()-starting_size < num_mols){
      //cout << " num mols " << num_mols << " parents " << parents.size() << endl;
      for (int i=0; i<tmp.size(); i++){
         if (selected[i] == 0){

            float prob = 0;
            if ( ga_niching ){
               prob = exp(-(tmp[i].fitness - fitness_ref)/(ga_mc_k*new_temp));
               //DELTE cout << "prob " << prob << endl;
            }
            else{
               prob = exp(-(tmp[i].current_score - fitness_ref)/(ga_mc_k*new_temp));
               //DELTE cout << "prob " << prob << endl;
            }

            if ( prob > 1 ){
                parents.push_back(tmp[i]);
                continue;
            }
            float rand_prob = ((float) rand() / RAND_MAX);
            //DELETE cout << " rand prob " << rand_prob << endl;

            if ( prob > rand_prob ){
               parents.push_back(tmp[i]);
            }
            if (parents.size() == num_mols){ break;}
         }
      }
      cycle_counter++;
      //DELETE cout << " cycle counter " << cycle_counter << endl;
      if ( cycle_counter == 4 ) { break;}
   }
   return;
} // end GA_Recomb::true_mc_metropolis()





            //////////////////////////////////////////
            //****MULTIOBJECTIVE FUNCTIONS *****//
            /////////////////////////////////////////




// +++++++++++++++++++++++++++++++++++++++++
// Function to perform niche fitness sharing method if using descriptor or mg 
// Numerical code for each scoring objective
void
GA_Recomb::niche_sharing( std::vector <DOCKMol> & mol,  Master_Score & score )
{
   Trace trace( "GA_Recomb::niche_sharing()" );
   //DELETE cout << "Inside niche_sharing" <<endl;
   
   // If vector is empty, return - CS 05/30/2018
   if (mol.size() == 0){
      return;
   }

   if ( mol.size() < 2 ){
      mol[0].fitness = -1;
   }
   else{
   // Reset the fitness of the molecules to zero
   for (int i=0; i<mol.size(); i++){
       mol[i].fitness = 0.00;
   }

   cout << "UPDATING FITNESS" <<endl;
   // Determine which descriptors are being using
   // STEP 1: If using descriptor score
   //DELETE cout << "DESC " << score.c_desc.use_score << " MG " << score.c_mg_nrg.use_score <<endl;
   if (score.c_desc.use_score == true){

      // Rank by each descriptor
      // A: GRID
      if (score.c_desc.desc_use_nrg){
         // Sort by grid energy
         sort(mol.begin(), mol.end(), grid_sort);
         // Change score range: 0-1
         score_scaling(mol, 0, score);
         // Assign crowding distance
         std::vector < std::vector <float> > dmatrix(mol.size(), std::vector<float>(mol.size())); 
         dmatrix = distance_matrix(mol, 0, score);
         float sigma = niche_radius(mol, dmatrix);
         sharing_function(mol, 0, dmatrix, sigma, score);
      }

      // B: MG
      if (score.c_desc.desc_use_mg_nrg){

         // 1: Sort & crowd by MG
         // Sort by mg energy
         sort(mol.begin(), mol.end(), mg_sort);
         // Change score range: 0-1
         score_scaling(mol, 1, score);
         // Assign crowding distance
         std::vector < std::vector <float> > dmatrix(mol.size(), std::vector<float>(mol.size())); 
         dmatrix = distance_matrix(mol, 1, score);
         float sigma = niche_radius(mol, dmatrix);
         sharing_function(mol, 1, dmatrix, sigma, score);
         // 2: Sort & crowd by fps 
         // Sort by fps energy
         sort(mol.begin(), mol.end(), fps_sort);
         // Change score range: 0-1
         score_scaling(mol, 3, score);
         // Assign crowding distance
         dmatrix = distance_matrix(mol, 3, score);
         sigma = niche_radius(mol, dmatrix);
         sharing_function(mol, 3, dmatrix, sigma, score);
      }
      
      // C: Cont Energy
      if (score.c_desc.desc_use_cont_nrg){
         // Sort by cont energy
         sort(mol.begin(), mol.end(), cont_sort);
         // Change score range: 0-1
         score_scaling(mol, 2, score);
         // Assign crowding distance
         std::vector < std::vector <float> > dmatrix(mol.size(), std::vector<float>(mol.size())); 
         dmatrix = distance_matrix(mol, 2, score);
         float sigma = niche_radius(mol, dmatrix);
         sharing_function(mol, 2, dmatrix, sigma, score);
      }

      // D: Footprint
      if (score.c_desc.desc_use_fps){
         // Sort by fps energy
         sort(mol.begin(), mol.end(), fps_sort);
         // Change score range: 0-1
         score_scaling(mol, 3, score);
         // Assign crowding distance
         std::vector < std::vector <float> > dmatrix(mol.size(), std::vector<float>(mol.size())); 
         dmatrix = distance_matrix(mol, 3, score);
         float sigma = niche_radius(mol, dmatrix);
         sharing_function(mol, 3, dmatrix, sigma, score);
      }

      // E: Pharmacophore
      if (score.c_desc.desc_use_ph4){
         // Sort by ph4 energy
         sort(mol.begin(), mol.end(), ph4_sort);
         // Change score range: 0-1
         score_scaling(mol, 4, score);
         // Assign crowding distance
         std::vector < std::vector <float> > dmatrix(mol.size(), std::vector<float>(mol.size())); 
         dmatrix = distance_matrix(mol, 4, score);
         float sigma = niche_radius(mol, dmatrix);
         sharing_function(mol, 4, dmatrix, sigma, score);
      }

      // F: Tanimoto
      if (score.c_desc.desc_use_tan){
         // Sort by tan energy
         sort(mol.begin(), mol.end(), tan_sort);
         // Change score range: 1-0
         score_scaling(mol, 5, score);
         // Assign crowding distance
         std::vector < std::vector <float> > dmatrix(mol.size(), std::vector<float>(mol.size())); 
         dmatrix = distance_matrix(mol, 5, score);
         float sigma = niche_radius(mol, dmatrix);
         sharing_function(mol, 5, dmatrix, sigma, score);
      }

      // G: Hungarian
      if (score.c_desc.desc_use_hun){
         // Sort by hun energy
         sort(mol.begin(), mol.end(), hun_sort);
         // Change score range: 0-1
         score_scaling(mol, 6, score);
         // Assign crowding distance
         std::vector < std::vector <float> > dmatrix(mol.size(), std::vector<float>(mol.size())); 
         dmatrix = distance_matrix(mol, 6, score);
         float sigma = niche_radius(mol, dmatrix);
         sharing_function(mol, 6, dmatrix, sigma, score);
      }

      // H: Volume
      if (score.c_desc.desc_use_volume){
         // Sort by vol energy
         sort(mol.begin(), mol.end(), vol_sort);
         // Change score range: -1-0
         score_scaling(mol, 7, score);
         // Assign crowding distance
         std::vector < std::vector <float> > dmatrix(mol.size(), std::vector<float>(mol.size())); 
         dmatrix = distance_matrix(mol, 7, score);
         float sigma = niche_radius(mol, dmatrix);
         sharing_function(mol, 7, dmatrix, sigma, score);
      }
   } 

   else if (score.c_mg_nrg.use_score == true){
      // Rank by each descriptor
      // A: MG
      if (score.c_desc.desc_use_mg_nrg){
         // Sort by mg energy
         sort(mol.begin(), mol.end(), mg_sort);
         // Change score range: 0-1
         score_scaling(mol, 1, score);
         // Assign crowding distance
         std::vector < std::vector <float> > dmatrix(mol.size(), std::vector<float>(mol.size())); 
         dmatrix = distance_matrix(mol, 1, score);
         float sigma = niche_radius(mol, dmatrix);
         sharing_function(mol, 1, dmatrix, sigma, score);
      }
      // B: Footprint
      if (score.c_desc.desc_use_fps){
         // Sort by fps energy
         sort(mol.begin(), mol.end(), fps_sort);
         // Change score range: 0-1
         score_scaling(mol, 3, score);
         // Assign crowding distance
         std::vector < std::vector <float> > dmatrix(mol.size(), std::vector<float>(mol.size())); 
         dmatrix = distance_matrix(mol, 3, score);
         float sigma = niche_radius(mol, dmatrix);
         sharing_function(mol, 3, dmatrix, sigma, score);
      }
   }
   }
   return;
} // end GA_Recomb::niche_sharing()




// +++++++++++++++++++++++++++++++++++++++++
// Define distance matrix for each molecule for one scoring function 
std::vector <std::vector <float> >
GA_Recomb::distance_matrix( std::vector <DOCKMol> & mol, int objective, Master_Score & score )
{
   Trace trace( "GA_Recomb::distance_matrix()" );
   // Define a matrix to determine whether a pair has already been used
   std::vector < std::vector <float> > dist_matrix(mol.size(), std::vector <float> (mol.size()));
   // Intialize the matrix to -1 so the same molecule can be distinguished from 0 distance
   for (int i=0; i<mol.size(); i++){
       for (int j=0; j<mol.size(); j++){
           dist_matrix[i][j] = -1;
       }
   }
       
   // Prepare for each objective
   // A: GRID
   if (objective == 0){
      cout << "For objective 0 (grid):" <<endl;
      /*DELETE for (int i=0; i<mol.size(); i++){
         cout <<"GRID score " << mol[i].score_nrg <<endl;
      }*/
      // Calculate Euc distance between each mol + save to matrix 
      // For each mol i (except last molecule)
      for (int i=0; i<mol.size()-1; i++){
          // For each mol, excluding itself
          for (int j=i+1; j<mol.size(); j++){
              dist_matrix[i][j] = dist_matrix[j][i] = fabs( (mol[i].score_nrg - mol[j].score_nrg) );
          }
      }
   }

   // B: MG
   else if (objective == 1){
      cout << "For objective 1 (mg):" <<endl;
      // Calculate Euc distance between each mol + save to matrix 
      // For each mol i (except last molecule)
      for (int i=0; i<mol.size()-1; i++){
          // For each mol, excluding itself
          for (int j=i+1; j<mol.size(); j++){
              dist_matrix[i][j] = dist_matrix[j][i] = fabs( (mol[i].score_mg_nrg - mol[j].score_mg_nrg) );
          }
      }
   }
   // C: Cont Energy
   else if (objective == 2){
      cout << "For objective 2 (cont):" <<endl;
      // Calculate Euc distance between each mol + save to matrix 
      // For each mol i (except last molecule)
      for (int i=0; i<mol.size()-1; i++){
          // For each mol, excluding itself
          for (int j=i+1; j<mol.size(); j++){
              dist_matrix[i][j] = dist_matrix[j][i] = fabs( (mol[i].score_cont_nrg - mol[j].score_cont_nrg) );
          }
      }
   }
   // D: Footprint
   else if (objective == 3){
      cout << "For objective 3 (fps):" <<endl;
      // Calculate Euc distance between each mol + save to matrix 
      // For each mol i (except last molecule)
      for (int i=0; i<mol.size()-1; i++){
          // For each mol, excluding itself
          for (int j=i+1; j<mol.size(); j++){
              dist_matrix[i][j] = dist_matrix[j][i] = fabs( (mol[i].score_fps - mol[j].score_fps) );
          }
      }
   }
   // E: Pharmacophore
   else if (objective == 4){
      cout << "For objective 4 (ph4):" <<endl;
      // Calculate Euc distance between each mol + save to matrix 
      // For each mol i (except last molecule)
      for (int i=0; i<mol.size()-1; i++){
          // For each mol, excluding itself
          for (int j=i+1; j<mol.size(); j++){
              dist_matrix[i][j] = dist_matrix[j][i] = fabs( (mol[i].score_ph4 - mol[j].score_ph4) );
          }
      }
   }
   // F: Tanimoto
   else if (objective == 5){
      cout << "For objective 5 (tan):" <<endl;
      // Calculate Euc distance between each mol + save to matrix 
      // For each mol i (except last molecule)
      for (int i=0; i<mol.size()-1; i++){
          // For each mol, excluding itself
          for (int j=i+1; j<mol.size(); j++){
              dist_matrix[i][j] = dist_matrix[j][i] = fabs( (mol[i].score_tan - mol[j].score_tan) );
          }
      }
   }
   // G: Hungarian
   if (objective == 6){
      cout << "For objective 6 (hun):" <<endl;
      // Calculate Euc distance between each mol + save to matrix 
      // For each mol i (except last molecule)
      for (int i=0; i<mol.size()-1; i++){
          // For each mol, excluding itself
          for (int j=i+1; j<mol.size(); j++){
              dist_matrix[i][j] = dist_matrix[j][i] = fabs( (mol[i].score_hun - mol[j].score_hun) );
          }
      }
   }
   // H: Volume
   if (objective == 7){
      cout << "For objective 7 (vol):" <<endl;
      // Calculate Euc distance between each mol + save to matrix 
      // For each mol i (except last molecule)
      for (int i=0; i<mol.size()-1; i++){
          // For each mol, excluding itself
          for (int j=i+1; j<mol.size(); j++){
              dist_matrix[i][j] = dist_matrix[j][i] = fabs( (mol[i].score_vol - mol[j].score_vol) );
          }
       }
   }

   // Print dist_matrixi DELETE
   /*for (int i=0; i<mol.size(); i++){
       for (int j=0; j<mol.size(); j++){
           cout << dist_matrix[i][j] << " ";
       } cout <<endl;
   }cout <<endl;*/
   
   return dist_matrix;
}// end GA_Recomb::distance_matrix()




// +++++++++++++++++++++++++++++++++++++++++
// Define the niche radius - threshold of similarity to be part of a niche 
// Sigma = average of min abs (distance) per molecule
float
GA_Recomb::niche_radius( std::vector <DOCKMol> & mol,  std::vector <std::vector <float> > dmatrix )
{
   Trace trace( "GA_Recomb::niche_radius()" );
   //DELETE cout <<"Inside niche radius"<<endl;
   // Sum of min distances
   float total_min_dist = 0;

   // The min distance for each molecule should be between the current molecule and either
   // molecule +1 or molecule -1, since the molecules are sorted by score
   // For each molecule
   for (int i=0; i<mol.size(); i++){
       // Placeholder for min distance
       float min = 0;

       // For the first and last molecule
       if ( i == 0){
          // Min distancce is between molecule 0 + 1
          min = dmatrix[i][i+1];
       }
       else if ( i == (mol.size()-1) ){
          min = dmatrix[i][i-1];
       }
       // For the remainder id which distance is smaller
       else {
         if ( dmatrix[i][i+1] < dmatrix[i][i-1] ){
            min = dmatrix[i][i+1];
         }
         else{ 
            min = dmatrix[i][i-1]; 
         }
       }
       // Update total min
       total_min_dist += min;
   }

   // Calculare sigma
   float sigma = total_min_dist/mol.size();
   //DELETE cout <<"sigma " << sigma <<endl;
   return sigma;
}//end GA_Recomb::niche_radius




// +++++++++++++++++++++++++++++++++++++++++
// Calculates the similarity between a molecule and the remaining populatin 
// For this function, the best fitness is the abs value (largest number)
void
GA_Recomb::sharing_function( std::vector <DOCKMol> & mol, int objective, std::vector <std::vector <float> > dmatrix, 
                             float sigma, Master_Score & score )
{
   Trace trace( "GA_Recomb::sharing_function()" );
   for (int i=0; i<mol.size(); i++){
       // Niche count for each molecule
       float niche_count = 0.00;
       int num_below_sigma = 0;
       for (int j=0; j<mol.size(); j++){
           if ( i != j ){
             // If distance < sigma
             if ( dmatrix[i][j] < sigma ){
                // Compute similarity
                //cout << dmatrix[i][j]/sigma <<endl;
                niche_count += 1 - ( dmatrix[i][j] / sigma );
                num_below_sigma++;
                //DELETE cout << "niche_count " << niche_count << endl;
             }
          }
      }

      // If there are no similar molecules, do not perform
      // Would divide by 0
      //DELETE cout << "num_below_sigma " << num_below_sigma <<endl;

      // Calculate new shared fitness
      // A: GRID
      if (objective == 0){
         if (num_below_sigma == 0){
            // Update mol.fitness
            mol[i].fitness += mol[i].score_nrg;
            continue;
         }
         //cout << "For objective 0 (grid):" <<endl;
         //DELETE cout << "mol " << i << " fitness initial " << mol[i].fitness << " score " << mol[i].score_nrg<<endl;
         float new_fitness = mol[i].score_nrg/niche_count;
         
         if ( new_fitness > mol[i].score_nrg ){
            mol[i].fitness += new_fitness;
         }
         else{
            mol[i].fitness += mol[i].score_nrg;
         }
         //DELETE cout << "mol " << i << " new fitness " << new_fitness <<" final fitness " << mol[i].fitness <<endl;
      }

      // B: MG
      else if (objective == 1){
         if (num_below_sigma == 0){
            // Update mol.fitness
         cout << "For objective 1 (mg): nothing similar" <<endl;
         cout << "mol " << i << " fitness initial " << mol[i].fitness << " score " << mol[i].score_mg_nrg<<endl;
            mol[i].fitness += mol[i].score_mg_nrg;
         cout << "mol " << i << " final fitness " << mol[i].fitness <<endl;
            continue;
         }
         cout << "For objective 1 (mg):" <<endl;
         cout << "mol " << i << " fitness initial " << mol[i].fitness << " score " << mol[i].score_mg_nrg<<endl;
         float new_fitness = mol[i].score_mg_nrg/niche_count;
         if ( new_fitness > mol[i].score_mg_nrg ){
            mol[i].fitness += new_fitness;
         }
         else{
            mol[i].fitness += mol[i].score_mg_nrg;
         }
         cout << "mol " << i << " new fitness " << new_fitness <<" final fitness " << mol[i].fitness <<endl;
      }

      // C: Cont Energy
      else if (objective == 2){
         if (num_below_sigma == 0){
            // Update mol.fitness
            mol[i].fitness += mol[i].score_cont_nrg;
            continue;
         }
         //cout << "For objective 2 (cont):" <<endl;
         //DELETE cout << "mol " << i << " fitness initial " << mol[i].fitness <<endl;
         float new_fitness = mol[i].score_cont_nrg/niche_count;
         if ( new_fitness > mol[i].score_cont_nrg ){
            mol[i].fitness += new_fitness;
         }
         else{
            mol[i].fitness += mol[i].score_cont_nrg;
         }
         //DELETE cout << "mol " << i << " new fitness " << new_fitness <<" final fitness " << mol[i].fitness <<endl;
      }

      // D: Footprint
      else if (objective == 3){
         if (num_below_sigma == 0){
            // Update mol.fitness
         cout << "For objective 3 (fps): nothing similar " <<endl;
         cout << "mol " << i << " fitness initial " << mol[i].fitness << " score " << mol[i].score_fps <<endl;
            mol[i].fitness += mol[i].score_fps;
         cout << "mol " << i << " final fitness " << mol[i].fitness <<endl;
            continue;
         }
         cout << "For objective 3 (fps):" <<endl;
         cout << "mol " << i << " fitness initial " << mol[i].fitness << " score " << mol[i].score_fps <<endl;
         float new_fitness= mol[i].score_fps/niche_count;
         if ( new_fitness > mol[i].score_fps ){
            mol[i].fitness += new_fitness;
         }
         else{
            mol[i].fitness += mol[i].score_fps;
         }
         cout << "mol " << i << " new fitness " << new_fitness <<" final fitness " << mol[i].fitness <<endl;
      }

      // E: Pharmacophore
      else if (objective == 4){
         if (num_below_sigma == 0){
         cout << "For objective 4 (ph4): nothing similar " <<endl;
         cout << "mol " << i << " fitness initial " << mol[i].fitness << " score " << mol[i].score_ph4<<endl;
            // Update mol.fitness
            mol[i].fitness += mol[i].score_ph4;
         cout << "mol " << i << " final fitness " << mol[i].fitness <<endl;
            continue;
         }
         cout << "For objective 4 (ph4):" <<endl;
         cout << "mol " << i << " fitness initial " << mol[i].fitness << " score " << mol[i].score_ph4<<endl;
         float new_fitness = mol[i].score_ph4/niche_count;
         if ( new_fitness > mol[i].score_ph4 ){
            mol[i].fitness += new_fitness;
         }
         else{
            mol[i].fitness += mol[i].score_ph4;
         }
         cout << "mol " << i << " new fitness " << new_fitness <<" final fitness " << mol[i].fitness <<endl;
      }

      // F: Tanimoto
      else if (objective == 5){
         if (num_below_sigma == 0){
            // Update mol.fitness
            mol[i].fitness += mol[i].score_tan;
            continue;
         }
         //cout << "For objective 5 (tan):" <<endl;
         //DELETE cout << "mol " << i << " fitness initial " << mol[i].fitness <<endl;
         float new_fitness = mol[i].score_tan/niche_count;
         if ( new_fitness > mol[i].score_tan ){
            mol[i].fitness += new_fitness;
         } 
         else{
            mol[i].fitness += mol[i].score_tan;
         }
         //DELETE cout << "mol " << i << " new fitness " << new_fitness <<" final fitness " << mol[i].fitness <<endl;
      }

      // G: Hungarian
      if (objective == 6){
         if (num_below_sigma == 0){
            // Update mol.fitness
            mol[i].fitness += mol[i].score_hun;
            continue;
         }
         //cout << "For objective 6 (hun):" <<endl;
         //DELETE cout << "mol " << i << " fitness initial " << mol[i].fitness <<endl;
         float new_fitness= mol[i].score_hun/niche_count;
         if ( new_fitness > mol[i].score_hun ){
            mol[i].fitness += new_fitness;
         }
         else{
            mol[i].fitness += mol[i].score_hun;
         }
         //DELETE cout << "mol " << i << " new fitness " << new_fitness <<" final fitness " << mol[i].fitness <<endl;
      }

      // H: Volume
      if (objective == 7){
         if (num_below_sigma == 0){
         cout << "For objective 7 (vol): nothing similar " <<endl;
         cout << "mol " << i << " fitness initial " << mol[i].fitness << "score " << mol[i].score_vol<<endl;
            // Update mol.fitness
            mol[i].fitness += mol[i].score_vol;
         cout << "mol " << i << " final fitness " << mol[i].fitness <<endl;
            continue;
         }
         cout << "For objective 7 (vol):" <<endl;
         cout << "mol " << i << " fitness initial " << mol[i].fitness << "score " << mol[i].score_vol<<endl;
         float new_fitness= mol[i].score_vol/niche_count;
         if ( new_fitness > mol[i].score_vol ){
            mol[i].fitness += new_fitness;
         }
         else{
            mol[i].fitness += mol[i].score_vol;
         }
         cout << "mol " << i << " new fitness " << new_fitness <<" final fitness " << mol[i].fitness <<endl;
      }
   }
   return;
}//end GA_Recomb::sharing_function()




// +++++++++++++++++++++++++++++++++++++++++
// Score normalization function - scales scores between -1 - 0
// For the shared fitness equation the absolute value of highest number is the most fit
// so all of the scores are converted from range 0 - 1 to -1 - 0
// Send this function a presorted list
void
GA_Recomb::score_scaling( std::vector <DOCKMol> & mol, int objective, Master_Score & score )
{
   Trace trace( "GA_Recomb::score_scaling()" );
   // Identify the min and max score for the ensemble for current scoring type   
   float min_score = 0.00;
   float max_score = 0.00;

   // Calculate new shared fitness
   // A: GRID
   if (objective == 0){
      //cout << "For objective 0 (grid):" <<endl;
      // Set initial score as min and max
      min_score = mol[0].score_nrg;
      max_score = mol[mol.size()-1].score_nrg;

      // Change the range of all scores
      for (int i=0; i<mol.size(); i++){
          //DELETE cout << "mol " << i << " score " << mol[i].score_nrg;
          mol[i].score_nrg = (( mol[i].score_nrg - min_score )/( max_score - min_score )) - 1;
          //DELETE cout << " new score " << mol[i].score_nrg <<endl;
      }
   }

   // B: MG
   else if (objective == 1){
      //cout << "For objective 1 (mg):" <<endl;
      // Set initial score as min and max
      min_score = mol[0].score_mg_nrg;
      max_score = mol[mol.size()-1].score_mg_nrg;

      // Change the range of all scores
      for (int i=0; i<mol.size(); i++){
          //DELETE cout << "mol " << i << " score " << mol[i].score_mg_nrg;
          cout << "mol " << i << " score_mg " << mol[i].score_mg_nrg << " min " << min_score << " max " << max_score<< "score - min " << mol[i].score_mg_nrg - min_score << " denom " << max_score - min_score;
          mol[i].score_mg_nrg = (( mol[i].score_mg_nrg - min_score )/( max_score - min_score )) - 1;
          cout << " new score " << mol[i].score_mg_nrg <<endl;
      }
   }
   // C: Cont Energy
   else if (objective == 2){
      //cout << "For objective 2 (cont):" <<endl;
      // Set initial score as min and max
      min_score = mol[0].score_cont_nrg;
      max_score = mol[mol.size()-1].score_cont_nrg;

      // Change the range of all scores
      for (int i=0; i<mol.size(); i++){
          //DELETE cout << "mol " << i << " score " << mol[i].score_cont_nrg<<endl;
          mol[i].score_cont_nrg = (( mol[i].score_cont_nrg - min_score )/( max_score - min_score )) - 1;
          //DELETE cout << "mol " << i << " new score " << mol[i].score_cont_nrg <<endl;
      }
   }
   // D: Footprint
   else if (objective == 3){
      //cout << "For objective 3 (fps):" <<endl;
      // Set initial score as min and max
      min_score = mol[0].score_fps;
      max_score = mol[mol.size()-1].score_fps;

      // Change the range of all scores
      for (int i=0; i<mol.size(); i++){
          //DELETE cout << "mol " << i << " score " << mol[i].score_fps;
          cout << "mol " << i << " score_fps " << mol[i].score_fps << " min " << min_score << " max " << max_score<< "score - min " << mol[i].score_fps - min_score << " denom " << max_score - min_score;
          mol[i].score_fps = (( mol[i].score_fps - min_score )/( max_score - min_score )) - 1;
          //DELETE cout << " new score " << mol[i].score_fps <<endl;
          cout << " new score " << mol[i].score_fps <<endl;
      }
   }
   // E: Pharmacophore
   else if (objective == 4){
      //cout << "For objective 4 (ph4):" <<endl;
      // Set initial score as min and max
      min_score = mol[0].score_ph4;
      max_score = mol[mol.size()-1].score_ph4;

      // Change the range of all scores
      for (int i=0; i<mol.size(); i++){
          //DELETE cout << "mol " << i << " score " << mol[i].score_ph4;
          cout << "mol " << i << " score_ph4 " << mol[i].score_ph4 << " min " << min_score << " max " << max_score<< "score - min " << mol[i].score_ph4 - min_score << " denom " << max_score - min_score;
          mol[i].score_ph4 = (( mol[i].score_ph4 - min_score )/( max_score - min_score )) - 1;
          //DELETE cout << " new score " << mol[i].score_ph4 <<endl;
          cout << " new score " << mol[i].score_ph4 <<endl;
      }
   }
   // F: Tanimoto
   else if (objective == 5){
      //cout << "For objective 5 (tan):" <<endl;
      // Set initial score as min and max
      max_score = mol[0].score_tan;
      min_score = mol[mol.size()-1].score_tan;
 
      // This scoring function is already between 0 - 1 but 1 is best
      // Make all scores negative
      for (int i=0; i<mol.size(); i++){
          //DELETE cout << "mol " << i << " score " << mol[i].score_tan<<endl;
          //mol[i].score_tan = -mol[i].score_tan;
          mol[i].score_tan = -(( mol[i].score_tan - min_score )/( max_score - min_score ));
          //DELETE cout << "mol " << i << " new score " << mol[i].score_tan <<endl;
      }
   }
   // G: Hungarian
   if (objective == 6){
      //cout << "For objective 6 (hun):" <<endl;
      // Set initial score as min and max
      min_score = mol[0].score_hun;
      max_score = mol[mol.size()-1].score_nrg;

      // Change the range of all scores
      for (int i=0; i<mol.size(); i++){
          //DELETE cout << "mol " << i << " score " << mol[i].score_hun<<endl;
          mol[i].score_hun = (( mol[i].score_hun - min_score )/( max_score - min_score )) - 1;
          //DELETE cout << "mol " << i << " new score " << mol[i].score_hun <<endl;
      }
   }
   // H: Volume
   if (objective == 7){
      //cout << "For objective 7 (vol):" <<endl;

      // Set initial score as min and max
      max_score = mol[0].score_vol;
      min_score = mol[mol.size()-1].score_vol;

      // This scoring function is already between 0 - 1 but 1 is best
      // Make all scores negative
      for (int i=0; i<mol.size(); i++){
          cout << "mol " << i << " score_vol " << mol[i].score_vol << " min " << min_score << " max " << max_score<< "score - min " << mol[i].score_vol - min_score << " denom " << max_score - min_score;
          //mol[i].score_vol = -mol[i].score_vol;
          mol[i].score_vol = -(( mol[i].score_vol - min_score )/( max_score - min_score ));
          //DELETE cout << " new score " << mol[i].score_vol <<endl;
          cout << " new score " << mol[i].score_vol <<endl;
      }
   }
   return;
}// end GA_Recomb::score_scaling()




// +++++++++++++++++++++++++++++++++++++++++
// Function to perform niche crowding method if using descriptor or mg 
// Numerical code for each scoring objective
void
GA_Recomb::niche_crowding( std::vector <DOCKMol> & mol,  Master_Score & score )
{
   Trace trace( "GA_Recomb::niche_crowding()" );

   // If offspring is empty leave niching function - 2018.06.05 LEP
   if (mol.size() == 0){
       return;
   }
   // Reset the rank and crowding distance of the parents
   for (int i=0; i<mol.size(); i++){
       mol[i].rank = 0;
       mol[i].crowding_dist = 0.00;
   }

   // Determine which descriptors are being using
   // STEP 1: If using descriptor score
   if (score.c_desc.use_score == true){

      // Rank by each descriptor
      // A: GRID
      if (score.c_desc.desc_use_nrg){
         // Sort by grid energy
         sort(mol.begin(), mol.end(), grid_sort);
         // Increment rank 
         increment_rank(mol);
         // Assign crowding distance
         crowding_dist(mol, 0, score);
      }

      // B: MG
      if (score.c_desc.desc_use_mg_nrg){

         // 1: Sort & crowd by MG
         // Sort by mg energy
         sort(mol.begin(), mol.end(), mg_sort);
         // Increment rank 
         increment_rank(mol);
         // Assign crowding distance
         crowding_dist(mol, 1, score);

         // 2: Sort & crowd by fps 
         // Sort by fps energy
         sort(mol.begin(), mol.end(), fps_sort);
         // Increment rank 
         increment_rank(mol);
         // Assign crowding distance
         crowding_dist(mol, 3, score);
      }
      
      // C: Cont Energy
      if (score.c_desc.desc_use_cont_nrg){
         // Sort by cont energy
         sort(mol.begin(), mol.end(), cont_sort);
         // Increment rank 
         increment_rank(mol);
         // Assign crowding distance
         crowding_dist(mol, 2, score);
      }

      // D: Footprint
      if (score.c_desc.desc_use_fps){
         // Sort by fps energy
         sort(mol.begin(), mol.end(), fps_sort);
         // Increment rank 
         increment_rank(mol);
         // Assign crowding distance
         crowding_dist(mol, 3, score);
      }

      // E: Pharmacophore
      if (score.c_desc.desc_use_ph4){
         // Sort by ph4 energy
         sort(mol.begin(), mol.end(), ph4_sort);
         // Increment rank 
         increment_rank(mol);
         // Assign crowding distance
         crowding_dist(mol, 4, score);
      }

      // F: Tanimoto
      if (score.c_desc.desc_use_tan){
         // Sort by tan energy
         sort(mol.begin(), mol.end(), tan_sort);
         // Increment rank 
         increment_rank(mol);
         // Assign crowding distance
         crowding_dist(mol, 5, score);
      }

      // G: Hungarian
      if (score.c_desc.desc_use_hun){
         // Sort by hun energy
         sort(mol.begin(), mol.end(), hun_sort);
         // Increment rank 
         increment_rank(mol);
         // Assign crowding distance
         crowding_dist(mol, 6, score);
      }

      // H: Volume
      if (score.c_desc.desc_use_volume){
         // Sort by vol energy
         sort(mol.begin(), mol.end(), vol_sort);
         // Increment rank 
         increment_rank(mol);
         // Assign crowding distance
         crowding_dist(mol, 7, score);
      }
   } 

   else if (score.c_mg_nrg.use_score == true){
      // Rank by each descriptor
      // A: MG
      if (score.c_desc.desc_use_mg_nrg){
         // Sort by mg energy
         sort(mol.begin(), mol.end(), mg_sort);
         // Increment rank 
         increment_rank(mol);
         // Assign crowding distance
         crowding_dist(mol, 1, score);
      }
      // B: Footprint
      if (score.c_desc.desc_use_fps){
         // Sort by fps energy
         sort(mol.begin(), mol.end(), fps_sort);
         // Increment rank 
         increment_rank(mol);
         // Assign crowding distance
         crowding_dist(mol, 3, score);
      }
   }
   return;
} // end GA_Recomb::niche_crowding()




// +++++++++++++++++++++++++++++++++++++++++
// Function to rank molecules if using descriptor or mg 
void
GA_Recomb::increment_rank( std::vector <DOCKMol> & mol )
{  
   Trace trace( "GA_Recomb::increment_rank()" );
   // For each mol
   for (int i=0; i<mol.size(); i++){
      // Increment rank based on position
      mol[i].rank += i;
   }

   return;
} // end GA_Recomb::increment_rank()



// +++++++++++++++++++++++++++++++++++++++++
// Measures crowding distance for each descriptor
void
GA_Recomb::crowding_dist( std::vector <DOCKMol> & mol, int objective,  Master_Score & score )
{
   Trace trace( "GA_Recomb::crowding_dist()" );
   // Prepare for each objective
   // A: GRID
   if (objective == 0){
      //cout << "For objective 0 (grid):" <<endl;
      // Assign end points
      mol[0].crowding_dist += 2*fabs( (mol[0].score_nrg - mol[1].score_nrg) / 
                                   (mol[0].score_nrg - mol[mol.size()-1].score_nrg) );

      //DELETE cout <<"Mol 0 w/ rank: " << mol[0].rank << " has dist " << mol[0].crowding_dist <<endl;

      mol[mol.size()-1].crowding_dist += 2*fabs( (mol[mol.size()-2].score_nrg - mol[mol.size()-1].score_nrg) / 
                                              (mol[0].score_nrg - mol[mol.size()-1].score_nrg) );

      //DELETE cout <<"Mol " << mol.size() << " w/ rank: " << mol[mol.size()-1].rank << " has dist " << mol[mol.size()-1].crowding_dist <<endl;

      // For other molecules
      for (int i=1; i<mol.size()-1; i++){
          mol[i].crowding_dist += fabs( (mol[i+1].score_nrg - mol[i-1].score_nrg)/
                                       (mol[0].score_nrg - mol[mol.size()-1].score_nrg) );
      //DELETE cout <<"Mol " << i << " w/ rank: " << mol[i].rank << " has dist " << mol[i].crowding_dist <<endl;
      }
   }

   // B: MG
   else if (objective == 1){
      // Assign end points
      //cout << "For objective 1 (mg):" <<endl;
      mol[0].crowding_dist += 2*fabs( (mol[0].score_mg_nrg - mol[1].score_mg_nrg) / 
                                   (mol[0].score_mg_nrg - mol[mol.size()-1].score_mg_nrg) );
  
      //DELETE cout <<"Mol 0 w/ rank: " << mol[0].rank << " has dist " << mol[0].crowding_dist <<endl;
  
      mol[mol.size()-1].crowding_dist += 2*fabs( (mol[mol.size()-2].score_mg_nrg - mol[mol.size()-1].score_mg_nrg) / 
                                            (mol[0].score_mg_nrg - mol[mol.size()-1].score_mg_nrg) );
      //DELETE cout <<"Mol " << mol.size() << " w/ rank: " << mol[mol.size()-1].rank << " has dist " << mol[mol.size()-1].crowding_dist <<endl;

      // For other molecules
      for (int i=1; i<mol.size()-1; i++){
          mol[i].crowding_dist += fabs( (mol[i+1].score_mg_nrg - mol[i-1].score_mg_nrg)/
                                       (mol[0].score_mg_nrg - mol[mol.size()-1].score_mg_nrg) );
      //DELETE cout <<"Mol " << i << " w/ rank: " << mol[i].rank << " has dist " << mol[i].crowding_dist <<endl;
      }

   }
   // C: Cont Energy
   else if (objective == 2){
      // Assign end points
      mol[0].crowding_dist += 2*fabs( (mol[0].score_cont_nrg - mol[1].score_cont_nrg) / 
                                   (mol[0].score_cont_nrg - mol[mol.size()-1].score_cont_nrg) );
      mol[mol.size()-1].crowding_dist += 2*fabs( (mol[mol.size()-2].score_cont_nrg - mol[mol.size()-1].score_cont_nrg) / 
                                              (mol[0].score_cont_nrg - mol[mol.size()-1].score_cont_nrg) );

      // For other molecules
      for (int i=1; i<mol.size()-1; i++){
          mol[i].crowding_dist += fabs( (mol[i+1].score_cont_nrg - mol[i-1].score_cont_nrg)/
                                       (mol[0].score_cont_nrg - mol[mol.size()-1].score_cont_nrg) );
      }

   }
   // D: Footprint
   else if (objective == 3){
      //cout << "For objective 3 (fps):" <<endl;
      // Assign end points
      mol[0].crowding_dist += 2*fabs( (mol[0].score_fps - mol[1].score_fps) / 
                                   (mol[0].score_fps - mol[mol.size()-1].score_fps) );
      //DELETE cout <<"Mol 0 w/ rank: " << mol[0].rank << " has dist " << mol[0].crowding_dist <<endl;
      mol[mol.size()-1].crowding_dist += 2*fabs( (mol[mol.size()-2].score_fps - mol[mol.size()-1].score_fps) / 
                                              (mol[0].score_fps - mol[mol.size()-1].score_fps) );
      //DELETE cout <<"Mol " << mol.size() << " w/ rank: " << mol[mol.size()-1].rank << " has dist " << mol[mol.size()-1].crowding_dist <<endl;

      // For other molecules
      for (int i=1; i<mol.size()-1; i++){
          mol[i].crowding_dist += fabs( (mol[i+1].score_fps - mol[i-1].score_fps)/
                                       (mol[0].score_fps - mol[mol.size()-1].score_fps) );
      //DELETE cout <<"Mol " << i << " w/ rank: " << mol[i].rank << " has dist " << mol[i].crowding_dist <<endl;
      }

   }
   // E: Pharmacophore
   else if (objective == 4){
      // Assign end points
      mol[0].crowding_dist += 2*fabs( (mol[0].score_ph4 - mol[1].score_ph4) / 
                                   (mol[0].score_ph4 - mol[mol.size()-1].score_ph4) );
      mol[mol.size()-1].crowding_dist += 2*fabs( (mol[mol.size()-2].score_ph4 - mol[mol.size()-1].score_ph4) / 
                                            (mol[0].score_ph4 - mol[mol.size()-1].score_ph4) );

      // For other molecules
      for (int i=1; i<mol.size()-1; i++){
          mol[i].crowding_dist += fabs( (mol[i+1].score_ph4 - mol[i-1].score_ph4)/
                                       (mol[0].score_ph4 - mol[mol.size()-1].score_ph4) );
      }

   }
   // F: Tanimoto
   else if (objective == 5){
      // Assign end points
      mol[0].crowding_dist += 2*fabs( (mol[0].score_tan - mol[1].score_tan) / 
                                   (mol[0].score_tan - mol[mol.size()-1].score_tan) );
      mol[mol.size()-1].crowding_dist += 2*fabs( (mol[mol.size()-2].score_tan - mol[mol.size()-1].score_tan) / 
                                            (mol[0].score_tan - mol[mol.size()-1].score_tan) );

      // For other molecules
      for (int i=1; i<mol.size()-1; i++){
          mol[i].crowding_dist += fabs( (mol[i+1].score_tan - mol[i-1].score_tan)/
                                       (mol[0].score_tan - mol[mol.size()-1].score_tan) );
      }

   }
   // G: Hungarian
   if (objective == 6){
      // Assign end points
      mol[0].crowding_dist += 2*fabs( (mol[0].score_hun - mol[1].score_hun) / 
                                   (mol[0].score_hun - mol[mol.size()-1].score_hun) );
      mol[mol.size()-1].crowding_dist += 2*fabs( (mol[mol.size()-2].score_hun - mol[mol.size()-1].score_hun) / 
                                              (mol[0].score_hun - mol[mol.size()-1].score_hun) );

      // For other molecules
      for (int i=1; i<mol.size()-1; i++){
          mol[i].crowding_dist += fabs( (mol[i+1].score_hun - mol[i-1].score_hun)/
                                       (mol[0].score_hun - mol[mol.size()-1].score_hun) );
      }

   }
   // H: Volume
   if (objective == 7){
      // Assign end points
      mol[0].crowding_dist += 2*fabs( (mol[0].score_vol - mol[1].score_vol) / 
                                   (mol[0].score_vol - mol[mol.size()-1].score_vol) );
      mol[mol.size()-1].crowding_dist += 2*fabs( (mol[mol.size()-2].score_vol - mol[mol.size()-1].score_vol) / 
                                              (mol[0].score_vol - mol[mol.size()-1].score_vol) );

      // For other molecules
      for (int i=1; i<mol.size()-1; i++){
          mol[i].crowding_dist += fabs( (mol[i+1].score_vol - mol[i-1].score_vol)/
                                       (mol[0].score_vol - mol[mol.size()-1].score_vol) );
      }
   
   }
   
   return;
} // end GA_Recomb::crowding_dist()






            //////////////////////////////////////////
            //****MIS FUNCTIONS *****//
            /////////////////////////////////////////




// +++++++++++++++++++++++++++++++++++++++++
// Deactivate all of the atoms AND bonds in a given DOCKMol - called in sample_minimized_torsions 
void
GA_Recomb::deactivate_mol( DOCKMol & mol )
{
   Trace trace( "GA_Recomb::deactivate_mol()" );
   // Iterate through all atoms and set atom_active_flag to true
   for ( int i=0; i<mol.num_atoms; i++ ){
       mol.atom_active_flags[i] = false;
   }

   // Iterate through all bonds and set bond_active_flag to true
   for ( int i=0; i<mol.num_bonds; i++ ){
       mol.bond_active_flags[i] = false;
   }

   return;

} // end GA_Recomb::deactivate_mol()




// +++++++++++++++++++++++++++++++++++++++++
// Activate all of the atoms AND bonds in a given DOCKMol - called in sample_minimized_torsions 
void
GA_Recomb::activate_mol( DOCKMol & mol )
{
   Trace trace( "GA_Recomb::activate_mol()" );
   // Iterate through all atoms and set atom_active_flag to true
   for ( int i=0; i<mol.num_atoms; i++ ){
       mol.atom_active_flags[i] = true;
   }

   // Iterate through all bonds and set bond_active_flag to true
   for ( int i=0; i<mol.num_bonds; i++ ){
       mol.bond_active_flags[i] = true;
   }

   return;

} // end GA_Recomb::activate_mol()




// ++++++++++++++++++++++++++++++++++++++++
// Activate an entire vector
void
GA_Recomb::activate_vector( std::vector <DOCKMol> & mol_vec )
{
   Trace trace( "GA_Recomb::activate_vector()" );
   // Activate parents before amber_typer
   for ( int i=0; i<mol_vec.size(); i++ ){
       for ( int j=0; j<mol_vec[i].num_bonds; j++ ){
           mol_vec[i].bond_active_flags[j]=true;
       }
   }

   for ( int i=0; i<mol_vec.size(); i++ ){
       for ( int j=0; j<mol_vec[i].num_atoms; j++ ){
           mol_vec[i].atom_active_flags[j]=true;
       }
   }
   return;
} // end GA_Recomb::activate_vector()



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by energy score
bool
mol_sort(DOCKMol a, DOCKMol b)
{
   Trace trace( "GA_Recomb::mol_sort()" );
   return (a.current_score < b.current_score);
}// end GA_Recomb::mol_sort()



/*// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by energy score
bool
fragment_sort( const Fragment & a, const Fragment & b )
{
    return (a.mol.current_score < b.mol.current_score);
}// end GA_Reomb::fragment_sort()*/



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by grid score
bool
grid_sort(DOCKMol a, DOCKMol b)
{
   Trace trace( "GA_Recomb::grid_sort()" );
   return (a.score_nrg < b.score_nrg);
}// end GA_Recomb::grid_sort()



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by mg score
bool
mg_sort(DOCKMol a, DOCKMol b)
{
   Trace trace( "GA_Recomb::mg_sort()" );
   return (a.score_mg_nrg < b.score_mg_nrg);
}// end GA_Recomb::mg_sort()



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by cont score
bool
cont_sort(DOCKMol a, DOCKMol b)
{
   Trace trace( "GA_Recomb::cont_sort()" );
   return (a.score_cont_nrg < b.score_cont_nrg);
}// end GA_Recomb::cont_sort()



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by fps score
bool
fps_sort(DOCKMol a, DOCKMol b)
{
   Trace trace( "GA_Recomb::fps_sort()" );
   return (a.score_fps < b.score_fps);
}// end GA_Recomb::fps_sort()



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by ph4 score
bool
ph4_sort(DOCKMol a, DOCKMol b)
{
   Trace trace( "GA_Recomb::ph4_sort()" );
   return (a.score_ph4 < b.score_ph4);
}// end GA_Recomb::ph4_sort()



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by tan score
bool
tan_sort(DOCKMol a, DOCKMol b)
{
   Trace trace( "GA_Recomb::tan_sort()" );
   return (a.score_tan > b.score_tan);
}// end GA_Recomb::tan_sort()



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by hun score
bool
hun_sort(DOCKMol a, DOCKMol b)
{
   Trace trace( "GA_Recomb::hun_sort()" );
   return (a.score_hun < b.score_hun);
}// end GA_Recomb::hun_sort()



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by vol score
bool
vol_sort(DOCKMol a, DOCKMol b)
{
   Trace trace( "GA_Recomb::vol_sort()" );
   return (a.score_vol > b.score_vol);
}// end GA_Recomb::vol_sort()



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by rank score
bool
rank_sort(DOCKMol a, DOCKMol b)
{
   Trace trace( "GA_Recomb::rank_sort()" );
   return (a.rank < b.rank);
}// end GA_Recomb::rank_sort()



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by rank score
bool
fitness_sort(DOCKMol a, DOCKMol b)
{
   Trace trace( "GA_Recomb::fitness_sort()" );
   return (a.fitness < b.fitness);
}// end GA_Recomb::rank_sort()



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by crowding_dist score
bool
crowd_sort(DOCKMol a, DOCKMol b)
{
   Trace trace( "GA_Recomb::crowd_sort()" );
   return (a.crowding_dist >  b.crowding_dist);
}// end GA_Recomb::crowd_sort()



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by crowding_dist score
bool
compare_rank(DOCKMol a, DOCKMol b)
{
   Trace trace( "GA_Recomb::compare_rank()" );
   return ((a.rank < b.rank) || 
          ((a.rank == b.rank) && (a.crowding_dist > b.crowding_dist)));
}// end GA_Recomb::compare_rank()

// +++++++++++++++++++++++++++++++++++++++++
// Function to renumber atoms to have unique atom numbering
void
GA_Recomb::renumber_atom_numbering( DOCKMol & mol , bool dna_flag )
{
    int Hcounter = 0;
    int Ccounter = 0;
    int Ncounter = 0;
    int Ocounter = 0;
    int Scounter = 0;
    int Pcounter = 0;
    int Fcounter = 0;
    int Clcounter = 0;
    int Brcounter = 0;
   
    for (int i = 0; i < mol.num_atoms; ++i){

        if ( mol.atom_types[i] == "H" )
        {     
              Hcounter++;
              if ( mol.atom_names[i] != "DNM" || !dna_flag ){
                stringstream stuff;
                stuff << "H"<< Hcounter;
                mol.atom_names[i] = stuff.str();
                stuff.clear();
              }
        }

        else if ( mol.atom_types[i] == "C.3" || mol.atom_types[i] == "C.2" || mol.atom_types[i] == "C.1" || mol.atom_types[i] == "C.ar" || mol.atom_types[i]  == "C.cat" )
        {     
              Ccounter++;
              if ( mol.atom_names[i] != "DNM" || !dna_flag ){
                stringstream stuff;
                stuff << "C"<< Ccounter;
                mol.atom_names[i] = stuff.str();
                stuff.clear();
              }

        }

        else if ( mol.atom_types[i] == "N.4" || mol.atom_types[i] == "N.3" || mol.atom_types[i] == "N.2" || mol.atom_types[i] == "N.1" || mol.atom_types[i]  == "N.ar" ||
              mol.atom_types[i] == "N.am" || mol.atom_types[i] == "N.pl3" )
        {     
              Ncounter++; 
              if ( mol.atom_names[i] != "DNM" || !dna_flag ){
                stringstream stuff;
                stuff << "N" << Ncounter;
                mol.atom_names[i] = stuff.str();
                stuff.clear();
              }
        }

        else if ( mol.atom_types[i] == "O.3" || mol.atom_types[i] == "O.2" || mol.atom_types[i] == "O.co2" )
        {     
              Ocounter++; 
              if ( mol.atom_names[i] != "DNM" || !dna_flag ){
                stringstream stuff;
                stuff << "O" << Ocounter;
                mol.atom_names[i] = stuff.str();
                stuff.clear();
              }
        }

        else if ( mol.atom_types[i] == "S.3" || mol.atom_types[i] == "S.2" || mol.atom_types[i] == "S.O" || mol.atom_types[i] == "S.o" || mol.atom_types[i] == "S.O2" ||
              mol.atom_types[i] == "S.o2" )
        {     
              Scounter++; 
              if ( mol.atom_names[i] != "DNM" || !dna_flag ){
                stringstream stuff;
                stuff << "S" << Scounter;
                mol.atom_names[i] = stuff.str();
                stuff.clear();
              }
        }

        else if ( mol.atom_types[i] == "P.3" )
        {    
             Pcounter++; 
             if ( mol.atom_names[i] != "DNM" || !dna_flag ){
                stringstream stuff;
                stuff << "P" << Pcounter;
                mol.atom_names[i] = stuff.str();
                stuff.clear();
             }
        }

        else if ( mol.atom_types[i] == "F" )
        {    
             Fcounter++;
             if ( mol.atom_names[i] != "DNM" || !dna_flag ){
                stringstream stuff;
                stuff << "F" << Fcounter;
                mol.atom_names[i] = stuff.str();
                stuff.clear();
             }
        }

        else if ( mol.atom_types[i] == "Cl" )
        {    
             Clcounter++;
             if ( mol.atom_names[i] != "DNM" || !dna_flag ){
                stringstream stuff;
                stuff << "Cl" << Clcounter;
                mol.atom_names[i] = stuff.str();
                stuff.clear();
             }
        }

        else if ( mol.atom_types[i] == "Br" )
        {    Brcounter++; 
             if ( mol.atom_names[i] != "DNM" || !dna_flag ){
                stringstream stuff;
                stuff << "Br" << Brcounter;
                mol.atom_names[i] = stuff.str();
                stuff.clear();
             }
        }

        else{
            cout << "Atom renumbering: atom type not recognized - " << mol.atom_types[i]<< endl;
        }
    }
} // end GA_Recomb::renumbering_atom_numbering 

            //////////////////////////////////////////
            //****DEBUGGING FUNCTIONS *****//
            /////////////////////////////////////////



//++++++++++++++++++++++++++++++++++++++++++
//Check functions
void
GA_Recomb::is_active_vec(std::vector <DOCKMol> & mol_vec)
{
   Trace trace( "GA_Recomb::is_active_vec()" );
    for (int i=0; i<mol_vec.size(); i++){
        for (int j=0; j<mol_vec[i].num_bonds; j++){
            if(mol_vec[i].bond_active_flags[j]==true){
               cout <<"Active bonds (#) " << j <<endl;
            }

            if(mol_vec[i].bond_keep_flags[j]==true){
              cout <<"Keep bonds (#) " << j <<endl;
            }
        }
    }

    for (int i=0; i<mol_vec.size(); i++){
        for (int j=0; j<mol_vec[i].num_atoms; j++){
            if(mol_vec[i].atom_active_flags[j]==true){
               cout <<"Active atoms (#) " << j <<endl;
            }
        }
    }

}// end GA_Recomb::is_active()



//++++++++++++++++++++++++++++++++++++++++++
//Check functions

void
GA_Recomb::is_active(DOCKMol & mol)
{
   Trace trace( "GA_Recomb::is_active()" );
        for (int j=0; j<mol.num_bonds; j++){
            if(mol.bond_active_flags[j]==true){
               cout <<"Active bonds (#) " << j <<endl;
            }

            if(mol.bond_keep_flags[j]==true){
              cout <<"Keep bonds (#) " << j <<endl;
            }
        }

        for (int j=0; j<mol.num_atoms; j++){
            if(mol.atom_active_flags[j]==true){
               cout <<"Active atoms (#) " << j <<endl;
            }
        }

}// end GA_Recomb::is_active()



// +++++++++++++++++++++++++++++++++++++++++
// Debugging function to print the TorEnv_GA data structure
void
GA_Recomb::print_torenv( std::vector <Tor_Env>  tmp_torenv_vector )
{
   Trace trace( "GA_Recomb::print_torenv()" );
   int total = 0;
   cout <<endl;

    for (int i=0; i<tmp_torenv_vector.size(); i++){
        cout <<tmp_torenv_vector[i].origin_env <<":\n ";

        for (int j=0; j<tmp_torenv_vector[i].target_envs.size(); j++){

            cout <<"\t\t" <<tmp_torenv_vector[i].target_envs[j] <<"\t" <<tmp_torenv_vector[i].target_freqs[j] <<"\n";
            //cout <<tmp_torenv_vector[i].origin_env <<"-" <<tmp_torenv_vector[i].target_envs[j] <<"\n";

        }

        cout <<"\t\tTotal: " <<tmp_torenv_vector[i].target_envs.size() <<"\n";
        total = total+tmp_torenv_vector[i].target_envs.size();
    }

    cout <<"Final Total = " <<total <<"\n";

    return;

} // end DN_GA_Build::print_torenv()



// +++++++++++++++++++++++++++++++++++++++++++++++
double
GA_Recomb::time_seconds()
{
   Trace trace( "GA_Recomb::time_seconds()" );
   time_t          t;
   if (static_cast < time_t > (-1) == time(&t)) {
       cout << "Error from time function!  Elapsed time is erroneous." << endl;
   }
   return static_cast < double >(t);
}

            //////////////////////////////////////////
            //****OLD/ARCHIVE FUNCTIONS*****//
            /////////////////////////////////////////




// +++++++++++++++++++++++++++++++++++++++++
// Prepare and score parents - the program will exit if no scoring function was used
void
GA_Recomb::score_parents( std::vector <DOCKMol> & tmp, Master_Score & score,  AMBER_TYPER & typer)
{
   for (int i=0; i<tmp.size(); i++){
       activate_mol(tmp[i]);
       //cout << "Entering Amber typer" << endl;
       tmp[i].prepare_molecule();
       typer.skip_verbose_flag = true;
       typer.prepare_molecule(tmp[i], true, score.use_chem, score.use_ph4, score.use_volume);

       // Declare a temporary fingerprint object and compute atom environments (necessary for Gasteiger)
       Fingerprint temp_finger;
       for ( int j=0; j<tmp[i].num_atoms; j++ ){
           tmp[i].atom_envs[j] = temp_finger.return_environment( tmp[i], j );
       }
       float total_charges = 0.00;
       total_charges = compute_gast_charges( tmp[i] );


       /*if (score.use_score) {
          if (score.use_primary_score){
             score.compute_primary_score(tmp[i]);

             // compute internal energy for pruning
             tmp[i].internal_energy = score.primary_score->compute_ligand_internal_energy(tmp[i]);
             cout << "#### Parent " << i << "  Rescored to Fitness: " << tmp[i].current_score <<endl;
          }
      }*/
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
             score.primary_score->initialize_internal_energy(tmp[i]);
          }
          score.compute_primary_score(tmp[i]);
          cout << "Parent " << i << "  Rescored to Fitness: " << tmp[i].current_score <<endl;
      }
      else{
          cout << "No scoring function was found! Please re-run the GA experiment with a scoring function." <<endl;
          exit(0);
      }
   }
   return;
} //end GA_Recomb::score_parents()




// +++++++++++++++++++++++++++++++++++++++++++
// Eliminate redundate parents - within one vector
void
GA_Recomb::unique_parents(std::vector <DOCKMol> & tmp, Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & typer)
{

    Hungarian_RMSD h;
    pair <double, int> result;
    vector <DOCKMol> p;

    // If the two mols do not have the same number of atoms, rmsd is -1000
    //double torsion_rmsd_cutoff = 2.0;//MARK: user-specified?
    
    // These children need to be prepared with amber typer before rmsd calculation
    for ( int i=0; i<tmp.size(); i++ ){
        //typer.prepare_molecule(pruned_children[i], true, false, false);
        tmp[i].used = false;
    }

    // For each molecule in pruned_children
    for ( int i=0; i<tmp.size(); i++ ){
        // If it has not been marked used
        if ( !tmp[i].used ){
           //cout << " child " << i << " num atoms " << pruned_children[i].num_atoms << endl;
           //cout << " atom types i " << pruned_children[i].atom_types[0] << " " << pruned_children[i].atom_types[1] << endl;
           // For each molecule listed below it
           for ( int j=i+1; j<tmp.size(); j++ ){
               // If the below molecule has not been marked used
               if ( !tmp[j].used ){
                  //cout << " child " << j << " num atoms " << pruned_children[j].num_atoms << endl;
                  //cout << " atom types j " << pruned_children[j].atom_types[0] << " " << pruned_children[j].atom_types[1] << endl;
                  // Calculate and save hunarian rmsd
                  result = h.calc_Hungarian_RMSD_dissimilar( tmp[i], tmp[j] );
                  //rmsd = h.calc_Hungarian_RMSD( pruned_children[i], pruned_children[j] );
                  /*cout << "calc H rmsd = " << result.first <<endl;
                  cout << "num unmatched = " << result.second << endl;
                  cout <<"for child#= " << i << " with score " <<pruned_children[i].current_score <<endl;
                  cout <<"for child#= " << j << " with score " <<pruned_children[j].current_score <<endl;*/
                  // If the two confomers are very similar, mark the one with higher energy as used
                  // If the two molecules do not have the same number of atoms than rmsd is -1000
                  //if ( rmsd < torsion_rmsd_cutoff & rmsd != -1000 ){
                  if (result.first < ga_heur_matched_rmsd && result.second < ga_heur_unmatched_num) {
                      tmp[j].used = true;
                  }
               }
            }
         }
      }
    // At this point, only cluster heads will not be flagged 'used', add those to scored_generation
    // MARK - may want to sort by score before pruning also note that this rmsd does not take into consideration
    // the position of the molecule in the binding site so energy rank could be better
    for (int i=0; i<tmp.size(); i++){
        if (!tmp[i].used){
           p.push_back(tmp[i]);
        }
    }

    tmp.clear();
    for (int i=0; i<p.size(); i++){
        tmp.push_back(p[i]);
    }

    p.clear();
    

    return;
}//end GA_Recomb::unique_parents()



// ++++++++++++++++++++++++++++++++++++++++++++
// Create a unique name with the parents, generation number, and location in file
// If mol from all new, name is dn_gennum_loc
void
GA_Recomb::naming_function( DOCKMol & mol, int gen, int loc )
{

   //cout << "inside NAMING FUNCTION" <<endl;
   // Place holder for Name list
   vector <string> list;
   
   list.push_back(mol.subst_names[0]);

   // Search through the substrate names
   for ( int i=1; i<mol.num_atoms; i++ ){
       // If the next substrate name is equal to the previous
       //cout <<" name " << mol.subst_names[i] <<endl;
       if ( (mol.subst_names[i] != mol.subst_names[i-1])){
          // Check to see if it is already in the list
          size_t found = mol.subst_names[i].find_first_of("p");
          //cout << "NAMING FUNCTION " <<found <<endl;
 
          // If the atoms are from a parent and not in the list already
          if ( found == 0 && ( find ( list.begin(), list.end(), mol.subst_names[i] ) == list.end() )){
             // Add the parent name to the list
             list.push_back(mol.subst_names[i]);
          }
       }
   }


   // Add the items in list to the subst name
   ostringstream new_title;
   
   // Add gen first  
   if ( gen < 10 ){
       new_title << ga_name_identifier << "_g000" << gen ;
   }
   else if ( gen < 100 ){
       new_title << ga_name_identifier << "_g00" << gen ;
   }
   else if ( gen < 1000 ){
       new_title << ga_name_identifier << "_g0" << gen ;
   }
   else {
       new_title << ga_name_identifier << "_g" << gen ;
   }

   //cout << "new_title: " << new_title << endl;
   //mol.title = new_title.str();
   
   // Add the items in list to new_energy as a string
   /*for (int i=0; i<list.size(); i++){

       //Check to make sure that the first in the list is not <0>
       if (i == 0){
          size_t found = list[i].find_first_of("p");
          // If the atoms are from a parent and not in the list already
          if ( found == 0 ){
             new_title << list[i];
          }
       }
       else{
           new_title << "_" << list[i];
       }
   }


   // If list is empty, all new molecule
   if (list.size() == 0){
      new_title << "dn";  
   }*/
   //Add location
   if ( loc < 10 ){
       new_title << "_i000" << loc;
   }
   else if ( loc < 100 ){
       new_title << "_i00" << loc;
   }
   else if ( loc < 1000 ){
       new_title << "_i" << loc;
   }
   else {
       new_title << "_i" << loc;
   }

   //new_title << "_r" << loc;
   /* Add generation and location
   new_title << "_g" << gen << "_" << loc;
   //cout << "new_title: " << new_title << endl;*/
   mol.title = new_title.str();

   // Clear new_energy
   new_title.str("");
   list.clear();

   return;
}//end GA_Recomb::naming_function()


//JDB - calculate pairwise distances between all atoms in a molecule
void GA_Recomb::calc_pairwise_distance(DOCKMol & a)
{

    //Declare RMSD and total number of heavy atoms
    float distance = 0.0;
    int heavy_total = 0;
    bool lauren_var = false;

    //iterate through every atom
    for (int i=0; i<a.num_atoms; i++){
        for (int j=i+1; j<a.num_atoms; j++){
            distance = 0.0;
            //calculates the distance between two atoms           
            distance += sqrt((a.x[i]-a.x[j]) * (a.x[i]-a.x[j]) +
                        (a.y[i]-a.y[j]) * (a.y[i]-a.y[j]) +
                        (a.z[i]-a.z[j]) * (a.z[i]-a.z[j]));
            //if the distance is below a hardcoded cutoff, output an error
            if (distance < 0.5){
                cout << "GA_DEBUG::SMALL BOND SMALL BOND SMALL BOND" << endl;
                cout << "NAME: " << a.title << endl;
                //some_molecules_go_here.push_back(a);
                lauren_var = true;
 
            }
        }
    }
    //if (lauren_var) { some_molecules_go_here.push_back(a);}

}
// +++++++++++++++++++++++++++++++++++++++++++++++
// Calculate tanimoto pairwise LEP
void GA_Recomb::pairwise_tanimoto_calc(vector <DOCKMol> & mol)
{
    //prep some files
    ostringstream fout_txt_name;
    fout_txt_name << ga_output_prefix << "_tan_pairwise_matrix.csv";
    fstream fout_crap;
    fout_crap.open(fout_txt_name.str().c_str(), fstream::out|fstream::app );
    
    // for clustering TODO
    /*ostringstream fout_clus_txt_name;
    fout_clus_txt_name << ga_output_prefix << "_best_first_cluster_tan.csv";
    fstream fout_crap1;
    fout_crap1.open(fout_clus_txt_name.str().c_str(), fstream::out|fstream::app );
    */

    //declare fingerprint obj
    Fingerprint finger;
    float tanimoto;

    for (int i=0; i<mol.size(); i++){
        mol[i].used = false;
    }

    //initialize some var/vects
    vector<float> list_of_tanimotos;
    INTVec tmp_list_of_indeces;
    INTVec list_of_indeces;
    INTVec tmp_cluster;
    INTVec list_of_mol;
    int clust_counter=0;

    for (int i=0; i<mol.size(); i++){
        list_of_mol.push_back(i);
        for (int j=i+1; j<mol.size(); j++){
            tanimoto = finger.compute_tanimoto(mol[i], mol[j]);

	    fout_crap << tanimoto << ",";

            //if (i ==1){
            //    list_of_tanimotos.push_back(tanimoto);
	    //	list_of_mol.push_back(j);
            //}
            if (tanimoto == 1){
                cout << "tanimoto = 1" << endl;
		mol[j].used = true;
            }
        }
   }
    // if tan sim prune turned on push only those without tan=1
    if (ga_tan_similarity_prune){
        for (int i=0; i<mol.size(); i++){
            if( mol[i].used == false ){
                prunemol.push_back(mol[i]);
            }
        }
    }
    // start of clustering by tanimoto TODO
    /*if (ga_tan_best_first_clustering){
        if (list_of_mol.size()> 1){
        for (int i=0; i<list_of_mol.size(); i ++){
            if (list_of_tanimotos[i] >= 0.9){
                tmp_cluster.push_back(list_of_mol[i]);
            }else{
                tmp_list_of_indeces.push_back(list_of_mol[i]);
            }
        }
        for  (a=0, a<tmp_cluster.size();a++){
            fout_crap1 << mol[tmp_cluster[a]].title << "," << clust_counter << endl;
        }

        clust_counter++;
        tmp_cluster.clear();
        list_of_tanimotos.clear();
        list_of_mol.clear();

        for (int i=0; i<tmp_list_of_indeces.size(); i ++){

            list_of_mol.push_back(tmp_list_of_indeces[i]);
            for (int j=i+1; j<tmp_list_of_indeces.size(); j++){
                tanimoto = finger.compute_tanimoto(mol[tmp_list_of_indeces[i]], mol[tmp_list_of_indeces[j]]);
                list_of_tanimotos.push_back(tanimoto);
            }
        }
        }//if size greater than 1
    }else{
        list_of_tanimotos.clear();
        list_of_mol.clear();
    }
    if (tmp_cluster.size()>0){

    }*/
    fout_crap.close();
    //fout_crap1.close();
    fout_txt_name.clear();

    return;
}
// +++++++++++++++++++++++++++++++++++++++++++++++
// calculate HMS pairwise LEP
void GA_Recomb::pairwise_hms_calc(vector <DOCKMol> & mol)
{
    //prep the outfiles
    ostringstream fout_txt_name;
    fout_txt_name << ga_output_prefix << "_hms_pairwise_matrix.txt";
    fstream fout_crap;
    fout_crap.open(fout_txt_name.str().c_str(), fstream::out|fstream::app );

    //declare some stuff
    Hungarian_RMSD hms;
    std::pair <double, int>  hms_num;
    int hun_ref_heavy_atoms;
    float temp_hun_score;

    //loop through mols
    for (int i=0; i<mol.size(); i++){

	//init var at 0
        hun_ref_heavy_atoms=0;

	//loop through all atoms to count heavy atoms and increment var
        for (int k=0; k<mol[i].num_atoms; k++){
            if (mol[i].atom_types[k].compare("H") != 0 && mol[i].atom_types[k].compare("Du") != 0){
                hun_ref_heavy_atoms++;
            }
        }
	
	// loop through and calc rmsd of matching atoms
        for (int j=i+1; j<mol.size(); j++){
            hms_num = hms.calc_Hungarian_RMSD_dissimilar(mol[j], mol[i]);
            temp_hun_score = (-5) * (hun_ref_heavy_atoms - hms_num.second) / hun_ref_heavy_atoms + hms_num.first;

            //spit out pairiwse hms numbers	
            fout_crap << temp_hun_score << ",";
        }
    }
    
    //close and clear
    fout_crap.close();
    fout_txt_name.clear();

    return;
}


// +++++++++++++++++++++++++++++++++++++++++++++++
// Populate donotmutate bools LEP
void
GA_Recomb::set_DNM_bools( DOCKMol & mol )
{
    for (int i = 0; i < mol.num_atoms; i++){
        if (mol.atom_names[i] == "DNM"){
            mol.atom_dnm_flag[i] = true;
            if (!mol.mol_dnm_flag){
                mol.mol_dnm_flag=true;
            }
        }
        else{
            mol.atom_dnm_flag[i] = false;
        }
    }
    return;
}
// +++++++++++++++++++++++++++++++++++++++++++++++
//Check through the whole parent list for DNM - sets the boolean for the run
//and each atom's boolean.
void
GA_Recomb::set_DNM_bools_start( std::vector <DOCKMol> & parents)
{
    dnm_enabled=false;
    //loop through all parents
    for (int par_id=0; par_id<parents.size(); par_id++){
        //all atoms in parent
        for (int atom_ind=0; atom_ind<parents[par_id].num_atoms; atom_ind++) {
            if ( parents[par_id].atom_names[atom_ind] == "DNM"){ //DNM Case
                parents[par_id].atom_dnm_flag[atom_ind] = true;
                // sets DNM status for run to true on first instance
                if (!dnm_enabled){ 
                    dnm_enabled=true;
                }
                // sets DNM status for this molecule to true on first instance
                if (!parents[par_id].mol_dnm_flag){ 
                    parents[par_id].mol_dnm_flag = true;
                }
            // All other cases, atom isn't DNM
            } else { 
                parents[par_id].atom_dnm_flag[atom_ind] = false;
            }
        }
    }

}
// +++++++++++++++++++++++++++++++++++++++++++++++


