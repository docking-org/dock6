#include "conf_gen_dn_ga.h"
#include "conf_gen_ga.h"
#include "fingerprint.h"
#include "hungarian.h"
#include "gasteiger.h"
#include <iomanip>
#include "trace.h"
//#include <cstdlib>
#include <list>
#include <time.h>
using namespace std;

class Bump_Filter;
class Fingerprint;
class AMBER_TYPER;
class Master_Score;

// static member initializers
// These parameters are consistent with those in Base_Score()
const string DN_GA_Build::DELIMITER    = "########## ";
const int    DN_GA_Build::FLOAT_WIDTH  = 20;
const int    DN_GA_Build::STRING_WIDTH = 17 + 19;


//
// +++++++++++++++++++++++++++++++++++++++++
// Some constructors and destructors

TorEnv_GA::TorEnv_GA(){
    origin_env = "";
}
TorEnv_GA::~TorEnv_GA(){
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

DN_GA_Build::DN_GA_Build(){
    // Flag used to save & print mutated molecules for GA - CDS
    dn_ga_flag = false;
    dn_ga_gen = 0;
    dn_ga_mut_attempt = 0;
    tmp_first_check = false;
    tmp_second_check = false; //temporary fix for uninitialized value, set to true to be conservative, takes extra time but checks both direction sof the torenv (?)
    iso_picks=0;
 dn_ga_verbose = 0;
}
DN_GA_Build::~DN_GA_Build(){
    scaffolds.clear();
    linkers.clear();
    sidechains.clear();
    anchors.clear();
    scaf_link_sid.clear();
    dn_ga_verbose = 0;
}



// +++++++++++++++++++++++++++++++++++++++++
// Read parameters from the dock.in file (called in master_conf.cpp)
void
DN_GA_Build::input_parameters( Parameter_Reader & parm )
{
    cout <<endl <<"De Novo Build Parameters" <<endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

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

    // Ask the user for a torsion environment table of allowable environments
    // (The only time you might NOT want this is if you are doing simple_build. Otherwise it really
    // is useful in preventing some garbage from being created.)
    //dn_use_torenv_table = (parm.query_param("dn_use_torenv_table", "yes", "yes no") == "yes");
    // LEP 2019.03.13 hardcode torsion table on 
    dn_use_torenv_table = true;
    if (dn_use_torenv_table){
        dn_torenv_table = parm.query_param("dn_torenv_table", "fraglib_torenv.dat");
        dn_use_roulette = parm.query_param("dn_use_roulette", "yes", "yes no") == "yes";
    }
    
    denovo_name = parm.query_param( "dn_name_identifier",
                                                  "denovo" );

    // Specify a specific dn sampling method
    // ex = exhaustive (slow)
    // rand = random (user specified limit)
    // graph = sample from graph of related fragments
    dn_sampling_method = parm.query_param("dn_sampling_method", "graph", "ex | rand | graph");
    dn_sampling_method_ex = false;
    dn_sampling_method_rand = false;
    dn_sampling_method_graph = false;

    if (dn_sampling_method.compare("ex") == 0){
        dn_sampling_method_ex = true;
    } else if (dn_sampling_method.compare("rand") == 0){
        dn_sampling_method_rand = true;
    } else if (dn_sampling_method.compare("graph") == 0){
        dn_sampling_method_graph = true;
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

    // User specifies a hard upper score cutoff, similar to pruning_conformer_score_cutoff in 
    // anchor-and-grow
    dn_pruning_conformer_score_cutoff = atof(parm.query_param("dn_pruning_conformer_score_cutoff",
                                                              "100.0").c_str());
    dn_pruning_conformer_score_cutoff_begin = dn_pruning_conformer_score_cutoff;
    dn_pruning_conformer_score_scaling_factor = atof(parm.query_param("dn_pruning_conformer_score_scaling_factor", "1.0").c_str());

    // User specifies a cutoff for pruning torsions - equivalent to the anchor and grow heuristic
    dn_pruning_clustering_cutoff = atof(parm.query_param("dn_pruning_clustering_cutoff",
                                                         "100.0").c_str());
    // User can specify an upper bound for molecular weight 
    dn_constraint_mol_wt = atof(parm.query_param("dn_constraint_mol_wt", "550.0").c_str());

    // User can specify whether they want a hard or soft MW cutoff - SMT
    dn_MW_cutoff_type = parm.query_param("dn_mol_wt_cutoff_type", "soft", "hard | soft");
    dn_MW_cutoff_type_hard = false;
    dn_MW_cutoff_type_soft = false;
  
    if (dn_MW_cutoff_type.compare("hard") == 0) {
        dn_MW_cutoff_type_hard = true;
    } else if (dn_MW_cutoff_type.compare("soft") == 0) {
        dn_MW_cutoff_type_soft = true;
    } else {
        cout <<"You chose...poorly." <<endl;
        exit(0);
    }

    //User can specify the std deviation of the MW distribtuion - SMT
    if (dn_MW_cutoff_type_soft){
        dn_MW_std_dev = atof(parm.query_param("dn_mol_wt_std_dev", "50.0").c_str());
    }
 
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

    // Write checkpoint files at each layer, useful to monitor progress and/or restart jobs
    dn_write_checkpoints = (parm.query_param("dn_write_checkpoints", "yes", "yes no") == "yes");

    // Write the molecules that are pruned at each layer, useful to see what we are throwing away
    dn_write_prune_dump = (parm.query_param("dn_write_prune_dump", "no", "yes no") == "yes");

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

} // end DN_GA_Build::input_parameters()



// +++++++++++++++++++++++++++++++++++++++++
// The initialize function tells the user the filenames it is reading from, and it also populates 
// the fragment library vectors. (called in master_conf.cpp)
void
DN_GA_Build::initialize()
{
    Trace trace("DN_GA_Build::initialize()");
    cout <<endl <<"Initializing De Novo Growth Routines..." <<endl;

    cout <<" Reading the scaffold library from " <<dn_fraglib_scaffold_file <<"...";
    read_library( scaffolds, dn_fraglib_scaffold_file );
    cout <<"Done (#=" <<scaffolds.size() <<")" <<endl;

    cout <<" Reading the linker library from " <<dn_fraglib_linker_file <<"...";
    read_library( linkers, dn_fraglib_linker_file );
    cout <<"Done (#=" <<linkers.size() <<")" <<endl;
       
    cout <<" Reading the sidechain library from " <<dn_fraglib_sidechain_file <<"...";
    read_library( sidechains, dn_fraglib_sidechain_file );
    cout <<"Done (#=" <<sidechains.size() <<")" <<endl;

    // We don't need to include rigids in de novo growth, but this can be uncommented for testing
    /*
    cout <<" Reading the rigid library from " <<dn_fraglib_rigid_file <<"...";
    read_library( rigids, dn_fraglib_rigid_file );
    cout <<"Done (#=" <<rigids.size() <<")" <<endl;
    */

    if (dn_user_specified_anchor){
        cout <<" Reading the anchor library from " <<dn_fraglib_anchor_file <<"...";
        read_library( anchors, dn_fraglib_anchor_file );
        cout <<"Done (#=" <<anchors.size() <<")" <<endl;
    }
    if (dn_use_torenv_table){
        cout <<" Reading the torenv table from " <<dn_torenv_table <<"...";
        read_torenv_table(dn_torenv_table);
    }
    if (dn_ga_flag==false){
        if (dn_use_roulette)
        {
            cout << " Generating the roulette wheel from the torenv table ... ";
            generate_roulette();
        }
    }
    // We used to do this for just graph, but now we do it for all sampling functions for simplicity
    if (dn_sampling_method_graph || dn_sampling_method_rand || dn_sampling_method_ex){

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

    // This function will find what fragments are similar to what other fragments in each of the
    // libraries, storing that information as a FragGraph object
    if (dn_sampling_method_graph){

        prepare_fragment_graph( scaf_link_sid, scaf_link_sid_graph );

        // To visualize or inspect graph, uncomment this. Note: exits program on finishing.
        //print_fraggraph();
    }
    return;

} // end DN_GA_Build::initialize()



// +++++++++++++++++++++++++++++++++++++++++
// Get the internal energy parameters from Master_Conformer_Search
void
DN_GA_Build::initialize_internal_energy_parms( bool uie, int rep_exp, int att_exp, float diel, float iec )
{
    Trace trace("DN_GA_Build::initialize_internal_energy_parms()");
    // These just need to be defined before we try to call the minimizer

    use_internal_energy = uie;         // boolean, use internal energy in ligand minimization?
    ie_rep_exp          = rep_exp;     // int. energy vdw repulsive exponent (default 12)
    ie_att_exp          = att_exp;     // int. energy vdw attractive exponent (6)
    ie_diel             = diel;        // int. energy dielectric (4.0)
  //internal_energy_cutoff = iec;
    ie_cutoff           = iec;         // BCF internal energy cutoff 
    return;

} // end DN_GA_Build::initialize_internal_energy_parms()



// +++++++++++++++++++++++++++++++++++++++++
// Read fragment libraries from file and save in appropriate vector of <Fragments>
void
DN_GA_Build::read_library( vector <Fragment> & frag_vec, string filename )
{
    // Create temporary DOCKMol object
    DOCKMol tmp_frag_mol;

    // Create filehandle for reading mol2 library
    ifstream fin_frags;

    // Open the filehandle for reading a fragment library
    fin_frags.open(filename.c_str());
    if (fin_frags.fail()){
        cout <<"Could not open " <<filename <<". Program will terminate." <<endl;
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

} // end DN_GA_Build::read_library()



// +++++++++++++++++++++++++++++++++++++++++
// This function reads the table of allowable environments and stores them in a data structure 
void
DN_GA_Build::read_torenv_table( string torenv_table ) 
{
    // Open input filestream given the filename from parameter file 
    ifstream fin_torenv_table;
    fin_torenv_table.open(torenv_table.c_str());

    if (fin_torenv_table.fail()){
        cout <<"Could not open " <<torenv_table <<". Program will terminate." <<endl;
        exit(0);
    }

    // Declare some temporary variables
    TorEnv_GA temp_torenv;
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

                // Add it to the vector of TorEnv_GA data structures
                counter += temp_torenv.target_envs.size();
                torenv_vector.push_back(temp_torenv);
                temp_torenv.origin_env.clear();
                temp_torenv.target_envs.clear();
            }

            // Then begin populating a new TorEnv_GA
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

} // end DN_GA_Build::read_torenv_table()

//JDB
//Comparator functions for vector of vector sorting
//order_torsions sorts least to greatest by the second position of each vector
bool ga_roulette_sort_order_torsions(const vector<double>& vect1, const vector<double>& vect2)
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
bool ga_roulette_sort_order_frequencies(const vector<double>& vect1, const vector<double>& vect2)
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


//JDB - function to generate a roulette wheel from the provided torsion table
void
DN_GA_Build::generate_roulette () {
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
    sort(freq_ordering.begin(), freq_ordering.end(), ga_roulette_sort_order_frequencies);
    
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
            duplicate_frequency_collector = double(freq_percent[i-1]) + (double(duplicate_pos_collector.size()) * freq_percent[i]);
            
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
    sort(freq_ordering.begin(), freq_ordering.end(), ga_roulette_sort_order_torsions);
    
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
} // end DN_GA_Build::generate_roulette()


// +++++++++++++++++++++++++++++++++++++++++
// Iterate through a vector of fragments, compute pair-wise Tanimotos, and remember rank ordered
// lists of similarity using the FragGraph data structure.
void
DN_GA_Build::prepare_fragment_graph( vector <Fragment> & frag_vec, vector <FragGraph> & fraggraph_vec )
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
        sort(temp_fraggraph.tanvec.begin(), temp_fraggraph.tanvec.end(), ga_fgpair_sort);

        // Record the ranks in the rank vector
        // (start at 1 to avoid recording the self-self Tanimoto)
        for (int j=1; j<temp_fraggraph.tanvec.size(); j++){
            temp_fraggraph.rankvec.push_back( temp_fraggraph.tanvec[j].second );
        }

        // Add this object to the vector of FragGraph
        fraggraph_vec.push_back(temp_fraggraph);
    }

    return;

} // end DN_GA_Build::prepare_fragment_graph()



// +++++++++++++++++++++++++++++++++++++++++
// This function is the main de novo engine for building molecules. (called in dock.cpp)
void
DN_GA_Build::build_molecules( Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & typer,
                           Orient & orient )
{   
    Trace trace("DN_GA_Build::build_molecules()");
    if (dn_ga_flag == false){
       cout <<endl <<endl <<"##### Entering the main de novo build engine #####" <<endl <<endl;
    }
    molecule_counter = 0;    // this is an estimate of how many molecules are sampled
    growth_tree_index = 0;   // in case we are writing growth trees

    if (verbose){ cout <<"verbose is turned on" <<endl; }

    if (orient.orient_ligand == false && dn_user_specified_anchor == false){
        cout <<endl
             <<"# Note: No anchors were specified AND orient is turned off." <<endl
             <<"# DOCK now assumes that all scaffolds, linkers, and sidechains" <<endl
             <<"# provided are already oriented to the binding site." <<endl <<endl;
    }

    // fout_molecules is the filestream for complete molecules
    ostringstream fout_molecules_name;
    fout_molecules_name <<dn_output_prefix <<".denovo_build.mol2";
    fstream fout_molecules;
    //LEP - 2018.06.20 
    if (dn_ga_flag == false){
        fout_molecules.open (fout_molecules_name.str().c_str(), fstream::out|fstream::app);
    }
    fout_molecules_name.clear();
    

    // Declare some fragment vectors for growing
    vector <Fragment> root;          // starting fragments, returned each layer
    vector <Fragment> layer;         // vector for completing a layer
    vector <Fragment> next_layer;    // to be moved to the next layer
    vector <Fragment> growing;       // temporary vector for torsions / minimizing
    vector <Fragment> prune_dump;    // not as bad as it sounds


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
        typer.skip_verbose_flag = true;
        typer.prepare_molecule(anchors[i].mol, true, score.use_chem, score.use_ph4, score.use_volume);
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
    //sort(anchors.begin(), anchors.end(), ga_size_sort);
    frag_sort(anchors, ga_size_sort);

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
    anchors.erase (anchors.begin()+dn_unique_anchors, anchors.end());


    ////////////////////////////////////////////////////////////////////////////
    // Begin Main Loop                                                        //
    ////////////////////////////////////////////////////////////////////////////

    //simplex.initialize();
    // Initiate a new round of de novo growth for each anchor
    for (int a=0; a<anchors.size(); a++){

        // Reset some values here that may have been scaled for the previous anchor
        if (dn_sampling_method_graph){
           dn_graph_temperature = dn_temp_begin;
        }

        dn_pruning_conformer_score_cutoff = dn_pruning_conformer_score_cutoff_begin;
        
        // Reset the random seed so that it is the same for every anchor
        // This behavior is quasi-consistent with A&G
        //simplex.initialize();

       // Clear the non-bonded pairlist() from the previous anchor
        score.primary_score->nb_int.clear();
        trace.note("nb array clear");
    
        // If orient_ligand is set to yes
        if (orient.orient_ligand){
    
            // Orient one anchor and return the results on the vector root
            cout <<"##### Now orienting anchor #" <<(a+1) <<endl;

           // if (verbose){
                // This outputs what anchor is being chosen CPC
            //    cout<<"##### Name of anchor chosen is " <<anchors[a].mol.energy<<endl;
            //}    

            orient_fragments(root, anchors[a], score, simplex, typer, orient);
            //cout << "Orient Anchor Charge Debug" <<  
        } else {

            // Just put one anchor onto root - it should already be in the binding site
            if (dn_ga_flag == false){
               cout <<"##### Now growing from anchor #" <<(a+1) <<endl;
            }

            root.push_back(anchors[a]);
        }
    
        // counter for controlling the next while loop
        int counter = 0;
    
        // last_layer feature to help control sampling
        bool last_layer = false;
    
        // while there are still molecules in root && for as many layers as you want to grow out
        // (each iteration of this while loop is equivalent to growing out one layer)
        while (counter < dn_max_grow_layers){
    
            if (counter == (dn_max_grow_layers-1)) {last_layer = true;}
            if (dn_ga_flag == false){
               if (root.size() == 0){ cout <<"Root is empty, continuing." <<endl <<endl; break;} 
               cout <<"  ##### Entering layer of growth #" <<(counter+1) <<endl;
               cout <<"  Root size at the beginning of this layer = " <<root.size() <<endl;
            }
 
            // for each fragment in root
            for (int i=0; i<root.size(); i++){
    
                // In case the fragment has no attachment points and it is the beginning of growth,
                // prepare it now just so there will actually be a score
                if (counter == 0 && root[i].aps.size() == 0){ 
                    activate_mol(root[i].mol);
                    root[i].mol.prepare_molecule();

                    typer.skip_verbose_flag = true;
                    typer.prepare_molecule(root[i].mol, true, score.use_chem, score.use_ph4, score.use_volume);
                    score.compute_primary_score(root[i].mol);
                    root[i].scaffolds_this_layer = 0;
                }
    
                // Set scaffolds this layer to 0 before we start growth
                root[i].scaffolds_this_layer = 0;
               
                // Copy the current frag from root into layer
                layer.push_back(root[i]);
    
                // Count the number of attachment points in layer[0]
                int temp = layer[0].aps.size();
    
                // For each of those attachment points
                // (do the last attachment point first, then increment down. this will have the effect
                // of filling every attachment point on the fragment before going to the next layer)
                for (int j=temp-1; j>-1; j--){
    
                    // For all of the fragments in layer
                    for (int k=0; k<layer.size(); k++){
                        // Copy on the old dockmol for the growth tree
                        if (dn_write_growth_trees){
                            layer[k].frag_growth_tree.push_back(layer[k].mol);
                        }
//TODO
// in some of these different sampling methods, add additional print statements that are activated by a verbose flag on the command line
                        // remember how many things are on growing before sampling
                        int growing_size_before = growing.size();
    
                        // Sampling method = exhaustive
                        if (dn_sampling_method_ex){
                            sample_fraglib_exhaustive( layer[k], j, scaf_link_sid, growing, 
                                                       score, simplex, typer, last_layer );
    
                        // Sampling method = random
                        } else if (dn_sampling_method_rand){
                            sample_fraglib_rand( layer[k], j, scaf_link_sid, growing,
                                                 score, simplex, typer, last_layer );
    
                        // Sampling method = graph
                        } else if (dn_sampling_method_graph){
                            sample_fraglib_graph( layer[k], j, scaf_link_sid, scaf_link_sid_graph,
                                                  growing, score, simplex, typer, last_layer );
    
                        } else {
                            cout <<"You chose...poorly" <<endl;
                            exit(0);
                        }

//TODO
// dummy_to_H(layer[k], layer[k].aps[j].heavy_atom, layer[k].aps[j].dummy_atom);  growing.push_back(layer[k]);

                    } // end for each fragment in layer


                    // clear layer and get ready to copy growing onto layer
                    layer.clear();
    
                    // sort growing on the mol object primary score
                    cout << "Sorting molecules made after denovo" << endl;
                    //sort(growing.begin(), growing.end(), ga_fragment_sort);
		    frag_sort(growing, ga_fragment_sort);
                    molecule_counter += growing.size();
                cout << endl << "Molecules made after denovo: " << growing.size() << endl << endl;

                    /*ostringstream fout_prefix_name;
                    fout_prefix_name << "zzz.dnmols_after_randsampl_gen" << dn_ga_gen+1 << ".mol2";
                    fstream fout_molecules2;
                    fout_molecules2.open ( fout_prefix_name.str().c_str(), fstream::out|fstream::app );
                    for( int i=0; i< growing.size(); i++ ){
                        fout_molecules2 <<DELIMITER<<setw(STRING_WIDTH)<<"Frag_String:"<<setw(FLOAT_WIDTH)<<growing[i].mol.energy <<endl;
                        fout_molecules2 <<DELIMITER<<setw(STRING_WIDTH)<<"Generation:"<<setw(FLOAT_WIDTH)<< dn_ga_gen+1 <<endl;
                        fout_molecules2 <<DELIMITER<<setw(STRING_WIDTH)<<"Mutation_Attempt:"<<setw(FLOAT_WIDTH)<< dn_ga_mut_attempt <<endl;
                        fout_molecules2 << growing[i].mol.current_data << endl;
                        Write_Mol2( growing[i].mol, fout_molecules2 );
                    }
                    fout_molecules2.close();*/

//TODO these descriptors are now computed in sample_minimized torsions and filtering is done then
//     doing it again here may be redundant if that information is copied when you are creating
//     the torsions
                    // copy at most dn_max_layer_size back onto layer - at this point the loop will
                    // continue adding fragments to the remaining available attachment points to fill
                    // a layer
                    //cout << "growing size: " << growing.size() << endl;
                    for (int x=0; x<growing.size(); x++){
    
                        calc_mol_wt(growing[x].mol);
                        calc_rot_bonds(growing[x].mol);
                        calc_formal_charge(growing[x].mol);
                        //calc_num_HA_HD(growing[i].mol);

// write and implement calc_chiral_centers() function

                        if (dn_MW_cutoff_type_hard){
                        // Only copy things that are within the user-defined property constraints - SMT
                        	if (growing[x].mol.mol_wt <= dn_constraint_mol_wt  && 
                           	    growing[x].mol.rot_bonds <= dn_constraint_rot_bon  &&
                           	    fabs(growing[x].mol.formal_charge) <= (fabs(dn_constraint_formal_charge)+0.1) )
				            {   //cout << "Molecules passed MW, RB, and FC" << endl;

    
                            // And only copy until you reach user specified max_layer_size
                           	    if (layer.size() < dn_max_layer_size){
                               	        layer.push_back(growing[x]);
                           	    } else {
                              	        break;
                           	    }
    
                       	     }else{
                                 cout << "Did not pass MW, RB, and/or FC" << endl;
                             }
			            }
		        // Only copy things where MW probability beats random number
                        else if (dn_MW_cutoff_type_soft)
			{
                                if (growing[x].mol.mol_wt <= dn_constraint_mol_wt  &&
                                    growing[x].mol.rot_bonds <= dn_constraint_rot_bon  &&
                                    fabs(growing[x].mol.formal_charge) <= (fabs(dn_constraint_formal_charge)+0.1) ){
                                    
                            
                            // And only copy until you reach user specified max_layer_size     
                                    if (layer.size() < dn_max_layer_size)
				    {
                                        layer.push_back(growing[x]);
                                     }
           
                                } else if (growing[x].mol.mol_wt > dn_constraint_mol_wt){
                                    rando_num = ((double) rand() / (RAND_MAX));
                                    //cout <<"Max Random number that could be generated is: " <<(RAND_MAX) <<endl;
                                    if (Parameter_Reader::verbosity_level() > 1) {
				        cout <<"Random # between 0 and 1 is: " <<(rando_num) <<endl;
				    }
                                          
                                    //var_one_MW = 1 / (dn_MW_std_dev * sqrt(2 * M_PI));
                                    var_two_MW = growing[x].mol.mol_wt - dn_constraint_mol_wt;
                                    var_three_MW = (var_two_MW*var_two_MW);
                                    var_four_MW = -1 * var_three_MW;
                                    var_five_MW = (dn_MW_std_dev*dn_MW_std_dev);
                                    var_six_MW = 2 * var_five_MW;
                                    var_seven_MW = var_four_MW / var_six_MW;
                                    var_eight_MW = exp (var_seven_MW);
				    if (Parameter_Reader::verbosity_level() > 1) {
                                        cout <<"Molecular weight of molecule in question is :" 
                                             <<(growing[x].mol.mol_wt) <<endl;
                                        cout <<"Probability of accepting molecule based on Norm. Dist. is :" 
                                             <<(var_eight_MW) <<endl;
                                    }
                                    if (rando_num <= var_eight_MW &&
                                        growing[x].mol.rot_bonds <= dn_constraint_rot_bon &&
                                        fabs(growing[x].mol.formal_charge) <= (fabs(dn_constraint_formal_charge)+0.1)) {
                                        if (layer.size() < dn_max_layer_size){
                                            layer.push_back(growing[x]);
                                         }
                                    } else {
                                        break;
                                    }
                                }
                        }
                    }
		    //SMT - ends here
                    // clear growing to prepare for the next round of fragments
                    growing.clear();
    
                } // end for each attachment point in the fragment
                if (dn_ga_flag == false){ 
                   cout <<"    From Root[" <<i <<"], created " <<layer.size() 
                        <<" new fragments with an added layer" <<endl;
                }

                // put all the contents of layer onto next_layer - this clears up layer for the next
                // fragment in root, and later we will sort / prune all the contents of next_layer to
                // pick fragments to put back onto root
                for (int x=0; x<layer.size(); x++){
                    next_layer.push_back(layer[x]);
                }
                cout << "Number of molecules in next_layer: " << next_layer.size() << endl;
    
                // clear the contents of layer
                layer.clear();
    
            } // end for each fragment in root
    
//TODO
// add an option to GA through all the new fragments you just built at the appropriate place in the code below (no additions or deletions, just cross-overs with occasional mutations) 
    
            // at this point, we have built one complete layer - clear root and sort / prune the
            // contents of next_layer to figure out what frags will move onto the next layer of growth
            root.clear();
            cout << "Sorting denovo molecules based on score." << endl;
            // choose the sort function based on input parameters
            // sort(next_layer.begin(), next_layer.end(), ga_fragment_sort);
	    frag_sort(next_layer, ga_fragment_sort);
            if (dn_ga_flag == false){ 
               cout <<"      From all frags in root, kept a total of " <<next_layer.size() 
                    <<" new fragments with an added layer" <<endl;
            }

            // once next_layer is sorted, remove some redundancy based on the Hungarian RMSD heuristic
            cout << "Entering Hungarian RMSD heuristic pruning." << endl; 
            prune_h_rmsd(next_layer);
            if (dn_ga_flag == false){
               cout <<"      After pruning by the RMSD heuristic, that number goes down to " 
                    <<next_layer.size() <<endl;
            }
    
            // prepare to calculate some ensemble properties for printing to stdout
            float avg_mol_wt = 0.0;
            float avg_rot_bonds = 0.0;
            float avg_formal_charge = 0.0;
            float avg_num_HA = 0.0;
            float avg_num_HD = 0.0;
    
            if (dn_ga_flag == false){
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
            }
    
            // If for some reason the next_layer size was 0, you don't want to get an accidental
            // divide by 0 error here
            if (dn_ga_flag == false){
               if (next_layer.size() != 0){
                   avg_mol_wt /= next_layer.size();
                   cout <<"      From this ensemble, the average molecular weight is " <<avg_mol_wt <<endl; 
                   avg_rot_bonds /= next_layer.size();
                   cout <<"      From this ensemble, the average number of rotatable bonds is " <<avg_rot_bonds <<endl; 
                   avg_formal_charge /= next_layer.size();
                   cout <<"      From this ensemble, the average formal_charge is " <<avg_formal_charge <<endl; 
                  // avg_num_HA /= next_layer.size();
                  // cout <<"      From this ensemble, the average hbond acceptors is " <<avg_num_HA <<endl;
                  // avg_num_HD /= next_layer.size();
                  // cout <<"      From this ensemble, the average hbond donors is " <<avg_num_HA <<endl;

               }
            }

            // prepare to write some molecules to file
            int write_to_file_counter=0;
    
            // Go through next_layer one by one. If there are still dummy atoms, push it back onto root
            // (provided you don't exceed dn_max_root_size). If there are no dummy atoms, it is a 
            // complete molecule, so write it to file.
            
            // If using GA, save a subset of molecules - CDS
            if (dn_ga_flag && (next_layer.size() != 0)){
               // If at last layer, return all molecules (-1 since the counter has not been incremented
               if (counter == dn_max_grow_layers-1){
                  for ( int j=0; j<next_layer.size(); j++ ){
                     tmp_mutants.push_back(next_layer[j]);
                  }
               }
               // Otherwise keep the best score intermediate
               else{
                 tmp_mutants.push_back(next_layer[0]);
               }
               cout << "Number of molecules being sent to GA: " << tmp_mutants.size() << endl; 
            }


            for (int i=0; i<next_layer.size(); i++){
    
                // If there is a dummy remaining...
                if (dummy_in_mol(next_layer[i].mol)){

                     // ...put it back in root for the next layer
                     if (root.size() < dn_max_root_size){
                         root.push_back(next_layer[i]);
                     }

                     // also move up to NNNN onto prune_dump, cap if you can, and print to file
                     if (prune_dump.size() < 10000){
                         prune_dump.push_back(next_layer[i]);
                     }

                // Else if there is no dummy remaining...
                } else {
                    if (dn_ga_flag == false){
                       // LEP - save unique name to dockmol object title (instead of energy)  
                       ostringstream new_comment;
		       new_comment << denovo_name<< "_layer"<<counter+1<<"_"<< write_to_file_counter;
		       next_layer[i].mol.title = new_comment.str();
		       new_comment.clear();
		       // write mols to file
                       //fout_molecules <<next_layer[i].mol.current_data;
		       //LEP added in unique denovo name and scientific notation to floating point for output
		               fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Name:" << setw(FLOAT_WIDTH) << next_layer[i].mol.title <<endl;
                       fout_molecules <<DELIMITER<<setw(STRING_WIDTH)<<"Molecular_Weight:"<<setw(FLOAT_WIDTH)<<next_layer[i].mol.mol_wt <<endl;
                       fout_molecules <<DELIMITER<<setw(STRING_WIDTH)<<"DOCK_Rotatable_Bonds:"<<setw(FLOAT_WIDTH)<<next_layer[i].mol.rot_bonds <<endl;
                       fout_molecules << fixed <<DELIMITER<<setw(STRING_WIDTH)<<"Formal_Charge:"<<setw(FLOAT_WIDTH)<<next_layer[i].mol.formal_charge <<endl;
                       fout_molecules << fixed <<DELIMITER<<setw(STRING_WIDTH)<<"HBond_Acceptors:"<<setw(FLOAT_WIDTH)<<next_layer[i].mol.hb_acceptors <<endl;
                       fout_molecules << fixed <<DELIMITER<<setw(STRING_WIDTH)<<"HBond_Donors:"<<setw(FLOAT_WIDTH)<<next_layer[i].mol.hb_donors <<endl;
                       fout_molecules <<DELIMITER<<setw(STRING_WIDTH)<<"Layer_Completed:"<<setw(FLOAT_WIDTH)<<counter+1 <<endl;
                       fout_molecules <<next_layer[i].mol.current_data;
                       fout_molecules <<DELIMITER<<setw(STRING_WIDTH)<<"Frag_String:"<<setw(FLOAT_WIDTH)<<next_layer[i].mol.energy << "\n"<<endl;
                       Write_Mol2(next_layer[i].mol, fout_molecules);
                       write_to_file_counter++;
		       //new_comment.clear();
                    }

                    //Refinement does not display fragment string name due to spacing issues - Fixed LEP

                    // If growth_trees are turned on, write those here as well
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
            }
    
            // If dn_write_prune_dump is turned on, do that here.
            if (dn_write_prune_dump){

                // Filename for filehandle
                ostringstream fout_prune_dump_name;
                fout_prune_dump_name <<dn_output_prefix <<".anchor_" <<(a+1) <<".prune_dump_layer_" <<counter+1 <<".mol2";
                fstream fout_prune_dump;
    
                // Open the filehandle
                fout_prune_dump.open (fout_prune_dump_name.str().c_str(), fstream::out|fstream::app);

                // Check the molecules and cap as appropriate
                for (int i=0; i<prune_dump.size(); i++){

                    bool prune_dump_flag = true;
                 
                    // Iterate over each attachment point
                    for (int j=0; j<prune_dump[i].aps.size(); j++){

                        string bond_type;

                        // Loop over all bonds in frag1
                        for (int k=0; k<prune_dump[i].mol.num_bonds; k++){
                            if (prune_dump[i].mol.bonds_origin_atom[k] == prune_dump[i].aps[j].dummy_atom ||
                                prune_dump[i].mol.bonds_target_atom[k] == prune_dump[i].aps[j].dummy_atom   ){
                                if (prune_dump[i].mol.bonds_origin_atom[k] == prune_dump[i].aps[j].heavy_atom ||
                                    prune_dump[i].mol.bonds_target_atom[k] == prune_dump[i].aps[j].heavy_atom   ){

                                    // This is the bond you're looking for...
                                    bond_type = prune_dump[i].mol.bond_types[k];
                                    break;
                                }
                            }
                        }

                        // Only cap it if the dummy is on a single bond 
                        if ( bond_type == "1" ){

                            dummy_to_H(prune_dump[i], prune_dump[i].aps[j].heavy_atom, prune_dump[i].aps[j].dummy_atom);
                            prune_dump[i].mol.mol_wt += ATOMIC_WEIGHT_H;

                        // Otherwise, ignore this molecule
                        } else {
                            prune_dump_flag = false;
                            break;
                        }
                    }
                    // LEP - run prune_dump mols through multi torenv check to not make garbage
                    if (valid_torenv_multi(prune_dump[i])){

                        // If there was no problem capping the molecules, then write to file
                        if (prune_dump_flag){
                            calc_mol_wt(prune_dump[i].mol);
                            calc_rot_bonds(prune_dump[i].mol);
                            calc_formal_charge(prune_dump[i].mol);

                            //TODO there is an extra space being printed after current_data
                            fout_prune_dump << prune_dump[i].mol.current_data <<endl;
                            fout_prune_dump <<DELIMITER<<setw(STRING_WIDTH)<<"Molecular_Weight:"<<setw(FLOAT_WIDTH)<<prune_dump[i].mol.mol_wt <<endl;
                            fout_prune_dump <<DELIMITER<<setw(STRING_WIDTH)<<"DOCK_Rotatable_Bonds:"<<setw(FLOAT_WIDTH)<<prune_dump[i].mol.rot_bonds <<endl;
                            fout_prune_dump <<DELIMITER<<setw(STRING_WIDTH)<<"Formal_Charge:"<<setw(FLOAT_WIDTH)<<prune_dump[i].mol.formal_charge <<endl;
                            fout_prune_dump <<DELIMITER<<setw(STRING_WIDTH)<<"Layer_Completed:"<<setw(FLOAT_WIDTH)<<counter+1 <<"\n" <<endl;
                            Write_Mol2( prune_dump[i].mol, fout_prune_dump );
                            write_to_file_counter++;
                        }
                    }
                }
    
                // Clear some memory
                fout_prune_dump_name.clear();
                fout_prune_dump.close();
            }
    
            // clear the prune dump vector for the next round
            prune_dump.clear();
            if (dn_ga_flag == false){
               cout <<"      After writing " <<write_to_file_counter <<" of those to file, " <<root.size() 
                    <<" are returned to root for the next iteration" <<endl <<endl;
            } 
    
            // Print out all the contents of root, this will be useful for restarts CDS 
            if ((dn_write_checkpoints || counter==(dn_max_grow_layers-1))&& dn_ga_flag == false){
    
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
 
        } // end main while loop
    
    
        // clear root for the next round of anchor orients
        root.clear();
    
    } // end for each anchor
    if (dn_ga_flag == false) {
       cout <<"There are no more anchors. Exiting." <<endl <<endl;
    }

    ////////////////////////////////////////////////////////////////////////////
    // End Main Loop                                                          //
    ////////////////////////////////////////////////////////////////////////////


    // clean up any filehandles
    fout_molecules.close();

    return;

} // end DN_GA_Build::build_molecules()



// +++++++++++++++++++++++++++++++++++++++++
// This function will combine things without too much thinking
void
DN_GA_Build::simple_build( Master_Score & score, Simplex_Minimizer & simplex, AMBER_TYPER & typer )
{
//TODO
// Needs an overhaul

    cout <<endl <<endl <<"##### Entering simple build" <<endl;
    fstream fout_molecules;
    fout_molecules.open ("final.mol2", fstream::out|fstream::app);
    vector <Fragment> root;
    vector <Fragment> layer;
    vector <Fragment> next_layer;
    vector <Fragment> growing;
    int counter = 0;
    if (dn_user_specified_anchor){
        root.push_back(anchors[0]);
    }
    while (counter < dn_max_grow_layers){
        if (root.size() == 0){ cout <<"Root is empty, exiting." <<endl; break;}
        cout <<"##### Entering layer of growth #" <<counter+1 <<endl;
        cout <<"Root size at the beginning of this layer = " <<root.size() <<endl;
        for (int i=0; i<root.size(); i++){
            activate_mol(root[i].mol);
            layer.push_back(root[i]);
            int temp = layer[0].aps.size();
            for (int j=temp-1; j>-1; j--){
                for (int k=0; k<layer.size(); k++){
                    for (int l=0; l<sidechains.size(); l++){
                        for (int m=0; m<sidechains[l].aps.size(); m++){
                            Fragment frag_combined = combine_fragments( layer[k],
                                                                        layer[k].aps[j].dummy_atom,
                                                                        layer[k].aps[j].heavy_atom,
                                                                        sidechains[l],
                                                                        sidechains[l].aps[m].dummy_atom,
                                                                        sidechains[l].aps[m].heavy_atom );
                            activate_mol(frag_combined.mol);
                            growing.push_back(frag_combined);
                        }
                    }
                } 
                layer.clear();
                for (int x=0; x<growing.size(); x++){
                    if (layer.size() < dn_max_layer_size){
                        layer.push_back(growing[x]);
                    } else {
                        break;
                    }
                }
                growing.clear();
            } 
            cout <<"  From Root[" <<i <<"], created " <<layer.size()
                 <<" fragments with a new layer" <<endl;
            for (int x=0; x<layer.size(); x++){
                next_layer.push_back(layer[x]);
            }
            layer.clear();
        } 
        root.clear();
        cout <<"    From all frags in root, kept a total of " <<next_layer.size() 
             <<" fragments with a new layer" <<endl;
        for (int i=0; i<next_layer.size(); i++){
           if (dummy_in_mol(next_layer[i].mol)){
                if (root.size() < dn_max_root_size){
                     root.push_back(next_layer[i]);
                }
           } else {
                fout_molecules <<next_layer[i].mol.current_data << endl;
                Write_Mol2(next_layer[i].mol, fout_molecules);
           }
        }
        cout <<"    After writing some (or none) of those out to file, " <<root.size() 
             <<" are returned to root for the next iteration" <<endl <<endl;
        next_layer.clear();
        counter++;
    } 
    root.clear();
    fout_molecules.close();
    return;

} // end DN_GA_Build::simple_build()



// +++++++++++++++++++++++++++++++++++++++++
// This function will orient a fragment and return filtered / clustered / pruned results on the
// root vector
void
DN_GA_Build::orient_fragments( vector <Fragment> & root, Fragment & frag, Master_Score & score,
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
    //sort(orients.begin(), orients.end(), ga_fragment_sort);
    frag_sort(orients, ga_fragment_sort);

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

} // end DN_GA_Build::orient_fragments()



// +++++++++++++++++++++++++++++++++++++++++
// Sample the contents of a given fragment library exhaustively
void
DN_GA_Build::sample_fraglib_exhaustive( Fragment & layer_frag, int j, vector <Fragment> & fraglib,
                                     vector <Fragment> & growing_ref, Master_Score & score, 
                                     Simplex_Minimizer & simplex, AMBER_TYPER & typer, bool last_layer)
{

    // If the fragment library is empty, return without doing anything
    if (fraglib.size() == 0) { return; }

    // For every fragment in the fragment library (i.e. sidechains, linkers, scaffolds)
    for (int l=0; l<fraglib.size(); l++){

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
            if (compare_dummy_bonds(layer_frag, layer_frag.aps[j].dummy_atom,
                                    layer_frag.aps[j].heavy_atom, fraglib[l],
                                    fraglib[l].aps[m].dummy_atom, fraglib[l].aps[m].heavy_atom )){

                // If it is compatible, combine the two fragments to make a new, larger fragment
                Fragment frag_combined = combine_fragments(layer_frag, layer_frag.aps[j].dummy_atom,
                                                           layer_frag.aps[j].heavy_atom, fraglib[l],
							   fraglib[l].aps[m].dummy_atom,
                                                           fraglib[l].aps[m].heavy_atom);

                // check growing_ref to see if it gets bigger
                int growing_ref_size_before = growing_ref.size();        

                // If the torenv is turned on...
                if (dn_use_torenv_table){
                    bool valid_torenv_bool = false;
                    // First check to see if the newest torsion environment is allowable
                    // (note that there may be some dummy atoms included in the environment, and
                    // even though it passes a check now, it may not pass the check later)
                    if (dn_ga_flag==false){ 
                        if(dn_use_roulette){
                            valid_torenv_bool=roulette_valid_torenv(frag_combined);
                        } else {
                            valid_torenv_bool=valid_torenv(frag_combined);
                        }
                    }                     

                    if (valid_torenv_bool){
                        // Second check to see whether all of the old torsion environments are
                        // still compatible. This is necessary because of the dummy as wild card
                        // thing
                        
                        if (dn_ga_flag){
                        // CDS - Prepare the molecule using the amber_typer to assign bond types to each bond
                           frag_combined.mol.prepare_molecule();
                           typer.skip_verbose_flag = true;
                           typer.prepare_molecule(frag_combined.mol, true, score.use_chem, score.use_ph4, score.use_volume);
                           GA_Recomb c_ga;                                                    
                           c_ga.prepare_torenv_indices(frag_combined);
                        }

                        if (valid_torenv_multi(frag_combined)){

                            // If both checks pass, then this is a valid molecule. Sample torsions.
                            sample_minimized_torsions(frag_combined, growing_ref, score, simplex, typer);
                        }
                    }

                // If the torenv is not turned on, just skip straight to sampling torsions
                } else {
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
 
} // end DN_GA_Build::sample_fraglib_exhaustive()



// +++++++++++++++++++++++++++++++++++++++++
// Sample the contents of a given fragment library randomly
void
DN_GA_Build::sample_fraglib_rand( Fragment & layer_frag, int j, vector <Fragment> & fraglib,
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

    // This loop will iterate dn_num_random_picks times, each time trying to add exactly one new
    // frgament. The counter increments no matter whether the fragment is kept or not.
    while (choice < dn_num_random_picks){

        // The passes int counts every time we either keep or skip a fragment from the libary, so
        // once its size reaches the fraglib.size(), we have looked at everything
        if (passes >= fraglib.size()){ return; }

        // Choose a random starting position, fraglib[dnrn]
        // (rand has been seeded if and only if the minimizer is turned on)
        dnrn=rand() % fraglib.size();
        cout << endl <<  "DN random number: "<< dnrn << endl; 
        //LEP - verbose dn frag selection random seed output
        cout << endl << "De novo fragment selection: " << fraglib[dnrn].mol.energy << endl; 
        // Keep track of things we try so we don't try the same thing twice
        if (picked.size() > 0){
            if (picked.size() >= fraglib.size()){
                cout << "Number of selections equals size of fraglib - exiting fraglib sampling." << endl;
               picked.clear();
               return; 
            }
            for (int p=0; p<picked.size(); p++){
                if (dnrn == picked[p]){ 
                    cout << "Tried this fragment before." << endl;
                    continue;
                }
            }
        }
        picked.push_back(dnrn);
        

        // (1) Would create too many scaffolds per layer
        if ( fraglib[dnrn].aps.size() > 2 &&
             layer_frag.scaffolds_this_layer >= dn_max_scaffolds_per_layer ){ passes++;  
             cout << "Does not pass scaffolds per layer cutoff." << endl;
             continue;
        }

        // (2) Would give too many current attachment points
        int temp_aps = fraglib[dnrn].aps.size() + layer_frag.aps.size() - 2;
        if (temp_aps > dn_max_current_aps) { passes++; 
            cout << "Does not pass attachment points per layer cutoff." << endl;
            continue;
        }

        // (3) Rotatable bond cutoff
        // BCF assuming rot_bon does not take into account attachment points
        // will have at least one new rot bond per attachment point, on each fragment
        // except when combined two aps will become one rot bond.
        // Floor for remaining rot_bonds if all remaining aps are capped with sidechains.
        int temp_rbs = layer_frag.mol.rot_bonds + temp_aps;
        if (temp_rbs > dn_constraint_rot_bon) { passes++; 
            cout << "Does not pass rotatable bond cutoff." << endl;
            continue;
         }

        // (4) Molecular weight cutoff (MW computed for layer frag at last step, MW precomputed for all frags)
        //int temp_mw = fraglib[dnrn].mol.mol_wt + layer_frag.mol.mol_wt;
        //if (temp_mw > dn_constraint_mol_wt) { passes++; continue; }

        // (5) Formal Charge Cutoff
        int temp_fc = fraglib[dnrn].mol.formal_charge + layer_frag.mol.formal_charge;
        if (fabs(float(temp_fc)) > (fabs(float(dn_constraint_formal_charge))+0.1)) { passes++; 
            cout << "Does not pass formal charge cutoff." << endl;
            continue;
        }
        // (6) Torsion environment prune is done later depending on user input


        // If this is the last layer, we do not want to sample anything with more than 1 att point
        if ( last_layer && fraglib[dnrn].aps.size() > 1 ){
            passes++;
            //continue;
            cout << "Do not want to sample anything with more than 1 attachment point." << endl;
            continue;
        }

        if (iso_picks>0){
            std::stringstream newName;
            newName << "h_"<< fraglib[dnrn].mol.energy;
            for (int s =0; s<fraglib[dnrn].mol.num_atoms; s++){
                fraglib[dnrn].mol.subst_names[s] = newName.str();
            }
        }

        cout << "Number of attachment points of frag: " << fraglib[dnrn].aps.size() << endl;
        // For each attachement point on that fragment
        for (int m=0; m<fraglib[dnrn].aps.size(); m++){
            //cout << "number of attachment points of frag: " << fraglib[dnrn].aps.size() << endl;
            //cout << "Trying attachment point #" << m << endl; //LEP - 03.24.18

            // Check to see if the bond order is compatible
            cout << "Checking bond order compatibility..." << endl;
            if (compare_dummy_bonds(layer_frag, layer_frag.aps[j].dummy_atom, 
                                    layer_frag.aps[j].heavy_atom, fraglib[dnrn],
                                    fraglib[dnrn].aps[m].dummy_atom,
                                    fraglib[dnrn].aps[m].heavy_atom )){
                
                cout << "Passed bond order compatibility check." << endl;
                // If it is compatible, combine the two fragments to make a new, larger fragment
                Fragment frag_combined = combine_fragments(layer_frag, layer_frag.aps[j].dummy_atom,
                                                           layer_frag.aps[j].heavy_atom,
                                                           fraglib[dnrn],
                                                           fraglib[dnrn].aps[m].dummy_atom,
                                                           fraglib[dnrn].aps[m].heavy_atom);
                if (frag_combined.iso_targeted_AP >= 0 
                    || frag_combined.iso_frag_att_num_APs >= 0){
                    std::cout << "ERROR: iso_targeted_AP is not -1" << std::endl;
                    std::cout << "ERROR: iso_frag_att_num_APs is not -1" << std::endl;
                    exit(1);
                }

                // get the attachemnt point that was targeted
                frag_combined.iso_targeted_AP = j;
 
                // get the num of AP for the frag that was recently attached
                // for the isoswap_GA
                frag_combined.iso_frag_att_num_APs = fraglib[dnrn].aps.size();
  
                // get the index of the fragment used
                frag_combined.iso_head_attached_ind = dnrn;

                // check growing_ref to see if it gets bigger
                int growing_ref_size_before = growing_ref.size();        

                // If the torenv is turned on...
                if (dn_use_torenv_table){
                    bool valid_torenv_bool=false;
                    // First check to see if the newest torsion environment is allowable
                    // (note that there may be some dummy atoms included in the environment, and
                    // even though it passes a check now, it may not pass the check later)
                    if(dn_use_roulette){
                        valid_torenv_bool=roulette_valid_torenv(frag_combined);
                    } else {
                        cout << "Checking single torenv..." << endl;
                        valid_torenv_bool=valid_torenv(frag_combined);
                    }
                    if (valid_torenv_bool){
                        cout << "Passed single torenv check." << endl;

                        // Second check to see whether all of the old torsion environments are
                        // still compatible. This is necessary because of the dummy as wild card
                        // thing
                        if (dn_ga_flag){
                        // CDS - Prepare the molecule using the amber_typer to assign bond types to each bond
                           frag_combined.mol.prepare_molecule();
                           typer.skip_verbose_flag = true;
                           typer.prepare_molecule(frag_combined.mol, true, score.use_chem, score.use_ph4, score.use_volume);
                           GA_Recomb c_ga;                                                    
                           c_ga.prepare_torenv_indices(frag_combined);
                        }
                        cout << "Checking multi torenv..." << endl;
                        if (valid_torenv_multi(frag_combined)){
                            cout << "Passed multi torenv check" << endl;
                            // If both checks pass, then this is a valid molecule. Sample torsions.
                            sample_minimized_torsions(frag_combined, growing_ref, score, simplex, typer);
                            //cout << "size of growing ref: " << growing_ref.size() << endl;
                        }
                        else{
                            cout << "Did not pass multi torenv check." << endl;
                        }    
                    }
                    else{
                        cout << "Did not pass single torenv check." << endl;
                    }

                // If the torenv is not turned on, just skip straight to sampling torsions
                } else {
                    sample_minimized_torsions(frag_combined, growing_ref, score, simplex, typer);
                }

                // If we got this far, and if the fragment we are looking at is a scaffold,
                // and if growing_ref increased in size, then we added a scaffold this layer
                if ( (fraglib[dnrn].aps.size() > 2) && 
                     (growing_ref.size() > growing_ref_size_before)){ 

                    layer_frag.scaffolds_this_layer += 1;
                }
            }else{
                cout << "Did not pass bond order compatibility check." << endl;
            }
        }

        // Increment choice whether a new frag was generated or not
        choice++;
        passes++;

    } // ... end While N choices

    picked.clear();
    return;

} // end DN_GA_Build::sample_fraglib_rand()



// +++++++++++++++++++++++++++++++++++++++++
// Sample the contents of a given fragment library following the FragGraph
void
DN_GA_Build::sample_fraglib_graph( Fragment & layer_frag, int j, vector <Fragment> & fraglib,
                                vector <FragGraph> & fraggraph_vec, vector <Fragment> & growing_ref,
                                Master_Score & score, Simplex_Minimizer & simplex,
                                AMBER_TYPER & typer, bool last_layer )
{
    //cout << "-SAMPLE_FRAGLIB_GRAPH-" << endl;
    // If the fragment library is empty, return without doing anything
    if (fraglib.size() == 0) { return; }

    //BCF Calc MW RB and FC for pruning later
    calc_mol_wt(layer_frag.mol);
    calc_rot_bonds(layer_frag.mol);
    calc_formal_charge(layer_frag.mol);

    // If this is the last layer, we do not want to sample anything with more than 1 att point
    if (last_layer){
        for (int f=0; f<fraglib.size(); f++){
            if (fraglib[f].aps.size() > 1) { fraggraph_vec[f].visited = true; }
        }
    }


    // For a certain number of random starting positions
    int i = 0;
    while (i<dn_graph_max_picks){
        bool counted     = false; // BCF 4-19 do not want to count multiple aps of the same fragment as different tries
        // Choose a random starting position, fraglib[dnrn]
        // (rand has been seeded if and only if the minimizer is turned on)
        int dnrn=rand() % fraglib.size();

        // Check to see if it has been visited
        if (fraggraph_vec[dnrn].visited){

            int temp_counter = 1;

            // If it has, pick a new starting point, keep picking until finding a valid one
            while ( temp_counter < fraglib.size() ){

                // This check will stop incrementing dnrn when it reaches the end of fraglib,
                // then start over at zero
                dnrn++;
                if ( dnrn > (fraglib.size()-1) ){ dnrn=0; }

                // If it finds one that is not yet visited, continue with that dnrn
                if ( !fraggraph_vec[dnrn].visited ){ break; }
                else { temp_counter++; }

                // If there are no un-visited fragments, reset the flags and return
                if (temp_counter == fraglib.size()){
                    for (int r=0; r<fraggraph_vec.size(); r++){ fraggraph_vec[r].visited = false; }
                    //cout << "visited all" << endl;
                    return;
                }
            }
        }

        // Pair the starting position with the current score (1/19/16 from layer_frag), add it to a
        // sampling queue for a breadth-first search
        vector < pair <int, float> > sampling_queue;
        pair <int, float> starting_point = make_pair( dnrn, layer_frag.mol.current_score );
        sampling_queue.push_back( starting_point );

        // Limit the amount of sampling performed
        int tries = 0;
        int max_tries = pow(float(dn_graph_breadth), dn_graph_depth);

        // Keep searching along additional fragment lines that improve the score
        float compare_score = starting_point.second;
        bool improvement_flag = false;

        // While there are still things to try in the sampling_queue
        while (sampling_queue.size() > 0 && tries < max_tries){

            // pull off the first value, delete it from the queue
            pair <int, float> this_dnrn;
            this_dnrn = sampling_queue.front();
            sampling_queue.erase( sampling_queue.begin() );

            // Some things can be added to the queue more than once - e.g. they are already in the
            // queue, but not yet visited, so they are added again.
            if (fraggraph_vec[this_dnrn.first].visited == true){ continue; }

            // Mark this position as visited in the FragGraph, we are now "trying" a fragment
            fraggraph_vec[this_dnrn.first].visited = true;

            /***********************************************
 *            Now going through all the reasons to not add a fragment
 *           */

            // (1) Would create too many scaffolds per layer
            if ( fraglib[this_dnrn.first].aps.size() > 2 && 
                 layer_frag.scaffolds_this_layer >= dn_max_scaffolds_per_layer ){ continue; }

            // (2) Would give too many current attachment points  
            int temp_aps = fraglib[this_dnrn.first].aps.size() + layer_frag.aps.size() - 2;
            if (temp_aps > dn_max_current_aps) { continue; }

            // (3) Rotatable bond cutoff   
            // BCF assuming rot_bon does not take into account attachment points
            // will have at least one new rot bond per attachment point, on each fragment
            // except when combined two aps will become one rot bond.
            // Floor for remaining rot_bonds if all remaining aps are capped with sidechains.
            int temp_rbs = layer_frag.mol.rot_bonds + temp_aps;
            //cout << layer_frag.mol.rot_bonds << endl;
            if (temp_rbs > dn_constraint_rot_bon) { continue; }
            
            // (4) Molecular weight cutoff (MW computed for layer frag at last step, MW precomputed for all frags)
            //int temp_mw = fraglib[this_dnrn.first].mol.mol_wt + layer_frag.mol.mol_wt;
            //if (temp_mw > dn_constraint_mol_wt) { continue; } 

            // (5) Formal Charge Cutoff
            int temp_fc = fraglib[this_dnrn.first].mol.formal_charge + layer_frag.mol.formal_charge;
            if (fabs(float(temp_fc)) > (fabs(float(dn_constraint_formal_charge))+0.1)) { continue; }
            // (6) Torsion environment prune is done later depending on user input
           

            // Print statements to visualize which frags have been visited in the graph
            //for (int m=0; m<fraggraph_vec.size(); m++){ cout <<fraggraph_vec[m].visited; } cout <<endl;

            // For each attachement point on that fragment
            counted = false; //want to reset for each unique mol but not for diff aps on same mol
            for (int m=0; m<fraglib[this_dnrn.first].aps.size(); m++){
        
                // Check to see if the bond order is compatible
                if (compare_dummy_bonds(layer_frag, layer_frag.aps[j].dummy_atom,
                                        layer_frag.aps[j].heavy_atom, fraglib[this_dnrn.first],
                                        fraglib[this_dnrn.first].aps[m].dummy_atom,
                                        fraglib[this_dnrn.first].aps[m].heavy_atom )){
        
                    // If so, attach it
                    Fragment frag_combined = combine_fragments(layer_frag, layer_frag.aps[j].dummy_atom,
                                                               layer_frag.aps[j].heavy_atom,
                                                               fraglib[this_dnrn.first],
                                                               fraglib[this_dnrn.first].aps[m].dummy_atom,
                                                               fraglib[this_dnrn.first].aps[m].heavy_atom);

                    if (frag_combined.iso_targeted_AP >= 0 
                        || frag_combined.iso_frag_att_num_APs >= 0){
                        std::cout << "ERROR: iso_targeted_AP is not -1" << std::endl;
                        std::cout << "ERROR: iso_frag_att_num_APs is not -1" << std::endl;
                        exit(1);
                    }

                    // get the attachemnt point that was targeted
                    frag_combined.iso_targeted_AP = j;
 
                    // get the num of AP for the frag that was recently attached
                    // for the isoswap_GA
                    frag_combined.iso_frag_att_num_APs = fraglib[dnrn].aps.size();
  
                    // get the index of the fragment used
                    frag_combined.iso_head_attached_ind = dnrn;

                    // check growing_ref to see if it gets bigger
                    int growing_ref_size_before = growing_ref.size();        

                    // If the torenv is turned on...
                    if (dn_use_torenv_table){
                        bool valid_torenv_bool=false;
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
                            if (dn_ga_flag){
                               // Prepare the molecule using the amber_typer to assign bond types to each bond
                               frag_combined.mol.prepare_molecule();
 
                               typer.skip_verbose_flag = true;
                               typer.prepare_molecule(frag_combined.mol, true, score.use_chem, score.use_ph4, score.use_volume);
                               GA_Recomb c_ga;
                               c_ga.prepare_torenv_indices(frag_combined);                                                                            }
                            if (valid_torenv_multi(frag_combined)){

                                // If both checks pass, then this is a valid molecule. Sample torsions.
                                sample_minimized_torsions(frag_combined, growing_ref, score, simplex, typer);
                                if (!counted ) { //BCF 
                                    i++; //increment random pick iterator
                                    tries++; //increment graph pick iterator
                                    counted = true; //don't want to count multiple aps on one frag
                                    //cout << i << endl;
                                }
                            }
                        }

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

} // end DN_GA_Build::sample_fraglib_graph()



// +++++++++++++++++++++++++++++++++++++++++
// Given two fragments and connection point data, combine them into one and return it
Fragment
DN_GA_Build::combine_fragments( Fragment & frag1, int dummy1, int heavy1,
                             Fragment frag2, int dummy2, int heavy2 )
{
    double finalmat[3][3];
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
        //finalmat[3][3] = { { -1, 0, 0}, {0, -1, 0}, {0, 0, -1} };
        finalmat[0][0]=-1; finalmat[0][1]=0;  finalmat[0][2]=0;
        finalmat[1][0]=0;  finalmat[1][1]=-1; finalmat[1][2]=0;
        finalmat[2][0]=0;  finalmat[2][1]=0;  finalmat[2][2]=-1;

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
        //double finalmat[3][3];

        finalmat[0][0]=1;  finalmat[0][1]=0;  finalmat[0][2]=0;
        finalmat[1][0]=0;  finalmat[1][1]=1;  finalmat[1][2]=0;
        finalmat[2][0]=0;  finalmat[2][1]=0;  finalmat[2][2]=1;

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

    // Save the orientation mat for isoswap
    for(int i=0;i<3;i++){
        for( int j=0;j<3;j++){
            frag_combined.iso_ori_mat[i][j] = finalmat[i][j];
        }
    }  

    return frag_combined;

} // end DN_GA_Build::combine_fragments()



// +++++++++++++++++++++++++++++++++++++++++
// Return a covalent radius depending on atom type 
float
DN_GA_Build::calc_cov_radius( string atom )
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
               <<" in DN_GA_Build::calc_cov_radius()"  <<endl; 
          return 0.71; }

} // end DN_GA_Build::calc_cov_radius()



// +++++++++++++++++++++++++++++++++++++++++
// Attach two fragments that are already overlaid at an attachment point
Fragment
DN_GA_Build::attach( Fragment & frag1, int dummy1, int heavy1, 
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
    //ostringstream new_name;
    //new_name <<frag1.mol.energy <<"--" <<frag2.mol.energy;
    //combined_mol.energy = new_name.str();
    //if (verbose) cout << "Trying fragment: " << frag2.mol.energy << endl; //LEP - 03.24.18
    //new_name.clear();

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

    // Clear some vectors
    tmp_dummy_atoms.clear();
    tmp_heavy_atoms.clear();

    //BCF new fragment needs to remember how many scaffolds were added this layer
    //it will be reset at the start of a new layer
    final.scaffolds_this_layer = frag1.scaffolds_this_layer;

    return final;

} // end DN_GA_Build::attach()



// +++++++++++++++++++++++++++++++++++++++++
// Check new atom environments against table of allowable environments
bool
DN_GA_Build::valid_torenv( Fragment & frag )
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

    if (verbose) { cout <<"VT: input torenv is " <<tmp_pair.first <<"-" <<tmp_pair.second <<endl;}
    if (verbose) { cout <<"VT: Starting a new check" <<endl; }
    if (verbose) { cout <<"VT:   tmp_pair.first.size() = " <<tmp_pair.first.size() <<endl; }
    // For every char in the string
    for (int i=0; i<tmp_pair.first.size(); i++){
        // If a dummy atom is encountered (Y in fingerprint.cpp)
        if (verbose) { cout <<"VT:     tmp_pair.first[i] =" <<tmp_pair.first[i] <<endl; }

        if (tmp_pair.first[i] == 'Y'){
            if (verbose) { cout <<"VT:       found a match...breaking" <<endl; }
            // Then make this flag true and break the loop
            wc_invoked = true;
            break;
        }
    }

    if (verbose) { cout <<"VT:   tmp_pair.second.size() = " <<tmp_pair.second.size() <<endl; }

    // Iterate over the second string in the same way if needed
    if (!wc_invoked) {
        for (int i=0; i<tmp_pair.second.size(); i++){
            if (verbose) { cout <<"VT:     tmp_pair.second[i] =" <<tmp_pair.second[i] <<endl; }
            if (tmp_pair.second[i] == 'Y'){
            if (verbose) { cout <<"VT:       found a match...breaking" <<endl; }
                wc_invoked = true;
                break;
            }
        }
    }

    if (verbose) { cout <<"VT:         After both checks, wc_invoked = " <<wc_invoked <<endl <<endl; }


    // The returned tmp_pair is alphabetized, as is the torenv vector, so you will only ever have
    // to check tmp_pair.first against the origin_env
    for (int i=0; i<torenv_vector.size(); i++){

        if (verbose) {cout <<"VT: Checking to see if " <<tmp_pair.first <<" is in the fragment library " <<endl;}

        // If the first environment of the pair finds a match in the torenv_vector...
        if (compare_atom_environments(torenv_vector[i].origin_env, tmp_pair.first)){

            if (verbose){ cout <<"VT: Yep it's in there, looking for a match... " <<endl;}

            // ...then start checking to see if the second environment in the pair has a 
            // corresponding match
            for (int j=0; j<torenv_vector[i].target_envs.size(); j++){

                if (verbose){ cout <<"VT:     Now checking if it is connected to a " <<tmp_pair.second <<endl;}

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

        if (verbose) {cout <<"VT: Checking to see if " <<tmp_pair.second <<" is in the fragment library " <<endl;}

        if (compare_atom_environments(torenv_vector[i].origin_env, tmp_pair.second)){

            if (verbose){ cout <<"VT: Yep it's in there, looking for a match... " <<endl;}

            for (int j=0; j<torenv_vector[i].target_envs.size(); j++){

                if (verbose){ cout <<"VT:     Now checking if it is connected to a " <<tmp_pair.first <<endl;}

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

} // end DN_GA_Build::valid_torenv()




// +++++++++++++++++++++++++++++++++++++++++
// Check new atom environments against table of allowable environments
// and accept based on roulette probability
bool
DN_GA_Build::roulette_valid_torenv( Fragment & frag )
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
    
    if (verbose) { cout <<"VT: input torenv is " <<tmp_pair.first <<"-" <<tmp_pair.second <<endl;}
    if (verbose) { cout <<"VT: Starting a new check" <<endl; }
    if (verbose) { cout <<"VT:   tmp_pair.first.size() = " <<tmp_pair.first.size() <<endl; }
    // For every char in the string
    for (int i=0; i<tmp_pair.first.size(); i++){
        // If a dummy atom is encountered (Y in fingerprint.cpp)
        if (verbose) { cout <<"VT:     tmp_pair.first[i] =" <<tmp_pair.first[i] <<endl; }
        
        if (tmp_pair.first[i] == 'Y'){
            if (verbose) { cout <<"VT:       found a match...breaking" <<endl; }
            // Then make this flag true and break the loop
            wc_invoked = true;
            break;
        }
    }
    
    if (verbose) { cout <<"VT:   tmp_pair.second.size() = " <<tmp_pair.second.size() <<endl; }
    
    // Iterate over the second string in the same way if needed
    if (!wc_invoked) {
        for (int i=0; i<tmp_pair.second.size(); i++){
            if (verbose) { cout <<"VT:     tmp_pair.second[i] =" <<tmp_pair.second[i] <<endl; }
            if (tmp_pair.second[i] == 'Y'){
                if (verbose) { cout <<"VT:       found a match...breaking" <<endl; }
                wc_invoked = true;
                break;
            }
        }
    }
    
    if (verbose) { cout <<"VT:         After both checks, wc_invoked = " <<wc_invoked <<endl <<endl; }
    
    
    // The returned tmp_pair is alphabetized, as is the torenv vector, so you will only ever have
    // to check tmp_pair.first against the origin_env
    bool first_tors_check=false, second_tors_check=false;
    tmp_first_check=false;
    tmp_second_check=false;
    for (int i=0; i<torenv_vector.size(); i++){
        
        if (verbose) {cout <<"VT: Checking to see if " <<tmp_pair.first <<" is in the fragment library " <<endl;}
        
        // If the first environment of the pair finds a match in the torenv_vector...
        tmp_first_check=true;
        if (compare_atom_environments(torenv_vector[i].origin_env, tmp_pair.first)){
            
            if (verbose){ cout <<"VT: Yep it's in there, looking for a match... " <<endl;}
            
            // ...then start checking to see if the second environment in the pair has a
            // corresponding match
            tmp_first_check=false;
            tmp_second_check=true;
            for (int j=0; j<torenv_vector[i].target_envs.size(); j++){
                
                if (verbose){ cout <<"VT:     Now checking if it is connected to a " <<tmp_pair.second <<endl;}
                
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
            
            if (verbose) {cout <<"VT: Checking to see if " <<tmp_pair.second <<" is in the fragment library " <<endl;}
            
            tmp_second_check=true;
            if (compare_atom_environments(torenv_vector[i].origin_env, tmp_pair.second)){
                
                if (verbose){ cout <<"VT: Yep it's in there, looking for a match... " <<endl;}
                
                for (int j=0; j<torenv_vector[i].target_envs.size(); j++){
                    
                    if (verbose){ cout <<"VT:     Now checking if it is connected to a " <<tmp_pair.first <<endl;}
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
            roulette_compare_environment=tmp_pair.second + "-" + tmp_pair.first;
        }
    }
    
    //initialize variables - frequency is from the table, compare_frequency is a random number
    double roulette_frequency; 
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
     0     1 ^ 2 3
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
    
} // end DN_GA_Build::roulette_valid_torenv()



// +++++++++++++++++++++++++++++++++++++++++
// Check new atom environments against table of allowable environments
bool
DN_GA_Build::valid_torenv_multi( Fragment & frag )
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
            if (verbose){ cout <<"VTM: Checking to see if " <<tmp_pair.first <<" is in the fragment library " <<endl;}
    
            // If the first environment of the pair finds a match in the torenv_vector...
            if (compare_atom_environments(torenv_vector[j].origin_env, tmp_pair.first)){
    
                if (verbose){ cout <<"VTM: Yep it's in there, looking for a match... " <<endl; }
    
                // ...then start checking to see if the second environment in the pair has a 
                // corresponding match
                for (int k=0; k<torenv_vector[j].target_envs.size(); k++){
    
                    if (verbose){ cout <<"VTM:     Now checking if it is connected to a " <<tmp_pair.second <<endl; }
    
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
                if (verbose){ cout <<"VTM: Checking to see if " <<tmp_pair.second <<" is in the fragment library " <<endl;}

                if (compare_atom_environments(torenv_vector[j].origin_env, tmp_pair.second)){
 
                    if (verbose){ cout <<"VTM: Yep it's in there, looking for a match... " <<endl; }

                    for (int k=0; k<torenv_vector[j].target_envs.size(); k++){

                        if (verbose){ cout <<"VTM:     Now checking if it is connected to a " <<tmp_pair.first <<endl; }

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

} // end DN_GA_Build::valid_torenv_multi()



// +++++++++++++++++++++++++++++++++++++++++
// Function to compare two atom environments assuming that dummy atom Y is a wildcard
bool
DN_GA_Build::compare_atom_environments( string s1, string s2 )
{
    // s1 is the reference string from the vector of torenvs
    // s2 is a candidate that we are checking

    if (verbose) {cout <<"CAE: Comparing s1 = " <<s1 <<"\t and s2 = " <<s2 <<endl; }

    // If they are exactly the same, then it's a match
    if (s1 == s2){
        if (verbose) { cout <<"CAE: Matched s1 and s2: Strings are identical" <<endl <<endl; }
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


    if (verbose) { cout <<"CAE: Comparing s1_new = " <<s1_new <<"\t and s2_new = " <<s2_new <<endl; }

    // Compare contents of s1 and s2
    for (int i=0; i<s1_new.size(); i+=2){
        if ( (s1_new.c_str()[i] == s2_new.c_str()[count]) && 
             (s1_new.c_str()[i+1] == s2_new.c_str()[count+1]) ){
            if (verbose){
                cout <<"CAE: Matched  s1_new[" <<i <<"] = " <<s1_new.c_str()[i] 
                     <<" and s2_new[" <<count <<"] = " <<s2_new.c_str()[count] 
                     <<endl;
                cout <<"CAE: Matched  s1_new[" <<i+1 <<"] = " <<s1_new.c_str()[i+1] 
                     <<" and s2_new[" <<count+1 <<"] = " <<s2_new.c_str()[count+1] 
                     <<endl;
            }
            count+=2;
        } else if (count_dummie > 0){
            if (verbose) {cout <<"CAE: Matched  s1_new[" <<i <<"] = " <<s1_new.c_str()[i] <<" and dummy" <<endl;}
            if (verbose) {cout <<"CAE: Matched  s1_new[" <<i+1 <<"] = " <<s1_new.c_str()[i+1] <<" and dummy" <<endl;}
            count_dummie--;
        } else {
            if (verbose) {cout <<"CAE: Did not match s1 =\t" <<s1 <<"\tand s2 = " <<s2 <<endl <<endl;}
            return false;
        }
    }

    if (verbose){ cout <<"CAE: Matched s1 =\t" <<s1 <<"\tand s2 = " <<s2 <<endl <<endl;}
    
    //JDB logic check for which variable is being passed, and appropriate variable assignment
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
    return true;

} // end DN_GA_Build::compare_atom_environments()



// +++++++++++++++++++++++++++++++++++++++++
// Check to see if two dummy bond connections are the same bond order
bool
DN_GA_Build::compare_dummy_bonds( Fragment frag1, int dummy1, int heavy1, 
                               Fragment & frag2, int dummy2, int heavy2 )
{
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

} // end DN_GA_Build::compare_dummy_bonds();



// +++++++++++++++++++++++++++++++++++++++++
// Given a fragment and an empty vector of fragments, sample torsions for the most recent addition
// and return the results in the vector.
void
DN_GA_Build::sample_minimized_torsions( Fragment & frag1, vector <Fragment> & list_of_frags, 
                                     Master_Score & score, Simplex_Minimizer & simplex, 
                                     AMBER_TYPER & typer, std::pair<bool,float> tors_spec
 )
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

    typer.skip_verbose_flag = true;
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
    if (!tors_spec.first){
        frag_torsion_drive(frag1, tmp_list);
    } else{
       frag_torsion_drive(frag1, tmp_list, tors_spec);
    }


    std::vector <DOCKMol> v_frag_writeout {};

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
   
        // for verbose out
        if ( dn_ga_verbose > 1 ){
            v_frag_writeout.push_back(tmp_list[i].mol);
        }

        // Flexible minimization
        simplex.minimize_flexible_growth(tmp_list[i].mol, score, bond_tors_vectors);

        // Compute internal energy and primary score, store it in the dockmol
        if (score.use_primary_score) {
            valid_orient = score.compute_primary_score(tmp_list[i].mol);
            //cout << "Is this a valid orient? - " << valid_orient << endl;
        } else {
            valid_orient = true;
            tmp_list[i].mol.current_score = 0;
            tmp_list[i].mol.internal_energy = 0;
        }

        // If the energy or if the internal energy is greater than the cutoff, erase that conformer
        if ( tmp_list[i].mol.current_score > dn_pruning_conformer_score_cutoff  || 
             !valid_orient ){
            cout << "Did not pass internal energy cutoff" << endl;
            tmp_list[i].used = true;

// BCF part
            //tmp_list.erase(tmp_list.begin()+i);
            //i--;
        } else {
         
              if (use_internal_energy) {
                      if ( tmp_list[i].mol.internal_energy > ie_cutoff ) {
                          cout << "Did not pass internal energy cutoff." << endl;
                          tmp_list[i].used = true;
                          //tmp_list.erase(tmp_list.begin()+i);
                          //i--;
                      }
              }
         }
// end BCF part
        bond_tors_vectors.clear();
    }


    if ( dn_ga_verbose > 1 && iso_picks > 0 ){
        ostringstream fout_molecules_name;
        std::ostringstream gen;

        gen << std::setw(4) << std::setfill('0') << dn_ga_gen+1;
        fout_molecules_name << dn_output_prefix << ".isoswap_addtn" << gen.str() <<".mol2";
        
        print_molecules(
            fout_molecules_name.str(),
            v_frag_writeout,
            dn_ga_gen,
            ""
        );
    }

    v_frag_writeout.clear();

    // Cluster post-minimization torsions by RMSD and remove those under a certain cutoff 
    float rmsd;
    cout << "Molecules prior to clustering: "<< tmp_list.size() << endl;
    cout << "Entering clustering of post-min torsions by RMSD" << endl;
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
                    //cout << "RMSD: " << rmsd << endl;

                    // If the rmsd is below the specified cut-off, they are in the same 'cluster'
                    if ((float) j / rmsd > dn_pruning_clustering_cutoff) {
                                cout << "Molecules clustered." << endl;
                                tmp_list[j].used = true;
                    }
                }        
            }
        }
    }

//TODO
// the above comparison is checking heavy atom RMSD for the WHOLE MOLECULE. In A&G, they only check heavy atom RMSD for the newly-added segment. We need to reconcile that

    int pos_counter = 0;
    // At this point, only cluster heads will not be flagged 'used', add those to list_of_frags
    for (int i=0; i<tmp_list.size(); i++){
        if (!tmp_list[i].used){
            list_of_frags.push_back(tmp_list[i]);
            pos_counter++;
        }
    }

    // Clear the memory of tmp_list
    tmp_list.clear();

    cout << "Molecules retained: " << pos_counter << endl;
    
    return;

} // end DN_GA_Build::sample_minimized_torsions()



// +++++++++++++++++++++++++++++++++++++++++
// Compute the RMSD between two fragments considering all atoms except H
// (equivalent to heavy atoms + dummy atoms)
float
DN_GA_Build::calc_fragment_rmsd( Fragment & a, Fragment & b )
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

} // end DN_GA_Build::calc_fragment_rmsd()



// +++++++++++++++++++++++++++++++++++++++++
// Populate a list of frags with different torsions. In A&G, this is called segment_torsion_drive
void
DN_GA_Build::frag_torsion_drive( Fragment & new_frag, 
                                 vector <Fragment> & frag_list, 
                                 std::pair<bool,float> tors_spec
                                 )
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
             << "DN_GA_Build::frag_torsion_drive" <<endl;
        return;
    }

    // Make sure there is no overlap in atom assignment - This might happen for adding new 
    // fragments with only two atoms (Du-X, for example). In that case, no torsions need to be
    // sampled.
    if (atom1 == atom3 || atom4 == atom2) {
        cout <<"Warning: atoms were assigned in a weird way in DN_GA_Build::frag_torsion_drive" <<endl;
        cout <<"Perhaps a fragment was added that does not need to sample torsions." <<endl;
        return;
    }


    if ( tors_spec.first ) {
        new_frag.iso_tors_turned = tors_spec.second;
        new_frag.mol.set_torsion(atom1, atom2, atom3, atom4, tors_spec.second);
        frag_list.push_back( new_frag ); 
        return;
    } else {
        // Number of torsions comes from amber_bt_torsion_total (returns a different number for a given
        // bond type)
        num_torsions = new_frag.mol.amber_bt_torsion_total[ new_frag.mol.num_bonds-1 ];

        // Create torsions at appropriate intervals
        for (int i=0; i<num_torsions; i++){

            // Figure out the angle for this torsion and convert it to radians
            new_angle = new_frag.mol.amber_bt_torsions[ new_frag.mol.num_bonds-1 ][i];
            new_angle = (PI / 180) * new_angle;

            if (verbose){
                cout <<"FTD: amber_bt_torsions[" <<i <<"] = " 
                     <<new_frag.mol.amber_bt_torsions[new_frag.mol.num_bonds-1][i] <<endl;
                cout <<"FTD: new_angle = " <<new_angle <<endl;
            }

            // Push the fragment onto the end of the frag_list vector
            frag_list.push_back(new_frag);

            // Then adjust the torsions of that fragment according to the angle at this iteration
            frag_list[frag_list.size()-1].mol.set_torsion(atom1, atom2, atom3, atom4, new_angle);
 
            // save the torsioned turned for isoswap
            frag_list[frag_list.size()-1].iso_tors_turned = new_angle;
        }
    }

    return;

} // end DN_GA_Build::frag_torsion_drive()



// +++++++++++++++++++++++++++++++++++++++++
// Calculate the internal energy of a fragment
void
DN_GA_Build::prepare_internal_energy( Fragment & frag1, Master_Score & score )
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

} // end DN_GA_Build::prepare_internal_energy()



// +++++++++++++++++++++++++++++++++++++++++
// Function to prune fragments by Hungarian RMSD using a generic best-first clustering
// algorithm
void
DN_GA_Build::prune_h_rmsd( vector <Fragment> & list_of_frags )
{ 
    cout << "Number of molecules before Hungarian RMSD pruning: " << list_of_frags.size() << endl;
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
    cout << "Number of molecules after Hungarian RMSD pruning: " << list_of_frags.size()<< endl;
    temp_vec.clear();

    return;

} // end DN_GA_Build::prune_h_rmsd()



// +++++++++++++++++++++++++++++++++++++++++
// Calculate the molecular weight of the dockmol object of a fragment
// (modified from amber_typer.cpp)
void
DN_GA_Build::calc_mol_wt( DOCKMol & mol )
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
                   <<" in DN_GA_Build::calc_mol_wt()" <<endl; }

    }

    // Assign the molecular weight to the mol object
    mol.mol_wt = mw;

    return;

} // end DN_GA_Build::calc_mol_wt()



// +++++++++++++++++++++++++++++++++++++++++
// Given a fragment, populate the rot_bonds field
void
DN_GA_Build::calc_rot_bonds( DOCKMol & mol )
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

} // end DN_GA_Build::calc_rot_bonds();



// +++++++++++++++++++++++++++++++++++++++++
// Given a fragment, populate the formal_charge field
void
DN_GA_Build::calc_formal_charge( DOCKMol & mol )
{
    float charge = 0.0;

    // Iterate over all atoms, find the partial charge
    for (int i=0; i<mol.num_atoms; i++){
        charge += mol.charges[i];
    }

    // Assign it directly to the referenced mol object
    mol.formal_charge = charge;

    return;

} // end DN_GA_Build::calc_formal_charge();

// +++++++++++++++++++++++++++++++++++++++++
// LEP: Given a mol, calculate the hydrogen acceptors and donor
void
DN_GA_Build::calc_num_HA_HD( DOCKMol & mol )
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

} //end DN_GA_Build::calc_num_HA_HD();


// +++++++++++++++++++++++++++++++++++++++++
// Activate all of the atoms AND bonds in a given DOCKMol
void
DN_GA_Build::activate_mol( DOCKMol & mol )
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

} // end DN_GA_Build::activate_mol()



// +++++++++++++++++++++++++++++++++++++++++
// Check to see if there are still unsatisfied attachment points in the fragment
bool
DN_GA_Build::dummy_in_mol( DOCKMol & mol )
{
    // For every atom, if one of them is a 'Du', return true
    for (int i=0; i<mol.num_atoms; i++){
        if (mol.atom_types[i].compare("Du") == 0){
            return true; 
        }
    }

    return false;

} // end DN_GA_Build::dummy_in_mol()



// +++++++++++++++++++++++++++++++++++++++++
// Given a specific dummy atom, change it to a hydrogen
void
DN_GA_Build::dummy_to_H( Fragment & frag, int heavy, int dummy )
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

} // end DN_GA_Build::dummy_to_H()



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Accessory, debugging, or unused functions                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by energy score
bool 
ga_fragment_sort( const Fragment & a, const Fragment & b )
{
    return (a.mol.current_score < b.mol.current_score);
}



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting FragGraph Tanimotos
bool 
ga_fgpair_sort( const pair<float,int> & a, const pair<float,int> & b )
{
    return (a.first > b.first);
}



// +++++++++++++++++++++++++++++++++++++++++
// Comparison function for sorting fragment vectors by number of heavy atoms
// (more heavy atoms comes first), then consider number of att points
bool 
ga_size_sort( const Fragment & frag1, const Fragment & frag2 )
{
    if (frag1.mol.heavy_atoms != frag2.mol.heavy_atoms)
        return (frag1.mol.heavy_atoms > frag2.mol.heavy_atoms);
    else
        return (frag1.aps.size() > frag2.aps.size());
}



// +++++++++++++++++++++++++++++++++++++++++
// Debugging function to print the TorEnv_GA data structure
void
DN_GA_Build::print_torenv( vector <TorEnv_GA> tmp_torenv_vector )
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

} // end DN_GA_Build::print_torenv()



// +++++++++++++++++++++++++++++++++++++++++
// Debugging function to print the FragGraph data structure
void
DN_GA_Build::print_fraggraph( )
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
        fout.open (fout_name.str().c_str(), fstream::out|fstream::app);
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
}

void DN_GA_Build::frag_sort(std::vector<Fragment> & vec_frag, std::function<bool(const Fragment&, const Fragment&)> func){

    std::list<Fragment> list_growing(vec_frag.begin(), vec_frag.end());
    list_growing.sort(func);

    vec_frag.clear();
    for (Fragment tmp_frag : list_growing){
        vec_frag.push_back(tmp_frag);
    }

}

void DN_GA_Build::frag_sort(std::vector< std::pair <Fragment, int> > & vec_frag,
                         std::function<bool(const std::pair <Fragment, int> &, const std::pair <Fragment, int> &)> func){

    std::list< pair < Fragment, int >> list_growing(vec_frag.begin(), vec_frag.end());
    list_growing.sort(func);
    for (pair<Fragment,int> tmp_frag : list_growing){
        vec_frag.push_back(tmp_frag);
    }

}


void
DN_GA_Build::print_molecules(std::string fout_molecules_name, std::vector <DOCKMol> & parents, int gen, std::string prefix){
    fstream fout_molecules;
    fout_molecules.open (fout_molecules_name, fstream::out|fstream::app);
    for (int j=0; j<parents.size(); j++){
        activate_mol(parents[j]);
    
        // calculate some ensemble properties
        calc_mol_wt( parents[j] );
        calc_rot_bonds( parents[j] );
        calc_formal_charge( parents[j] );
        calc_num_HA_HD( parents[j] );

        // Only Update name if a new molecule
        // BTB - if limit max change is on then the molecule is already named in PDM
        if (parents[j].parent == 0){
            naming_function(parents[j], gen+1, j, prefix);
        }

        // if you compartamaliz res nums, change the res id where
        // every time there is a change in subst id change the 
        // atom_res_num by incrementing by one. 
        // This is PURELY for Chimera visual analysis. Nothing else

        std::vector <std::string> l_atom_res_subst {};
        //l_atom_res_subst= parents[j].subst_names;
        for (int z=0; z < parents[j].num_atoms; z++){
            l_atom_res_subst.push_back(parents[j].subst_names[z]);
        }

        std::vector <std::string> new_l_atom_nums {};

        std::string string_holder = "";
        int subst_id = 0;
        for ( std::string res_n : l_atom_res_subst ){
            if ( string_holder == "" || string_holder != res_n){
                 subst_id++;
                 string_holder = res_n;
                 new_l_atom_nums.push_back( std::to_string( subst_id ) );
            } else{ new_l_atom_nums.push_back( std::to_string( subst_id ) ); }
        }
        for ( int z=0; z < parents[j].num_atoms; z++ ){
            parents[j].atom_residue_numbers[z] = new_l_atom_nums[z];
        }


        // sort and make the list of substance names unique 
        // to write out to the header
        std::sort(l_atom_res_subst.begin(), l_atom_res_subst.end());
        auto uniq_subt_name  =
            std::unique( l_atom_res_subst.begin(), l_atom_res_subst.end());
        l_atom_res_subst.erase(uniq_subt_name, l_atom_res_subst.end());

        std::stringstream frag_string_newname;
        for ( std::string str: l_atom_res_subst){
            frag_string_newname << str << "**";
        }

        fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Name:"
                       << setw(FLOAT_WIDTH) << parents[j].title <<endl;
        fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Generation:"
                       << setw(FLOAT_WIDTH) << gen+1 <<endl;
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

        fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Type:"
                       << setw(FLOAT_WIDTH) << parents[j].energy <<endl;

        fout_molecules << DELIMITER << setw(STRING_WIDTH) << "GA_Frag_String:"
                       << setw(FLOAT_WIDTH)<<" " + frag_string_newname.str() << std::endl;
        //current _data may not be initialized on write out
        fout_molecules << parents[j].current_data << endl;
        Write_Mol2(parents[j], fout_molecules);

        // Reset the atom_residue_nums back to an array of "1"
        for ( int z=0; z < parents[j].num_atoms; z++ ){
            parents[j].atom_residue_numbers[z] = "1";
        }
    }
    fout_molecules.close();

}


// ++++++++++++++++++++++++++++++++++++++++++++
// Create a unique name with the parents, generation number, and location in file
// If mol from all new, name is dn_gennum_loc
void
DN_GA_Build::naming_function( DOCKMol & mol, int gen, int loc, std::string prefix)
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

   //BTB 2022.10.24 - provide additional context for the molecule in the name
   new_title << prefix;
   
   // Add gen first  
   if ( gen < 10 ){ 
       new_title << dn_ga_name_identifier << "_g000" << gen ; 
   }
   else if ( gen < 100 ){
       new_title << dn_ga_name_identifier << "_g00" << gen ; 
   }
   else if ( gen < 1000 ){
       new_title << dn_ga_name_identifier << "_g0" << gen ; 
   }
   else {
       new_title << dn_ga_name_identifier << "_g" << gen ; 
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
   if (mol.energy != PARENT_TAG){
       mol.title = new_title.str();   
   }

   // Clear new_energy
   new_title.str("");
   list.clear();

   return;
}//end GA_Recomb::naming_function()
