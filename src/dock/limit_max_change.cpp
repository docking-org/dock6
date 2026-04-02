#include "limit_max_change.h"
#include "fingerprint.h"
#include "hungarian.h"
#include "trace.h"

// +++++++++++++++++++++++++++++++++++++++++
// static member initializers

// These are the same as those in Base_Score.
const string LimitMaxChange::DELIMITER    = "########## ";
const int    LimitMaxChange::FLOAT_WIDTH  = 20;
const int    LimitMaxChange::STRING_WIDTH = 17 + 19;

AMBER_TYPER LimitMaxChange::c_typer; //typer definition

// +++++++++++++++++++++++++++++++++++++++++
// Some constructors and destructors
LimitMaxChange::LimitMaxChange(){
    ramping_enabled = false;
    use_ramp = false;
    //Delta Max functions
    // (cutoff percent, threshold, name, invert comparisons, function address)
    limit_max_change_push_new_func_parameter(1, -2, "HMS", false, limit_max_change_HMS );
    limit_max_change_push_new_func_parameter(1, 40, "mol_wt", false, limit_max_change_mw  );
    limit_max_change_push_new_func_parameter(1, 0.7, "tan", true, limit_max_change_tanimoto );
    limit_max_change_push_new_func_parameter(1, 0.75, "VOS", true, limit_max_change_VOS );
    limit_max_change_push_new_func_parameter(1, 0.5, "FMS", true, limit_max_change_FMS );
}

LimitMaxChange::~LimitMaxChange(){
    limit_max_change_functions.clear();
    limit_max_change_possible_scoring_functions.clear();
    limit_max_change_possible_scoring_function_names.clear();
    limit_max_change_function_names.clear();
    limit_max_change_possible_thresholds.clear();
    limit_max_change_thresholds.clear();
    limit_max_change_possible_cutoffs.clear();
    limit_max_change_percent_match.clear();
    limit_max_change_possible_inversions.clear();
    limit_max_change_comparison_inversions.clear();
}

// +++++++++++++++++++++++++++++++++++++++++
// SETTERS & GETTERS
// file output prefix setter
void 
LimitMaxChange::set_output_prefix(std::string prefix){
    output_prefix = prefix;
}
bool
LimitMaxChange::get_ramping(){
    return ramping_enabled;
}

void LimitMaxChange::set_c_typer(AMBER_TYPER & typer){
    LimitMaxChange::c_typer = typer; //set the c_typer associated w/ delta max class equal to the parameter c_typer
        
    LimitMaxChange::c_typer.atom_model = typer.atom_model;
    LimitMaxChange::c_typer.vdw_defn_file =typer.vdw_defn_file;
    LimitMaxChange::c_typer.flex_defn_file = typer.flex_defn_file;
    LimitMaxChange::c_typer.flex_drive_tbl = typer.flex_drive_tbl;
    LimitMaxChange::c_typer.chem_defn_file = typer.chem_defn_file;
    LimitMaxChange::c_typer.ph4_defn_file = typer.ph4_defn_file;
    

    bool read_vdw = true;
    bool use_chem = true;
    bool use_ph4 = true;
    bool use_volume = true;
    bool read_gb_parm = false;
    //cout << "initialize c_typer " << endl;
    LimitMaxChange::c_typer.initialize(read_vdw, read_gb_parm,
                       use_chem, use_ph4, use_volume);
    
}

// +++++++++++++++++++++++++++++++++++++++++
// Read parameters from the dock.in file
void
LimitMaxChange::input_parameters( Parameter_Reader & parm )
{
    Trace trace( "Entering LimitMaxChange::input_parameters()" );
    cout << "\nLimit Max Change Parameters\n";
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    // BLOCK 1: INPUT FILES
    for (int k = 0; k < limit_max_change_possible_scoring_functions.size(); k++){
        std::string query = "delta_max_use_function_";
        query += limit_max_change_possible_scoring_function_names[k];
        if (parm.query_param(query,"yes", " no | yes") == "yes") {
            limit_max_change_functions.push_back(limit_max_change_possible_scoring_functions[k]);
            query = "limit_change_minimum_match_fraction_";
            query += limit_max_change_possible_scoring_function_names[k];
            limit_max_change_percent_match.push_back(atof(parm.query_param(query, to_string(limit_max_change_possible_cutoffs[k])).c_str()) );
            query = "limit_change_cutoff_";
            query += limit_max_change_possible_scoring_function_names[k];
            limit_max_change_thresholds.push_back(atof(parm.query_param(query, to_string(limit_max_change_possible_thresholds[k])).c_str()));
            limit_max_change_function_names.push_back(limit_max_change_possible_scoring_function_names[k]);
            limit_max_change_comparison_inversions.push_back(limit_max_change_possible_inversions[k]);
        }
    }

    ramping_enabled = (parm.query_param("limit_change_enable_ramping","yes", " no | yes") == "yes"); //determines if ramping is enabled
    limit_max_change_ceiling = (parm.query_param("limit_change_ceiling","yes", " no | yes") == "yes"); //determines if the percent of parents to matched is cieling or floored
    limit_max_change_possible_scoring_functions.clear();
    limit_max_change_possible_cutoffs.clear();
    limit_max_change_possible_thresholds.clear();
    limit_max_change_possible_scoring_function_names.clear();
    limit_max_change_possible_inversions.clear();
    return;
}  //end LimitMaxChange::input_parameters




// +++++++++++++++++++++++++++++++++++++++++++++++
// Function to initialize a delta_max function
void 
LimitMaxChange::limit_max_change_push_new_func_parameter(float cutoff, 
		float threshold, 
		string name, 
		bool invert_comparison, 
		std::function <float(DOCKMol &, DOCKMol &)> func){
    Trace trace( "Entering LimitMaxChange::limit_max_change_push_new_func_parameter()" );
    using namespace std::placeholders; 
    limit_max_change_possible_scoring_functions.push_back(func);
    limit_max_change_possible_scoring_function_names.push_back(name);
    limit_max_change_possible_thresholds.push_back(threshold);
    limit_max_change_possible_cutoffs.push_back(cutoff);
    limit_max_change_possible_inversions.push_back(invert_comparison);

    //ramping stuff
    limit_max_change_ramp_cutoff.push_back(threshold);
    limit_max_change_ramp_percent_match.push_back(cutoff); //currently unused
} //end LimitMaxChange::limit_max_change_push_new_func_parameter

// +++++++++++++++++++++++++++++++++++++++++
// ALL POSSIBLE delta_max FUNCTIONS
//takes MW difference of the two moles and absolute values
float
LimitMaxChange::limit_max_change_mw(DOCKMol & mol, DOCKMol & MW_ref_mol){
    Trace trace( "Entering LimitMaxChange::limit_max_change_mw()" );
    int local_verbosity = 1;
    float MW_diff = mol.mol_wt - MW_ref_mol.mol_wt;
    if (MW_diff < 0) {
        MW_diff = -MW_diff;
    }
    return MW_diff;
}

// +++++++++++++++++++++++++++++++++++++++++++++++
// tanimoto similarity
float 
LimitMaxChange::limit_max_change_tanimoto(DOCKMol & mol, DOCKMol & tan_ref_mol){
    Trace trace( "Entering LimitMaxChange::limit_max_change_tanimoto" );
    int local_verbosity = 1;
    Fingerprint finger;
    return finger.compute_tanimoto(mol, tan_ref_mol);
}
// +++++++++++++++++++++++++++++++++++++++++++++++


// +++++++++++++++++++++++++++++++++++++++++++++++
//Volume overlap Score 
float
LimitMaxChange::limit_max_change_VOS(DOCKMol & mol, DOCKMol & vos_ref_mol){
    Trace trace( "Entering LimitMaxChange::limit_max_change_VOS" );
    int local_verbosity = 1;
    Volume_Score VOS;
    
    float old_score = mol.current_score;
    VOS.compute_score(mol, vos_ref_mol, c_typer);
    float vol_score = mol.current_score;
    cout << VOS.output_score_summary(mol) << endl;
    mol.current_score = old_score;
    return vol_score;
}
// +++++++++++++++++++++++++++++++++++++++++++++++

// +++++++++++++++++++++++++++++++++++++++++++++++
//Pharmcophore Matching Similarity
float
LimitMaxChange::limit_max_change_FMS(DOCKMol & mol, DOCKMol & FMS_ref_mol){
    Trace trace( "Entering LimitMaxChange::limit_max_change_FMS" );
    int local_verbosity = 1;
    Ph4_Score FMS;
    float old_score = mol.current_score;
    std::string ph4_compare_type = "overlap";
    FMS.compute_score(mol, FMS_ref_mol, ph4_compare_type, c_typer);
    //return the score
    float FMS_score = mol.current_score;
    cout << FMS.output_score_summary(mol) << endl;
    mol.current_score = old_score;
    return FMS_score;
}
// +++++++++++++++++++++++++++++++++++++++++++++++

// +++++++++++++++++++++++++++++++++++++++++++++++
//Footprint Similarity Score 
float
LimitMaxChange::limit_max_change_FPS(DOCKMol & mol, DOCKMol & FPS_ref_mol){
    //Trace trace( "Entering LimitMaxChange::limit_max_change_FPS" );
    //int local_verbosity = 1;
    //Footprint_Similarity_Score FPS;
    //FPS.receptor = receptor;
    //FPS.fps_foot_compare_type = "Pearson";
    //FPS.volume_ref_mol = DOCKMol(FPS_ref_mol);

    //DOCKMol mol_cop = DOCKMol(mol);
    //FPS.compute_score(mol_cop);
    //float FPS_score = mol_cop.current_score;
    //return FPS_score;
    return 0;
}
// +++++++++++++++++++++++++++++++++++++++++++++++

// +++++++++++++++++++++++++++++++++++++++++++++++
// Hungarian Matching Similarity Score
float
LimitMaxChange::limit_max_change_HMS( DOCKMol & mol, DOCKMol & hun_ref_mol){
    Trace trace( "Entering LimitMaxChange::limit_max_change_HMS()" );
    //init local vars
    int local_verbosity = 1;
    Hungarian_RMSD hms;
    //init var at 0
    int hun_ref_heavy_atoms=0;
    std::pair <double, int>  hms_num;
    float temp_hun_score;
    int name_assigned = 0;

    limit_max_change_activate_mol(mol);
    limit_max_change_activate_mol(hun_ref_mol);

    hms_num = hms.calc_Hungarian_RMSD_dissimilar(mol, hun_ref_mol);
    //loop through all atoms to count heavy atoms and increment var
    for (int k=0; k<hun_ref_mol.num_atoms; k++){
        if (hun_ref_mol.atom_types[k].compare("H") != 0 && hun_ref_mol.atom_types[k].compare("Du") != 0){
            hun_ref_heavy_atoms++;
        }
    }

    // calc rmsd of matching atoms
    temp_hun_score = (-5.0) * (hun_ref_heavy_atoms - hms_num.second) / hun_ref_heavy_atoms + hms_num.first;
    
    if (local_verbosity){
        cout << " number of ref hvy atoms: " << hun_ref_heavy_atoms << " num unmatched atoms total: " << hms_num.second << " num of lig hvy atoms: " << mol.heavy_atoms << " matched rmsd: " << hms_num.first << endl; 
    }

    limit_max_change_deactivate_mol(mol);
    limit_max_change_deactivate_mol(hun_ref_mol);

    return temp_hun_score;
}//end LimitMaxChange::limit_max_change_HMS()
// +++++++++++++++++++++++++++++++++++++++++++++++


// +++++++++++++++++++++++++++++++++++++++++++++++
// DEBUGGING FUNCTION
//this function is for dev purposes only. It is meant to help with debugging of the prune_divergent_molecules function
void LimitMaxChange::print_PDM_comparisons(vector<float> & limit_max_change_scores, std::vector <DOCKMol> & parents, DOCKMol & child, int c_index, int score_index, std::vector <bool> limit_max_change_evals, int min_parents){
    Trace trace( "Entering LimitMaxChange::print_PDM_comparisons()" );
    int using_title = 0;
    std::string  score_name = limit_max_change_function_names[score_index];
    std::string child_name;
    if (using_title){
         child_name = child.title;
    }
    else{
        stringstream ss;
        ss << "child index: "<< c_index << ":";
        child_name = ss.str();
    }
    int matched = 0;
    int total = limit_max_change_evals.size();
    for (int index = 0; index < total; index++){
        if (limit_max_change_evals[index]) {
            matched += 1;  
        }
    }
    cout << "Child " << c_index << ": " << endl; 
    cout << matched << " out of " << total << " parents with a " << score_name << " cutoff of " << limit_max_change_thresholds[score_index] << endl;
    std::string result = "";
    if (matched < min_parents){
        result = "rejected";
    }
    else{
        result = "accepted";
    }
    cout << matched << " is >= " << min_parents << " required, child " << result << endl;

/*    if (matched < min_parents){
        cout << child_name << " REJECTED on " << score_name << endl;
    }
    else{
        cout << child_name << " SURVIVED on " << score_name << endl;
    }
*/
    // Output for scores greater than or equal to the threshold
    cout << score_name << " matched: ";
    for (int i = 0; i < parents.size(); i++){
        if (limit_max_change_evals[i]) {
            cout << "[" << parents[i].title << ", " << limit_max_change_scores[i] << "] ";
        }
    }

    // Output for scores less than the threshold
    cout << endl << score_name << " unmatched" << ": ";
    for (int i = 0; i < parents.size(); i++){
        if (!limit_max_change_evals[i]) {
            cout << "[" << parents[i].title << ", " << limit_max_change_scores[i] << "] ";
        }
    }

    cout << endl;
}//end LimitMaxChange::print_PDM_comparisons()
// +++++++++++++++++++++++++++++++++++++++++++++++

// +++++++++++++++++++++++++++++++++++++++++++++++
// PUBLIC FUNCTION TO CALL MAIN FUNCTION IN OTHER MODULES
void 
LimitMaxChange::prune_divergent_molecules(std::vector <DOCKMol> & scored_generation, std::vector <DOCKMol> & divergent_children, std::vector <DOCKMol> & old_parents, int gen_num){
      prune_divergent_molecules_main(scored_generation, 
		      divergent_children, 
		      old_parents, 
		      limit_max_change_functions, 
		      limit_max_change_percent_match, 
		      limit_max_change_thresholds, 
		      limit_max_change_function_names, 
		      limit_max_change_comparison_inversions, 
		      gen_num);
}

// +++++++++++++++++++++++++++++++++++++++++++++++
// MAIN LIMIT_MAX_CHANGE_FUNCTION
// this function reads in a vector of children, compares them to a vector of parents, and then modifies the reference of children by removing children that fell too far from the hungarian tree. These children that fell too far will be stored in pruned_children for future print out if desired. 
void
LimitMaxChange::prune_divergent_molecules_main(
		std::vector <DOCKMol> & children,
	       	std::vector <DOCKMol> & pruned_children, 
		std::vector <DOCKMol> & parents, 
		std::vector <std::function <float (DOCKMol &, DOCKMol &)>> limit_max_change_score_functions,
	       	std::vector <float> & percent_match, 
		std::vector <float> & thresholds, 
		std::vector <std::string> & score_names, 
		std::vector <bool> & conversions, 
		int generation){

    Trace trace( "Entering LimitMaxChange::prune_divergent_molecules_main()" );
    //local variables for dev purposes
    int local_verbosity = 2;

    //parameters and calculated variables
    vector <int> cutoff;

    float match_val;
    //casting is necessary 
    for (int k = 0; k < limit_max_change_score_functions.size(); k++){
	if (limit_max_change_ceiling) {
	    match_val = ceil(int(percent_match[k] * float(parents.size())) );
	}
	else {
	    match_val = floor(int(percent_match[k] * float(parents.size())) );
	}
        cutoff.push_back(match_val);
    
        if (percent_match[k] == 0){
            cutoff[k] = 0;
        }
    }

    //loop control variables
    int i = 0;
    int count = 0;
    vector <int> passed; //stores the number of passes a child gets for each score type
    vector <vector <bool> > eval; //stores the individual parent/child match values, will store a false if the child deviates too much from the parent, true if it matches; stores for each score type
    bool mol_passed; //this is set to true if the mol passess delta max
    float cur_limit_max_change_score; 

    //variables to store properties 
    vector <vector <float> >  limit_max_change_scores; // stores the comparisons from delta max calculations
    vector <vector <vector <float> > > limit_max_change_scores_passed_mols; // stores the comparisons of passed moles
    vector <vector <vector <float> > > limit_max_change_scores_failed_mols; // stores the comparisons of failed moles
    //variables to store/calculate ramping values
    vector <vector <float> > super_greatest_miss_values; //stores collections of worst scores for rejected children
    vector <float> greatest_miss_values; //stores the furthest value from the threshold (what the threshold would have to be for the mol to have passed) for each failed mol (inner vector) for each score (outer vector)

    //initialize vectors
    // need to double check whether this is necessary
    for (int k = 0; k<limit_max_change_score_functions.size(); k++){
        vector <float> filler;
	vector <bool> eval_filler;
        limit_max_change_scores.push_back(filler);
        passed.push_back(0);
	eval.push_back(eval_filler);
    }

    if (local_verbosity >= 2){
        cout << endl << DELIMITER << " Entering main PDM (Prune Divergent Molecules/Delta Max) function for generation " << generation << " " << DELIMITER << endl; 
	for (int k = 0; k<limit_max_change_score_functions.size(); k++){
            cout << "Score type: " << score_names[k] << endl; 
	    cout << "Fraction Parents Required to pass: " << percent_match[k] << endl;
	    cout << "Score to match a parent: " << thresholds[k] << " " << score_names[k] << endl;
	    cout << "Number of parents in population: " << parents.size() << endl;
	    float unadjusted_val = percent_match[k] * float(parents.size());
	    cout << "Number of parents required to match: " << unadjusted_val << "(" << cutoff[k] << " when rounded)";
	    cout << " based on " << parents.size() << " * " << percent_match[k] << " parents "<< endl; 
        }  
    }
    if (local_verbosity >= 4){
        for (int i = 0; i < children.size(); i++){
            cout << "child info:" << "[child index: " << i << ",child type: " << children[i].energy << "]" << endl;
        }

        for (int j = 0; j < parents.size(); j++){
            cout << "parent info:" << "[parent title: " << parents[j].title << ",parent type: " << parents[j].energy << "]" << endl;
        }
    }

    //beginning of main PDM loom
    if (local_verbosity >= 4){
        cout << "debug function vec size: "<< limit_max_change_score_functions.size() << endl;
        cout << "debug children size: "<< children.size() << endl;
        cout << "debug parent size: "<< parents.size() << endl << endl;
    }
    if (local_verbosity >= 2){
       cout << "Printing children to parent PDM comparisons: " << endl << endl;
    }
    while (i < children.size()){
	
        mol_passed = true;
	
        for(int k = 0; k<limit_max_change_score_functions.size(); k++){
             passed[k] = 0;
	     eval[k].clear();
	     for (int j = 0; j < parents.size(); j++){
	         eval[k].push_back(false);
	     }

        }


	for (int k=0; k<limit_max_change_score_functions.size(); k++){
            for( int j = 0; j<parents.size(); j++){
                auto limit_max_change_func = limit_max_change_score_functions[k];
		if (local_verbosity >= 4){
		    cout << "i: " << i << " k : " << k << " j: " << j << endl;
		}
		cur_limit_max_change_score = limit_max_change_func(children[i], parents[j]);
                if ( (cur_limit_max_change_score < thresholds[k] && !conversions[k]) || (cur_limit_max_change_score > thresholds[k] && conversions[k]) ){
                    passed[k] += 1;
		    eval[k][j] = true;
                }


                limit_max_change_scores[k].push_back(cur_limit_max_change_score);
            }
        }

        for(int k = 0; k<limit_max_change_score_functions.size(); k++){
            if (passed[k] < cutoff[k]){
                mol_passed = false;
		float RampElement;
		if (conversions[k]){
		    RampElement = *std::max_element(std::begin(limit_max_change_scores[k]), std::end(limit_max_change_scores[k]) ); //get minimum value 
		}
		else{
		    RampElement = *std::min_element(std::begin(limit_max_change_scores[k]), std::end(limit_max_change_scores[k]) ); //get maximum value 
		}
                
		//push the threshold back the molecule would've passed with
		if (local_verbosity >= 6){
		    cout << "ramp value: " << RampElement << endl;
		}
		greatest_miss_values.push_back(RampElement);
                //break;
            }
	    else{
		//if the molecule passed with the default threshold pushback the default threshold
		if (local_verbosity >= 6){
		    cout << "default threshold: " << thresholds[k] << endl;
		}
	        greatest_miss_values.push_back(thresholds[k]);
	    }
        }

	//temporarily name the child
        std::string child_name;
        stringstream ss;
        ss << "child index: " << count;
        child_name = ss.str();
        
        if (local_verbosity > 0){
            std::string result_text;
	    cout << DELIMITER << endl;
	    if (children[i].parent && local_verbosity >= 4){
                cout << "index " << i << " is a (mutated) parent molecule " << endl;
            } 
            for(int k = 0; k<limit_max_change_score_functions.size(); k++){
                    print_PDM_comparisons(limit_max_change_scores[k], parents, children[i], count, k, eval[k], cutoff[k]);
            }
	    //print out the overall result
            if (!mol_passed){
                result_text = "REJECTED ";
                }
            else{
                result_text = "ACCEPTED ";
            }
            cout << endl << "PDM() overall " << result_text << child_name << endl << DELIMITER << endl << endl;
        }
        //push child to appropriate vector if it didnt pass PDM
        if (!mol_passed){
	    limit_max_change_scores_failed_mols.push_back(limit_max_change_scores);
            pruned_children.push_back(children[i]);
            children.erase(children.begin()+i);
        }
        //increment i if the child passed
        else{
            i++;
        }
        count++;

	//save scores for passed child
        if (mol_passed){
            limit_max_change_scores_passed_mols.push_back(limit_max_change_scores); 
        }
	else{
	    super_greatest_miss_values.push_back(greatest_miss_values);
	}
	for (int k = 0; k < limit_max_change_scores.size(); k++){
	    limit_max_change_scores[k].clear();
	}
        greatest_miss_values.clear();	

    }//END WHILE LOOP


    //IDENTIFY THE LEAST AMOUNT OF RAMPING REQUIRED TO ALLOW A FAILED MOLECULE TO PASS
    //calculate sum square of the difference of the failed scores and the thresholds (if it passed for a score this will be zero), and then sum these for each score to come up with a molecule score, smaller scores are better
    float min_molecule_score = 99999999; //need to change this later to avoid using a hardcoded intial value
    int min_molecule_index = 0;
    for (int i = 0; i < super_greatest_miss_values.size(); i++){
        float molecule_score = 0;
	//index variable k is used here because we are iterating the different delta max functions used
	for (int k = 0; k < super_greatest_miss_values[i].size(); k++){
	    molecule_score += pow(super_greatest_miss_values[i][k] -  thresholds[k], 2);
	}
	if (molecule_score < min_molecule_score){
	    min_molecule_score = molecule_score;
	    min_molecule_index = i;
	}
    }



    //take the molecule with the smallest molecule score and use it to determine ramping values for the next generation, should they be needed/enabled
    if (super_greatest_miss_values.size() > 0){
        for (int k = 0; k < super_greatest_miss_values[min_molecule_index].size(); k++){
             limit_max_change_ramp_cutoff[k] = super_greatest_miss_values[min_molecule_index][k]; 
        }
    }

    //just readd the closest molecule to passing to the scored generation when ramping is enabled
    if (ramping_enabled && children.size() == 0 && pruned_children.size() > 0){
	cout << "ramping molecule" << endl;
        children.push_back(pruned_children[min_molecule_index]);
	limit_max_change_scores_passed_mols.push_back(limit_max_change_scores_failed_mols[min_molecule_index]);
    }


    //Final IO && deallocation
    //deallocate vectors
    for (int i = 0; i < super_greatest_miss_values.size(); i++){
        super_greatest_miss_values[i].clear();
    }

    //WRITE OUT ALL SCORES FOR PASSED MOLECULE
    // generate file name and stream
    std::ostringstream PDM_fout_molecules_name;
    std::ostringstream gen;
    gen << std::setw(4) << std::setfill('0') << generation+1;
    PDM_fout_molecules_name << output_prefix << ".delta_max" << gen.str() <<".dat";
   
    cout << "foutmolecules name: " << PDM_fout_molecules_name.str() << endl;
    std::fstream PDM_fout_molecules;

    PDM_fout_molecules.open(PDM_fout_molecules_name.str(), fstream::out);



    //Write out delta max scores for passed molecules
    int compare_value = children.size();
    for (int i=0; i<children.size(); i++){
         for (int k = 0; k< limit_max_change_scores_passed_mols[i].size(); k++){
         PDM_fout_molecules << DELIMITER << setw(STRING_WIDTH) << "Molecule:" << setw(FLOAT_WIDTH) << gen.str() << "_i" << i+1 << endl;
             for (int j=0; j<limit_max_change_scores_passed_mols[i][k].size(); j++){
                  PDM_fout_molecules << DELIMITER << setw(STRING_WIDTH) << "limit_max_change " << score_names[k] << " comparison " << j << " : " << setw(FLOAT_WIDTH) << limit_max_change_scores_passed_mols[i][k][j] <<endl;
             }
             limit_max_change_scores_passed_mols[i][k].clear();

         }
         limit_max_change_scores_passed_mols[i].clear();
    }
 
    passed.clear();
    eval.clear(); //this is a vector of vectors, need to ensure that the inner vectors are cleared and not creating mem leaks
    limit_max_change_scores.clear();
    limit_max_change_scores_passed_mols.clear();
    limit_max_change_scores_failed_mols.clear();
    super_greatest_miss_values.clear();
    PDM_fout_molecules.close(); 
    PDM_fout_molecules_name.str(std::string());
    gen.str(std::string());
    cout << DELIMITER << " Exiting PDM/Delta Max " << DELIMITER << endl << endl << endl << "-----------------------------------" << endl;
    

}//end LimitMaxChange::prune_divergent_molecules()
// +++++++++++++++++++++++++++++++++++++++++++++++

// +++++++++++++++++++++++++++++++++++++++++++++++
// SOME UTILITY FUNCTIONS FOR ACTIVATING/DEACTIVATING MOLS
void
LimitMaxChange::limit_max_change_deactivate_mol( DOCKMol & mol )
{
   Trace trace( "Entering LimitMaxChange::deactivate_mol()" );
   // Iterate through all atoms and set atom_active_flag to true
   for ( int i=0; i<mol.num_atoms; i++ ){
       mol.atom_active_flags[i] = false;
   }

   // Iterate through all bonds and set bond_active_flag to true
   for ( int i=0; i<mol.num_bonds; i++ ){
       mol.bond_active_flags[i] = false;
   }

   return;

} // end LimitMaxChange::deactivate_mol()




// +++++++++++++++++++++++++++++++++++++++++
// Activate all of the atoms AND bonds in a given DOCKMol - called in sample_minimized_torsions
void
LimitMaxChange::limit_max_change_activate_mol( DOCKMol & mol )
{
   Trace trace( "Entering limit_max_change::activate_mol()" );
   // Iterate through all atoms and set atom_active_flag to true
   for ( int i=0; i<mol.num_atoms; i++ ){
       mol.atom_active_flags[i] = true;
   }

   // Iterate through all bonds and set bond_active_flag to true
   for ( int i=0; i<mol.num_bonds; i++ ){
       mol.bond_active_flags[i] = true;
   }

   return;

} // end LimitMaxChange::activate_mol()
 
// +++++++++++++++++++++++++++++++++++++++++++++++ 
// END OF FILE
