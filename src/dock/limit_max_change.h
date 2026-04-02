#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <functional> //for limit_max_change

#include "amber_typer.h"
#include "dockmol.h"
#include "master_score.h"
#include "utils.h"

#ifdef BUILD_DOCK_WITH_RDKIT
#include "rdtyper.h"
class RDTYPER;
#endif

using namespace std;

// ++++++++++++++++++++++++++++++++++++
// CLASS DECLARATION
class LimitMaxChange {
    public:    
	LimitMaxChange(); //constructor
	~LimitMaxChange(); //destructor

        // ++++++++++++++++++++++++++++++++++++
        // PUBLIC MEMBER VARIABLES
        // use these for consistent formatting in output_score_summary
        static const std::string DELIMITER;
        static const int         FLOAT_WIDTH;
        static const int         STRING_WIDTH;                         // standard output header

        // ++++++++++++++++++++++++++++++++++++
        // PUBLIC MEMBER FUNCTIONS
        void                     prune_divergent_molecules(
                                       std::vector <DOCKMol> & children,
                                       std::vector <DOCKMol> & pruned_children,
                                       std::vector <DOCKMol> & parents,
                                       int generations); // A public version to be called from inside other .cpp files to perform limit_max_change/delta_max
        void                     input_parameters( Parameter_Reader & parm );                          // initialize limit_max_change input parameters                                                                                                                                                                                                                                         
	// ++++++++++++++++++++++++++++++++++++
	// SETTERS
	static  void                    set_c_typer(AMBER_TYPER & c_typer);
	void                     set_output_prefix(std::string prefix);


	// GETTERS
	bool                     get_ramping();

    private: 
        // ++++++++++++++++++++++++++++++++++++
        // MEMBER VARIABLES
	static AMBER_TYPER              c_typer; //typer, provided from calling superclass, needed for FMS 
	std::string              output_prefix = "";         // prefix for file writouts

        vector <float>           limit_max_change_cutoff;             // BTB - HMS cutoff for max change
        vector <float>           limit_max_change_percent_match;            // BTB - what % of parents a child must match to be retained in "prune_divergent_molecules" function 
        vector <float>           limit_max_change_thresholds;               // BTB - threshold for each max change function (i.e. -2, -1, or  0 for hms)     
        vector <bool>            limit_max_change_comparison_inversions;
        vector <float>           limit_max_change_possible_thresholds;
        vector <float>           limit_max_change_possible_cutoffs;         vector <bool>            limit_max_change_possible_inversions;
        vector <std::function    <float(DOCKMol &, DOCKMol &)> >           limit_max_change_functions;                // BTB - which functions used for max change
        vector <std::string>     limit_max_change_function_names;           // BTB - names of max change functions
        vector <std::function    <float(DOCKMol &, DOCKMol &)> >           limit_max_change_possible_scoring_functions;  // BTB - a comprehensive list of all potential delta max functions, important for the input file decision tree
        vector <std::string>     limit_max_change_possible_scoring_function_names; // BTB - corresponding names of each function in the input file decision tree
        bool                     limit_max_change_ceiling; //whether to floor/cieling the number of parents required to match for dmax

        //RAMPING VARIABLES
	bool                     ramping_enabled; // if we should ever ramp
	bool                     use_ramp; // if ramp values should be used for the current generation
	vector <float>           limit_max_change_ramp_cutoff; // the new threshold to use when ramping is on
        vector <float>           limit_max_change_ramp_percent_match; // The new percentage of parents needed to meet when ramping is on

        // ++++++++++++++++++++++++++++++++++++
        // PRIVATE MEMBER Functions
        void                    limit_max_change_push_new_func_parameter(
                                      float cutoff,
                                      float threshold,
                                      string name,
                                      bool invert_comparison,
                                      std::function <float(DOCKMol &, DOCKMol &)> func);              // initialize input parameters and defaults for delta max functions
        
        void                   prune_divergent_molecules_main(
                                      std::vector <DOCKMol> & children,
                                      std::vector <DOCKMol> & pruned_children,
                                      std::vector <DOCKMol> & parents,
                                      std::vector <std::function <float (DOCKMol &, DOCKMol &)>> limit_max_change_score_functions,
                                      std::vector <float> & percent_match,
                                      std::vector <float> & thresholds,
                                      std::vector <std::string> & score_names,
                                      std::vector <bool> & conversions,
                                      int generation);                                                // limit_max_change main body, this function calls each delta max function specified in the input file and uses it to identify molecules which deviate too far from another set of molecules (i.e. children to parents in the GA)
        
        
        void                    print_PDM_comparisons(vector<float> & limit_max_change_scores, std::vector <DOCKMol> & parents, DOCKMol & child, int c_index, int score_index, std::vector <bool> limit_max_change_evals, int min_num_parents); //BTB - for debugging
        
        // ++++++++++++++++++++++++++++++++++++
        //DELTA MAX FUNCTIONS
                                        //child is reference, parent is mol
        static  float                   limit_max_change_HMS( DOCKMol & parent, DOCKMol & child); //BTB this needs to be static to push back onto a vector, HMS for limit_max_change
        static  float                   limit_max_change_VOS(DOCKMol & mol, DOCKMol & vos_ref_mol); // BTB - VOS
        static  float                   limit_max_change_mw(DOCKMol & mol, DOCKMol & MW_ref_mol); // Compare via molecular weight
        static  float                   limit_max_change_tanimoto(DOCKMol & mol, DOCKMol & tan_ref_mol); // Compare via fingerprint similarity
	static  float                   limit_max_change_FPS(DOCKMol&, DOCKMol&); // Compare via FPA (NYI)
	static  float                   limit_max_change_FMS(DOCKMol&, DOCKMol&); // Compare via Pharmacophore Score (NYI)

        //DELTA MAX FUNCTION SUBROUTINES (i.e. performing initialization that would normally be done during input_parameter collection and using the chem defn/ph4 defn files)
            
            
        // ++++++++++++++++++++++++++++++++++++
        //UTILITY FUNCTIONS
        static  void                   limit_max_change_deactivate_mol ( DOCKMol &);                                          // turn off active atom/bond flags
        static  void                   limit_max_change_activate_mol ( DOCKMol &);                                            // activate atom/bond flags

};
