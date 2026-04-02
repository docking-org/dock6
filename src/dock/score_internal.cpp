#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <sstream>
#include "amber_typer.h"
#include "base_score.h"
#include "dockmol.h"
#include "grid.h" 
//#include "score_descriptor.h"
#include "score_internal.h"
#include "utils.h"

#include <math.h>

using namespace std;

// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++

// +++++++++++++++++++++++++++++++++++++++++

Internal_Energy_Score::Internal_Energy_Score()
    :   att_exp( 0 ),
        rep_exp( 0 )
{
}

// +++++++++++++++++++++++++++++++++++++++++
Internal_Energy_Score::~Internal_Energy_Score()
{
}

void
Internal_Energy_Score::close()
{
    // clear contance of c++ strings 
    //receptor_filename.clear();
}

// +++++++++++++++++++++++++++++++++++++++++
void
Internal_Energy_Score::input_parameters(Parameter_Reader & parm,
                                          bool & primary_score,
                                          bool & secondary_score)
{
    string          tmp;
   // string          tmp2;

    use_primary_score = false;
    use_secondary_score = false;


    // intailize parameters.

    cout << "\nInternal Energy Score Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    if (!primary_score) {
        tmp = parm.query_param("internal_energy_score_primary", "no", "yes no");
        if (tmp == "yes")
            use_primary_score = true;
        else
            use_primary_score = false;

        primary_score = use_primary_score;
    }

    if (!secondary_score) {
        //tmp = parm.query_param("hbond_score_secondary", "no", "yes no");
        //if (tmp == "yes")
            //use_secondary_score = true;
        //else
            use_secondary_score = false;

        secondary_score = use_secondary_score;
    }

    ie_att_exp = 6;
    ie_diel = 4.0;
    if (use_primary_score || use_secondary_score){
        ie_rep_exp = atoi(parm.query_param("internal_energy_rep_exp", "12").c_str());
        use_score = true;
    }
    else
        use_score = false;
}

// +++++++++++++++++++++++++++++++++++++++++
void
Internal_Energy_Score::initialize(AMBER_TYPER & typer)
{

    if (use_score) {
        //use_internal_energy = true;

        init_vdw_energy(typer, att_exp, rep_exp);
        //initialize_internal_energy(

    }
}


// +++++++++++++++++++++++++++++++++++++++++
bool
Internal_Energy_Score::compute_score(DOCKMol & mol)
{
      bool temp_int = use_internal_energy;
      use_internal_energy = true; // it must be true to calculate internal energy
      mol.current_score = compute_ligand_internal_energy(mol); 
      mol.current_data = output_score_summary(mol);

      use_internal_energy = temp_int; // change it back

      return true;
}

// +++++++++++++++++++++++++++++++++++++++++
string
Internal_Energy_Score::output_score_summary(DOCKMol & mol)
{
    ostringstream text;
//    cout << "inside output_score_summary" << endl;
    if (use_score) {
        text << DELIMITER << setw(STRING_WIDTH) << "Internal_Score:"
             << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl;
    }
    return text.str();
}


// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++
