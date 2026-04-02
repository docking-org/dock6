#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <sstream>
#include "amber_typer.h"
#include "base_score.h"
#include "dockmol.h"
#include "grid.h" 
#include "score.h"
#include "utils.h"

using namespace std;


// +++++++++++++++++++++++++++++++++++++++++
Bump_Filter::Bump_Filter()
{
    bump_grid = NULL;
}

// +++++++++++++++++++++++++++++++++++++++++
Bump_Filter::~Bump_Filter()
{
    delete bump_grid;
}

// +++++++++++++++++++++++++++++++++++++++++
//collect user defined parameters
void
Bump_Filter::input_parameters(Parameter_Reader & parm)
{
    cout << "\nBump Filter Parameters\n"
         << "--------------------------------------------------------"
            "----------------------------------"
         << endl;

    bump_filter = parm.query_param("bump_filter", "no", "yes no") == "yes";

    if (bump_filter) {
        
        grid_file_name = parm.query_param( "bump_grid_prefix", "grid" );

        anchor_bump_max = atoi(parm.query_param("max_bumps_anchor", "2").c_str());
        if (anchor_bump_max <= 0) {
            cout << "ERROR:  Parameter must be integer greater than zero."
                    "  Program will terminate."
                 << endl;
            exit(0);
        }

        growth_bump_max = atoi(parm.query_param("max_bumps_growth", "2").c_str());
        if (growth_bump_max <= 0) {
            cout << "ERROR:  Parameter must be integer greater than zero."
                    "  Program will terminate."
                 << endl;
            exit(0);
        }
    }
}

// +++++++++++++++++++++++++++++++++++++++++
//call to read in data from pregenerated Bump Grid
//and calculate location of grid in space
void
Bump_Filter::initialize()
{
    if (bump_filter) {
        cout << "Initializing Bump Filter Routines..." << endl;
        // get an instance of the bump grid merely to read the bump grid file.
        bump_grid = new Bump_Grid(); 
        bump_grid->get_instance(grid_file_name);
        use_score = true;
    } else
        use_score = false;
}

// +++++++++++++++++++++++++++++++++++++++++
//call to calculate number of bumps for anchor orientation
//if number of bumps exceeds user defined cutoff, do
//not allow anchor to continue
bool
Bump_Filter::check_anchor_bumps(DOCKMol & mol, bool more_orients)
{
    if(!more_orients)
        return true;

    if (bump_filter) {
        int num_bumps = get_bump_score(mol);
        return (num_bumps <= anchor_bump_max) && (num_bumps != -1);
    } else
        return true;
}
 
// +++++++++++++++++++++++++++++++++++++++++
//call to calculate number of bumps for grown conformation
//if number of bumps exceeds user defined cutoff, do
//not allow conformation to continue
bool
Bump_Filter::check_growth_bumps(DOCKMol & mol)
{
    if (bump_filter) {
        int num_bumps = get_bump_score(mol);
        return (num_bumps <= growth_bump_max) && (num_bumps != -1);
    } else
        return true;
}

// +++++++++++++++++++++++++++++++++++++++++
//calculate number of bumps
int
Bump_Filter::get_bump_score(DOCKMol & mol)
{
    int             bump_count = 0;
    int             vdw_threshold = 0;

    for (int atom = 0; atom < mol.num_atoms; atom++) {
        if (mol.atom_active_flags[atom] && mol.amber_at_heavy_flag[atom]) {
            if (bump_grid->is_inside_grid_box(mol.x[atom], mol.y[atom], mol.z[atom])) {
                bump_grid->find_grid_neighbors(mol.x[atom], mol.y[atom], mol.z[atom]);
                vdw_threshold = NINT(10 * mol.amber_at_radius[atom]);

                if (vdw_threshold >= (int) bump_grid->bump[bump_grid->nearest_neighbor])
                    bump_count++;
            } else
                return -1;  //atom is outside the grid box
        }
    }

    return (bump_count);
}


// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++
Energy_Score::Energy_Score()
{
    energy_grid = NULL;
    grid_lig_efficiency = false; //this gets called when it is not initialized, so this will prevent any issues
}

// +++++++++++++++++++++++++++++++++++++++++
Energy_Score::~Energy_Score()
{
    delete energy_grid;
}

// +++++++++++++++++++++++++++++++++++++++++
void
Energy_Score::input_parameters(Parameter_Reader & parm, bool & primary_score,
                               bool & secondary_score)
{
    use_primary_score = false;
    use_secondary_score = false;

    cout << "\nGrid Score Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    string          tmp;

    if (!primary_score) {
        tmp = parm.query_param("grid_score_primary", "yes", "yes no");
        use_primary_score = tmp == "yes";
        primary_score = use_primary_score;
    }

    if (!secondary_score) {
        //tmp = parm.query_param("grid_score_secondary", "no", "yes no");
        tmp = "no";
	use_secondary_score = (tmp == "yes");
        secondary_score = use_secondary_score;
    }

    if (use_primary_score || use_secondary_score)
        use_score = true;
    else
        use_score = false;

    if (use_score) {
        rep_radius_scale = atof(parm.query_param("grid_score_rep_rad_scale", "1").c_str());
        if (rep_radius_scale <= 0.0) {
                cout << "ERROR:  Parameter must be a float greater than zero."
                        " Program will terminate."
                     << endl;
                exit(0);
        }
        vdw_scale = atof(parm.query_param("grid_score_vdw_scale", "1").c_str());
        if (vdw_scale <= 0.0) {
            bool off = parm.query_param( "grid_score_turn_off_vdw", "yes",
                                         "yes no") == "yes";
            if (!off) {
                cout << "ERROR:  Parameter must be a float greater than zero."
                        " Program will terminate."
                     << endl;
                exit(0);
            }
        }
        es_scale = atof(parm.query_param("grid_score_es_scale", "1").c_str());
        if (es_scale <= 0.0) {
            bool off = parm.query_param( "grid_score_turn_off_es", "yes",
                                         "yes no") == "yes";
            if (!off) {
                cout << "ERROR:  Parameter must be a float greater than zero."
                        " Program will terminate."
                     << endl;
                exit(0);
            }
        }
        grid_lig_efficiency = parm.query_param("grid_lig_efficiency", "no", "yes no") == "yes";
        grid_file_name = parm.query_param( "grid_score_grid_prefix", "grid" );
    }
}

// +++++++++++++++++++++++++++++++++++++++++
void
Energy_Score::initialize(AMBER_TYPER & typer)
{

    if (use_score) {
        cout << "Initializing Grid Score Routines..." << endl;
        energy_grid = new Energy_Grid(); 
        energy_grid->get_instance(grid_file_name);
        init_vdw_energy(typer, energy_grid->att_exp, energy_grid->rep_exp);
    } 
}
// +++++++++++++++++++++++++++++++++++++++++
bool
Energy_Score::compute_score(DOCKMol & mol)
{
    int atom;

    if (use_score) {

        // check to see if molecule is inside grid box
        for (atom = 0; atom < mol.num_atoms; atom++) {
            // is an active atom and is inside the grid
            if (mol.atom_active_flags[atom] && !energy_grid->is_inside_grid_box(mol.x[atom], mol.y[atom], mol.z[atom])) {
                mol.current_score = -MIN_FLOAT;  // arbitrarily large score
                mol.current_data = "ERROR:  Conformation could not be scored."
                    "\nConformation not completely within grid box.\n";
                return false;
            }
        }

        float es_val = 0.0;
        float vdw_val = 0.0;
        float total = 0.0;
        int num_active_heavy = 0;
        for (atom = 0; atom < mol.num_atoms; atom++) {

            if (mol.atom_active_flags[atom]) {

                energy_grid->find_grid_neighbors(mol.x[atom], mol.y[atom], mol.z[atom]);
                //cout << vdwA[mol.amber_at_id[atom]] << " " << vdwB[mol.amber_at_id[atom]] << endl;

                vdw_val +=
                    ((vdwA[mol.amber_at_id[atom]] * energy_grid->interpolate(energy_grid->avdw)) -
                     (vdwB[mol.amber_at_id[atom]] * energy_grid->interpolate(energy_grid->bvdw))) *
                    vdw_scale;
                es_val += mol.charges[atom] * energy_grid->interpolate(energy_grid->es) * es_scale;
                if (mol.atom_types[atom] != "H"){
                    num_active_heavy +=1;
                }
            }
        }

        total = vdw_val + es_val;
        grid_total = total;
        vdw_component = vdw_val;
        es_component = es_val;

        total_eff = total / num_active_heavy;
        vdw_eff = vdw_component / num_active_heavy;
        es_eff = es_component / num_active_heavy;

        if ( !grid_lig_efficiency) {
            mol.current_score = total;
        } else {
            mol.current_score = total_eff;
        }

        //JDB This is assigned above
        //mol.current_score = total;
        mol.current_data = output_score_summary(mol);

    }

    return true;
}

// +++++++++++++++++++++++++++++++++++++++++
string
Energy_Score::output_score_summary(DOCKMol & mol)
{
    ostringstream text;

    if (use_score) {
        
        if (grid_lig_efficiency) {
            text << DELIMITER << setw(STRING_WIDTH) << "Grid_Efficiency_Score:" //Put in underscore to be consistent. Yuchen 10/24/2016
                 << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "Grid_vdw_efficiency:" 
                 << setw(FLOAT_WIDTH) << fixed << vdw_eff << endl
                 << DELIMITER << setw(STRING_WIDTH) << "Grid_es_efficiency:"  
                 << setw(FLOAT_WIDTH) << fixed << es_eff  << endl
                 << DELIMITER << setw(STRING_WIDTH) << "Grid_Score:"  
                 << setw(FLOAT_WIDTH) << fixed << grid_total  << endl
                 << DELIMITER << setw(STRING_WIDTH) << "Grid_vdw_energy:" 
                 << setw(FLOAT_WIDTH) << fixed << vdw_component << endl
                 << DELIMITER << setw(STRING_WIDTH) << "Grid_es_energy:"  
                 << setw(FLOAT_WIDTH) << fixed << es_component  << endl;
        } else { 
            text << DELIMITER << setw(STRING_WIDTH) << "Grid_Score:" //Put in underscore to be consistent. Yuchen 10/24/2016
                << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl;
            text << DELIMITER << setw(STRING_WIDTH) << "Grid_vdw_energy:" 
                << setw(FLOAT_WIDTH) << fixed << vdw_component << endl
                << DELIMITER << setw(STRING_WIDTH) << "Grid_es_energy:"  
                << setw(FLOAT_WIDTH) << fixed << es_component  << endl;
        }
        // Compute lig internal energy with segments (sudipto 12-12-08)
        //float int_vdw_att, int_vdw_rep, int_es;
        //compute_ligand_internal_energy(mol, int_vdw_att, int_vdw_rep, int_es);

        if (use_internal_energy) 
            text << DELIMITER << setw(STRING_WIDTH) << "Internal_energy_repulsive:" //Put in underscore to be consistent. Yuchen 10/24/2016
                 << setw(FLOAT_WIDTH) << fixed << mol.internal_energy << endl;

        // To get the correct internal energy, it MUST be calculated before calling the primary score

    }
    return text.str();
}


// +++++++++++++++++++++++++++++++++++++++++
Continuous_Energy_Score::Continuous_Energy_Score()
{
}

// +++++++++++++++++++++++++++++++++++++++++
Continuous_Energy_Score::~Continuous_Energy_Score()
{
}

// +++++++++++++++++++++++++++++++++++++++++
void
Continuous_Energy_Score::input_parameters(Parameter_Reader & parm,
                                          bool & primary_score,
                                          bool & secondary_score)
{
    string          tmp;

    use_primary_score = false;
    use_secondary_score = false;

    cout << "\nContinuous Energy Score Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    if (!primary_score) {
        tmp = parm.query_param("continuous_score_primary", "no", "yes no");
        if (tmp == "yes")
            use_primary_score = true;
        else
            use_primary_score = false;

        primary_score = use_primary_score;
    }

    if (!secondary_score) {
        //tmp = parm.query_param("continuous_score_secondary", "no", "yes no");
        tmp = "no";
        if (tmp == "yes")
            use_secondary_score = false;
        else
            use_secondary_score = false;

        secondary_score = use_secondary_score;
    }

    if (use_primary_score || use_secondary_score)
        use_score = true;
    else
        use_score = false;

    if (use_score) {
        receptor_filename =
            parm.query_param("cont_score_rec_filename", "receptor.mol2");

        att_exp = atoi(parm.query_param("cont_score_att_exp", "6").c_str());
        if (att_exp <= 0) {
            cout <<
                "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }
        rep_exp = atoi(parm.query_param("cont_score_rep_exp", "12").c_str());
        if (rep_exp <= 0) {
            cout <<
                "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }
        rep_radius_scale = atof(parm.query_param("cont_score_rep_rad_scale", "1").c_str());
        if (rep_radius_scale <= 0.0) {
                cout <<
                    "ERROR:  Parameter must be a float greater than zero. Program will terminate."
                    << endl;
                exit(0);
        }

        //Ask whether we should use a distance dependent dieleectric

        tmp = parm.query_param("cont_score_use_dist_dep_dielectric", "yes", "yes no");
        if (tmp == "yes")
            use_ddd = true;
        else
            use_ddd = false;

        //Use default dielectric = 4.0 for distance dependent dielectric
        //else use d=1.0 assuming gas phase

        char default_diel[4];
        if (use_ddd) strcpy(default_diel, "4.0");
        else strcpy(default_diel,"1.0");

        diel_screen =
            atof(parm.query_param("cont_score_dielectric", default_diel).c_str());
        if (diel_screen <= 0.0) {
            cout <<
                "ERROR: Parameter must be a float greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }

        vdw_scale = atoi(parm.query_param("cont_score_vdw_scale", "1").c_str());
        if (vdw_scale <= 0) {
            bool            off;
            off =
                (parm.query_param("cont_score_turn_off_vdw", "yes", "yes no") ==
                 "yes") ? true : false;
            if (!off) {
                cout <<
                    "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }
        }
        es_scale = atoi(parm.query_param("cont_score_es_scale", "1").c_str());
        if (es_scale <= 0) {
            bool            off;
            off =
                (parm.query_param("cont_score_turn_off_es", "yes", "yes no") ==
                 "yes") ? true : false;
            if (!off) {
                cout <<
                    "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }
        }

    }
}

// +++++++++++++++++++++++++++++++++++++++++
void
Continuous_Energy_Score::initialize(AMBER_TYPER & typer)
{
    ifstream        rec_file;
    bool            read_vdw,
                    use_chem,
                    use_ph4,
                    use_volume;

    if (use_score) {

        init_vdw_energy(typer, att_exp, rep_exp);

        rec_file.open(receptor_filename.c_str());

        if (rec_file.fail()) {
            cout << "Error Opening Receptor File!" << endl;
            exit(0);
        }

        if (!Read_Mol2(receptor, rec_file, false, false, false)) {
            cout << "Error Reading Receptor Molecule!" << endl;
            exit(0);
        }

        rec_file.close();

        read_vdw = true;
        use_chem = false;
        use_ph4  = false;
        use_volume = false;
        typer.prepare_molecule(receptor, read_vdw, use_chem, use_ph4, use_volume);

    }
}

// +++++++++++++++++++++++++++++++++++++++++
bool
Continuous_Energy_Score::compute_score(DOCKMol & mol)
{
    //cout << "Continuous_Energy_Score::compute_score Enter" << endl;
    float           vdw_val,
                    es_val,
                    total;
    int             i,
                    j;
    float           dist;

    total = vdw_val = es_val = 0.0;

    if (use_score == 1) {

        for (i = 0; i < mol.num_atoms; i++) {

            for (j = 0; j < receptor.num_atoms; j++) {

                if (mol.atom_active_flags[i]) {

                    dist =
                        sqrt(((mol.x[i] - receptor.x[j])*(mol.x[i] - receptor.x[j])) +
                             ((mol.y[i] - receptor.y[j])*(mol.y[i] - receptor.y[j])) + 
                             ((mol.z[i] - receptor.z[j])*(mol.z[i] - receptor.z[j])));

                    // for covalent docking we use a dummy atoms that likely overlaps with a receptor atom.
                    // if the distance is 0.0 then check if it is a dummy atom.
                    if (dist < 0.0001 and mol.atom_types[i] == "Du"){
                       // the dummy should have 0.0 charge and 0.0 vdw radius.
                       // skipe to next atom
                       cout << "Dummy atom too close. skip to next atom. " << endl;
                       continue;
                    }

                    vdw_val +=
                        (((vdwA[mol.amber_at_id[i]] *
                           vdwA[receptor.amber_at_id[j]]) / pow(dist,
                                                                rep_exp)) -
                         ((vdwB[mol.amber_at_id[i]] *
                           vdwB[receptor.amber_at_id[j]]) / pow(dist,
                                                                att_exp))) *
                        vdw_scale;

                    if (use_ddd) es_val +=
                        ((332 * mol.charges[i] * receptor.charges[j]) /
                         ((dist*dist) * diel_screen)) * es_scale;
                    else es_val +=
                        ((332 * mol.charges[i] * receptor.charges[j]) /
                         (dist * diel_screen)) * es_scale;
                }
            }
        }

        total = vdw_val + es_val;

        vdw_component = vdw_val;
        es_component = es_val;

        mol.current_score = total;
        mol.current_data = output_score_summary(mol);

        // add_score_to_mol(mol, "Energy Score", total);

    }
    //cout << "Continuous_Energy_Score::compute_score Exit" << endl;

    return true;
}

// +++++++++++++++++++++++++++++++++++++++++
string
Continuous_Energy_Score::output_score_summary(DOCKMol & mol)
{
    ostringstream text;

    if (use_score) {
        text << DELIMITER << setw(STRING_WIDTH) << "Continuous_Score:"
             << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl
             << DELIMITER << setw(STRING_WIDTH) << "Continuous_vdw_energy:"
             << setw(FLOAT_WIDTH) << fixed << vdw_component << endl
             << DELIMITER << setw(STRING_WIDTH) << "Continuous_es_energy:"
             << setw(FLOAT_WIDTH) << fixed << es_component << endl;

        if (use_internal_energy)
            text << DELIMITER << setw(STRING_WIDTH) << "Internal_energy_repulsive:"
                 << setw(FLOAT_WIDTH) << fixed << mol.internal_energy << endl;

    }
    return text.str();
}


// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++
Contact_Score::Contact_Score()
{
    bump_grid    = NULL;
    contact_grid = NULL;
}

// +++++++++++++++++++++++++++++++++++++++++
Contact_Score::~Contact_Score()
{
    delete bump_grid;
    delete contact_grid;
}

// +++++++++++++++++++++++++++++++++++++++++
void
Contact_Score::input_parameters(Parameter_Reader & parm, bool & primary_score,
                                bool & secondary_score)
{
    string          tmp;

    use_primary_score = false;
    use_secondary_score = false;

    cout << "\nContact Score Parameters\n"
            "--------------------------------------------------------"
            "----------------------------------"
         << endl;

    if (!primary_score) {
        tmp = parm.query_param("contact_score_primary", "no", "yes no");
        use_primary_score = (tmp == "yes");
        primary_score = use_primary_score;
    }

    if (!secondary_score) {
        //tmp = parm.query_param("contact_score_secondary", "no", "yes no");
        tmp = "no";
	use_secondary_score = (tmp == "yes");
        secondary_score = use_secondary_score;
    }

    if (use_primary_score || use_secondary_score)
        contact_score = 1;
    else
        contact_score = 0;

    if (contact_score == 1) {
        use_score = true;
        cutoff_distance = atof
            (parm.query_param("contact_score_cutoff_distance", "4.5").c_str());
        if (cutoff_distance <= 0.0) {
            cout << "ERROR:  Parameter must be a float greater than zero."
                    "  Program will terminate."
                 << endl;
            exit(0);
        }
        clash_overlap = atof
            (parm.query_param("contact_score_clash_overlap", "0.75").c_str());
        if (clash_overlap <= 0.0) {
            cout << "ERROR:  Parameter must be a float greater than zero."
                    "  Program will terminate."
                 << endl;
            exit(0);
        }
        clash_penalty = atof
            (parm.query_param("contact_score_clash_penalty", "50").c_str());
        if (clash_penalty <= 0.0) {
            cout << "ERROR:  Parameter must be a float greater than zero."
                    "  Program will terminate."
                 << endl;
            exit(0);
        }
        grid_file_name = parm.query_param( "contact_score_grid_prefix", "grid" );
    }
}

// +++++++++++++++++++++++++++++++++++++++++
void
Contact_Score::initialize(AMBER_TYPER & typer)
{

    if (contact_score == 1) {
        cout << "Initializing Contact Score Routines..." << endl;
        contact_grid = new Contact_Grid(); 
        contact_grid->get_instance(grid_file_name);
        bump_grid = new Bump_Grid(); 
        bump_grid->get_instance(grid_file_name);
        use_score = true;
    } else
        use_score = false;
}

// +++++++++++++++++++++++++++++++++++++++++
// Note that there is a seg fault when bump_filter and Contact_Score are used together. Trent Balius 2012
bool
Contact_Score::compute_score(DOCKMol & mol)
{
    int atom;

    if (contact_score == 1) {

        for (atom = 0; atom < mol.num_atoms; atom++) {
            // is an active atom and is inside the grid
            if ( mol.atom_active_flags[atom] && !contact_grid->is_inside_grid_box(mol.x[atom], mol.y[atom], mol.z[atom])) {
                mol.current_score = -MIN_FLOAT;  // arbitrarily large score
                mol.current_data = "ERROR:  Conformation could not be scored."
                    "\nConformation not completely within grid box.\n";
                return false;
            }
        }

        float cnt_score = 0.0;
        int vdw_threshold = 0;
        for (atom = 0; atom < mol.num_atoms; atom++) {

            if (mol.amber_at_heavy_flag[atom]) {

                if (mol.atom_active_flags[atom]) {

                    contact_grid->find_grid_neighbors(mol.x[atom], mol.y[atom], mol.z[atom]);

                    vdw_threshold = NINT(10 * mol.amber_at_radius[atom]);

                    if (vdw_threshold >= (int) bump_grid->bump[contact_grid->nearest_neighbor]) {
                        cnt_score += clash_penalty;
                    } else {
                        cnt_score += (float) contact_grid->cnt[contact_grid->nearest_neighbor];
                    }
                }
            }
        }

        mol.current_score = cnt_score;
        mol.current_data = output_score_summary(cnt_score);

    }

    return true;
}

// +++++++++++++++++++++++++++++++++++++++++
string
Contact_Score::output_score_summary(float score)
{
    ostringstream text;

    if (use_score) {
        text << DELIMITER << setw(STRING_WIDTH) << "Contact_Score:"
             << setw(FLOAT_WIDTH) << fixed << score << endl
        ;
    }
    return text.str();
}
