#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <sstream>
#include "amber_typer.h"
#include "base_score.h"
#include "dockmol.h"
#include "grid.h" 
#include "score_gist.h"
#include "utils.h"

using namespace std;


// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++
GIST_Score::GIST_Score()
{
    //sigma2 = 0.0;
    gist_grid = NULL;
    gist_H_grid = NULL;
}

// +++++++++++++++++++++++++++++++++++++++++
GIST_Score::~GIST_Score()
{
    delete gist_grid;
    delete gist_H_grid;
}

// +++++++++++++++++++++++++++++++++++++++++
void
GIST_Score::input_parameters_main(Parameter_Reader & parm, string parm_head)
{
        att_exp = atoi(parm.query_param(parm_head+"att_exp", "6").c_str());
        if (att_exp <= 0) {
            cout <<
                "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }
        rep_exp = atoi(parm.query_param(parm_head+"rep_exp", "12").c_str());
        if (rep_exp <= 0) {
            cout <<
                "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }
        gist_scale = atof(parm.query_param(parm_head+"gist_scale", "-1").c_str());
        if (gist_scale == 0.0) {
             cout << gist_scale << endl;
             cout << "ERROR:  Parameter must not be zero."
                     " Program will terminate."
                  << endl;
             exit(0);
        }
        //if (gist_scale <= 0.0) {
        //     cout << "ERROR:  Parameter must be a float greater than zero."
        //             " Program will terminate."
        //          << endl;
        //     exit(0);
        //}
        gist_type = parm.query_param(parm_head+"gist_type", "trilinear", "trilinear displace blurry_displace").c_str();

        //if (gist_type == "blurry_displace") {
        //    sigma2 = pow(atof(parm.query_param(parm_head+"sigma", "2.0").c_str()),2.0); //square it. 
        //}

        grid_file_name = parm.query_param( parm_head+"grid_file", "grid.dx" );
        if (gist_type == "trilinear") {
            grid_H_file_name = parm.query_param( parm_head+"Hydrogen_grid_file", "grid_h.dx" );
        }
}

// +++++++++++++++++++++++++++++++++++++++++
void
GIST_Score::input_parameters(Parameter_Reader & parm, bool & primary_score,
                               bool & secondary_score)
{
    use_primary_score = false;
    use_secondary_score = false;

    cout << "\nGIST Score Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;
    use_primary_score = false;
    use_secondary_score = false;
    string          tmp;

    if (!primary_score) {
        tmp = parm.query_param("gist_score_primary", "no", "yes no");
        use_primary_score = tmp == "yes";
        primary_score = use_primary_score;
    }
    
    //LEP - removed secondary score in gist 
    //if (!secondary_score) {
    //    tmp = parm.query_param("gist_score_secondary", "no", "yes no");
    //    use_secondary_score = tmp == "yes";
    //    secondary_score = use_secondary_score;
    //}

    secondary_score = use_secondary_score; // false

    if (use_primary_score || use_secondary_score)
        use_score = true;
    else
        use_score = false;

    if (use_score) {
        input_parameters_main(parm, "gist_score_");
   }
}

// +++++++++++++++++++++++++++++++++++++++++
void
GIST_Score::initialize(AMBER_TYPER & typer)
{

    if (use_score) {
        cout << "Initializing GIST Grid Score Routines..." << endl;
        gist_grid = new GIST_Grid(); 
        gist_grid->get_instance(grid_file_name);
        gist_grid->att_exp = att_exp;
        gist_grid->rep_exp = rep_exp;
        init_vdw_energy(typer, gist_grid->att_exp, gist_grid->rep_exp);
        if (gist_type == "trilinear") {
           cout << "Initializing GIST H Grid..." << endl;
           gist_H_grid->got_the_grid = false;
           gist_H_grid = new GIST_Grid(); 
           gist_H_grid->get_instance(grid_H_file_name);
           init_vdw_energy(typer, gist_H_grid->att_exp, gist_H_grid->rep_exp);
        }
    } 
}
// +++++++++++++++++++++++++++++++++++++++++
bool
GIST_Score::compute_score(DOCKMol & mol)
{
    int atom;

    if (use_score) {

        // check to see if molecule is inside grid box
        for (atom = 0; atom < mol.num_atoms; atom++) {
            // is an active atom and is inside the grid
            if (mol.atom_active_flags[atom] && !gist_grid->is_inside_grid_box(mol.x[atom], mol.y[atom], mol.z[atom])) {
                mol.current_score = -MIN_FLOAT;  // arbitrarily large score
                mol.current_data = "ERROR:  Conformation could not be scored."
                    "\nConformation not completely within grid box.\n";
                return false;
            }
            if (gist_type == "trilinear") {
                if (mol.atom_active_flags[atom] && !gist_H_grid->is_inside_grid_box(mol.x[atom], mol.y[atom], mol.z[atom])) {
                    mol.current_score = -MIN_FLOAT;  // arbitrarily large score
                    mol.current_data = "ERROR:  Conformation could not be scored."
                        "\nConformation not completely within grid box.\n";
                    return false;
                }
            }
        }
        float gist_total = 0.0;
        if (gist_type == "trilinear") {
            for (atom = 0; atom < mol.num_atoms; atom++) {
                  if (mol.atom_active_flags[atom]) {
                    if (mol.atom_types[atom] == "H"){
                       //cout << mol.atom_types[atom] << atom << endl;
                       gist_H_grid->find_grid_neighbors(mol.x[atom], mol.y[atom], mol.z[atom]);
                       //gist_total = gist_total + (gist_H_grid->interpolate(gist_H_grid->gist));
                       gist_total = gist_total + (gist_H_grid->interpolate(gist_H_grid->gist));
                       //cout << gist_total << endl;
                    }
                    else { //heavy atoms
                       //cout << mol.atom_types[atom] << atom << endl;
                       gist_grid->find_grid_neighbors(mol.x[atom], mol.y[atom], mol.z[atom]);
                       gist_total = gist_total + (gist_grid->interpolate(gist_grid->gist));
                       //cout << gist_total << endl;
                    }
                  } 
            }
            gist_total = gist_grid->vol*gist_total; //
        } 
        else if (gist_type == "displace") { // perform displacement
            bool *displace;
            displace = new bool[gist_grid->size];
            for (atom = 0; atom < mol.num_atoms; atom++) {
                  if (mol.atom_active_flags[atom]) {
                      gist_grid->find_grid_neighbors(mol.x[atom], mol.y[atom], mol.z[atom]);
//                     radius = (sqrt(2.0)* sra(atom_vdwtype(atomindex)) /
//                          &               srb(atom_vdwtype(atomindex)))**(1.0/3.0) / 2
                      float radius = pow((sqrt(2.0)*vdwA[mol.amber_at_id[atom]]/vdwB[mol.amber_at_id[atom]]),(1.0/3.0)) / 2.0;
                      //cout << mol.amber_at_id[atom] << " " << mol.amber_at_id[atom] << endl;
                      //cout << vdwA[mol.amber_at_id[atom]] << " " << vdwB[mol.amber_at_id[atom]] << endl;
                      cout << mol.atom_types[atom] << " radius = "<< radius << endl;
                      float gist_atom = gist_grid->atomic_displacement(mol.x[atom], mol.y[atom], mol.z[atom], radius, displace);
                      cout << gist_atom << endl;
                      gist_total = gist_total + gist_atom; 
                  }
            }
            //gist_grid->write_gist_grid("dis_gist.dx", displace); // this is for debuging
            gist_total = gist_grid->vol*gist_total;
            //exit(0);//
        }
        else if (gist_type == "blurry_displace") {
            for (atom = 0; atom < mol.num_atoms; atom++) {
                  if (mol.atom_active_flags[atom]) {
                      gist_grid->find_grid_neighbors(mol.x[atom], mol.y[atom], mol.z[atom]);
                      // consider moving this to so we percompute radius. 
                      float radius = pow((sqrt(2.0)*vdwA[mol.amber_at_id[atom]]/vdwB[mol.amber_at_id[atom]]),(1.0/3.0)) / 2.0;
                      //cout << mol.amber_at_id[atom] << " " << mol.amber_at_id[atom] << endl;
                      //cout << vdwA[mol.amber_at_id[atom]] << " " << vdwB[mol.amber_at_id[atom]] << endl;
                      cout << mol.atom_types[atom] << " radius = "<< radius << endl;
                      float sigma2 = pow(radius/2.0,2.0);
                      float gist_atom = gist_grid->atomic_blurry_displacement(mol.x[atom], mol.y[atom], mol.z[atom], radius, sigma2);
                      cout << gist_atom << endl;
                      gist_total = gist_total + gist_atom; 
                  }
            }
            //gist_grid->write_gist_grid("dis_gist.dx", displace); // this is for debuging
            gist_total = gist_grid->vol*gist_total;
            //exit(0);//

        }
/*
        float es_val = 0.0;
        float vdw_val = 0.0;
        for (atom = 0; atom < mol.num_atoms; atom++) {

            if (mol.atom_active_flags[atom]) {

                energy_grid->find_grid_neighbors(mol.x[atom], mol.y[atom], mol.z[atom]);

                vdw_val +=
                    ((vdwA[mol.amber_at_id[atom]] * energy_grid->interpolate(energy_grid->avdw)) -
                     (vdwB[mol.amber_at_id[atom]] * energy_grid->interpolate(energy_grid->bvdw))) *
                    vdw_scale;
                es_val += mol.charges[atom] * energy_grid->interpolate(energy_grid->es) * es_scale;
            }
        }

        total = vdw_val + es_val;

        vdw_component = vdw_val;
        es_component = es_val;
*/
        //cout << "gist_total = "<< gist_total << endl;
        mol.current_score = gist_scale*gist_total;
        mol.current_data = output_score_summary(mol);

    }

    return true;
}

// +++++++++++++++++++++++++++++++++++++++++
string
GIST_Score::output_score_summary(DOCKMol & mol)
{
    ostringstream text;

    if (use_score) {

        text << DELIMITER << setw(STRING_WIDTH) << "GIST_Score:" //Put in underscore to be consistent. Yuchen 10/24/2016
             << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl;
/*        text << DELIMITER << setw(STRING_WIDTH) << "Grid_vdw_energy:" 
             << setw(FLOAT_WIDTH) << fixed << vdw_component << endl
             << DELIMITER << setw(STRING_WIDTH) << "Grid_es_energy:"  
             << setw(FLOAT_WIDTH) << fixed << es_component  << endl;
*/
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


