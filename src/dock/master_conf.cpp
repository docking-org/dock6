#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits.h>
#include <sstream>
#include <time.h>

#include "master_conf.h"
#include "master_score.h"
#include "simplex.h"
#include "trace.h"

class Bump_Filter;
class Parameter_Reader;

using namespace std;


// +++++++++++++++++++++++++++++++++++++++++
void
Master_Conformer_Search::input_parameters(Parameter_Reader & parm)
{
    Trace trace( "Master_Conformer_Search::input_parameters" );
    flexible_ligand = false;
    genetic_algorithm = false;
    denovo_design = false;



    // cout << "\nFlexible Ligand Parameters" << endl;
    // cout <<
    //     "------------------------------------------------------------------------------------------"
    //     << endl;
    //
    //
    // flexible_ligand =
    //     (parm.query_param("flexible_ligand", "yes", "yes no") == "yes") ? true : false;
    //
    // method = 0;                 // assume rigid ligand. Method 0 = No anchor & grow
    //
    // if (flexible_ligand) {
    //
    //     method = 1;   // perform Anchor & Grow
    //     c_ag_conf.input_parameters(parm);
    //     c_ag_conf.initialize_internal_energy_parms(use_internal_energy, ie_rep_exp, ie_att_exp, ie_diel); 
    //
    //             /**** Currently Disabled ****
    //             // present the HDB option
    //             if(method == 0) {
    //                     c_hdb_conf.input_parameters(parm);
    // 
    //                     if(c_hdb_conf.flexible_ligand)
    //                             method = 2;
    //             }
    //             ****************************/
    // 
    //      trace.note( "Performing anchor and grow docking, or fixed anchor docking." );
    //
    // } else
    //      trace.note( "Performing rigid docking, single point calc, "
    //                  "or minimization without orientation." );


    // WJA 08/22/2011 - Does this thing need an overhaul? - instead of 'Flexible Ligand Parameters',
    // maybe should ask for a type of sampling function:
    //     (0=rigid, 1=flex, 2=denovo, 3=hdb, 4=ga, etc.)

    cout << "\nConformational Sampling Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    // Other search methods would include HDB and possibly a Genetic Algorithm
    //conformer_search_type = parm.query_param("conformer_search_type", "flex", "rigid | flex | denovo | covalent | genetic | HDB");
    conformer_search_type = parm.query_param("conformer_search_type", "flex", "rigid | flex | denovo | covalent | genetic ");

    if (conformer_search_type.compare("rigid") == 0) {
        method = 0;
        trace.note("Performing rigid docking, single point calc, or minimization without orientation.");
    } 
    else if (conformer_search_type.compare("flex") == 0) {
        method = 1;
        flexible_ligand = true;
        c_ag_conf.input_parameters(parm);
        //c_ag_conf.initialize_internal_energy_parms(use_internal_energy, ie_rep_exp, ie_att_exp, ie_diel); 
        trace.note("Performing anchor and grow docking, or fixed anchor docking.");

    }
    else if (conformer_search_type.compare("denovo") == 0) {
        method = 2;
        flexible_ligand = true;
        denovo_design = true;
        c_dn_build.input_parameters(parm);
        //c_dn_build.initialize_internal_energy_parms(use_internal_energy, ie_rep_exp, ie_att_exp, ie_diel); 
        trace.note("Performing de novo molecule growth.");

    }
    else if (conformer_search_type.compare("covalent") == 0) {
        method = 3;
        flexible_ligand = true;
        c_cg_conf.input_parameters(parm);
        //c_dn_build.initialize_internal_energy_parms(use_internal_energy, ie_rep_exp, ie_att_exp, ie_diel);
        trace.note("Performing covalent molecule growth.");
        //
    }
    else if (conformer_search_type.compare("genetic") == 0) {
         method = 4;
         flexible_ligand = true;
         genetic_algorithm = true;
         c_ga_recomb.input_parameters(parm);
         //c_dn_build.input_parameters(parm);
         //c_ga_recomb.initialize_internal_energy_parms(use_internal_energy, ie_rep_exp, ie_att_exp, ie_diel);
         trace.note("Performing genetic algorithm.");

    }
    else if (conformer_search_type.compare("HDB") == 0) {
         cout << "HDB is turned off" << endl;
         exit(0); 
         method = 5;
         flexible_ligand = true;
         c_hdb_conf.input_parameters(parm);
         //c_dn_build.input_parameters(parm);
         //c_ga_recomb.initialize_internal_energy_parms(use_internal_energy, ie_rep_exp, ie_att_exp, ie_diel);
         trace.note("Performing genetic algorithm.");
       
    }
         

    else {
        cout <<"The type of conformational sampling you have chosen does not make sense\n";
        method = -1;
    }


    cout << "\nInternal Energy Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    use_internal_energy = (parm.query_param("use_internal_energy", "yes", "yes no") == "yes");

    if (use_internal_energy) {
        ie_rep_exp = atoi(parm.query_param("internal_energy_rep_exp", "12").c_str());
        if (ie_rep_exp <= 0) {
            cout << "ERROR: internal_energy_rep_exp needs to be positive" << endl;
            exit(0);
        }
        ie_cutoff = atof(parm.query_param("internal_energy_cutoff", "100.0").c_str());
        if (ie_cutoff <= 0.0) {
            cout << "ERROR: internal_energy_cutoff needs to be positive" << endl;
            exit(0);
        }

        ie_att_exp = 6;
        ie_diel = 4.0;
        cout <<"Note: Internal energy only includes repulsive VDW for growth and/or minimization."
             << endl;

        if (method == 1)
            c_ag_conf.initialize_internal_energy_parms(use_internal_energy, ie_rep_exp, ie_att_exp, ie_diel, ie_cutoff); 
        else if (method == 2)
            c_dn_build.initialize_internal_energy_parms(use_internal_energy, ie_rep_exp, ie_att_exp, ie_diel, ie_cutoff); 
        else if (method == 3)
            c_cg_conf.initialize_internal_energy_parms(use_internal_energy, ie_rep_exp, ie_att_exp, ie_diel, ie_cutoff); 
        else if (method == 4)
            c_ga_recomb.initialize_internal_energy_parms(use_internal_energy, ie_rep_exp, ie_att_exp, ie_diel, ie_cutoff);

    }
    else {
            c_ag_conf.initialize_internal_energy_null(use_internal_energy);
   }

}


// +++++++++++++++++++++++++++++++++++++++++
void
Master_Conformer_Search::initialize()
{
    Trace trace( "Master_Conformer_Search::initialize" );

    switch (method) {

    case 0: //Rigid
        break;

    case 1: //Flex
        trace.note( "doing Flex docking" );
        c_ag_conf.initialize();
        c_ag_conf.count_conf_num = 0;  // trent balius Dec. 15, 2008
        break;

    case 2: //De novo
        trace.note( "doing de novo growth" );
        c_dn_build.initialize();
        break;
    case 3: //Covalent Flex
        trace.note( "doing covalent docking" );
        c_cg_conf.initialize();
        c_cg_conf.count_conf_num = 0;  
        break;

    // WJA 08/22/2011 - de novo is the new case 2, HDB is now case 3
    // CDS 10/2014 - GA is the new case 3, HDB is now case 4
    case 4: //Genetic algorithm
         trace.note( "doing genetic algorithm" );
         c_ga_recomb.initialize();
         //c_dn_build.initialize();
         break;


        // case 5:
        // c_hdb_conf.initialize();
        // break;

        // Note to whoever implemented HDB Search:
        // if you have trouble making your HDB search to work after we
        // included internal energy during ag_search, please contact
        // Trent Balius & Sudipto Mukherjee (Stony Brook) for help

    }

}


// +++++++++++++++++++++++++++++++++++++++++
void
Master_Conformer_Search::prepare_molecule(DOCKMol & mol)
{
    //cout << "Master_Conformer_Search::prepare_molecule(DOCKMol & mol)" << endl;
    switch (method) {

    case 0:  //Rigid
        more_anchors = true;
        last_conformer = true;
        copy_molecule(orig, mol); //trent 2009-02-12
        break;

    case 1:  //Flex
        c_ag_conf.prepare_molecule(mol);
        break;
    case 3:  //covalent
        c_cg_conf.prepare_molecule(mol);
        break;

    // WJA 08/22/2011 - not sure if we need to do this in dock.cpp yet
    // case 2:  //De Novo
    //    c_dn_build.prepare_molecule(mol);
    //    break;

        // case 4:
        // c_hdb_conf.prepare_molecule(mol);
        // break;

    }

}


// +++++++++++++++++++++++++++++++++++++++++
bool
Master_Conformer_Search::next_anchor(DOCKMol & mol)
{
    bool            return_val;

    switch (method) {

    case 0://Rigid
        if (more_anchors) {
            more_anchors = false;
            return_val = true;
        } else
            return_val = false;
        break;

    case 1://Flex
        return_val = c_ag_conf.next_anchor(mol);
        break;

    case 3://covalent
        return_val = c_cg_conf.next_anchor(mol);
        more_anchors = false;
        return_val = false;
        break;

        // case 4:
        // return_val = c_hdb_conf.next_anchor(mol);
        // break;

    default:
        throw;
        break;
    }

    return return_val;
}


// +++++++++++++++++++++++++++++++++++++++++
bool
Master_Conformer_Search::submit_anchor_orientation(DOCKMol & mol,
                                                   bool more_orients)
{
    bool            return_val;

    switch (method) {

    case 0:  // if rigid
        if (more_orients) {
            more_anchors = true;
            last_conformer = true;
            return_val = true;
        } else if (last_conformer) {
            more_anchors = true;
            last_conformer = false;
            return_val = true;
        } else
            return_val = false;
        break;

    case 1: // if flex
        return_val =
            c_ag_conf.submit_anchor_orientation(mol, more_orients);
        break;
    case 3: // if covalent
        return_val =
            c_cg_conf.submit_anchor_orientation(mol, more_orients);
        break;

        // case 2:
        // return_val = c_hdb_conf.submit_anchor_orientation(mol, more_orients, 
        // score);
        // break;

    default:
        throw;
        break;
    }

    return return_val;
}


// +++++++++++++++++++++++++++++++++++++++++
void
Master_Conformer_Search::grow_periphery(Master_Score & score,
                                        Simplex_Minimizer & simplex, Bump_Filter & bump)
{

    // sudipto & trent: this command prevents segfault when no scoring function is specified
//    if (!score.use_score){
//       cout << "Currently a scoring function must be specified" << endl;
//       return;
//    }

    // always pass this to base_score (even if score.use_score is false)
    score.primary_score->method = method;
 
    switch (method) {

    case 0: //Rigid 
        //cout << "Master_Conformer_Search::grow_periphery -> Rigid docking" << endl;
        // initialization of internal energy 
        // if the flag use_internal_energy is true.

        // we should move the internal energy initialization 
        // because it really should be called only once for each
        // new ligand during rigid docking
        if (use_internal_energy) { // initializes internal energy ones for every achor
           score.primary_score->ie_att_exp = ie_att_exp;
           score.primary_score->ie_rep_exp = ie_rep_exp;
           score.primary_score->ie_diel = ie_diel;
           score.primary_score->use_internal_energy = true;

           // need to use DOCKMol with radii and segments assigned
           // it does not matter which atoms are labeled active
           score.primary_score->initialize_internal_energy(orig);
           initialize_once = false;
        } 
        
        break;

    case 1: // if flex call anchor and grow
        // initialize of internal energy is performed in grow_periphery
        c_ag_conf.grow_periphery(score, simplex, bump);
        break;

    case 3: //covalent
        c_cg_conf.grow_periphery(score, simplex, bump);
        break;

        // case 2:
        // c_hdb_conf.grow_periphery(score, bump);
        // break;

    }

}


// +++++++++++++++++++++++++++++++++++++++++
bool
Master_Conformer_Search::next_conformer(DOCKMol & mol)
{
    bool            return_val;

    switch (method) {

    case 0:  // if Rigid
        if (more_anchors) {
            more_anchors = false;
            last_conformer = true;
            return_val = true;
        } else
            return_val = false;
        break;

    case 1: // if flex
        return_val = c_ag_conf.next_conformer(mol);
        break;

    case 3: // if covalent
        return_val = c_cg_conf.next_conformer(mol);
        break;

        // case 4:
        // return_val = c_hdb_conf.next_conformer(mol);
        // break;

    default:
        throw;
        break;
    }

    return return_val;
}


