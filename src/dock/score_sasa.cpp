#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <time.h>
#include "amber_typer.h"
#include "sasa.h"
#include "score_sasa.h"
#include "score_solvent.h"

using namespace std;

/************************************************/
// SASA Descriptor Score 
// Written by Yulin Huang, Trent Balius and Sudipto Mukherjee
// Rizzo Group, Stony Brook University
/************************************************/

// SASA Descriptor Score has been renamed to SASA Score
// in order to prevent confusion with Descriptor Score
// renaming it will also make incorportating it into
// score_descriptor.cpp much cleaner



/************************************************/
SASA_score::SASA_score()
{
}

/************************************************/
SASA_score::~SASA_score()
{
}

/************************************************/
void
SASA_score::input_parameters(Parameter_Reader & parm, bool & primary_score,
                             bool & secondary_score)
{
    string          tmp;

    use_primary_score = false;
    use_secondary_score = false;

    cout << "\nSASA Score Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    if (!primary_score) {
        tmp = parm.query_param("SASA_score_primary", "no", "yes no");
        if (tmp == "yes")
            use_primary_score = true;
        else
            use_primary_score = false;

        primary_score = use_primary_score;
    }

    secondary_score = false;
    //if (!secondary_score) {
    //    //tmp = parm.query_param("SASA_score_secondary", "no", "yes no");
    //    tmp == "no";
    //    if (tmp == "yes")
    //        use_secondary_score = true;
    //    else
    //        use_secondary_score = false;

    //    secondary_score = use_secondary_score;
    //}

    if (use_primary_score || use_secondary_score)
        use_score = true;
    else
        use_score = false;

    if (use_score) {
        rec_filename =
            parm.query_param("SASA_score_rec_filename",
                             "receptor.mol2");
    }
    verbose = 0 != parm.verbosity_level();
}

/************************************************/
void
SASA_score::initialize(AMBER_TYPER & typer)
{
    ifstream        rec_file;

    use_score = true;

    if (use_score) {
        cout << "Initializing SASA Score Routines..." << endl;


        // read in receptor file
        rec_file.open(rec_filename.c_str());
        if (rec_file.fail()) {
            cout << "Error Opening Receptor File!" << endl;
            exit(0);
        }

        if (!Read_Mol2(receptor, rec_file, false, false, false)) {
            cout << "Error Reading Receptor Molecule!" << endl;
            exit(0);
        }
        rec_file.close();

        bool read_vdw = true;
        bool use_chem = false;
        bool use_ph4  = false;
        bool use_volume = false;
        typer.prepare_molecule(receptor, read_vdw, use_chem, use_ph4, use_volume);

        // calculate receptor alone

        prepare_receptor();
    }

}

/************************************************/
bool
SASA_score::compute_score(DOCKMol & mol)
{

    if (use_score) {

        prepare_ligand(mol);
        prepare_complex(mol);

        // percent of the ligand that is solvent exposed
        percent_lig_exposed          = (com_sasa_lig_tot / lig_sasa_tot)* 100;
        // the percent of the hydrophobic portion of ligand which is buried
        percent_lig_buried_is_phobic = ((lig_sasa_phobic - com_sasa_lig_phobic) / lig_sasa_phobic) * 100;
        // the percent of the hydrophobic portion of ligand which is exposed
        // double-check that we do not divide by zero and output -nan
        if (lig_sasa_phobic==0){
          percent_lig_phobic_exposed=lig_sasa_phobic;
        }
        else {
          percent_lig_phobic_exposed = (com_sasa_lig_phobic / lig_sasa_phobic) * 100;
        }

        // the percent of the hydrophilic portion of ligand which is exposed
        // double-check that we do not divide by zero and output -nan
        if (lig_sasa_philic==0){
          percent_lig_philic_exposed=lig_sasa_philic;
        }
        else {
        percent_lig_philic_exposed  = (com_sasa_lig_philic / lig_sasa_philic) * 100;
        }

        // the percent of the "other" portion of ligand which is exposed
        // double-check that we do not divide by zero and output -nan
        if (lig_sasa_other ==0){ 
          percent_lig_other_exposed=lig_sasa_other;
        }
        else {
          percent_lig_other_exposed  = (com_sasa_lig_other / lig_sasa_other) * 100;
        }

        // the percent of the hydrophobic portion of the receptor which is buried
        percent_rec_buried_is_phobic = ((rec_sasa_phobic - com_sasa_rec_phobic) / rec_sasa_phobic) * 100;

        // the percent of buried ligand which is hydrophobic
        percent_buried_lig_phobic = ((lig_sasa_phobic - com_sasa_lig_phobic) / (lig_sasa_tot - com_sasa_lig_tot)) * 100;
        // the percent of buried receptor which is hydrophobic
        percent_buried_rec_phobic = ((rec_sasa_phobic - com_sasa_rec_phobic) / (rec_sasa_tot - com_sasa_rec_tot)) * 100;

        //total_score = -1.0 * (percent_lig_exposed + percent_lig_buried_is_phobic + percent_rec_buried_is_phobic);
        total_score = percent_lig_exposed; 
/*
                    rec_sasa_tot        * parm1 +
                    rec_sasa_phobic     * parm1 +
                    rec_sasa_philic     * parm +
                    rec_sasa_other      * parm +
                    rec_sasa_tot_lig    * parm +
                    rec_sasa_phobic_lig * parm +
                    rec_sasa_phobic_lig * parm +
                    rec_sasa_other_lig  * parm +
                    rec_sasa_tot_rec    * parm +
                    rec_sasa_philic_rec * parm +
                    rec_sasa_philic_rec * parm +
                    rec_sasa_other_rec  * parm +

                    lig_sasa_tot        * parm +
                    lig_sasa_phobic     * parm +
                    lig_sasa_philic     * parm +
                    lig_sasa_other      * parm +
                    lig_sasa_tot_lig    * parm +
                    lig_sasa_phobic_lig * parm +
                    lig_sasa_phobic_lig * parm +
                    lig_sasa_other_lig  * parm +
                    lig_sasa_tot_rec    * parm +
                    lig_sasa_philic_rec * parm +
                    lig_sasa_philic_rec * parm +
                    lig_sasa_other_rec  * parm +

                    com_sasa_tot        * parm +
                    com_sasa_phobic     * parm +
                    com_sasa_philic     * parm +
                    com_sasa_other      * parm +
                    com_sasa_tot_lig    * parm +
                    com_sasa_phobic_lig * parm +
                    com_sasa_phobic_lig * parm +
                    com_sasa_other_lig  * parm +
                    com_sasa_tot_rec    * parm +
                    com_sasa_philic_rec * parm +
                    com_sasa_philic_rec * parm + 
                    com_sasa_other_rec  * parm ;

*/

        mol.current_score = total_score;
        mol.current_data = output_score_summary(mol);

    }

    return true;
}

/************************************************/
string
SASA_score::output_score_summary(DOCKMol & mol)
{
    // SASA decsriptor score

    ostringstream text;

    if (use_score) {
        text << setprecision(3)
             << DELIMITER << setw(STRING_WIDTH) << "SASA_Score:"
             << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl
             << DELIMITER << setw(STRING_WIDTH) << "SASA_com_lig_sasa_tot:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_lig_tot     << endl
             << DELIMITER << setw(STRING_WIDTH) << "SASA_lig_sasa_tot:"
             << setw(FLOAT_WIDTH) << fixed << lig_sasa_tot         << endl
             << DELIMITER << setw(STRING_WIDTH) << "SASA_%_phobic_lig_exposed:"
             << setw(FLOAT_WIDTH) << fixed << percent_lig_phobic_exposed      << endl
             << DELIMITER << setw(STRING_WIDTH) << "SASA_com_lig_sasa_phobic:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_lig_phobic  << endl
             << DELIMITER << setw(STRING_WIDTH) << "SASA_lig_sasa_phobic:"
             << setw(FLOAT_WIDTH) << fixed << lig_sasa_phobic      << endl
             << DELIMITER << setw(STRING_WIDTH) << "SASA_%_philic_lig_exposed:"
             << setw(FLOAT_WIDTH) << fixed << percent_lig_philic_exposed      << endl
             << DELIMITER << setw(STRING_WIDTH) << "SASA_com_lig_sasa_philic:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_lig_philic  << endl
             << DELIMITER << setw(STRING_WIDTH) << "SASA_lig_sasa_philic:"
             << setw(FLOAT_WIDTH) << fixed << lig_sasa_philic      << endl
             << DELIMITER << setw(STRING_WIDTH) << "SASA_%_other_lig_exposed:"
             << setw(FLOAT_WIDTH) << fixed << percent_lig_other_exposed      << endl
             << DELIMITER << setw(STRING_WIDTH) << "SASA_com_lig_sasa_other:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_lig_other  << endl
             << DELIMITER << setw(STRING_WIDTH) << "SASA_lig_sasa_other:"
             << setw(FLOAT_WIDTH) << fixed << lig_sasa_other      << endl
              ;



/*          Below is the original header info printed to the mol2 file when SASA
            score was used. This section has been retained in case users wanted to be able to still 
            print this information by simply uncommenting the code. 

             << DELIMITER << setw(STRING_WIDTH) << "SASA_Score:"
             << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl
             << "###############################################" << endl
             << DELIMITER << setw(STRING_WIDTH) << "%_lig_exposed:"
             << setw(FLOAT_WIDTH) << fixed << percent_lig_exposed  << endl
             << DELIMITER << setw(STRING_WIDTH) << "%_phobic_lig_buried:"
             << setw(FLOAT_WIDTH) << fixed << percent_lig_buried_is_phobic      << endl
             << DELIMITER << setw(STRING_WIDTH) << "%_phobic_rec_buried:"
             << setw(FLOAT_WIDTH) << fixed << percent_rec_buried_is_phobic      << endl
             << DELIMITER << setw(STRING_WIDTH) << "%_lig_phobic_buried/tot_buried:"
             << setw(FLOAT_WIDTH) << fixed << percent_buried_lig_phobic      << endl
             << DELIMITER << setw(STRING_WIDTH) << "%_rec_phobic_buried/tot_buried:"
             << setw(FLOAT_WIDTH) << fixed << percent_buried_rec_phobic      << endl
             *<< "###############################################" << endl
             << DELIMITER << setw(STRING_WIDTH) << "rec_sasa_tot:"
             << setw(FLOAT_WIDTH) << fixed << rec_sasa_tot         << endl
             << DELIMITER << setw(STRING_WIDTH) << "rec_sasa_phobic:"
             << setw(FLOAT_WIDTH) << fixed << rec_sasa_phobic      << endl
             << DELIMITER << setw(STRING_WIDTH) << "rec_sasa_philic:"
             << setw(FLOAT_WIDTH) << fixed << rec_sasa_philic      << endl
             << DELIMITER << setw(STRING_WIDTH) << "rec_sasa_other:"
             << setw(FLOAT_WIDTH) << fixed << rec_sasa_other       << endl
             << DELIMITER << setw(STRING_WIDTH) << "lig_sasa_tot:"
             << setw(FLOAT_WIDTH) << fixed << lig_sasa_tot         << endl
             << DELIMITER << setw(STRING_WIDTH) << "lig_sasa_phobic:"
             << setw(FLOAT_WIDTH) << fixed << lig_sasa_phobic      << endl
             << DELIMITER << setw(STRING_WIDTH) << "lig_sasa_philic:"
             << setw(FLOAT_WIDTH) << fixed << lig_sasa_philic      << endl
             << DELIMITER << setw(STRING_WIDTH) << "lig_sasa_other:"
             << setw(FLOAT_WIDTH) << fixed << lig_sasa_other       << endl
             << DELIMITER << setw(STRING_WIDTH) << "com_sasa_tot:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_tot         << endl
             << DELIMITER << setw(STRING_WIDTH) << "com_sasa_phobic:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_phobic      << endl
             << DELIMITER << setw(STRING_WIDTH) << "com_sasa_philic:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_philic      << endl
             << DELIMITER << setw(STRING_WIDTH) << "com_sasa_other:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_other       << endl
             << DELIMITER << setw(STRING_WIDTH) << "com_sasa_tot_lig:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_lig_tot     << endl
             << DELIMITER << setw(STRING_WIDTH) << "com_sasa_phobic_lig:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_lig_phobic  << endl
             << DELIMITER << setw(STRING_WIDTH) << "com_sasa_philic_lig:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_lig_philic  << endl
             << DELIMITER << setw(STRING_WIDTH) << "com_sasa_other_lig:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_lig_other   << endl
             << DELIMITER << setw(STRING_WIDTH) << "com_sasa_tot_rec:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_rec_tot     << endl
             << DELIMITER << setw(STRING_WIDTH) << "com_sasa_phobic_rec:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_rec_phobic  << endl
             << DELIMITER << setw(STRING_WIDTH) << "com_sasa_philiic_rec:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_rec_philic  << endl
             << DELIMITER << setw(STRING_WIDTH) << "com_sasa_other_rec:"
             << setw(FLOAT_WIDTH) << fixed << com_sasa_rec_other   << endl
              ; */

        if (use_internal_energy)
            text << DELIMITER << setw(STRING_WIDTH) << "Internal_energy_repulsive:"
                 << setw(FLOAT_WIDTH) << fixed << mol.internal_energy << endl;

    }
    return text.str();
}

/************************************************/
void
SASA_score::prepare_receptor()
{

    int             i,
                    j;

    if (use_score) {

        short_rec.allocate_short_arrays(receptor.num_atoms);

        // collect information about receptor
        short_rec.num_atoms = receptor.num_atoms;
        for (i = 0; i < receptor.num_atoms; i++) {
            if (receptor.atom_active_flags[i]) {

                short_rec.atom_names[i] = receptor.atom_names[i];
                short_rec.atom_types[i] = receptor.atom_types[i];

                short_rec.coord[i].v[0] = receptor.x[i];
                short_rec.coord[i].v[1] = receptor.y[i];
                short_rec.coord[i].v[2] = receptor.z[i];

                short_rec.charges[i] = receptor.charges[i];
                short_rec.vdw_radius[i] = receptor.amber_at_radius[i];
                short_rec.gb_radius[i] = receptor.gb_hawkins_radius[i];
                short_rec.gb_scale[i] = receptor.gb_hawkins_scale[i];
            }
        }

        // calculate gpol for receptor
        if (verbose) {
            cout << "\n----------------------------------------------------\n";
            cout << "RECEPTOR\n";
        }

        // calculate solvent accessible surface area of receptor

        // Measure SASA time
        clock_t         start;
        clock_t         stop;

        start = clock();
        sasa           *s = new sasa();
        //float rec_sasa = s->getSASA(short_rec,receptor.num_atoms);

        //com_sasa = getSASA(short_com,receptor.num_atoms);

        // these are not passed to the class SASA_score
        float sasa_lig_tot, sasa_lig_philic, sasa_lig_phobic, sasa_lig_other,
              sasa_rec_tot, sasa_rec_philic, sasa_rec_phobic, sasa_rec_other;

        s->getSASA(short_rec, receptor.num_atoms,
                                   rec_sasa_tot, rec_sasa_philic, rec_sasa_phobic, rec_sasa_other,
                                   sasa_lig_tot, sasa_lig_philic, sasa_lig_phobic, sasa_lig_other,
                                   sasa_rec_tot, sasa_rec_philic, sasa_rec_phobic, sasa_rec_other );

        stop = clock();

//        gnpol_rec = (rec_sasa * 0.00542) + 0.92;
        delete          s;
        if (verbose) {
            // cout.precision(3);
            cout << "SASA=" << rec_sasa_tot << endl
                << "\tTime=" << ((long) stop -
                                 (long) start) /
                (float) CLOCKS_PER_SEC << "s" << endl;
        }

    }
}

/************************************************/
bool
SASA_score::prepare_ligand(DOCKMol & ligand)
{
    int             i,
                    j;
    vector < float >tmp;

    if (use_score) {

        short_lig.allocate_short_arrays(ligand.num_atoms);

        // collect information about ligand
        short_lig.num_atoms = ligand.num_atoms;
        for (i = 0; i < ligand.num_atoms; i++) {
            if (ligand.atom_active_flags[i]) {

                short_lig.atom_names[i] = ligand.atom_names[i];
                short_lig.atom_types[i] = ligand.atom_types[i];

                short_lig.coord[i].v[0] = ligand.x[i];
                short_lig.coord[i].v[1] = ligand.y[i];
                short_lig.coord[i].v[2] = ligand.z[i];

                short_lig.charges[i] = ligand.charges[i];
                short_lig.vdw_radius[i] = ligand.amber_at_radius[i];
                short_lig.gb_radius[i] = ligand.gb_hawkins_radius[i];
                short_lig.gb_scale[i] = ligand.gb_hawkins_scale[i];
            }
        }

        // allocate born and distance matrices
        // calculate solvent accessible surface area for ligand

        clock_t         start;
        clock_t         stop;

        start = clock();
        sasa           *s = new sasa();
        //lig_sasa = s->getSASA(short_lig,0);

        //com_sasa = getSASA(short_com,receptor.num_atoms);
        // these are not passed to the class SASA_score
        float sasa_lig_tot, sasa_lig_philic, sasa_lig_phobic, sasa_lig_other,
              sasa_rec_tot, sasa_rec_philic, sasa_rec_phobic, sasa_rec_other;

        s->getSASA(short_lig, receptor.num_atoms, 
                                   lig_sasa_tot, lig_sasa_philic, lig_sasa_phobic, lig_sasa_other,
                                   sasa_lig_tot, sasa_lig_philic, sasa_lig_phobic, sasa_lig_other,
                                   sasa_rec_tot, sasa_rec_philic, sasa_rec_phobic, sasa_rec_other);
        stop = clock();

        //gnpol_lig = (lig_sasa * 0.00542) + 0.92;
        delete          s;
        if (verbose) {
            cout << "SASA=" << lig_sasa_tot << endl
                << "\tTime=" << ((long) stop -
                                 (long) start) /
                (float) CLOCKS_PER_SEC << "s" << endl;
        }
    }

    return true;

}

/************************************************/
void
SASA_score::prepare_complex(DOCKMol & ligand)
{
    int             i,
                    j;
    float           dist;
    vector < float >tmp;

    if (use_score) {

        short_com.allocate_short_arrays(ligand.num_atoms + receptor.num_atoms);

        // collect information about complex
        short_com.num_atoms = ligand.num_atoms + receptor.num_atoms;
        // receptor portion
        for (i = 0; i < receptor.num_atoms; i++) {
            if (receptor.atom_active_flags[i]) {

                short_com.atom_names[i] = receptor.atom_names[i];
                short_com.atom_types[i] = receptor.atom_types[i];

                short_com.coord[i].v[0] = receptor.x[i];
                short_com.coord[i].v[1] = receptor.y[i];
                short_com.coord[i].v[2] = receptor.z[i];

                short_com.charges[i] = receptor.charges[i];
                short_com.vdw_radius[i] = receptor.amber_at_radius[i];
                short_com.gb_radius[i] = receptor.gb_hawkins_radius[i];
                short_com.gb_scale[i] = receptor.gb_hawkins_scale[i];

            }
        }
        // ligand portion
        for (i = 0; i < ligand.num_atoms; i++) {
            if (ligand.atom_active_flags[i]) {

                short_com.atom_names[i + receptor.num_atoms] =
                    ligand.atom_names[i];
                short_com.atom_types[i + receptor.num_atoms] =
                    ligand.atom_types[i];

                short_com.coord[i + receptor.num_atoms].v[0] = ligand.x[i];
                short_com.coord[i + receptor.num_atoms].v[1] = ligand.y[i];
                short_com.coord[i + receptor.num_atoms].v[2] = ligand.z[i];

                short_com.charges[i + receptor.num_atoms] = ligand.charges[i];
                short_com.vdw_radius[i + receptor.num_atoms] =
                    ligand.amber_at_radius[i];
                short_com.gb_radius[i + receptor.num_atoms] =
                    ligand.gb_hawkins_radius[i];
                short_com.gb_scale[i + receptor.num_atoms] =
                    ligand.gb_hawkins_scale[i];


            }
        }

        // calculate solvent accessible surface area for complex

        clock_t         start;
        clock_t         stop;

        start = clock();
        sasa           *s = new sasa();

        //com_sasa = getSASA(short_com,receptor.num_atoms);
        s->getSASA(short_com, receptor.num_atoms,
                                   com_sasa_tot,     com_sasa_philic,     com_sasa_phobic,     com_sasa_other,
                                   com_sasa_lig_tot, com_sasa_lig_philic, com_sasa_lig_phobic, com_sasa_lig_other,
                                   com_sasa_rec_tot, com_sasa_rec_philic, com_sasa_rec_phobic, com_sasa_rec_other); 

        stop = clock();

        //gnpol_com = (com_sasa * 0.00542) + 0.92;
        
        delete          s;
        if (verbose) {
            cout << "SASA=" << com_sasa_tot << endl
                << "\tTime=" << ((long) stop -
                                 (long) start) /
                (float) CLOCKS_PER_SEC << "s" << endl;
        }
    }
}


