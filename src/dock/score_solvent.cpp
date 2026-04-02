#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <time.h>
#include "amber_typer.h"
#include "sasa.h"
#include "score_solvent.h"
#include "trace.h"


using namespace std;


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 
// Abbreviated DOCKMol--needed to loop over complexes in solvation scoring.
//

/*************************************/
ShortDOCKMol::ShortDOCKMol()
{

    coord = NULL;
    charges = NULL;
    atom_names = NULL;
    atom_types = NULL;

    vdw_radius = NULL;

    gb_radius = NULL;
    gb_scale = NULL;
}

/*************************************/
ShortDOCKMol::~ShortDOCKMol()
{

    delete[]coord;
    delete[]charges;
    delete[]atom_names;
    delete[]atom_types;

    delete[]vdw_radius;

    delete[]gb_radius;
    delete[]gb_scale;
}

/*************************************/
void
ShortDOCKMol::allocate_short_arrays(int natoms)
{

    num_atoms = natoms;

    delete[]coord;
    coord = new POINT[num_atoms];

    delete[]charges;
    charges = new float[num_atoms];

    delete[]atom_names;
    atom_names = new string[num_atoms];

    delete[]atom_types;
    atom_types = new string[num_atoms];

    delete[]vdw_radius;
    vdw_radius = new float[num_atoms];

    delete[]gb_radius;
    gb_radius = new float[num_atoms];

    delete[]gb_scale;
    gb_scale = new float[num_atoms];

    // assign init values
    for (int i = 0; i < num_atoms; i++) {

        atom_names[i] = "";
        atom_types[i] = "";
        charges[i] = 0.0;

        vdw_radius[i] = 0.0;

        gb_radius[i] = 0.0;
        gb_scale[i] = 0.0;
    }
}


/************************************************/
// Copy constructor to convert DockMol to ShortDockMol object
ShortDOCKMol::ShortDOCKMol(DOCKMol const& mol) {

    num_atoms = mol.num_atoms;

    //copy constructor needs to allocate memory
    coord = new POINT[num_atoms];
    charges = new float[num_atoms];
    atom_names = new string[num_atoms];
    atom_types = new string[num_atoms];
    vdw_radius = new float[num_atoms];
    gb_radius = new float[num_atoms];
    gb_scale = new float[num_atoms];

    for (int i = 0; i < mol.num_atoms; i++) {

        if (mol.atom_active_flags[i]) {
            atom_names[i] = mol.atom_names[i];
            atom_types[i] = mol.atom_types[i];
            
            //cout << "###" << i << "###" << mol.atom_types[i] << endl;

            coord[i].v[0] = mol.x[i];
            coord[i].v[1] = mol.y[i];
            coord[i].v[2] = mol.z[i];

            charges[i] = mol.charges[i];
            vdw_radius[i] = mol.amber_at_radius[i];
            gb_radius[i] = mol.gb_hawkins_radius[i];
            gb_scale[i] = mol.gb_hawkins_scale[i];
        }
    }
}


/************************************************/
// Combine two Dockmol objects to create a complex
ShortDOCKMol::ShortDOCKMol(ShortDOCKMol const& mol1,ShortDOCKMol const& mol2) {
    num_atoms = mol1.num_atoms + mol2.num_atoms;

    //copy constructor needs to allocate memory
    coord = new POINT[num_atoms];
    charges = new float[num_atoms];
    atom_names = new string[num_atoms];
    atom_types = new string[num_atoms];
    vdw_radius = new float[num_atoms];
    gb_radius = new float[num_atoms];
    gb_scale = new float[num_atoms];

    int count = 0;
    for (int i = 0; i < mol1.num_atoms; i++) {

            atom_names[count] = mol1.atom_names[i];
            atom_types[count] = mol1.atom_types[i];

            coord[count].v[0] = mol1.coord[i].v[0];  //x
            coord[count].v[1] = mol1.coord[i].v[1];  //y
            coord[count].v[2] = mol1.coord[i].v[2];  //z

            charges[count] = mol1.charges[i];
            vdw_radius[count] = mol1.vdw_radius[i];
            gb_radius[count] = mol1.gb_radius[i];
            gb_scale[count] = mol1.gb_scale[i];
            count++;
    }

    //count is now at i+1, we can start appending
    for (int i = 0; i < mol2.num_atoms; i++) {

            atom_names[count] = mol2.atom_names[i];
            atom_types[count] = mol2.atom_types[i];

            coord[count].v[0] = mol2.coord[i].v[0];
            coord[count].v[1] = mol2.coord[i].v[1];
            coord[count].v[2] = mol2.coord[i].v[2];

            charges[count] = mol2.charges[i];
            vdw_radius[count] = mol2.vdw_radius[i];
            gb_radius[count] = mol2.gb_radius[i];
            gb_scale[count] = mol2.gb_scale[i];
            count++;
    }
    //assert( count == num_atoms );
}

/************************************************/
// this is copy constructor that takes two DOCKMols and 
// makes a ShortDOCKMol
/************************************************/
ShortDOCKMol::ShortDOCKMol(DOCKMol const& mol1,DOCKMol const& mol2) {
   ShortDOCKMol smol1(mol1); //make a shortdockmol object 
   ShortDOCKMol smol2(mol2); //make a shortdockmol object 
   ShortDOCKMol smol_combined(smol1,smol2);
   for (int i = 0; i < smol_combined.num_atoms; i++) {
        
            atom_names[i] = smol_combined.atom_names[i];
            atom_types[i] = smol_combined.atom_types[i];

            coord[i].v[0] = smol_combined.coord[i].v[0];
            coord[i].v[1] = smol_combined.coord[i].v[1];
            coord[i].v[2] = smol_combined.coord[i].v[2];

            charges[i] = smol_combined.charges[i];
            vdw_radius[i] = smol_combined.vdw_radius[i];
            gb_radius[i] = smol_combined.gb_radius[i];
            gb_scale[i] = smol_combined.gb_scale[i];
    }
}

/************************************************/
/************************************************/
/************************************************/
// XIAOQIN ZOU GBSA
/************************************************/
float
square_distance(float a[3], float b[3])
{
    return
        (a[0] - b[0]) * (a[0] - b[0]) +
        (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]);
}

/************************************************/
GB_Pairwise::GB_Pairwise()
{
    gb_grid = NULL;
    sa_grid = NULL;
}

/************************************************/
GB_Pairwise::~GB_Pairwise()
{
    delete gb_grid;
    delete sa_grid;
}

/************************************************/
void
GB_Pairwise::input_parameters(Parameter_Reader & parm, bool & primary_score,
                              bool & secondary_score)
{
    string          tmp;

    cout << "\nZou GB/SA Score Parameters\n"
            "--------------------------------------------------------"
            "----------------------------------"
         << endl;

    use_primary_score = false;
    use_secondary_score = false;

    if (!primary_score) {
        tmp = parm.query_param("gbsa_zou_score_primary", "no", "yes no");
        use_primary_score = (tmp == "yes");
        primary_score = use_primary_score;
    }

    if (!secondary_score) {
        //tmp = parm.query_param("gbsa_zou_score_secondary", "no", "yes no");
        tmp = "no";
	use_secondary_score = (tmp == "yes");
        secondary_score = use_secondary_score;
    }

    use_score = (use_primary_score || use_secondary_score);

    if (use_score) {
        gb_grid_prefix = parm.query_param("gbsa_zou_gb_grid_prefix", "gb_grid");
        sa_grid_prefix = parm.query_param("gbsa_zou_sa_grid_prefix", "sa_grid");
        vdw_grid_prefix = parm.query_param("gbsa_zou_vdw_grid_prefix", "grid");
        screen_file = parm.query_param("gbsa_zou_screen_file", "screen.in");
        tmp = parm.query_param("gbsa_zou_solvent_dielectric", "78.300003");
        sscanf(tmp.c_str(), "%f", &solvent_dielectric);
    }

    verbose = 0 != parm.verbosity_level();
}

/************************************************/
void
GB_Pairwise::initialize(AMBER_TYPER & typer)
{
    if (use_score) {
        cout << "Initializing Zou GB/SA Score Routines..." << endl;
        
        vdw_score.use_score = true;
        vdw_score.es_scale = 0;
        vdw_score.vdw_scale = 1;
        vdw_score.rep_radius_scale = 1;
        vdw_score.grid_file_name = vdw_grid_prefix;
        vdw_score.initialize(typer);
        
        //gb_grid = GB_Grid :: get_instance(gb_grid_prefix);
        //sa_grid = SA_Grid :: get_instance(sa_grid_prefix);
        gb_grid = new GB_Grid();
        gb_grid->get_instance(gb_grid_prefix);
        sa_grid = new SA_Grid();
        sa_grid->get_instance(sa_grid_prefix);

        gb_initialized = false;

        cout << "Done Initializing Zou GB/SA Score Routines." << endl;
    }
}

/************************************************/
bool
GB_Pairwise::compute_score(DOCKMol & mol)
{
    float           total;

    if (use_score) {
        //cout << "I AM HERE"<<endl;
        get_gb_solvation_score(mol);
        //cout << "I AM HERE"<<endl;
        get_sa_solvation_score(mol);
        //cout << "I AM HERE"<<endl;

        if (vdw_score.compute_score(mol))
            vdw_component = vdw_score.vdw_component;
        else
            return false;

        total =
            gb_component + 0.6 * vdw_component - 0.015 * sa_component +
            0.025 * del_SAS_hp;

        mol.current_score = total;
        mol.current_data = output_score_summary(mol);

        return true;
    } else
        return false;
}

/************************************************/
string
GB_Pairwise::output_score_summary(DOCKMol & mol)
{
    ostringstream text;

    if (use_score) {
        text << DELIMITER << setw(STRING_WIDTH) << "Zou_GBSA_Score:"
             << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl
             << DELIMITER << setw(STRING_WIDTH) << "Zou_GBSA_vdw_energy:"
             << setw(FLOAT_WIDTH) << fixed << 0.6 * vdw_component << endl
             << DELIMITER << setw(STRING_WIDTH) << "Zou_GBSA_gb_energy:"
             << setw(FLOAT_WIDTH) << fixed << gb_component << endl
             << DELIMITER << setw(STRING_WIDTH) << "Zou_GBSA_sa_energy:"
             << setw(FLOAT_WIDTH) << fixed << -0.015 * sa_component << endl
             << DELIMITER << setw(STRING_WIDTH) << "Zou_GBSA_delta_SAS_hp_energy:"
             << setw(FLOAT_WIDTH) << fixed << 0.025 * del_SAS_hp << endl;

        if (use_internal_energy)
            text << DELIMITER << setw(STRING_WIDTH) << "Internal_energy_repulsive:"
                 << setw(FLOAT_WIDTH) << fixed << mol.internal_energy << endl;

    }
    return text.str();
}
/************************************************/
float
GB_Pairwise::get_gb_solvation_score(DOCKMol & mol)
{
    float           factor,
                    sub_total,
                    solv_lp,
                    es,
                    temp;
    float           crd[3];
    float           rvdweff0,
                    sum;
    float           r2,
                    rij2,
                    inv_aij2,
                    dij,
                    f_gb;
    int             nearpt[3];
    int             i,
                    j;
    float           INVA_NEG = (float) 0.1;
    float           ri,
                    sj,
                    sj2,
                    rij,
                    RLij,
                    RUij,
                    RLij2,
                    RUij2,
                    hij;
    float           screen[25];
    float           ligand_solvation,
                    ligand_solvation_partial;
    float           receptor_solvation,
                    receptor_solvation_partial;
    float           ligand_desolvation,
                    receptor_desolvation,
                    screened_es;
    FILE           *fp1;

    if ( ! use_score )
        return 0;

    gb_initialized = false;

    // read in the screen.in file
    fp1 = fopen(screen_file.c_str(), "r");

    if (fp1 == NULL) {
        cout << "\n\nCould not open " << screen_file <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    for (i = 0; i < 25; i++) {
        fscanf(fp1, "%f", &screen[i]);
    }
    fclose(fp1);

    // set the solvent dielectric factor (hardcoded for now)
    factor = -166.0 * (1.0 - 1.0 / solvent_dielectric);
    sub_total = 0.0;

    // Conditional- initialize GB if necessary
    if (!gb_initialized) {
        gb_initialized = true;

        // Clear & resize ligand GB arrays
        gb_lig_coords.clear();
        gb_lig_inv_a.clear();
        gb_lig_radius_eff.clear();
        gb_lig_sum_lig.clear();
        gb_lig_flag.clear();

        gb_lig_coords.resize(mol.num_atoms);
        gb_lig_inv_a.resize(mol.num_atoms);
        gb_lig_radius_eff.resize(mol.num_atoms);
        gb_lig_sum_lig.resize(mol.num_atoms);
        gb_lig_flag.resize(mol.num_atoms);

        // Initialize ligand GB arrays
        for (i = 0; i < mol.num_atoms; i++) {
            gb_lig_coords[i].crd[0] = 0.0;
            gb_lig_coords[i].crd[1] = 0.0;
            gb_lig_coords[i].crd[2] = 0.0;
            gb_lig_inv_a[i] = 0.0;
            gb_lig_radius_eff[i] = 0.0;
            gb_lig_sum_lig[i] = 0.0;
            gb_lig_flag[i] = 0;
        }

        // ID ligand atoms->grid points
        for (i = 0; i < mol.num_atoms; i++) {
            if (mol.amber_at_well_depth[i] != 0.0) {

                // Label atoms that fall outside the grid
                crd[0] = mol.x[i] - gb_grid->gb_grid_origin[0];
                crd[1] = mol.y[i] - gb_grid->gb_grid_origin[1];
                crd[2] = mol.z[i] - gb_grid->gb_grid_origin[2];

                for (j = 0; j < 3; j++) {
                    nearpt[j] = NINT(crd[j] / gb_grid->gb_grid_spacing);
                    if ((nearpt[j] < 0) | (nearpt[j] >= gb_grid->gb_grid_span[j])) {
                        gb_lig_flag[i] = 1;
                    }
                    gb_lig_coords[i].crd[j] = crd[j];
                }

                rvdweff0 = mol.amber_at_radius[i] + VDWOFF;
                gb_lig_radius_eff[i] = rvdweff0;

            }
        }

        // Find solvation energy for pure ligand
        for (i = 0; (i < mol.num_atoms) && (gb_lig_flag[i] == 0); i++) {
            if (mol.amber_at_well_depth[i] != 0.0) {

                rvdweff0 = gb_lig_radius_eff[i];
                sum = 0.0;
                ri = rvdweff0 - VDWOFF;

                // Pairwise calculation of Born Radius
                for (j = 0;
                     (j < mol.num_atoms) && (gb_lig_flag[j] == 0)
                     && (j != i); j++) {
                    if (mol.amber_at_well_depth[j] != 0.0) {
                        rvdweff0 = gb_lig_radius_eff[j];
                        sj = (rvdweff0 -
                              VDWOFF) * screen[mol.amber_at_id[j] + 1];
                        sj2 = sj * sj;
                        r2 = square_distance(gb_lig_coords[i].crd,
                                             gb_lig_coords[j].crd);
                        rij = sqrt(r2);
                        if (ri < rij + sj) {
                            RLij = 1. / MAX(ri, rij - sj);
                            RUij = 1. / (rij + sj);
                            RLij2 = RLij * RLij;
                            RUij2 = RUij * RUij;
                            hij = RLij - RUij + 0.25 * (sj2 / rij - rij) *
                                 (RLij2 - RUij2) - 0.5 / rij * log(RLij / RUij);
                        } else {
                            hij = 0;
                        }

                        sum += hij;
                    }
                }           // End pairwise born radii

                gb_lig_sum_lig[i] = sum;
                // jpcb v108 p5455 eq 21:
                gb_lig_inv_a[i] = 1.0 / (ri) - 0.5 * sum;
                if (gb_lig_inv_a[i] < 0) {
                    if (verbose)
                        cout << "Pure ligand inv_a error: " << i << "\t" <<
                            gb_lig_inv_a[i] << endl;
                    gb_lig_inv_a[i] = INVA_NEG;
                }
            }
        }

        ligand_solvation = 0.0;

        for (i = 0; (i < mol.num_atoms) && (gb_lig_flag[i] == 0); i++) {
            for (j = 0; (j < mol.num_atoms) && (gb_lig_flag[j] == 0); j++) {
                rij2 =
                    square_distance(gb_lig_coords[i].crd,
                                    gb_lig_coords[j].crd);
                if (rij2 <= gb_grid->gb_grid_cutsq) {
                    inv_aij2 = gb_lig_inv_a[i] * gb_lig_inv_a[j];
                    dij = rij2 * inv_aij2 / gb_grid->gb_grid_f_scale;
                    f_gb = sqrt(rij2 + exp(-dij) / inv_aij2);
                    ligand_solvation +=
                        factor * mol.charges[i] * mol.charges[j] / f_gb;
                }
            }
        }                   // End Ligand solvation energy

        // Calculate effect of receptor atoms on ligand atom Born radii
        for (i = 0; (i < mol.num_atoms) && (gb_lig_flag[i] == 0); i++) {
            if (mol.amber_at_well_depth[i] != 0.0) {
                rvdweff0 = gb_lig_radius_eff[i];
                sum = gb_lig_sum_lig[i];
                ri = rvdweff0 - VDWOFF;

                for (j = 0; j < gb_grid->gb_rec_total_atoms; j++) {
                    rvdweff0 = gb_grid->gb_rec_radius_eff[j];
                    sj = (rvdweff0 - VDWOFF) * gb_grid->gb_rec_receptor_screen[j];
                    sj2 = sj * sj;
                    r2 = square_distance(gb_lig_coords[i].crd,
                                         gb_grid->gb_rec_coords[j].crd);
                    rij = sqrt(r2);
                    if (ri < rij + sj) {
                        RLij = 1. / MAX(ri, rij - sj);
                        RUij = 1. / (rij + sj);
                        RLij2 = RLij * RLij;
                        RUij2 = RUij * RUij;
                        hij = RLij - RUij + 0.25 * (sj2 / rij - rij) *
                              (RLij2 - RUij2) - 0.5 / rij * log(RLij / RUij);
                    } else {
                        hij = 0;
                    }

                    sum += hij;
                }           // End receptor loop

                gb_lig_inv_a[i] = 1.0 / (ri) - 0.5 * sum;

                if (gb_lig_inv_a[i] < 0) {
                    if (verbose)
                        cout << "Receptor on ligand inv_a error: " << i <<
                            "\t" << gb_lig_inv_a[i] << endl;
                    gb_lig_inv_a[i] = INVA_NEG;
                }
            }
        }                   // End receptor atom -> ligand born radii
                            // calculation

        // Calculate effect of ligand atoms on receptor atom born radii
        for (i = 0; i < gb_grid->gb_rec_total_atoms; i++) {
            rvdweff0 = gb_grid->gb_rec_radius_eff[i];
            sum = 0.0;

            ri = rvdweff0 - VDWOFF;
            for (j = 0; (j < mol.num_atoms) && (gb_lig_flag[j] == 0); j++) {
                rvdweff0 = gb_lig_radius_eff[j];
                sj = (rvdweff0 - VDWOFF) * screen[mol.amber_at_id[j] + 1];
                sj2 = sj * sj;
                r2 = square_distance(gb_grid->gb_rec_coords[i].crd,
                                     gb_lig_coords[j].crd);
                rij = sqrt(r2);

                if (ri < rij + sj) {
                    RLij = 1. / MAX(ri, rij - sj);
                    RUij = 1. / (rij + sj);
                    RLij2 = RLij * RLij;
                    RUij2 = RUij * RUij;
                    hij = RLij - RUij + 0.25 * (sj2 / rij - rij) *
                          (RLij2 - RUij2) - 0.5 / rij * log(RLij / RUij);
                } else {
                    hij = 0;
                }

                sum += hij;
            }

            gb_grid->gb_rec_inv_a[i] = gb_grid->gb_rec_inv_a_rec[i] - sum * 0.5;

            if (gb_grid->gb_rec_inv_a[i] < 0) {
                if (verbose)
                    cout << "Ligand on receptor inva error: " << i << "\t"
                        << gb_grid->gb_rec_inv_a[i] << endl;
                gb_grid->gb_rec_inv_a[i] = INVA_NEG;
            }

        }                   // End ligand atom -> receptor born radii
                            // calculation

        // Polarization solvation energy of L-P complex
        // ll term in solv_complex
        ligand_solvation_partial = 0.0;
        for (i = 0; (i < mol.num_atoms) && (gb_lig_flag[i] == 0); i++) {
            for (j = 0; (j < mol.num_atoms) && (gb_lig_flag[j] == 0); j++) {
                rij2 =
                    square_distance(gb_lig_coords[i].crd,
                                    gb_lig_coords[j].crd);

                if (rij2 <= gb_grid->gb_grid_cutsq) {
                    inv_aij2 = gb_lig_inv_a[i] * gb_lig_inv_a[j];
                    dij = rij2 * inv_aij2 / gb_grid->gb_grid_f_scale;
                    f_gb = sqrt(rij2 + exp(-dij) / inv_aij2);
                    ligand_solvation_partial +=
                        factor * mol.charges[i] * mol.charges[j] / f_gb;
                }
            }
        }

        // pp term in solv_complex
        receptor_solvation_partial = 0.0;
        for (i = 0; i < gb_grid->gb_rec_total_atoms; i++) {
            for (j = 0; j < gb_grid->gb_rec_total_atoms; j++) {
                rij2 =
                    square_distance(gb_grid->gb_rec_coords[i].crd,
                                    gb_grid->gb_rec_coords[j].crd);

                if (rij2 <= gb_grid->gb_grid_cutsq) {
                    inv_aij2 = gb_grid->gb_rec_inv_a[i] * gb_grid->gb_rec_inv_a[j];
                    dij = rij2 * inv_aij2 / gb_grid->gb_grid_f_scale;
                    f_gb = sqrt(rij2 + exp(-dij) / inv_aij2);

                    receptor_solvation_partial +=
                        factor * gb_grid->gb_rec_charge[i] * gb_grid->gb_rec_charge[j] / f_gb;

                }
            }
        }

        receptor_solvation = gb_grid->gb_rec_solv_rec;

    } else {                // end initialization routine

        for (i = 0; (i < mol.num_atoms) && (gb_lig_flag[i] == 0); i++) {
            if (mol.amber_at_well_depth[i] != 0.0) {
                for (j = 0; j < 3; j++) {
                    crd[j] = gb_lig_coords[i].crd[j] - gb_grid->gb_grid_origin[j];
                    nearpt[j] = NINT(crd[j] / gb_grid->gb_grid_spacing);

                    if ((nearpt[j] < 0) | (nearpt[j] >= gb_grid->gb_grid_span[j])) {
                        gb_lig_flag[i] = 1;
                    }

                    gb_lig_coords[i].crd[j] = crd[j];
                }
            }
        }

    }                       // end all initialization routines

    // find the atomic contribution to screened es interactions
    solv_lp = 0.0;
    es = 0.0;

    for (i = 0; (i < mol.num_atoms) && (gb_lig_flag[i] == 0); i++) {
        temp = 0.0;

        for (j = 0; j < gb_grid->gb_rec_total_atoms; j++) {
            rij2 =
                square_distance(gb_lig_coords[i].crd, gb_grid->gb_rec_coords[j].crd);

            if (rij2 < gb_grid->gb_grid_cutsq) {
                inv_aij2 = gb_lig_inv_a[i] * gb_grid->gb_rec_inv_a[j];
                dij = rij2 * inv_aij2 / gb_grid->gb_grid_f_scale;
                f_gb = sqrt(rij2 + exp(-dij) / inv_aij2);
                solv_lp +=
                    2 * factor * mol.charges[i] * gb_grid->gb_rec_charge[j] / f_gb;
                es += 332 * mol.charges[i] * gb_grid->gb_rec_charge[j] / sqrt(rij2);
                temp +=
                    2 * factor * mol.charges[i] * gb_grid->gb_rec_charge[j] / f_gb +
                    332 * mol.charges[i] * gb_grid->gb_rec_charge[j] / sqrt(rij2);
            }
        }
    }

    ligand_desolvation = ligand_solvation_partial - ligand_solvation;
    receptor_desolvation = receptor_solvation_partial - receptor_solvation;
    screened_es = solv_lp + es;
    gb_component = ligand_desolvation + receptor_desolvation + screened_es;

    if (verbose) {
        // setw applies only to the next emission, but fixed applies to all
        ios_base::fmtflags save_flags = cout.flags();
        cout << "\nGB Calculation\n"
             << "--------------------"
             << "\nligand desolvation:\t\t"
             << setw(FLOAT_WIDTH) << fixed << ligand_desolvation
             << "\nreceptor desolvation:\t\t"
             << setw(FLOAT_WIDTH) << fixed << receptor_desolvation
             << "\nscreened electrostatics:\t"
             << setw(FLOAT_WIDTH) << fixed << screened_es
             << "\ntotal:\t\t\t\t"
             << setw(FLOAT_WIDTH) << fixed << gb_component << endl;
        cout.flags( save_flags );  // restore flags
    }

    return gb_component;
}

/************************************************/
float
GB_Pairwise::get_sa_solvation_score(DOCKMol & mol)
{
    int             i,
                    j,
                    k,
                    l,
                    n,
                    i1,
                    mark;
    float           crd[3],
                    crd0[3],
                    crd1[3],
                    r2;
    int             vdwn,
                    nsum;
    int             nearpt[3];
    int             neighbor_num[MAX_ATOM_LIG],
                    neighbors[MAX_ATOM_LIG][MAX_ATOM_LIG];
    int             neighbor_num1[MAX_ATOM_REC+MAX_ATOM_LIG],
                    neighbors1[MAX_ATOM_REC+MAX_ATOM_LIG][MAX_ATOM_LIG];
    int             nsas_lig,
                    delta_nsas_rec,
                    delta_nsas_lig;
    int             nsas_hp_lig,
                    delta_nsas_hp_rec,
                    delta_nsas_hp_lig;
    int             nsas_pol_lig,
                    delta_nsas_pol_rec,
                    delta_nsas_pol_lig;
    float           score[2],
                    factor;

    // conider making neighbor_num, neighbors, neighbor_num1, neighbors1 allocated dynamicly
    // int* neighbor_num1;
    // int** neighbors1;
    // neighbor_num1 = new [MAX_ATOM_LIG] int

    if (use_score) {
        // clear all ligand parameters
        sa_grid_mark_sas_lig.clear();
        sa_lig_coords.clear();
        sa_lig_radius_eff.clear();

        // resize all ligand parameters
        sa_grid_mark_sas_lig.resize(mol.num_atoms);
        sa_lig_coords.resize(mol.num_atoms);
        sa_lig_radius_eff.resize(mol.num_atoms);

        // initialize ligand parameters
        for (i = 0; i < mol.num_atoms; i++)
            for (j = 0; j < 2000; j++)
                sa_grid_mark_sas_lig[i].vals[j] = 0;


        for (i = 0; i < mol.num_atoms; i++) {
            sa_lig_coords[i].crd[0] = 0.0;
            sa_lig_coords[i].crd[1] = 0.0;
            sa_lig_coords[i].crd[2] = 0.0;
            sa_lig_radius_eff[i] = 0.0;
            neighbor_num[i] = 0;
            neighbor_num1[i] = 0;
            for (j = 0; j < 150; j++) {
                neighbors[i][j] = 0;
                neighbors1[i][j] = 0;
            }
        }

        //cout << "I AM HERE in SA (1)"<<endl;
        //cout << "mol.num_atoms = " << mol.num_atoms <<endl;
        //cout << "sa_grid->sa_rec_total_atoms = " << sa_grid->sa_rec_total_atoms <<endl;
        if (MAX_ATOM_REC< sa_grid->sa_rec_total_atoms){
            cout << "Error in get_sa_solvation_score :: MAX_ATOM_REC< sa_grid->sa_rec_total_atoms :: "<< MAX_ATOM_REC << "<" << sa_grid->sa_rec_total_atoms << endl;
        } 
        for (i = 0; i < sa_grid->sa_rec_total_atoms; i++) {
            n = i + mol.num_atoms;
            neighbor_num1[n] = 0;
            for (j = 0; j < 150; j++)
                neighbors1[n][j] = 0;
        }

        //cout << "I AM HERE in SA (2)"<<endl;
        // Pure ligand atoms
        for (i = 0; i < mol.num_atoms; i++) {
            if (mol.amber_at_well_depth[i] != 0.0) {

                crd[0] = mol.x[i] - sa_grid->sa_grid_origin[0];
                crd[1] = mol.y[i] - sa_grid->sa_grid_origin[1];
                crd[2] = mol.z[i] - sa_grid->sa_grid_origin[2];

                for (j = 0; j < 3; j++) {
                    nearpt[j] = NINT(crd[j] / sa_grid->sa_grid_spacing);
                    sa_lig_coords[i].crd[j] = crd[j];
                }

                r2 = SQR(mol.amber_at_radius[i] + sa_grid->sa_grid_r_probe);
                sa_lig_radius_eff[i] = r2;
            }
        }

        // find neighboring atoms
        for (i = 0; i < mol.num_atoms; i++) {
            if (mol.amber_at_well_depth[i] != 0.0) {

                for (k = 0; k < 3; k++)
                    crd[k] = sa_lig_coords[i].crd[k];

                for (j = i + 1; j < mol.num_atoms; j++) {
                    if (mol.amber_at_well_depth[j] != 0.0) {

                        for (k = 0; k < 3; k++)
                            crd1[k] = sa_lig_coords[j].crd[k];

                        r2 = square_distance(crd, crd1);
                        if (r2 < sa_grid->sa_grid_r2_cutoff) {
                            neighbor_num[i]++;
                            neighbor_num[j]++;
                            neighbors[i][neighbor_num[i]] = j;
                            neighbors[j][neighbor_num[j]] = i;
                        }
                    }
                }

                for (j = 0; j < sa_grid->sa_rec_total_atoms; j++) {
                    for (k = 0; k < 3; k++)
                        crd1[k] = sa_grid->sa_rec_coords[j].crd[k];

                    n = j + mol.num_atoms;
                    r2 = square_distance(crd, crd1);
                    if (r2 < sa_grid->sa_grid_r2_cutoff) {
                        neighbor_num1[i]++;
                        neighbor_num1[n]++;
                        neighbors1[i][neighbor_num1[i]] = n;
                        neighbors1[n][neighbor_num1[n]] = i;
                    }
                }

            }
        }

        // Find SASA for pure ligand
        nsas_lig = 0;
        nsas_hp_lig = 0;
        nsas_pol_lig = 0;

        for (i = 0; i < mol.num_atoms; i++) {
            if (mol.amber_at_well_depth[i] != 0.0) {
                vdwn = mol.amber_at_id[i];
                nsum = 0;

                for (k = 0; k < 3; k++)
                    crd0[k] = sa_lig_coords[i].crd[k];

                for (l = 0; l < sa_grid->sa_rec_nsphgrid[vdwn]; l++) {
                    for (k = 0; k < 3; k++)
                        crd[k] = crd0[k] + sa_grid->sa_grid_sphgrid_crd[vdwn][l].crd[k];

                    sa_grid_mark_sas_lig[i].vals[l] = 0;
                    mark = 0;

                    for (j = 1; j <= neighbor_num[i]; j++) {
                        n = neighbors[i][j];

                        for (k = 0; k < 3; k++)
                            crd1[k] = sa_lig_coords[n].crd[k];

                        r2 = square_distance(crd, crd1);
                        if (r2 < sa_lig_radius_eff[n]) {
                            mark = 1;
                            break;
                        }
                    }

                    if (mark == 0) {
                        sa_grid_mark_sas_lig[i].vals[l] = 1;
                        nsas_lig++;
                        nsum++;
                    }
                }

                if (sa_grid->sa_rec_nhp[vdwn] == 1)
                    nsas_hp_lig += nsum;
                else
                    nsas_pol_lig += nsum;
            }
        }


        // Effect of receptor atoms on SAS of ligand atoms
        delta_nsas_lig = 0;
        delta_nsas_hp_lig = 0;
        delta_nsas_pol_lig = 0;

        for (i = 0; i < mol.num_atoms; i++) {
            if (mol.amber_at_well_depth[i] != 0.0) {
                nsum = 0;
                vdwn = mol.amber_at_id[i];

                for (k = 0; k < 3; k++)
                    crd0[k] = sa_lig_coords[i].crd[k];

                for (l = 0; l < sa_grid->sa_rec_nsphgrid[vdwn]; l++) {
                    if (sa_grid_mark_sas_lig[i].vals[l] == 1) {

                        for (k = 0; k < 3; k++)
                            crd[k] =
                                crd0[k] + sa_grid->sa_grid_sphgrid_crd[vdwn][l].crd[k];

                        for (j = 1; j <= neighbor_num1[i]; j++) {
                            n = neighbors1[i][j] - mol.num_atoms;

                            for (k = 0; k < 3; k++)
                                crd1[k] = sa_grid->sa_rec_coords[n].crd[k];

                            r2 = square_distance(crd, crd1);
                            if (r2 < sa_grid->sa_rec_radius_eff[n]) {
                                delta_nsas_lig++;
                                nsum++;
                                break;
                            }
                        }
                    }
                }

                if (sa_grid->sa_rec_nhp[vdwn] == 1)
                    delta_nsas_hp_lig += nsum;
                else
                    delta_nsas_pol_lig += nsum;
            }
        }

        // Effect of ligand atoms on SAS of receptor atoms
        delta_nsas_rec = 0;
        delta_nsas_hp_rec = 0;
        delta_nsas_pol_rec = 0;

        for (i = 0; i < sa_grid->sa_rec_total_atoms; i++) {
            nsum = 0;
            vdwn = sa_grid->sa_rec_vdwn_rec[i] - 1;

            for (k = 0; k < 3; k++)
                crd0[k] = sa_grid->sa_rec_coords[i].crd[k];

            i1 = i + mol.num_atoms;

            for (l = 0; l < sa_grid->sa_rec_nsphgrid[vdwn]; l++) {
                if (sa_grid->sa_grid_mark_sas[i][l] == 1) {
                    for (k = 0; k < 3; k++)
                        crd[k] = crd0[k] + sa_grid->sa_grid_sphgrid_crd[vdwn][l].crd[k];

                    for (j = 1; j <= neighbor_num1[i1]; j++) {
                        n = neighbors1[i1][j];

                        for (k = 0; k < 3; k++)
                            crd1[k] = sa_lig_coords[n].crd[k];

                        r2 = square_distance(crd, crd1);
                        if (r2 < sa_lig_radius_eff[n]) {
                            delta_nsas_rec++;
                            nsum++;
                            break;
                        }

                    }
                }
            }

            if (sa_grid->sa_rec_nhp[vdwn] == 1)
                delta_nsas_hp_rec += nsum;
            else
                delta_nsas_pol_rec += nsum;
        }

        factor = SQR(sa_grid->sa_grid_sspacing);
        score[0] = -(float) (delta_nsas_rec) * factor;
        score[1] = -(float) (delta_nsas_lig) * factor;

        sa_component = score[0] + score[1];
        del_SAS_hp = (-(delta_nsas_hp_rec + delta_nsas_hp_lig) * factor);

        if (verbose) {
            ios_base::fmtflags save_flags = cout.flags();
            cout << "\nSA Calculation:\n"
                 << "--------------------\n"
                 << fixed
                 << "SAS_rec = " << (sa_grid->sa_rec_nsas * factor)
                 << ", delta_SAS_rec = " << score[0] << "\n"
                 << "SAS_lig = " << (nsas_lig * factor)
                 << ", delta_SAS_lig = " << score[1] << "\n"
                 << "delta_SAS_pol = "
                 << (-(delta_nsas_pol_rec + delta_nsas_pol_lig) * factor)
                 << ", delta_SAS_hp = "
                 << (-(delta_nsas_hp_rec + delta_nsas_hp_lig) * factor) << "\n"
                 << "total = " << (score[0] + score[1]) << endl << endl;
            cout.flags( save_flags );  // restore flags
        }

        return sa_component;
    } else
        return 0;

}

/************************************************/
/************************************************/
/************************************************/
// HAWKINS, CRAMER, TRUHLAR GBSA--GB PORTION
/************************************************/
float
distance2(POINT coord1, POINT coord2)
{

    float           sum = 0.0;
    for (int i = 0; i < 3; i++) {
        sum += ((coord1.v[i] - coord2.v[i])*(coord1.v[i] - coord2.v[i]));
    }
    return sum;
}

/************************************************/
float
born_radii_calc(float Rkkp, float Pkp, float Pk)
{
    float           Lkkp,
                    Ukkp,
                    iLkkp,
                    iLkkp2,
                    iUkkp,
                    iUkkp2;
    float           sum1,
                    sum2,
                    sum3,
                    sum4;

    // Calculate Lkkp
    Lkkp = 0.0;

    if ((Rkkp + Pkp) <= Pk)
        Lkkp = 1.0;
    if ((Rkkp - Pkp) <= Pk)
        Lkkp = Pk;
    if (Pk < (Rkkp + Pkp))
        Lkkp = Pk;
    if (Pk <= (Rkkp - Pkp))
        Lkkp = (Rkkp - Pkp);

    // calculate Ukkp
    Ukkp = 0.0;

    if ((Rkkp + Pkp) <= Pk)
        Ukkp = 1.0;
    if (Pk < (Rkkp + Pkp))
        Ukkp = (Rkkp + Pkp);

    // compute some inverse numbers
    iLkkp = 1.0 / Lkkp;
    iLkkp2 = (iLkkp*iLkkp);

    iUkkp = 1.0 / Ukkp;
    iUkkp2 = (iUkkp*iUkkp);

    sum1 = iLkkp - iUkkp;
    sum2 = ((Rkkp / 4.0) * (iUkkp2 - iLkkp2));
    sum3 = ((1.0 / (2.0 * Rkkp)) * (log(Lkkp / Ukkp)));
    sum4 = (((Pkp * Pkp) / (4.0 * Rkkp)) * (iLkkp2 - iUkkp2));

    return sum1 + sum2 + sum3 + sum4;
}

/************************************************/
GB_Hawkins::GB_Hawkins()
{

    born_radius_rec = NULL;
    rec_dist_mat = NULL;
    born_radius_lig = NULL;
    lig_dist_mat = NULL;
    born_radius_com = NULL;
    com_dist_mat = NULL;

}

/************************************************/
GB_Hawkins::~GB_Hawkins()
{

    delete[]born_radius_rec;
    delete[]rec_dist_mat;
    delete[]born_radius_lig;
    delete[]lig_dist_mat;
    delete[]born_radius_com;
    delete[]com_dist_mat;
}

/************************************************/
void
GB_Hawkins::born_array_calc(ShortDOCKMol & mol, float offset, float *dist_mat,
                            float *born_array)
{

    int             i,
                    j;
    float           sum,
                    tmp_dist;
    float           Rkkp,
                    ialphak,
                    Pk,
                    Pkp;

    for (i = 0; i < mol.num_atoms; i++) {

        sum = 0.0;
        Pk = mol.gb_radius[i] - offset;

        for (j = 0; j < mol.num_atoms; j++) {

            if (i != j) {

                // calc atom pair distance & init variables
                tmp_dist = distance2(mol.coord[i], mol.coord[j]);
                dist_mat[i * mol.num_atoms + j] = tmp_dist;

                Rkkp = sqrt(tmp_dist);
                Pkp = mol.gb_scale[j] * (mol.gb_radius[j] - offset);

                sum = sum + born_radii_calc(Rkkp, Pkp, Pk);
            }
        }                       // end second atom loop

        sum = sum / 2.0;
        ialphak = (1.0 / Pk) - sum;
        born_array[i] = (1.0 / ialphak);

        if (born_array[i] <= 0.0) {
            born_array[i] = -MIN_FLOAT;  // an arbitrarily large number
        }

    }                           // end first atom loop

}

/************************************************/
float
GB_Hawkins::gpol_calc(ShortDOCKMol & mol, float *dist_mat, float *born,
                      bool salt, float kappa, float espout, float gfac)
{

    int             i,
                    j;

    float           dist,
                    totgpol,
                    gpself,
                    gpshld,
                    alphaij,
                    alphaij2;
    float           Dij,
                    fgb,
                    engshld,
                    engself,
                    charge_sum,
                    screen;

    // calculate gpol
    totgpol = 0.0;
    gpself = 0.0;
    gpshld = 0.0;

    for (i = 0; i < mol.num_atoms; i++) {

        for (j = i + 1; j < mol.num_atoms; j++) {

            if (i != j) {

                dist = dist_mat[i * mol.num_atoms + j];

                alphaij2 = born[i] * born[j];
                alphaij = sqrt(alphaij2);

                Dij = dist / (4.0 * alphaij2);
                fgb = sqrt(dist + (alphaij2 * exp(-Dij)));

                if (!salt)
                    engshld = (gfac * mol.charges[i] * mol.charges[j]) / fgb;
                else {
                    charge_sum = (mol.charges[i] * mol.charges[j]) / fgb;
                    screen = (exp(-kappa * fgb)) / espout;
                    engshld = -166 * (charge_sum - (screen * charge_sum));
                }
                gpshld = gpshld + engshld;


            }                   // end calculation inner loop

        }                       // end second atom loop

        if (!salt)
            engself = (gfac * mol.charges[i] * mol.charges[i]) / born[i];
        else {
            charge_sum = (mol.charges[i] * mol.charges[i]) / born[i];
            screen = (exp(-kappa * born[i])) / espout;
            engself = -166 * (charge_sum - (screen * charge_sum));
        }
        gpself = gpself + engself;

    }                           // end first atom loop

    gpshld = gpshld * 2;
    totgpol = gpself + gpshld;

    if (verbose)
        cout << "\tGPSELF=" << gpself << "\tGPSHLD=" << gpshld << "\t";

    return totgpol;

}

/************************************************/
void
GB_Hawkins::input_parameters(Parameter_Reader & parm, bool & primary_score,
                             bool & secondary_score)
{
    string          tmp;

    use_primary_score = false;
    use_secondary_score = false;

    cout << "\nHawkins GB/SA Score Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    if (!primary_score) {
        tmp = parm.query_param("gbsa_hawkins_score_primary", "no", "yes no");
        if (tmp == "yes")
            use_primary_score = true;
        else
            use_primary_score = false;

        primary_score = use_primary_score;
    }

    if (secondary_score) { 
        cout << "Warning. Secondary Score is no longer uses. set secondary_score to false. " << endl;
        secondary_score = false;
    }
    //if (!secondary_score) {
        //tmp = parm.query_param("gbsa_hawkins_score_secondary", "no", "yes no");
        //tmp == "no";
        //if (tmp == "yes")
        //    use_secondary_score = true;
        //else
        //    use_secondary_score = false;

        //secondary_score = use_secondary_score;
    //}

    if (use_primary_score || use_secondary_score)
        use_score = true;
    else
        use_score = false;

    if (use_score) {
        rec_filename =
            parm.query_param("gbsa_hawkins_score_rec_filename",
                             "receptor.mol2");

        espout =
            atof(parm.
                 query_param("gbsa_hawkins_score_solvent_dielectric",
                             "78.5").c_str());
        if (espout <= 0.0) {
            cout <<
                "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }
        salt_screen =
            (parm.query_param("gbsa_hawkins_use_salt_screen", "no", "yes no") ==
             "yes") ? true : false;
        if (salt_screen) {
            salt_conc =
                atof(parm.query_param("gbsa_hawkins_score_salt_conc(M)", "0.0").
                     c_str());
            if (salt_conc <= 0.0) {
                cout <<
                    "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }
        }
        gb_offset =
            atof(parm.query_param("gbsa_hawkins_score_gb_offset", "0.09").
                 c_str());
        if (gb_offset <= 0.0) {
            cout <<
                "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }
        tmp =
            parm.query_param("gbsa_hawkins_score_cont_vdw_and_es", "yes",
                             "yes no");
        if (tmp == "yes") {
            cont_vdw = true;
            att_exp =
                atoi(parm.query_param("gbsa_hawkins_score_vdw_att_exp", "6").
                     c_str());
            if (att_exp <= 0) {
                cout <<
                    "ERROR:  Parameter must be an integer greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }
            rep_exp =
                atoi(parm.query_param("gbsa_hawkins_score_vdw_rep_exp", "12").
                     c_str());
            if (rep_exp <= 0) {
                cout <<
                    "ERROR:  Parameter must be an integer greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }
            rep_radius_scale = atof(parm.query_param("grid_score_rep_rad_scale", "1").c_str());   
        } else {
            cont_vdw = false;
            vdw_score.grid_file_name =
                parm.query_param("gbsa_hawkins_score_grid_prefix", "grid");
            cout << "   WARNING:  Grid file should have been prepared with" <<
                endl;
            cout <<
                "   \"distance_dielectric = no\" and \"dielectric_factor = 1\""
                << endl;
        }

    }
    verbose = 0 != parm.verbosity_level();

}

/************************************************/
void
GB_Hawkins::initialize(AMBER_TYPER & typer)
{
    ifstream        rec_file;
    bool            read_vdw,
                    flexible,
                    use_chem,
                    use_ph4,
                    use_volume;

    use_score = true;

    if (use_score) {
        cout << "Initializing Hawkins GB/SA Score Routines..." << endl;

        // calculate scaling factors for gb
        if (salt_screen) {
            kappa = sqrt((8.3943 * salt_conc) / espout);
            kappa = 0.73 * kappa;
        } else
            gfac = -166.0 * (1.0 - (1.0 / espout));

        // initialize VDW for continuous vdw
        if (cont_vdw) {
            init_vdw_energy(typer, att_exp, rep_exp);
        } else {                // initialize VDW for grided vdw
            vdw_score.use_score = true;
            vdw_score.es_scale = 1;
            vdw_score.vdw_scale = 1;
            vdw_score.rep_radius_scale = 1;
            vdw_score.initialize(typer);
        }

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

        read_vdw = true;
        use_chem = false;
        use_ph4  = false;
        use_volume = false;
        typer.prepare_molecule(receptor, read_vdw, use_chem, use_ph4, use_volume);

        // calculate receptor alone

        prepare_receptor();
    }

}

/************************************************/
bool
GB_Hawkins::compute_score(DOCKMol & mol)
{

    if (use_score) {

        prepare_ligand(mol);
        prepare_complex(mol);

        // Print full GBSA score
        double          dgnpol = gnpol_com - (gnpol_rec + gnpol_lig);
        double          dgpol = gpol_com - (gpol_rec + gpol_lig);

        if (!cont_vdw) {
            if (vdw_score.compute_score(mol)) {
                
                vdw_component = vdw_score.vdw_component;
                es_component = vdw_score.es_component;
            } else {
                return false;
            }
        }

        total_score = vdw_component + es_component + dgnpol + dgpol;

        mol.current_score = total_score;
        mol.current_data = output_score_summary(mol);

        //cout << mol.current_data<<endl;

    }

    return true;
}

/************************************************/
string
GB_Hawkins::output_score_summary(DOCKMol & mol)
{
    // Print full GBSA score
    double          dgnpol = gnpol_com - (gnpol_rec + gnpol_lig);
    double          dgpol = gpol_com - (gpol_rec + gpol_lig);

    ostringstream text;

    if (use_score) {
        text << DELIMITER << setw(STRING_WIDTH) << "Hawkins_GBSA_Score:"
             << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl
             << DELIMITER << setw(STRING_WIDTH) << "Hawkins_GBSA_vdw_energy:"
             << setw(FLOAT_WIDTH) << fixed << vdw_component << endl
             << DELIMITER << setw(STRING_WIDTH) << "Hawkins_GBSA_es_energy:"
             << setw(FLOAT_WIDTH) << fixed << es_component << endl
             << DELIMITER << setw(STRING_WIDTH) << "Hawkins_GBSA_gb_energy:"
             << setw(FLOAT_WIDTH) << fixed << dgpol << endl
             << DELIMITER << setw(STRING_WIDTH) << "Hawkins_GBSA_sa_energy:"
             << setw(FLOAT_WIDTH) << fixed << dgnpol << endl ;

        if (use_internal_energy)
            text << DELIMITER << setw(STRING_WIDTH) << "Internal_energy_repulsive:"
                 << setw(FLOAT_WIDTH) << fixed << mol.internal_energy << endl;

    }
    return text.str();
}

/************************************************/
void
GB_Hawkins::prepare_receptor()
{

    int             i,
                    j;

    gpol_rec = 0.0;

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

        // allocate born and distance matrices
        delete[]born_radius_rec;
        born_radius_rec = new float[short_rec.num_atoms];

        delete[]rec_dist_mat;
        rec_dist_mat = new float[short_rec.num_atoms * short_rec.num_atoms];

        // calculate born radii
        born_array_calc(short_rec, gb_offset, rec_dist_mat, born_radius_rec);

        // calculate gpol for receptor
        if (verbose) {
            cout << "\n----------------------------------------------------\n";
            cout << "RECEPTOR\n";
        }
        gpol_rec =
            gpol_calc(short_rec, rec_dist_mat, born_radius_rec, salt_screen,
                      kappa, espout, gfac);

        // calculate solvent accessible surface area of receptor

        // Measure SASA time
        clock_t         start;
        clock_t         stop;

        start = clock();
        sasa           *s = new sasa();
        //rec_sasa = s->getSASA(short_rec);
        rec_sasa = s->getSASA(short_rec,receptor.num_atoms);
        stop = clock();

        gnpol_rec = (rec_sasa * 0.00542) + 0.92;
        delete          s;
        if (verbose) {
            // cout.precision(3);
            cout << "SASA=" << rec_sasa << endl
                << "\tTotal:\tGB=" << gpol_rec << "\tSA=" << gnpol_rec
                << "\tTime=" << ((long) stop -
                                 (long) start) /
                (float) CLOCKS_PER_SEC << "s" << endl;
        }

    }
}

/************************************************/
bool
GB_Hawkins::prepare_ligand(DOCKMol & ligand)
{
    int             i,
                    j;
    vector < float >tmp;

    gpol_lig = 0.0;

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
        delete[]born_radius_lig;
        born_radius_lig = new float[short_lig.num_atoms];

        delete[]lig_dist_mat;
        lig_dist_mat = NULL;
        lig_dist_mat = new float[short_lig.num_atoms * short_lig.num_atoms];

        // calculate born radii
        born_array_calc(short_lig, gb_offset, lig_dist_mat, born_radius_lig);

        // calculate gpol for ligand
        if (verbose)
            cout << "LIGAND\n";
        gpol_lig =
            gpol_calc(short_lig, lig_dist_mat, born_radius_lig, salt_screen,
                      kappa, espout, gfac);

        // calculate solvent accessible surface area for ligand

        clock_t         start;
        clock_t         stop;

        start = clock();
        sasa           *s = new sasa();
        //lig_sasa = s->getSASA(short_lig);
        lig_sasa = s->getSASA(short_lig,0);
        stop = clock();

        gnpol_lig = (lig_sasa * 0.00542) + 0.92;
        delete          s;
        if (verbose) {
            cout << "SASA=" << lig_sasa << endl
                << "\tTotal:\tGB=" << gpol_lig << "\tSA=" << gnpol_lig
                << "\tTime=" << ((long) stop -
                                 (long) start) /
                (float) CLOCKS_PER_SEC << "s" << endl;
        }
    }

    return true;

}

/************************************************/
void
GB_Hawkins::prepare_complex(DOCKMol & ligand)
{
    int             i,
                    j;
    float           dist;
    vector < float >tmp;

    gpol_com = 0.0;
    vdw_component = 0.0;
    es_component = 0.0;

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

                // calculate vdw and electrostatic values for continuous
                if (cont_vdw) {
                    for (j = 0; j < ligand.num_atoms; j++) {
                        dist =
                            sqrt(((ligand.x[j] - receptor.x[i])*(ligand.x[j] - receptor.x[i])) +
                                 ((ligand.y[j] - receptor.y[i])*(ligand.y[j] - receptor.y[i])) + 
                                 ((ligand.z[j] - receptor.z[i])*(ligand.z[j] - receptor.z[i])));
                        vdw_component +=
                            ((vdwA[ligand.amber_at_id[j]] *
                              vdwA[receptor.amber_at_id[i]]) /
                             pow(dist,
                                 rep_exp)) -
                            ((vdwB[ligand.amber_at_id[j]] *
                              vdwB[receptor.amber_at_id[i]]) /
                             pow(dist, att_exp));
                        es_component +=
                            ((332 * receptor.charges[i] * ligand.charges[j]) /
                             dist);
                        // for gbsa, gas phase electrostatics should be used
                        // as in q1*q2/r not q1*q2/r^2
                    }
                }
            }
        }
        // ligand portion
        for (i = 0; i < ligand.num_atoms; i++) {
            if (ligand.atom_active_flags[i]) {

                short_com.atom_names[i + receptor.num_atoms] =
                    ligand.atom_names[i];

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



        // allocate born and distance matrices
        delete[]born_radius_com;
        born_radius_com = new float[short_com.num_atoms];

        delete[]com_dist_mat;
        com_dist_mat = new float[short_com.num_atoms * short_com.num_atoms];

        // calculate born radii
        born_array_calc(short_com, gb_offset, com_dist_mat, born_radius_com);

        // calculate gpol for complex
        if (verbose)
            cout << "COMPLEX\n";
        gpol_com =
            gpol_calc(short_com, com_dist_mat, born_radius_com, salt_screen,
                      kappa, espout, gfac);

        // calculate solvent accessible surface area for complex

        clock_t         start;
        clock_t         stop;

        start = clock();
        sasa           *s = new sasa();
        //com_sasa = s->getSASA(short_com);
        com_sasa = s->getSASA(short_com,receptor.num_atoms);
        stop = clock();

        gnpol_com = (com_sasa * 0.00542) + 0.92;
        delete          s;
        if (verbose) {
            cout << "SASA=" << com_sasa << endl
                << "\tTotal:\tGB=" << gpol_com << "\tSA=" << gnpol_com
                << "\tTime=" << ((long) stop -
                                 (long) start) /
                (float) CLOCKS_PER_SEC << "s" << endl;
        }
    }
}

