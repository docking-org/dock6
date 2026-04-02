// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 
// score_amber.cpp
// 
// implementation of class Amber_Score
// 
// This class defines continuum scoring for the AMBER energy.
// It is based on the other derived classes of Base_Score.
// The AMBER functionality actually relies on NAB code.
// 
// The initial version was written in 2004 by Scott Brozell while he was in
// the Dave Case group at TSRI with cooperation from Demetri Moustakas while
// he was in the Irwin "Tack" Kuntz group at UCSF.  During late 2005 the
// score24 protocol of Devleena (Mazumder) Shivakumar, a member of the Case
// group, was incorporated.  In addition in the first half of 2006,
// many minor modifications were made by Scott,
// with financial support from Tack, by Devleena, and with input from
// Terry (Downing) Lang, a member of the Kuntz group.
// In the last half of 2006, Scott and Terry added the distance movable region.
// In 2007 and 2008, incremental improvements and corrections were made by
// Scott while working on the DOCK6 paper with Terry at UCB.
// Amber_Score modifications made by Scott, mainly in early 2009 and early
// 2011, involved preparation and useability, not C++ coding changes.
// 
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 
// 
// This software is copyrighted, 2004-2018,
// by Scott R. Brozell and David A. Case. 
// 
// The authors hereby grant permission to use, copy, modify, and re-distribute
// this software and its documentation for any purpose, provided
// that existing copyright notices are retained in all copies and that this
// notice is included verbatim in any distributions. No written agreement,
// license, or royalty fee is required for any of the authorized uses.
// Modifications to this software may be distributed provided that
// the nature of the modifications are clearly indicated.
// 
// IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY
// FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
// ARISING OUT OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY
// DERIVATIVES THEREOF, EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
// 
// THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES,
// INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE
// IS PROVIDED ON AN "AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE
// NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
// MODIFICATIONS.
// 
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <assert.h>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>

#include "utils.h"
#include "dockmol.h"
#include "amber_typer.h"
#include "base_score.h"
#include "score_amber.h"
#include "sphere.h"

using namespace std;

// C linkage convention for NAB.
extern "C" {
    #include "nabcode.h"
    FILE *nabout;         // non-fatal-error, non-lex emissions from NAB.
    extern FILE *mmoout;  // lex emissions from NAB.
}

const Real AMBER_ELECTROSTATIC = 18.2223;
const char final_pdb_extension[] = ".final_pose.amber.pdb";
const char restart_extension[] = ".final_pose.amber.restart";

// +++++++++++++++++++++++++++++++++++++++++
Amber_Score::Amber_Score()
{
    // input_parameters and initialize should be called here, but
    // DOCK does not use constructors to create well defined objects.
    // Each class has an initialize member which is explicitly called.
    // DOCK could be redesigned to use C++ better, but for now we comply.

    complex = NULL;
    ligand = NULL;
    receptor = NULL;

    // the default for NAB output is stdout, and -o overrides stdout.
    nabout = stdout;
}


// +++++++++++++++++++++++++++++++++++++++++
Amber_Score::~Amber_Score()
{
    // freemolecule( NULL ) is safe.
    freemolecule(complex);
    freemolecule(ligand);
    freemolecule(receptor);
}


// +++++++++++++++++++++++++++++++++++++++++
void
Amber_Score::input_parameters(Parameter_Reader & parm, bool & primary_score,
                              bool & secondary_score)
{
    cout << "\nAmber Score Parameters\n"
         << "--------------------------------------------------------"
            "----------------------------------"
         << endl;

    use_primary_score = false;
    use_secondary_score = false;

    string          yesorno;
    if (!primary_score) {
        yesorno = parm.query_param("amber_score_primary", "no", "yes no");
        use_primary_score = (yesorno == "yes");
        primary_score = use_primary_score;
    }
    if (!secondary_score) {
        //yesorno = parm.query_param("amber_score_secondary", "no", "yes no");
	yesorno = "no";
        use_secondary_score = (yesorno == "yes");
        // temporarily disable amber_score_secondary until next release
        if ( use_secondary_score ) {
          cout<< "TEMPORARILY DEPRECATED OPTION!"
              << endl
              << "  Amber score as a secondary score is disabled."
              << endl
              << "  The recommended protocol is to perform two DOCK runs with"
              << endl
              << "  the second run specifying amber_score as the primary_score."
              << endl;
          exit(0);
        use_secondary_score = false;
        secondary_score = use_secondary_score;
        }
    }
    use_score = (use_primary_score || use_secondary_score);

    if (use_score) {
        receptor_file_prefix = parm.query_param(
                                   "amber_score_receptor_file_prefix", "rec");
        read_movable_region(parm);
        switch (movable) {
        case DISTANCE:
        case EVERYTHING:
        case LIGAND:
        case NAB_ATOM_EXPRESSION:
            istringstream(parm.query_param(
                "amber_score_minimization_rmsgrad", "0.01")) >>
                conv_criterion_rmsgrad;
            istringstream(parm.query_param(
                "amber_score_before_md_minimization_cycles", "100")) >>
                num_premd_min_cycles;
            if ( ! ( istringstream(parm.query_param("amber_score_md_steps",
                                                "3000")) >> num_md_steps ) ) {
                cout << "Error:  Invalid input for "
                     << "amber_score_md_steps !  "
                     << "Program will terminate."
                     << endl ;
                exit(0);
            }
            istringstream(parm.query_param(
                "amber_score_after_md_minimization_cycles", "100")) >>
                num_postmd_min_cycles;
        break;
        case NOTHING:
            num_premd_min_cycles = 0;
            num_md_steps = 0;
            num_postmd_min_cycles = 0;
        break;
        default:
            throw;
        break;
        }

        gb_model = parm.query_param("amber_score_gb_model", "5", "1 2 5");
        nonbonded_cutoff = parm.query_param("amber_score_nonbonded_cutoff",
                                            "18.0");
        temperature = parm.query_param("amber_score_temperature", "300.0");
        yesorno = parm.query_param("amber_score_abort_on_unprepped_ligand",
                                   "yes", "yes no");
        skip_unprepped = (yesorno == "no");
        verbose = 0 != parm.verbosity_level();   // -v is for extra scoring info
    }
}


// +++++++++++++++++++++++++++++++++++++++++
void
Amber_Score::read_movable_region(Parameter_Reader & parm)
{

    const string distance = "distance";
    const string everything = "everything";
    const string ligand = "ligand";
    const string nab_atom_expression = "nab_atom_expression";
    const string nothing = "nothing";
    const string space = " ";
    const string legal = everything + space + distance + space + ligand
                         + space + nab_atom_expression + space + nothing;
    string mov = parm.query_param("amber_score_movable_region", ligand, legal);
    if (mov == distance) {
        movable = DISTANCE;
        istringstream(parm.query_param("amber_score_movable_distance_cutoff",
                                       "3.0")) >> distance_cutoff;
        Active_Site_Spheres :: set_sphere_file_name( parm );
    } else if (mov == everything) {
        movable = EVERYTHING;
    } else if (mov == ligand) {
        movable = LIGAND;
    } else if (mov == nab_atom_expression) {
        movable = NAB_ATOM_EXPRESSION;
        // see NAB manual, pages 19-20 53-54 for NAB atom expressions.
        // Science-wise there should be no default;
        // but to be consistent with query_param we use a simple one.
        receptor_move_atomexpr =
            parm.query_param("amber_score_receptor_movable_atom_expr", "::");
        complex_move_atomexpr =
            parm.query_param("amber_score_complex_movable_atom_expr", "::");
    } else if (mov == nothing) {
        movable = NOTHING;
    } else {
        // programmer error since parm.query_param validates input
        throw;
    }
}


// +++++++++++++++++++++++++++++++++++++++++
void
Amber_Score::initialize(AMBER_TYPER & typer)
{
    if ( ! verbose ) {
        // send all NAB output to the black hole.
        nabout = fopen( "/dev/null", "a" );
        mmoout = nabout;
    }

    if (verbose) {
        cout<< "Initializing Amber_Score...\n\n"
            << "Reading the receptor input files.\n" ;
    }

    receptor_pdb = receptor_file_prefix + ".amber.pdb";
    receptor_prmtop = receptor_file_prefix + ".prmtop";
    receptor = getpdb(const_cast< char *>(receptor_pdb.c_str()), NULL);
    readparm(receptor, const_cast< char *>(receptor_prmtop.c_str()));
    num_receptor_atoms = receptor->m_prm->Natom;

    if (verbose) {
        cout<< "Number of receptor strands is " << receptor->m_nstrands << "\n"
            << "Number of receptor residues is " << receptor->m_prm->Nres <<"\n"
            << "Number of receptor atoms is " << num_receptor_atoms << "\n" ;
    }

    receptor_xyz.reserve(3 * num_receptor_atoms);
    setxyz_from_mol(&receptor, NULL, &receptor_xyz[0]);

    // Lipinski's rule of 500 Daltons
    int upperish_ligand_atoms = (500 + 100) / ((12 + 1) / 2);
    ligand_xyz.reserve(3 * upperish_ligand_atoms);
    complex_xyz.reserve(3 * (upperish_ligand_atoms + num_receptor_atoms));
    gradient_xyz.reserve(3 * (upperish_ligand_atoms + num_receptor_atoms));
    velocity_xyz.reserve(3 * (upperish_ligand_atoms + num_receptor_atoms));

    // Initialize the Molecular Mechanics options.
    // mm_options for minimizations
    nab_minimization_options = "diel=C, gbsa=0, nsnb=99999, cut="
        + nonbonded_cutoff + ", gb=" + gb_model;

    // mm_options for md's
    nab_md_options = "diel=C, dt=0.001, gamma_ln=2, rattle=0, cut="
        + nonbonded_cutoff
        + ", gb=" + gb_model
        + ", tempi=" + temperature + ", temp0=" + temperature;

    // mm_options for energy calculations on the final structures
    nab_energy_options = "cut=999.0, gbsa=1";

    if (verbose) {
        ostringstream min_frac;
        min_frac << 1;
        nab_minimization_options += ", ntpr=" + min_frac.str();
        nab_mme_initial_iteration = -1;
        ostringstream md_frac;
        md_frac << 1;
        nab_md_options += ", ntwx=" + md_frac.str();
        nab_md_options += ", ntpr=" + md_frac.str();
    } else {
        int const frac = 1;
        ostringstream min_frac;
        min_frac << max(num_premd_min_cycles / frac, 1);
        nab_minimization_options += ", ntpr=" + min_frac.str();
        nab_mme_initial_iteration = 1;
        ostringstream md_frac;
        md_frac << max(num_md_steps / frac, 1);
        nab_md_options += ", ntwx=" + md_frac.str();
        nab_md_options += ", ntpr=" + md_frac.str();
    }

    // initialize the NAB movable region
    switch (movable) {
    case DISTANCE: {
        int n = create_atomexpr_from_sphere_receptor_distance();
        cout<< "\n" << n << " receptor residues are within "
            << distance_cutoff << " Angstroms of the active site.\n" ;
        if (verbose) {
            cout<< "The distance based atom expression for the receptor is "
                << receptor_move_atomexpr
                << "\n" 
                << "The ligand must be the first strand.\n"
                << "The distance based atom expression for the complex is "
                << complex_move_atomexpr
                << "\n" ;
        }
    break; }
    case EVERYTHING:
        receptor_move_atomexpr = "::";
        complex_move_atomexpr = "::";
    break;
    case LIGAND: {
        receptor_move_atomexpr = "::ZZZZ";
        // assume the ligand is the first strand
        complex_move_atomexpr = "1::";
    break; }
    case NAB_ATOM_EXPRESSION: {
        // receptor_move_atomexpr has been initialized with the user input.
        // assume the ligand is the first strand
        // instead of creating the complex move atomexpr,
        // complex_move_atomexpr has been initialized with the user input.
        // this is only partly for programming simplicity; there might be
        // some utility in having complex_move_atomexpr not just
        // receptor_move_atomexpr plus a completely movable ligand.
    break; }
    case NOTHING:
        receptor_move_atomexpr = "::ZZZZ";
        complex_move_atomexpr = "::ZZZZ";
    break;
    default:
        throw;
    break;
    }

    char restraint_atomexpr[] = "::ZZZZ";        // no restraints
    string nab_allfree_atoms = "::";
    Real dummy[2];
    int dimensionr = 3 * receptor->m_prm->Natom;
    Real expected_drop_on_first_iter = 10.0;
    Real fret = 0.0;
    char* mv_a_e = const_cast< char* >( receptor_move_atomexpr.c_str() );
    Real rmsgrad = conv_criterion_rmsgrad;
    int status;
    int random_number = -1;
    if (verbose) {
        cout<< "\n************************************** \n"
            << "Computing the receptor energy:\n"
            << endl ;
    }
    // readparm( receptor, const_cast<char *>( receptor_prmtop.c_str() ) );
    switch (movable) {
    case DISTANCE:
    case EVERYTHING:
    case NAB_ATOM_EXPRESSION: {
    // the receptor protocol: minimize, md, minimize, and compute the energy 
        mm_options(const_cast< char *>(nab_minimization_options.c_str()));
        mme_init(receptor, mv_a_e, restraint_atomexpr, dummy, NULL);
        status = conjgrad(&receptor_xyz[0], &dimensionr, &fret, mme, &rmsgrad,
                          &expected_drop_on_first_iter, &num_premd_min_cycles);
        if (verbose) {
            cout<< ConjGradReturnCode( status );
        }
        // readparm( receptor, const_cast<char *>( receptor_prmtop.c_str() ) );
        rand2(&random_number);
        if (verbose) {
            cout<< "The random number is " << random_number << "\n" ;
        }
        mm_options(const_cast< char *>(nab_md_options.c_str()));
        status = md(dimensionr, num_md_steps, &receptor_xyz[0],
                    &gradient_xyz[0], &velocity_xyz[0], mme);
        // readparm( receptor, const_cast<char *>( receptor_prmtop.c_str() ) );
        mm_options(const_cast< char *>(nab_minimization_options.c_str()));
        status = conjgrad(&receptor_xyz[0], &dimensionr, &fret, mme, &rmsgrad,
                          &expected_drop_on_first_iter, &num_postmd_min_cycles);
        if (verbose) {
            cout<< ConjGradReturnCode( status );
        }
        readparm(receptor, const_cast< char *>(receptor_prmtop.c_str()));
        mm_options(const_cast< char *>(nab_energy_options.c_str()));
        mme_init(receptor, const_cast< char *>(nab_allfree_atoms.c_str()),
                 restraint_atomexpr, dummy, NULL);
        receptor_energy =
            mme(&receptor_xyz[0], &gradient_xyz[0], &nab_mme_initial_iteration);
        setmol_from_xyz(&receptor, NULL, &receptor_xyz[0]);
        string receptor_final_pdb = receptor_file_prefix + final_pdb_extension;
        putpdb(const_cast< char* >(receptor_final_pdb.c_str()), receptor,NULL);
        if (verbose) {
            string receptor_restart( receptor_file_prefix + restart_extension );
            char* fn = const_cast< char *>( receptor_restart.c_str() );
            putxv( & fn, & fn, & num_receptor_atoms, & fret, & receptor_xyz[0],
                   & velocity_xyz[0] );
        }
    break; }
    case LIGAND:
    case NOTHING:
    // the receptor protocol: compute the energy 
        mm_options(const_cast< char *>(nab_minimization_options.c_str()));
        mm_options(const_cast< char *>(nab_energy_options.c_str()));
        mme_init(receptor, const_cast< char *>(nab_allfree_atoms.c_str()),
                 restraint_atomexpr, dummy, NULL);
        receptor_energy =
            mme(&receptor_xyz[0], &gradient_xyz[0], &nab_mme_initial_iteration);
    break;
    default:
        throw;
    break;
    }

    if (verbose) {
        cout<< "Done computing the receptor energy.\n"
            << "Done initializing Amber_Score.\n"
            << endl ;
    }
}


// +++++++++++++++++++++++++++++++++++++++++
bool
Amber_Score::compute_score(DOCKMol & mol)
{
    if (use_score) {
        bool successful_import = import_dockmol_ligand(mol);
        if ( ! successful_import ) {
            mol.current_score = -MIN_FLOAT;  // arbitrarily large score
            mol.current_data = "Not scored because of a file I/O problem.\n";
            return false;
        }
        mol.current_score = calculate_amber_energy();
        mol.current_data = output_score_summary(mol.current_score);
        // send dock the ligand coordinates from the final complex structure.
        mol.setxyz( &complex_xyz[ 0 ] );
        // send dock the ligand charges from the NAB ligand molecule.
        mol.setcharges( ligand->m_prm->Charges, 1.0/AMBER_ELECTROSTATIC );
        freemolecule( ligand );
        ligand = NULL;
        freemolecule( complex );
        complex = NULL;
    }

    return true;
}


// +++++++++++++++++++++++++++++++++++++++++
string
Amber_Score::output_score_summary(float current_score)
{
    ostringstream text;

    if (use_score) {
        text << DELIMITER << setw(STRING_WIDTH) << "Amber_Score:"
             << setw(FLOAT_WIDTH) << fixed << current_score << endl
             << DELIMITER << setw(STRING_WIDTH) << "Amber_complex_energy:"
             << setw(FLOAT_WIDTH) << fixed << complex_energy << endl
             << DELIMITER << setw(STRING_WIDTH) << "Amber_receptor_energy:"
             << setw(FLOAT_WIDTH) << fixed << -receptor_energy << endl
             << DELIMITER << setw(STRING_WIDTH) << "Amber_ligand_energy:"
             << setw(FLOAT_WIDTH) << fixed << -ligand_energy << endl
        ;
    }
    return text.str();
}


// +++++++++++++++++++++++++++++++++++++++++
// Read the ligand and complex files that are external to DOCK
// and otherwise prepare for the score calculation.
// Normally, returns true on successful importing or terminates
// on any error; but if skipping of unprepped ligands has been
// requested then returns false on file opening errors.
bool
Amber_Score::import_dockmol_ligand(DOCKMol & lig)
{
    lig_file_prefix = lig.amber_score_ligand_id;
    string lig_inpcrd = lig_file_prefix + ".inpcrd";
    string lig_pdb = lig_file_prefix + ".amber.pdb";
    string lig_prmtop = lig_file_prefix + ".prmtop";
    complex_file_prefix = receptor_file_prefix + "." + lig_file_prefix;
    string complex_inpcrd = complex_file_prefix + ".inpcrd";
    string complex_pdb = complex_file_prefix + ".amber.pdb";
    complex_prmtop = complex_file_prefix + ".prmtop";

    if (verbose) {
        cout<< "\nReading the ligand " << lig_file_prefix << " input files.\n"
            << endl ;
    }
    if (skip_unprepped) {
        FILE *fp = fopen( const_cast< char *>(lig_pdb.c_str()), "r" );
        if( NULL == fp ) {
            perror( const_cast< char *>(("dock6: " + lig_pdb).c_str()) );
            cout<< "\nWarning cannot open " << lig_pdb
                << "\nSkipping this ligand as requested.\n" ;
            return false;
        }
        fclose( fp );
    }

    ligand = getpdb(const_cast< char *>(lig_pdb.c_str()), NULL);
    readparm(ligand, const_cast< char *>(lig_prmtop.c_str()));
    if (verbose) {
        cout<< "Number of ligand strands is " << ligand->m_nstrands << "\n"
            << "Number of ligand residues is " << ligand->m_prm->Nres <<"\n"
            << "Number of ligand atoms is " << ligand->m_prm->Natom << "\n" ;
    }
    assert( lig.num_atoms == ligand->m_prm->Natom );

    if (verbose) {
        cout<< "\nReading the complex " << complex_file_prefix
            << " input files.\n" ;
    }
    complex = getpdb(const_cast< char *>(complex_pdb.c_str()), NULL);
    readparm(complex, const_cast< char *>(complex_prmtop.c_str()));
    if (verbose) {
        cout<< "Number of complex strands is " << complex->m_nstrands << "\n"
            << "Number of complex residues is " << complex->m_prm->Nres <<"\n"
            << "Number of complex atoms is " << complex->m_prm->Natom << "\n" ;
    }
    assert( complex->m_prm->Natom == receptor->m_prm->Natom +
        ligand->m_prm->Natom );

    if (ligand_xyz.capacity() < 3 * lig.num_atoms) {
        ligand_xyz.reserve(3 * (lig.num_atoms));
        complex_xyz.reserve(3 * (num_receptor_atoms + lig.num_atoms));
        gradient_xyz.reserve(3 * (num_receptor_atoms + lig.num_atoms));
        velocity_xyz.reserve(3 * (num_receptor_atoms + lig.num_atoms));
    }
    setxyz_from_mol(&complex, NULL, &complex_xyz[0]);
    setxyz_from_mol(&ligand, NULL, &ligand_xyz[0]);

    if (false && verbose) {
        // TODO allocate only once
        MOLECULE_T*
            m = getpdb(const_cast< char *>(receptor_pdb.c_str()), NULL);
        setmol_from_xyz(&m, NULL, &complex_xyz[0]);
        // TODO need a dynamic filename
        //putpdb("importedligand.pdb", m, NULL);
        freemolecule(m);
    }

    return true;
}


// +++++++++++++++++++++++++++++++++++++++++
float
Amber_Score::calculate_amber_energy()
{
    char restraint_atomexpr[] = "::ZZZZ";        // no restraints
    string nab_allfree_atoms = "::";
    Real dummy[2];
    int dimensionc = 3 * complex->m_prm->Natom;
    int dimensionl = 3 * ligand->m_prm->Natom;
    Real expected_drop_on_first_iter = 10.0;
    Real fret;
    char* mv_a_e = const_cast< char *>(complex_move_atomexpr.c_str());
    Real rmsgrad = conv_criterion_rmsgrad;
    int status;
    // comments by Devleena: reset the random number generator. Useful when
    // scoring multiple ligands in the same DOCK run.
    int random_number = -1;

    if (verbose) {
        cout<< "\n************************************** \n"
            << "Computing the ligand energy:\n"
            << endl ;
    }
    switch (movable) {
    case DISTANCE:
    case EVERYTHING:
    case LIGAND:
    case NAB_ATOM_EXPRESSION: {
    // the ligand protocol: minimize, md, minimize, and compute the energy 
        // readparm( ligand, const_cast<char *>( ligand_prmtop.c_str() ) );
        mm_options(const_cast< char *>(nab_minimization_options.c_str()));
        mme_init(ligand, const_cast< char *>(nab_allfree_atoms.c_str()),
                 restraint_atomexpr, dummy, NULL);
        status = conjgrad(&ligand_xyz[0], &dimensionl, &fret, mme, &rmsgrad,
                          &expected_drop_on_first_iter, &num_premd_min_cycles);
        if (verbose) {
            cout<< ConjGradReturnCode( status );
        }
        // readparm( ligand, const_cast<char *>( ligand_prmtop.c_str() ) );
        rand2(&random_number);
        if (verbose) {
            cout<< "The random number is " << random_number << "\n" ;
        }
        mm_options(const_cast< char *>(nab_md_options.c_str()));
        status = md(dimensionl, num_md_steps, &ligand_xyz[0], &gradient_xyz[0],
                    &velocity_xyz[0], mme);
        // readparm( ligand, const_cast<char *>( ligand_prmtop.c_str() ) );
        mm_options(const_cast< char *>(nab_minimization_options.c_str()));
        status = conjgrad(&ligand_xyz[0], &dimensionl, &fret, mme, &rmsgrad,
                          &expected_drop_on_first_iter, &num_postmd_min_cycles);
        if (verbose) {
            cout<< ConjGradReturnCode( status );
        }
        // readparm( ligand, const_cast<char *>( ligand_prmtop.c_str() ) );
        mm_options(const_cast< char *>(nab_energy_options.c_str()));
        mme_init(ligand, const_cast< char *>(nab_allfree_atoms.c_str()),
                 restraint_atomexpr, dummy, NULL);
        ligand_energy =
            mme(&ligand_xyz[0], &gradient_xyz[0], &nab_mme_initial_iteration);
        setmol_from_xyz(&ligand, NULL, &ligand_xyz[0]);
        string lig_final_pdb = lig_file_prefix + final_pdb_extension;
        putpdb(const_cast< char* >(lig_final_pdb.c_str()), ligand, NULL);
        if (verbose) {
            string ligand_restart( lig_file_prefix + restart_extension );
            char* fn = const_cast< char *>( ligand_restart.c_str() );
            int num_ligand_atoms = ligand->m_prm->Natom;
            putxv( & fn, & fn, & num_ligand_atoms, & fret, & ligand_xyz[0],
                   & velocity_xyz[0] );
        }
    break; }
    case NOTHING:
    // the ligand protocol: compute the energy 
        mme_init(ligand, const_cast< char *>(nab_allfree_atoms.c_str()),
                 restraint_atomexpr, dummy, NULL);
        ligand_energy =
            mme(&ligand_xyz[0], &gradient_xyz[0], &nab_mme_initial_iteration);
    break;
    default:
        throw;
    break;
    }

    if (verbose) {
        cout<< "Done computing the ligand energy.\n\n"
            << endl ;
    }

    if (verbose) {
        cout<< "************************************** \n"
            << "Computing the complex energy:\n"
            << endl ;
    }
    switch (movable) {
    case DISTANCE:
    case EVERYTHING:
    case LIGAND:
    case NAB_ATOM_EXPRESSION: {
    // the complex protocol: minimize, md, minimize, and compute the energy 
        // readparm( complex, const_cast<char *>( complex_prmtop.c_str() ) );
        mm_options(const_cast< char *>(nab_minimization_options.c_str()));
        mme_init(complex, mv_a_e, restraint_atomexpr, dummy, NULL);
        status = conjgrad(&complex_xyz[0], &dimensionc, &fret, mme, &rmsgrad,
                          &expected_drop_on_first_iter, &num_premd_min_cycles);
        if (verbose) {
            cout<< ConjGradReturnCode( status );
        }
        // readparm( complex, const_cast<char *>( complex_prmtop.c_str() ) );
        rand2(&random_number);
        if (verbose) {
            cout<< "The random number is " << random_number << "\n" ;
        }
        mm_options(const_cast< char *>(nab_md_options.c_str()));
        status = md(dimensionc, num_md_steps, &complex_xyz[0], &gradient_xyz[0],
                    &velocity_xyz[0], mme);
        // readparm( complex, const_cast<char *>( complex_prmtop.c_str() ) );
        mm_options(const_cast< char *>(nab_minimization_options.c_str()));
        status = conjgrad(&complex_xyz[0], &dimensionc, &fret, mme, &rmsgrad,
                          &expected_drop_on_first_iter, &num_postmd_min_cycles);
        if (verbose) {
            cout<< ConjGradReturnCode( status );
        }
        readparm(complex, const_cast< char *>(complex_prmtop.c_str()));
        mm_options(const_cast< char *>(nab_energy_options.c_str()));
        mme_init(complex, const_cast< char *>(nab_allfree_atoms.c_str()),
                 restraint_atomexpr, dummy, NULL);
        complex_energy =
            mme(&complex_xyz[0], &gradient_xyz[0], &nab_mme_initial_iteration);
        setmol_from_xyz ( &complex, NULL, &complex_xyz[0] );
        // comments by Devleena: Write out the final pose of the complex
        string complex_final_pdb = complex_file_prefix + final_pdb_extension;
        putpdb(const_cast< char* >(complex_final_pdb.c_str()), complex, NULL);
        if (verbose) {
            string complex_restart( complex_file_prefix + restart_extension );
            char* fn = const_cast< char *>( complex_restart.c_str() );
            int num_complex_atoms = complex->m_prm->Natom;
            putxv( & fn, & fn, & num_complex_atoms, & fret, & complex_xyz[0],
                   & velocity_xyz[0] );
        }
    break; }
    case NOTHING:
    // the complex protocol: compute the energy 
        mme_init(complex, const_cast< char *>(nab_allfree_atoms.c_str()),
                 restraint_atomexpr, dummy, NULL);
        complex_energy =
            mme(&complex_xyz[0], &gradient_xyz[0], &nab_mme_initial_iteration);
    break;
    default:
        throw;
    break;
    }

    if (verbose) {
        cout<< "Done computing the complex energy.\n\n"
            << endl ;
    }

    if (verbose) {
        cout << "Energy of the complex " << complex_file_prefix << " is "
             << setprecision(12) << setw(18) << complex_energy << "\n"
             << "Energy of the receptor " << receptor_file_prefix << " is "
             << setprecision(12) << setw(18) << receptor_energy << "\n"
             << "Energy of the ligand " << lig_file_prefix << " is "
             << setprecision(12) << setw(18) << ligand_energy << "\n"
             << endl ;
    }

    return complex_energy - receptor_energy - ligand_energy;
}

// +++++++++++++++++++++++++++++++++++++++++
// Select residues by ligand-receptor distance.
// The ligand is represented by the active site sphere list.
// The receptor is represented by the NAB receptor MOLECULE_T.
// If any atom in a residue is within Amber_Score :: distance_cutoff
// Angstroms of the ligand then the whole residue is selected.
// The NAB atom expressions for the receptor and the complex are created.
// Return the number of selected residues.
int
Amber_Score::create_atomexpr_from_sphere_receptor_distance( )
{
    SphereVec active_site = Active_Site_Spheres :: get_instance();
    Active_Site_Spheres :: const_iterator i;
    typedef const char * Strand_Name;
    typedef int Strand_Sequence_Number;
    typedef int Residue_Sequence_Number;
    typedef map< Strand_Sequence_Number, set< Residue_Sequence_Number> >
        ResiduesPerStrand;
    ResiduesPerStrand residue_numbers;
    MOLECULE_T **m = & receptor;
    STRAND_T *sp;
    RESIDUE_T *res;
    ATOM_T *ap;
    int r, a, n;

    n = 0;
    // NAB idiom to loop over all atoms
    for( sp = (*m)->m_strands; sp; sp = sp->s_next ){
        set< Residue_Sequence_Number> rsnset;
        for( r = 0; r < sp->s_nresidues; r++ ){
            res = sp->s_residues[ r ];
            for( a = 0; a < res->r_natoms; a++ ){
                ap = &res->r_atoms[ a ];
                // STL idiom to loop over all spheres
                for ( i = active_site.begin(); i != active_site.end(); ++ i ) {
                    XYZCRD axyz = XYZCRD( ap->a_pos[ 0 ], ap->a_pos[ 1 ],
                        ap->a_pos[ 2 ] );
                    if ( i->crds.distance_squared( axyz ) <
                        distance_cutoff*distance_cutoff ) {
                        rsnset.insert( r + 1 );
                        ++ n;
                        goto NEXT_RESIDUE;  // exit sphere and atom loops
                    }
                }
            }
NEXT_RESIDUE:
        ;  // null statement
        }
        if ( ! rsnset.empty() ) {
            // strands generated from LEaP prmtops have positive integral names
            Strand_Sequence_Number ssn = atoi( sp->s_strandname ) ;
            if ( ssn <= 0 ) {
                cout << "Error:  Invalid strand sequence number !  "
                     << "Contact dock-fans@docking.org \n"
                     << "Program will terminate."
                     << endl ;
                exit(0);
            }
            residue_numbers.insert( ResiduesPerStrand::value_type(
                ssn, rsnset ) );
        }
    }

    ostringstream rmae;
    ResiduesPerStrand :: const_iterator j;
    for ( j = residue_numbers.begin(); j != residue_numbers.end(); ++j ) {
        rmae << j->first;  // add the strand label
        rmae << ":";       // terminate the strand section
        copy( residue_numbers[ j->first ].begin(), 
            residue_numbers[ j->first ].end(),
            ostream_iterator< Residue_Sequence_Number >(rmae, ",") );
        rmae << ":";  // terminate the residue labels section
        rmae << "|";  // prepare for next strand with the alteration token
    }
    if ( residue_numbers.empty() ) {
        receptor_move_atomexpr = "::ZZZZ";
    }
    else {
        receptor_move_atomexpr = rmae.str();
        // remove the trailing alteration token
        receptor_move_atomexpr.erase( receptor_move_atomexpr.size() - 1, 1 );
    }

    ostringstream mae;
    // add all ligand atoms
    // assume the ligand is the first strand
    mae << "1::";
    for ( j = residue_numbers.begin(); j != residue_numbers.end(); ++j ) {
        mae << "|";
        mae << j->first + 1 ;  // add the strand label
        mae << ":";            // terminate the strand section
        copy( residue_numbers[ j->first ].begin(), 
            residue_numbers[ j->first ].end(),
            ostream_iterator< Residue_Sequence_Number >(mae, ",") );
        mae << ":";  // terminate the residue labels section
    }
    complex_move_atomexpr = mae.str();

    return n;
}


// +++++++++++++++++++++++++++++++++++++++++
// validate and store the return code from NAB's conjgrad.
ConjGradReturnCode::ConjGradReturnCode( int cod )
: code( cod )
{
    assert( code > 0 || ( code <= -1 && code >= -4 ) );
}

// +++++++++++++++++++++++++++++++++++++++++
int
ConjGradReturnCode::getcode( )
const
{
    return code;
}

// +++++++++++++++++++++++++++++++++++++++++
std::ostream &
operator<<( std::ostream & os, ConjGradReturnCode const & cg )
{
    if ( cg.getcode() > 0 ) {
        cout<< "Minimization converged in " 
            << cg.getcode()
            << " iterations.\n" ;
    } else if ( cg.getcode() == -1 ) {
        cout<< "Minimization terminated due to a bad line search.\n" ;
    } else if ( cg.getcode() == -2 ) {
        cout<< "Minimization terminated due to an uphill search direction.\n" ;
    } else if ( cg.getcode() == -3 ) {
        cout<< "Minimization reached the maximum number of iterations.\n" ;
    } else if ( cg.getcode() == -4 ) {
        cout<< "Minimization terminated due to a nondecreasing function "
               "value.\n" ;
    }
    return os ;
}

