// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// dock.cpp
//
// definition of the main function
//
// The pre version 6.4 history is documented in the manual.
// The version 6.4 changes were primarily from SUNY Stony Brook:
// Code modifications made by the group of
// Robert C Rizzo (Stony Brook University) with
// Sudipto Mukherjee (Stony Brook University), and
// Trent E Balius (Stony Brook University);
// as well as Demetri Moustakas.
// They include:
// rep_vdw ligand clash filter, lig_int_nrg in mol2,
// multi-anchor growth tree,
// torsion pre-minimizer during growth of active atoms bugfix,
// clear anchor_positions array in next_anchor(), question tree update,
// DM's orienting bugfix (updated), last anchor correction, and MW & Charge.
//
// The lists of other past and present authors is in the manual.
//
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//
// This software is copyrighted, 2004-2023,
// by the DOCK Developers.
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

#include <time.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdint.h>
#include <map>
#include "stdlib.h"
#include "string.h"
#include "amber_typer.h"
#include "master_conf.h"
#include "library_file.h"
#include "master_score.h"
#include "orient.h"
#include "simplex.h"
#include "trace.h"
#include "utils.h"
#include "version.h"
#include "filter.h"
//#include <gzstream.h>
#include "gzstream/gzstream.h"

#ifndef __APPLE__
#include <sys/sysinfo.h>
#endif


// Temporary support on cygwin for calculate_simulation_time vs wall_clock_seconds
// https://stackoverflow.com/questions/71324957/cygwin-c-standard-library-does-not-support-time-utc-and-timespec-get
#if __cplusplus > 199711L
  #define TIME_PRECISION
#endif
#ifdef __CYGWIN__
  #undef TIME_PRECISION
#endif

void            dock_new_handler();
void            print_header( bool USE_MPI, int processes );
//LEP - Preprocesor directive determines if the compiler is >CPP11 for certain functionality
double          wall_clock_seconds();
#ifdef TIME_PRECISION
inline double   calculate_simulation_time(timespec t_start, timespec t_end);
#endif
using namespace std;


/************************************************/
int
main(int argc, char **argv)
{

#ifndef __APPLE__
    // set up memory info
    struct sysinfo memInfo;
    sysinfo (&memInfo);
#endif


    // set the function that will be called if new fails.
    std::set_new_handler(dock_new_handler);

    // synchronize C++ stream io and C stdio.
    ios::sync_with_stdio();

#ifdef TRACE
    Trace::traceOn();
#else
    Trace::traceOff();
#endif
    Trace trace( "::main" );

    // DOCK Classes
    // DOCK does not use constructors, but each class has an initialize member.
    // Note that the compiler-generated default constructors effectively
    // do nothing; so that these objects should be considered uninitialized
    // until their initialize members are called below.
    Parameter_Reader         c_parm;
    Library_File             c_library;
    Orient                   c_orient;
    Master_Conformer_Search  c_master_conf;
    Bump_Filter              c_bmp_score;
    Simplex_Minimizer        c_simplex;
    Master_Score             c_master_score;
    AMBER_TYPER              c_typer;
    Filter                   c_filter;
    DOCKMol                  mol;

    // Local vars
    ofstream        outfile;
    streambuf      *original_cout_buffer;
    char            fname[500];

    bool            USE_MPI;
    // Preprocessor macro MPI is controlled by the platform configuration
    // which is specified during installation.
#ifdef BUILD_DOCK_WITH_MPI
    USE_MPI = true;
#else
    USE_MPI = false;
#endif

    // Initialize mpi and determine if proper number of processors has been called.
    if (USE_MPI){
        // mpi_init requires pointers to argc and argv.
        USE_MPI = c_library.initialize_mpi(&argc, &argv);
    }

    // Direct the output to a file or stdout
    if (check_commandline_argument(argv, argc, "-o") != -1) {
        if (USE_MPI) {
            if (c_library.rank > 0) {
                sprintf(fname, "%s.%d",
                        parse_commandline_argument(argv, argc, "-o").c_str(),
                        c_library.rank);
            } else {
                sprintf(fname, "%s",
                        parse_commandline_argument(argv, argc, "-o").c_str());
            }
        } else {
            sprintf(fname, "%s",
                    parse_commandline_argument(argv, argc, "-o").c_str());
        }
        // Redirect non error output to the -o file.
        // C stdio.
        freopen( fname, "a", stdout );
        // C++ stream io.
        outfile.open(fname);
        original_cout_buffer = cout.rdbuf();
        cout.rdbuf(outfile.rdbuf());

    } else if (USE_MPI) {
        cout << "Error: DOCK must be run with the -o outfile option under MPI"
             << endl;
        c_library.finalize_mpi();
        exit(0);
    }


    // /////////////////////////////////////////////////////////////////////////
    // Begin timing
#ifndef TIME_PRECISION
    double          start_time = wall_clock_seconds();
    //int             wall_clock_nseconds();
#endif
#ifdef TIME_PRECISION
    struct timespec start_time;
    timespec_get(&start_time, TIME_UTC);
#endif
    print_header( USE_MPI, c_library.comm_size );

    // Read input parameters
    c_parm.initialize(argc, argv);
    c_master_conf.input_parameters(c_parm);
    if (c_master_conf.method == 0 || c_master_conf.method == 1 ||  c_master_conf.method == 3 )
    //if (c_master_conf.method == 0 || c_master_conf.method == 1 ||  c_master_conf.method == 3  || c_master_conf.method == 5)
        c_library.input_parameters_input(c_parm);
    c_filter.input_parameters(c_parm);   //dbfilter code
    c_orient.input_parameters(c_parm);
    c_bmp_score.input_parameters(c_parm);
    c_master_score.input_parameters(c_parm);
    c_simplex.input_parameters(c_parm, c_master_conf.flexible_ligand, c_master_conf.genetic_algorithm, c_master_conf.denovo_design, c_master_score);

    // check parms compatablity
    if (c_master_conf.method == 3){ // if covalent some minimizer parameters are not compatable. 
         if (c_simplex.use_min_rigid_anchor){
             cout << "min anchor must be turned off for covalent" << endl;
             c_simplex.use_min_rigid_anchor = false;
             exit(0);
         }
         if (c_simplex.flex_min_max_iterations > 0.0){
             cout << "flex_min_max_iterations must be set to 0.0 for covalent" << endl;
             c_simplex.flex_min_max_iterations = 0.0;
             exit(0);
         }
    } 
     //Put a check to make sure denovo/ga and MPI are not specified together
     //Currently using de novo/ga under MPI causes weird things to happen
    if (USE_MPI && c_master_conf.method ==2){
        cout << "\nError: DOCK does not support the DOCK_DN  option under MPI"
             << endl;
        c_library.finalize_mpi();
        exit(0);
    }
    if (USE_MPI && c_master_conf.method ==3){
        cout << "\nError: DOCK does not support the GA option under MPI"
             << endl;
        c_library.finalize_mpi();
        exit(0);
    }
    if (USE_MPI && c_master_conf.method ==5){
        cout << "\nError: DOCK does not support the HDB option under MPI, for now"
             << endl;
        c_library.finalize_mpi();
        exit(0);
    }


    if (c_master_conf.flexible_ligand   || c_simplex.minimize_ligand 
        || c_filter.use_database_filter || c_orient.orient_ligand 
        || c_bmp_score.bump_filter      || c_library.calc_rmsd)
        c_master_score.read_vdw = true;

    // we should also read in the vdw.defn if we are calculating xlogp in the 
    // database filter i.e. if c_filter.use_database_filter is true. For now,
    // the xlogp code has been removed, so this is no longer needed
    
    //added this if statement and line so that we have the information from these parm files available for delta_max, the else portion was what was there originally
    if ( (c_master_conf.method==4 && c_master_conf.c_ga_recomb.use_limit_max_change) ){
        bool use_volume = true;
        bool use_ph4 = true;
        bool use_chemical_matching = true;
        bool use_vdw = true;
        c_master_score.use_ph4 = use_ph4;
        c_master_score.use_volume = use_volume; 
        c_orient.use_chemical_matching = use_chemical_matching;
        c_master_score.read_vdw = use_vdw;
        c_typer.input_parameters(c_parm, use_vdw, use_chemical_matching, use_ph4, use_volume);
        cout << "MODIFYING C_TYPER" << endl;
    }
    else{
        c_typer.input_parameters(c_parm, c_master_score.read_vdw,
                             c_orient.use_chemical_matching, c_master_score.use_ph4, c_master_score.use_volume); 
    }
    //if (c_master_conf.method == 0 || c_master_conf.method == 1 ||  c_master_conf.method == 3 )
    //if (c_master_conf.method == 0 || c_master_conf.method == 1 ||  c_master_conf.method == 3  || c_master_conf.method == 5)
    if (c_master_conf.method == 0 || c_master_conf.method == 1 ||  c_master_conf.method == 3  || c_master_conf.method == 5)
        c_library.input_parameters_output(c_parm, c_master_score, USE_MPI);
    mol.input_parameters(c_parm);

    if (!USE_MPI)
        c_parm.write_params();

    // Exit the program if parameterization fails //
    if (!c_parm.parameter_input_successful()) {
        if (USE_MPI)
            c_library.finalize_mpi();
        if (outfile.is_open()) {
            // return non error output to standard output.
            cout.rdbuf(original_cout_buffer);
            outfile.close();
            fclose(stdout);  // C stdio.
        }
        return 1;               // non-zero return indicates error
    }
    // /////////////////////////////////////////////////////////////////////////
    // Initialization routines
    c_master_conf.initialize();
    c_typer.initialize(c_master_score.read_vdw, c_master_score.read_gb_parm,
                       c_orient.use_chemical_matching, c_master_score.use_ph4, c_master_score.use_volume);
    //if (c_master_conf.method == 0 || c_master_conf.method == 1 ||  c_master_conf.method == 3  || c_master_conf.method == 5)
    if (c_master_conf.method == 0 || c_master_conf.method == 1 ||  c_master_conf.method == 3 )
        c_library.initialize(argc, argv, USE_MPI);
    c_orient.initialize(argc, argv);
    c_bmp_score.initialize();
    c_master_score.initialize_all(c_typer, argc, argv);  //removed reference to flex_min_add_internal: not used any more
    c_simplex.initialize();

#ifdef BUILD_DOCK_WITH_MPI
    if ((USE_MPI) && (c_library.rank > 0)) {
        cout << "DOCK is currently running on ";
        char klient[ MPI_MAX_PROCESSOR_NAME ];
        int length;
        MPI_Get_processor_name(klient, &length);
        for ( int i = 0; i < length; ++i )
            cout << klient[i];
        cout << endl;
    }
#endif

    //fstream main_dock_loop_anchors;
    //main_dock_loop_anchors.open ("unmin_anchors.mol2", fstream::out|fstream::app);

    // /////////////////////////////////////////////////////////////////////////
    // Main loop
    //if ligand has been docked
    //      write out and read in next ligand
    //if entire docking is complete
    //      perform final analysis and write out
    //else
    //      read in ligand

    // If you are doing de novo growth, enter this function
    if (c_master_conf.method == 2) {
        //if (c_master_conf.c_dn_build.simple_build_flag)
        //    c_master_conf.c_dn_build.simple_build(c_master_score, c_simplex, c_typer);
        //else
            c_master_conf.c_dn_build.build_molecules(c_master_score, c_simplex, c_typer, c_orient);

    } else if (c_master_conf.method == 3) { // this is covalent
    while (c_library.get_mol(mol,false, USE_MPI, c_master_score.amber, c_typer, c_master_score, c_simplex)) { 
        // If MPI is used this is done on the compute nodes.
        // filtering must be done here because it needs all prep for docking
        // before the mols can be eliminated.
 
        // keep track of time for individual molecules
        double          mol_start_time = wall_clock_seconds();
        //int             mol_start_ntime = wall_clock_nseconds();

        //seed random number generator
        c_simplex.initialize();

        //parse ligand atoms into child lists
        mol.prepare_molecule();
        mol.flag_write_solvation = c_library.write_solv_mol2;

        //label ligand atoms with proper vdw, bond, and chem
        //types
        c_typer.prepare_molecule(mol, c_master_score.read_vdw,
                                 c_orient.use_chemical_matching, c_master_score.use_ph4, c_master_score.use_volume);

        //parse ligand into rigid and flexible portions
        c_master_conf.prepare_molecule(mol);
 
        // Writing fragment libraries for denovo is done inside c_master_conf.prepare_molecule, 
        // so write some frags then continue to the next molecule
        if (c_master_conf.c_ag_conf.write_fragment_libraries){ continue; }

        // set this bool to true of every new ligand 
        // this is for rigid sampling
        // this ensures that the internal energy is only 
        // initialized once per ligand.
        c_master_conf.initialize_once = true;

        // Database Filter code
        //filter molecules by descriptors
        if (c_filter.use_database_filter)
        {
            //calculate descriptors for the ligand
            // c_filter.calc_descriptors(mol);
            // descriptors are now computed & printed in amber_typer.cpp

            //print out the descriptors
            // cout << c_filter.get_descriptors(mol);

            if (c_filter.fails_filter(mol))  {
                 //ligand failed the filter
                 //move to the next ligand to be docked
                 cout << "\n" "-----------------------------------" "\n";
                 cout << "Molecule: " << mol.title << "\n\n";
                 cout << mol.current_data << endl;
                 continue;
            }
        }

        // we want to run it once to store the orignal
        c_master_conf.next_anchor(mol);
        c_orient.match_ligand_covalent(mol,1.8,1.8,PI,false,false,true); // set the original molecule to the start
        c_master_conf.c_cg_conf.current_anchor = 0;  // reset to start at one.  
        //transform the ligand to covalent bond
        //c_orient.match_ligand_covalent(mol,c_master_conf.c_cg_conf.bondlength,c_master_conf.c_cg_conf.bondlength2,c_master_conf.c_cg_conf.angleval*(PI/180.0));
        //c_master_conf.next_anchor(mol);
        //cout << c_master_conf.c_cg_conf.bondlength << " " << c_master_conf.c_cg_conf.blstep << " " << c_master_conf.c_cg_conf.blstop << endl;

        float bl_val,bl_val2,aval;
        //float bl_val  = c_master_conf.c_cg_conf.bondlength;
        bl_val  = c_master_conf.c_cg_conf.bondlength;
        while ( bl_val <= c_master_conf.c_cg_conf.blstop) {
            //cout << "for debuging ... bl_val  = " << bl_val << endl;
            //float bl_val2 = c_master_conf.c_cg_conf.bondlength2;
            bl_val2 = c_master_conf.c_cg_conf.bondlength2;
            while ( bl_val2 <= c_master_conf.c_cg_conf.blstop2) {
                //cout << "                 bl_val2 = " << bl_val2 << endl;
                //float aval    = c_master_conf.c_cg_conf.angleval;
                aval    = c_master_conf.c_cg_conf.angleval;
                while ( aval <= c_master_conf.c_cg_conf.avstop) {
                    // keep track of time for individual placements
                    double          mol_start_time = wall_clock_seconds();

                    c_simplex.initialize();

                    //cout << "for debuging ... bl_val  = " << bl_val << endl;
                    //cout << "                 bl_val2 = " << bl_val2 << endl;
                    //cout << "                 aval    = " << aval << endl;

                    c_master_conf.next_anchor(mol);
                    // ajust = true, orient = true, frist = false
                    c_orient.match_ligand_covalent(mol,bl_val,bl_val2,aval*(PI/180.0),true,true,false);

                    float angle = 0.0;
                    //float inc_angle = (10.0 * (PI/180.0));
                    //float inc_angle = (2.0*PI/c_master_conf.c_cg_conf.num_sample_angles); 
                    float inc_angle = (float(c_master_conf.c_cg_conf.dihideral_step) * (PI/180.0)); 
                    float max_angle = 2.0*PI;
                    //while (angle < 2*PI){ // 360 degrees == 2*PI radians
                    //for (angle =0.0; angle < max_angle; angle=angle+inc_angle){ // 360 degrees == 2*PI radians
                    for (angle =0.0; angle - max_angle < -0.001; angle=angle+inc_angle){ // 360 degrees == 2*PI radians
                        //cout << angle << endl;
                        //cout << max_angle * 180.0 / PI << " " << angle * 180.0 / PI << endl;
                        //cout << max_angle << " " << angle << endl;
                        //cout << (350.0*(PI/180.0)) << endl;
                        // rotate about the covalent bond
                        bool tmpflag=true; // more_orients
                        //if (fabs(angle - (350.0*(PI/180.0))) < 0.001) { 
                        if (fabs(angle - (max_angle-inc_angle)) < 0.001) { 
                        //if ((fabs( bl_val - (c_master_conf.c_cg_conf.blstop - c_master_conf.c_cg_conf.blstep))<0.001) and 
                        //   (fabs( bl_val2 - (c_master_conf.c_cg_conf.blstop2 - c_master_conf.c_cg_conf.blstep2))<0.001) and 
                        //   (fabs( aval - (c_master_conf.c_cg_conf.avstop - c_master_conf.c_cg_conf.avstep))<0.001) and  
                        //   (fabs(angle - (max_angle-inc_angle)) < 0.001) ){ 
                             tmpflag=false;
                             //cout << "last angle" << endl;
                        }
                        c_orient.new_next_orientation_covalent(mol,angle);
                            //perform bump check on anchor and filter if fails
                            if (c_bmp_score.check_anchor_bumps(mol, tmpflag)) {
                            //cout << "Entering check_anchor_bumps" << endl; 
                                //Write_Mol2(mol, main_dock_loop_anchors);
                                //anchor minimization for flexible docking is moved to grow_periphery()
                                //score orientation and write out.
                                //submit_orientation returns false only if there is a
                                //problem with the score.
                                //no ligand output is produced outside of this if.
                                if(c_library.submit_orientation(mol, c_master_score, c_orient.orient_ligand) ||
                                   ! tmpflag ){
                                // score the ligands. if outside of grid, returns false. there was an issue if the 
                                // last orient was false then the anchors were not submitted to the growth routine.
                                    //cout << "Entering submit_orientation" << endl; 
                                    //rank orientations but only keep number user cutoff
                                    if(c_master_conf.submit_anchor_orientation(mol, tmpflag)){ 
                                        //cout << "Entering submit_anchor_orientation" << endl; 
                                        //add mol to (orientation) anchor_positions array 
                                        //prune anchors, then perform growth, minimization, and pruning
                                        //until molecule is fully grown
                                        c_master_conf.grow_periphery(c_master_score, c_simplex, c_bmp_score);
                                        //for the fully grown conformations, if there are conformations
                                        //remaining
                                        while (c_master_conf.next_conformer(mol)) {
                                                //cout << "Entering while loop next_conformer" << endl; 
                    
                                                //minimize the final pose
                                                c_simplex.minimize_final_pose(mol, c_master_score, c_typer);
                    
                                                //calculate score and internal
                                                //c_master_score.compute_primary_score(mol); 
                                                //cout << mol.score_text_data << endl;
                                                // add info about angle and bond length to mol2 header
                                                    int FLOAT_WIDTH = 20;
                                                    int STRING_WIDTH = 17 + 19;
                                                    string DELIMITER    = "########## ";
                                               
                                                    ostringstream text;
                                                    //ostringstream blbl2adha ;
                                                    ostringstream blbl2a ;
                                                    //blbl2adha << bl_val << ';' << bl_val2 << ';' << aval << ';' << angle*180/PI ;
                                                    blbl2a << bl_val << ';' << bl_val2 << ';' << aval  ;
                                                    //text << DELIMITER << setw(STRING_WIDTH) << "bl;bl2;a;dha:"      << setw(FLOAT_WIDTH) << fixed << blbl2adha.str()   << endl;
                                                    text << DELIMITER << setw(STRING_WIDTH) << "bl;bl2;a:"      << setw(FLOAT_WIDTH) << fixed << blbl2a.str()   << endl;
                                                    mol.current_data = mol.current_data + text.str();
                                                    //cout << mol.current_data << endl;
                                               
                                                //
                    
                                                //add best scoring pose to list for
                                                //ranking and further analysis
                                                c_library.submit_scored_pose(mol, c_master_score, c_simplex);
                                        }
                                        //write out list of final conformations
                                        c_library.submit_conformations(c_master_score);
                                    }
                                }
                            }
                        //angle = angle + (10.0 * (PI/180.0)); // 10 degree incraments
                    }

                    //cout << "\nI AM HERE" << endl;
            
        //print error messages if orienting or growth has failed
        if (c_library.num_orients == 0) {
            double          mol_stop_time = wall_clock_seconds();
//          int             mol_stop_ntime = wall_clock_nseconds();

            cout << "\n" "-----------------------------------" "\n";
            cout << "Molecule: " << mol.title << "\n\n";
            cout << " Elapsed time:\t" << fixed << setprecision(3)
                 << mol_stop_time - mol_start_time << " seconds\n\n";
            //cout << " Elapsed time:\t" << mol_stop_ntime - mol_start_ntime
            //     << " nseconds\n\n";
            cout << " ERROR:  Could not find a valid orientation." << endl;
            cout << "         (For rigid docking the sought orientation is "
                    "the whole ligand;\n"
                    "         for flexible docking the sought orientation is "
                    "the anchor.)\n"
                    "         Confirm that all spheres are inside the grid box"
                    "         \nand that the grid box is big enough to contain "
                    "an orientation.\n";
        } else if (c_library.num_confs == 0) {
            double          mol_stop_time = wall_clock_seconds();
            // int             mol_stop_ntime = wall_clock_nseconds();
            cout << "\n" "-----------------------------------" "\n";
            cout << "Molecule: " << mol.title << "\n\n";
            cout << " Elapsed time:\t" << fixed << setprecision(3)
                 << mol_stop_time - mol_start_time << " seconds\n\n";
            //cout << " Elapsed time:\t" << mol_stop_ntime - mol_start_ntime
            //     << " nseconds\n\n";
            cout << " ERROR:  Could not complete growth." << endl;
            cout << "         Confirm that the grid box is large enough to "
                    "contain the ligand,\n"
                    "         and try increasing max_orientations.\n";
        }
        // End individual molecule timing
        if (c_library.num_orients > 0 && c_library.num_confs > 0) {
            double          mol_stop_time = wall_clock_seconds();
            //int             mol_stop_ntime = wall_clock_nseconds();
            cout << "\n" "-----------------------------------" "\n";
            cout << "Molecule: " << mol.title << "\n\n";
            cout << " Elapsed time for docking:\t" << fixed << setprecision(3)
                 << mol_stop_time - mol_start_time << " seconds\n\n";
            //cout << " Elapsed time for docking:\t" << mol_stop_ntime - mol_start_ntime
            //     << " nseconds\n\n";

        }
                        
                    aval    = aval + c_master_conf.c_cg_conf.avstep;
                }
                bl_val2 = bl_val2 + c_master_conf.c_cg_conf.blstep2;
            }
            bl_val  = bl_val + c_master_conf.c_cg_conf.blstep;

                    //cout << "for debuging ... bl_val  = " << bl_val << endl;
                    //cout << "                 bl_val2 = " << bl_val2 << endl;
                    //cout << "blstop2 = "  << c_master_conf.c_cg_conf.blstop2 << endl;
                    //cout << "                 aval    = " << aval << endl;
        }

    }
  
    // If you are doing genetic algorithm, enter this function
    } else if (c_master_conf.method == 4) {
          c_master_conf.c_ga_recomb.max_breeding(c_master_score, c_simplex, c_typer, c_orient);
    
    // If you are doing hierarchical database (HDB) searching enter this function
    } else if (c_master_conf.method == 5) {
          // check if gist is used in combination with HDB search
          // full displacement gist is incompatable.
          if (c_master_score.c_gist.use_primary_score){
             if (c_master_score.c_gist.gist_type == "displace"){
                 cout << " Error... full displacement gist is incompatible with HDB search.\n\tYou can use blurry_displacement or trilinear instead.\n\tExiting the program. " << endl;
                 exit(0);
             }  
          } 
          c_library.db2flag = true;
          c_library.initialize(argc, argv, USE_MPI);
          c_library.initialize_input(); //  TEB 2024/04/11, input from JB 
          // mimic what is done in DOCK3.7 for now this will only work with serial program.
          HDB_Mol db2_data;
          DOCKMol mol_ac; //# all_coords_rigid_seg;
          //igzstream  db2_stream;
          string filename = c_master_conf.c_hdb_conf.db2filename;
          //bool flag_skip_broken = c_master_conf.c_hdb_conf.skip_broken; // if true do not read in broken sets (branches) if false read-in them. 

          string filenamedb2;

          cout << "opening file: " << filename << endl;

          vector<string> list;

          //c_library.read_sdifile("./split_database_file", list);
          c_library.read_sdifile(filename, list);
 
          cout << list.size()<<endl;
          for (int l=0; l<list.size();l++){

               filenamedb2 = list[l];

               
               cout << "opening :" << filenamedb2 <<endl;
               igzstream  db2_stream;
               db2_stream.clear();
               db2_stream.open(filenamedb2.c_str());

               if (db2_stream.fail()) {
                  cout << "\n\nCould not open " << filenamedb2 <<
                          " for reading.  Program will terminate." << endl << endl;
                  exit(0);
               }


               //while (c_library.read_hierarchy_db2( c_master_conf.c_hdb_conf.db2filename, db2_data)){//
               while (c_library.read_hierarchy_db2( db2_stream, db2_data)){
                  double          mol_start_time = wall_clock_seconds();
                  //DOCKMol rigid; //# anchor;
              
                  //mol.clear_molecule();  // this is done in c_hdb_conf.create_mol when allocation happens
                  //mol_ac.clear_molecule();
                  //c_library.read_hierarchy_db2( c_master_conf.c_hdb_conf.db2filename, db2_data); 
                  //ifstream   db2_stream;
                  //c_orient.cliques.clear();
                   
                  c_simplex.initialize();
                  
                  c_master_conf.c_hdb_conf.all_poses.clear();
              
                  c_master_conf.c_hdb_conf.create_mol(mol,mol_ac,db2_data,0);
              
                  mol.prepare_molecule();
              
                  c_typer.prepare_molecule(mol, c_master_score.read_vdw,
                                         c_orient.use_chemical_matching, c_master_score.use_ph4, c_master_score.use_volume);
                  c_master_conf.prepare_molecule(mol);
              
                  c_typer.prepare_molecule(mol_ac, c_master_score.read_vdw,
                                         c_orient.use_chemical_matching, c_master_score.use_ph4, c_master_score.use_volume);
                  
                  //copy_molecule(rigid, mol);
                  //cout << "I AM HERE" << endl;
                  c_master_conf.c_hdb_conf.set_rigid_active(mol,db2_data);
                  c_orient.match_ligand(mol);
                  c_library.num_anchors = 1;
              
                  int num = 0;
                  while (c_orient.new_next_orientation(mol, false)){
                    c_master_conf.c_hdb_conf.search(c_master_score,c_orient,c_bmp_score,mol,mol_ac,db2_data,num); 
                    //c_master_conf.c_hdb_conf.set_rigid_active(mol,db2_data);
                    num++;
                  }
              
                  c_library.num_orients = num;
                  
                  //mol_ac.flag_write_solvation = true;
                  //mol.flag_write_solvation = true;
                  mol_ac.flag_write_solvation = c_library.write_solv_mol2;
                  mol.flag_write_solvation = c_library.write_solv_mol2;
                  // sort the poses by energy
                  int number_min = 0;
                  if (c_master_conf.c_hdb_conf.num_per_search < c_master_conf.c_hdb_conf.all_poses.size()) {
                     number_min= c_master_conf.c_hdb_conf.num_per_search;
                  } else{ // 
                     number_min = c_master_conf.c_hdb_conf.all_poses.size();
                  }

                  if (number_min == 0){
                      cout << "Warning... number_min == 0." << endl;
                      //break;
                  }
                  //cout << "Before sorting:  " << c_master_conf.c_hdb_conf.all_poses[0].first << endl;
                  //
                  //int *top_X = new int [number_min]; // list of index for top x poses. 
                  int top_X[number_min]; // list of index for top x poses. 
                  
                  //c_library.sort_top_X_mol(c_master_conf.c_hdb_conf.all_poses, number_min, top_X);
                  if (! c_library.sort_top_X_mol(c_master_conf.c_hdb_conf.all_poses, number_min, top_X)) { 
                      
                      cout << "Warning... no poses found with min score." << endl;
                      //break;
                  }
                  //exit(0);
                  //sort(c_master_conf.c_hdb_conf.all_poses.begin(), c_master_conf.c_hdb_conf.all_poses.end(), less_than_pair);
                  //cout << "Affer sorting:  " << c_master_conf.c_hdb_conf.all_poses[0].first << endl;
                  
                  //Read in db2 file. 
                  //search db2 file and score segments
                  //score viable poses and write out the top scoring poses
                  //
                  //c_simplex.initialize();

                  //for (int i = 0; i < number_min; i++) {  
                  for (int ii = 0; ii < number_min; ii++) {  
                  //for (int i = 0; i < c_master_conf.c_hdb_conf.all_poses.size(); i++) {  
                         int i = top_X[ii];
                         //cout << i << endl;
                         //mol.flag_write_solvation = true;
                         mol.flag_write_solvation = c_library.write_solv_mol2;
                         //cout << "debug" << i << " " << c_master_conf.c_hdb_conf.all_poses[i].first <<endl;
                         copy_molecule(mol, c_master_conf.c_hdb_conf.all_poses[i].second);
                         //copy_crds(mol, c_master_conf.c_hdb_conf.all_poses[i].second);
                         //c_master_score.compute_primary_score(mol); 
                         //cout << i << " : score before min = " << mol.current_score << endl;  
                         //c_typer.prepare_molecule(mol, c_master_score.read_vdw,
                         //                c_orient.use_chemical_matching, c_master_score.use_ph4, c_master_score.use_volume);
                         //c_master_conf.prepare_molecule(mol);
                         c_simplex.minimize_final_pose(mol, c_master_score, c_typer);
                         c_master_score.compute_primary_score(mol); 
                         //cout << i << " : score affter min = " << mol.current_score << endl;  
                         //mol.score_text_data = mol.score_text_data + mol.hdb_data;
                         mol.current_data = mol.current_data + mol.hdb_data;
                         //add best scoring pose to list for
                         //ranking and further analysis
                         //mol.flag_write_solvation = true;
                         mol.flag_write_solvation = c_library.write_solv_mol2;
                         c_library.submit_scored_pose(mol, c_master_score, c_simplex);
                         //write out list of final conformations
                  }

                  //delete top_X;
                  //delete[] top_X;
                  c_library.submit_conformations(c_master_score);
               
                  //cout << "debug all_poses.size() = " << c_master_conf.c_hdb_conf.all_poses.size() << endl;
                  if (c_library.num_orients == 0) {
                      double          mol_stop_time = wall_clock_seconds();
                      cout << "\n" "-----------------------------------" "\n";
                      cout << "Molecule: " << mol.title << "\n\n";
                      cout << " Elapsed time:\t" << mol_stop_time - mol_start_time
                           << " seconds\n\n";
                      cout << " ERROR:  Could not find a valid orientation." << endl;
                      cout << "         (For rigid docking the sought orientation is "
                              "the whole ligand;\n"
                              "         for flexible docking the sought orientation is "
                              "the anchor.)\n"
                              "         Confirm that all spheres are inside the grid box"
                              "         \nand that the grid box is big enough to contain "
                              "an orientation.\n";
                  } else if (c_library.num_confs == 0) {
                      double          mol_stop_time = wall_clock_seconds();
                      cout << "\n" "-----------------------------------" "\n";
                      cout << "Molecule: " << mol.title << "\n\n";
                      cout << " Elapsed time:\t" << mol_stop_time - mol_start_time
                           << " seconds\n\n";
                      cout << " ERROR:  Could not complete growth." << endl;
                      cout << "         Confirm that the grid box is large enough to "
                              "contain the ligand,\n"
                              "         and try increasing max_orientations.\n";
                  }
                  // End individual molecule timing
                  if (c_library.num_orients > 0 && c_library.num_confs > 0) {
                      double          mol_stop_time = wall_clock_seconds();
                      cout << "\n" "-----------------------------------" "\n";
                      cout << "Molecule: " << mol.title << "\n\n";
                      cout << " Elapsed time for docking:\t"
                           << mol_stop_time - mol_start_time << " seconds\n\n";
                  }
                  //c_library.write_scored_poses(USE_MPI, c_master_score);
                  //c_library.sort_write(USE_MPI, c_master_score, c_simplex);
                  //c_library.sort_write(USE_FILT,USE_MPI, c_master_score, c_simplex);
                  c_library.sort_write(false,USE_MPI, c_master_score, c_simplex);
               }
               db2_stream.clear();
               db2_stream.close();
          } // list (of files)
/**/

    // Else if you are doing flexible, rigid, or fixed anchor docking, enter here
    } else {

    while (c_library.get_mol(mol,c_filter.use_database_filter, USE_MPI, c_master_score.amber, c_typer, c_master_score, c_simplex)) { 
        // If MPI is used this is done on the compute nodes.
        // filtering must be done here because it needs all prep for docking
        // before the mols can be eliminated.
 
        // keep track of time for individual molecules
#ifndef TIME_PRECISION
        double          mol_start_time = wall_clock_seconds();
#endif
#ifdef TIME_PRECISION
        struct timespec mol_start_time;
        timespec_get(&mol_start_time, TIME_UTC);
#endif
        //seed random number generator
        c_simplex.initialize();

        //parse ligand atoms into child lists
        mol.prepare_molecule();
        mol.flag_write_solvation = c_library.write_solv_mol2;

        //label ligand atoms with proper vdw, bond, and chem
        //types
        c_typer.prepare_molecule(mol, c_master_score.read_vdw,
                                 c_orient.use_chemical_matching, c_master_score.use_ph4, c_master_score.use_volume);

        c_master_conf.initialize_once = true;  // this is need for initializing the internal energy fuction for rigid docking and minimization
        //parse ligand into rigid and flexible portions
        c_master_conf.prepare_molecule(mol);
 
        // Writing fragment libraries for denovo is done inside c_master_conf.prepare_molecule, 
        // so write some frags then continue to the next molecule
        if (c_master_conf.c_ag_conf.write_fragment_libraries){ continue; }


        // Database Filter code
        //filter molecules by descriptors
        if (c_filter.use_database_filter)
        {
            //calculate RDKIT-related descriptors for the ligand
            c_filter.calc_descriptors(mol);
            //other descriptors are now computed & printed in amber_typer.cpp
            //print out the descriptors
            //cout << c_filter.get_descriptors(mol);
            #ifdef BUILD_DOCK_WITH_RDKIT
            if (c_filter.fails_filter(mol))  {
                 //ligand failed the filter 
                 mol.fails_filt = true;
                 //move to the next ligand to be docked
                 //cout << "\n" "-----------------------------------" "\n";
                 //cout << "Molecule: " << mol.title << "\n\n";
                 cout << mol.current_data << endl;
                 //continue;
            }else{
                mol.fails_filt = false;
            } 
            #else
            if (c_filter.fails_filter(mol))  {
                 //ligand failed the filter
                 //move to the next ligand to be docked
                 //cout << "\n" "-----------------------------------" "\n";
                 //cout << "Molecule: " << mol.title << "\n\n";
                 cout << mol.current_data << endl;
                 continue;
            }
            #endif
        }

        // sudipto & trent Dec 09, 2008
        // possible problem: Anchors are being generated inside this loop
        // Solution: Generate anchors outside this loop

        //while there is still another anchor fragment to be docked
        // when no anchor & grow is done, this while loop executes only once 
        while (c_master_conf.next_anchor(mol)) {
            trace.note( "Entering while loop next_anchor" );

            c_library.num_anchors++;

            //generate entire list of atom center-sphere center matches
            c_orient.match_ligand(mol);
            
            //transform the ligand to match
            // called only once for singlepoint or minimization only
            while (c_orient.new_next_orientation(mol)) {
                trace.note( "Entering while loop new_next_orientation" );

                //perform bump check on anchor and filter if fails
                trace.boolean("Perform bump_check::", c_bmp_score.check_anchor_bumps(mol, c_orient.more_orientations()));
                if (c_bmp_score.check_anchor_bumps(mol, c_orient.more_orientations())) {
                    trace.note( "Entering check_anchor_bumps if stmt" );

                    //Write_Mol2(mol, main_dock_loop_anchors);

                    //anchor minimization for flexible docking is moved to grow_periphery()

                    //score orientation and write out.
                    //submit_orientation returns false only if there is a
                    //problem with the score.
                    //no ligand output is produced outside of this if.

                    //this ensures that non-bonded pairlist is properly cleared
                    //before scoring orients from a new ligand during flex

                    //if (c_master_conf.method == 1 || c_master_conf.method == 2 || c_master_conf.method == 3 ||  c_master_conf.method == 4){
                    if (c_master_conf.method == 1 ){
                    //if (c_master_conf.method != 0){ // For Rigid, method == 0
                          c_master_score.primary_score->nb_int.clear();
                          trace.note("Clearing the Flex/non-bonded pairlist");
                    }
                    //this ensures that non-bonded pairlist is properly cleared
                    // and rebuilt before scoring orients from a new ligand 
                    //during rigid sampling
                    else {
                          c_master_score.primary_score->nb_int.clear();
                          trace.note("Clearing the Rigid non-bonded pairlist");
                          c_master_conf.initialize_once = true; // reset for internal energy for rigid docking. 
                          c_master_conf.grow_periphery(c_master_score, c_simplex, c_bmp_score);
                    }


                    if(c_library.submit_orientation(mol, c_master_score, c_orient.orient_ligand) ||
                       ! c_orient.more_orientations() ){
                    // score the ligands. if outside of grid, returns false. there was an issue if the 
                    // last orient was false then the anchors were not submitted to the growth routine.
                        trace.note( "Entering submit_orientation if stmt" );

                        //rank orientations but only keep number user cutoff
                        if(c_master_conf.submit_anchor_orientation(mol, c_orient.more_orientations())){ 
                            trace.note( "Entering submit_anchor_orientation if stmt" );
                            //add mol to (orientation) anchor_positions array 
                            //prune anchors, then perform growth, minimization, and pruning
                            //until molecule is fully grown
                            c_master_conf.grow_periphery(c_master_score, c_simplex, c_bmp_score);

                            //for the fully grown conformations, if there are conformations
                            //remaining
                            while (c_master_conf.next_conformer(mol)) {
                                trace.note( "Entering while loop next_conformer" );

                                //minimize the final pose
                                c_simplex.minimize_final_pose(mol, c_master_score, c_typer);

                                //calculate score and internal
                                // c_master_score.compute_primary_score(mol); 

                                //add best scoring pose to list for
                                //ranking and further analysis
                                c_library.submit_scored_pose(mol, c_master_score, c_simplex);
                            }

                        }
                    }
                }
            }
        }
        //write out list of final conformations
        c_library.submit_conformations(c_master_score);
        
#ifdef TIME_PRECISION
        //print error messages if orienting or growth has failed
            if (c_library.num_orients == 0) {
                //double          mol_stop_time = wall_clock_seconds();
                struct timespec mol_stop_time;
                timespec_get(&mol_stop_time, TIME_UTC);
                double          time_per_mol = calculate_simulation_time(mol_start_time, mol_stop_time);
                cout << "\n" "-----------------------------------" "\n";
                cout << "Molecule: " << mol.title << "\n\n";
                cout << " Elapsed time:\t";
                cout << fixed << setprecision(3) << time_per_mol << " seconds\n\n";
                cout << " ERROR:  Could not find a valid orientation." << endl;
                cout << "         (For rigid docking the sought orientation is "
                        "the whole ligand;\n"
                        "         for flexible docking the sought orientation is "
                        "the anchor.)\n"
                        "         Confirm that all spheres are inside the grid box"
                        "         \nand that the grid box is big enough to contain "
                        "an orientation.\n";
            } else if (c_library.num_confs == 0) {
            //double          mol_stop_time = wall_clock_seconds();
                struct timespec mol_stop_time;
                timespec_get(&mol_stop_time, TIME_UTC);
                double          time_per_mol = calculate_simulation_time(mol_start_time, mol_stop_time);
                cout << "\n" "-----------------------------------" "\n";
                cout << "Molecule: " << mol.title << "\n\n";
                cout << " Elapsed time:\t";
                cout << fixed << setprecision(3) << time_per_mol << " seconds\n\n";
                cout << " ERROR:  Could not complete growth." << endl;
                cout << "         Confirm that the grid box is large enough to "
                        "contain the ligand,\n"
                        "         and try increasing max_orientations.\n";
            }
            // End individual molecule timing
            if (!c_filter.use_database_filter){
                if (c_library.num_orients > 0 && c_library.num_confs > 0) {
                //double          mol_stop_time = wall_clock_seconds();
                    struct timespec mol_stop_time;
                    timespec_get(&mol_stop_time, TIME_UTC);
                    double          time_per_mol = calculate_simulation_time(mol_start_time, mol_stop_time);
                    cout << "\n" "-----------------------------------" "\n";
                    cout << "Molecule: " << mol.title << "\n\n";
                    cout << " Elapsed time for docking:\t";
                    cout << fixed << setprecision(3) << time_per_mol << " seconds\n\n";
                }
            } else {
                if (c_library.num_orients > 0 && c_library.num_confs > 0) {
                    struct timespec mol_stop_time;
                    timespec_get(&mol_stop_time, TIME_UTC);
                    double          time_per_mol = calculate_simulation_time(mol_start_time, mol_stop_time);
                    
                    cout << " Elapsed time for docking:\t";
                    cout << fixed << setprecision(3) << time_per_mol << " seconds\n\n"; 
                }
            }
             
#endif
#ifndef TIME_PRECISION
            //print error messages if orienting or growth has failed
            if (c_library.num_orients == 0) {
                double          mol_stop_time = wall_clock_seconds();

                cout << "\n" "-----------------------------------" "\n";
                cout << "Molecule: " << mol.title << "\n\n";
                cout << " Elapsed time:\t" << fixed << setprecision(3)
                     << mol_stop_time - mol_start_time << " seconds\n\n";
                cout << " ERROR:  Could not find a valid orientation." << endl;
                cout << "         (For rigid docking the sought orientation is "
                        "the whole ligand;\n"
                        "         for flexible docking the sought orientation is "
                        "the anchor.)\n"
                        "         Confirm that all spheres are inside the grid box"
                        "         \nand that the grid box is big enough to contain "
                        "an orientation.\n";
            } else if (c_library.num_confs == 0) {
                double          mol_stop_time = wall_clock_seconds();
                cout << "\n" "-----------------------------------" "\n";
                cout << "Molecule: " << mol.title << "\n\n";
                cout << " Elapsed time:\t" << fixed << setprecision(3)
                     << mol_stop_time - mol_start_time << " seconds\n\n";
                cout << " ERROR:  Could not complete growth." << endl;
                cout << "         Confirm that the grid box is large enough to "
                        "contain the ligand,\n"
                        "         and try increasing max_orientations.\n";
            }
            // End individual molecule timing
            if (c_library.num_orients > 0 && c_library.num_confs > 0) {
                double          mol_stop_time = wall_clock_seconds();
                cout << "\n" "----------------------------------" "\n";
                cout << "Molecule: " << mol.title << "\n\n";
                cout << " Elapsed time for docking:\t" << fixed << setprecision(3)
                     << mol_stop_time - mol_start_time << " seconds\n\n";

            }
#endif

    }
    }
    //main_dock_loop_anchors.close();
    // /////////////////////////////////////////////////////////////////////////


    // Write fragment libraries to file here
    if (c_master_conf.c_ag_conf.write_fragment_libraries)
        //c_master_conf.c_ag_conf.write_unique_fragments();
        c_master_conf.c_ag_conf.write_unique_fragments();
    // Rescore library using secondary scoring
    if (c_master_score.use_secondary_score) {
        c_library.secondary_rescore_poses(c_master_score, c_simplex);
        c_library.submit_secondary_pose();
    }
    
    int end_process_mols = 0;
    if ((!USE_MPI) || (c_library.rank == 0)){
        if (c_master_conf.method == 2){
            end_process_mols = c_master_conf.c_dn_build.molecule_counter;
            cout << "\n\n" << c_master_conf.c_dn_build.molecule_counter
                 << " Attachments Processed" << endl;
        } else if (! (c_master_conf.method == 3)){ 
            end_process_mols = c_library.total_mols - c_library.initial_skip;
            cout << "\n\n" << c_library.total_mols - c_library.initial_skip
                 << " Molecules Processed" << endl;
        }
    } else {
        if (c_master_conf.method == 2){
            end_process_mols = c_master_conf.c_dn_build.molecule_counter;
            cout << "\n\n" << c_master_conf.c_dn_build.molecule_counter
                 << " Attachments Processed" << endl;
        } else if (! (c_master_conf.method == 3)){
            end_process_mols = c_library.completed;
            cout << "\n\n" << c_library.completed
                 << " Molecules Processed" << endl;
        }
    }    

    // End timing
#ifndef TIME_PRECISION
    double          stop_time = wall_clock_seconds();
    double          total_time = stop_time - start_time;
    cout << "Total elapsed time:\t";
    cout << fixed << setprecision(3) << total_time << " seconds\n";
#endif
#ifdef TIME_PRECISION
    struct timespec stop_time;
    timespec_get(&stop_time, TIME_UTC);
    double          total_time = calculate_simulation_time(start_time, stop_time); 
    cout << "Total elapsed time:\t";
    cout << fixed << setprecision(3) << total_time << " seconds\n";
#endif
  
    float mol_per_sec = end_process_mols / total_time;

    if (c_master_conf.method == 2) {
        cout << "Number of attachments per second:\t"
             << fixed << setprecision(3) << mol_per_sec << "\n";
    } else if (! (c_master_conf.method == 3)) {
        cout << "Number of molecules per second:\t"
             << fixed << setprecision(3) << mol_per_sec << "\n";
    }    
    
#ifndef __APPLE__

    long long virtmem = getVirtValue();
    std::cout << "Virtual memory used for this process: " << virtmem
              << " kilobytes\n";

    long long physmem = getPhysValue();
    std::cout << "Physical memory used for this process: " << physmem
              << " kilobytes" << std::endl;
#endif


    if (USE_MPI)
        c_library.finalize_mpi();

    if (outfile.is_open()) {
        // return non error output to standard output.
        cout.rdbuf(original_cout_buffer);
        outfile.close();
        fclose(stdout);  // C stdio.
    }

    return 0;                   // zero return indicates success

}

/************************************************/
void
dock_new_handler()
{

    // Define the default behavior for the failure of new.

    cout << "Error: memory exhausted!" << endl
        << "  If this occurs during grid reading then a likely cause\n"
        << "  is a grid that is too large.  Some machines have a very\n"
        << "  restrictive policy on allocating available resources:\n"
        << "  increase the datasize, stacksize, and memoryuse\n"
        << "  using the limit, ulimit, or unlimit commands;\n"
        << "  for many Linuxes this command sequence will work:\n"
        << "  limit; unlimit; limit\n"
        << "  Otherwise read the man pages or apply trial and error\n"
        << "  to find the apt use of these commands.\n"
        << endl;

    cerr << "Error: memory exhausted!" << endl
        << "  If this occurs during grid reading then a likely cause\n"
        << "  is a grid that is too large.  Some machines have a very\n"
        << "  restrictive policy on allocating available resources:\n"
        << "  increase the datasize, stacksize, and memoryuse\n"
        << "  using the limit, ulimit, or unlimit commands;\n"
        << "  for many Linuxes this command sequence will work:\n"
        << "  limit; unlimit; limit\n"
        << "  Otherwise read the man pages or apply trial and error\n"
        << "  to find the apt use of these commands.\n"
        << endl;

    // Throw so that recovery or special error notification can be performed.
    throw           bad_alloc();
}

/************************************************/
void
print_header( bool USE_MPI, int processes )
{
    cout << "\n\n\n--------------------------------------\n"
         << DOCK_VERSION
         << "\n\nReleased " << DOCK_RELEASE_DATE
         << "\nCopyright UCSF"
         << "\n--------------------------------------\n";
    if (USE_MPI) {
        cout << "Parallel dock running "
             << processes
             << " MPI processes"
             << "\n--------------------------------------\n";
    }
    cout << endl;
}

/************************************************/
//make sure c++11 is available - full functionality not guaranteed
#ifdef TIME_PRECISION

double
calculate_simulation_time(timespec t_start, timespec t_end)
{
/***********************************************************************
timespec returns time as a number of seconds, timespec.tv_sec,       
and a number of nanoseconds, timespec.tv_nsec. The correct moment    
in time equals to             
timespec.tv_sec + static_cast<double>(timespec.tv_nsec) / 1000000000.0
***********************************************************************/
    double sec_diff = t_end.tv_sec - t_start.tv_sec;
    long nsec_diff = t_end.tv_nsec - t_start.tv_nsec;
    if (nsec_diff < 0){
        sec_diff = sec_diff - 1.0;
        nsec_diff = 1000000000 + nsec_diff;
    }
    double sec_fraction = static_cast<double>(nsec_diff) / 1000000000.0;
    double simulation_time = sec_diff + sec_fraction;
    return simulation_time;
}
#endif
//Just in case they don't have c++11
////#ifndef TIME_PRECISION
double
wall_clock_seconds()
{
    time_t          t;
    if (static_cast < time_t > (-1) == time(&t)) {
        cout << "Error from time function!  Elapsed time is erroneous." << endl;
    }
    return static_cast < double >(t);
}
/****************************************************/
/*int
wall_clock_nseconds()
{
    struct timespec ts;
    timespec_get(&ts, TIME_UTC);
    return static_cast < int >(ts.tv_nsec);
}*/
////#endif
