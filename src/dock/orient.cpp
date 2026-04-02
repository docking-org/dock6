#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include "orient.h"

using namespace std;


/*** matfit stuff ***/
#define SMALL  1.0e-20
#define SMALSN 1.0e-10
#define ABS(x)   (((x)<0)   ? (-(x)) : (x))

/************************************************************************/
static void
minimized_fit(double umat[3][3], double rm[3][3])
{
    // Sudipto: could someone add comments to explain what each varible
    // here describes, and more details about the algorithm used?

    double          rot[3][3],
                    turmat[3][3],
                    c[3][3],
                    coup[3],
                    dir[3],
                    step[3],
                    v[3],
                    rtsum,
                    rtsump,
                    rsum,
                    stp,
                    stcoup,
                    ud,
                    tr,
                    ta,
                    cs,
                    sn,
                    ac,
                    delta,
                    deltap,
                    gfac,
                    cle,
                    clep;
    int             i,j,k,l,m,  //loop counters
                    jmax,
                    ncyc,
                    nsteep,
                    nrem;

    /*
     * Rotate repeatedly to reduce couple about initial direction to zero.
     * Clear the rotation matrix
     */
    for (l = 0; l < 3; l++) {
        for (m = 0; m < 3; m++)
            rot[l][m] = 0.0;
        rot[l][l] = 1.0;
    }

    /*
     * Copy vmat[][] (sp) into umat[][] (dp) 
     */
    jmax = 30;
    rtsum = umat[0][0] + umat[1][1] + umat[2][2];
    delta = 0.0;

    for (ncyc = 0; ncyc < jmax; ncyc++) {
        /*
         * Modified CG. For first and every NSTEEP cycles, set previous step as 
         * zero and do an SD step 
         */
        nsteep = 3;
        nrem = ncyc - nsteep * (int) (ncyc / nsteep);

        if (!nrem) {
            for (i = 0; i < 3; i++)
                step[i] = 0.0;
            clep = 1.0;
        }

        /*
         * Couple 
         */
        coup[0] = umat[1][2] - umat[2][1];
        coup[1] = umat[2][0] - umat[0][2];
        coup[2] = umat[0][1] - umat[1][0];
        cle = sqrt(coup[0] * coup[0] + coup[1] * coup[1] + coup[2] * coup[2]);

        /*
         * Gradient vector is now -coup 
         */
        gfac = (cle / clep) * (cle / clep);

        /*
         * Value of rtsum from previous step 
         */
        rtsump = rtsum;
        deltap = delta;
        clep = cle;
        if (cle < SMALL){
            //cout <<ncyc<< ": cle < SMALL:" << cle << "<" << SMALL << endl;
            break;
        }

        /*
         * Step vector conjugate to previous 
         */
        stp = 0.0;
        for (i = 0; i < 3; i++) {
            step[i] = coup[i] + step[i] * gfac;
            stp += (step[i] * step[i]);
        }
        stp = 1.0 / sqrt(stp);

        /*
         * Normalised step 
         */
        for (i = 0; i < 3; i++)
            dir[i] = stp * step[i];

        /*
         * Couple resolved along step direction 
         */
        stcoup = coup[0] * dir[0] + coup[1] * dir[1] + coup[2] * dir[2];

        /*
         * Component of UMAT along direction 
         */
        ud = 0.0;
        for (l = 0; l < 3; l++)
            for (m = 0; m < 3; m++)
                ud += umat[l][m] * dir[l] * dir[m];


        tr = umat[0][0] + umat[1][1] + umat[2][2] - ud;
        ta = sqrt(tr * tr + stcoup * stcoup);
        cs = tr / ta;
        sn = stcoup / ta;

        /*
         * If cs<0 then posiiton is unstable, so don't stop 
         */
        if ((cs > 0.0) && (ABS(sn) < SMALSN)){
            //cout <<ncyc <<": (cs > 0.0) && (ABS(sn) < SMALSN):"<< cs << ">"<< 0.0 << " && " << ABS(sn) << "<" << SMALSN << endl;
            break;
        }

        /*
         * Turn matrix for correcting rotation:
         * 
         * Symmetric part 
         */
        ac = 1.0 - cs;
        for (l = 0; l < 3; l++) {
            v[l] = ac * dir[l];
            for (m = 0; m < 3; m++)
                turmat[l][m] = v[l] * dir[m];
            turmat[l][l] += cs;
            v[l] = dir[l] * sn;
        }

        /*
         * Asymmetric part 
         */
        turmat[0][1] -= v[2];
        turmat[1][2] -= v[0];
        turmat[2][0] -= v[1];
        turmat[1][0] += v[2];
        turmat[2][1] += v[0];
        turmat[0][2] += v[1];

        /*
         * Update total rotation matrix 
         */
        for (l = 0; l < 3; l++) {
            for (m = 0; m < 3; m++) {
                c[l][m] = 0.0;
                for (k = 0; k < 3; k++)
                    c[l][m] += turmat[l][k] * rot[k][m];
            }
        }

        for (l = 0; l < 3; l++)
            for (m = 0; m < 3; m++)
                rot[l][m] = c[l][m];

        /*
         * Update umat tensor 
         */
        for (l = 0; l < 3; l++)
            for (m = 0; m < 3; m++) {
                c[l][m] = 0.0;
                for (k = 0; k < 3; k++)
                    c[l][m] += turmat[l][k] * umat[k][m];
            }

        for (l = 0; l < 3; l++)
            for (m = 0; m < 3; m++)
                umat[l][m] = c[l][m];

        rtsum = umat[0][0] + umat[1][1] + umat[2][2];
        delta = rtsum - rtsump;

        /*
         * If no improvement in this cycle then stop 
         */
        if (ABS(delta) < SMALL){
            //cout <<ncyc<< ": ABS(delta) < SMALL: " << ABS(delta) << "<" << SMALL << endl;
            break;
        }

        /*
         * Next cycle 
         */
    }

    rsum = rtsum;

    /*
     * Copy rotation matrix for output 
     */

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            rm[i][j] = rot[i][j];       // can be transposed
}
//
/*************************************************************************/
// x1 is the ligand atoms, x2 is the sphere atoms.
// rm should rotate x1 onto x2
int
compute_rot_matrix(XYZVec & x1, XYZVec & x2, double rm[3][3], int n)
{
    int             i,
                    j;
    double          umat[3][3];


    if (n < 2) {
        return (0);
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            umat[i][j] = 0.0;
    }

    for (j = 0; j < n; j++) {
        umat[0][0] += x1[j].x * x2[j].x;
        umat[1][0] += x1[j].y * x2[j].x;
        umat[2][0] += x1[j].z * x2[j].x;

        umat[0][1] += x1[j].x * x2[j].y;
        umat[1][1] += x1[j].y * x2[j].y;
        umat[2][1] += x1[j].z * x2[j].y;

        umat[0][2] += x1[j].x * x2[j].z;
        umat[1][2] += x1[j].y * x2[j].z;
        umat[2][2] += x1[j].z * x2[j].z;
    }

    /*
    // calculate the distance/residual between spheres and orients // for debuging TEB 2021
    double tempresid = 0; 
    for (j = 0; j < n; j++) {
         double tempdist2 = 0;
         tempdist2 =+ (x1[j].x - x2[j].x) * (x1[j].x - x2[j].x);
         tempdist2 =+ (x1[j].y - x2[j].y) * (x1[j].y - x2[j].y);
         tempdist2 =+ (x1[j].z - x2[j].z) * (x1[j].z - x2[j].z);
         tempresid =+ sqrt(tempdist2);
    }
    cout << "calculate the distance/residual between spheres and orients (before minimized_fit): " << tempresid << endl;
    */

    minimized_fit(umat, rm);

    /*
    // calculate the distance/residual between spheres and orients // for debuging TEB 2021
    XYZCRD P; 
    tempresid = 0; 
    for (j = 0; j < n; j++) {
         double tempdist2 = 0;
         P.x = rm[0][0] * x1[j].x + rm[0][1] * x1[j].y + rm[0][2] * x1[j].z;
         P.y = rm[1][0] * x1[j].x + rm[1][1] * x1[j].y + rm[1][2] * x1[j].z;
         P.z = rm[2][0] * x1[j].x + rm[2][1] * x1[j].y + rm[2][2] * x1[j].z;
         tempdist2 =+ (P.x - x2[j].x) * (P.x - x2[j].x);
         tempdist2 =+ (P.y - x2[j].y) * (P.y - x2[j].y);
         tempdist2 =+ (P.z - x2[j].z) * (P.z - x2[j].z);
         tempresid =+ sqrt(tempdist2);
    }
    cout << "calculate the distance/residual between spheres and orients (after minimized_fit): " << tempresid << endl;
    */

    return (1);
}


/************************************************/
/************************************************/
/************************************************/
Orient::Orient()
{

    sph_dist_mat = NULL;
    lig_dist_mat = NULL;
    residual_mat = NULL;

}
/************************************************/
Orient::~Orient()
{

    delete[]sph_dist_mat;
    delete[]lig_dist_mat;
    delete[]residual_mat;

}
/************************************************/
void
Orient::input_parameters(Parameter_Reader & parm)
{
    cout << "\nOrient Ligand Parameters\n";
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    use_chemical_matching = false;
    orient_ligand = parm.query_param("orient_ligand", "yes", "yes no") == "yes";

    if (orient_ligand) {
        automated_matching = 
            parm.query_param("automated_matching", "yes", "yes no") == "yes";

        if (automated_matching) {
            tolerance = 0.25;
            //dist_min = 2.0;
            dist_min  = 0.0;
            min_nodes = 3;
            max_nodes = 10;
        } else {
            automated_matching = 
                 parm.query_param("automated_matching_iteration", "no", "yes no") == "yes";
            tolerance =
                atof(parm.query_param("distance_tolerance", "0.25").c_str());
            if (tolerance <= 0.0) {
                cout <<
                    "ERROR: Parameter must be a float greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }
            dist_min =
                atof(parm.query_param("distance_minimum", "0.5").c_str());
            //if (dist_min <= 0.0) {
            if (dist_min < 0.0) {
                cout <<
                    "ERROR: Parameter must be a float greater than or equal to zero.  Program will terminate."
                    << endl;
                exit(0);
            }
            min_nodes = atoi(parm.query_param("nodes_minimum", "3").c_str());
            if (min_nodes <= 0) {
                cout <<
                    "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }
            max_nodes = atoi(parm.query_param("nodes_maximum", "10").c_str());
            if (max_nodes <= 0) {
                cout <<
                    "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }

        }

        Active_Site_Spheres :: set_sphere_file_name( parm );

        max_orients = atoi(parm.query_param("max_orientations", "1000").c_str());
        if (max_orients <= 0) {
            cout <<
                "ERROR: Parameter must be a integer greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }

        critical_points =
            parm.query_param("critical_points", "no", "yes no") == "yes";

        use_chemical_matching =
            parm.query_param("chemical_matching", "no", "yes no") == "yes";
        if (use_chemical_matching) {
            chem_match_tbl_fname =
                parm.query_param("chem_match_tbl", "chem_match.tbl");
        }

        use_ligand_spheres =
            parm.query_param("use_ligand_spheres", "no", "yes no") == "yes";
        if (use_ligand_spheres) {
            lig_sphere_filename =
                parm.query_param("ligand_sphere_file", "ligand.sph");
        }

        verbose = 0 != parm.verbosity_level();   // -v is for verbose flag
        //TODO: Print more info about orienting progress on verbose
    }
}

/************************************************/
void
Orient::initialize(int argc, char **argv)
{

    if (orient_ligand) {
        cout << "Initializing Orienting Routines...\n";

        orig_tolerance = tolerance;
        prepare_receptor();

        // check validity of manual matching
        if (max_nodes < min_nodes) {
            cout << endl <<
                "ERROR:  Invalid range of nodes.  Program will terminate." <<
                endl;
            exit(0);
        }
    }
}

/************************************************/
void
Orient::prepare_receptor()
{
    get_spheres();
    calculate_sphere_distance_matrix();
    if (use_chemical_matching)
        read_chem_match_tbl();
}

/************************************************/
void
Orient::clean_up()
{
    current_clique = 0; //TEB 2012-02-12
    delete[]lig_dist_mat;
    delete[]residual_mat;
    lig_dist_mat = NULL;
    residual_mat = NULL;
}

/************************************************/
void
Orient::get_spheres()
{
    bool            found_cluster;
    CRITICAL_CLUSTER tmp_cluster;

    spheres = Active_Site_Spheres :: get_instance();
    num_spheres = spheres.size();

    if (verbose) cout << "Read in " << num_spheres 
                      << " spheres for orienting." << endl;

    // process critical points
    if (critical_points) {
        receptor_critical_clusters.clear();

        // loop over spheres
        for (int i = 0; i < spheres.size(); i++) {

            // if sphere belongs to a critical cluster
            if (spheres[i].critical_cluster > 0) {

                found_cluster = false;

                for (int j = 0; j < receptor_critical_clusters.size(); j++) {
                    // if the cluster has been identified previously
                    if (spheres[i].critical_cluster ==
                        receptor_critical_clusters[j].index) {
                        // add the sphere to the existing cluster
                        receptor_critical_clusters[j].spheres.push_back(i);
                        found_cluster = true;
                        break;
                    }
                }

                // else create new cluster entry
                if (!found_cluster) {
                    tmp_cluster.index = spheres[i].critical_cluster;
                    tmp_cluster.spheres.clear();
                    tmp_cluster.spheres.push_back(i);
                    receptor_critical_clusters.push_back(tmp_cluster);
                }
            }
        }
    }
}

/************************************************/
void
Orient::calculate_sphere_distance_matrix()
{
    delete[]sph_dist_mat; 
    sph_dist_mat = new double[num_spheres * num_spheres];

    for (int i = 0; i < num_spheres; i++)
        for (int j = 0; j < num_spheres; j++) {
            sph_dist_mat[num_spheres * i + j] = spheres[i].distance(spheres[j]);
        }
}

/************************************************/
// Function will use active heavy atoms or dummy atoms as spheres
void
Orient::get_centers(DOCKMol & mol)
{
    centers.clear();
    num_centers = 0;

    if (!use_ligand_spheres) {

        Sphere tmp;
        for (int atom = 0; atom < mol.num_atoms; atom++) {
            // CDS-09/26/16: added dummy atom
            if (mol.amber_at_heavy_flag[atom] || mol.atom_types[atom] == "Du") {
                if (mol.atom_active_flags[atom]) {
                    //cout << "atom "<< atom << "is a center point" << endl;
                    tmp.crds.x = mol.x[atom];
                    tmp.crds.y = mol.y[atom];
                    tmp.crds.z = mol.z[atom];
                    tmp.radius = 0.0;
                    tmp.surface_point_i = 0;
                    tmp.surface_point_j = 0;
                    tmp.critical_cluster = 0;
                    if (use_chemical_matching){
                        cout << "atom = " << atom << "; color = " << mol.chem_types[atom] << endl;
                        tmp.color = mol.chem_types[atom];
                    }
                    centers.push_back(tmp);
                    num_centers++;
                }
            }
        }

    } else {
        //read in ligand spheres 
        num_centers = read_spheres( lig_sphere_filename, centers );
    }
    // print the number of anchor heavy atoms
    if (verbose) cout << "Orienting " << num_centers 
                      << " anchor heavy atom centers" << endl;
}

/************************************************/
// Function for ligand to ligand matching, which will prepare the reference ligand spheres
// similar to get_spheres().  Added by CDS - 09/24/16
// If Critical Points is true, it will make the dummy spheres critical only! 
void
Orient::get_lig_reference_spheres(DOCKMol & mol)
{
   orig_tolerance = tolerance;
   // Use the ligand reference function to identify the heavy atoms/spheres
   // Clear spheres
   spheres.clear();
   num_spheres = 0;

   if (!use_ligand_spheres) {

      Sphere tmp;
      int critical_cluster = 1; // place holder to update critical cluster
      for (int atom = 0; atom < mol.num_atoms; atom++) {
          // CDS-09/26/16: added dummy atom
          if (mol.amber_at_heavy_flag[atom] || mol.atom_types[atom] == "Du") {
              if (mol.atom_active_flags[atom]) {
                 tmp.crds.x = mol.x[atom];
                 tmp.crds.y = mol.y[atom];
                 tmp.crds.z = mol.z[atom];
                 tmp.radius = 0.0;
                 tmp.surface_point_i = 0;
                 tmp.surface_point_j = 0;
                 // If there is a dummy atom, make it a critical cluster
                 if ( mol.atom_types[atom] == "Du" ){
                    tmp.critical_cluster = critical_cluster;
                    critical_cluster++;
                 }
                 else{
                    tmp.critical_cluster = 0;
                 }
                 if (use_chemical_matching)
                    tmp.color = mol.chem_types[atom];
                 spheres.push_back(tmp);
                 num_spheres++;
              }
          }
      }

    } else {
        //read in ligand spheres 
        num_spheres = read_spheres( lig_sphere_filename, centers );
    }
    // print the number of anchor heavy atoms
    if (verbose) cout << "Orienting " << num_spheres 
                      << " anchor heavy atom centers" << endl;



   // The code below is from get_spheres()
   bool            found_cluster;
   CRITICAL_CLUSTER tmp_cluster;

   if (verbose) cout << "Read in " << num_spheres 
                     << " spheres for orienting." << endl;

   // process critical points
   if (critical_points) {
       receptor_critical_clusters.clear();

       // loop over spheres
       for (int i = 0; i < spheres.size(); i++) {

           // if sphere belongs to a critical cluster
           if (spheres[i].critical_cluster > 0) {

               found_cluster = false;

               for (int j = 0; j < receptor_critical_clusters.size(); j++) {
                   // if the cluster has been identified previously
                   if (spheres[i].critical_cluster ==
                       receptor_critical_clusters[j].index) {
                       // add the sphere to the existing cluster
                       receptor_critical_clusters[j].spheres.push_back(i);
                       found_cluster = true;
                       break;
                   }
               }

               // else create new cluster entry
               if (!found_cluster) {
                   tmp_cluster.index = spheres[i].critical_cluster;
                   tmp_cluster.spheres.clear();
                   tmp_cluster.spheres.push_back(i);
                   receptor_critical_clusters.push_back(tmp_cluster);
               }
           }
       }
   }
   calculate_sphere_distance_matrix();
   if (use_chemical_matching)
      read_chem_match_tbl();
}

/************************************************/
void
Orient::calculate_ligand_distance_matrix()
{

    lig_dist_mat = new double[num_centers * num_centers];

    for (int i = 0; i < num_centers; i++)
        for (int j = 0; j < num_centers; j++) {
            lig_dist_mat[num_centers * i + j] = centers[i].distance(centers[j]);
        }
}


/************************************************/
// This fuction is called in the main dock loop.
// This function will rotate/translate the anchor to the covalent bond. 
// 
// this code is adapted from denovo.
//// +++++++++++++++++++++++++++++++++++++++++
//// Given two fragments and connection point data, combine them into one and return it
//Fragment
//DN_Build::combine_fragments( Fragment & frag1, int dummy1, int heavy1,
//                             Fragment frag2, int dummy2, int heavy2 )
//{
void
//Orient::match_ligand_covalent(DOCKMol & mol)
//Orient::match_ligand_covalent(DOCKMol & mol, float bondlength, float bondlength2, float angle)
Orient::match_ligand_covalent(DOCKMol & mol, float bondlength, float bondlength2, float angle, bool flag_adjust_bonds, bool flag_orient, bool flag_first)
{
    SCOREMol        tmp_mol;
    float           mol_score;
    //int             insert_point;
    //int             i;

    int             itmp;

    // if (verbose) { 
    if (verbose and not flag_first) { 
         cout << "-----------------------------------" << endl 
              << "COVALENT ORIENTING INFO :" << endl << endl;
    }
    if (not flag_first){  //  if it is not the frist time, we want to copy the already oreinted molecule into the mol // this will help reproduceablity.
       copy_molecule( mol, original);
    } 
    //LEP - debug flag of orienting 
    //if (Parameter_Reader::verbosity_level() > 0) {
        //ofstream myfile;
        //myfile.open ("debug.mol2");
        //Write_Mol2(mol, myfile);
    //}

    //cout << "CG_Conformer_Search::submit_anchor_orientation" << endl;
    mol_score = mol.current_score;

    int dummy1 = mol.dummy1;
    int dummy2 = mol.dummy2;
    int atomtag = mol.atomtag;

    //SphereVec recsph; 
    //
    //string file = "rec.sph";

    //itmp = read_spheres( file, recsph );
    // rec residue is in the correct position - the objective is to translate / rotate lig
    // so that the bond from dummy2->dummy1 is overlapping the bond from sphere2->sphere1

/*
    cout << dummy1 << " " 
         << dummy2 << " " 
         << atomtag << endl;
    cout << "dummy atoms from the ligand\n"
         << "DUMMY 1: \n"
         << " " << mol.x[dummy1] 
         << " " << mol.y[dummy1] 
         << " " << mol.z[dummy1] 
         << " " << mol.atom_names[dummy1] << endl;

    cout << "DUMMY2: \n" 
         << " " << mol.x[dummy2] 
         << " " << mol.y[dummy2] 
         << " " << mol.z[dummy2] 
         << " " << mol.atom_names[dummy2] << endl;
    cout << "ATOM for dhideral sampling: \n" 
         << " " << mol.x[atomtag] 
         << " " << mol.y[atomtag] 
         << " " << mol.z[atomtag] 
         << " " << mol.atom_names[atomtag] << endl;
*/

    if ( flag_adjust_bonds ){

    // Step 1. Adjust the bond length of covalent bond.
    //
    //    D2--D1--A1(this is the frist atom of the ligand)
    //

    //    Step 1.1 Frist, we will adjust the bond lenth for D1--A1 
    if (verbose) { 
      cout << "input bond lengths: " << endl;
      cout << "bondlength D1--A1 = " << sqrt(pow(mol.x[atomtag] - mol.x[dummy1],2.0) + pow( mol.y[atomtag] - mol.y[dummy1],2.0) +  pow(mol.z[atomtag] - mol.z[dummy1],2.0)) << endl;
      cout << "bondlength D2--D1 = " << sqrt(pow(mol.x[dummy1] - mol.x[dummy2],2.0)  + pow( mol.y[dummy1] - mol.y[dummy2],2.0)  +  pow(mol.z[dummy1] - mol.z[dummy2],2.0)) << endl;
    }

    if (bondlength2 > 0.0 and bondlength > 0.0 ) {  // if the bond lenth is positive then ajust the length. if you just adjusted D1-A1 and not D2-D1 then D2-D2 might be wierd (distorted): so, if you adjust D1-A1 you must adjust D2-D1.    
        float new_rad = bondlength2;
        // Calculate the x-y-z components of the bond vector (bond_vec)
        DOCKVector bond_vec;
        bond_vec.x = mol.x[atomtag] - mol.x[dummy1];
        bond_vec.y = mol.y[atomtag] - mol.y[dummy1];
        bond_vec.z = mol.z[atomtag] - mol.z[dummy1];
       
        // Normalize the bond vector then multiply each component by new_rad so that it is the desired
        // length
        bond_vec = bond_vec.normalize_vector();
        bond_vec.x *= new_rad;
        bond_vec.y *= new_rad;
        bond_vec.z *= new_rad;
        // Change the coordinates of the dummy atom so that the bond length is correct
        mol.x[dummy1] = mol.x[atomtag] - bond_vec.x;
        mol.y[dummy1] = mol.y[atomtag] - bond_vec.y;
        mol.z[dummy1] = mol.z[atomtag] - bond_vec.z;
        //Write_Mol2(mol, myfile);
    } // else (if the bond is negative) do not adjust the bond. 
 
    //    Step 1.1 second we will adjust the bond lenth for D2--D1 

    // Calculate the desired bond length and remember as 'new_rad'
    //float new_rad = calc_cov_radius(frag1.mol.atom_types[heavy1]) +
    //                calc_cov_radius(frag2.mol.atom_types[heavy2]);
    //float new_rad = 1.4;
    if (bondlength > 0.0 ) {  // if the bond lenth is positive then ajust the length.  (D2-D1 can be adjusted regardless if D1-A1 is adjusted.)
        float new_rad = bondlength;
        
        // Calculate the x-y-z components of the bond vector (bond_vec)
        DOCKVector bond_vec;
        bond_vec.x = mol.x[dummy1] - mol.x[dummy2];
        bond_vec.y = mol.y[dummy1] - mol.y[dummy2];
        bond_vec.z = mol.z[dummy1] - mol.z[dummy2];
        
        // Normalize the bond vector then multiply each component by new_rad so that it is the desired
        // length
        bond_vec = bond_vec.normalize_vector();
        bond_vec.x *= new_rad;
        bond_vec.y *= new_rad;
        bond_vec.z *= new_rad;
        // Change the coordinates of the dummy atom so that the bond length is correct
        mol.x[dummy2] = mol.x[dummy1] - bond_vec.x;
        mol.y[dummy2] = mol.y[dummy1] - bond_vec.y;
        mol.z[dummy2] = mol.z[dummy1] - bond_vec.z;
        //Write_Mol2(mol, myfile);
    } // else (if the bond is negative) do not adjust the bond. 
    if (verbose) { 
      cout << "ajusted bond lengths: " << endl;
      cout << "bondlength D1--A1 = " << sqrt(pow(mol.x[atomtag] - mol.x[dummy1],2.0) + pow( mol.y[atomtag] - mol.y[dummy1],2.0) +  pow(mol.z[atomtag] - mol.z[dummy1],2.0)) << endl;
      cout << "bondlength D2--D1 = " << sqrt(pow(mol.x[dummy1] - mol.x[dummy2],2.0)  + pow( mol.y[dummy1] - mol.y[dummy2],2.0)  +  pow(mol.z[dummy1] - mol.z[dummy2],2.0)) << endl;
    }

    // Step 1.5 Adjust angle.  
    // if angle is negative then leave angle unchanged. 

    if (angle < 0.0){
        cout << " set angle is not used. because value is negative.  " << endl;
    } else{
        mol.set_angle(atomtag,dummy1,dummy2,angle);
    } 

/*
    cout << "dummy atoms from the ligand\n"
         << "DUMMY 1: \n"
         << " " << mol.x[dummy1] 
         << " " << mol.y[dummy1] 
         << " " << mol.z[dummy1] 
         << " " << mol.atom_names[dummy1] << endl;

    cout << "DUMMY2: \n" 
         << " " << mol.x[dummy2] 
         << " " << mol.y[dummy2] 
         << " " << mol.z[dummy2] 
         << " " << mol.atom_names[dummy2] << endl;

    cout << "ATOM for dhideral sampling: \n" 
         << " " << mol.x[atomtag] 
         << " " << mol.y[atomtag] 
         << " " << mol.z[atomtag] 
         << " " << mol.atom_names[atomtag] << endl;
*/

    } // if flag_adjust_bonds


    if ( flag_orient ){
    // Step 2. Translate dummy2 of frag2 to the origin

    // Figure out what translation is required to move the dummy atom to the origin
    DOCKVector trans1;
    trans1.x = -mol.x[dummy2];
    trans1.y = -mol.y[dummy2];
    trans1.z = -mol.z[dummy2];

    // Use the dockmol function to translate the fragment so the dummy atom is at the origin
    mol.translate_mol(trans1);
    //myfile << "####### I AM HERE" << endl;
    //Write_Mol2(mol, myfile);

    // Step 3. Calculate dot product to determine theta (theta = angle between vec1 and vec2)

    // vec1 = vector pointing from heavy1 to dummy1 in frag1
    DOCKVector vec1;
    vec1.x = spheres[0].crds.x - spheres[1].crds.x;
    vec1.y = spheres[0].crds.y - spheres[1].crds.y;
    vec1.z = spheres[0].crds.z - spheres[1].crds.z;

/*
    cout  << "\n\nSpheres representing the covalent residue" << endl;
    cout   <<  "Atom 1, place dummy 1 here:\n"; 
    cout   <<  spheres[0].crds.x 
    << " " <<  spheres[0].crds.y 
    << " " <<  spheres[0].crds.z << endl; 
    cout   <<  "Atom 2, place dummy 2 here:\n"; 
    cout   <<  spheres[1].crds.x 
    << " " <<  spheres[1].crds.y 
    << " " <<  spheres[1].crds.z << endl; 
    cout   <<  "Atom 3, for diherdral sampling\n"; 
    cout   <<  spheres[2].crds.x 
    << " " <<  spheres[2].crds.y 
    << " " <<  spheres[2].crds.z << endl; 
*/

    // vec2 = vector pointing from dummy2 to heavy2 in frag2 (dummy2 is at the origin)
    DOCKVector vec2;
    vec2.x = mol.x[dummy1];
    vec2.y = mol.y[dummy1];
    vec2.z = mol.z[dummy1];

    // Declare some variables
    float dot;          // dot product value of vec1 and vec2
    float vec1_magsq;   // vec1 magnitude-squared
    float vec2_magsq;   // vec2 magnitude-squared
    float cos_theta;    // cosine of theta
    float sin_theta;    // sine of theta

    // Compute the dot product using the function in utils.cpp
    dot = dot_prod(vec1, vec2);

    // Compute these magnitudes (squared)
    vec1_magsq = (vec1.x * vec1.x) + (vec1.y * vec1.y) + (vec1.z * vec1.z);
    vec2_magsq = (vec2.x * vec2.x) + (vec2.y * vec2.y) + (vec2.z * vec2.z);

    // Compute cosine and sine of theta (theta itself is not actually calculated)
    //if (vec1_magsq != 0.0 && vec2_magsq != 0.0 ){
    if (vec1_magsq > 0.0000001 && vec2_magsq > 0.0000001 ){
      cos_theta = dot / (sqrt (vec1_magsq * vec2_magsq));
    }
    else{ // 
      cos_theta = 1.0;
      //if dot >= 0.0{
      //   cos_theta = 1.0;
      //}else{
      //   cos_theta = -1.0;
      //}
    }

    // cos_theta should never be greater than one. 
    // if it is then it is a numeric issue. 
    if (cos_theta >= 1.0){
       cos_theta = 1.0;
       sin_theta = 0.0;
    }
    else if (cos_theta <= -1.0){ // likewise if it is less than negative one.  
       cos_theta = -1.0;
       sin_theta = 0.0;
    }
    else{
       sin_theta = sqrt (1 - (cos_theta * cos_theta));
    }

    //cout << "dot product info: " << dot << " " <<
    //        vec1_magsq << " " <<
    //        vec2_magsq << " " <<
    //        cos_theta << " " <<
    //        sin_theta <<endl;
    // dot product info: 3.60413 3.24744 4 1 -nan

    // Step 4. Rotate vec2 to be coincident with vec1

    // If cos_theta is -1, the vectors are parallel but in the opposite direction
    if (cos_theta == -1){

        // Declare the rotation matrix and rotate frag2
        double finalmat[3][3] = { { -1, 0, 0}, {0, -1, 0}, {0, 0, -1} };
        mol.rotate_mol(finalmat);
    }

    // If cos_theta is 1, vec1 and vec2 are already parallel - only translation is needed.
    // Otherwise, enter this loop and calculate out how to rotate frag2
    else if (cos_theta != 1) {

        // Calculate cross product of vec1 and vec2 to get U (function from utils.cpp)
        DOCKVector normalU = cross_prod(vec1, vec2);

        // Calculate cross product of vec2 and U to get ~W
        DOCKVector normalW = cross_prod(vec2, normalU);

        // Normalize the vectors
        vec2 = vec2.normalize_vector();
        normalU = normalU.normalize_vector();
        normalW = normalW.normalize_vector();


        // (1) Make coordinate rotation matrix, which rotates {e1, e2, e3} coordinate to
        // {normalW, vec2, normalU} coordinate
        float coorRot[3][3];
        coorRot[0][0] = normalW.x;  coorRot[0][1] = vec2.x;  coorRot[0][2] = normalU.x;
        coorRot[1][0] = normalW.y;  coorRot[1][1] = vec2.y;  coorRot[1][2] = normalU.y;
        coorRot[2][0] = normalW.z;  coorRot[2][1] = vec2.z;  coorRot[2][2] = normalU.z;



        // (3) Make inverse  matrix of coorRot matrix - since coorRot is an orthogonal matrix,
        // the inverse is its transpose, (coorRot)^T
        float invcoorRot[3][3];
        invcoorRot[0][0] = coorRot[0][0];  invcoorRot[0][1] = coorRot[1][0];
        invcoorRot[0][2] = coorRot[2][0];
       
        invcoorRot[1][0] = coorRot[0][1];  invcoorRot[1][1] = coorRot[1][1];
        invcoorRot[1][2] = coorRot[2][1];

        invcoorRot[2][0] = coorRot[0][2];  invcoorRot[2][1] = coorRot[1][2];
        invcoorRot[2][2] = coorRot[2][2];

      /************************ TEB, 2020/04/03
       * Why the angle is never negative? 
       *
       * It is because the rotation matrices [ norm_W ,  norm_v2 , norm_U] transposed  
       * includes norm_U which is determined by taking the cross product of the vectors 
       * v1 and v2.  if the angle between v1 and v2 is positive (counterclockwise) then 
       * the vector is pointed out of the plane, if the angle is negative then the vector
       * (norm_U) is pointing into the plane.  This vector (norm_U) is placed onto the 
       * positive z axis.  Thus, even, if the the cross-poduct is pointing into the plane
       * when it is transformed it is point out. 
       *
       ************************  
       */

      /*
        // check the sign of the angle, and sign.  TEB 2020/04/02
        float v1[3], v2[3];
        v1[0] = vec1.x; v2[0] = vec2.x; 
        v1[1] = vec1.y; v2[1] = vec2.y; 
        v1[2] = vec1.z; v2[2] = vec2.z; 
        
        float sign = check_neg_angle(v1,v2,invcoorRot);
        //float sign = check_neg_angle(v2,v1,invcoorRot);
        // sin(-theata) = -sin(theata)
        sin_theta = sign * sin_theta;
       */

        // (2) Make rotation matrix, which rotates vec2 theta angle on a plane of vec2 and normalW
        // to the direction of normalW
        float planeRot[3][3];
        planeRot[0][0] =  cos_theta;  planeRot[0][1] = sin_theta;  planeRot[0][2] = 0;
        planeRot[1][0] = -sin_theta;  planeRot[1][1] = cos_theta;  planeRot[1][2] = 0;
        planeRot[2][0] =          0;  planeRot[2][1] =         0;  planeRot[2][2] = 1;


        // (4) Multiply three matrices together:  [coorRot * planeRot * invcoorRot]
        float temp[3][3];
        double finalmat[3][3];

        // First multiply coorRot * planeRot, save as temp
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                temp[i][j] = 0.0;
                for (int k=0; k<3; k++){
                    temp[i][j] += coorRot[i][k]*planeRot[k][j];
                }
            }
        }

        // Then multiply temp * invcoorRot, save as finalmat
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++){
                finalmat[i][j] = 0.0;
                for (int k=0; k<3; k++){
                    finalmat[i][j] += temp[i][k]*invcoorRot[k][j];
                }
            }
        }

        // Rotate frag2 using finalmat[3][3]
        mol.rotate_mol(finalmat);
    }

    //Write_Mol2(mol, myfile);

    // Step 5. Translate frag2 to frag1

    // This is the translation vector to move dummy2 to heavy1 (dummy2 is at the origin)
    DOCKVector trans2;
    trans2.x = spheres[1].crds.x;
    trans2.y = spheres[1].crds.y;
    trans2.z = spheres[1].crds.z;

    // Use the dockmol function to translate frag2
    mol.translate_mol(trans2);
    //Write_Mol2(mol, myfile);
    //myfile.close();
    } // flag_orient
    if ( flag_first){ // copy the mol into original // we can do this without orienting it at all.  I thought that we only want to do this once after we orient it, but I think it is better to not orient it.   
       copy_molecule(original, mol);
    }
}

/************************************************/
void
//bool
Orient::new_next_orientation_covalent(DOCKMol & mol, float angle)
{

   /**************************************
 *     
 *                     ( Ligand )
 *                     o a4
 *                    /
 *                   /
 *       a2  o---(---o a3
 *          /    di
 *         /
 *     a1 o
 *   (Receptor)
 *
 *  a1 is sphere 3
 *  a2 is sphere 2 and dummy 2
 *  a3 is sphere 1 and dummy 1 (this is the attachement point
 *  a4 is is a ligand atom
 *  di is the dihideral to sample
 *
   ***************************************/ 

   
    int dummy1 = mol.dummy1;
    int dummy2 = mol.dummy2;
    int atomtag = mol.atomtag;
    //ofstream myfile;
    //myfile.open ("debug.mol2",ios::app);

    //float angle = 0.0;
//    while (angle < 2*PI){ // 360 degrees == 2*PI radians
      // myfile << "######## angle:  " << angle << endl;
      set_torsion(
      //spheres[2].crds.x, mol.x[dummy1],  mol.x[dummy2], mol.x[atomtag],
      //spheres[2].crds.y, mol.y[dummy1],  mol.y[dummy2], mol.y[atomtag],
      //spheres[2].crds.z, mol.z[dummy1],  mol.z[dummy2], mol.z[atomtag],  
      spheres[2].crds.x, mol.x[dummy2],  mol.x[dummy1], mol.x[atomtag],
      spheres[2].crds.y, mol.y[dummy2],  mol.y[dummy1], mol.y[atomtag],
      spheres[2].crds.z, mol.z[dummy2],  mol.z[dummy1], mol.z[atomtag],  
      angle, mol);
    //SCOREMol        tmp_mol;
      //anchor_positions.push_back(tmp_mol);
      //copy_molecule(anchor_positions[anchor_positions.size() - 1].
      //                    second, mol);
      // Write_Mol2(mol,myfile); 

//      angle = angle + (10.0 * (PI/180.0)); // 10 degree incraments
//    }
    //myfile.close();
//    return true;
}




/*************************************/
// This is copied from dockmol.cpp 
// This function seems to be using simple rotation matrices
// Why not use the fancy quternion stuff? :sudipto
// this function is is modifed from that which is in ?? 
// note that the angles is the new dihideral angle requested.
void
Orient::set_torsion(float x0, float x1, float x2, float x3, 
 float y0 , float y1, float y2, float y3, 
 float z0 , float z1, float z2, float z3, 
 float angle, DOCKMol & mol)
{
    //int             tor[4];
    //vector < int   >atoms;
    float           v1x,
                    v1y,
                    v1z,
                    v2x,
                    v2y,
                    v2z,
                    v3x,
                    v3y,
                    v3z;
    float           c1x,
                    c1y,
                    c1z,
                    c2x,
                    c2y,
                    c2z,
                    c3x,
                    c3y,
                    c3z;
    float           c1mag,
                    c2mag,
                    radang,
                    costheta,
                    m[9];
    float           nx,
                    ny,
                    nz,
                    mag,
                    rotang,
                    sn,
                    cs,
                    t,
                    tx,
                    ty,
                    tz;
    int             i,
                    j,
                    idx;


    // calculate the torsion angle
    v1x = x0 - x1;
    v2x = x1 - x2;
    v1y = y0 - y1;
    v2y = y1 - y2;
    v1z = z0 - z1;
    v2z = z1 - z2;
    v3x = x2 - x3;
    v3y = y2 - y3;
    v3z = z2 - z3;

    c1x = v1y * v2z - v1z * v2y;
    c2x = v2y * v3z - v2z * v3y;
    c1y = -v1x * v2z + v1z * v2x;
    c2y = -v2x * v3z + v2z * v3x;
    c1z = v1x * v2y - v1y * v2x;
    c2z = v2x * v3y - v2y * v3x;
    c3x = c1y * c2z - c1z * c2y;
    c3y = -c1x * c2z + c1z * c2x;
    c3z = c1x * c2y - c1y * c2x;

    c1mag = pow(c1x, 2) + pow(c1y, 2) + pow(c1z, 2);
    c2mag = pow(c2x, 2) + pow(c2y, 2) + pow(c2z, 2);

    if (c1mag * c2mag < 0.01)
        costheta = 1.0;         // avoid div by zero error
    else
        costheta = (c1x * c2x + c1y * c2y + c1z * c2z) / (sqrt(c1mag * c2mag));

    if (costheta < -0.999999)
        costheta = -0.999999f;
    if (costheta > 0.999999)
        costheta = 0.999999f;

    if ((v2x * c3x + v2y * c3y + v2z * c3z) > 0.0)
        radang = -acos(costheta);
    else
        radang = acos(costheta);

    // 
    // now we have the torsion angle (radang) - set up the rot matrix
    // 

    // find the difference between current and requested
    rotang = angle - radang;

    sn = sin(rotang);
    cs = cos(rotang);
    t = 1 - cs;

    // normalize the rotation vector
    mag = sqrt(pow(v2x, 2) + pow(v2y, 2) + pow(v2z, 2));
    nx = v2x / mag;
    ny = v2y / mag;
    nz = v2z / mag;

    // set up the rotation matrix
    m[0] = t * nx * nx + cs;
    m[1] = t * nx * ny + sn * nz;
    m[2] = t * nx * nz - sn * ny;
    m[3] = t * nx * ny - sn * nz;
    m[4] = t * ny * ny + cs;
    m[5] = t * ny * nz + sn * nx;
    m[6] = t * nx * nz + sn * ny;
    m[7] = t * ny * nz - sn * nx;
    m[8] = t * nz * nz + cs;

    // 
    // now the matrix is set - time to rotate the atoms
    // 
    tx = x1;
    ty = y1;
    tz = z1;

    //for (i = 0; i < mol.num_atoms; i++) {
        //j = atoms[i];
    for (j = 0; j < mol.num_atoms; j++) {

        // for(i=0;i<num_atoms;i++) {
        // j = i;

        // if(child_list[a2*num_atoms + i] == a3) { //////////////

        mol.x[j] -= tx;
        mol.y[j] -= ty;
        mol.z[j] -= tz;

        nx = mol.x[j] * m[0] + mol.y[j] * m[1] + mol.z[j] * m[2];
        ny = mol.x[j] * m[3] + mol.y[j] * m[4] + mol.z[j] * m[5];
        nz = mol.x[j] * m[6] + mol.y[j] * m[7] + mol.z[j] * m[8];

        mol.x[j] = nx;
        mol.y[j] = ny;
        mol.z[j] = nz;
        mol.x[j] += tx;
        mol.y[j] += ty;
        mol.z[j] += tz;

        // } /////////////////////////
    }

}



/************************************************/
// This fuction is called in the main dock loop.

void
Orient::match_ligand(DOCKMol & mol)
{
    int             s1,
                    s2,
                    c1,
                    c2;
    int             i,
                    j,
                    k;
    INTVec          cmt_sphere_idx,
                    cmt_center_idx;

    if (orient_ligand) {

        if (verbose) cout << "-----------------------------------" << endl 
             << "VERBOSE ORIENTING STATS :" << endl << endl;

        //cached_orient.clear_molecule();
        original.clear_molecule();
        centers.clear();
        clique_spheres.clear();
        clique_centers.clear();

        num_orients = 0;
        last_orient_flag = false;

        copy_molecule(original, mol);

        get_centers(original);

        // calc chem matching table if chem matching is used
        if (use_chemical_matching) {
            if (verbose) {
                cout << " using chemical matching " << endl;
                for (i = 0; i < spheres.size(); i++) {
                   cout << "sphere = " << i << "; color = " << spheres[i].color << endl;
                }
            } 
            // idx = center*num_spheres + sphere
            chem_match_align_tbl.clear();
            chem_match_align_tbl.resize((spheres.size() * centers.size()), 0);

            // assign chemical match type indices to the spheres
            cmt_sphere_idx.clear();
            cmt_sphere_idx.resize(spheres.size(), 0);

            // loop over spheres
            for (i = 0; i < spheres.size(); i++) {
                // loop over match types
                for (j = 0; j < chem_match_tbl_labels.size(); j++) {
                    // if types match, make the assignment
                    if (spheres[i].color == chem_match_tbl_labels[j]) {
                        cmt_sphere_idx[i] = j;
                    }
                }
            }

            // assign chemical match type indices to the centers
            cmt_center_idx.clear();
            // bug fix from revision 1.32 and 1.33; srb
            cmt_center_idx.resize( centers.size(), 0 );

            // loop over spheres
            for (i = 0; i < centers.size(); i++) {
                // loop over match types
                for (j = 0; j < chem_match_tbl_labels.size(); j++) {
                    // if types match, make the assignment
                    if (centers[i].color == chem_match_tbl_labels[j]) {
                        cmt_center_idx[i] = j;
                    }
                }
            }

            // generate table of legal sphere/center pairings
            for (i = 0; i < spheres.size(); i++) {
                for (j = 0; j < centers.size(); j++) {
                    chem_match_align_tbl[j * spheres.size() + i] =
                        chem_match_tbl_matrix[cmt_sphere_idx[i] *
                                              chem_match_tbl_labels.size() +
                                              cmt_center_idx[j]];
                }
            }

        }
        // end chemical matching code

        calculate_ligand_distance_matrix();

	// num_centers is the number of anchor heavy/dummy atoms
        num_nodes = num_spheres * num_centers;
        //cout << "Number of spheres for orienting: " << num_spheres << endl;
        // initialize matrix full of node-node residuals
        residual_mat = new double[num_nodes * num_nodes];

        k = 0;
        //cout << "match_ligand().num_nodes: "<< num_nodes << endl;
	// populate residual matrix
        for (i = 0; i < num_nodes; i++) {
            s1 = i % num_spheres;
            c1 = i / num_spheres;

            for (j = 0; j < num_nodes; j++) {
                s2 = j % num_spheres;
                c2 = j / num_spheres;

                //residual_mat[k] =
                //    fabs(sph_dist_mat[s1 * num_spheres + s2] -
                //         lig_dist_mat[c1 * num_centers + c2]);
                // if two spheres are too close together discard them. 
                if (sph_dist_mat[s1 * num_spheres + s2] < dist_min) {
                   residual_mat[k] = 10000.0;
                }
                else {
                   residual_mat[k] =
                      fabs(sph_dist_mat[s1 * num_spheres + s2] -
                           lig_dist_mat[c1 * num_centers + c2]);
                }
                k++;
            }
        }

        level = 0;
        am_iteration_num = 1;

        // perform clique detection
        id_all_cliques();

    } else
        last_orient_flag = false;
}

/************************************************/
void
Orient::calculate_translations()
{

    double sph_com_x = 0.0;
    double sph_com_y = 0.0;
    double sph_com_z = 0.0;
    double cen_com_x = 0.0;
    double cen_com_y = 0.0;
    double cen_com_z = 0.0;

    for (int i = 0; i < clique_size; i++) {
        cen_com_x += clique_centers[i].x;
        cen_com_y += clique_centers[i].y;
        cen_com_z += clique_centers[i].z;
        sph_com_x += clique_spheres[i].x;
        sph_com_y += clique_spheres[i].y;
        sph_com_z += clique_spheres[i].z;
    }

    sph_com_x = sph_com_x / clique_size;
    sph_com_y = sph_com_y / clique_size;
    sph_com_z = sph_com_z / clique_size;
    cen_com_x = cen_com_x / clique_size;
    cen_com_y = cen_com_y / clique_size;
    cen_com_z = cen_com_z / clique_size;

    spheres_com.x = sph_com_x;
    spheres_com.y = sph_com_y;
    spheres_com.z = sph_com_z;
    centers_com.x = cen_com_x;
    centers_com.y = cen_com_y;
    centers_com.z = cen_com_z;
}

/************************************************/
void
Orient::translate_clique_to_origin()
{

    for (int i = 0; i < clique_size; i++) {

        clique_centers[i].x = clique_centers[i].x - centers_com.x;
        clique_centers[i].y = clique_centers[i].y - centers_com.y;
        clique_centers[i].z = clique_centers[i].z - centers_com.z;

        clique_spheres[i].x = clique_spheres[i].x - spheres_com.x;
        clique_spheres[i].y = clique_spheres[i].y - spheres_com.y;
        clique_spheres[i].z = clique_spheres[i].z - spheres_com.z;
    }
}

/************************************************/
void
Orient::calculate_rotation()
{
    int n = clique_size;
    int i, j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            rotation_matrix[i][j] = 0.0;

    compute_rot_matrix(clique_centers, clique_spheres, rotation_matrix, n);
}

/************************************************/
bool
Orient::more_orientations()
{
    if (orient_ligand) {
        return ! last_orient_flag;
    } else
        return false;
}

/************************************************/
// called in match_ligand()

void
Orient::id_all_cliques()
{
    int             i,
                    j,
                    k,
                    size;
    int             next_cand,
                    index;
    // int notcount;
    CLIQUE          tmp_clique;
    double          *new_level_residuals;
    double           tmp_resid;

    int limit_cliques = 0;
    am_iteration_num = 1;
    cliques.clear();

    // init arrays
    size = num_centers * num_spheres;
    bool *new_candset = new bool[size * max_nodes];
    bool *new_notset = new bool[size * max_nodes];
    int *new_state = new int[max_nodes];
    new_level_residuals = new double[max_nodes];

    double loop_tolerance =  tolerance;

    if (verbose){
           cout <<"Sphere Center Matching Parameters:" << endl
                <<"   tolerance: "<<tolerance 
                <<"; dist_min: "<< dist_min
                <<"; min_nodes: "<<min_nodes
                <<"; max_nodes: "<<max_nodes
		<<"; max orients: "<<max_orients<< endl; //BTB`
    }
    // Loop over automated matching loop &&(am_iteration_num < 2)
    while ((cliques.size() <= max_orients)
           && ((automated_matching) || (am_iteration_num == 1))
           && (am_iteration_num < 10)) {

        loop_tolerance = am_iteration_num * orig_tolerance;

        for (i = 0; i < size * max_nodes; i++) {
            if (i < size)
                new_candset[i] = true;
            else
                new_candset[i] = false;

            new_notset[i] = false;

            if (i < max_nodes) {
                new_state[i] = -1;
                new_level_residuals[i] = 0.0;
            }
        }

        int new_level = 0;

        // main loop over levels
        while (new_level > -1) {
            // find next true in candset
            next_cand = -1;
            for (i = new_state[new_level] + 1; (i < size) && (next_cand == -1); i++) {
                index = new_level * size + i;

                if ((new_candset[index] == true) && (new_notset[index] == false)) {
                    next_cand = i;
		}
                // compute residuals
                if (next_cand != -1) {

                    // compute residuals
                    if (new_level > 0){
                        new_level_residuals[new_level] =
                            new_level_residuals[new_level - 1];
		    }
                    else if (new_level == 0){
                        new_level_residuals[new_level] = 0.0;
		    }
                    else { 
                        cout << "ERROR new_level is negative." << endl;
                    }
                    for (j = new_level - 1; j > -1; j--) {
                        k = next_cand * size + new_state[j];
                        // sum total resid method
                        new_level_residuals[new_level] += residual_mat[k];

                        // single max resid method
                        // if(residual_mat[k] > new_level_residuals[new_level])
                        // new_level_residuals[new_level] = residual_mat[k];
                    }

                    if (new_level_residuals[new_level] > loop_tolerance) {
                        next_cand = -1;
                        new_candset[index] = false;
                        new_notset[index] = true;
                    }

                }
                // end compute residuals
            }

            // if a candidate node is found at the current level
            if (next_cand != -1) {

                new_state[new_level] = next_cand;
                new_candset[new_level * size + new_state[new_level]] = false;

                new_level++;
                if (new_level < max_nodes)
                    new_state[new_level] = new_state[new_level - 1];
                if (new_level >= max_nodes) { 
                    // add state to cliques if tolerance is proper
                    if ((new_level_residuals[new_level - 1] >=
                         orig_tolerance * (am_iteration_num - 1))
                        && (new_level_residuals[new_level - 1] < loop_tolerance)) {
                        limit_cliques++;

                        tmp_clique.nodes.clear();
                        tmp_clique.residual =
                            new_level_residuals[new_level - 1];
                        for (i = 0; i < new_level; i++)
                            tmp_clique.nodes.push_back(new_state[i]);

                        if (check_clique_critical_points(tmp_clique))
                            if (check_clique_chemical_match(tmp_clique))
                                cliques.push_back(tmp_clique);

                    }
                    // end clique add code

                    // if you've reached the end of the tree and still have
                    // nodes to add
                    new_level--;

                    if (new_level >= 0)
                        new_notset[new_level * size + new_state[new_level]] =
                            true;

                } else {

                    // recompute candset and notset
                    for (i = 0; i < size; i++) {
                        index = new_state[new_level - 1] * size + i;
                        new_candset[new_level * size + i] =
                            ((new_candset[(new_level - 1) * size + i]));
                        new_notset[new_level * size + i] =
                            ((new_notset[(new_level - 1) * size + i]));

                        // compute residuals
                        if ((new_candset[new_level * size + i])
                            || (new_notset[new_level * size + i])) {

                            tmp_resid = new_level_residuals[new_level - 1];

                            for (j = new_level - 1; j > -1; j--) {
                                k = i * size + new_state[j];

                                // sum total resid method
                                tmp_resid += residual_mat[k];

                                // single max resid method
                                // if(residual_mat[k] > tmp_resid)
                                // tmp_resid = residual_mat[k];

                            }

                            if (tmp_resid > loop_tolerance) {
                                new_candset[new_level * size + i] = false;
                                new_notset[new_level * size + i] = false;
                            }
                        }
                        // End residual comp
                    }

                }

            } else {            // if no candidates are found
                //cout << "a";
                if (new_level >= min_nodes) {

                    // add state to cliques
                    if ((new_level_residuals[new_level - 1] >=
                         orig_tolerance * (am_iteration_num - 1))
                        && (new_level_residuals[new_level - 1] < loop_tolerance)) {

                        tmp_clique.nodes.clear();
                        tmp_clique.residual =
                            new_level_residuals[new_level - 1];
                        for (i = 0; i < new_level; i++)
                            tmp_clique.nodes.push_back(new_state[i]);

                        // perform critical point checking
                        if (check_clique_critical_points(tmp_clique))
                            if (check_clique_chemical_match(tmp_clique))
                                cliques.push_back(tmp_clique);

                    }
                    // end clique add code
                }

                new_state[new_level] = -1;
                new_level--;
            }

            // end if 100X the # of orients are found.  This is mostly for the
            // dense tree cases.
            if (limit_cliques >= 100 * max_orients) {
                cout << "Warning:  Match Search Truncated due to too many " <<
                    max_nodes << " cliques." << endl;
                break;
            }

        }

        am_iteration_num++;

    } // End automated matching loop

    // sort by residuals
    sort(cliques.begin(), cliques.end());

    if (verbose){
        cout.precision(4);
        cout << fixed;
        cout << "Num of cliques generated: " << cliques.size() << endl;
        if (cliques.size() > 0){
            cout << " Residual Info:"  << endl;
            cout << "   min residual:    " << cliques[0].residual << endl;
            cout << "   median residual: " << cliques[(int)(cliques.size() / 2)].residual << endl;
            cout << "   max residual:    " << cliques[cliques.size() - 1].residual << endl;
            double temp_resid_sum=0;
            double temp_resid2_sum=0;
            int max_node_size = 0;
            int min_node_size = 999999;
            double node_size_sum = 0;
            for (int i = 0; i < cliques.size(); i++){
                temp_resid_sum += cliques[i].residual;
                temp_resid2_sum += pow(cliques[i].residual,2);
                if (min_node_size > cliques[i].nodes.size())
                     min_node_size = cliques[i].nodes.size();
                if (max_node_size < cliques[i].nodes.size())
                     max_node_size = cliques[i].nodes.size();
                node_size_sum += cliques[i].nodes.size();
            }
            double mean = temp_resid_sum / cliques.size();
            double std  = sqrt(temp_resid2_sum/ cliques.size() - pow(mean,2));
            cout << "   mean residual:   " << mean << endl;
            cout << "   std residual:    "  << std  << endl;
            cout << " Node Sizes:"  << endl;
            cout << "   min nodes:    "  << min_node_size << endl;
            cout << "   max nodes:    "  << max_node_size << endl;
            cout << "   mean nodes:   " << node_size_sum / cliques.size() << endl;
            cout.unsetf ( ios_base::fixed  ); 
        }
             
    }

    current_clique = 0;

    // clean up arrays
    delete[]new_candset;
    delete[]new_notset;
    delete[]new_state;
    delete[]new_level_residuals;
    new_candset = NULL;
    new_notset = NULL;
    new_state = NULL;
    new_level_residuals = NULL;

}

/************************************************/
void
Orient::new_extract_coords_from_clique(CLIQUE & clique)
{
    // XYZCRD tmp;
    Sphere          tmp;
    int             index,
                    center,
                    sphere;
    int             i;

    clique_spheres.clear();
    clique_centers.clear();
    clique_size = clique.nodes.size();

    for (i = 0; i < clique_size; i++) {
        index = clique.nodes[i];
/*
        sphere = index / num_centers;
        center = index % num_centers;
*/
        // replace it with the proper sphere/center indexing
        sphere = index % num_spheres;
        center = index / num_spheres;
        // DTM - End removal of faulty code - 1/30/07

        tmp = spheres[sphere];
        clique_spheres.push_back(tmp.crds);

        tmp = centers[center];
        clique_centers.push_back(tmp.crds);
    }

}

/************************************************/
bool
Orient::check_clique_chemical_match(CLIQUE & clique)
{
    int             i;
    int             idx,
                    sphere_idx,
                    center_idx;

    if (use_chemical_matching) {

        // loop over clique nodes
        for (i = 0; i < clique.nodes.size(); i++) {

            // extract sphere & center indices
            idx = clique.nodes[i];
/*
            sphere_idx = idx / num_centers;
            center_idx = idx % num_centers;
*/
            // correct orienting: sudipto & DTM
            // replace it with the proper sphere/center indexing
            sphere_idx = idx % num_spheres;
            center_idx = idx / num_spheres;
            // DTM - End removal of faulty code - 1/30/07

            // check the chem_match_align_tbl - return false if match is
            // illegal
            idx = center_idx * num_spheres + sphere_idx;

            if (chem_match_align_tbl[idx] == 0)
                return false;

        }

        // else return true if match is legal
        return true;

    } else
        return true;

}

/************************************************/
bool
Orient::check_clique_critical_points(CLIQUE & clique)
{
    Sphere          tmp;
    int             index,
                    sphere;
    int             i,
                    j,
                    k;
    INTVec          tmp_clique_spheres;
    INTVec          hits;
    int             hit_sum;

    if (critical_points) {

        tmp_clique_spheres.clear();

        for (i = 0; i < clique.nodes.size(); i++) {
            index = clique.nodes[i];
            sphere = index % num_spheres;    // correct orienting: sudipto & DTM
            tmp_clique_spheres.push_back(sphere);
        }

        hits.clear();
        hits.resize(receptor_critical_clusters.size(), 0);

        // loop over the clusters
        for (i = 0; i < hits.size(); i++) {

            // loop over the cluster spheres
            for (j = 0; j < receptor_critical_clusters[i].spheres.size(); j++) {

                // loop over the clique spheres
                for (k = 0; k < tmp_clique_spheres.size(); k++) {

                    if (receptor_critical_clusters[i].spheres[j] ==
                        tmp_clique_spheres[k]) {
                        hits[i] = 1;
                        break;
                    }
                }

                if (hits[i] == 1)
                    break;
            }
        }

        hit_sum = 1;

        for (i = 0; i < hits.size(); i++) {
            hit_sum *= hits[i];
        }

        if (hit_sum == 0)
            return false;
        else
            return true;

    } else {
        return true;
    }
}


/************************************************/
// Called in main loop in dock.cpp.  Is a condition in while loop.
bool
Orient::new_next_orientation(DOCKMol & mol)
{
    bool all_atoms = true;
    return new_next_orientation(mol, all_atoms);
}

bool
//Orient::new_next_orientation(DOCKMol & mol)
Orient::new_next_orientation(DOCKMol & mol, bool all_atoms)
{
    //cout << "new_next_orientation" << endl;
    if (orient_ligand) {

        // in case no cliques could be found
        if (cliques.size() == 0) {
	    if (verbose) cout << "No orients found for current anchor" << endl;
            return false;
	}

        if (last_orient_flag) {
	    if (verbose) cout << "Current clique:" << current_clique << endl;
	    if (verbose) cout << "(last_orient_flag==true)Current clique:" << current_clique << endl;
            clean_up();
	    if (verbose) cout << "(last_orient_flag==true)Current clique:" << current_clique << endl;
            return false;
        }

        new_extract_coords_from_clique(cliques[current_clique]);
        calculate_translations();
        translate_clique_to_origin();
        calculate_rotation();

        copy_molecule(mol, original);
        //copy_crds(mol, original);
        //mol.translate_mol(-centers_com);
        //mol.rotate_mol(rotation_matrix);
        //mol.translate_mol(spheres_com);
        mol.translate_mol(-centers_com,all_atoms);
        mol.rotate_mol(rotation_matrix,all_atoms);
        mol.translate_mol(spheres_com,all_atoms);

        current_clique++;

        // DTM - change this line to ensure all cliques are examined as orients (stop skipping the last one) - 1/30/07
        if ((current_clique == max_orients) || (current_clique == cliques.size())){
	    if (verbose) cout << "Current clique:" << current_clique << endl;
            last_orient_flag = true;
        } else{
            last_orient_flag = false;
            
        }

        return true;

    } else {
        // code to return the input mol once, and then return false after that
        // (allow one pass through the loop)
        if (!last_orient_flag) {
            last_orient_flag = true;
            return true;
        } else {
            return false;
        }
    }
    return false;

}

bool
//Orient::new_next_orientation(DOCKMol & mol)
Orient::iso_new_next_orientation(Fragment & iso_frag, bool all_atoms)
{
    //cout << "new_next_orientation" << endl;
    if (orient_ligand) {

        // in case no cliques could be found
        if (cliques.size() == 0) {
	    if (verbose) cout << "No orients found for current anchor" << endl;
            return false;
	}

        if (last_orient_flag) {
	    if (verbose) cout << "Current clique:" << current_clique << endl;
	    if (verbose) cout << "(last_orient_flag==true)Current clique:" << current_clique << endl;
            clean_up();
	    if (verbose) cout << "(last_orient_flag==true)Current clique:" << current_clique << endl;
            return false;
        }

        new_extract_coords_from_clique(cliques[current_clique]);
        calculate_translations();
        translate_clique_to_origin();
        calculate_rotation();

        copy_molecule(iso_frag.mol, original);


        // grab rotataion matrx info for reference
        for(int i =0; i<3; i++){
            for(int j =0; j<3; j++){
                iso_frag.iso_ori_mat[i][j] = rotation_matrix[i][j];
            }
        }
        iso_frag.iso_centers_com = centers_com;
        iso_frag.iso_spheres_com = spheres_com;
 

        iso_frag.mol.translate_mol(-centers_com,all_atoms);
        iso_frag.mol.rotate_mol(rotation_matrix,all_atoms);
        iso_frag.mol.translate_mol(spheres_com,all_atoms);

        current_clique++;

        // DTM - change this line to ensure all cliques are examined as orients (stop skipping the last one) - 1/30/07
        if ((current_clique == max_orients) || (current_clique == cliques.size())){
	    if (verbose) cout << "Current clique:" << current_clique << endl;
            last_orient_flag = true;
        } else{
            last_orient_flag = false;
            
        }

        return true;

    } else {
        // code to return the input mol once, and then return false after that
        // (allow one pass through the loop)
        if (!last_orient_flag) {
            last_orient_flag = true;
            return true;
        } else {
            return false;
        }
    }
    return false;

}

/************************************************/
// Called in hdb search routine to orient each segement before scorring

bool
Orient::orientation_HDB(DOCKMol & mol, bool first)
{
    //cout << "orientation_HDB" << endl;
    if (orient_ligand) {
        //cout << "current_clique = " << current_clique << endl;
        if (first) { // only do this if it is the frist set
           new_extract_coords_from_clique(cliques[current_clique-1]);
           calculate_translations();
           translate_clique_to_origin();
           calculate_rotation();
        }
        //copy_molecule(mol, original);
        mol.translate_mol(-centers_com);
        mol.rotate_mol(rotation_matrix);
        mol.translate_mol(spheres_com);
    }
    return false;
}


/************************************************/
void
Orient::read_chem_match_tbl()
{
    FILE           *ifp;
    char            line[100],
                    tmp[100];
    string          tmp_string;
    int             i,
                    j,
                    z;

    if (use_chemical_matching) {

        // open match table file and loop over lines
        ifp = fopen(chem_match_tbl_fname.c_str(), "r");

        if (ifp == NULL) {
            cout << "\n\nCould not open " << chem_match_tbl_fname <<
                " for reading.  Program will terminate." << endl << endl;
            exit(0);
        }

        while (fgets(line, 100, ifp) != NULL) {

            // read in match table chemical labels
            if (!strncmp(line, "label", 5)) {
                sscanf(line, "%*s %s", tmp);
                tmp_string = tmp;
                chem_match_tbl_labels.push_back(tmp_string);
            }

            chem_match_tbl_matrix.
                resize((chem_match_tbl_labels.size() *
                        chem_match_tbl_labels.size()), 0);

            // read in match table matrix
            if (!strncmp(line, "table", 5)) {
                for (i = 0; i < chem_match_tbl_labels.size(); i++) {
                    for (j = 0; j <= i; j++) {
                        fscanf(ifp, "%d", &z);
                        chem_match_tbl_matrix[i * chem_match_tbl_labels.size() +
                                              j] = z;
                        chem_match_tbl_matrix[j * chem_match_tbl_labels.size() +
                                              i] = z;
                    }
                }
            }

        }

        fclose(ifp);
    }
    /*
     * for(i=0;i<chem_match_tbl_labels.size();i++) {
     * for(j=0;j<chem_match_tbl_labels.size();j++) { cout <<
     * chem_match_tbl_matrix[i*chem_match_tbl_labels.size() + j] << " "; } cout 
     * << endl; } 
     */
}
