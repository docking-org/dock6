#include "sasa.h"

////////////////////////////////////////////////////////////////////////////////
// / //////
// / Solvent Accessible Surface Area ////// 
// / //////
// / C++ implementation of icosasurf.f from amber8 ////// 
// / (C) Noel Carrascal, Sudipto Mukherjee, Robert C. Rizzo ////// 
// / SUNY Stony Brook //////
// / March 22, 2006 //////
// / //////
// / Modified by Yulin Huang and Trent Balius 2011-2012////// 
// / //////
////////////////////////////////////////////////////////////////////////////////


// / Noel Carrascal's implementation of icosasurf.f from amber8 in dock begins
// here
const int       ictri[20][3] =
    { {0, 4, 1}, {0, 9, 4}, {9, 5, 4}, {4, 5, 8}, {4, 8, 1},
{8, 10, 1}, {8, 3, 10}, {5, 3, 8}, {5, 2, 3}, {2, 7, 3},
{7, 10, 3}, {7, 6, 10}, {7, 11, 6}, {11, 0, 6}, {0, 1, 6},
{6, 1, 10}, {9, 0, 11}, {9, 11, 2}, {9, 2, 5}, {7, 2, 11}
};
const double    icosasym[12][3] =
    { {-0.525731112119133606, 0.0, 0.850650808352039932},
{0.525731112119133606, 0.0, 0.850650808352039932},
{-0.525731112119133606, 0.0, -0.850650808352039932},
{0.525731112119133606, 0.0, -0.850650808352039932},
{0.0, 0.850650808352039932, 0.525731112119133606},
{0.0, 0.850650808352039932, -0.525731112119133606},
{0.0, -0.850650808352039932, 0.525731112119133606},
{0.0, -0.850650808352039932, -0.525731112119133606},
{0.850650808352039932, 0.525731112119133606, 0.0},
{-0.850650808352039932, 0.525731112119133606, 0.0},
{0.850650808352039932, -0.525731112119133606, 0.0},
{-0.850650808352039932, -0.525731112119133606, 0.0}
};

// Max nof neighboring atoms including an icosahedra point
const int       iemax = 200;
int             ind = 0;

/*
 * bool icosa_point_exclusion(DOCKMol &,double[],double[],int,int
 * [],int*,int*,double,double,double); void icosa_init(int, int, double); void 
 * icosa_atom_environment(int *,int[],int *,int[]); double
 * icosa_sphere_approx(DOCKMol &,int, double[], double[],int *,int[]); void
 * gen_rot_mat(double, double [], double [][3]); void rot_point(double [],
 * double[][3]); double icosa_patch_approx(DOCKMol
 * &,int,double[],double[],double*,double*,double*,
 * bool*,bool*,bool*,int*,int*,int*,int*, int*,int*,int,int[],int); void
 * feed_decsasa(); 
 */

double
sasa::getSASA(ShortDOCKMol & molec, int recatomnum)
{
    // Wrapper function that calls getSASA(components)
    float sasa_tot;
    float sasa_hydrophilic;
    float sasa_hydrophobic;
    float sasa_other;
    float sasa_lig_tot;
    float sasa_lig_hydrophilic;
    float sasa_lig_hydrophobic;
    float sasa_lig_other;
    float sasa_rec_tot;
    float sasa_rec_hydrophilic;
    float sasa_rec_hydrophobic;
    float sasa_rec_other;
    //getSASA(molec,sasa_tot ,sasa_hydrophilic, sasa_hydrophobic, sasa_other);
    getSASA(molec, recatomnum, sasa_tot,sasa_hydrophilic, sasa_hydrophobic, sasa_other,
                               sasa_lig_tot, sasa_lig_hydrophilic, sasa_lig_hydrophobic, sasa_lig_other,
                               sasa_rec_tot ,sasa_rec_hydrophilic, sasa_rec_hydrophobic, sasa_rec_other);

/*
    cout <<"SASAtot     = " << sasa_tot         << "; " 
         <<"SASAphilic  = " << sasa_hydrophilic << "; " 
         <<"SASAphobic  = " << sasa_hydrophobic << "; " 
         <<"SASAother   = " << sasa_other     << endl;
    
    cout <<"SASA_ligtot     = " << sasa_lig_tot         << "; "
         <<"SASA_ligphilic  = " << sasa_lig_hydrophilic << "; "
         <<"SASA_ligphobic  = " << sasa_lig_hydrophobic << "; "
         <<"SASA_ligother   = " << sasa_lig_other     << endl;

    cout <<"SASA_rectot     = " << sasa_rec_tot         << "; "
         <<"SASA_recphilic  = " << sasa_rec_hydrophilic << "; "
         <<"SASA_recphobic  = " << sasa_rec_hydrophobic << "; "
         <<"SASA_recother   = " << sasa_rec_other     << endl;
    cout <<"###############################################" << endl;
*/

    //CHECK if the sasa compoenets add up
    //if (fabs(sasa_tot - (sasa_hydrophilic + sasa_hydrophobic + sasa_other)) > 0.01 )
    //    cout << "Warning: sasa_tot ["<< sasa_tot <<"]" << "!=" << "(sasa_hydrophilic + sasa_hydrophobic + sasa_other) [" 
    //         <<(sasa_hydrophilic + sasa_hydrophobic + sasa_other) << "]" <<" -- diff == " 
    //         << (sasa_tot - (sasa_hydrophilic + sasa_hydrophobic + sasa_other)) << endl; 
    return sasa_tot;
}

// Need a better atom typing scheme
// use the one in X log P code?? : sudipto
//TODO: Use amber typer to find out better types
// Note that H includes both polar and non-polar hydrogens

bool sasa::atom_is_hydrophobic(string atom_type){
    if (atom_type=="C.3"||atom_type=="C.2"||atom_type=="C.1"
      ||atom_type=="C.ar" ||atom_type=="C.2"||atom_type=="H"
      ||atom_type == "H.spc"||atom_type == "H.t3p"
      ||atom_type=="S.3"||atom_type=="S.2"||atom_type=="N.ar")
        return true;

    return false;
}

bool sasa::atom_is_hydrophilic(string atom_type){
    if (atom_type=="O.3"||atom_type=="O.2"||atom_type=="O.co2"
      ||atom_type=="O.spc" ||atom_type=="O.t3p" 
      ||atom_type=="F" ||atom_type=="Cl"||atom_type=="Br"||atom_type=="I"
      ||atom_type=="P.3"||atom_type == "Mn"
      ||atom_type== "Ca"|| atom_type== "Zn"||atom_type == "Mg"
      ||atom_type == "Na"||atom_type == "K"||atom_type == "Co" 
      ||atom_type=="N.3" ||atom_type=="N.2"||atom_type=="N.1"
      ||atom_type=="N.am" ||atom_type=="N.pl3"||atom_type=="N.4"
      ||atom_type == "C.cat"||atom_type == "S.O"||atom_type == "S.O2"
      ||atom_type == "C.cat"
      ) 
        return true;

    return false;
}

void
sasa::getSASA(ShortDOCKMol & molec, int recatomnum,
              float & sasa_tot,     float& sasa_hydrophilic,     float& sasa_hydrophobic,     float& sasa_other,
              float & sasa_lig_tot, float& sasa_lig_hydrophilic, float& sasa_lig_hydrophobic, float& sasa_lig_other,
              float & sasa_rec_tot, float& sasa_rec_hydrophilic, float& sasa_rec_hydrophobic, float& sasa_rec_other)
{
    sasa_tot = 0.0; sasa_hydrophobic = 0.0;
    sasa_hydrophilic = 0.0; sasa_other =0.0;

    sasa_lig_tot = 0.0; sasa_lig_hydrophobic = 0.0;
    sasa_lig_hydrophilic = 0.0; sasa_lig_other =0.0;

    sasa_rec_tot = 0.0; sasa_rec_hydrophobic = 0.0;
    sasa_rec_hydrophilic = 0.0; sasa_rec_other =0.0;

    double          xi,
                    yi,
                    zi;
    double          xj,
                    yj,
                    zj;

    // this can be change to a vector 

    int            *ineighbor = new int[molec.num_atoms * 200];
    int             ineighborpt = 0;
    int             icount;
    int             count;

    // consider making this vectors
    // vector <int> jj;
   // the array is already allocated dynamically, 
   // using a vector won't make it faster: sudipto 
    int            *jj = new int[molec.num_atoms];
    double         *r2x = new double[molec.num_atoms];

    double          dij;
    double          x[3];
    double         *const vdwrad = new double[molec.num_atoms];
    double          rgbmaxpsmax2 = 276.9428;

    // rgbmaxpsmax2 = GB radius max ... squared = (16.64)^2
    // why this arbitrary value: sudipto

    count = 0;

    icosa_init(2, 3, 0.0);

    for (int i = 0; i < molec.num_atoms; i++) {
        xi = molec.coord[i].v[0];
        yi = molec.coord[i].v[1];
        zi = molec.coord[i].v[2];
        x[0] = xi;
        x[1] = yi;
        x[2] = zi;
        icount = 0;

        for (int j = 0; j < molec.num_atoms; j++) {
            if (i != j) {
                xj = molec.coord[j].v[0];
                yj = molec.coord[j].v[1];
                zj = molec.coord[j].v[2];
                double          r2 =
                    (xi - xj) * (xi - xj) + (yi - yj) * (yi - yj) + (zi -
                                                                     zj) * (zi -
                                                                            zj);
                // cout << "i=" << i << " j=" << j << " r2=" << r2 << endl;

                if (r2 < rgbmaxpsmax2) {
                    jj[icount] = j;
                    r2x[icount] = r2;
                    icount++;
                }
            }
        }

        vdwrad[i] = molec.gb_radius[i] + 1.4;

        double          temp1,
                        temp2;
        for (int k = 0; k < icount; k++) {
            int             j = jj[k];
            dij = sqrt(r2x[k]);
            temp1 = molec.gb_radius[i] + 1.4;
            temp2 = molec.gb_radius[j] + 1.4;
            if ((temp1 + temp2) > dij) {
                ineighbor[count] = j;
                count++;
            } else {
            }
        }
        ineighbor[count] = -1;
        count++;
        /*
         * if(i == 0){ ineighborpt = 0; icosa_init(2, 3, 0.0); }
         */

        // ineighborpt is a interger value the defines the 
        // neightors 
        // is incremented (meny times) every time icosa_sphere_approx is called.
        // ineighborpt_(all) should be the same value to when the fuction is called
        // to t

        double sphere_approx = icosa_sphere_approx(molec, i, x, vdwrad, &ineighborpt, ineighbor);
        sasa_tot = sasa_tot + sphere_approx;

        if (atom_is_hydrophobic(molec.atom_types[i])){
          sasa_hydrophobic += sphere_approx;
        }
        else if (atom_is_hydrophilic(molec.atom_types[i])){ // not hydrophobic
          sasa_hydrophilic += sphere_approx;
        }
        else { // not hydrophobic or hydrophilic
          //cout << "###" << molec.atom_types[i] << "###" << endl;
          sasa_other += sphere_approx;
        }

        // caclculate just receptor portion of SASA.
        if (i < recatomnum ){
            sasa_rec_tot = sasa_rec_tot + sphere_approx;
            if (atom_is_hydrophobic(molec.atom_types[i])){
              sasa_rec_hydrophobic += sphere_approx;
            }
            else if (atom_is_hydrophilic(molec.atom_types[i])){ // not hydrophobic
              sasa_rec_hydrophilic += sphere_approx;
            }
            else { // not hydrophobic
              sasa_rec_other += sphere_approx;
            }
        }
        // calculate just ligand portion of SASA.
        // in the case of ligand recatomnum = 0.
        else {
            sasa_lig_tot = sasa_lig_tot + sphere_approx;
            if (atom_is_hydrophobic(molec.atom_types[i])){
              sasa_lig_hydrophobic += sphere_approx;
            }
            else if (atom_is_hydrophilic(molec.atom_types[i])){ // not hydrophobic
              sasa_lig_hydrophilic += sphere_approx;
            }
            else { // not hydrophobic
              sasa_lig_other += sphere_approx;
            }
        }
    }
    delete[]jj;
    jj = NULL;

    delete[]r2x;
    r2x = NULL;

    delete[]ineighbor;
    ineighbor = NULL;

    delete[]vdwrad;

    return;
}

double
sasa::icosa_sphere_approx(ShortDOCKMol & molec, int i, double x[],
                          double vdwrad[], int *ineighborpt, int ineighbor[])
{
    double          pi = 3.141592654;
    double          xi,
                    yi,
                    zi,
                    r,
                    total_sas,
                    sas,
                    sastmp;
    int             iatenvcnt = 10;
    int             iatenv[iemax];
    int             iatinclcnt[12];
    int             iatincl[12 * 200];
    double          v[12 * 3];
    bool            exc[12];
    float           sphere_approx;


    xi = x[0];
    yi = x[1];
    zi = x[2];
    r = vdwrad[i];
    total_sas = (4.0 * pi * r * r) / 20;
    sas = 0.0;

    icosa_atom_environment(ineighborpt, ineighbor, &iatenvcnt, iatenv);


    if (iatenvcnt > 0) {
        for (int m = 0; m < 12; m++) {
            v[3 * m] = icosa[m][0] * r;
            v[3 * m + 1] = icosa[m][1] * r;
            v[3 * m + 2] = icosa[m][2] * r;
            exc[m] = icosa_point_exclusion(molec, x, vdwrad, iatenvcnt, iatenv,
                                      &iatinclcnt[m], &iatincl[m * 200],
                                      xi + v[3 * m], yi + v[3 * m + 1],
                                      zi + v[3 * m + 2]);

        }

        for (int n = 0; n < 20; n++) {
            if (exc[ictri[n][0]] && exc[ictri[n][1]] && exc[ictri[n][2]]) {     // if 
                                                                                // statement 
                                                                                // checked 
                                                                                // with 
                                                                                // sander
                sas = sas + total_sas;
            } else {

                sastmp = icosa_patch_approx(molec, i, x, vdwrad,
                                            &v[3 * ictri[n][0]],
                                            &v[3 * ictri[n][1]],
                                            &v[3 * ictri[n][2]],
                                            &exc[ictri[n][0]],
                                            &exc[ictri[n][1]],
                                            &exc[ictri[n][2]],
                                            &iatinclcnt[ictri[n][0]],
                                            &iatincl[ictri[n][0] * iemax],
                                            &iatinclcnt[ictri[n][1]],
                                            &iatincl[ictri[n][1] * iemax],
                                            &iatinclcnt[ictri[n][2]],
                                            &iatincl[ictri[n][2] * iemax],
                                            iatenvcnt, iatenv, 20);
                sas = sas + sastmp;
            }
        }
    }

    sphere_approx = sas;
    return sphere_approx;
}

// icosa called
double
sasa::icosa_patch_approx(ShortDOCKMol & molec, int i, double x[],
                         double vdwrad[], double *v1, double *v2, double *v3,
                         bool * exc1, bool * exc2, bool * exc3,
                         int *iatinclcnt1, int *iatincl1, int *iatinclcnt2,
                         int *iatincl2, int *iatinclcnt3, int *iatincl3,
                         int iatenvcnt, int iatenv[], int nof_parts)
{
    double          norm_x,
                    norm_y,
                    norm_z,
                    normsum,
                    normi;
    double          v12[3],
                    v23[3],
                    v31[3];
    double          pi = 3.141592654;
    int             iatinclcnt12;
    int             iatinclcnt23;
    int             iatinclcnt31;
    int             iatincl12[iemax],
                    iatincl23[iemax],
                    iatincl31[iemax];
    bool            exc12,
                    exc23,
                    exc31;
    double          xi,
                    yi,
                    zi;
    double          r,
                    total_sas,
                    sas;
    int             nof_exc,
                    idecomp;
    xi = x[0];
    yi = x[1];
    zi = x[2];
    r = vdwrad[i];
    total_sas = (4.0 * pi * r * r) / nof_parts;
    if (nof_parts >= ipmax) {
        nof_exc = 0;
        if (*exc1) {
            nof_exc = nof_exc + 1;
        }                       /* else if(idecomp>0){ feed_decsasa(); } */
        if (*exc2) {
            nof_exc = nof_exc + 1;
        }                       /* else if(idecomp>0){ feed_decsasa(); } */
        if (*exc3) {
            nof_exc = nof_exc + 1;
        }                       /* else if(idecomp>0){ feed_decsasa(); } */
        sas = 0.0;
        if (nof_exc == 1) {
            sas = total_sas * 0.3333333333;
        } else if (nof_exc == 2) {
            sas = total_sas * 0.6666666666;
        } else if (nof_exc == 3) {
            sas = total_sas;
        }

    } else {
        norm_x = *v1 + *v2;
        v1 += 1;
        v2 += 1;
        norm_y = *v1 + *v2;
        v1 += 1;
        v2 += 1;
        norm_z = *v1 + *v2;
        normsum = (norm_x*norm_x) + (norm_y*norm_y) + (norm_z*norm_z);
        normi = 1 / sqrt(normsum);
        v12[0] = norm_x * normi * r;
        v12[1] = norm_y * normi * r;
        v12[2] = norm_z * normi * r;

        v2 -= 2;

        norm_x = *v2 + *v3;
        v2 += 1;
        v3 += 1;
        norm_y = *v2 + *v3;
        v2 += 1;
        v3 += 1;
        norm_z = *v2 + *v3;
        normsum = (norm_x*norm_x) + (norm_y*norm_y) + (norm_z*norm_z);
        normi = 1 / sqrt(normsum);
        v23[0] = norm_x * normi * r;
        v23[1] = norm_y * normi * r;
        v23[2] = norm_z * normi * r;

        v1 -= 2;
        v3 -= 2;

        norm_x = *v3 + *v1;
        v3 += 1;
        v1 += 1;
        norm_y = *v3 + *v1;
        v3 += 1;
        v1 += 1;
        norm_z = *v3 + *v1;

        v1 -= 2;
        v2 -= 2;
        v3 -= 2;
        normsum = (norm_x*norm_x) + (norm_y*norm_y) + (norm_z*norm_z);
        normi = 1 / sqrt(normsum);
        v31[0] = norm_x * normi * r;
        v31[1] = norm_y * normi * r;
        v31[2] = norm_z * normi * r;
        exc12 = icosa_point_exclusion(molec, x, vdwrad, iatenvcnt, iatenv,
                                      &iatinclcnt12, &iatincl12[0],
                                      xi + v12[0], yi + v12[1], zi + v12[2]);
        exc23 = icosa_point_exclusion(molec, x, vdwrad, iatenvcnt, iatenv,
                                      &iatinclcnt23, &iatincl23[0],
                                      xi + v23[0], yi + v23[1], zi + v23[2]);
        exc31 = icosa_point_exclusion(molec, x, vdwrad, iatenvcnt, iatenv,
                                      &iatinclcnt31, &iatincl31[0],
                                      xi + v31[0], yi + v31[1], zi + v31[2]);
        sas = 0.0;
        // triangel 1
        if (*exc1 && exc12 && exc31) {
            sas = sas + total_sas / 4.0;
        } else if ((nof_parts < ipmin) || *exc1 || exc12 || exc31) {
            sas = sas + icosa_patch_approx(molec, i, x, vdwrad, v1, v12, v31,
                                           exc1, &exc12, &exc31,
                                           iatinclcnt1, iatincl1,
                                           &iatinclcnt12, &iatincl12[0],
                                           &iatinclcnt31, &iatincl31[0],
                                           iatenvcnt, iatenv, 4 * nof_parts);
        }
        /*
         * else if(idecomp > 0){ feed_decsasa(); feed_decsasa();
         * feed_decsasa(); }
         */
        // triangel 2
        if (*exc2 && exc23 && exc12) {
            sas = sas + total_sas / 4.0;
        } else if ((nof_parts < ipmin) || *exc2 || exc23 || exc12) {
            sas = sas + icosa_patch_approx(molec, i, x, vdwrad, v2, v23, v12,
                                           exc2, &exc23, &exc12,
                                           iatinclcnt2, iatincl2,
                                           &iatinclcnt23, &iatincl23[0],
                                           &iatinclcnt12, &iatincl12[0],
                                           iatenvcnt, iatenv, 4 * nof_parts);
        }
        /*
         * else if(idecomp > 0){ feed_decsasa(); feed_decsasa();
         * feed_decsasa(); }
         */
        // triangel 3
        if (*exc3 && exc31 && exc23) {
            sas = sas + total_sas / 4.0;
        } else if ((nof_parts < ipmin) || *exc3 || exc31 || exc23) {
            sas = sas + icosa_patch_approx(molec, i, x, vdwrad, v3, v31, v23,
                                           exc3, &exc31, &exc23,
                                           iatinclcnt3, iatincl3,
                                           &iatinclcnt31, &iatincl31[0],
                                           &iatinclcnt23, &iatincl23[0],
                                           iatenvcnt, iatenv, 4 * nof_parts);
        }
        /*
         * else if(idecomp > 0){ feed_decsasa(); feed_decsasa();
         * feed_decsasa(); }
         */
        // triangel c
        if (exc12 && exc23 && exc31) {
            sas = sas + total_sas / 4.0;
        } else if ((nof_parts < ipmin) || exc12 || exc23 || exc31) {
            sas = sas + icosa_patch_approx(molec, i, x, vdwrad, v12, v23, v31,
                                           &exc12, &exc23, &exc31,
                                           &iatinclcnt12, &iatincl12[0],
                                           &iatinclcnt23, &iatincl23[0],
                                           &iatinclcnt31, &iatincl31[0],
                                           iatenvcnt, iatenv, 4 * nof_parts);
        }                       /* else if(idecomp > 0){ feed_decsasa();
                                 * feed_decsasa(); feed_decsasa(); } */
    }


    return sas;
}

/*
 * void sasa::feed_decsasa(){
 * 
 * }
 */


// this is the slowest function found in profiling
// called 37% if the time
bool
sasa::icosa_point_exclusion(ShortDOCKMol & molec, double x[], double vdwrad[],
                            int iatenvcnt, int iatenv[], int *iatinclcnt,
                            int *iatincl, double xp, double yp, double zp)
{
    int             j;
    double          xj,
                    yj,
                    zj,
                    r2,
                    vdw;
    bool            binside = true;
    *iatinclcnt = 0;
    for (int f = 0; f < iatenvcnt; f++) {
        j = iatenv[f];
        xj = molec.coord[j].v[0];
        yj = molec.coord[j].v[1];
        zj = molec.coord[j].v[2];
        //r2 = pow((xp - xj), 2) + pow((yp - yj), 2) + pow((zp - zj), 2);
        r2 = (xp - xj)*(xp - xj) + (yp - yj)*(yp - yj) + (zp - zj)*(zp - zj);
        vdw = molec.gb_radius[j] + 1.4;
        vdw = vdw * vdw;
        if (r2 < vdw) {
            binside = false;
            *iatincl = j;
            iatincl++;
            *iatinclcnt = *iatinclcnt + 1;
        }
    }
    return binside;
}

void
sasa::icosa_atom_environment(int *ineighborpt, int ineighbor[], int *iatenvcnt,
                             int iatenv[])
{
    *iatenvcnt = 0;
    while (ineighbor[*ineighborpt] != -1) {
        if (*iatenvcnt > iemax) {
            cout << "iatenvcnt > iemax in icosasurf" << endl;
            exit(EXIT_FAILURE);
        }
        iatenv[*iatenvcnt] = ineighbor[*ineighborpt];
        *ineighborpt = *ineighborpt + 1;
        *iatenvcnt = *iatenvcnt + 1;
    }
    *ineighborpt = *ineighborpt + 1;
}

void
sasa::icosa_init(int isasmin, int isasmax, double sasrad)
{
    double          angle = 11.5;
    double          axis[3] = { 1.0, 1.0, 1.0 };
    double          point[3];
    double          rotmat[3][3];
    gen_rot_mat(angle, axis, rotmat);
    for (int a = 0; a < 12; a++) {
        point[0] = icosasym[a][0];
        point[1] = icosasym[a][1];
        point[2] = icosasym[a][2];
        rot_point(point, rotmat);
        icosa[a][0] = point[0];
        icosa[a][1] = point[1];
        icosa[a][2] = point[2];
    }
    if (isasmin < 0 || isasmax < isasmin || sasrad < 0) {
        cout << "Wrong input parameters for icosa algorithm";
        cout << "isasmin: " << isasmin << endl;
        cout << "isasmax: " << isasmax << endl;
        cout << "sasrad: " << sasrad << endl;
        exit(EXIT_FAILURE);     // !call mexit(6,1)
    }
    ismin = isasmin;
    ismax = isasmax;
    srad = sasrad;

    ipmin = 20;
    for (int g = 0; g < ismin; g++) {
        ipmin = ipmin * 4;
    }
    ipmax = 20;
    for (int h = 0; h < ismax; h++) {
        ipmax = ipmax * 4;
    }
    // cout<<"ismin="<<ismin<<" ismax="<<ismax<<" srad="<<srad<<endl;
    // cout<<"ipmin = "<<ipmin<<" ipmax = "<<ipmax<<endl;
    // cout<<"icosa_init"<<endl;

}


void
sasa::gen_rot_mat(double angle, double axis[], double rotmat[][3])
{
    double          pi = 3.141592654;
    double          norm,
                    c,
                    s,
                    t;

    angle = (angle * pi) / 180.0;
    c = cos(angle);
    s = sin(angle);
    t = 1 - c;
    norm = 1 / (sqrt((axis[0]*axis[0]) + (axis[1]*axis[1]) + (axis[2]*axis[2])));
    axis[0] = axis[0] * norm;
    axis[1] = axis[1] * norm;
    axis[2] = axis[2] * norm;

    rotmat[0][0] = t * (axis[0]*axis[0]) + c;
    rotmat[1][0] = t * axis[0] * axis[1] - s * axis[2];
    rotmat[2][0] = t * axis[0] * axis[2] + s * axis[1];

    rotmat[0][1] = t * axis[0] * axis[1] + s * axis[2];
    rotmat[1][1] = t * (axis[1]*axis[1]) + c;
    rotmat[2][1] = t * axis[1] * axis[2] - s * axis[0];

    rotmat[0][2] = t * axis[0] * axis[2] - s * axis[1];
    rotmat[1][2] = t * axis[1] * axis[2] + s * axis[0];
    rotmat[2][2] = t * (axis[2]*axis[2]) + c;

}

void
sasa::rot_point(double point[], double rotmat[][3])
{
    double          ptmp[3];

    for (int a = 0; a < 3; a++) {
        ptmp[a] = 0.0;
        int             c = 0;
        for (int b = 0; b < 3; b++) {
            ptmp[a] = ptmp[a] + rotmat[b][a] * point[b];
            c++;
        }

    }

    for (int c = 0; c < 3; c++) {
        point[c] = ptmp[c];
    }
}

// ////////////////////////////////////////////////////////////////////////////////////
// / //////
// / Solvent Accessible Surface Area ////// 
// / //////
// / C++ implementation of icosasurf.f from amber8 ////// 
// / (C) Noel Carrascal, Sudipto Mukherjee, Robert C. Rizzo ////// 
// / SUNY Stony Brook //////
// / March 22, 2006 //////
// / //////
// ////////////////////////////////////////////////////////////////////////////////////
