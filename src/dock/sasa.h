////////////////////////////////////////////////////////////////////////////////
// / //////
// / Solvent Accessible Surface Area //////
// / //////
// / C++ implementation of icosasurf.f from amber8 //////
// / (C) Noel Carrascal, Sudipto Mukherjee, Robert C. Rizzo //////
// / SUNY Stony Brook //////
// / //////
////////////////////////////////////////////////////////////////////////////////
//
#ifndef SASA_H
#define SASA_H

#include <iostream>
#include "score_solvent.h"

using namespace std;

class ShortDOCKMol;

class sasa {
public:
    double         getSASA(ShortDOCKMol &, int );
    void           getSASA(ShortDOCKMol &, int,
                           float & ,float& , float& , float&, 
                           float & ,float& , float& , float&, 
                           float & ,float& , float& , float& );
private:
    double        icosa[12][3];
    int           ismin;
    int            ismax;
    int            ipmin;
    int            ipmax;
    double        arad;       // Atomic radii
    double          srad;       // Water radius
    bool icosa_point_exclusion(ShortDOCKMol &, double[], double[], int, int[],
                                int *, int *, double, double, double);

    void           icosa_init(int, int, double);
    void           icosa_atom_environment(int *, int[], int *, int[]);
    double         icosa_sphere_approx(ShortDOCKMol &, int, double[], double[],
                                        int *, int[]);
    void           gen_rot_mat(double, double[], double[][3]);
    void           rot_point(double[], double[][3]);
    double         icosa_patch_approx(ShortDOCKMol &, int, double[], double[],
                                       double *, double *, double *, bool *,
                                       bool *, bool *, int *, int *, int *,
                                       int *, int *, int *, int, int[], int);
    void           feed_decsasa();
    bool           atom_is_hydrophobic(string atom_type);
    bool           atom_is_hydrophilic(string atom_type);

};


#endif  // SASA_H

