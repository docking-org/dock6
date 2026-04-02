// definition of classes COORDS, Bump_Grid, Contact_Grid, Energy_Grid,
// GB_Grid, and SA_Grid

#ifndef GRID_H
#define GRID_H

#include <string>
#include <vector>
#include "base_grid.h"

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class           COORDS {
  public:
    float           crd[3];
};

typedef         std::vector < COORDS > CRDVec;


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class Bump_Grid is the bump grid from the grid program .bmp file.
// The bump grid is a community resource.
// This class uses the Singleton pattern:
// use Bump_Grid :: get_instance(std::string filename) to access the bump grid.
// But note that there is no enforcement of it being a singleton.
//
class       Bump_Grid:public Base_Grid {

public:
    
    Bump_Grid();
    void get_instance(std::string file_prefix);

private:

    void   read_bump_grid(std::string file_prefix);
    static bool got_the_grid;
};


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class Contact_Grid is the contact grid from the grid program .cnt file.
// The contact grid is a community resource.
// This class uses the Singleton pattern:
// use Contact_Grid :: get_instance(std::string filename) to access the grid.
// But note that there is no enforcement of it being a singleton.
//
class       Contact_Grid:public Base_Grid {

public:

    Contact_Grid();
    ~Contact_Grid();
    void get_instance(std::string filename);
    short int      *cnt; // in lieu of a better singleton implementation
                         // where cnt might be static.

private:

    static bool got_the_grid;
    void   read_contact_grid(std::string file_prefix);
};


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class Energy_Grid is the energy grid from the grid program .nrg file.
// The energy grid is a community resource.
// This class uses the Singleton pattern:
// use Energy_Grid :: get_instance(std::string filename)
// to access the energy grid.
// But note that there is no enforcement of it being a singleton.
//
class       Energy_Grid:public Base_Grid {

public:

    float  *avdw, *bvdw, *es;
    int     atom_model, att_exp, rep_exp;
    // 2011-01-07 -- trent e balius
    // move from private to public inorder to
    // reset got_the_grid so that more than 
    // one grid can be read in.
    static bool got_the_grid;

    // functions
    ~Energy_Grid();
    Energy_Grid();
    void   clear_grid();
    void get_instance(std::string filename);

private:

    void   read_energy_grid(std::string file_prefix);

    // functions

};


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class GIST_Grid is the receptor desolvation grid from the cpptraj program .dx file.
//
class       GIST_Grid:public Base_Grid {

public:

    //float  *avdw, *bvdw, *es;
    float  *gist;
    float   vol;  // volumn of a voxel
    //float   sigma2;  // for blurry_diplacement
    int     atom_model, att_exp, rep_exp;
    // 2011-01-07 -- trent e balius
    // move from private to public inorder to
    // reset got_the_grid so that more than
    // one grid can be read in.
    static bool got_the_grid;

    // functions
    ~GIST_Grid();
    GIST_Grid();
    void   clear_grid();
    void get_instance(std::string filename);
    float atomic_displacement(float, float, float, float, bool* & );
    float atomic_blurry_displacement(float, float, float, float, float );
    void write_gist_grid(std::string, bool* & );

private:

    void   read_gist_grid(std::string filename);

    // functions

};


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class GB_Grid is the generalized born grid from the nchemgrid_GB program 
// both the .bmp and .rec files are read
// The GB grid is a community resource.
// This class uses the Singleton pattern:
// use GB_Grid :: get_instance(std::string filename) to access the GB grid.
// But note that there is no enforcement of it being a singleton.
//
class       GB_Grid:public Base_Grid {

public:

    GB_Grid();
    void get_instance(std::string filename);
    int             gb_rec_total_atoms;
    float           gb_rec_solv_rec;
    std::vector < COORDS > gb_rec_coords;
    FLOATVec        gb_rec_radius_eff;
    FLOATVec        gb_rec_charge;
    FLOATVec        gb_rec_inv_a_rec;
    FLOATVec        gb_rec_inv_a;
    FLOATVec        gb_rec_receptor_screen;
    float           gb_grid_origin[3];
    float           gb_grid_spacing;
    float           gb_grid_grdcut;
    float           gb_grid_grdcuto;
    float           gb_grid_cutsq;
    float           gb_grid_f_scale;
    int             gb_grid_span[3];
    int             gb_grid_size;

private:

    static bool got_the_grid;
    void   read_gb_grid(std::string file_prefix);
};


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class SA_Grid is the sa grid from the nchemgrid_SA program 
// the .bmp, .sas, and .sasmark files are all read int
// The SA grid is a community resource.
// This class uses the Singleton pattern:
// use SA_Grid :: get_instance(std::string filename) to access the SA grid.
// But note that there is no enforcement of it being a singleton.
//
class       SA_Grid:public Base_Grid {

public:

    SA_Grid();
    void get_instance(std::string filename);
    int             sa_rec_total_atoms;
    float           sa_rec_solv_rec;
    std::vector < COORDS > sa_rec_coords;
    FLOATVec        sa_rec_radius_eff;
    INTVec          sa_rec_vdwn_rec;
    INTVec          sa_rec_nsphgrid;
    INTVec          sa_rec_nhp;
    int             sa_rec_nsas;
    int             sa_rec_nsas_hp;
    int             sa_rec_nsas_pol;
    float           sa_grid_origin[3];
    float           sa_grid_spacing;
    float           sa_grid_grdcut;
    float           sa_grid_grdcuto;
    float           sa_grid_cutsq;
    float           sa_grid_f_scale;
    int             sa_grid_span[3];
    int             sa_grid_size;
    int             sa_grid_nvtyp;
    std::vector < CRDVec > sa_grid_sphgrid_crd;
    INTVecVec       sa_grid_mark_sas;
    float           sa_grid_r_probe;
    float           sa_grid_sspacing;
    float           sa_grid_r2_cutoff;

private:

    static bool got_the_grid;
    void   read_sa_grid(std::string file_prefix);
};

#endif  // GRID_H
