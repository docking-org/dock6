//
#ifndef SCORE_SOLVENT_H
#define SCORE_SOLVENT_H 

#include <string>
#include <vector>
#include "base_score.h"
#include "grid.h"  // COORDS
#include "score.h"  // Energy_Score
class GB_Grid;
class SA_Grid;

#define VDWOFF -0.09
#define MAX_ATOM_REC 2000
#define MAX_ATOM_LIG 200

#define SQR(x) ((x)*(x))
#define MAX(x, y) ((x) > (y) ? (x) : (y))


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class           INT2000 {
  public:
    int             vals[2000];
};

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
typedef struct point_def {
    float           v[3];
} POINT;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class           ShortDOCKMol {

  public:

    int             num_atoms;
    POINT          *coord;
    float          *vdw_radius;
    float          *gb_radius;
    float          *gb_scale;
    float          *charges;
    std::string    *atom_names;
    std::string    *atom_types;

                    ShortDOCKMol();
                    ShortDOCKMol(DOCKMol const&);
                    ShortDOCKMol(ShortDOCKMol const&, ShortDOCKMol const&); //appends mol1+mol2
                    ShortDOCKMol(DOCKMol const& ,DOCKMol const& ); //convert to ShortDOCKMol
                                                                   //and appends mol1+mol2
                    virtual ~ ShortDOCKMol();

    void            allocate_short_arrays(int);

};

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class           GB_Pairwise:public Base_Score {

  public:

    // GB Variable Section ///////////////////////////////

    // Ligand data
    std::vector < COORDS > gb_lig_coords;
    FLOATVec        gb_lig_radius_eff;
    FLOATVec        gb_lig_sum_lig;
    FLOATVec        gb_lig_inv_a;
    INTVec          gb_lig_flag;

    bool            gb_initialized;
    // SA Variable Section //////////////////////////////
    // Ligand data
    std::vector < INT2000 > sa_grid_mark_sas_lig;
    std::vector < COORDS > sa_lig_coords;
    FLOATVec        sa_lig_radius_eff;                  

    // VDW calculation section
    Energy_Score    vdw_score;
    GB_Grid *       gb_grid;
    SA_Grid *       sa_grid;
    std::string     vdw_grid_prefix;

    // General Variables
    int             use_gbsa;
    std::string     gb_grid_prefix;
    std::string     sa_grid_prefix;
    std::string     screen_file;
    float           solvent_dielectric;
    bool            verbose;
    float           vdw_component;
    float           gb_component;
    float           sa_component;
    float           del_SAS_hp;

    // Functions ////////////////////////////////////////////

    GB_Pairwise();
    virtual ~ GB_Pairwise();

    void            input_parameters(Parameter_Reader & parm,
                                     bool & primary_score,
                                     bool & secondary_score);
    void            initialize(AMBER_TYPER & typer);
    bool            compute_score(DOCKMol & mol);
    std::string     output_score_summary(DOCKMol & mol);

    float           get_gb_solvation_score(DOCKMol & mol);
    float           get_sa_solvation_score(DOCKMol & mol);
};

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class           GB_Hawkins:public Base_Score {

  public:

    // parameters
    std::string     rec_filename;
    bool            salt_screen,
                    cont_vdw;
    float           gb_offset,
                    salt_conc,
                    espout;
    bool            verbose;

    // data structures for solvation
    DOCKMol         receptor;
    ShortDOCKMol    short_rec,
                    short_lig,
                    short_com;
    Energy_Score    vdw_score;

    // structures for GB
    float           kappa,
                    gfac;
    float           gpol_rec,
                    gpol_lig,
                    gpol_com;

    float          *born_radius_rec,
                   *born_radius_lig,
                   *born_radius_com;
    float          *rec_dist_mat,
                   *lig_dist_mat,
                   *com_dist_mat;

    // SASA Portion
    float           rec_sasa,
                    lig_sasa,
                    com_sasa;
    float           gnpol_rec,
                    gnpol_lig,
                    gnpol_com;
    // internal data for VDW & ES score
    float           vdw_component,
                    es_component;
    float           att_exp, rep_exp;
    float           total_score;

    // GB Portion
    void            born_array_calc(ShortDOCKMol &, float, float *, float *);
    float           gpol_calc(ShortDOCKMol &, float *, float *, bool, float,
                              float, float);

    // Energy Scoring
                    GB_Hawkins();
                    virtual ~ GB_Hawkins();
    void            input_parameters(Parameter_Reader & parm,
                                     bool & primary_score,
                                     bool & secondary_score);
    void            initialize(AMBER_TYPER &);

    void            prepare_receptor();
    bool            prepare_ligand(DOCKMol &);
    void            prepare_complex(DOCKMol &);

    bool            compute_score(DOCKMol &);
    std::string     output_score_summary(DOCKMol & mol);

};


#endif  // SCORE_SOLVENT_H

