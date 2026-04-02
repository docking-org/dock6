//
#ifndef SCORE_CHEMGRID_H
#define SCORE_CHEMGRID_H 

#include <string>
#include "base_grid.h"
#include "base_score.h"


#define VDWOFF -0.09
#define MAX_ATOM_REC 2000
#define MAX_ATOM_LIG 200

#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

/******************************************************************/
class           Shape_Filter:public Base_Score {

  public:

    FILE * shape_grid_in;
    std::string     shape_file_prefix;
    bool            use_shape_score;    // for distmap filter kxr

    // Shape Filter Globals kxr 
    int             xmax;
    int             ymax;
    int             zmax;
    int             grddiv;

    int             max_mismatch;       // corr. to bump_max kxr
    int             shape_filter;
    float           shape_score;

    void            input_parameters(Parameter_Reader & parm);
    void            initialize();

    void            read_shape_grid();
    bool            match_shape(DOCKMol & mol, bool bypass, bool shp_score);
    int             get_shape_score(DOCKMol & mol);

};
/*****************************************************************/
class           Chemgrid_Grid:public Base_Grid {

  public:

    //variables for grid functions
    FILE * grid_in;
    std::string     file_prefix;
    std::string     solv_file_prefix;
    std::string     hsolv_file_prefix;
    FILE           *solv_grid_in;
    FILE           *hsolv_grid_in;
    FILE           *rdsol_grid_in;
    int             dsize;      // For DelPhi grids kxr
    bool            hsolv_flag;  // boolian to use hydrogen desolv grid

    float           dspacing;   // DelPhi grids kxr
    float           phi_coords[3];
    float           phi_origin[3];      // DelPhi grid box origin kxr
    float           oldmid[3];  // DelPhi grid box center kxr
    int             dspan[3];   // kxr
    //int             solv_grd_ex[3]; // top corner of the desolv box
    float           solv_grd_ex[3]; // top corner of the desolv box
    int             solv_grd_N[3]; // grid dimensions 
    int             perang;


    XYZCRD          phi_corners[8];     // kxr
    int             phi_neighbors[8];   // kxr
    float           phi_cube_coords[3]; // kxr
    float           phi_int_coords[3];  // kxr
    int             phi_nearest_neighbor;       // kxr
    float           phi_x_min,
                    phi_x_max,
                    phi_y_min,
                    phi_y_max,
                    phi_z_min,
                    phi_z_max;  // kxr
    float           phi_space[3];
 
    int             rdxmax;
    int             rdymax;
    int             rdzmax;
    int             rdgrdiv;
    float          *avdw;
    float          *bvdw;
    float          *es;
    float          *dslx;
    float          *hdslx; // hydrogen desolvation grid, TEB 2020
    float          *dslb;
    int             atom_model;
    int             att_exp, rep_exp;

    float          *phi;

    void            read_phi_grid();
    void            read_chm_grid();
    void            read_odm_grid();
    void            read_rdsol_grid();
    void            read_solv_grid(bool);
        // DelPhi BKS interpolation data kxr
    void            phi_corner_coords();        // kxr
    bool            phi_box_boundaries(float x, float y, float z);      // kxr
    void            real_to_phi_coords(float x, float y, float z);
    float           interpolpot();
    void            find_phi_neighbors(float x, float y, float z);      // kxr
    int             find_phi_index(int x, int y, int z);
    float           interpolphi(float *pgrid);  // kxr
};

/*****************************************************************/
class           Chemgrid_Score:public Base_Score {

  public:

    Chemgrid_Grid *  chem_grid;
    Shape_Filter    shape;
    std::string     rdsol_file_prefix;
    std::string     atomic_contrib;
    bool            use_chemgrid_score; // for chemgrid scoring kxr
    bool            use_conf_entropy;   // for conf entropy from rot bonds kxr
//    bool            add_ligand_internal;// for adding ligand internal energy to score kx
    bool            use_solv_score;     // for solvent occlusion score kxr
    bool            redist_pos_desol;   // for redistributing positive desolvation
    bool            use_odm_score;      // for occlusion desolvation method kxr
    bool            total_ligand_dsolv; // for total of volume based dsolv kxr
    bool            use_recep_dsolv;    // for delphi calc. receptor desolvation kxr
    bool            write_atomic_energy;
    bool            use_delphi_score;   // for delphi elec. score kxr
    int             interpol_method;    // For choosing interpolation method
    bool            hsolv_flag;  // boolian to use hydrogen desolv grid

    
    // Receptor desolvation grid globals kxr 0306

    float           vdw_scale;
    float           es_scale;
    int             energy_score;
    bool            vdwA_allocated;
    bool            vdwB_allocated;
    float           vdw_component;
    float           es_component;
    float           rdsol_component;
    float           polsolv_component;
    float           apolsolv_component;
    float           conf_entropy_component;
    float           bulksolv_component;
    float           explsolv_component;
    float           shape_score;

                    Chemgrid_Score();
                    virtual ~ Chemgrid_Score();
    void            input_parameters(Parameter_Reader & parm,
                                     bool & primary_score,
                                     bool & secondary_score);
    void            input_parameters_main(Parameter_Reader & parm, std::string parm_head);
    void            initialize(AMBER_TYPER &);
    bool            compute_score(DOCKMol & mol);
    std::string     output_score_summary(float, float);



    // Solvent Occlusion grid BKS data kxr


    float           solv_score;
    float           psolv_score;
    float           apsolv_score;
    float           socc_component;
    bool            redist_positve_desolv(DOCKMol & mol);       // redist lig.
                                                                // desol BKS
                                                                // kxr
    bool            compute_ligand_desolvation(DOCKMol & mol, int atom_id);     // Lig. 
                                                                                // desolvation 

    bool            compute_ligand_desolvation_interpolate(DOCKMol & mol, int atom_id , float &, float &);     // Lig. 
                                                                                            // desolvation TEB add on 2020/02/11
    // Receptor desolvation kxr
    float           rdsol_score;
    float           solx_val;
    float           solb_val;
    bool            compute_rdsol_score(DOCKMol & mol, int atom_id);    // Rec. 
                                                                        // desolvation 
                                                                        // kxr
    bool            compute_odm_score(DOCKMol & mol, int atom_id);      // Occup. 
                                                                        // desolvation 
                                                                        // kxr
    // Shape score BKS kxr
    float           get_atom_shape_score(float x, float y, float z);
    float           compute_buried_fraction(DOCKMol & mol);
    bool            atom_buried(DOCKMol & mol, int atom_id);
    float           get_conf_entropy(DOCKMol & mol);


};

#endif  // SCORE_CHEMGRID_H

