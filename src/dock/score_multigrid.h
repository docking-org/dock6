//
#ifndef SCORE_MULTIGRID_H
#define SCORE_MULTIGRID_H 

#include <string>
#include "base_score.h"
#include "dockmol.h"
class AMBER_TYPER;
class Bump_Grid;
class Contact_Grid;
class Energy_Grid;
class Multigrid_Energy_Grid;
class Parameter_Reader;

using namespace std;

#define VDWOFF -0.09
#define MAX_ATOM_REC 2000
#define MAX_ATOM_LIG 200
//#define bltzmn_const 0.0019872041 // kcal/mol/K

#define SQR(x) ((x)*(x))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

/*****************************************************************/
class           Multigrid_Energy_Score:public Base_Score {

  public:

    int numgrids;

    //vector < Energy_Grid * >  energy_grids;
    //vector < string >    grid_file_names;

    // array of pointers.
    //Energy_Grid * energy_grids[];
    Energy_Grid * energy_grids;
    string * grid_file_names;

    // scaler values to define the energy function.
    float           vdw_scale;
    float           es_scale;

    float           vdw_cor_scale;
    float           vdw_euc_scale;
    float           vdw_norm_scale;
//    float           vdw_sum_scale;
    float           es_cor_scale;
    float           es_euc_scale;
    float           es_norm_scale;
//    float           es_sum_scale;
    float           bltzmn_Z;  //BCF  scalar value of partition function
    float           bltzmn_kt; // BCF kT
    float           bltzmn_temp; // BCF kT

    // array
    float           * vdw_pose_array;
    float           * es_pose_array;
    float           * vdw_ref_array;
    float           * es_ref_array;
    float           * mgweights_array;
    float           * bltzmn_weight_array; //BCF

    DOCKMol         footprint_reference;
    //bool            constant_footprint_ref;

    // multigrid inputs 
    bool            fp_mol;
    bool            fp_txt;
    bool            wt_txt;
    bool            bltzmn;
    bool            ir_ensemble;
    string          footprint_ref_file;
    string          footprint_txt_file;
    string          gridweight_txt_file;

    // types of footprints
    bool            use_euc;
    bool            use_norm;
    bool            use_cor;

    // footprint comparison (similarity) value    
    float           vdw_cor;
    float           vdw_euc;
    float           vdw_norm;
    float           vdw_sum;
    float           es_cor;
    float           es_euc;
    float           es_norm;
    float           es_sum;

    // functions
                    Multigrid_Energy_Score();
                    virtual ~ Multigrid_Energy_Score();
    void            input_parameters(Parameter_Reader & parm,
                                     bool & primary_score,
                                     bool & secondary_score);
    void            input_parameters_main(Parameter_Reader & parm, std::string parm_head);
    bool            read_mgfootprint_txt(istream & );
    bool            read_mgweights_txt(istream & );
    void            submit_wts( );

    void            initialize(AMBER_TYPER &);
    void            close();

    bool            compute_multigrid(DOCKMol & ,float *, float *);
    void            submit_footprint_reference(AMBER_TYPER &);
    bool            compute_score(DOCKMol & mol);
    bool            compute_ir_ensemble_score(DOCKMol & mol);
    string          output_score_summary(DOCKMol & mol);
    void            partitionFunc(float &, float *, float *, float *, int, float);  //BCF
    float           bltzmn_sum(float *x, float *wt,int size);
};


#endif // SCORE_MULTIGRID_H
