//
#ifndef SCORE_VOLUME_H
#define SCORE_VOLUME_H

#include <string>
#include "base_score.h"
#include "dockmol.h"
class AMBER_TYPER;
class Bump_Grid;
class Contact_Grid;
class Energy_Grid;
class Parameter_Reader;


#define VDWOFF -0.09
#define MAX_ATOM_REC 2000
#define MAX_ATOM_LIG 200

#define SQR(x) ((x)*(x))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

/*****************************************************************/
class           Cloud_Struct {

  public:

    std::vector<float>   x;
    std::vector<float>   y;
    std::vector<float>   z;
    std::vector<std::string> type;

    void                 initialize();
    void                 clear();
};

/*****************************************************************/
// To calculate Volume overlap without using input parameters
// This focuses on the analytical method only
class           Volume_Score_Comp {

    public:
        float           total_component;
        float           heavy_atom_component;
        float           positive_component;
        float           negative_component;
        float           hydrophobic_component;
        float           hydrophilic_component;
        float           VOS;
  
        // reset the scores to null
        void            clear();
        Volume_Score_Comp compute_score_analytical(DOCKMol & , DOCKMol &);
        Volume_Score_Comp();
        ~Volume_Score_Comp();
        Volume_Score_Comp( const Volume_Score_Comp & );
        void operator=(const Volume_Score_Comp &);

        

};

/*****************************************************************/
class           Volume_Score:public Base_Score {

  public:

/*
    Energy_Grid *   energy_grid;

    std::string     grid_file_name;
    float           vdw_scale;
    float           es_scale;

    float           vdw_component;
    float           es_component;
*/


    std::string     volume_ref_file;
    std::string     tmp;
    DOCKMol         volume_ref_mol;
    float           separation;
    int             steps;
    std::string     volume_compare_type;

    float           overlap_ref_ref;
    float           overlap_ref_ref_hvy;
    float           overlap_ref_ref_pho;
    float           overlap_ref_ref_phi;
    float           overlap_ref_ref_neg;
    float           overlap_ref_ref_pos;
    float           total_component;
    float           heavy_atom_component;
    float           positive_component;
    float           negative_component;
    float           hydrophobic_component;
    float           hydrophilic_component;

    bool            output_cloud;

    Volume_Score();
    virtual ~ Volume_Score();

    void            input_parameters(Parameter_Reader & parm,
                                     bool & primary_score,
                                     bool & secondary_score);
    void            input_parameters_main(Parameter_Reader & parm, std::string parm_head);
    void            initialize(AMBER_TYPER &);

    bool            Grid_method(DOCKMol & mol);
    bool            Analytical_method(DOCKMol & mol,float &,float &,float &,float &,float &,float &);
    bool            compute_score(DOCKMol & mol);
    bool            compute_score(DOCKMol & mol, DOCKMol & ref_mol, AMBER_TYPER & c_typer);
    void            write_cloud(Cloud_Struct &);
    std::string     output_score_summary(DOCKMol & mol);

  private:
    //subroutines
    void            initialize_reference_overlap_helper();
    //light weight initialize function (not sure to make this public or private yet)
    void            initialize(AMBER_TYPER &, DOCKMol &);

};
#endif // SCORE_VOLUME_H

