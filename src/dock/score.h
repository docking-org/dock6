//
#ifndef SCORE_H
#define SCORE_H 

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
class           Bump_Filter:public Base_Score {

  public:

    //---------------------------------------------//
    //CONSTRUCTOR AND DESTRUCTOR
    Bump_Filter();
    virtual ~ Bump_Filter();

    //---------------------------------------------//
    //VARIABLES
    Bump_Grid *     bump_grid;

    //flag to use bump filter
    bool            bump_filter;
    //name of grid file
    std::string     grid_file_name;
    //user defined parameter for maximum number of bumps for anchor
    int             anchor_bump_max;
    //user defined parameter for maximum number of bumps for growth
    int             growth_bump_max;

    //---------------------------------------------//
    //FUNCTIONS

    //user defined parameters
    void            input_parameters(Parameter_Reader & parm);
    //call to read and initialize grid
    void            initialize();
    //evaluate whether anchor should pass filter
    bool            check_anchor_bumps(DOCKMol & mol, bool);
    //evaluate whether conformation should pass filter
    bool            check_growth_bumps(DOCKMol & mol);
    //calculate bump score
    int             get_bump_score(DOCKMol & mol);

};

/*****************************************************************/
class           Energy_Score:public Base_Score {

  public:

    Energy_Grid *   energy_grid;

    std::string     grid_file_name;
    float           vdw_scale;
    float           es_scale;
    bool            grid_lig_efficiency;
    float           grid_total;
    float           vdw_component;
    float           es_component;
    float           vdw_eff;
    float           es_eff;
    float           total_eff;

                    Energy_Score();
                    virtual ~ Energy_Score();
    void            input_parameters(Parameter_Reader & parm,
                                     bool & primary_score,
                                     bool & secondary_score);
    void            initialize(AMBER_TYPER &);

    bool            compute_score(DOCKMol & mol);
    std::string     output_score_summary(DOCKMol & mol);

};

/*****************************************************************/
class           Continuous_Energy_Score:public Base_Score {

  public:

    std::string     receptor_filename;
    DOCKMol         receptor;

    float           vdw_scale;
    float           es_scale;
    float           rep_exp, att_exp;
    float           diel_screen;
    bool            use_ddd;

    float           vdw_component;
    float           es_component;

                    Continuous_Energy_Score();
                    virtual ~ Continuous_Energy_Score();
    void            input_parameters(Parameter_Reader & parm,
                                     bool & primary_score,
                                     bool & secondary_score);
    void            initialize(AMBER_TYPER &);
    bool            compute_score(DOCKMol &);
    std::string     output_score_summary(DOCKMol & mol);

};

/*****************************************************************/
class           Contact_Score:public Base_Score {

  public:

    Contact_Grid *  contact_grid;
    Bump_Grid *     bump_grid;
    std::string     grid_file_name;
    int             contact_score;
    float           cutoff_distance;
    float           clash_overlap;
    float           clash_penalty;

                    Contact_Score();
                    virtual ~ Contact_Score();
    void            input_parameters(Parameter_Reader & parm,
                                     bool & primary_score,
                                     bool & secondary_score);
    void            initialize(AMBER_TYPER &);
    bool            compute_score(DOCKMol & mol);
    std::string     output_score_summary(float);

};

#endif  // SCORE_H

