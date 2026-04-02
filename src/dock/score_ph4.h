//
#ifndef SCORE_PH4_H
#define SCORE_PH4_H

#include <string>
#include "base_score.h"
#include "dockmol.h"
class AMBER_TYPER;
class Bump_Grid;
class Contact_Grid;
class Energy_Grid;
class Parameter_Reader;


#define VDWOFF -0.09

#define SQR(x) ((x)*(x))
#define MAX(x, y) ((x) > (y) ? (x) : (y))


/*****************************************************************/
class           Ph4_Struct {

  public:

    int            num_features;

    float          *x;
    float          *y;
    float          *z;
    float          *v_x;
    float          *v_y;
    float          *v_z;
    float          *radius;
    int            *id;
    std::string    *ph4_types;
    bool           *contri_to_score;

                    Ph4_Struct();
                    virtual ~ Ph4_Struct();
    void            initialize();
    void            clear();



};

/*****************************************************************/
class           Ph4_Score:public Base_Score {

  public:
    //parameters read specified from dock.in file
    bool            use_primary_score;
    bool            use_secondary_score;
    bool            bool_ph4_ref_mol2;
    bool            bool_ph4_ref_txt;
    bool            bool_out_ref_ph4_mol2;
    bool            bool_out_ref_ph4_txt;
    bool            bool_out_cad_ph4;
    bool            bool_out_mat_ph4;

    std::string     constant_ph4_ref_file;

    std::string     ref_ph4_out_mol2_filename;
    std::string     ref_ph4_out_txt_filename;
    std::string     cad_ph4_out_filename;
    std::string     mat_ph4_out_filename;
      
    std::string     ph4_compare_type;

    bool            ph4_full_match;

    bool            ph4_specify_mismatch_num;
    bool            ph4_specify_mismatch_percent;
    int             ph4_mismatch_num;
    int             ph4_mismatch_percent;

    double          ph4_rate_k;
    float           ph4_dist_r;
    double          ph4_proj_cos;
    int             ph4_max_x;

    //reference dockmol and reference ph4
    DOCKMol         reference;
    Ph4_Struct      Ph4_ref;

    //definition of ph4
    int             *match_comp;
    double          *components;
    int             *match_num;
    double          *match_term;

    std::vector <std::string> list_comp;

                    Ph4_Score();
                    virtual ~ Ph4_Score();

    //these are the functions defined in Ph4_Score
    bool            find_ring_member(DOCKMol & mol,
                                         int & id,
                           std::vector <int> & visited,
                           std::vector <int> & ring_member);
    void            close();
    void            input_parameters(Parameter_Reader & parm,
                                                 bool & primary_score,
                                                 bool & secondary_score);    
    void            input_parameters_main(Parameter_Reader & parm, std::string parm_head);
    void            initialize(AMBER_TYPER &);//prepare the reference ph4
    void            initialize(AMBER_TYPER & typer, DOCKMol & ref_mol); //overloaded initialize for dmax and more modular usage

    void            refresh_file(std::string);
    void            submit_ph4_mol2(Ph4_Struct &, std::string, std::string, bool);
    void            submit_ph4_txt(Ph4_Struct &, std::string);
    bool            compute_ph4(DOCKMol &, Ph4_Struct &);//compute ph4 from the mol2
    int             calc_num_match(Ph4_Struct &, Ph4_Struct &);
    double          calc_ph4_similarity(Ph4_Struct &, char);

    bool            compute_score(DOCKMol &);//output to .out file
    bool            compute_score(DOCKMol & mol, DOCKMol & ref_mol, std::string ph4_compare_type, AMBER_TYPER & c_typer); //alternative score computation meant for when the reference may be changing a lot
    std::string     output_score_summary(DOCKMol &);//output to the mol2
  private:
   //initialization subroutines
   void              initialize_prepare_Ph4_components_helper(AMBER_TYPER & typer);
   void              initialize_read_Mol2_helper(AMBER_TYPER & typer); 
   void              initialize_prepare_Mol2_helper(AMBER_TYPER & typer, DOCKMol & reference);
   void              initialize_read_ph4_helper(AMBER_TYPER & typer);
   void              initialize_output_helper();
   void              initialize_define_score_coeffecients_helper();
};
#endif // SCORE_PH4_H
