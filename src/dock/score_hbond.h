//

#include <string>
#include "base_score.h"
#include "dockmol.h"
class AMBER_TYPER;
class Bump_Grid;
class Contact_Grid;
class Energy_Grid;
class Parameter_Reader;


#define VDWOFF -0.09
#define MAX_ATOM_REC 1000
#define MAX_ATOM_LIG 200

#define SQR(x) ((x)*(x))
#define MAX(x, y) ((x) > (y) ? (x) : (y))



/*****************************************************************/
class           Hbond_Energy_Score:public Base_Score {

  public:

    std::string     receptor_filename;
    DOCKMol         receptor;
    float           vdw_scale;
    float           es_scale;
    float           hb_scale;
    float           internal_scale;

    float           rep_exp, att_exp;
    float           diel_screen;
    bool            use_ddd;

    float           vdw_component;
    float           es_component;
    //hbond add
    int             hbond;
    // count unsatisfied interations
    int             rec_lig_hb_don,rec_lig_hb_acc; 
    int             lig_rec_hb_don,lig_rec_hb_acc;
    int             rec_rec_hb_don,rec_rec_hb_acc; 
    int             lig_lig_hb_don,lig_lig_hb_acc;
    int             rec_hb_don,rec_hb_acc; // total
    int             lig_hb_don,lig_hb_acc; // total

    float           acc_don_rec_distance  ; // this is a parameter for ditermining the h-bond acceptors/doners on the receptor. (distance from ligand)
    float           threshold ; // this is the distance used for defining an h-bond
    float           min_angle ; // this is the minimum angle permited for defining an b-bond


                    Hbond_Energy_Score();
                    virtual ~ Hbond_Energy_Score();


    void            input_parameters(Parameter_Reader & parm,
                                     bool & primary_score,
                                     bool & secondary_score);

    void            initialize(AMBER_TYPER &);
    bool            hbond_cal_func1(DOCKMol &, DOCKMol &, int, int ,double &, double &);
    bool            compute_hbonds(DOCKMol &);
    void            close();
    bool            compute_score(DOCKMol &);
    std::string     output_score_summary(DOCKMol &);
};


/*****************************************************************/
// auxiliary functions
//void Tokenizer(std::string, std::vector < std::string > &, char);
// Tokenizer was moved to utils 
//bool hbond_cal(DOCKMol &, DOCKMol &, int, int ,double &, double &);
//void print_foot(std::vector <std::string> *,std::vector <double> *,std::vector <double> *, std::string *);

