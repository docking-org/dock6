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
class           Internal_Energy_Score:public Base_Score {

  public:

    float           rep_exp, att_exp;
    float           diel_screen;
    bool            use_ddd;

    //hbond add
    // count unsatisfied interations

                    Internal_Energy_Score();
                    virtual ~ Internal_Energy_Score();


    void            input_parameters(Parameter_Reader & parm,
                                     bool & primary_score,
                                     bool & secondary_score);

    void            initialize(AMBER_TYPER &);
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

