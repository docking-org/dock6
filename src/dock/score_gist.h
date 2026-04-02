//

#include <string>
#include "base_score.h"
#include "dockmol.h"
class GIST_Grid;

/*****************************************************************/
class           GIST_Score:public Base_Score {

  public:

    GIST_Grid *     gist_grid;
    GIST_Grid *     gist_H_grid;

    std::string     grid_file_name;
    std::string     grid_H_file_name;
    float           gist_scale;
    //float           sigma2;  // used for blurry displacement. 
    std::string     gist_type;

    float           gist_component;
    float           att_exp;
    float           rep_exp;
    float           bgist_div;

                    GIST_Score();
                    virtual ~ GIST_Score();
    void            input_parameters(Parameter_Reader & parm,
                                     bool & primary_score,
                                     bool & secondary_score);
    void            input_parameters_main(Parameter_Reader & parm, std::string parm_head);
    void            initialize(AMBER_TYPER &);

    bool            compute_score(DOCKMol & mol);
    std::string     output_score_summary(DOCKMol & mol);

};


