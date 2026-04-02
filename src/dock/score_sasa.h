#include <string>
#include <vector>
#include "base_score.h"
#include "grid.h"  // COORDS
#include "score.h"  // Energy_Score
#include "score_solvent.h"  // Energy_Score
#include "sasa.h"  


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class           SASA_score:public Base_Score {

  public:

    // parameters
    std::string     rec_filename;

    // data structures for solvation
    DOCKMol         receptor;
    ShortDOCKMol    short_rec,
                    short_lig,
                    short_com;


    bool            verbose;

    // SASA Portion
    float           rec_sasa_tot,
                    rec_sasa_phobic,
                    rec_sasa_philic,
                    rec_sasa_other,
                    lig_sasa_tot,
                    lig_sasa_phobic,
                    lig_sasa_philic,
                    lig_sasa_other,
                    com_sasa_tot,
                    com_sasa_phobic,
                    com_sasa_philic,
                    com_sasa_other,
                    com_sasa_lig_tot,
                    com_sasa_lig_phobic,
                    com_sasa_lig_philic,
                    com_sasa_lig_other,
                    com_sasa_rec_tot,
                    com_sasa_rec_phobic,
                    com_sasa_rec_philic,
                    com_sasa_rec_other;

    float           percent_lig_exposed,         
                    percent_lig_buried_is_phobic,                    
                    percent_rec_buried_is_phobic,
                    percent_buried_lig_phobic,
                    percent_buried_rec_phobic,
                    percent_lig_philic_exposed,
                    percent_lig_other_exposed,
                    percent_lig_phobic_exposed;

    //float           parm1, parm2, parm3;

    float           total_score;

    // Energy Scoring
                    SASA_score();
                    virtual ~ SASA_score();
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

