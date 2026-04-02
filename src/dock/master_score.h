//
#ifndef MASTER_SCORE_H
#define MASTER_SCORE_H 

#include "score.h"
//#include "score_gist.h"
#include "score_multigrid.h"
#include "score_footprint.h"
#include "score_ph4.h"
#include "score_volume.h"
#include "score_hbond.h"
#include "score_internal.h"
#include "score_descriptor.h"
#include "score_amber.h"
#include "score_chemgrid.h"
#include "score_solvent.h"
#include "score_sasa.h"

class           Master_Score {

  public:
    Base_Score                 c_bas; //TEB 2010-03-31
    Energy_Score               c_nrg;
    GIST_Score                 c_gist; //TED 2018-07-24
    Multigrid_Energy_Score     c_mg_nrg;
    Chemgrid_Score             c_cmg;      // kxr & bks scoring 
    Continuous_Energy_Score    c_cont_nrg;
    Footprint_Similarity_Score c_fps;
    Ph4_Score                  c_ph4; //LINGLING 2015-01-15
    //Volume_Score               c_volume;
    Hbond_Energy_Score         c_hbond;
    Internal_Energy_Score      c_int;
    Descriptor_Energy_Score    c_desc;
    Contact_Score              c_cnt;
    GB_Pairwise                c_gbsa;
    GB_Hawkins                 c_gbsa_hawkins;
    SASA_score                 c_sasa;
    Amber_Score                c_amber;

    Base_Score     *primary_score;
    Base_Score     *secondary_score;

    bool            primary_score_found;
    bool            secondary_score_found;

    bool            use_primary_score;
    bool            use_secondary_score;
    bool            use_score;
    bool            read_gb_parm;
    bool            primary_min;
    bool            secondary_min;
    bool            read_vdw;
    bool            use_chem;
    bool            use_ph4;
    bool            use_volume;
    bool            amber;
    bool            use_nrg_grid;
    bool            use_cnt_grid;

    // ensemble
    bool            ir_ensemble;

                    Master_Score();
                    virtual ~ Master_Score();
    void            input_parameters(Parameter_Reader & parm);
    void            initialize_all(AMBER_TYPER &, int argc, char **argv);
    void            close_all();
    bool            compute_primary_score(DOCKMol & mol);
    bool            compute_secondary_score(DOCKMol & mol);

};

#endif  // MASTER_SCORE_H

