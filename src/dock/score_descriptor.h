#include <string>
#include "dockmol.h"
#include "base_score.h"
#include "score.h"
#include "score_multigrid.h"
#include "score_footprint.h"
#include "score_ph4.h"
#include "fingerprint.h"
#include "hungarian.h"
#include "score_volume.h"
#include "score_gist.h"
#include "score_chemgrid.h"

class AMBER_TYPER;
class Bump_Grid;
class Contact_Grid;
class Energy_Grid;
class GIST_Grid;
//class GIST_Score;
class Parameter_Reader;


#define VDWOFF -0.09
#define MAX_ATOM_REC 2000
#define MAX_ATOM_LIG 200
#define SQR(x) ((x)*(x))
#define MAX(x, y) ((x) > (y) ? (x) : (y))


/*****************************************************************/
class           Descriptor_Energy_Score:public Base_Score {

  public:

    // Input parameters
    bool            desc_use_nrg;
    bool            desc_use_mg_nrg;
    bool            desc_use_cont_nrg;
    bool            desc_use_fps;
    bool            desc_use_ph4;       // add by LINGLING
    //bool          desc_use_mmgbsa;    // add later, mutually exclusive at top level

    bool            desc_use_tan;
    bool            desc_use_hun;
    bool            desc_use_volume;    // add by YUCHEN
    bool            desc_use_gist;      // add by TEB and LEP
    bool            desc_use_cmg;       // add by TEB DOCK3.7 score
    //bool          desc_use_molprop;   // add later, additive with other properties
    //bool          desc_use_sasa;      // add later, additive with other properties


    // Score classes
    Energy_Score                desc_c_nrg;
    Multigrid_Energy_Score      desc_c_mg_nrg;
    Continuous_Energy_Score     desc_c_cont_nrg;
    Footprint_Similarity_Score  desc_c_fps;
    Ph4_Score                   desc_c_ph4;
    Fingerprint                 desc_c_fing;
    Hungarian_RMSD              desc_c_hun;
    Volume_Score                desc_c_vol;
    GIST_Score                  desc_c_gist;
    Chemgrid_Score              desc_c_cmg;      // TEB added


    // Other parameters (not in classes)
    std::string     desc_fing_ref_filename;
    int             desc_fing_depth;
    DOCKMol         fing_ref_mol;

    std::string              desc_hun_ref_filename;
    DOCKMol                  hun_ref_mol;
    int                      hun_ref_heavy_atoms;
    float                    desc_hun_matching_coeff;
    float                    desc_hun_rmsd_coeff;
    std::pair <double, int>  desc_hun_result;


    // Weights
    int             desc_weight_nrg;
    int             desc_weight_mg_nrg;
    int             desc_weight_cont_nrg;
    int             desc_weight_fps;
    int             desc_weight_ph4;
    int             desc_weight_tan;
    int             desc_weight_hun;
    int             desc_weight_volume;
    int             desc_weight_gist;
    int             desc_weight_cmg;

    // Temp scores
    float     temp_nrg_score,
              temp_mg_nrg_score,
              temp_cont_nrg_score,
              temp_fps_score,
              temp_ph4_score,
              temp_tan_score,
              temp_hun_score,
              temp_vol_score,
              temp_gist_score,
              temp_cmg_score,
              temp_desc_score;


/*
    float           vdw_scale;
    float           es_scale;
    float           hb_scale;
    float           fp_vdw_scale;
    float           fp_es_scale;
    float           fp_hb_scale;
    float           internal_scale;
    float           rot_bonds_scale;
    float           heavyatoms_scale;
    float           hb_don_scale;
    float           hb_acc_scale;
    float           molecular_wt_scale;
    float           formal_charge_scale;
    float           XlogP_scale;
*/


    // Functions
    Descriptor_Energy_Score();
    virtual ~ Descriptor_Energy_Score();

    void            close();
    void            input_parameters(Parameter_Reader & parm, bool & primary_score, bool & secondary_score);
    void            initialize(AMBER_TYPER &);
    bool            compute_score(DOCKMol &);
    std::string     output_score_summary(DOCKMol &);
  private:
    bool desc_score_valid = false; //This variable is default to false and will be turned to true if a descriptor score component score is selected for use
};

/*****************************************************************/

