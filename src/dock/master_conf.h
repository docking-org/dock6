#ifndef CONF_GEN_H
#define CONF_GEN_H 

#include <string>
#include <vector>
#include "dockmol.h"
#include "utils.h"   // INTVec
//#include "conf_gen_ag.h"
#include "conf_gen_cg.h"
#include "conf_gen_dn.h"
#include "conf_gen_ga.h"
#include "conf_gen_hdb.h"

class Bump_Filter;
class Master_Score;
class Parameter_Reader;
class Simplex_Minimizer;
class AG_Conformer_Search;
class DN_Build;


/********************************************************************/
class           Master_Conformer_Search {

  public:
    AG_Conformer_Search  c_ag_conf;
    CG_Conformer_Search  c_cg_conf;
    DN_Build             c_dn_build;          // object in de novo class
    GA_Recomb            c_ga_recomb;
    HDB_Conformer_Search  c_hdb_conf;

    bool            flexible_ligand;
    bool            genetic_algorithm;        // LEP need to remove some simpelx features
    bool            denovo_design;
    bool            last_conformer;
    int             method;
    std::string     conformer_search_type;    // for revamp of conf_gen.cpp
    bool            more_anchors;


    DOCKMol         orig;                     // original molecule struct

    bool            use_internal_energy;      // int nrg function superseded by func in base_score
    bool            initialize_once;          // this should be set true for every new ligand.

    // moved from AG_Conformer_Search 
    int             ie_att_exp;
    int             ie_rep_exp;
    float           ie_diel;
    float           ie_cutoff;

    void            input_parameters(Parameter_Reader & parm);
    void            initialize();
    void            prepare_molecule(DOCKMol &);
    bool            next_anchor(DOCKMol &);
    bool            submit_anchor_orientation(DOCKMol &, bool);
    void            grow_periphery(Master_Score &, Simplex_Minimizer &, Bump_Filter &);
    bool            next_conformer(DOCKMol &);

};

#endif  // CONF_GEN_H
