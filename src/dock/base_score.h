// definition of class Base_Score
// 
// This class defines the basic properties of all DOCK Score classes.
// Specific DOCK Score classes are to be derived from Base_Score.
// The characterization of a basic property is this equivalence:
// a property is basic to all Score classes if every Score class needs it
// and if a property is in Base_score then every Score class needs it and
// it is a basic property.
// 

#ifndef BASE_SCORE_H
#define BASE_SCORE_H

#include <string>
#include "utils.h" 
class AMBER_TYPER;
class DOCKMol;


class           Base_Score {

  public:

    // use these for consistent formatting in output_score_summary
    static const std::string DELIMITER;
    static const int         FLOAT_WIDTH;
    static const int         STRING_WIDTH;

    bool            use_score;
    bool            use_primary_score;
    bool            use_secondary_score;
    
    //general vdw parameters from vdw file
    float          *vdwA;
    float          *vdwB;
    float	   rep_radius_scale;

    //internal energy parameters
    int             method; // get from master_conformer_search 
                            // identify whether Rigid or Flex
    bool            use_internal_energy;  
    int             ie_att_exp;
    int             ie_rep_exp;
    float           ie_diel;
    float          *ie_vdwA;  //indexed by mol.num_atoms
    float          *ie_vdwB;  //indexed by mol.num_atoms

    std::vector <INTPair> nb_int; //intpair vector of non-bonded interactions

    Base_Score();
    virtual ~ Base_Score();

    // Generic Function Definitions 
    virtual bool    compute_score(DOCKMol &) {
        return true;
    };
    virtual std::string  output_score_summary(float) {
        return "No Score";
    };
 
    void            initialize_internal_energy(DOCKMol & mol);
    float           compute_ligand_internal_energy(DOCKMol & mol);     
    //float           compute_ligand_internal_energy(DOCKMol & mol, float& , float& , float& );     

    //float           compute_ligand_internal_energy_all_atom(DOCKMol & mol);
    //float           compute_ligand_internal_energy_all_atom(DOCKMol & mol, float& , float& , float&);
    void            init_vdw_energy(AMBER_TYPER &, float, float);

};

#endif  // BASE_SCORE_H

