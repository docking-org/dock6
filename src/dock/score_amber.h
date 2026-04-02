// definition of classes Amber_Score and ConjGradReturnCode
// 
// Amber_Score defines continuum scoring for the AMBER energy.
// It is a DOCK Score class and is thus derived from class Base_Score.
// 
// ConjGradReturnCode defines a translator from NAB's conjgrad
// exit code to human readable text.
// 
#ifndef SCORE_AMBER_H
#define SCORE_AMBER_H

#include <iosfwd>
#include <string>
#include <vector>
// C linkage convention for NAB.
extern "C" {
#include "nab/nab.h"
#include "nab/defreal.h"
}
typedef REAL_T  Real;
#include "base_score.h"
class AMBER_TYPER;
class ConjGradReturnCode;
class DOCKMol;
class Parameter_Reader;


class Amber_Score:public Base_Score {

  public:

    Amber_Score();
    virtual ~ Amber_Score();    // TODO ~Base_Score should be virtual 
    bool            compute_score(DOCKMol & mol);
    void            initialize(AMBER_TYPER &);
    void            input_parameters(Parameter_Reader & parm,
                                     bool & primary_score,
                                     bool & secondary_score);
    std::string     output_score_summary(float current_score);

  private:

    // prevent copying - do not implement
    Amber_Score(Amber_Score const &initializer);

    // prevent assignment - do not implement
    Amber_Score & operator=(Amber_Score const &rightHandSide);

    float               calculate_amber_energy();
    int                 create_atomexpr_from_sphere_receptor_distance();
    bool                import_dockmol_ligand(DOCKMol &);
    void                read_movable_region(Parameter_Reader & parm);

    enum MovableRegion {
        DISTANCE,
        EVERYTHING,
        LIGAND,
        NAB_ATOM_EXPRESSION,
        NOTHING,
    };

    std::string         complex_move_atomexpr;
    Real                conv_criterion_rmsgrad;
    Real                distance_cutoff;
    std::string         gb_model;
    std::string         ligands_frcmod;
    MovableRegion       movable;    // flexibility indicator
    std::string         nab_minimization_options;
    std::string         nab_md_options;
    std::string         nab_energy_options;
    int                 nab_mme_initial_iteration;
    std::string         nonbonded_cutoff;
    int                 num_premd_min_cycles;
    int                 num_md_steps;
    int                 num_postmd_min_cycles;
    std::string         receptor_move_atomexpr;
    bool                skip_unprepped;
    std::string         temperature;

    int                 num_receptor_atoms;
    MOLECULE_T         *receptor;   // the NAB representation of the receptor
    Real                receptor_energy;
    std::string         receptor_file_prefix;
    std::string         receptor_pdb;
    std::string         receptor_prmtop;
    std::vector< Real > receptor_xyz;

    MOLECULE_T         *complex;    // the NAB representation of the complex
    Real                complex_energy;
    std::string         complex_file_prefix;
    std::string         complex_pdb;
    std::string         complex_prmtop;
    std::vector< Real > complex_xyz;

    MOLECULE_T         *ligand;     // the NAB representation of the ligand
    Real                ligand_energy;
    std::string         lig_file_prefix;
    std::string         ligand_pdb;
    std::string         ligand_prmtop;
    std::vector< Real > ligand_xyz;

    bool                verbose;    // emit NAB energy breakdown
    std::vector< Real > gradient_xyz;       // unused output from NAB's mme
    std::vector< Real > velocity_xyz;       // unused output from NAB's md

};


std::ostream & operator<<( std::ostream & os, ConjGradReturnCode const & code );

class ConjGradReturnCode {

  public:

    ConjGradReturnCode( int code );
    int getcode( ) const;

  private:

    int code;
};

#endif  // SCORE_AMBER_H
