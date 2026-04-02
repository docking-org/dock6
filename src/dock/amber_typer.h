//
#ifndef AMBER_TYPER_H
#define AMBER_TYPER_H 

#include <string>
#include <vector>
#include "utils.h"  // STRING20
class DOCKMol;
class Parameter_Reader;


/***********************************/
class           ATOM_TYPE_NODE {

  public:
    char            type[6];
    int             include;
    int             next_total;
    int             vector_atom;
    int             multiplicity;
    float           weight;
    std::vector < ATOM_TYPE_NODE > next;
};

/***********************************/
class           ATOM_TYPE {

  public:
    char            name[100];
    char            atom_model;
    float           radius;
    float           well_depth;
    int             heavy_flag;
    int             valence;
    int             bump_id;
    float           gbradius;
    float           gbscale;
    std::vector < ATOM_TYPE_NODE > definitions;

                    ATOM_TYPE();
                    virtual ~ ATOM_TYPE();
};

/***********************************/
class           ATOM_TYPER {

  public:
    std::vector < ATOM_TYPE > types;
    std::vector < int   >atom_types;

    void            get_vdw_labels(std::string fname, bool read_gb_parm);
    int             assign_vdw_labels(DOCKMol &, int);
    // float vdw_radius_from_type(int type){return types[type].radius;};
    // float well_depth_from_type(int type){return types[type].well_depth;};

};

/***********************************/
class           BOND_TYPE {

  public:
    char            name[100];
    int             drive_id;
    int             minimize;
    int             torsion_total;
    std::vector < float >torsions;
    ATOM_TYPE_NODE  definition[2];
};

/***********************************/
class           BOND_TYPER {

  public:
    std::vector < BOND_TYPE > types;
    std::vector < int   >flex_ids;
    int             total_torsions;

    void            get_flex_search(std::string fname);
    void            get_flex_labels(std::string fname);
    void            apply_flex_labels(DOCKMol &);

    bool            is_rotor(int);
};

/***********************************/
class           CHEM_TYPE {

  public:
    STRING20 name;              /* Member name */
    std::vector < ATOM_TYPE_NODE > definitions;      /* Member definitions */
    int             definition_total;   /* Number of definitions */

};

/***********************************/
class           CHEM_TYPER {

  public:
    std::vector < CHEM_TYPE > types;
    std::vector < int   >chem_type_ids;      /* Label definitions */
    int             total;      /* Number of members */

    void            get_chem_labels(std::string fname);
    void            apply_chem_labels(DOCKMol &);
};

/***********************************/
class           PH4_TYPE {

  public:
    STRING20 name;              /* Member name */
    std::vector < ATOM_TYPE_NODE > definitions;      /* Member definitions */
    int             definition_total;   /* Number of definitions */

};

/***********************************/
class           PH4_TYPER {

  public:
    std::vector < PH4_TYPE > types;
    std::vector < int>       ph4_type_ids;      /* Label definitions */
    int             total;      /* Number of members */

    void            get_ph4_labels(std::string fname);
    void            apply_ph4_labels(DOCKMol &);
};

/***********************************/
class           AMBER_TYPER {

  public:

    // parameter file locations
    std::string     vdw_defn_file;
    std::string     flex_defn_file;
    std::string     flex_drive_tbl;
    std::string     chem_defn_file;
    std::string     ph4_defn_file;
    char            atom_model;
    int             verbose;
    bool	    skip_verbose_flag;
    // amber atom and bond typing classes
    ATOM_TYPER      atom_typer;
    BOND_TYPER      bond_typer;
    CHEM_TYPER      chem_typer;
    PH4_TYPER       ph4_typer;

    void            initialize(bool read_vdw, bool read_gb_parm, bool use_chem, bool use_ph4, bool use_volume);
    void            input_parameters(Parameter_Reader & parm, bool read_vdw,
                                     bool use_chem, bool use_ph4, bool use_volume);
    void            assign_hbond_labels(DOCKMol &);
    void            prepare_molecule(DOCKMol &, bool read_vdw, bool use_chem, bool use_ph4, bool use_volume);
    void            prepare_for_torsions( DOCKMol & );
    float           getMW (DOCKMol & mol);

};

#endif  // AMBER_TYPER_H

