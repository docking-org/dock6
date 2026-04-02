//
#ifndef ORIENT_H
#define ORIENT_H 

#include <string>
#include <vector>
#include "dockmol.h"
#include "sphere.h"
#include "utils.h"
class Parameter_Reader;


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class           CLIQUE {
  public:
    INTVec nodes;
    double           residual;

    bool            operator<(CLIQUE clique) const {
        return (residual < clique.residual);
    };
};

typedef         std::vector < CLIQUE > CLIQUEVec;


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class           CLIQUE_EDGE {

  public:
    int             spheres[2];
    int             centers[2];
    double           residual;

    bool            operator<(CLIQUE_EDGE clnode) const {
        return (residual < clnode.residual);
    };
    bool            operator>(CLIQUE_EDGE clnode) const {
        return (residual > clnode.residual);
    };
    void            operator=(CLIQUE_EDGE cn) {
        spheres[0] = cn.spheres[0];
        spheres[1] = cn.spheres[1];
        centers[0] = cn.centers[0];
        centers[1] = cn.centers[1];
        residual = cn.residual;
    };
};


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class           CRITICAL_CLUSTER {

  public:
    int             index;
    INTVec          spheres;
};

typedef         std::vector < CRITICAL_CLUSTER > CLUSTERVec;


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
class           Orient {

  public:
    // General data members
    bool            orient_ligand;      // flag to use orienting routines
    DOCKMol         original;           // Original molecule
//    DOCKMol         cached_orient;      // Cached orientation
    bool            last_orient_flag;   // flag to control iterations
                                        // (false=more orients, true=last
                                        // orient)
    int             num_orients;        // number of orientations produced
    bool            critical_points;    // flag for use of critical point
                                        // matching
    bool            use_ligand_spheres;
    std::string     lig_sphere_filename;

    // Sphere & Center data
    SphereVec       centers;
    SphereVec       spheres;
    int             num_spheres;
    int             num_centers;
    int             num_nodes;
    double         *sph_dist_mat;
    double         *lig_dist_mat;
    double         *residual_mat;
//    double        **residual_mat;
    int             residual_mat_size;
    std::vector < CLIQUE_EDGE > edges;

    // matching data structures

/*
    bool           *candset;
    bool           *notset;
    int            *state;
    double         *level_residuals;
*/

    int             level;
    XYZVec          clique_spheres;
    XYZVec          clique_centers;
    int             clique_size;
    DOCKVector      spheres_com;
    DOCKVector      centers_com;
    double          rotation_matrix[3][3];

    CLIQUEVec       cliques;
    int             current_clique;

    // File I/O members
    std::string     sphere_file_name;   // Sphere infile name

    // Matching parameters
    double          tolerance;
    double          orig_tolerance;
    double          dist_min;
    int             min_nodes;
    int             max_nodes;
    int             max_orients;
    bool            automated_matching;
    int             am_iteration_num;

    // chemical matching parameters
    bool            use_chemical_matching;
    std::string     chem_match_tbl_fname;       // chem match table fname
    std::vector < std::string > chem_match_tbl_labels;    // chem type names
    INTVec          chem_match_tbl_matrix;      // table of allowed chem
                                                // type/chem type matchings
    INTVec          chem_match_align_tbl;       // table of allowed sphere/lig
                                                // atom matchings

    // critical sphere parameters
    CLUSTERVec      receptor_critical_clusters;

    bool            verbose;

    // Member functions
    Orient();
    virtual ~ Orient();
    void            input_parameters(Parameter_Reader & parm);  // Read in
                                                                // parameters
    void            initialize(int argc, char **argv);  // initialize
    void            prepare_receptor(); // read in receptor params
    void            get_spheres();      // read in receptor spheres
    void            get_lig_reference_spheres(DOCKMol &);      // generate spheres
                                                               // from a ref ligand
    void            calculate_sphere_distance_matrix(); // calculate matrix of
                                                        // sphere-sphere
                                                        // diststances
    void            match_ligand(DOCKMol &);    // generate information needed
                                                // for matching
    void            get_centers(DOCKMol &);     // Generate ligand centers from 
                                                // heavy atom positions
    void            calculate_ligand_distance_matrix(); // generate matrix of
                                                        // center-center
                                                        // distances
    void            clean_up(); // frees allocated memory
    bool            next_orientation(DOCKMol &);        // calls for a new
                                                        // orientation, manages 
                                                        // control flow for
                                                        // last orient
    bool            generate_orientation(DOCKMol &);    // generates a new
                                                        // orientation
    void            extract_coords_from_clique();       // extracts the spheres 
                                                        // and centers from a
                                                        // clique
    void            calculate_translations();
    void            translate_clique_to_origin();
    void            calculate_rotation();
    bool            more_orientations();
    void            id_all_cliques();   // test fxn to find all cliques
    bool            new_next_orientation(DOCKMol &);
    bool            new_next_orientation(DOCKMol &, bool);
    bool            orientation_HDB(DOCKMol &, bool); // TEB for HDB
    //bool            orientation_HDB(DOCKMol &); // TEB for HDB
    void            new_extract_coords_from_clique(CLIQUE &);
    bool            check_clique_critical_points(CLIQUE &);
    bool            check_clique_chemical_match(CLIQUE &);
    void            read_chem_match_tbl();

// added for covalent docking.
    //void            match_ligand_covalent(DOCKMol & mol, float lenth, float lenth2, float angle);
    void            match_ligand_covalent(DOCKMol & mol, float lenth, float lenth2, float angle, bool flag_adjust_bond, bool flag_orient, bool flag_first);
    //bool            new_next_orientation_covalent(DOCKMol & mol, float angle);
    void            new_next_orientation_covalent(DOCKMol & mol, float angle);
    void            set_torsion(float x0, float x1, float x2, float x3,
                                float y0 , float y1, float y2, float y3,
                                float z0 , float z1, float z2, float z3,
                                float angle, DOCKMol & mol);
};

#endif  // ORIENT_H

