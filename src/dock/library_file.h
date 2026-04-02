//
#ifndef LIBRARY_FILE_H
#define LIBRARY_FILE_H 

#include <fstream>
#include <string>
#include <vector>
#include "base_mpi.h"
#include "dockmol.h"
#include "amber_typer.h"
#include "gzstream/gzstream.h"
class Master_Score;
class Parameter_Reader;
class Simplex_Minimizer;
class Filter;
class RDTYPER;


class           Library_File:public Base_MPI {

  public:

    Library_File();

    // use these for consistent formatting in output_score_summary
    static const std::string DELIMITER;
    static const int         FLOAT_WIDTH;
    static const int         STRING_WIDTH;

    // data to identify database type (mol2 or hierarchy)
    int             db_type;    // 1=MOL2 (*.mol2), 2=HDB (*.db)
    std::vector < DOCKMol > mol_branches;    // vector for the
                                             // branches of an HDB molecule
    bool            db2flag;

    // general lib class data members
    bool            use_secondary_score;
    std::ifstream   ligand_in;
    std::ofstream   ligand_out_orients;
    std::ofstream   ligand_out_scored;
    std::ofstream   ligand_out_confs; 
    std::ofstream   ligand_out_ranked; 
    #ifdef BUILD_DOCK_WITH_RDKIT 
    std::ofstream   ligand_out_rejected;
    #endif
    std::ofstream   ligand_out_footprint_scored;
    std::ofstream   ligand_out_footprint_confs;
    std::ofstream   ligand_out_footprint_ranked;
    std::ofstream   ligand_out_hbond_scored;
    std::ofstream   ligand_out_hbond_confs;
    std::ofstream   ligand_out_hbond_ranked;
    std::ofstream   secondary_ligand_out_confs;
    std::ofstream   secondary_ligand_out_scored;
    std::ofstream   secondary_ligand_out_ranked;

    bool            max_mol_limit;
    int             max_mols;
    bool            read_solvation,
                    read_color;
    bool            write_solv_mol2; // for writing solv section to mol2 file. 
    int             initial_skip;
    int             total_mols;
    std::string     input_file_name;
    std::string     output_file_prefix;
    std::string     output_file_orient;
    std::string     output_file_confs;
    std::string     output_file_scored;
    std::string     output_file_ranked;
    #ifdef BUILD_DOCK_WITH_RDKIT 
    std::string     output_file_rejected;
    #endif
    std::string     output_file_footprint_confs;
    std::string     output_file_footprint_scored;
    std::string     output_file_footprint_ranked;
    std::string     output_file_hbond_confs;
    std::string     output_file_hbond_scored;
    std::string     output_file_hbond_ranked;
    std::string     secondary_output_file_confs;
    std::string     secondary_output_file_scored;
    std::string     secondary_output_file_ranked;

    DOCKMol         rmsd_reference;
    bool            calc_rmsd;
    int             num_anchors;
    int             num_orients;
    int             num_confs;
    bool            constant_rmsd_ref;
    std::string     constant_rmsd_ref_file;

    //DOCKMol         bestmol;
    //float           bestscore;
    //std::string     bestdata;
    float           score_thres;  // worst energy to print out.

    std::vector < RANKMol > ranked_poses;
    int             high_position;
    float           high_score;
    int             num_scored_poses;
    int             num_secondary_scored_poses;

    bool            write_footprints;
    bool            write_hbonds;
    bool            write_orients;
    bool            write_conformers;
    bool            write_secondary_conformers;
    bool            write_scored_ligands;

    bool            read_success;

    bool            rank_ligands;
    std::vector < SCOREMol > ranked_list;
    int             max_ranked;
    bool            rank_secondary_ligands;
    int             max_secondary_ranked;
    bool            rank_ligand_output_override;
    float           rank_list_highest_energy;
    int             rank_list_highest_pos;
    int             completed;

    // final pose RMSD clustering
    bool            cluster_ranked_poses;
    INTVec          cluster_assignments;  // maps ranked_poses to cluster size;
    float           cluster_rmsd_threshold;
    int             num_clusterheads_rescore;

    // secondary pose rescoring
    std::vector < DOCKMol > poses_for_rescore;
    std::vector < SCOREMol > ranked_secondary_list;
    

    void            input_parameters_input(Parameter_Reader & parm);
    void            input_parameters_output(Parameter_Reader &, Master_Score &, bool);

    void            initialize_input();
    void            initialize(int, char **, bool);
    void            open_files();
    void            close_files();

    bool            read_mol(DOCKMol &, bool);
    void            write_mol(DOCKMol &, std::ofstream &);
    bool            get_mol(DOCKMol &,bool, bool, bool, AMBER_TYPER &, Master_Score &, Simplex_Minimizer &);

    void            write_footprint_file(DOCKMol &, Master_Score &, std::ofstream &);
    void            write_hbond_file(DOCKMol &, Master_Score &, std::ofstream &);
//    void            clear_fp_hb(DOCKMol &);

    bool            submit_orientation(DOCKMol &, Master_Score &, bool);
//    void            submit_footprint(Master_Score &);
    void            submit_conformation(Master_Score &);
    void            submit_secondary_conformation(Master_Score &, Simplex_Minimizer & min);
    void            submit_scored_pose(DOCKMol &, Master_Score &,
                                       Simplex_Minimizer & min);
    void            submit_secondary_pose();
    void            rescore_conformation(DOCKMol &, Master_Score &,
                                       Simplex_Minimizer & min);
    void            secondary_rescore_poses(Master_Score &,
                                       Simplex_Minimizer &);
//    void            write_scored_poses(bool, Master_Score &);
    void            write_scored_poses(bool,bool, Master_Score &, Simplex_Minimizer &);
    void            sort_write(bool, bool, Master_Score &, Simplex_Minimizer & );
    void            write_ranked_ligands(bool,Master_Score &);
    //void            old_submit_scored_pose(DOCKMol &, Master_Score &);
    //void            old_write_scored_poses();

    void            submit_rmsd_reference(DOCKMol &, AMBER_TYPER &);
    void            calculate_rmsd(DOCKMol &,double *);
    void            calculate_rmsd(DOCKMol &, DOCKMol &, double *);
    float           calculate_std_rmsd(DOCKMol & refmol, DOCKMol & mol);
    float           calculate_min_rmsd(DOCKMol &, DOCKMol &);
    float           calculate_min_rmsd_one_way(DOCKMol &, DOCKMol &);
    float           calculate_min_cor_rmsd(DOCKMol &, DOCKMol &);
    std::string     calc_rmsd_string(DOCKMol & );
    void            cluster_list();
    void            calc_num_HBA_HBD( DOCKMol & mol );

    // hierarchy DB functions
    std::string     hdb_atom_type_converter(int);
    bool            read_hierarchy_db(DOCKMol &, std::ifstream &);
    //bool            read_hierarchy_db(DOCKMol &, igzstream &);
    //bool            read_hierarchy_db2(DOCKMol & mol, ifstream & ifs)
    //bool            read_hierarchy_db2(std::string, HDB_Mol &);
    //bool            read_hierarchy_db2(std::ifstream &, HDB_Mol &);
    bool            read_hierarchy_db2(igzstream &, HDB_Mol & );
    bool            read_sdifile(std::string, std::vector < std::string> & );

    //void            sort_top_X_mol(DOCKMol &, int );
    //void            sort_top_X_mol(std::vector <SCOREMol> &, int );
    //void            sort_top_X_mol(std::vector <SCOREMol> &, int , int*);
    bool            sort_top_X_mol(std::vector <SCOREMol> &, int , int*);

};

#endif  // LIBRARY_FILE_H

