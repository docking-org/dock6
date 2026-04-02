//
#ifndef SCORE_FOOTPRINT_H
#define SCORE_FOOTPRINT_H

#include <string>
#include "base_score.h"
#include "dockmol.h"
class AMBER_TYPER;
class Bump_Grid;
class Contact_Grid;
class Energy_Grid;
class Parameter_Reader;


#define VDWOFF -0.09
#define MAX_ATOM_REC 1000
#define MAX_ATOM_LIG 200

#define SQR(x) ((x)*(x))
#define MAX(x, y) ((x) > (y) ? (x) : (y))


/*****************************************************************/

class  RANGE {

  public:

    std::string     fps_range_vdw;
    std::string     fps_range_es;
    std::string     fps_range_hb;

  RANGE();

  RANGE & operator=( const RANGE & ranges) {
      fps_range_vdw = ranges.fps_range_vdw;
      fps_range_es  = ranges.fps_range_es;
      fps_range_hb  = ranges.fps_range_hb;
      return (*this);
  };
};

/*****************************************************************/
class           Footprint_Similarity_Score:public Base_Score {

  public:

    std::string     receptor_filename;
    DOCKMol         receptor;
    float           fp_vdw_scale;
    float           fp_es_scale;
    float           fp_hb_scale;

    float           rep_exp, att_exp;
    float           diel_screen;
    bool            use_ddd;

    float           vdw_component;
    float           es_component;
    int             hbond;

    double          vdw_foot_dist;
    double          es_foot_dist;
    double          hbond_foot_dist;
    std::string     fps_foot_compare_type;

    bool            bool_footref_mol2;
    bool            bool_footref_txt;
    bool            fps_foot_comp_all_residue;
    bool            fps_normalize_foot;
    
    RANGE           ref_ranges;
    RANGE           pose_ranges;
    std::string     fps_fp_info;

    // these are only used if fps_foot_comp_all_residue is false, i.e. not all residues are used.
    bool            fps_foot_specify_a_range; 
    bool            fps_foot_specify_a_threshold;
    bool            fps_use_remainder; 

    double          fps_vdw_threshold;
    double          fps_es_threshold;
    double          fps_hb_threshold;

    int             fps_vdw_num_resid;
    int             fps_es_num_resid;
    int             fps_hb_num_resid;

    DOCKMol         footprint_reference;
    bool            constant_footprint_ref;
    std::string     constant_footprint_ref_file;


    Footprint_Similarity_Score();
    virtual ~ Footprint_Similarity_Score();

    // ADD IN footprint reader the reads in a file in the 
    // same format as the output and stores the footprints
    // will read in h-bond, vdw, and es  
    // footprints are calculated outside dockrun.
   
    // one can pass a reference molecule to calculate the
    // reference footprints with in dockrun. 
    // footprint reference
    void            submit_footprint_reference(AMBER_TYPER &);
    bool            compute_footprint(DOCKMol &);
    double          calc_fp_similarity(DOCKMol &,char,char);
    double          Euclidean_distance(std::vector <double> *,std::vector <double> *); 
    double          Pearson_correlation(std::vector <double> *,std::vector <double> *); 
    // stop footprint

    void            input_parameters(Parameter_Reader & parm,
                                     bool & primary_score,
                                     bool & secondary_score);
    void            input_parameters_main(Parameter_Reader & parm, std::string parm_head);

    void            initialize(AMBER_TYPER &);
    void            close();
    bool            compute_score(DOCKMol &);
    std::string     output_score_summary(DOCKMol &);
    bool            on_range(int,std::string);
    RANGE           range_satisfying_threshold_all(std::vector <FOOTPRINT_ELEMENT>);
    std::string     union_of_ranges(std::string,std::string);
};


/*****************************************************************/
// auxiliary functions
bool read_footprint_txt(DOCKMol &, std::istream &);
bool hbond_cal(DOCKMol &, DOCKMol &, int, int ,double &, double &);
void range_satisfying_threshold(std::string &, double, std::vector <double>);
double Norm2(const std::vector <double> *);
void Normalization(std::vector <double> *);
void print_range(std::string, std::string, std::string, std::string, std::string *);
void print_foot(std::vector <std::string> *,std::vector <double> *,std::vector <double> *, std::string *);


#endif // SCORE_FOOTPRINT_H
