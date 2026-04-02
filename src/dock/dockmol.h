#ifndef DOCKMOL_H
#define DOCKMOL_H

#include <iosfwd>
#include <string>
#include <utility>  // pair
#include <vector>
#include "utils.h"

#ifdef BUILD_DOCK_WITH_RDKIT
// RDKit additions to DOCK6
// by Guilherme D. R. Matos
// February 19, 2020
// RDKit Modules
// Necessary for RDKit use
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/Bond.h>
#include <DataStructs/ExplicitBitVect.h>
#endif

#ifdef BUILD_DOCK_WITH_RDKIT
// GDRM, Feb 19, 2020
// Function that finds atomic numbers (Z, number of protons)
int                     find_atomic_number( std::string atom_type );
// Function that selects RDKit-compatible bond types given the DOCKMol-stored data
RDKit::Bond::BondType   select_bond_type( std::string bond_order );
// Cleanup methods that allow RDKit to be used alongside DOCK6
void                    fixNitroSubstructureAndCharge(RDKit::RWMol &res, unsigned int atIdx);
void                    fixPhosphateSubstructureAndCharge(RDKit::RWMol &res, unsigned int atIdx);
void                    fixProtonatedAmidineSubstructureAndCharge(RDKit::RWMol &res, unsigned int atIdx);
void                    guessFormalCharges(RDKit::RWMol &res);
unsigned int            chkNoHNeighbNOx(RDKit::RWMol &res, RDKit::ROMol::ADJ_ITER atIdxIt, int &toModIdx);
bool                    cleanUpSubstructures(RDKit::RWMol &res);
#endif

// sudipto & trent: 15 Sep 2009
// we should break up the DockMol class into two components,
// a static unchanging part (like the AMBER parm file)
// and a dynamic part with the coords, score, active atom list, etc.
// this would save a lot of memory esp. during growth
// also save time from having to copy dockmol objects

// possibility: make just the dynamic class containing a reference
// to a DockMol object that stores all the static parameters

typedef         std::vector < float >TOR_LIST;

/**************************************/

class           FOOTPRINT_ELEMENT {

  public:
    std::string     resname;
    int             resid;
    double          vdw;
    double          es;
    int             hb;

    FOOTPRINT_ELEMENT();

};


/***********************************/
//LINGLING
class           PH4_ELEMENT {

  public:
    std::string     ph4_type;
    int             ph4_id;
    double          x_crd;
    double          y_crd;
    double          z_crd;
    double          v_x;
    double          v_y;
    double          v_z;
    double          ph4_radius;

    PH4_ELEMENT();

};


/***********************************/
class           DOCKMol {

  public:

    // Molecule General Info
    std::string     title;
    std::string     mol_info_line;
    std::string     comment1;
    std::string     comment2;
    std::string     comment3;
    std::string     score_text_data;
    std::string     energy = ""; //default value for energy
    std::string     simplex_text;
    std::string     mol_data;   // misc data
    std::string     hdb_data;   // hdb data to store internal energy 

    bool            bad_molecule; // errors in atom section of mol2 file
    

    // Molecule Size and Property Info
    int             num_atoms;
    int             num_bonds;
    int             num_residues; // this is the number of residues in this molecule
                                  // this is not the number in the receptor unless
                                  // the molecule is the receptor
    float           total_dsol; // kxr
    bool            flag_write_solvation; //teb
    
/*
    // use functions in filter class to assign these values.

    // Molecule Descriptors for database filter
    // these values are assigned by call in dock.cpp
    // to filter class the function assign these values
    // not though in that class this is sloppy and needs
    // to be cleaned up.

    unsigned int    rot_bonds;
    unsigned int    heavy_atoms;
    unsigned int    hb_donors;
    unsigned int    hb_acceptors;
    float           mol_wt;
    float           formal_charge;
    float           xlogp;
*/

    // array allocation variable
    bool            arrays_allocated;

    // atom information
    std::string    *atom_data;  // misc data

    float          *x;
    float          *y;
    float          *z;


    bool *flag_acceptor; // this is a h-bond acceptor
    bool *flag_donator; // this is true for polar hydrogen
    int  *acc_heavy_atomid; // atom id of acceptor that is connected to polar h
                           // equal to zero if not polar hydrogen

    float          *charges;
    std::string    *atom_types;
    std::string    *atom_names;

    // GA flags: for pruning conformers - CS
    bool           used;
    bool           parent = 0; // 0 if offspring, 1 if parent, default value of 0

    //footprint info jwu
    std::string    *atom_number;
    std::string    *atom_residue_numbers;
    std::string    *subst_names; // residue name


    std::string    *atom_color; // kxr205
    float          *atom_psol;  // kxr205
    float          *atom_apsol; // kxr205

    // DTM - 11-12-08 - array to store rigid segment info
    int                        *atom_segment_ids;
    int                        *bond_segment_ids; // CDS 04/09/15

    // covalent docking atom
    int            dummy1 ;                    // attachment point 
    int            dummy2 ;                    // attachment point , vector
    int            atomtag ;

    // bond information
    int            *bonds_origin_atom;
    int            *bonds_target_atom;
    std::string    *bond_types;

    // ring info
    bool           *atom_ring_flags;
    bool           *bond_ring_flags;

    // activation info
    int             num_active_atoms;
    int             num_active_bonds;
    bool           *atom_active_flags;
    bool           *bond_active_flags;
    bool           *bond_keep_flags; //CDS 08/2014
    bool           *atom_dnm_flag; // LEP
    bool            mol_dnm_flag; // JDB

    // AMBER atom type info
    bool            amber_at_assigned;
    int            *amber_at_id;
    float          *amber_at_radius;
    float          *amber_at_well_depth;
    int            *amber_at_heavy_flag;
    int            *amber_at_valence;
    int            *amber_at_bump_id;

    // Molecule Descriptors for database filter
    unsigned int    rot_bonds;
    unsigned int    heavy_atoms;
    float           mol_wt;
    float           formal_charge;
    unsigned int    hb_donors;
    unsigned int    hb_acceptors;
    float           xlogp;


    // Colored sphere Type info
    bool            chem_types_assigned;
    std::string    *chem_types;

    // Pharmacophore Type info
    bool            ph4_types_assigned;
    std::string    *ph4_types;

    // AMBER bond type info
    bool            amber_bt_assigned;
    int            *amber_bt_id;
    int            *amber_bt_minimize;
    int            *amber_bt_torsion_total;
    FLOATVec       *amber_bt_torsions;

    // GB radius & scale factor
    float          *gb_hawkins_radius;
    float          *gb_hawkins_scale;

    // Multi-Grid Ensemble
    int             grid_num;

    // Amber Score Identity
    std::string     amber_score_ligand_id;

    // Arrays obtained from within xlogp.cpp
    // Removed by sudipto and trent until xlogp works again
    int             number_of_H_Donors;
    int             number_of_H_Acceptors;
    int             *H_bond_donor;
    int             *H_bond_acceptor;

    // General scoring information
    float           current_score;
    float           internal_energy;
    std::string     current_data;  //components of score from scoring function
    std::string     primary_data;
    std::string     hbond_text_data;
    float           vdw_comp;
    float           es_comp;


    // Save specific descriptor score information - CS 06-05-16
    float            score_nrg;     // Energy Score
    float            score_mg_nrg;  // MG Score
    float            score_cont_nrg;// Continuous Score
    float            score_fps;     // Footprint Score
    float            score_ph4;     // Ph4 Score
    float            score_tan;     // Fingerprint Score
    float            score_hun;     // Hungarian Score
    float            score_vol;     // Volume Score


    // Parameters for Genetic Algorithm ranking - CS 06-06-16
    int             rank;
    float           fitness;
    float           crowding_dist;


    // Atom environments for gasteiger calculation
    std::string     *atom_envs;

    // for footprints score
    // std::vector < double > footprint_vdw;
    // std::vector < double > footprint_es;
    // std::vector < int >    footprint_hb;
    std::vector <FOOTPRINT_ELEMENT> footprints;
    // for ph4 score, LINGLING
    std::vector <PH4_ELEMENT> ph4;
    // Neighbor lists
    INTVec         *neighbor_list;
    // bool *ie_neighbor_list;
    INTVec         *atom_child_list;    // list of children atoms for each bond
                                        // in mol (2x#bonds- for
                                        // directionality)
    //atom_keys created ans stored for tanimoto calculation.  YZ
    std::vector <std::string> atom_keys_0;
    std::vector <std::string> atom_keys_1;
    std::vector <std::string> atom_keys_2;

    // Internal Energy members
    // bool use_internal_energy;
    // int att_exp;
    // int rep_exp;
    // float dielectric;

    // General DOCK Functions
    bool           Write_Mol2(DOCKMol &, std::ostream &);
    //bool           initialize_from_mol2(std::istream &);
    //bool           Write_Footprint(DOCKMol &, std::ostream &);
    void           input_parameters(Parameter_Reader & parm);
    void           initialize();
    void           allocate_arrays(int, int,int=1000); //jwu
    void           clear_molecule();
    std::vector < int > get_atom_neighbors(int);
    std::vector < int > get_bond_neighbors(int);
    int            get_bond(int, int);
    void           id_ring_atoms_bonds();
    void           find_rings(std::vector < int >, std::vector < int >,
                           std::vector < bool > &, std::vector < bool > &, int);
    bool           bond_is_rotor(int index) {return amber_bt_id[index] != -1;};
    std::vector < int > get_atom_children(int, int);
    float          get_torsion(int, int, int, int);
    void           set_torsion(int, int, int, int, float);
    void           set_angle(int , int , int , float );
    void           translate_mol(const DOCKVector & vec);
    void           translate_mol(const DOCKVector & vec, bool);
    void           rotate_mol(double mat[9], bool);
    //void           rotate_mol(double mat[9]);
    void           rotate_mol(double mat[3][3]);
    void           rotate_mol(double mat[3][3], bool);
    bool           atoms_are_one_three(int, int);
    bool           atoms_are_one_four(int, int);
    void           prepare_molecule();
    void           setcharges(const double* charges, double units_factor = 1.0);
    void           setxyz( const double * xyz );
    void           cout_information(); // jdb debug function

    #ifdef BUILD_DOCK_WITH_RDKIT

    // Descriptors and related stuff
    bool           fails_filt;
    int            num_stereocenters; // number of stereocenters
    int            num_spiro_atoms; // number of spiro centers
    double         clogp; // Wildman-Crippen log P
    double         tpsa;  // topological polar surface area
    int            num_arom_rings; // number of aromatic rings
    int            num_alip_rings; // number of aliphatic rings
    int            num_sat_rings;  // number of saturated rings
    std::string    smiles;         // SMILES string
    double         qed_score; // QED, druglikeness
    double         sa_score;  // Synthetic Accessibility
    double         esol; // Solubility
    int            pns;
    std::vector<std::string>    pns_name;   
    boost::dynamic_bitset<> MACCS;
    //for MPI, MACCS_size is required to get mpi to work. 
    int            MACCS_size;
    //ExplicitBitVect *MACCS;


    // Descriptor-driven de novo variables
    bool           fail_clogp;
    bool           fail_esol;
    bool           fail_tpsa;
    bool           fail_qed;
    bool           fail_sa;
    bool           fail_stereo;
    bool           fail_pains;
    //int            dropped_at_layer;

    // RDKit-related methods
    RDKit::ROMol   DOCKMol_to_ROMol(bool create_smiles);
    RDKit::Atom    *createAtom(unsigned int atomIdx, RDGeom::Point3D & pos, bool noImplicit);
    RDKit::Bond    *createBond(unsigned int bondIdx);
    void           addAtomsToRWMol(RDKit::RWMol & tmp);
    void           addBondsToRWMol(RDKit::RWMol & tmp);
    bool           isNitro(std::string atom_type, int idx);
    //bool           isPhosphate(std::string atom_type, int idx);
    void           createRDKitDescriptors( bool create_smiles );

    #endif

    void           operator=(const DOCKMol & mol);

    DOCKMol();
    DOCKMol(const DOCKMol &);
    ~DOCKMol();

};

/**************************************************/
class           RANKMol {
  public:
    DOCKMol mol;
    float           score;
    std::string     data;

    void            clear_molecule() {
        mol.clear_molecule();
        score = 0.0;
        data = "";
    };
    int             operator<(RANKMol x) const {
        return (score < x.score);
    };

};

/********************************************************************/
class           BREADTH_SEARCH {

  public:

    INTVec atoms;
    INTVec          nbrs;
    INTVec          nbrs_next;

    int             get_search_radius(DOCKMol &, int, int);

};

typedef         std::pair < float, DOCKMol > SCOREMol;
int             less_than_pair(SCOREMol a, SCOREMol b);

bool            Read_Mol2(DOCKMol &, std::istream &, bool, bool, bool);
bool            Read_Mol2_retain(DOCKMol &, std::istream &, bool, bool, bool);
bool            Write_Mol2(DOCKMol &, std::ostream &);
void            copy_molecule(DOCKMol &, const DOCKMol &);
void            copy_molecule_shallow(DOCKMol & , const DOCKMol & );
void            copy_crds(DOCKMol &, DOCKMol &);
void            transform(DOCKMol &, float rmat[3][3], DOCKVector, DOCKVector);

/*****************************************************************/


/********************************************************************
 *  * db2 file has atom section that defines the atom type its charge and desolvation value
 *  */
class           HDB_atom {
  public:
     int            atom_num;
     std::string    atom_name;
     std::string    atom_syb_type;
     int            atom_seg_num;
     int            atom_dock_type;
     float          q  ;
     float          p_desolv  ;
     float          a_desolv  ;
     float          sasa  ;

   //void HDB_atom::initialize(int,std::string,std::string,int,int,float,float,float,float)
     void           initialize(int,std::string,std::string,int,int,float,float,float,float);
     void           print_atom();
     HDB_atom();
     ~HDB_atom();
};

/********************************************************************
 *  * db2 file has bond section this defineds conectivity
 *  */
class           HDB_bond {
  public:
     int            bond_num;
     int            atom1_num;
     int            atom2_num;
     std::string    bond_name; // sybyl bond type

     void           initialize(int,int,int,std::string);
     HDB_bond();
     ~HDB_bond();
};


/********************************************************************
 * */
class           HDB_coordinate {
  public:
     int            coord_num;
     int            atom_num;
     int            seg_num;
     float          x  ;
     float          y  ;
     float          z  ;

     void           initialize(int, int, int, float, float, float);
     void           print_coord();
     HDB_coordinate();
     ~HDB_coordinate();
};

/********************************************************************
 * */
class           HDB_segment {
  public:
     int seg_num;
     int start_coord;
     int stop_coord;

     void           initialize(int, int, int);
     HDB_segment();
     ~HDB_segment();
};

/********************************************************************
 * */
class           HDB_conformer {

  public:
    int conf_num;
    int num_of_seg;
    bool broken;
    float internal_energy1;
    float internal_energy2;
    int * list_of_seg;
    //std::vector <int> list_of_seg;


     //void           initialize(int, int, float, int * );
     void           initialize(int, int, bool, float, float );
     HDB_conformer();
     ~HDB_conformer();
};


/********************************************************************/
class           HDB_Mol{

  public:
    char *           name;
    int              num_of_atoms;
    int              num_of_bonds;
    int              num_of_coords;
    int              num_of_rigid;
    int              num_of_segs;  
    int              num_of_confs; //branches
    // heirarchy data
    // db2 info when reading in
    HDB_atom *       atoms;
    HDB_bond *       bonds;
    HDB_coordinate * coords;
    HDB_coordinate * rigid; // like the anchor, this is the portion for orienting,  this is the unmoving portion of the molecule.  
    HDB_segment *    segs;
    HDB_conformer *  confs;
    // store most of db2 info in total_mol for scoring.
    //
    //void           initialize();
    HDB_Mol();
    ~HDB_Mol();
    void           clear_molecule();
};

#endif  // DOCKMOL_H

