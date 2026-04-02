#ifndef _RDTYPER_H_
#define _RDTYPER_H_

#include <map>
#include <tuple>
#include <vector>
#include <string>
#include "dockmol.h"
#include "utils.h"
#include <GraphMol/GraphMol.h>

class DOCKMol;
//class Parameter_Reader;

// Auxiliary functions
std::tuple<double, double, double, double, double, double, double> get_adsParameters(std::string);
int getStructuralALERTS (RDKit::ROMol &);
int getStructuralALERTS (RDKit::ROMol &);
double  ADS(double, std::string);

// Class, atributes and methods
class RDTYPER {
public:
    //// Parameter file path
    //std::string sa_fraglib_path;
    
    //// Attributes
    //double clogp;
    //double tpsa;
    //double qed_score;
    //double sa_score;
    //double esol;
    //int num_arom_rings;
    //int num_alip_rings;
    //int num_sat_rings;
    //int num_spiro_atoms;
    //int num_stereocenters;
    //std::string smiles;

    //// Constructors and destructor
    //RDTYPER(); // default
    //RDTYPER(const RDTYPER &); // copy
    //~RDTYPER(); // destructor                         

    // Methods
    //void initialize_parameters(std::string);
    double calculate_qed_score(RDKit::ROMol &, double, double, std::vector <double>); 
    double calculate_sa_score(RDKit::ROMol &, int, int, std::map<unsigned int, double> &); 
    double calculate_esol(RDKit::ROMol &, double, bool);
    int get_pns(DOCKMol &, RDKit::ROMol &,std::map<std::string,std::string> &);
    int get_num_pns(RDKit::ROMol &, std::vector<std::string>&);
    void calculate_descriptors(DOCKMol &, std::map<unsigned int, double> &, bool, std::map<std::string, std::string> &);
    void calculate_descriptors_drive(DOCKMol &, std::map<unsigned int, double> , bool, std::map<std::string, std::string> );
    //void calculate_descriptors(DOCKMol &, std::map<unsigned int, double> &, bool);
    //void assign_descriptors(DOCKMol &);

};



#endif  // _RDTYPER_H_

