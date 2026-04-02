#ifndef FILTER_H
#define FILTER_H

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include "dockmol.h"

using namespace std;

class Library_File;
// Header file for the database filter

class Filter
{
	public:
	
	bool use_database_filter;

    // use these for consistent formatting in output_score_summary
    static const std::string DELIMITER;
    static const int         FLOAT_WIDTH;
    static const int         STRING_WIDTH;


	// Descriptor ranges
	unsigned int max_heavy_atoms;
	unsigned int min_heavy_atoms;
	unsigned int max_rot_bonds;
	unsigned int min_rot_bonds;
	unsigned int max_hb_donors;
	unsigned int min_hb_donors;
	unsigned int max_hb_acceptors;
	unsigned int min_hb_acceptors;
	float max_formal_charge;
	float min_formal_charge;
	float min_xlogp;
	float max_xlogp;
	float min_molwt;
	float max_molwt;

    // Descriptors -- Variables added because MOL2 files were printed without
    // descriptor values. Values were calculated but not passed to DOCKMol in
    // the main loop.
    //unsigned int numHbDonors;
    //unsigned int numHbAcceptors;
    //float formalCharge;
    //float filterXlogP;
    //float filterMolWt;

    #ifdef BUILD_DOCK_WITH_RDKIT

    // Library paths
    string sa_fraglib_path;
    string PAINS_path;

    // important data for SA and PAINS calculations
    std::map<unsigned int, double> fragMap;
    std::map<std::string, std::string> PAINSMap;

    // Descriptor ranges
    float min_clogp;
    float max_clogp;
    float min_logs;
    float max_logs;
    float min_tpsa;
    float max_tpsa;
    float min_qed;
    float max_qed;
    float min_sa;
    float max_sa;
    float max_pns;
    //float tpsa;         // topological polar surface area calculated by RDKit
    unsigned int max_number_stereocenters; // calculated using RDKit
    unsigned int min_number_stereocenters; // calculated using RDKit
    unsigned int max_number_spiro_centers; // calculated using RDKit
    unsigned int min_number_spiro_centers; // calculated using RDKit

    #endif

	void input_parameters(Parameter_Reader & parm);
        void calc_descriptors(DOCKMol & mol);

	bool fails_filter(DOCKMol & mol);
	string get_descriptors(DOCKMol & mol);
    string output_score_summary(DOCKMol & mol);
	
};


#endif // FILTER_H
