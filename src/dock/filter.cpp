#include "filter.h"
#include "library_file.h"
#include <fstream>
#include <map>
//#include "xlogp.h"
//#include "fingerprint.h"

// RDKit stuff
#ifdef BUILD_DOCK_WITH_RDKIT
#include "rdtyper.h"
#include <DataStructs/ExplicitBitVect.h>
//#include <GraphMol/Descriptors/MolDescriptors.h>
#endif

//// These are the same as those in Base_Score.
//const string Filter::DELIMITER    = "########## ";
//const int    Filter::FLOAT_WIDTH  = 20;
//const int    Filter::STRING_WIDTH = 17 + 19;

using namespace std;



void Filter::input_parameters(Parameter_Reader & parm)
{
    cout << "\nDatabase Filter Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    use_database_filter = parm.query_param("use_database_filter", "no", "yes no") == "yes";

    if (use_database_filter)
        {
        max_heavy_atoms = atoi(parm.query_param( "dbfilter_max_heavy_atoms", "999" ).c_str());
        min_heavy_atoms = atoi(parm.query_param( "dbfilter_min_heavy_atoms", "0" ).c_str());
        max_rot_bonds = atoi(parm.query_param( "dbfilter_max_rot_bonds", "999" ).c_str());
        min_rot_bonds = atoi(parm.query_param( "dbfilter_min_rot_bonds", "0" ).c_str());
        max_hb_donors = atoi(parm.query_param( "dbfilter_max_hb_donors", "999" ).c_str());
        min_hb_donors = atoi(parm.query_param( "dbfilter_min_hb_donors", "0" ).c_str());
        max_hb_acceptors = atoi(parm.query_param( "dbfilter_max_hb_acceptors", "999" ).c_str());
        min_hb_acceptors = atoi(parm.query_param( "dbfilter_min_hb_acceptors", "0" ).c_str());
        max_molwt = atof(parm.query_param( "dbfilter_max_molwt", "9999.0" ).c_str());
        min_molwt = atof(parm.query_param( "dbfilter_min_molwt", "0.0" ).c_str());
        max_formal_charge = atof(parm.query_param( "dbfilter_max_formal_charge", "10.0" ).c_str());
        min_formal_charge = atof(parm.query_param( "dbfilter_min_formal_charge", "-10.0" ).c_str());
        //max_xlogp = atof(parm.query_param( "dbfilter_max_xlogp", "20.0" ).c_str());
        //min_xlogp = atof(parm.query_param( "dbfilter_min_xlogp", "-20.0" ).c_str());
        #ifdef BUILD_DOCK_WITH_RDKIT
        max_number_stereocenters = atoi(parm.query_param("dbfilter_max_stereocenters","6").c_str());
        min_number_stereocenters = atoi(parm.query_param("dbfilter_min_stereocenters","0").c_str());
        max_number_spiro_centers = atoi(parm.query_param("dbfilter_max_spiro_centers","6").c_str());
        min_number_spiro_centers = atoi(parm.query_param("dbfilter_min_spiro_centers","0").c_str());
        max_clogp = atof(parm.query_param("dbfilter_max_clogp","20.0").c_str());
        min_clogp = atof(parm.query_param("dbfilter_min_clogp","-20.0").c_str());
        max_logs = atof(parm.query_param("dbfilter_max_logs","20.0").c_str());
        min_logs = atof(parm.query_param("dbfilter_min_logs","-20.0").c_str());
        max_qed = atof(parm.query_param("dbfilter_max_qed","1.0").c_str());
        min_qed = atof(parm.query_param("dbfilter_min_qed","0.0").c_str());
        max_sa = atof(parm.query_param("dbfilter_max_sa","10.0").c_str());
        min_sa = atof(parm.query_param("dbfilter_min_sa","1.0").c_str());
        max_pns = atoi(parm.query_param("dbfilter_max_pns","100").c_str());
        sa_fraglib_path = parm.query_param("filter_sa_fraglib_path","sa_fraglib.dat");
        PAINS_path = parm.query_param("filter_PAINS_path","pains_table.dat");

        // Import fragMap
        if (fragMap.empty() == true) {
            std::ifstream fin(sa_fraglib_path);
            double key;
            double val;
            while (fin >> key >> val) {
                fragMap[key] = val;
            }
        }
        
        // Import PAINSMap 
        std::vector<std::string> PAINStmp;
        if (PAINStmp.empty() == true) {
            std::ifstream fin(PAINS_path);
            std::string tmp_string;
            while (fin >> tmp_string) {
                PAINStmp.push_back(tmp_string);
            }
            fin.close();
        }
        int numofvectors = PAINStmp.size() - 1;
        for (int i{0}; i < numofvectors; ) {
            PAINSMap[PAINStmp[i]] = " " + PAINStmp[i + 1];
            i += 2;
        }
        PAINStmp.clear();
          
        #endif

        // Turn these on once someone gets a chance to test it
        //Fingerprint              c_finger;
        //c_finger.input_parameters(parm); // if using database filter, see if user also wants
                                           // to filter by tanimoto to a ref
    }

}


/*string Filter::output_score_summary(DOCKMol & mol)
{
    ostringstream text;
    #ifdef BUILD_DOCK_WITH_RDKIT
    text << DELIMITER << setw(STRING_WIDTH) << "CLOGP: "
         << setw(FLOAT_WIDTH) << fixed << mol.clogp << endl
         << DELIMITER << setw(STRING_WIDTH) << "XLOGP: "
         << setw(FLOAT_WIDTH) << fixed << mol.xlogp << endl
         << DELIMITER << setw(STRING_WIDTH) << "TPSA: "
         << setw(FLOAT_WIDTH) << fixed << mol.tpsa << endl
         << DELIMITER << setw(STRING_WIDTH) << "Stereocenters: "
         << setw(FLOAT_WIDTH) << fixed << mol.number_stereocenters  << endl
         << DELIMITER << setw(STRING_WIDTH) << "Spiro_centers: "
         << setw(FLOAT_WIDTH) << fixed << mol.number_spiro_centers << endl;
    #else
    text << DELIMITER << setw(STRING_WIDTH) << "XLOGP: "
         << setw(FLOAT_WIDTH) << fixed << mol.xlogp << endl;
    #endif
    return text.str();
}*/


#ifdef BUILD_DOCK_WITH_RDKIT
void Filter::calc_descriptors(DOCKMol & mol)
{
    // Declare xlogp object
    //xlogp filterXlogP {};
    // Calculate xlogp and store in mol object - this changes atom types and charges
    // DO NOT uncomment until code is fixed - LEP
    //#ifdef BUILD_DOCK_WITH_RDKIT
    //int counter{1};
    //for (auto const& [key, val] : PAINSMap) {
    //    cout << counter << " " << key << " == " << val << endl;
    //    ++counter;
    //}


    RDTYPER rdprops;
    Library_File lib;

    if (!mol.bad_molecule){
        try {

            // Calculate all RDKit-related descriptors
            rdprops.calculate_descriptors( mol, fragMap, true, PAINSMap );
            
            // Calculate DOCK6 descriptors
            mol.heavy_atoms = 0;
            mol.rot_bonds = 0;
            mol.hb_acceptors = 0;
            mol.hb_donors = 0;

            // Count the number of heavy atoms
            for (int i = 0; i < mol.num_atoms; i++){
                 if (mol.amber_at_heavy_flag[i]) { mol.heavy_atoms++; }
            }
            // Count the number of rotatable bonds
            for (int i = 0; i < mol.num_bonds; i++){
                if (mol.bond_is_rotor(i)) { mol.rot_bonds++; }
            }
            lib.calc_num_HBA_HBD(mol);

        } catch (...) {
            mol.heavy_atoms = 0;
            mol.rot_bonds = 0;
            mol.hb_acceptors = 0;
            mol.hb_donors = 0;

            // Count the number of heavy atoms
            for (int i = 0; i < mol.num_atoms; i++){
                 if (mol.amber_at_heavy_flag[i]) { mol.heavy_atoms++; }
            }

            // Count the number of rotatable bonds
            for (int i = 0; i < mol.num_bonds; i++){
                if (mol.bond_is_rotor(i)) { mol.rot_bonds++; }
            }
            lib.calc_num_HBA_HBD(mol);
        }
    } else {
        // Print to error file
        cerr << "Molecule " << mol.title << " could not be used by RDKit.";
        cerr << " Bad MOL2 file. Inspection advised." << endl;
    }

    string    message;
    message = get_descriptors(mol);
    cout << message << endl;
}
#else
void Filter::calc_descriptors(DOCKMol & mol)
{
    // Declare xlogp object
    //xlogp filterXlogP {};
    // Calculate xlogp and store in mol object - this changes atom types and charges
    // DO NOT uncomment until code is fixed - LEP
    Library_File lib;
    mol.heavy_atoms = 0;
    mol.rot_bonds = 0;
    mol.hb_acceptors = 0;
    mol.hb_donors = 0;

    // Count the number of heavy atoms
    for (int i = 0; i < mol.num_atoms; i++)
         if (mol.amber_at_heavy_flag[i]) mol.heavy_atoms++;

    // Count the number of rotatable bonds
    for (int i = 0; i < mol.num_bonds; i++)
         if (mol.bond_is_rotor(i)) mol.rot_bonds++;

    // Calculate the number of HBA and HBD
    lib.calc_num_HBA_HBD(mol);

    //mol.current_data = output_score_summary(mol);
    string    message;
    message = get_descriptors(mol);
    cout << message << endl; 
}
#endif


// Checks if mol passes the filter with the max and min limits
bool Filter::fails_filter(DOCKMol & mol)
{
        bool fails = false;
        ostringstream strstrm;

        // Is it a bad input file?
        if (mol.bad_molecule)
        {
            strstrm << "Filter Rejected: Bad MOL2 File generated this molecule." << endl;
            fails = true;
        }

        // Rotatable Bonds
        if (mol.rot_bonds > max_rot_bonds)
        {
                strstrm << "Filter Rejected: Too many rotatable bonds: "
                                << mol.rot_bonds << " > " << max_rot_bonds << endl;
                fails = true;
        } else if (mol.rot_bonds < min_rot_bonds)
        {
                strstrm << "Filter Rejected: Not enough rotatable bonds: "
                                << mol.rot_bonds << " < " << min_rot_bonds << endl;
                fails = true;
        }


        // Heavy atoms
        if (mol.heavy_atoms > max_heavy_atoms)
        {
                strstrm << "Filter Rejected: Too many heavy atoms: "
                                << mol.heavy_atoms << " > " << max_heavy_atoms << endl;
                fails = true;
        } else if (mol.heavy_atoms < min_heavy_atoms)
        {
                strstrm << "Filter Rejected: Not enough heavy atoms: "
                                << mol.heavy_atoms << " < " << min_heavy_atoms << endl;
                fails = true;
        }

        // Hydrogen Bond Donors
        if (mol.hb_donors > max_hb_donors)
        {
                strstrm << "Filter Rejected: Too many hydrogen bond donors: "
                                << mol.hb_donors << endl;
                fails = true;
        } else if (mol.hb_donors < min_hb_donors)
        {
                strstrm << "Filter Rejected: Not enough hydrogen bond donors: "
                                << mol.hb_donors << endl;
                fails = true;
        }

        // Hydrogen Bond Acceptors
        if (mol.hb_acceptors > max_hb_acceptors)
        {
                strstrm << "Filter Rejected: Too many hydrogen bond acceptors: "
                                << mol.hb_acceptors << endl;
                fails = true;
        } else if (mol.hb_acceptors < min_hb_acceptors)
        {
                strstrm << "Filter Rejected: Not enough hydrogen bond acceptors: "
                                << mol.hb_acceptors << endl;
                fails = true;
        }

        // Molecular Weight
	// LEP setprecsion - ideally fixed should go before setprecision
        if (mol.mol_wt > max_molwt)
        {
                strstrm << "Filter Rejected: Molecular Weight too high: "
                                << setprecision(2) << fixed
                                << mol.mol_wt << " > " << max_molwt << endl;
                fails = true;
        } else if (mol.mol_wt < min_molwt)
        {
                strstrm << "Filter Rejected: Molecular Weight too low: "
                                << setprecision(2) << fixed
                                << mol.mol_wt << " < " << min_molwt << endl;
                fails = true;
        }

        // Formal Charge
        if ( (mol.formal_charge - max_formal_charge) > 0.1 )
        {
                strstrm << "Filter Rejected: Formal charge too high: "
                                << setprecision(2) << fixed
                                << mol.formal_charge << " > " << max_formal_charge
                                << endl;
                fails = true;
        } else if ( (mol.formal_charge - min_formal_charge) < -0.1 )
        {
                strstrm << "Filter Rejected: Formal charge too low: "
                                << setprecision(2) << fixed
                                << mol.formal_charge << " < " << min_formal_charge
                                << endl;
                fails = true;
        }
        //// XlogP
        //// XlogP code has been removed from the filter because it crashed for
        //// molecules with hexavalent sulphur atoms
        //if (mol.xlogp > max_xlogp)
        //{
        //        strstrm << "Filter Rejected: xlogp is too high: "
        //                        << setprecision(2) << fixed
        //                        << mol.xlogp << endl;
        //        fails = true;
        //} else if (mol.xlogp < min_xlogp)
        //{
        //        strstrm << "Filter Rejected: xlogp is too low: "
        //                        << setprecision(2) << fixed
        //                        << mol.xlogp << endl;
        //        fails = true;
        //}

        //-------------------------------------------------------------
        // GDRM Feb 25, 2020
        #ifdef BUILD_DOCK_WITH_RDKIT
        // If it is a bad molecule, no RDKit::RWMol object will be created
        if (mol.bad_molecule)
        {
                strstrm << "No RDKit-related properties were calculated." << endl;
                fails = true;
        } else {
        // number of stereocenters
            if (mol.num_stereocenters > max_number_stereocenters)
            {
                    strstrm << "Filter Rejected: too many stereocenters: "
                                << mol.num_stereocenters << endl;
                    fails = true;

            } else if (mol.num_stereocenters < min_number_stereocenters)
            {
                    strstrm << "Filter Rejected: too few stereocenters: "
                                << mol.num_stereocenters << endl;
                    fails = true;
            }

            // number of spiro centers
            if (mol.num_spiro_atoms > max_number_spiro_centers)
            {
                    strstrm << "Filter Rejected: too many spiro atoms: "
                                << mol.num_spiro_atoms << endl;
                    fails = true;

            } else if (mol.num_spiro_atoms < min_number_spiro_centers)
            {
                    strstrm << "Filter Rejected: too few spiro atoms: "
                                << mol.num_spiro_atoms << endl;
                    fails = true;
            }

            // CLogP (Wildman-Crippen log P)
            if (mol.clogp > max_clogp)
            {
                    strstrm << "Filter Rejected: clogp is too high: "
                                << setprecision(2) << fixed
                                << mol.clogp << endl;
                    fails = true;
            } else if (mol.clogp < min_clogp)
            {
                    strstrm << "Filter Rejected: clogp is too low: "
                                << setprecision(2) << fixed
                                << mol.clogp << endl;
                    fails = true;
            }
        
            // ESOL (Estimated SOLubility calculator for log S)
            if (mol.esol > max_logs)
            {
                    strstrm << "Filter Rejected: logS is too high: "
                                << setprecision(2) << fixed
                                << mol.esol << endl;
                    fails = true;
            } else if (mol.esol < min_logs)
            {
                    strstrm << "Filter Rejected: logS is too low: "
                                << setprecision(2) << fixed
                                << mol.esol << endl;
                    fails = true;
            }

            // QED (druglikeness)
            if (mol.qed_score > max_qed)
            {
                    strstrm << "Filter Rejected: QED is too high: "
                                << setprecision(2) << fixed
                                << mol.qed_score << endl;
                    fails = true;
            } else if (mol.qed_score < min_qed)
            {
                    strstrm << "Filter Rejected: QED is too low: "
                                << setprecision(2) << fixed
                                << mol.qed_score << endl;
                    fails = true;
            }

            // SYNTHA (a.k.a. SA_Score)
            if (mol.sa_score > max_sa)
            {
                    strstrm << "Filter Rejected: SYNTHA is too high: "
                                << setprecision(2) << fixed
                                << mol.sa_score << endl;
                    fails = true;
            } else if (mol.sa_score < min_sa)
            {
                    strstrm << "Filter Rejected: SYNTHA is too low: "
                                << setprecision(2) << fixed
                                << mol.sa_score << endl;
                    fails = true;
            }

            // Number of PAINS groups
            if (mol.pns > max_pns)
            {
                    strstrm << "Filter Rejected: #PAINS is too high: "
                                << setprecision(2) << fixed
                                << mol.pns << endl;
                    fails = true;
            }


        }
        #endif

        //cout << strstrm.str();

        // save the reason for filtering ligand so that MPI docking does not print
        // "Could not score conformer" whenever something is filtered
        if (fails) mol.current_data = strstrm.str();
        if (fails) mol.primary_data = strstrm.str();
        return fails;
}



string Filter::get_descriptors(DOCKMol & mol)
{
        ostringstream strstrm;
        //preparing pains names
        #ifdef BUILD_DOCK_WITH_RDKIT 
        std::vector<std::string> vecpnstmp{};
        vecpnstmp = mol.pns_name;
        std::string molpns_name = std::accumulate(vecpnstmp.begin(), vecpnstmp.end(), std::string("")); 
        
        // Return a formatted string with all the descriptors to 
        // print in the dock output
        //ostringstream strstrm;
        
        strstrm         << "\n" "-----------------------------------" << endl
                        << " Molecule:         \t" << mol.title << endl
                        << " Rotatable_Bonds:  \t" << mol.rot_bonds << endl
                        << " Heavy_Atoms:      \t" << mol.heavy_atoms << endl
                        << " Molecular_Weight: \t" << mol.mol_wt << endl
                        << " Formal_Charge:    \t" << mol.formal_charge << endl
                        << " HBond_donors:     \t" << mol.hb_donors << endl
                        << " HBond_acceptors:  \t" << mol.hb_acceptors << endl
                        << " RD_num_arom_rings:\t" << mol.num_arom_rings << endl
                        << " RD_num_alip_rings:\t" << mol.num_alip_rings << endl
                        << " RD_num_sat_rings: \t" << mol.num_sat_rings << endl
                        << " RD_Stereocenters: \t" << mol.num_stereocenters << endl
                        << " RD_Spiro_atoms:   \t" << mol.num_spiro_atoms << endl
                        << setprecision(2) << fixed
                        //<< " xLog_P:          \t" << mol.xlogp << endl
                        << " RD_LogP:          \t" << mol.clogp << endl
                        << " RD_TPSA:          \t" << mol.tpsa << endl
                        << " RD_SYNTHA:        \t" << mol.sa_score << endl
                        << " RD_QED:           \t" << mol.qed_score << endl
                        << " RD_LogS:          \t" << mol.esol << endl
                        << " RD_PAINS:         \t" << mol.pns << endl
                        << " RD_PAINS_Names:   \t" << molpns_name << endl 
                        << " RD_SMILES:        \t" << mol.smiles << endl
                        << " RD_MACCS:         \t" << mol.MACCS << endl 
                        << endl; 
        #else
               
        strstrm         << "\n" "-----------------------------------" << endl
                        << " Molecule:        \t" << mol.title << endl
                        << " Rotatable_Bonds: \t" << mol.rot_bonds << endl
                        << " Heavy_Atoms:     \t" << mol.heavy_atoms << endl
                        << " HBond_donors:    \t" << mol.hb_donors << endl
                        << " HBond_acceptors: \t" << mol.hb_acceptors << endl
                        << setprecision(2) << fixed
                        << " Molecular_Weight:\t" << mol.mol_wt << endl
                        << " Formal_Charge:   \t" << mol.formal_charge << endl
                        //<< " xLog_P:          \t" << mol.xlogp << endl
                        << endl;
        #endif

        return strstrm.str();
}
