#include "xlogp.h"
using namespace std;

//////////////////////////////////////////////////////////////////
//                                                              //
// XlogP and other molecular descriptors                        //
// Noel Carrascal, Sudipto Mukherjee, Robert C. Rizzo           //
// Edited and modified by Trent E Balius                        //
//     J. Chem. Inf. Comput. Sci. 1997, 37, 615-621             //
//////////////////////////////////////////////////////////////////

// There are 83 atom types in XLOGP
// Contribution of each atom type to the total xlogp value
const double XLogP_AtomConribution[] =  
{   9999.9, 0.484, 0.168, -0.181, 0.358,
0.009, -0.344, -0.439, 0.051, -0.138,
-0.417, -0.454, -0.378, 0.223, -0.598,
-0.396, -0.699, -0.362, 0.395, 0.236,
-0.166, 1.726, 0.098, -0.108, 1.637,
1.774, 0.281, 0.142, 0.715, 0.302,
-0.064, 0.079, 0.200, 0.869, 0.316,
0.054, 0.347, 0.046, -0.399, -0.029,
-0.330, 0.397, 0.068, 0.327, -2.057,
0.218, -0.582, -0.449, -0.774, 0.040,
-0.381, 0.443, -0.117, -2.052, -1.716,
0.321, -0.921, -0.704, 0.119, 1.192,
0.434, 0.587, 0.668, -0.791, -0.212,
0.016, 0.752, 1.071, 0.964, -1.817,
-1.214, -0.778, 0.493, 1.010, 1.187,
1.489, -0.802, -0.256, 1.626, 0.077,
0.264, -1.02
};

// Make empty constructor
xlogp::xlogp() 
{

}


// Populate number_of_H_Donors[] and number_of_H_Acceptors[]
// arrays in the Dockmol object
void xlogp::get_H_BOND_Don_Acc (DOCKMol & recept, string s)
{
        //cout << "xlogp::get_H_BOND_Don_Acc" << endl;
	recept.number_of_H_Donors = 0;
	recept.number_of_H_Acceptors = 0;
	
        for (int i = 0; i < recept.num_atoms; i++)
	{
                if ( isAnOxygen(recept.atom_types[i]))
		{
			// atom i = oxygen, we found a hbond acceptor 
                        recept.H_bond_acceptor[recept.number_of_H_Acceptors] = i;
			recept.number_of_H_Acceptors++;
			
			// check for any polar hydrogens bound this this oxygen
			for (int j = 0; j < recept.num_bonds; j++)
			{
				if (recept.bonds_origin_atom[j] == i)
				{
					if (recept.atom_types[recept.bonds_target_atom[j]] == "H" ||
						recept.atom_types[recept.bonds_target_atom[j]] == "H.spc"
						|| recept.atom_types[recept.bonds_target_atom[j]] ==
						"H.t3p")
					{
                        // found a polar hydrogen O-H
						recept.H_bond_donor[recept.number_of_H_Donors] =
							recept.bonds_target_atom[j];
						recept.number_of_H_Donors++;
					}
				}
				
				// Why are we doing this a second time??
				if (recept.bonds_target_atom[j] == i)
				{
					//    cout<<"Went inside second if of oxygen"<<endl;
					if (recept.atom_types[recept.bonds_origin_atom[j]] == "H" ||
						recept.atom_types[recept.bonds_origin_atom[j]] == "H.spc"
						|| recept.atom_types[recept.bonds_origin_atom[j]] ==
						"H.t3p")
					{
						recept.H_bond_donor[recept.number_of_H_Donors] =
							recept.bonds_origin_atom[j];
						recept.number_of_H_Donors++;
					}
				}
			}

		}
                else if ( isANitrogen(recept.atom_types[i]) )
		{
            // atom i = nitrogen, one more hbond acceptor
			recept.H_bond_acceptor[recept.number_of_H_Acceptors] = i;
			recept.number_of_H_Acceptors++;

			for (int j = 0; j < recept.num_bonds; j++)
			{
				if (recept.bonds_origin_atom[j] == i)
				{
                    // found a polar hydrogen N-H
					if (recept.atom_types[recept.bonds_target_atom[j]] == "H" ||
						recept.atom_types[recept.bonds_target_atom[j]] == "H.spc"
						|| recept.atom_types[recept.bonds_target_atom[j]] ==
						"H.t3p")
					{
						recept.H_bond_donor[recept.number_of_H_Donors] =
							recept.bonds_target_atom[j];
						recept.number_of_H_Donors++;
					}
				}
				if (recept.bonds_target_atom[j] == i)
				{
					//        cout<<"Went inside second if of nitrogen"<<endl;
					if (recept.atom_types[recept.bonds_origin_atom[j]] == "H" ||
						recept.atom_types[recept.bonds_origin_atom[j]] == "H.spc"
						|| recept.atom_types[recept.bonds_origin_atom[j]] ==
						"H.t3p")
					{
						recept.H_bond_donor[recept.number_of_H_Donors] =
							recept.bonds_origin_atom[j];
						recept.number_of_H_Donors++;
					}
				}
			}
		}
                else if ( isASulfur( recept.atom_types[i] ) )
		{
			// atom i = sulphur, lipinsky's rules does not acceptor sulphur as a hbond aceptor
			recept.H_bond_acceptor[recept.number_of_H_Donors] = i;
			recept.number_of_H_Acceptors++;
			for (int j = 0; j < recept.num_bonds; j++)
			{
				if (recept.bonds_origin_atom[j] == i)
				{
					if (recept.atom_types[recept.bonds_target_atom[j]] == "H" ||
						recept.atom_types[recept.bonds_target_atom[j]] == "H.spc"
						|| recept.atom_types[recept.bonds_target_atom[j]] ==
						"H.t3p")
					{
						recept.H_bond_donor[recept.number_of_H_Donors] =
							recept.bonds_target_atom[j];
						recept.number_of_H_Donors++;
					}
				}
				if (recept.bonds_target_atom[j] == i)
				{
					//     cout<<"Went inside second if of sulfure"<<endl;
					if (recept.atom_types[recept.bonds_origin_atom[j]] == "H" ||
						recept.atom_types[recept.bonds_origin_atom[j]] == "H.spc"
						|| recept.atom_types[recept.bonds_origin_atom[j]] ==
						"H.t3p")
					{
						recept.H_bond_donor[recept.number_of_H_Donors] =
							recept.bonds_origin_atom[j];
						recept.number_of_H_Donors++;
					}
				}
			}
		}
	}
/*	cout << "Results for " << s << " " << recept.title << endl;
	cout << s << "_#_atoms = " << recept.num_atoms << endl;
	cout << s << "_#_H-donors = " << recept.number_of_H_Donors << endl;
	cout << s << "_#_H-acceptors = " << recept.number_of_H_Acceptors << endl; */
}


void xlogp::getInterMolec_H_Bonds (DOCKMol & lig, DOCKMol & rec)
{
	int Hydrogen_Bonds = 0;
	double distance = 0;
	double xc, yc, zc;
	int k, m;
	// int                             *H_bond_donor;
	// int                             *H_bond_acceptor;
	for (k = 0; k < lig.number_of_H_Donors; k++)
	{
		for (m = 0; m < rec.number_of_H_Acceptors; m++)
		{
			xc = lig.x[lig.H_bond_donor[k]] - rec.x[rec.H_bond_acceptor[m]];
			yc = lig.y[lig.H_bond_donor[k]] - rec.y[rec.H_bond_acceptor[m]];
			zc = lig.z[lig.H_bond_donor[k]] - rec.z[rec.H_bond_acceptor[m]];
			distance = (xc * xc) + (yc * yc) + (zc * zc);
			//  cout<<"LigDon("<<lig.H_bond_donor[k]<<")="<<lig.atom_types[k]<<", RecAcc("<<rec.H_bond_acceptor[m]<<")="<<rec.atom_types[k]<<" Dist="<<distance<<endl;
			if (distance < 6.25)
			{
				//     cout<<"        LigDon("<<lig.H_bond_donor[k]<<")="<<lig.atom_types[k]<<", RecAcc("<<rec.H_bond_acceptor[m]<<")="<<rec.atom_types[k]<<" Dist="<<distance<<endl;
				Hydrogen_Bonds++;
			}
		}
	}
	for (k = 0; k < lig.number_of_H_Acceptors; k++)
	{
		for (m = 0; m < rec.number_of_H_Donors; m++)
		{
			xc = lig.x[lig.H_bond_acceptor[k]] - rec.x[rec.H_bond_donor[m]];
			yc = lig.y[lig.H_bond_acceptor[k]] - rec.y[rec.H_bond_donor[m]];
			zc = lig.z[lig.H_bond_acceptor[k]] - rec.z[rec.H_bond_donor[m]];
			distance = (xc * xc) + (yc * yc) + (zc * zc);
			//   cout<<"LigAcc("<<lig.H_bond_acceptor[k]<<")="<<lig.atom_types[k]<<", RecDon=("<<rec.H_bond_donor[m]<<")="<<rec.atom_types[k]<<" Dist="<<distance<<endl;
			if (distance < 6.25)
			{
				//       cout<<"        LigAcc("<<lig.H_bond_acceptor[k]<<")="<<lig.atom_types[k]<<", RecDon=("<<rec.H_bond_donor[m]<<")="<<rec.atom_types[k]<<" Dist="<<distance<<endl;
				Hydrogen_Bonds++;
			}
		}
	}
	cout << "#_intermolec_H-bonds = " << Hydrogen_Bonds << endl << endl;
}

bool xlogp::isONSP (string s)
{
	if (isAnOxygen (s))
		return true;

	else if (isANitrogen (s))
		return true;

	else if (isASulfur (s))
		return true;

	else if (s == "P.3")
		return true;

	else
		return false;
}

bool xlogp::isONSP_HAL (string s)
{
        if (isONSP (s))
                return true;
        if (isHalogen (s))
                return true;
	else
		return false;
}

bool xlogp::isAromatic (string s)
{
	if (s == "C.ar" || s == "N.ar")
		return true;

	else
		return false;
}

bool xlogp::isHalogen (string s)
{
	if (s == "F" || s == "Cl" || s == "Br" || s == "I" || s == "Hal"){
                //cout << "Halogen " << s << endl;
		return true;

	}else
		return false;
}

bool xlogp::isPiBonded (string s)
{
	if (s == "C.ar" || s == "C.2" || s == "C.1" || s == "N.ar" || s == "N.am"
		|| s == "N.pl3" || s == "N.2" || s == "N.1" || s == "O.2" || s == "S.2")
		return true;

	else
		return false;
}

double xlogp::getAtomMW (string s)
{
	if (isAnOxygen (s))
		return 15.9994;

	else if (isASulfur (s))
		return 32.066;

	else if (isANitrogen (s))
		return 14.0067;

	else if (isACarbon (s) || s == "C.ar" || s == "C.cat")
		return 12.011;

	else if (s == "F")
		return 18.998403;

	else if (s == "Cl")
		return 35.453;

	else if (s == "Br")
		return 79.904;

	else if (s == "I")
		return 126.9045;

	else if (s == "H" || s == "H.spc" || s == "H.t3p")
		return 1.00794;

	else if (s == "LP")  // Lone Pair
		return 0;

	else if (s == "Li")
		return 6.941;

	else if (s == "Na")
		return 22.98977;

	else if (s == "Mg")
		return 24.305;

	else if (s == "Al")
		return 26.98154;

	else if (s == "Si")
		return 28.0855;

	else if (s == "K")
		return 39.0983;

	else if (s == "Ca")
		return 40.078;

	else if (s == "Cr.th" || s == "Cr.oh")
		return 51.996;

	else if (s == "Mn")
		return 54.9380;

	else if (s == "Fe")
		return 55.847;

	else if (s == "Co.oh")
		return 58.9332;

	else if (s == "Cu")
		return 63.546;

	else if (s == "Zn")
		return 65.39;

	else if (s == "Se")
		return 78.96;

	else if (s == "Mo")
		return 95.94;

	else if (s == "Sn")
		return 118.710;
	else return 0;   // fail-safe: atom not found
}


int xlogp::nHdonnors (string s, string s2, int k)
{
	// Check if atom type 's' is a polar hydrogen
	// OH or NH
	int h = 0;
	if (s == "H" || s == "H.spc" || s == "H.t3p")
	{
		if (isANitrogen (s2) || isAnOxygen (s2))
		{
			h = 1;
		}
	}
	return h;
}


bool xlogp::isACarbon (string s)
{
	if (s == "C.1" || s == "C.2" || s == "C.3" || s == "C.ar" || s == "C.cat")
		return true;

	else
		return false;
}

bool xlogp::isAnOxygen (string s)
{
	if (s == "O.3" || s == "O.2" || s == "O.co2" || s == "O.spc" || s == "O.t3p")
		return true;

	else
		return false;
}

bool xlogp::isASulfur (string s)
{
	if (s == "S.3" || s == "S.2" || s == "S.o" || s == "S.o2" || s == "S.O" || s == "S.O2")
return true;

	else
		return false;
}

bool xlogp::isANitrogen (string s)
{
	if (s == "N.1" || s == "N.2" || s == "N.3" || s == "N.ar" || s == "N.am"
		|| s == "N.pl3" || s == "N.4")
		return true;

	else
		return false;
}


void xlogp::getXLogP (DOCKMol & molec)
{
        ///cout << "xlogp::getXLogP" << endl;
	int atype;
	int nghbCnt = 0;
	//string atomsTable[molec.num_atoms][molec.num_atoms+1];
        //vector < vector <string> >  atomsTable(molec.num_atoms, vector <string> (molec.num_atoms+1,"") );
        vector < vector <string> >  atomsTable;
        atomsTable.resize(molec.num_atoms);
        for (int i = 0; i < molec.num_atoms; i++)
             atomsTable[i].resize(molec.num_atoms+1);


	int atomsTableInt[molec.num_atoms][molec.num_atoms+1];

        // one atom might be bound to all other atoms eg CH4, or SF6.
        // by adding on addtional column these in shurs that there is
        // in all rows a nn or -1.

	double atomsDistance[molec.num_atoms][molec.num_atoms];

	int nBonds[molec.num_atoms];
	int nHydrogens[molec.num_atoms];
	//string bondTypes[molec.num_atoms][molec.num_bonds + 1];
        vector < vector <string> > bondTypes;
        bondTypes.resize(molec.num_atoms);
        for (int i = 0; i<molec.num_atoms; i++){
             bondTypes[i].resize(molec.num_bonds+1);
        }

        //cout << "num_atoms = " << molec.num_atoms << endl;
        //cout << "num_bonds = " << molec.num_bonds << endl;

	int Hdonor[molec.num_atoms];
	int HdonorNumberHydrogens[molec.num_atoms];
	int Hacceptor[molec.num_atoms];
	int charge[molec.num_atoms];
	int N4number[molec.num_atoms];

	int adjacencyMatrix[molec.num_atoms][molec.num_atoms];

	int minDistance[molec.num_atoms];
	int HbondPairs[molec.num_atoms][2];

	int bondsBetween[molec.num_atoms][molec.num_atoms];

	// Assign HBond Donors and Acceptors
	get_H_BOND_Don_Acc (molec, "Lig");

	for (int i = 0; i < molec.num_atoms; i++)
	{
                //cout << "start :: atomsTableInt[0][1] =  " << atomsTableInt[0][1] << endl;
		Hdonor[i] = -1;           // Table of HBond donors
		Hacceptor[i] = -1;        // Table of Hbond acceptors
		HdonorNumberHydrogens[i] = 0;     // # of polar hydrogens attached to atom[i]
		charge[i] = 0;            // Assume all atoms are neutral initially
		N4number[i] = -1;
		minDistance[i] = -1;
		HbondPairs[i][0] = -1;
		HbondPairs[i][1] = -1;
		for (int ii = 0; ii < molec.num_atoms+1; ii++)
		{
			adjacencyMatrix[i][ii] = -9;
			bondsBetween[i][ii] = -1;
			atomsTable[i][ii] = "nn";
			atomsTableInt[i][ii] = -1;
			if (ii < (molec.num_bonds + 1))
				bondTypes[i][ii] = "nn";
			if (i == 0)
			{
				nBonds[ii] = 0;
				nHydrogens[ii] = 0;
			}
			if (i == ii)
			{
				atomsDistance[i][ii] = 0;
			}
			else if (ii < molec.num_atoms )
			{
				atomsDistance[i][ii] =
					((molec.x[i] - molec.x[ii]) * (molec.x[i] - molec.x[ii])) +
					((molec.y[i] - molec.y[ii]) * (molec.y[i] - molec.y[ii])) +
					((molec.z[i] - molec.z[ii]) * (molec.z[i] - molec.z[ii]));
			}
		}
                //cout << "middle :: atomsTableInt[0][1] =  " << atomsTableInt[0][1] << endl;

		// Find out which atom is connected to which
		// atomsTable[i] has a list of atom types bonded to atom[i]
		for (int j = 0; j < molec.num_bonds; j++)
		{
			if (molec.bonds_origin_atom[j] == i)
			{
				atomsTable[i][0] = molec.atom_types[i];
				nBonds[i] = nBonds[i] + 1;
				atomsTable[i][nBonds[i]] =
					molec.atom_types[molec.bonds_target_atom[j]];
				atomsTableInt[i][nBonds[i]] = molec.bonds_target_atom[j];
                                //cout << " (1) atomsTableInt[" << i << "][" << nBonds[i] << "] = " <<atomsTableInt[i][nBonds[i]] << endl;
				if (molec.atom_types[molec.bonds_target_atom[j]] == "H"
					|| molec.atom_types[molec.bonds_target_atom[j]] == "H.spc"
					|| molec.atom_types[molec.bonds_target_atom[j]] == "H.t3p")
				{
					nHydrogens[i] = nHydrogens[i] + 1;
				}
				bondTypes[i][nBonds[i]] = molec.bond_types[j];
                                //cout << " (1) bondTypes[" << i << "][" << nBonds[i]<< "] = " << bondTypes[i][nBonds[i]] << endl;
			}
			if (molec.bonds_target_atom[j] == i)
			{
				atomsTable[i][0] = molec.atom_types[i];
				nBonds[i] = nBonds[i] + 1;
				atomsTable[i][nBonds[i]] =
					molec.atom_types[molec.bonds_origin_atom[j]];
				atomsTableInt[i][nBonds[i]] = molec.bonds_origin_atom[j];
                                //cout << " (2) atomsTableInt[" << i << "][" << nBonds[i] << "] = " <<atomsTableInt[i][nBonds[i]] << endl;
				if (molec.atom_types[molec.bonds_origin_atom[j]] == "H"
					|| molec.atom_types[molec.bonds_origin_atom[j]] == "H.spc"
					|| molec.atom_types[molec.bonds_origin_atom[j]] == "H.t3p")
				{
					nHydrogens[i] = nHydrogens[i] + 1;
				}               
                                // Bug fixed: substituted target in last
				// two conditionals Not tested yet.
				bondTypes[i][nBonds[i]] = molec.bond_types[j];
                                //cout << " (2) bondTypes[" << i << "][" << nBonds[i]<< "] = " << bondTypes[i][nBonds[i]] << endl;
			}
		}
                
                //cout << "end::atomsTableInt[0][1] =  " << atomsTableInt[0][1] << endl;
	}

/*
        cout << "MATRIX" << endl;
        for (int i = 0; i<  molec.num_atoms; i++){
           for (int ii = 0; ii<  molec.num_atoms+1; ii++){

                 cout << atomsTable[i][ii] << ","; 
           }
           cout << endl;
        }
        cout << "END MATRIX" << endl;
*/

        // do we need to do this? trent
	//check and correct hybridization
	// End result: Calculate array charge[i]
	double numberBonds;
	int n = 1;
	double bondType = 0;
	int surroundingAromaticCarbons;
	for (int m = 0; m < molec.num_atoms; m++)
	{
		numberBonds = 0;
		if (isACarbon (atomsTable[m][0]))
		{
			n = 1;
			surroundingAromaticCarbons = 0;

			// atomsTable is initialized to 'nn' 
			while (atomsTable[m][n] != "nn" && n < (molec.num_bonds + 1))
			{
				if (bondTypes[m][n] == "3")
				{
					bondType = 3;
				}

				else if (bondTypes[m][n] == "2")
				{
					bondType = 2;
				}

				else if (bondTypes[m][n] == "1")
				{
					bondType = 1;
				}

				else if (bondTypes[m][n] == "ar")
				{
					surroundingAromaticCarbons++;
					if (surroundingAromaticCarbons == 1)
					{
						bondType = 1.5;
					}
					else if (surroundingAromaticCarbons == 2)
					{
						bondType = 1.5;
					}
					else if (surroundingAromaticCarbons == 3)
					{
						bondType = 1;
					}
				}
				//else if (bondTypes[m][n] == "am")
				//{
				//	bondType = 1;
				//}
                                else{
                                        cout << "alternative bondtype:" << bondTypes[m][n] << endl;
					bondType = 1;
                                }

				numberBonds = numberBonds + bondType;
				n++;
			}
			if (numberBonds == 3)// i replaced the inequality  with an equality, trent 2011.02.18
			{
				charge[m] = -1;
			}
			else if (numberBonds == 4)
			{
				// So this carbon is neutral
				charge[m] = 0;
			}
			else if (numberBonds == 5)
			{
				charge[m] = 1;
			}
			else // cech all
			{
                                cout << "Warning: Molecule has a bad Carbon: " << m+1 << endl;
                                cout << "numberBonds = " << numberBonds << endl; 
				charge[m] = 0;
                                //exit(1);
			}
		}
		else if (isANitrogen (atomsTable[m][0]))
		{
			n = 1;
			numberBonds = 0;
			while (atomsTable[m][n] != "nn" && n < (molec.num_bonds+1))
			{
				if (bondTypes[m][n] == "3")
				{
					bondType = 3;
				}

				else if (bondTypes[m][n] == "2")
				{
					bondType = 2;
				}

				else if (bondTypes[m][n] == "1")
				{
					bondType = 1;
				}

				else if (bondTypes[m][n] == "ar")
				{
					bondType = 1.5;
				}
				//else if (bondTypes[m][n] == "am")
				//{
				//	bondType = 1;
				//}
                                else{
                                        cout << "alternative bondtype:" << bondTypes[m][n] << endl;
					bondType = 1;
                                }
				numberBonds = numberBonds + bondType;
				n++;
			}
			if (numberBonds == 2)
			{
				charge[m] = -1;
			}
			else if (numberBonds == 3)
			{
				// This nitrogen is neutral
				charge[m] = 0;
			}
			else if (numberBonds == 4)
			{
				charge[m] = 1;
			}
                        else { // 
                                cout << "Warning: Molecule has a bad Nitrogen :" << m+1 << endl;
                                cout << "numberBonds = " << numberBonds << endl; 
				charge[m] = 0;
                                //exit(1);
                        }
		}
		else if (isAnOxygen (atomsTable[m][0]))
		{
			n = 1;
			numberBonds = 0;      //
			surroundingAromaticCarbons = 0;

                        //cout << "###  "<< atomsTable[m][0];
                        //cout << " and "<< atomsTable[m][n] << endl;

			while (atomsTable[m][n] != "nn" && n < (molec.num_bonds+1))
			{
				if (bondTypes[m][n] == "3")
				{
					bondType = 3;
				}

				else if (bondTypes[m][n] == "2")
				{
					bondType = 2;
				}

				else if (bondTypes[m][n] == "1")
				{
					bondType = 1;
				}
                                else if (bondTypes[m][n] == "ar")
                                {
                                        surroundingAromaticCarbons++;
                                        if (surroundingAromaticCarbons == 1)
                                            bondType = 2;
                                        else if (surroundingAromaticCarbons == 2)
                                            bondType = 1;
                                        else {
                                             cout << "more than 2 ar bonds with bad Oxygen :" << m+1 << endl;
                                             //exit (1);
                                        }
                                }
                                else{
                                        cout << "alternative bondtype:" << bondTypes[m][n] << endl;
					bondType = 1;
                                }
				numberBonds = numberBonds + bondType;
				n++;
			}
			if (numberBonds == 1)
			{
				charge[m] = -1;
			}
			else if (numberBonds == 2)
			{
				// This oxygen is neutral
				charge[m] = 0;
			}
			else if (numberBonds == 3)
			{
				charge[m] = 1;
			}
			else
			{
                                cout << "Warning: Molecule has a bad Oxygen :" << m+1 << endl;
                                cout << "numberBonds = " << numberBonds << endl; 
				charge[m] = 0;
                                //exit(1);
			}
		}
                else if ( isASulfur (atomsTable[m][0]))
                {
                        n = 1;
                        numberBonds = 0;      //
                        while (atomsTable[m][n] != "nn" && n < (molec.num_bonds+1))
                        {
                                if (bondTypes[m][n] == "3")
                                {
                                        bondType = 3;
                                }

                                else if (bondTypes[m][n] == "2")
                                {
                                        bondType = 2;
                                }

                                else if (bondTypes[m][n] == "1")
                                {
                                        bondType = 1;
                                }

                                else if (bondTypes[m][n] == "ar")
                                {
                                        surroundingAromaticCarbons++;
                                        if (surroundingAromaticCarbons == 1)
                                            bondType = 2;
                                        else if (surroundingAromaticCarbons == 2)
                                            bondType = 1;
                                        else {
                                             cout << "more than 2 ar bonds with bad Sulfur :" << m+1 << endl;
                                             //exit (1);
                                        }
                                }
                                else{
                                        cout << "alternative bondtype:" << bondTypes[m][n] << endl;
					bondType = 1;
                                }

                                numberBonds = numberBonds + bondType;
                                n++;
                        }
                        if (numberBonds == 1)
                        {
                                charge[m] = -1;
                        }
                        else if (numberBonds == 2)
                        {
                                // This sulfur is neutral
				charge[m] = 0;
                        }
                        else if (numberBonds == 3)
                        {
                                 charge[m] = 1;
                        }
                        else if (numberBonds == 4)
                        {
                                // This sulfur is neutral
				charge[m] = 0;
                        }
                        else if (numberBonds == 6)
                        {
                                // This sulfur is neutral
				charge[m] = 0;
                        }
                        // note the sulfur also may have oxidation of +8
                        // we deliberately do not take this into acount.
                        else
                        {
                                cout << "Error: Molecule has a bad sulfur :" << m+1 << endl;
                                cout << "numberBonds = " << numberBonds << endl; 
				charge[m] = 0;
                                //exit(1);
                        }
                }

	}

	// We do not correct protonation states for sulphur because
	// sulphur is too complicated



	// Now that we have the connectivity and protonation states corrected
	// we start going over every atom in the molecule and assign atom
	// types according to xlogp rules. xlogp requires that all atoms are 
	// neutral
	int temp = 1;
	int temp2 = 0;
	float XlogP = 0;
	float MWT = 0;
	bool isCNgroup = false;
	bool isNO2group = false;
	bool isN4inAminoAcid = false;
	bool checkIfAminoAcid = true;
	int FNumberGeminal = 0;
	int XNumberGeminal = 0;
	bool isAliphaticORAromatic = true;
	int numberHdonnors = 0;
	int numberHacceptors = 0;
	int numOfAliphaticAromaticContributions = 0;
	int N4Count = 0;

	// Calculate xlogP contribution of each atom
	// in the molecule
	for (int k = 0; k < molec.num_atoms; k++)
	{
		temp = 1;                 // Current atom in molec
		temp2 = 0;                // Heavy atom counter for nonH[]
		vector <string> nonH(10);           // Heavy atom type
		int nonH_INT[10];          // Corresponding Heavy atom number
		
		// We are assuming that no atom has more than four 
		// heavy atoms connected to it : Inceasing array size to 10 for testing

		// Make a list of all heavy atoms bonded to
		// current atom[k]
		while (temp < molec.num_atoms)
		{
			// 'nn' means we have reached the end of the atomstable list 
                        // for atom[k] 
                        // check if we reached the end of list atomsTable[k][temp]
			if (atomsTable[k][temp] != "nn")     
			{
				if (atomsTable[k][temp] != "H")  // this is a heavy atom
				{
					nonH[temp2] = atomsTable[k][temp];
					nonH_INT[temp2] = atomsTableInt[k][temp];
					temp2++;
				}
			}
			temp++;
		}                       //end while 

		// Calculate molecular weight of whole molecule
		MWT = MWT + getAtomMW (atomsTable[k][0]);

		// Assign xlogp contribution for carbon
		if (isACarbon (atomsTable[k][0]))
		{
			if (charge[k] == 0)
			{
				// If it is an sp3 carbon
				if (atomsTable[k][0] == "C.3")
				{
					if (nHydrogens[k] == 4)
					{
						// Our molecule is Methane
						atomsTableInt[k][0] = 1;
						XlogP = XlogP + XLogP_AtomConribution[1];
						numOfAliphaticAromaticContributions++;
					}
					else if (nHydrogens[k] == 3)  // We have CH3R or CH3X 
					{
						if (nonH[0] == "C.ar" || nonH[0] == "C.2"
							|| nonH[0] == "C.1")
						{
							// CH3R (pi!=0) in xlogp paper line #2
							atomsTableInt[k][0] = 2;
							XlogP = XlogP + XLogP_AtomConribution[2];
							numOfAliphaticAromaticContributions++;
						}
						else if (isONSP_HAL (nonH[0]))
						{
							// CH3X
							// isONSP_HAL = oxygen, nitrogen, sulphur phosphorus and halogens
							atomsTableInt[k][0] = 3;
							XlogP = XlogP + XLogP_AtomConribution[3];
						}
						else if (nonH[0] == "C.3" || nonH[0] == "C.cat")
						{
							// CH3R,  R is non aromatic (aliphatic)
							atomsTableInt[k][0] = 1;
							XlogP = XlogP + XLogP_AtomConribution[1];
							numOfAliphaticAromaticContributions++;
						}
					}
					else if (nHydrogens[k] == 2)
					{
						if ((isACarbon (nonH[0]) && !isACarbon (nonH[1]))
							|| (!isACarbon (nonH[0]) && isACarbon (nonH[1])))
						{
							atomsTableInt[k][0] = 6;
							XlogP = XlogP + XLogP_AtomConribution[6];
							numOfAliphaticAromaticContributions++;
						}
						else if (!isACarbon (nonH[0]) && !isACarbon (nonH[1]))
						{
							atomsTableInt[k][0] = 7;
							XlogP = XlogP + XLogP_AtomConribution[7];
							numOfAliphaticAromaticContributions++;
						}
						else if ((isACarbon (nonH[0])) && (isACarbon (nonH[1])))
						{
							if ((nonH[0] == "C.ar") || (nonH[1] == "C.ar") ||
								(nonH[0] == "C.2") || (nonH[1] == "C.2") ||
								(nonH[0] == "C.1") || (nonH[1] == "C.1"))
							{
								atomsTableInt[k][0] = 5;
								XlogP = XlogP + XLogP_AtomConribution[5];
								numOfAliphaticAromaticContributions++;
							}
							else
							{
								atomsTableInt[k][0] = 4;  // carbon is 1 unless an aromatic is found then it is 2
								XlogP = XlogP + XLogP_AtomConribution[4];
								numOfAliphaticAromaticContributions++;
							}
						}
						if (nonH[0] == "F" || nonH[1] == "F")
						{
							int temp = 0;
							if (isHalogen (nonH[0]))
								temp++;
							if (isHalogen (nonH[1]))
								temp++;
							if (temp == 1)
							{
							}
							else if (temp == 2)
							{
								FNumberGeminal = FNumberGeminal + (temp - 1);
							}
						}
						else
						{
							int temp = 0;
							if (isHalogen (nonH[0]))
								temp++;
							if (isHalogen (nonH[1]))
								temp++;
							if (temp == 1)
							{
							}
							else if (temp == 2)
							{
								XNumberGeminal = XNumberGeminal + (temp - 1);
							}
						}
					}
					else if (nHydrogens[k] == 1)
					{
						if (isACarbon (nonH[0]) && isACarbon (nonH[1])
							&& isACarbon (nonH[2]))
						{
							if ((nonH[0] == "C.ar" || nonH[1] == "C.ar"
								|| nonH[2] == "C.ar") || (nonH[0] == "C.2"
								|| nonH[1] == "C.2"
								|| nonH[2] == "C.2")
								|| (nonH[0] == "C.1" || nonH[1] == "C.1"
								|| nonH[2] == "C.1"))
							{
								atomsTableInt[k][0] = 9;
								XlogP = XlogP + XLogP_AtomConribution[9];
								numOfAliphaticAromaticContributions++;
							}
							else
							{
								atomsTableInt[k][0] = 8;
								XlogP = XlogP + XLogP_AtomConribution[8];
								numOfAliphaticAromaticContributions++;
							}
						}
						else
							if ((isACarbon (nonH[0]) && isONSP_HAL (nonH[1])
								&& isONSP_HAL (nonH[2])) || (isONSP_HAL (nonH[0])
								&&
								isACarbon (nonH[1])
								&&
								isONSP_HAL (nonH[2]))
								|| (isONSP_HAL (nonH[0]) && isONSP_HAL (nonH[1])
								&& isACarbon (nonH[2])))
							{
								atomsTableInt[k][0] = 11;
								XlogP = XlogP + XLogP_AtomConribution[11];
							}
							else if (isONSP_HAL (nonH[0]) && isONSP_HAL (nonH[1])
								&& isONSP_HAL (nonH[2]))
							{
								atomsTableInt[k][0] = 11;
								XlogP = XlogP + XLogP_AtomConribution[11];
							}
							else
								if ((isACarbon (nonH[0]) && isACarbon (nonH[1])
									&& isONSP_HAL (nonH[2])) || (isACarbon (nonH[0])
									&&
									isONSP_HAL (nonH[1])
									&&
									isACarbon (nonH[2]))
									|| (isONSP_HAL (nonH[0]) && isACarbon (nonH[1])
									&& isACarbon (nonH[2])))
								{
									atomsTableInt[k][0] = 10;
									XlogP = XlogP + XLogP_AtomConribution[10];
								}
								if (nonH[0] == "F" || nonH[1] == "F" || nonH[2] == "F")
								{
									int temp = 0;
									if (isHalogen (nonH[0]))
										temp++;
									if (isHalogen (nonH[1]))
										temp++;
									if (isHalogen (nonH[2]))
										temp++;
									if (temp == 1)
									{
									}
									else if (temp == 2)
									{
										FNumberGeminal = FNumberGeminal + (temp - 1);
									}
									else if (temp == 3)
									{
										FNumberGeminal = FNumberGeminal + temp;
									}
								}
								else
								{
									int temp = 0;
									if (isHalogen (nonH[0]))
										temp++;
									if (isHalogen (nonH[1]))
										temp++;
									if (isHalogen (nonH[2]))
										temp++;
									if (temp == 1)
									{
									}
									else if (temp == 2)
									{
										XNumberGeminal = XNumberGeminal + (temp - 1);
									}
									else if (temp == 3)
									{
										XNumberGeminal = XNumberGeminal + temp;
									}
								}
					}
					else if (nHydrogens[k] == 0)
					{
						if ((!isACarbon (nonH[0]) && isACarbon (nonH[1])
							&& isACarbon (nonH[2]) && isACarbon (nonH[3]))
							|| (isACarbon (nonH[0]) && !isACarbon (nonH[1])
							&& isACarbon (nonH[2]) && isACarbon (nonH[3]))
							|| (isACarbon (nonH[0]) && isACarbon (nonH[1])
							&& !isACarbon (nonH[2]) && isACarbon (nonH[3]))
							|| (isACarbon (nonH[0]) && isACarbon (nonH[1])
							&& isACarbon (nonH[2]) && !isACarbon (nonH[3])))
						{
							atomsTableInt[k][0] = 14;
							XlogP = XlogP + XLogP_AtomConribution[14];
						}
						else
							if ((isACarbon (nonH[0]) && isACarbon (nonH[1])
								&& !isACarbon (nonH[2]) && !isACarbon (nonH[3]))
								|| (isACarbon (nonH[0]) && !isACarbon (nonH[1])
								&& isACarbon (nonH[2]) && !isACarbon (nonH[3]))
								|| (isACarbon (nonH[0]) && !isACarbon (nonH[1])
								&& !isACarbon (nonH[2]) && isACarbon (nonH[3]))
								|| (!isACarbon (nonH[0]) && isACarbon (nonH[1])
								&& isACarbon (nonH[2]) && !isACarbon (nonH[3]))
								|| (!isACarbon (nonH[0]) && isACarbon (nonH[1])
								&& !isACarbon (nonH[2]) && isACarbon (nonH[3]))
								|| (!isACarbon (nonH[0]) && !isACarbon (nonH[1])
								&& isACarbon (nonH[2]) && isACarbon (nonH[3])))
							{
								atomsTableInt[k][0] = 15;
								XlogP = XlogP + XLogP_AtomConribution[15];
							}
							else
								if ((isACarbon (nonH[0]) && !isACarbon (nonH[1])
									&& !isACarbon (nonH[2]) && !isACarbon (nonH[3]))
									|| (!isACarbon (nonH[0]) && isACarbon (nonH[1])
									&& !isACarbon (nonH[2])
									&& !isACarbon (nonH[3]))
									|| (!isACarbon (nonH[0]) && !isACarbon (nonH[1])
									&& isACarbon (nonH[2]) && !isACarbon (nonH[3]))
									|| (!isACarbon (nonH[0]) && !isACarbon (nonH[1])
									&& !isACarbon (nonH[2]) && isACarbon (nonH[3])))
								{
									atomsTableInt[k][0] = 16;
									XlogP = XlogP + XLogP_AtomConribution[16];
								}
								else
									if ((!isACarbon (nonH[0]) && !isACarbon (nonH[1])
										&& !isACarbon (nonH[2]) && !isACarbon (nonH[3])))
									{
										atomsTableInt[k][0] = 17;
										XlogP = XlogP + XLogP_AtomConribution[17];
									}
									else
									{

										//After checking if there a Halogen(X) around the carbon, we check if there are 2 or more
										//aromatic carbons surrounding this carbon. If there are 2 or more, then this carbon is part of
										//a conjugated system. if there are not two aromatic carbons at least, then this carbon is not
										//part of a conjugated system. This does not check for conjugated systems that are not ring,
										//such as in amino acids where there are three carbons and the carbons in the ends do not meet
										//the criteria for this part of the code.  ??? I forgot What I meant by this.
										int numberOfAromaticNeighbors = 0;
										if ((nonH[0] == "N.ar") || (nonH[0] == "C.ar"))
											numberOfAromaticNeighbors++;
										if ((nonH[1] == "N.ar") || (nonH[1] == "C.ar"))
											numberOfAromaticNeighbors++;
										if ((nonH[2] == "N.ar") || (nonH[2] == "C.ar"))
											numberOfAromaticNeighbors++;
										if ((nonH[3] == "N.ar") || (nonH[3] == "C.ar"))
											numberOfAromaticNeighbors++;
										if (numberOfAromaticNeighbors < 2)
										{
											atomsTableInt[k][0] = 12; // carbon is 1 unless an aromatic is found then it is 2
											XlogP = XlogP + XLogP_AtomConribution[12];
											numOfAliphaticAromaticContributions++;
										}
										else
										{
											atomsTableInt[k][0] = 13; // carbon is 1 unless an aromatic is found then it is 2
											XlogP = XlogP + XLogP_AtomConribution[13];
											numOfAliphaticAromaticContributions++;
										}
									}
									if (nonH[0] == "F" || nonH[1] == "F" || nonH[2] == "F"
										|| nonH[3] == "F")
									{
										int temp = 0;
										if (isHalogen (nonH[0]))
											temp++;
										if (isHalogen (nonH[1]))
											temp++;
										if (isHalogen (nonH[2]))
											temp++;
										if (isHalogen (nonH[3]))
											temp++;
										if (temp == 1)
										{
										}
										else if (temp == 2)
										{
											FNumberGeminal = FNumberGeminal + (temp - 1);
										}
										else if (temp == 3)
										{
											FNumberGeminal = FNumberGeminal + temp;
										}
										else if (temp == 4)
										{
											FNumberGeminal = FNumberGeminal + (temp + 2);     //Number of pair when a carbon has 4 halogens
										}
									}
									else
									{
										int temp = 0;
										if (isHalogen (nonH[0]))
											temp++;
										if (isHalogen (nonH[1]))
											temp++;
										if (isHalogen (nonH[2]))
											temp++;
										if (isHalogen (nonH[3]))
											temp++;
										if (temp == 1)
										{
										}
										else if (temp == 2)
										{
											XNumberGeminal = XNumberGeminal + (temp - 1);
										}
										else if (temp == 3)
										{
											XNumberGeminal = XNumberGeminal + temp;
										}
										else if (temp == 4)
										{
											XNumberGeminal = XNumberGeminal + (temp + 2);     //Number of pair when a carbon has 4 halogens
										}
									}
					}
					else
					{
					}
				}
				else if (atomsTable[k][0] == "C.2")
				{
					if (nHydrogens[k] == 2)
					{
						atomsTableInt[k][0] = 18;
						XlogP = XlogP + XLogP_AtomConribution[18];
						numOfAliphaticAromaticContributions++;
					}
					else if (nHydrogens[k] == 1)
					{
						if (isACarbon (nonH[0]) && isACarbon (nonH[1]))
						{
							atomsTableInt[k][0] = 19;
							XlogP = XlogP + XLogP_AtomConribution[19];
							numOfAliphaticAromaticContributions++;
						}
						else if (!isACarbon (nonH[0]) && isACarbon (nonH[1]))
						{
							if (atomsTable[k][1] == "H")
							{
								if (isONSP_HAL (atomsTable[k][2])
									&& bondTypes[k][2] == "2")
								{       //bondtpes is off by 1 compared to nonH
									atomsTableInt[k][0] = 21;
									XlogP = XlogP + XLogP_AtomConribution[21];
								}
								else if (isONSP_HAL (atomsTable[k][3])
									&& bondTypes[k][3] == "2")
								{
									atomsTableInt[k][0] = 21;
									XlogP = XlogP + XLogP_AtomConribution[21];
								}
								else
								{
									atomsTableInt[k][0] = 20;
									XlogP = XlogP + XLogP_AtomConribution[20];
								}
							}
							else if (atomsTable[k][2] == "H")
							{
								if (isONSP_HAL (atomsTable[k][1])
									&& bondTypes[k][1] == "2")
								{       //bondtpes is off by 1 compared to nonH
									atomsTableInt[k][0] = 21;
									XlogP = XlogP + XLogP_AtomConribution[21];
								}
								else if (isONSP_HAL (atomsTable[k][3])
									&& bondTypes[k][3] == "2")
								{
									atomsTableInt[k][0] = 21;
									XlogP = XlogP + XLogP_AtomConribution[21];
								}
								else
								{
									atomsTableInt[k][0] = 20;
									XlogP = XlogP + XLogP_AtomConribution[20];
								}
							}
							else if (atomsTable[k][3] == "H")
							{
								if (isONSP_HAL (atomsTable[k][1])
									&& bondTypes[k][1] == "2")
								{       //bondtpes is off by 1 compared to nonH
									atomsTableInt[k][0] = 21;
									XlogP = XlogP + XLogP_AtomConribution[21];
								}
								else if (isONSP_HAL (atomsTable[k][2])
									&& bondTypes[k][2] == "2")
								{
									atomsTableInt[k][0] = 21;
									XlogP = XlogP + XLogP_AtomConribution[21];
								}
								else
								{
									atomsTableInt[k][0] = 20;
									XlogP = XlogP + XLogP_AtomConribution[20];
								}
							}
						}
						else if (isACarbon (nonH[0]) && !isACarbon (nonH[1]))
						{
							if (atomsTable[k][1] == "H")
							{
								if (isONSP_HAL (atomsTable[k][2])
									&& bondTypes[k][2] == "2")
								{       //bondtpes is off by 1 compared to nonH
									atomsTableInt[k][0] = 21;
									XlogP = XlogP + XLogP_AtomConribution[21];
								}
								else if (isONSP_HAL (atomsTable[k][3])
									&& bondTypes[k][3] == "2")
								{
									atomsTableInt[k][0] = 21;
									XlogP = XlogP + XLogP_AtomConribution[21];
								}
								else
								{
									atomsTableInt[k][0] = 20;
									XlogP = XlogP + XLogP_AtomConribution[20];
								}
							}
							else if (atomsTable[k][2] == "H")
							{
								if (isONSP_HAL (atomsTable[k][1])
									&& bondTypes[k][1] == "2")
								{       //bondtpes is off by 1 compared to nonH
									atomsTableInt[k][0] = 21;
									XlogP = XlogP + XLogP_AtomConribution[21];
								}
								else if (isONSP_HAL (atomsTable[k][3])
									&& bondTypes[k][3] == "2")
								{
									atomsTableInt[k][0] = 21;
									XlogP = XlogP + XLogP_AtomConribution[21];
								}
								else
								{
									atomsTableInt[k][0] = 20;
									XlogP = XlogP + XLogP_AtomConribution[20];
								}
							}
							else if (atomsTable[k][3] == "H")
							{
								if (isONSP_HAL (atomsTable[k][1])
									&& bondTypes[k][1] == "2")
								{       //bondtpes is off by 1 compared to nonH
									atomsTableInt[k][0] = 21;
									XlogP = XlogP + XLogP_AtomConribution[21];
								}
								else if (isONSP_HAL (atomsTable[k][2])
									&& bondTypes[k][2] == "2")
								{
									atomsTableInt[k][0] = 21;
									XlogP = XlogP + XLogP_AtomConribution[21];
								}
								else
								{
									atomsTableInt[k][0] = 20;
									XlogP = XlogP + XLogP_AtomConribution[20];
								}
							}
						}
						else if (!isACarbon (nonH[0]) && !isACarbon (nonH[1]))
						{
							atomsTableInt[k][0] = 21;
							XlogP = XlogP + XLogP_AtomConribution[21];
						}
					}
					else if (nHydrogens[k] == 0)
					{
						if (isACarbon (nonH[0]) && isACarbon (nonH[1])
							&& isACarbon (nonH[2]))
						{       //All Carbons
							atomsTableInt[k][0] = 22;
							XlogP = XlogP + XLogP_AtomConribution[22];
							numOfAliphaticAromaticContributions++;
						}
						else
							if ((isACarbon (nonH[0]) && isONSP_HAL (nonH[1])
								&& isONSP_HAL (nonH[2])) || (isONSP_HAL (nonH[0])
								&&
								isACarbon (nonH[1])
								&&
								isONSP_HAL (nonH[2]))
								|| (isONSP_HAL (nonH[0]) && isONSP_HAL (nonH[1])
								&& isACarbon (nonH[2])))
							{
								if (isACarbon (nonH[0]))
								{   // nonH and bondTypes are off by one carefule
									if (bondTypes[k][2] == "2"
										|| bondTypes[k][3] == "2")
									{
										atomsTableInt[k][0] = 24;
										XlogP = XlogP + XLogP_AtomConribution[24];
									}
									else
									{
										atomsTableInt[k][0] = 23;
										XlogP = XlogP + XLogP_AtomConribution[23];
									}
								}
								else if (isACarbon (nonH[1]))
								{
									if (bondTypes[k][1] == "2"
										|| bondTypes[k][3] == "2")
									{
										atomsTableInt[k][0] = 24;
										XlogP = XlogP + XLogP_AtomConribution[24];
									}
									else
									{
										atomsTableInt[k][0] = 23;
										XlogP = XlogP + XLogP_AtomConribution[23];
									}
								}
								else if (isACarbon (nonH[2]))
								{
									if (bondTypes[k][1] == "2"
										|| bondTypes[k][2] == "2")
									{
										atomsTableInt[k][0] = 24;
										XlogP = XlogP + XLogP_AtomConribution[24];
									}
									else
									{
										atomsTableInt[k][0] = 23;
										XlogP = XlogP + XLogP_AtomConribution[23];
									}
								}
							}
							else
								if ((isACarbon (nonH[0]) && isACarbon (nonH[1])
									&& isONSP_HAL (nonH[2])) || (isACarbon (nonH[0])
									&&
									isONSP_HAL (nonH[1])
									&&
									isACarbon (nonH[2]))
									|| (isONSP_HAL (nonH[0]) && isACarbon (nonH[1])
									&& isACarbon (nonH[2])))
								{
									if (isONSP_HAL (nonH[2]))
									{   // nonH and bondTypes are off by one carefule
										if (bondTypes[k][1] == "2"
											|| bondTypes[k][2] == "2")
										{
											atomsTableInt[k][0] = 23;
											XlogP = XlogP + XLogP_AtomConribution[23];
										}
										else
										{
											atomsTableInt[k][0] = 24;
											XlogP = XlogP + XLogP_AtomConribution[24];
										}
									}
									else if (isONSP_HAL (nonH[1]))
									{
										if (bondTypes[k][1] == "2"
											|| bondTypes[k][3] == "2")
										{
											atomsTableInt[k][0] = 23;
											XlogP = XlogP + XLogP_AtomConribution[23];
										}
										else
										{
											atomsTableInt[k][0] = 24;
											XlogP = XlogP + XLogP_AtomConribution[24];
										}
									}
									else if (isONSP_HAL (nonH[0]))
									{
										if (bondTypes[k][2] == "2"
											|| bondTypes[k][3] == "2")
										{
											atomsTableInt[k][0] = 23;
											XlogP = XlogP + XLogP_AtomConribution[23];
										}
										else
										{
											atomsTableInt[k][0] = 24;
											XlogP = XlogP + XLogP_AtomConribution[24];
										}
									}
								}
								else if (isONSP_HAL (nonH[0]) && isONSP_HAL (nonH[1])
									&& isONSP_HAL (nonH[2]))
								{
									atomsTableInt[k][0] = 25;
									XlogP = XlogP + XLogP_AtomConribution[25];
								}
								if (nonH[0] == "F" || nonH[1] == "F" || nonH[2] == "F")
								{
									int temp = 0;
									if (isHalogen (nonH[0]))
										temp++;
									if (isHalogen (nonH[1]))
										temp++;
									if (temp == 1)
									{
									}
									else if (temp == 2)
									{
										FNumberGeminal = FNumberGeminal + (temp - 1);
									}
								}
								else
								{
									int temp = 0;
									if (isHalogen (nonH[0]))
										temp++;
									if (isHalogen (nonH[1]))
										temp++;
									if (isHalogen (nonH[2]))
										temp++;
									if (temp == 1)
									{
									}
									else if (temp == 2)
									{
										XNumberGeminal = XNumberGeminal + (temp - 1);
									}
									else if (temp == 3)
									{
										XNumberGeminal = XNumberGeminal;
									}
								}
					}
					else
					{
						atomsTableInt[k][0] = 18; // carbon is 1 unless an aromatic is found then it is 2
						XlogP = XlogP + XLogP_AtomConribution[18];
						numOfAliphaticAromaticContributions++;
					}
				}
				else if (atomsTable[k][0] == "C.ar")
				{
					if (nHydrogens[k] == 1)
					{
						if (isACarbon (nonH[0]) && isACarbon (nonH[1]))
						{
							atomsTableInt[k][0] = 26;
							XlogP = XlogP + XLogP_AtomConribution[26];
						}
						else
						{
							if (isACarbon (nonH[0]) || isACarbon (nonH[1]))
							{
								atomsTableInt[k][0] = 27;
								XlogP = XlogP + XLogP_AtomConribution[27];
							}
							else
							{
								atomsTableInt[k][0] = 28;
								XlogP = XlogP + XLogP_AtomConribution[28];
							}
						}
					}
					else
					{
						if (isACarbon (nonH[0]) && isACarbon (nonH[1])
							&& isACarbon (nonH[2]))
						{
							atomsTableInt[k][0] = 29;
							XlogP = XlogP + XLogP_AtomConribution[29];
						}
						else
						{
							if ((isACarbon (nonH[0]) && isACarbon (nonH[1])))
							{
								if (nonH[0] == "C.ar" && nonH[1] == "C.ar")
								{
									atomsTableInt[k][0] = 30;
									XlogP = XlogP + XLogP_AtomConribution[30];
								}
								else
								{
									atomsTableInt[k][0] = 31;
									XlogP = XlogP + XLogP_AtomConribution[31];
								}
							}
							else if (isACarbon (nonH[1]) && isACarbon (nonH[2]))
							{
								if (nonH[1] == "C.ar" && nonH[2] == "C.ar")
								{
									atomsTableInt[k][0] = 30;
									XlogP = XlogP + XLogP_AtomConribution[30];
								}
								else
								{
									atomsTableInt[k][0] = 31;
									XlogP = XlogP + XLogP_AtomConribution[31];
								}
							}
							else if ((isACarbon (nonH[0]) && isACarbon (nonH[2])))
							{
								if (nonH[0] == "C.ar" && nonH[2] == "C.ar")
								{
									atomsTableInt[k][0] = 30;
									XlogP = XlogP + XLogP_AtomConribution[30];
								}
								else
								{
									atomsTableInt[k][0] = 31;
									XlogP = XlogP + XLogP_AtomConribution[31];
								}
							}
							else
								if ((isACarbon (nonH[0]) && !isACarbon (nonH[1])
									&& !isACarbon (nonH[2]))
									|| (isACarbon (nonH[1]) && !isACarbon (nonH[0])
									&& !isACarbon (nonH[2]))
									|| (isACarbon (nonH[2]) && !isACarbon (nonH[1])
									&& !isACarbon (nonH[0])))
								{
									atomsTableInt[k][0] = 32;
									XlogP = XlogP + XLogP_AtomConribution[32];
								}
								else
								{
									if ((isONSP (nonH[0]) && isONSP (nonH[1]))
										|| (isONSP (nonH[0]) && isONSP (nonH[2]))
										|| (isONSP (nonH[1]) && isONSP (nonH[2])))
									{
										atomsTableInt[k][0] = 33;
										XlogP = XlogP + XLogP_AtomConribution[33];
									}
									else
									{
										atomsTableInt[k][0] = 34;     // Note: this also accounts for connected to two Aromatic carbons and an X not just one as in the paper.                                                                                                                
										XlogP = XlogP + XLogP_AtomConribution[34];
									}
								}
						}
					}
				}
				else if (atomsTable[k][0] == "C.1")
				{
					if (nHydrogens[k] == 1)
					{
						if (nonH[0] == "N.1" && !isCNgroup)
						{
							atomsTableInt[k][0] = 77;
							XlogP = XlogP + XLogP_AtomConribution[77];
							isCNgroup = true;
						}
						else if (isACarbon (nonH[0]))
						{
							atomsTableInt[k][0] = 35;
							XlogP = XlogP + XLogP_AtomConribution[35];
						}
					}
					else
					{
						if ((nonH[0] == "N.1" || nonH[1] == "N.1") && !isCNgroup)
						{
							atomsTableInt[k][0] = 77;
							XlogP = XlogP + XLogP_AtomConribution[77];
							isCNgroup = true;
						}
						else
						{
							atomsTableInt[k][0] = 36;
							XlogP = XlogP + XLogP_AtomConribution[36];
							if ((isANitrogen (nonH[0]) && isASulfur (nonH[1]))
								|| (isANitrogen (nonH[1]) && isASulfur (nonH[0])))
							{
								XlogP = XlogP + 5.148;

								//For the group -NCS, atom contributions are added base on atom types. A correction is done when
								// a SP 1 carbon together with a N and a S is found. 
							}
						}
					}
				}
				else if (atomsTable[k][0] == "C.cat")
				{
				}
			}
			else
			{                   //if(charge!=0)
			}
		}
		else if (isANitrogen (atomsTable[k][0]))
		{
			isAliphaticORAromatic = false;
			if (charge[k] == 0)
			{
                                //cout << "I AM HERE (0)" << endl;
				if (atomsTable[k][0] == "N.4")
				{
				}
				else if (atomsTable[k][0] == "N.3")
				{
					if (nHydrogens[k] == 2)
					{
						if (isACarbon (nonH[0]))
						{
							if (isPiBonded (nonH[0]))
							{
								atomsTableInt[k][0] = 47;
								XlogP = XlogP + XLogP_AtomConribution[47];
							}
							else
							{
								atomsTableInt[k][0] = 46;
								XlogP = XlogP + XLogP_AtomConribution[46];
							}
						}
						else
						{
							atomsTableInt[k][0] = 48;
							XlogP = XlogP + XLogP_AtomConribution[48];
						}
						Hacceptor[numberHacceptors] = k;
						//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
						//molec.number_of_H_Acceptors++;
						numberHacceptors++;
						Hdonor[numberHdonnors] = k;
						//molec.H_bond_donor[molec.number_of_H_Donors] = k;
						//molec.number_of_H_Donors++;
						HdonorNumberHydrogens[numberHdonnors]++;
						HdonorNumberHydrogens[numberHdonnors]++;
						numberHdonnors++;
					}
					else if (nHydrogens[k] == 1)
					{
						if (isACarbon (nonH[0]) && isACarbon (nonH[1]))
						{
							atomsTableInt[k][0] = 49;
							XlogP = XlogP + XLogP_AtomConribution[49];
						}
						else
						{
							atomsTableInt[k][0] = 50;
							XlogP = XlogP + XLogP_AtomConribution[50];
						}
						if (nonH[0] == "F" || nonH[1] == "F")
						{
							int temp = 0;
							if (isHalogen (nonH[0]))
								temp++;
							if (isHalogen (nonH[1]))
								temp++;
							if (temp == 1)
							{
							}
							else if (temp == 2)
							{
								FNumberGeminal = FNumberGeminal + (temp - 1);
							}
						}
						else
						{
							int temp = 0;
							if (isHalogen (nonH[0]))
								temp++;
							if (isHalogen (nonH[1]))
								temp++;
							if (temp == 1)
							{
							}
							else if (temp == 2)
							{
								XNumberGeminal = XNumberGeminal + (temp - 1);
							}
						}
						Hacceptor[numberHacceptors] = k;
						//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
						//molec.number_of_H_Acceptors++;
						numberHacceptors++;
						Hdonor[numberHdonnors] = k;
						//molec.H_bond_donor[molec.number_of_H_Donors] = k;
						//molec.number_of_H_Donors++;
						HdonorNumberHydrogens[numberHdonnors]++;
						numberHdonnors++;
					}
					else
					{
						if (isACarbon (nonH[0]) && isACarbon (nonH[1])
							&& isACarbon (nonH[2]))
						{
							atomsTableInt[k][0] = 51;
							XlogP = XlogP + XLogP_AtomConribution[51];
						}
						else
						{
							atomsTableInt[k][0] = 52;
							XlogP = XlogP + XLogP_AtomConribution[52];
						}
						if (nonH[0] == "F" || nonH[1] == "F" || nonH[2] == "F")
						{
							int temp = 0;
							if (isHalogen (nonH[0]))
								temp++;
							if (isHalogen (nonH[1]))
								temp++;
							if (isHalogen (nonH[2]))
								temp++;
							if (temp == 1)
							{
							}
							else if (temp == 2)
							{
								FNumberGeminal = FNumberGeminal + (temp - 1);
							}
							else if (temp == 3)
							{
								FNumberGeminal = FNumberGeminal + temp;
							}
						}
						else
						{
							int temp = 0;
							if (isHalogen (nonH[0]))
								temp++;
							if (isHalogen (nonH[1]))
								temp++;
							if (isHalogen (nonH[2]))
								temp++;
							if (temp == 1)
							{
							}
							else if (temp == 2)
							{
								XNumberGeminal = XNumberGeminal + (temp - 1);
							}
							else if (temp == 3)
							{
								XNumberGeminal = XNumberGeminal + temp;
							}
						}
						Hacceptor[numberHacceptors] = k;
						//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
						//molec.number_of_H_Acceptors++;
						numberHacceptors++;
					}
				}
				else if (atomsTable[k][0] == "N.2")
				{
					if (nHydrogens[k] == 1)
					{
						atomsTableInt[k][0] = 53;
						XlogP = XlogP + XLogP_AtomConribution[53];
						Hacceptor[numberHacceptors] = k;
						//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
						//molec.number_of_H_Acceptors++;
						numberHacceptors++;
						Hdonor[numberHdonnors] = k;
						//molec.H_bond_donor[molec.number_of_H_Donors] = k;
						//molec.number_of_H_Donors++;
						HdonorNumberHydrogens[numberHdonnors]++;
						numberHdonnors++;
					}
					else if (nHydrogens[k] == 0)
					{
						if (temp2 == 2)
						{       //two neighbors only
							if (isACarbon (nonH[0]) && isACarbon (nonH[1]))
							{
								atomsTableInt[k][0] = 53;
								XlogP = XlogP + XLogP_AtomConribution[53];
							}
							else if (isONSP_HAL (nonH[0]) && isACarbon (nonH[1]))
							{
								if (bondTypes[k][1] == "2")
								{
									if (atomsTable[k][1] == "O.2")
									{   //Adjustment for -ON group
										atomsTableInt[k][0] = 79;
										XlogP =
											XlogP + XLogP_AtomConribution[79] -
											0.218;
									}
									else
									{
										atomsTableInt[k][0] = 55;
										XlogP = XlogP + XLogP_AtomConribution[55];
									}
								}
								else if (bondTypes[k][1] == "1")
								{
									atomsTableInt[k][0] = 54;
									XlogP = XlogP + XLogP_AtomConribution[54];
								}
							}
							else if (isACarbon (nonH[0]) && isONSP_HAL (nonH[1]))
							{
								if (bondTypes[k][2] == "2")
								{
									if (atomsTable[k][2] == "O.2")
									{   //Adjustment for -ON group
										atomsTableInt[k][0] = 79;
										XlogP =
											XlogP + XLogP_AtomConribution[79] -
											0.218;
									}
									else
									{
										atomsTableInt[k][0] = 55;
										XlogP = XlogP + XLogP_AtomConribution[55];
									}
								}
								else if (bondTypes[k][2] == "1")
								{
									atomsTableInt[k][0] = 54;
									XlogP = XlogP + XLogP_AtomConribution[54];
								}
							}
							else if (isONSP_HAL (nonH[0]) && isONSP_HAL (nonH[1]))
							{
								if (nonH[0] == "O.2" || nonH[1] == "O.2")
								{
									atomsTableInt[k][0] = 79;
									XlogP =
										XlogP + XLogP_AtomConribution[79] - 0.218;
								}
								else
								{
									atomsTableInt[k][0] = 56;
									XlogP = XlogP + XLogP_AtomConribution[56];
								}
							}
						}
						else if (temp2 == 3)
						{       //three neighbors
							if ((isAnOxygen (nonH[0]) && isAnOxygen (nonH[1]))
								|| (isAnOxygen (nonH[1]) && isAnOxygen (nonH[2]))
								|| (isAnOxygen (nonH[0]) && isAnOxygen (nonH[2])))
							{
								atomsTableInt[k][0] = 80;
								XlogP = XlogP + XLogP_AtomConribution[80] - 0.436;        // to compensate for oxygens added
							}
						}
						Hacceptor[numberHacceptors] = k;
						//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
						//molec.number_of_H_Acceptors++;
						numberHacceptors++;
					}
				}
				else if (atomsTable[k][0] == "N.ar")
				{
					atomsTableInt[k][0] = 57;
					XlogP = XlogP + XLogP_AtomConribution[57];
				}
				else if (atomsTable[k][0] == "N.pl3")
				{
					if (nHydrogens[k] == 2)
					{
						if (isACarbon (nonH[0]))
						{
							if (isPiBonded (nonH[0]))
							{
								atomsTableInt[k][0] = 47;
								XlogP = XlogP + XLogP_AtomConribution[47];
							}
							else
							{
								atomsTableInt[k][0] = 46;
								XlogP = XlogP + XLogP_AtomConribution[46];
							}
						}
						else
						{
							atomsTableInt[k][0] = 48;
							XlogP = XlogP + XLogP_AtomConribution[48];
						}
						Hacceptor[numberHacceptors] = k;
						//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
						//molec.number_of_H_Acceptors++;
						numberHacceptors++;
						Hdonor[numberHdonnors] = k;
						//molec.H_bond_donor[molec.number_of_H_Donors] = k;
						//molec.number_of_H_Donors++;
						HdonorNumberHydrogens[numberHdonnors]++;
						HdonorNumberHydrogens[numberHdonnors]++;
						numberHdonnors++;
					}
					else if (nHydrogens[k] == 1)
					{
						int firstNeighbor = nonH_INT[0];
						int secondNeighbor = nonH_INT[1];
						int firstNeighbor_neighbors[3];
						int secondNeighbor_neighbors[3];
						bool isFNNACarC2[3];
						bool isSNNACarC2[3];
						bool isNRing = false;
						for (int j = 0; j < 3; j++)
						{
							firstNeighbor_neighbors[j] =
								atomsTableInt[firstNeighbor][j + 1];
							secondNeighbor_neighbors[j] =
								atomsTableInt[secondNeighbor][j + 1];
						} for (int i = 0; i < 3; i++)
						{
							if (secondNeighbor_neighbors[i] == -1
								|| firstNeighbor_neighbors[i] == -1)
							{
							}

							else
							{
								if (atomsTable[firstNeighbor_neighbors[i]][0] !=
									"H")

								{
									isFNNACarC2[i] = true;
								}
								else
								{
									isFNNACarC2[i] = false;
								}
								if (atomsTable[secondNeighbor_neighbors[i]][0] !=
									"H")

								{
									isSNNACarC2[i] = true;
								}
								else
								{
									isSNNACarC2[i] = false;
								}
							}
						}
						for (int i = 0; i < 3; i++)
						{
							for (int j = 0; j < 3; j++)
							{
								if (isFNNACarC2[i] && isSNNACarC2[j])
								{
									for (int K = 0; K < 3; K++)
									{
										for (int L = 0; L < 3; L++)
										{
											if (secondNeighbor_neighbors[j] == -1
												|| firstNeighbor_neighbors[i] ==
												-1)
											{
											}

											else
											{
												if (firstNeighbor_neighbors[i] ==
													atomsTableInt
													[secondNeighbor_neighbors[j]]
												[L])

												{
													isNRing = true;
												}
											}
										}
									}
								}
							}
						}
						if (isNRing)
						{
							atomsTableInt[k][0] = 60;
							XlogP = XlogP + XLogP_AtomConribution[60];
						}
						else
						{
							if (isACarbon (nonH[0]) && isACarbon (nonH[1]))
							{
								atomsTableInt[k][0] = 58;
								XlogP = XlogP + XLogP_AtomConribution[58];
							}
							else if ((isACarbon (nonH[0]) && isONSP_HAL (nonH[1]))
								|| (isONSP_HAL (nonH[0])
								&& isACarbon (nonH[1]))
								|| (isONSP_HAL (nonH[0])
								&& isONSP_HAL (nonH[1])))
							{
								atomsTableInt[k][0] = 59;
								XlogP = XlogP + XLogP_AtomConribution[59];
							}
						}
						Hacceptor[numberHacceptors] = k;
						//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
						//molec.number_of_H_Acceptors++;
						numberHacceptors++;
						Hdonor[numberHdonnors] = k;
						//molec.H_bond_donor[molec.number_of_H_Donors] = k;
						//molec.number_of_H_Donors++;
						HdonorNumberHydrogens[numberHdonnors]++;
						numberHdonnors++;
					}
					else
					{
						int firstNeighbor = nonH_INT[0];
						int secondNeighbor = nonH_INT[1];
						int thirdNeighbor = nonH_INT[2];
						int firstNeighbor_neighbors[3];
						int secondNeighbor_neighbors[3];
						int thirdNeighbor_neighbors[3];
						bool isFNN_N[3];
						bool isSNN_N[3];
						bool isTNN_N[3];
						bool isNRing = false;

						if (isNRing)
						{
							atomsTableInt[k][0] = 62;
							XlogP = XlogP + XLogP_AtomConribution[62];
						}
						else
						{
							atomsTableInt[k][0] = 61;
							XlogP = XlogP + XLogP_AtomConribution[61];
						}
						Hacceptor[numberHacceptors] = k;
						//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
						//molec.number_of_H_Acceptors++;
						numberHacceptors++;
					}
				}
				else if (atomsTable[k][0] == "N.am")
				{
                                //cout << "I AM HERE (1)" << endl;
					isAliphaticORAromatic = false;
					if (nHydrogens[k] == 2)
					{
						Hdonor[numberHdonnors] = k;
						//molec.H_bond_donor[molec.number_of_H_Donors] = k;
						//molec.number_of_H_Donors++;
						HdonorNumberHydrogens[numberHdonnors]++;
						HdonorNumberHydrogens[numberHdonnors]++;
						numberHdonnors++;
						atomsTableInt[k][0] = 63;
						XlogP = XlogP + XLogP_AtomConribution[63];
					}
					else if (nHydrogens[k] == 1)
					{
						Hdonor[numberHdonnors] = k;
						//molec.H_bond_donor[molec.number_of_H_Donors] = k;
						//molec.number_of_H_Donors++;
						HdonorNumberHydrogens[numberHdonnors]++;
						numberHdonnors++;
						atomsTableInt[k][0] = 64;
						XlogP = XlogP + XLogP_AtomConribution[64];
					}
					else
					{
						atomsTableInt[k][0] = 65;
						XlogP = XlogP + XLogP_AtomConribution[65];
					}
					Hacceptor[numberHacceptors] = k;
					//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
					//molec.number_of_H_Acceptors++;
					numberHacceptors++;
				}
				else if (atomsTable[k][0] == "N.1")
				{
					isAliphaticORAromatic = false;
					if (isACarbon (nonH[0]) && !isCNgroup)
					{
						XlogP = XlogP + XLogP_AtomConribution[77];
						isCNgroup = true;
					}
					else
					{
					}
				}
			}
			else if (charge[k] == -1)
			{
				Hacceptor[numberHacceptors] = k;
				//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
				//molec.number_of_H_Acceptors++;
				numberHacceptors++;
				if (nHydrogens[k] == 1 || nHydrogens[k] == 0)
				{
					Hdonor[numberHdonnors] = k;
					//molec.H_bond_donor[molec.number_of_H_Donors] = k;
					//molec.number_of_H_Donors++;
					HdonorNumberHydrogens[numberHdonnors]++;
					numberHdonnors++;
				}
			}
			else if (charge[k] == 1)
			{
				if (atomsTable[k][0] == "N.4")
				{
					isN4inAminoAcid = true;
					N4number[N4Count] = k;
					N4Count++;
					if (nHydrogens[k] == 3)
					{
						if (nonH[0] == "C.ar" || nonH[0] == "C.2"
							|| nonH[0] == "C.1" || nonH[0] == "N.ar"
							|| nonH[0] == "N.2" || nonH[0] == "N.1"
							|| nonH[0] == "N.am" || nonH[0] == "N.pl2"
							|| nonH[0] == "O.2" || nonH[0] == "O.co2")
						{
							XlogP = XlogP - 0.046 - 0.449;
						}
						else
						{
							XlogP = XlogP - 0.046 - 0.582;
						}
						Hacceptor[numberHacceptors] = k;
						//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
						//molec.number_of_H_Acceptors++;
						numberHacceptors++;
						Hdonor[numberHdonnors] = k;
						//molec.H_bond_donor[molec.number_of_H_Donors] = k;
						//molec.number_of_H_Donors++;
						HdonorNumberHydrogens[numberHdonnors]++;
						numberHdonnors++;
					}
					else if (nHydrogens[k] == 2)
					{
						Hacceptor[numberHacceptors] = k;
						//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
						//molec.number_of_H_Acceptors++;
						numberHacceptors++;
						Hdonor[numberHdonnors] = k;
						//molec.H_bond_donor[molec.number_of_H_Donors] = k;
						//molec.number_of_H_Donors++;
						HdonorNumberHydrogens[numberHdonnors]++;
						HdonorNumberHydrogens[numberHdonnors]++;
						numberHdonnors++;

						//Do if needed
					}
					else if (nHydrogens[k] == 1)
					{
						XlogP = XlogP - 0.046 + 0.443;
						Hacceptor[numberHacceptors] = k;
						//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
						//molec.number_of_H_Acceptors++;
						HdonorNumberHydrogens[numberHdonnors]++;
						numberHacceptors++;
					}
					else if (nHydrogens[k] == 0)
					{

						//Do if needed
					}
				}
				else if (atomsTable[k][0] == "N.2")
				{
					if ((isAnOxygen (nonH[0]) && isAnOxygen (nonH[1])) ||
						(isAnOxygen (nonH[1]) && isAnOxygen (nonH[2])) ||
						(isAnOxygen (nonH[0]) && isAnOxygen (nonH[2])))
					{
						atomsTableInt[k][0] = 80;
						XlogP = XlogP + XLogP_AtomConribution[80] - 0.218;        // to compensate for oxygens added
					}
					Hacceptor[numberHacceptors] = k;
					//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
					//molec.number_of_H_Acceptors++;
					numberHacceptors++;
					if (nHydrogens[k] == 2)
					{
						Hdonor[numberHdonnors] = k;
						//molec.H_bond_donor[molec.number_of_H_Donors] = k;
						//molec.number_of_H_Donors++;
						HdonorNumberHydrogens[numberHdonnors]++;
						numberHdonnors++;
					}
				}
			}
		}
		else if (isAnOxygen (atomsTable[k][0]))
		{
                        //cout << "isAnOxygen:" << atomsTable[k][0] << endl;
			isAliphaticORAromatic = false;
			if (charge[k] == 0)
			{
                                //cout << "isAnOxygen:" << atomsTable[k][0] << endl;
				if (atomsTable[k][0] == "O.3")
				{
					if (nHydrogens[k] == 1)
					{
						Hdonor[numberHdonnors] = k;
						//molec.H_bond_donor[molec.number_of_H_Donors] = k;
						//molec.number_of_H_Donors++;
						HdonorNumberHydrogens[numberHdonnors]++;
						numberHdonnors++;
						if (isACarbon (nonH[0]))
						{
							if (isPiBonded (nonH[0]))
							{
								atomsTableInt[k][0] = 39;
								XlogP = XlogP + XLogP_AtomConribution[39];
							}
							else
							{
								atomsTableInt[k][0] = 38;
								XlogP = XlogP + XLogP_AtomConribution[38];
							}
						}
						else
						{
							atomsTableInt[k][0] = 40;
							XlogP = XlogP + XLogP_AtomConribution[40];
						}
					}
					else
					{
						if (isACarbon (nonH[0]) && isACarbon (nonH[1]))
						{
							if (((nonH[0] == "C.ar") || (nonH[0] == "C.2")) &&    // Need to add C.3?? 
								((nonH[1] == "C.ar") || (nonH[1] == "C.2")))
							{
								int firstNeighbor = nonH_INT[0];
								int secondNeighbor = nonH_INT[1];
								int firstNeighbor_neighbors[3];
								int firstNeighbor_neighbor_neighbors[3][3];
								int secondNeighbor_neighbors[3];
								int secondNeighbor_neighbor_neighbors[3][3];
								bool isFNNACarC2[3];
								bool isSNNACarC2[3];
								bool isFuranRing = false;
								for (int j = 0; j < 3; j++)
								{
									firstNeighbor_neighbors[j] =
										atomsTableInt[firstNeighbor][j + 1];
									secondNeighbor_neighbors[j] =
										atomsTableInt[secondNeighbor][j + 1];
								} for (int i = 0; i < 3; i++)
								{
									if ((atomsTable[firstNeighbor_neighbors[i]]
									[0] == "C.ar")
										||
										(atomsTable[firstNeighbor_neighbors[i]]
									[0] == "C.2"))

									{
										isFNNACarC2[i] = true;
									}
									else
									{
										isFNNACarC2[i] = false;
									}
									if ((atomsTable[secondNeighbor_neighbors[i]]
									[0] == "C.ar")
										||
										(atomsTable[secondNeighbor_neighbors[i]]
									[0] == "C.2"))

									{
										isSNNACarC2[i] = true;
									}
									else
									{
										isSNNACarC2[i] = false;
									}
								}
								for (int i = 0; i < 3; i++)
								{
									for (int j = 0; j < 3; j++)
									{
										if (isFNNACarC2[i] && isSNNACarC2[j])
										{
											for (int K = 0; K < 3; K++)
											{
												for (int L = 0; L < 3; L++)
												{
													if (atomsTableInt
														[firstNeighbor_neighbors
														[i]][K] ==
														atomsTableInt
														[secondNeighbor_neighbors
														[j]][L])

													{
														isFuranRing = true;
													}
												}
											}
										}
									}
								}
								if (isFuranRing)
								{
									atomsTableInt[k][0] = 43;
									XlogP = XlogP + XLogP_AtomConribution[43];
								}
								else
								{
									atomsTableInt[k][0] = 41;
									XlogP = XlogP + XLogP_AtomConribution[41];
								}
							}
							else if ((nonH[0] == "C.3") || (nonH[1] == "C.3"))
							{
								atomsTableInt[k][0] = 41;
								XlogP = XlogP + XLogP_AtomConribution[41];
							}
						}
						else
						{       //correction for deprotonated Sp3 Oxygen                                                                                        
							if (atomsTable[k][2] == "nn")
							{
								if (nonH[0] == "C.ar" || nonH[0] == "C.2"
									|| nonH[0] == "C.1")
								{
									atomsTableInt[k][0] = 47;
									XlogP = XlogP + XLogP_AtomConribution[47];
								}
								else
								{
									atomsTableInt[k][0] = 46;
									XlogP = XlogP + XLogP_AtomConribution[46];
								}
							}
							else
							{
								atomsTableInt[k][0] = 42;
								XlogP = XlogP + XLogP_AtomConribution[42];
							}
						}
					}
					Hacceptor[numberHacceptors] = k;
					//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
					//molec.number_of_H_Acceptors++;
					numberHacceptors++;
				}
				else if (atomsTable[k][0] == "O.2")
				{
                                        //cout << "I AM HERE (1)" << endl;
					if (!isNO2group)
					{
						if (temp2 == 2)
						{       //if O.2 has 2 molecules attached as in furan// FIX
                                                        //cout << "I AM HERE (2)" << endl;
							atomsTableInt[k][0] = 43;
							XlogP = XlogP + XLogP_AtomConribution[43];
						}
						else if (isONSP_HAL (nonH[0]))
						{
                                                        //cout << "I AM HERE (3): XLogP_AtomConribution:" << XLogP_AtomConribution[45] <<endl;
							atomsTableInt[k][0] = 45;
							XlogP = XlogP + XLogP_AtomConribution[45];
						}
						else
						{
                                                        //cout << "I AM HERE (4)" << endl;
							atomsTableInt[k][0] = 44;
							XlogP = XlogP + XLogP_AtomConribution[44];
						}
					}
                                        else // trent added 2011.02.20
                                        {// this will need to be cheeked 
                                                 //cout << "I AM HERE (5)" << endl;
                                                 atomsTableInt[k][0] = 44;
                                                 XlogP = XlogP + XLogP_AtomConribution[44];
                                        } 
					Hacceptor[numberHacceptors] = k;
					//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
					//molec.number_of_H_Acceptors++;
					numberHacceptors++;
				}
				else if (atomsTable[k][0] == "O.co2")
				               // || atomsTable[k][0] == "O.spc" || atomsTable[k][0] == "O.t3p")                                      
                                { 
					atomsTableInt[k][0] = 81;
					XlogP = XlogP + XLogP_AtomConribution[81];
					Hacceptor[numberHacceptors] = k;
					//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
					//molec.number_of_H_Acceptors++;
					numberHacceptors++;
				}
			}
			else if (charge[k] == -1)
			{
				if (atomsTable[k][0] == "O.3")
				{
					if (isPiBonded (nonH[0]) || isPiBonded (nonH[1])
						|| isPiBonded (nonH[2]))
					{
						atomsTableInt[k][0] = 39;
						XlogP = XlogP + XLogP_AtomConribution[39] + 0.046;
					}
					else
					{
						atomsTableInt[k][0] = 38;
						XlogP = XlogP + XLogP_AtomConribution[38] + 0.046;
					}
					Hacceptor[numberHacceptors] = k;
					//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
					//molec.number_of_H_Acceptors++;
					numberHacceptors++;
					Hdonor[numberHdonnors] = k;
					//molec.H_bond_donor[molec.number_of_H_Donors] = k;
					//molec.number_of_H_Donors++;
					HdonorNumberHydrogens[numberHdonnors]++;
					numberHdonnors++;
				}
				else if (atomsTable[k][0] == "O.2")
				{
					if (nonH[0] == "N.2")
					{
						if ((atomsTable[atomsTableInt[nonH_INT[0]][1]][0] ==
							"O.2"
							&& atomsTable[atomsTableInt[nonH_INT[0]][2]][0] ==
							"O.2")
							|| (atomsTable[atomsTableInt[nonH_INT[0]][1]][0] ==
							"O.2"
							&& atomsTable[atomsTableInt[nonH_INT[0]][3]][0]
						== "O.2")
							|| (atomsTable[atomsTableInt[nonH_INT[0]][2]][0] ==
							"O.2"
							&& atomsTable[atomsTableInt[nonH_INT[0]][3]][0]
						== "O.2"))
						{
						}
						else
						{
							Hdonor[numberHdonnors] = k;
							//molec.H_bond_donor[molec.number_of_H_Donors] = k;
							//molec.number_of_H_Donors++;
							HdonorNumberHydrogens[numberHdonnors]++;
							numberHdonnors++;
						}
					}
					Hacceptor[numberHacceptors] = k;
					//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
					//molec.number_of_H_Acceptors++;
					numberHacceptors++;
				}
				else if (atomsTable[k][0] == "O.co2")
               // || atomsTable[k][0] == "O.spc" || atomsTable[k][0] == "O.t3p")
				{
					atomsTableInt[k][0] = 81;
					XlogP = XlogP + XLogP_AtomConribution[81];
					Hacceptor[numberHacceptors] = k;
					//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
					//molec.number_of_H_Acceptors++;
					numberHacceptors++;
					Hdonor[numberHdonnors] = k;
					//molec.H_bond_donor[molec.number_of_H_Donors] = k;
					//molec.number_of_H_Donors++;
					HdonorNumberHydrogens[numberHdonnors]++;
					numberHdonnors++;
				}
			}
			else
			{
			}
		}
		else if (atomsTable[k][0] == "H")
		{
			atomsTableInt[k][0] = 37;
			XlogP = XlogP + XLogP_AtomConribution[37];
		}
		else if (atomsTable[k][0] == "S.3")
		{
                        //cout << "I AM HERE :" << atomsTable[k][0] << endl;
			isAliphaticORAromatic = false;
			if (charge[k] == 0)
			{
				temp = 1;
				temp2 = 0;
				vector <string>  nonH(10); //
				int nonH_INT[10];          // sulfurs may have more than 4 bonds
				while (atomsTable[k][temp] != "nn")
				{
                                        //cout << "atomsTable[k][" << temp << "]::"  << atomsTable[k][temp] << endl;
					if (atomsTable[k][temp] != "H")
					{
						nonH[temp2] = atomsTable[k][temp];
						nonH_INT[temp2] = atomsTableInt[k][temp];
						temp2++;
					}
					temp++;
				}               //end while
                                //cout << "number_of_neighbor = " << temp << endl;
				if (nHydrogens[k] == 1)
				{
					Hdonor[numberHdonnors] = k;
					//molec.H_bond_donor[molec.number_of_H_Donors] = k;
					//molec.number_of_H_Donors++;
					HdonorNumberHydrogens[numberHdonnors]++;
					numberHdonnors++;
					atomsTableInt[k][0] = 66;
					XlogP = XlogP + XLogP_AtomConribution[66];
				}
				else
				{
					if (((nonH[0] == "C.ar") || (nonH[0] == "C.2")) &&    // Need to add C.3??
						((nonH[1] == "C.ar") || (nonH[1] == "C.2")))
					{
						int firstNeighbor = nonH_INT[0];
						int secondNeighbor = nonH_INT[1];
						int firstNeighbor_neighbors[3];
						int firstNeighbor_neighbor_neighbors[3][3];
						int secondNeighbor_neighbors[3];
						int secondNeighbor_neighbor_neighbors[3][3];
						bool isFNNAromatic[3];
						bool isSNNAromatic[3];
						bool isSRing = false;
						for (int j = 0; j < 3; j++)
						{
							firstNeighbor_neighbors[j] =
								atomsTableInt[firstNeighbor][j + 1];
							secondNeighbor_neighbors[j] =
								atomsTableInt[secondNeighbor][j + 1];
						} for (int i = 0; i < 3; i++)
						{
							if ((atomsTable[firstNeighbor_neighbors[i]][0] ==
								"C.ar")
								|| (atomsTable[firstNeighbor_neighbors[i]][0] ==
								"C.2")
								|| (atomsTable[firstNeighbor_neighbors[i]][0] ==
								"N.2"))

							{
								isFNNAromatic[i] = true;
							}
							else
							{
								isFNNAromatic[i] = false;
							}
							if ((atomsTable[secondNeighbor_neighbors[i]][0] ==
								"C.ar")
								|| (atomsTable[secondNeighbor_neighbors[i]][0]
							== "C.2")
								|| (atomsTable[firstNeighbor_neighbors[i]][0] ==
								"N.2"))

							{
								isSNNAromatic[i] = true;
							}
							else
							{
								isSNNAromatic[i] = false;
							}
						}
						for (int i = 0; i < 3; i++)
						{
							for (int j = 0; j < 3; j++)
							{
								if (isFNNAromatic[i] && isSNNAromatic[j])
								{
									for (int K = 0; K < 3; K++)
									{
										for (int L = 0; L < 3; L++)
										{
											if (firstNeighbor_neighbors[i] ==
												atomsTableInt
												[secondNeighbor_neighbors[j]][L])

											{
												isSRing = true;
											}
										}
									}
								}
							}
						}
						if (isSRing)
						{
							atomsTableInt[k][0] = 68;
							XlogP = XlogP + XLogP_AtomConribution[68];
						}
						else
						{
							atomsTableInt[k][0] = 67;
							XlogP = XlogP + XLogP_AtomConribution[67];
						}
					}
					else
					{
                                                //cout << "I AM HERE (1)" << endl;
						atomsTableInt[k][0] = 67;
						XlogP = XlogP + XLogP_AtomConribution[67];
					}
				}
			}
			else if (charge[k] == -1)
			{
				Hdonor[numberHdonnors] = k;
				//molec.H_bond_donor[molec.number_of_H_Donors] = k;
				//molec.number_of_H_Donors++;
				HdonorNumberHydrogens[numberHdonnors]++;
				numberHdonnors++;
				atomsTableInt[k][0] = 66;
				XlogP = XlogP + XLogP_AtomConribution[66];
			}
			Hacceptor[numberHacceptors] = k;
			//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
			//molec.number_of_H_Acceptors++;
			numberHacceptors++;
		}
		else if (atomsTable[k][0] == "S.2")
		{
			isAliphaticORAromatic = false;
			temp = 1;
			temp2 = 0;
			vector <string> nonH(4);
			int nonH_INT[4];
			while (atomsTable[k][temp] != "nn")
			{
				if (atomsTable[k][temp] != "H")
				{
					nonH[temp2] = atomsTable[k][temp];
					nonH_INT[temp2] = atomsTableInt[k][temp];
					temp2++;
				}
				temp++;
			}                   //end while
			if (nHydrogens[k] == 1)
			{
				Hdonor[numberHdonnors] = k;
				//molec.H_bond_donor[molec.number_of_H_Donors] = k;
				//molec.number_of_H_Donors++;
				HdonorNumberHydrogens[numberHdonnors]++;
				numberHdonnors++;
				atomsTableInt[k][0] = 66;
				XlogP = XlogP + XLogP_AtomConribution[66];
			}
			else
			{
				if (((nonH[0] == "C.ar") || (nonH[0] == "C.2")) &&        // Need to add C.3??
					((nonH[1] == "C.ar") || (nonH[1] == "C.2")))
				{
					int firstNeighbor = nonH_INT[0];
					int secondNeighbor = nonH_INT[1];
					int firstNeighbor_neighbors[3];
					int firstNeighbor_neighbor_neighbors[3][3];
					int secondNeighbor_neighbors[3];
					int secondNeighbor_neighbor_neighbors[3][3];
					bool isFNNACarC2[3];
					bool isSNNACarC2[3];
					bool isSRing = false;
					for (int j = 0; j < 3; j++)
					{
						firstNeighbor_neighbors[j] =
							atomsTableInt[firstNeighbor][j + 1];
						secondNeighbor_neighbors[j] =
							atomsTableInt[secondNeighbor][j + 1];
					} for (int i = 0; i < 3; i++)
					{
						if ((atomsTable[firstNeighbor_neighbors[i]][0] == "C.ar")
							|| (atomsTable[firstNeighbor_neighbors[i]][0] ==
							"C.2"))

						{
							isFNNACarC2[i] = true;
						}
						else
						{
							isFNNACarC2[i] = false;
						}
						if ((atomsTable[secondNeighbor_neighbors[i]][0] ==
							"C.ar")
							|| (atomsTable[secondNeighbor_neighbors[i]][0] ==
							"C.2"))

						{
							isSNNACarC2[i] = true;
						}
						else
						{
							isSNNACarC2[i] = false;
						}
					}
					for (int i = 0; i < 3; i++)
					{
						for (int j = 0; j < 3; j++)
						{
							if (isFNNACarC2[i] && isSNNACarC2[j])
							{
								for (int K = 0; K < 3; K++)
								{
									for (int L = 0; L < 3; L++)
									{
										if (firstNeighbor_neighbors[i] ==
											atomsTableInt
											[secondNeighbor_neighbors[j]][L])

										{
											isSRing = true;
										}
									}
								}
							}
						}
					}
					if (isSRing)
					{
						atomsTableInt[k][0] = 68;
						XlogP = XlogP + XLogP_AtomConribution[68];
					}
					else
					{
						atomsTableInt[k][0] = 69;
						XlogP = XlogP + XLogP_AtomConribution[69];
					}
				}
				else
				{
					atomsTableInt[k][0] = 69;
					XlogP = XlogP + XLogP_AtomConribution[69];
				}
			}
			Hacceptor[numberHacceptors] = k;
			//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
			//molec.number_of_H_Acceptors++;
			numberHacceptors++;
		}
		else if (atomsTable[k][0] == "S.o")
		{
			isAliphaticORAromatic = false;
			atomsTableInt[k][0] = 70;
			XlogP = XlogP + XLogP_AtomConribution[70];
			Hacceptor[numberHacceptors] = k;
			//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
			//molec.number_of_H_Acceptors++;
			numberHacceptors++;
		}
		else if (atomsTable[k][0] == "S.o2")
		{
			isAliphaticORAromatic = false;
			atomsTableInt[k][0] = 71;
			XlogP = XlogP + XLogP_AtomConribution[71];
			Hacceptor[numberHacceptors] = k;
			//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
			//molec.number_of_H_Acceptors++;
			numberHacceptors++;
		}
		else if (atomsTable[k][0] == "F")
		{
			isAliphaticORAromatic = false;
			atomsTableInt[k][0] = 72;
			XlogP = XlogP + XLogP_AtomConribution[72];
			Hacceptor[numberHacceptors] = k;
			//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
			//molec.number_of_H_Acceptors++;
			numberHacceptors++;
		}
		else if (atomsTable[k][0] == "Cl")
		{
			isAliphaticORAromatic = false;
			atomsTableInt[k][0] = 73;
			XlogP = XlogP + XLogP_AtomConribution[73];
		}
		else if (atomsTable[k][0] == "Br")
		{
                        //cout << "I AM HERE (3)" << endl;
			isAliphaticORAromatic = false;
			atomsTableInt[k][0] = 74;
			XlogP = XlogP + XLogP_AtomConribution[74];
		}
		else if (atomsTable[k][0] == "I")
		{
			isAliphaticORAromatic = false;
			atomsTableInt[k][0] = 75;
			XlogP = XlogP + XLogP_AtomConribution[75];
		}
		else if (atomsTable[k][0] == "P.3")
		{
			Hacceptor[numberHacceptors] = k;
			//molec.H_bond_acceptor[molec.number_of_H_Acceptors] = k;
			//molec.number_of_H_Acceptors++;
			numberHacceptors++;
			isAliphaticORAromatic = false;
			atomsTableInt[k][0] = 76;
			XlogP = XlogP + XLogP_AtomConribution[76];
		}
		else
		{
			atomsTableInt[k][0] = -999;
		}
		if (isN4inAminoAcid && checkIfAminoAcid)
		{
			int C3number;
			int C2number;
			int Cco2 = 0;
			bool N4isConnectedToC3 = false;
			bool N4C3isConnectedToC2 = false;
			bool C2isConnectedToCarboxylate = false;
			int count = 0;
			while (N4number[count] != -1)
			{
				for (int i = 1; i < 5; i++)
				{
					if (atomsTable[N4number[count]][i] == "C.3")
					{
						C3number = atomsTableInt[N4number[count]][i];
						N4isConnectedToC3 = true;
					}
				}
				if (N4isConnectedToC3)
				{
					for (int i = 1; i < 5; i++)
					{
						if (atomsTable[C3number][i] == "C.2")
						{
							C2number = atomsTableInt[C3number][i];
							N4C3isConnectedToC2 = true;
						}
					}
				}
				if (N4C3isConnectedToC2)
				{
					for (int i = 1; i < 4; i++)
					{
						if (atomsTable[C2number][i] == "O.co2")
						{
							Cco2++;
							if (Cco2 == 2)
								C2isConnectedToCarboxylate = true;
						}
					}
				}
				if (N4isConnectedToC3 && N4C3isConnectedToC2
					&& C2isConnectedToCarboxylate)
				{
					XlogP = XlogP - 2.27;
					checkIfAminoAcid = false;
				}
				N4isConnectedToC3 = false;
				N4C3isConnectedToC2 = false;
				C2isConnectedToCarboxylate = false;
				count++;
			}
		}
	}
        //cout << "atomsTableInt[0][1] =  " <<atomsTableInt[0][1] << endl;

        //cout << "####### HERE I AM ()" << endl;
	//for (int d = 0; d < molec.num_atoms; d++){
	//   for (int dd = 0; dd < molec.num_atoms; dd++)
        //        cout << atomsTableInt[d][dd] << ","; 
        //   cout << endl;  
        //}


	float formalCharge = 0;
	for (int d = 0; d < molec.num_atoms; d++)
	{
		//for (int dd = 1; atomsTableInt[d][dd] != -1; dd++)
		for (int dd = 1; dd < molec.num_atoms; dd++)
		{
                        //cout << "d = " << d << "; dd = " << dd << "; atomsTableInt = " << atomsTableInt[d][dd] << endl;
			adjacencyMatrix[d][atomsTableInt[d][dd]] = 1;
			adjacencyMatrix[atomsTableInt[d][dd]][d] = 1;
		} 
                formalCharge = formalCharge + molec.charges[d];
		adjacencyMatrix[d][d] = 0;
	}
        int numHbonds = 0;
	int sepBonds = 0;
	int currentAtom = 0;
	bool found = false;
    molec.hb_donors = numberHdonnors;
    molec.hb_acceptors = numberHacceptors;
    //    cout << "-----------------------------------" << endl;
	//cout << "Molecule Name:\t" << molec.title << endl;
	//cout << "HBond_Donors(OH+NH):\t" << numberHdonnors << endl;
	//cout << "HBond_Acceptors(O+N):\t" << numberHacceptors << endl;

	/*	getBondsBetweenAllAtoms(&adjacencyMatrix[0][0],&bondsBetween[0][0],1,1,0,molec.num_atoms,&space);
	cout<<"Adjacency Matrix"<<endl;
	for(int d = 0; d < molec.num_atoms; d++){		
	for(int dd = 0; dd < molec.num_atoms; dd++){			
	if(adjacencyMatrix[d][dd] < 0){
	cout<<" "<<adjacencyMatrix[d][dd];				
	}else if(adjacencyMatrix[d][dd] > -1 && adjacencyMatrix[d][dd] < 10){
	cout<<"  "<<adjacencyMatrix[d][dd];
	}else if(adjacencyMatrix[d][dd] > 10){
	cout<<" "<<adjacencyMatrix[d][dd];
	}
	}cout<<endl;
	}cout<<endl<<endl;
	for(int d = 0; d < molec.num_atoms; d++){		
	for(int dd = 0; dd < molec.num_atoms; dd++){
	if(bondsBetween[d][dd] < 0){
	cout<<" "<<bondsBetween[d][dd];				
	}else if(bondsBetween[d][dd] > -1 && bondsBetween[d][dd] < 10){
	cout<<"  "<<bondsBetween[d][dd];
	}else if(bondsBetween[d][dd] > 10){
	cout<<" "<<bondsBetween[d][dd];
	}
	}cout<<endl;
	}cout<<endl;
	*/
	int numbHbondPairs = 0;
	for (int i = 0; i < numberHdonnors; i++)
	{
		if (Hdonor[i] != -1)
		{
			for (int j = 0; j < numberHacceptors; j++)
			{

				//If donnors and acceptors form a 6 member ring then they form an H-bond
				//cout<<"Hdonor["<<i<<"] = "<<Hdonor[i]<<" Hacceptor["<<j<<"] = "<<Hacceptor[j]<<endl;
				if (Hacceptor[j] != -1)
				{
					if (Hdonor[i] != Hacceptor[j])
					{
						getBondsBetweenAllAtoms (&adjacencyMatrix[0][0],
							&bondsBetween[0][0], Hdonor[i],
							Hdonor[i], 0, molec.num_atoms);
						if (bondsBetween[Hdonor[i]][Hacceptor[j]] == 4)
						{
							sepBonds = bondsBetween[Hdonor[i]][Hacceptor[j]];
						}
						if (sepBonds == 4)
						{
							if (HdonorNumberHydrogens[i] > 0)
							{
								if (atomsDistance[Hdonor[i]][Hacceptor[j]] < 15.5)
								{
									HdonorNumberHydrogens[i]--;
									numHbonds++;
									HbondPairs[numbHbondPairs][0] = Hdonor[i];
									HbondPairs[numbHbondPairs][1] = Hacceptor[j];
									numbHbondPairs++;
								}
							}
						}
						sepBonds = 0;
					}
				}
			}
		}
	}

	//eliminate redundand H-bonds
	//string nonH_donor[molec.num_atoms][4];
        vector < vector <string> > nonH_donor;
        nonH_donor.resize(molec.num_atoms); 
        for (int i = 0; i<molec.num_atoms; i++){
             nonH_donor[i].resize(4);
        }

	int nonH_INT_donor[molec.num_atoms][4];

	int numb_nonH_donor[molec.num_atoms];

	//string nonH_acceptor[molec.num_atoms][4];
        vector < vector <string> > nonH_acceptor;
        nonH_acceptor.resize(molec.num_atoms);
        for (int i = 0; i<molec.num_atoms; i++){
             nonH_acceptor[i].resize(4);
        }

	int nonH_INT_acceptor[molec.num_atoms][4];

	int numb_nonH_acceptor[molec.num_atoms];
	for (int mm = 0; mm < molec.num_atoms; mm++)
	{
		nonH_donor[mm][0] = nonH_donor[mm][1] = nonH_donor[mm][2] =
			nonH_donor[mm][3] = "nn";
		nonH_INT_donor[mm][0] = nonH_INT_donor[mm][1] = nonH_INT_donor[mm][2] =
			nonH_INT_donor[mm][3] = 0;
		numb_nonH_donor[mm] = 0;
		nonH_acceptor[mm][0] = nonH_acceptor[mm][1] = nonH_acceptor[mm][2] =
			nonH_acceptor[mm][3] = "nn";
		nonH_INT_acceptor[mm][0] = nonH_INT_acceptor[mm][1] =
			nonH_INT_acceptor[mm][2] = nonH_INT_acceptor[mm][3] = 0;
		numb_nonH_acceptor[mm] = 0;
	} for (int i = 0; i < numbHbondPairs; i++)
	{                           //HbondPairs[i][0] != -1; i++){//donor
		int temp = 1;
		int temp2 = 0;
		int temp3 = 1;
		int temp4 = 0;
		while (temp < molec.num_atoms)
		{
			if (atomsTable[HbondPairs[i][0]][temp] != "nn")
			{
				if (atomsTable[HbondPairs[i][0]][temp] != "H")
				{
					nonH_donor[i][temp2] = atomsTable[HbondPairs[i][0]][temp];
					nonH_INT_donor[i][temp2] =
						atomsTableInt[HbondPairs[i][0]][temp];
					numb_nonH_donor[i]++;
					temp2++;
				}
			}
			if (atomsTable[HbondPairs[i][1]][temp3] != "nn")
			{
				if (atomsTable[HbondPairs[i][1]][temp3] != "H")
				{
					nonH_acceptor[i][temp4] = atomsTable[HbondPairs[i][1]][temp3];
					nonH_INT_acceptor[i][temp4] =
						atomsTableInt[HbondPairs[i][1]][temp3];
					numb_nonH_acceptor[i]++;
					temp4++;
				}
			}
			temp++;
			temp3++;
		}                       //end while
	}
	if (numHbonds > 1)
	{
		for (int hh = 0; HbondPairs[hh][0] != -1; hh++)
		{

			//      cout<<"hhh "<<nonH_INT_donor[hh][0]<<" "<<nonH_INT_donor[hh][1]<<" "<<nonH_INT_donor[hh][2]<<" "<<nonH_INT_donor[hh][3]<<endl;
			//      cout<<"hhh "<<nonH_INT_acceptor[hh][0]<<" "<<nonH_INT_acceptor[hh][1]<<" "<<nonH_INT_acceptor[hh][2]<<" "<<nonH_INT_acceptor[hh][3]<<endl;
			for (int hhh = hh + 1; HbondPairs[hhh][0] != -1; hhh++)
			{
				if (numb_nonH_donor[hh] == 2)
				{
					if (numb_nonH_acceptor[hhh] == 3)
					{
						if (nonH_INT_donor[hh][0] == nonH_INT_donor[hhh][0])
						{
							if (nonH_INT_acceptor[hh][0] ==
								nonH_INT_acceptor[hhh][0])
							{
								numHbonds--;
							}
							else if (nonH_INT_acceptor[hh][1] ==
								nonH_INT_acceptor[hhh][1])
							{
								numHbonds--;
							}
							else if (nonH_INT_acceptor[hh][2] ==
								nonH_INT_acceptor[hhh][2])
							{
								numHbonds--;
							}
						}
						else if (nonH_INT_donor[hh][1] == nonH_INT_donor[hhh][1])
						{
							if (nonH_INT_acceptor[hh][0] ==
								nonH_INT_acceptor[hhh][0])
							{
								numHbonds--;
							}
							else if (nonH_INT_acceptor[hh][1] ==
								nonH_INT_acceptor[hhh][1])
							{
								numHbonds--;
							}
							else if (nonH_INT_acceptor[hh][2] ==
								nonH_INT_acceptor[hhh][2])
							{
								numHbonds--;
							}
						}
					}
					else if (numb_nonH_acceptor[hhh] == 2)
					{
						if (nonH_INT_donor[hh][0] == nonH_INT_donor[hhh][0])
						{
							if (nonH_INT_acceptor[hh][0] ==
								nonH_INT_acceptor[hhh][0])
							{
								numHbonds--;
							}
							else if (nonH_INT_acceptor[hh][1] ==
								nonH_INT_acceptor[hhh][1])
							{
								numHbonds--;
							}
						}
						else if (nonH_INT_donor[hh][1] == nonH_INT_donor[hhh][1])
						{
							if (nonH_INT_acceptor[hh][0] ==
								nonH_INT_acceptor[hhh][0])
							{
								numHbonds--;
							}
							else if (nonH_INT_acceptor[hh][1] ==
								nonH_INT_acceptor[hhh][1])
							{
								numHbonds--;
							}
						}
					}
					else if (numb_nonH_acceptor[hhh] == 1)
					{
						if ((nonH_INT_donor[hh][0] == nonH_INT_donor[hhh][0])
							&& (nonH_INT_acceptor[hh][0] ==
							nonH_INT_acceptor[hhh][0]))
						{
							numHbonds--;
						}
						else if ((nonH_INT_donor[hh][1] == nonH_INT_donor[hhh][1])
							&& (nonH_INT_acceptor[hh][0] ==
							nonH_INT_acceptor[hhh][0]))
						{
							numHbonds--;
						}
					}
				}
				else if (numb_nonH_donor[hh] == 1)
				{
					if (numb_nonH_acceptor[hhh] == 3)
					{
						if (nonH_INT_donor[hh][0] == nonH_INT_donor[hhh][0])
						{
							if (nonH_INT_acceptor[hh][0] ==
								nonH_INT_acceptor[hhh][0])
							{
								numHbonds--;
							}
							else if (nonH_INT_acceptor[hh][1] ==
								nonH_INT_acceptor[hhh][1])
							{
								numHbonds--;
							}
							else if (nonH_INT_acceptor[hh][2] ==
								nonH_INT_acceptor[hhh][2])
							{
								numHbonds--;
							}
						}
					}
					else if (numb_nonH_acceptor[hhh] == 2)
					{
						if (nonH_INT_donor[hh][0] == nonH_INT_donor[hhh][0])
						{
							if (nonH_INT_acceptor[hh][0] ==
								nonH_INT_acceptor[hhh][0])
							{
								numHbonds--;
							}
							else if (nonH_INT_acceptor[hh][1] ==
								nonH_INT_acceptor[hhh][1])
							{
								numHbonds--;
							}
						}
					}
					else if (numb_nonH_acceptor[hhh] == 1)
					{
						if ((nonH_INT_donor[hh][0] == nonH_INT_donor[hhh][0])
							&& (nonH_INT_acceptor[hh][0] ==
							nonH_INT_acceptor[hhh][0]))
						{
							numHbonds--;
						}
						else
							if ((nonH_INT_donor[hh][0] == nonH_INT_acceptor[hhh][0])
								|| (nonH_INT_acceptor[hh][0] ==
								nonH_INT_donor[hhh][0]))
							{
								numHbonds--;
							}
					}
				}
			}
		}
	}
	if (isAliphaticORAromatic == true)
	{
		XlogP = XlogP + numOfAliphaticAromaticContributions * 0.19;

		//      cout<<"  Hydrophobic Carbons = "<<numOfAliphaticAromaticContributions;
	}
	XlogP = XlogP + (0.6 * numHbonds);
	XlogP = XlogP + (0.08 * FNumberGeminal) + (-0.26 * XNumberGeminal);

	//      cout<<"   Geminal F = "<<FNumberGeminal<<endl;
	int rot_bonds = 0;
	for (int i = 0; i < molec.num_bonds; i++)
	{

		// cout<<"molec.amber_bt_torsion_total["<<i<<"]="<<molec.amber_bt_torsion_total[i]<<endl;
		if (molec.amber_bt_torsion_total[i] > 0)
			rot_bonds++;
	}

	// Output figures to two decimals places with zeros
	//cout << fixed << showpoint << setprecision(2)     
	//     << "Molecular_Weight:\t" << MWT << endl
	//     << "LogP(xlogP):\t\t" << XlogP << endl
	//     << "#_IntraMol_H-Bonds:\t" << numHbonds << endl
	//     << "Formal_Ligand_Charge:\t" << formalCharge << endl
	//     << "#_Rotatable_Bonds:\t" << rot_bonds << endl << endl;


/*
	cout << "Results for " << molec.title << endl;
	cout << "_#_atoms = " << molec.num_atoms << endl;
	cout << "_#_H-donors = " << molec.number_of_H_Donors << endl;
	cout << "_#_H-acceptors = " << molec.number_of_H_Acceptors << endl;
*/

	int hatom, hacceptor;

	//      for(int g = 0; g < molec.num_atoms; g++){
	/*
	if(atomsTable[g][0] == "H" || atomsTable[g][0] == "H.spc" || atomsTable[g][0] == "H.t3p"){
	if(isANitrogen(atomsTable[g][1]) || isAnOxygen(atomsTable[g][1]) || isASulfur(atomsTable[g][1])){
	for(int gg = 0; gg < molec.num_atoms; gg++){
	if(gg != g){
	if(isANitrogen(atomsTable[gg][0]) || isAnOxygen(atomsTable[gg][0]) || isASulfur(atomsTable[gg][0])){
	hatom = ((molec.x[g]-molec.x[gg])*(molec.x[g]-molec.x[gg]))+
	((molec.y[g]-molec.y[gg])*(molec.y[g]-molec.y[gg]))+
	((molec.z[g]-molec.z[gg])*(molec.z[g]-molec.z[gg]));
	if(hatom < 6.25){
	//   cout<<"Hydrogen bond between"<<atomsTable[g][0]<<" number="<<g<<" and "<<atomsTable[gg][0]<<" number "<<gg<<endl;
	}
	}
	}
	}
	}
	}    
	for(int h = 0; h < molec.num_atoms; h++){ 
	cout<<atomsTable[g][h]<<"   ";
	}cout<<endl;
	for(int t = 0; t < molec.num_atoms; t++){ 
	//cout<<atomsTableInt[g][t]<<"   ";
	}cout<<endl;
	//for(int t = 0; t < molec.num_atoms; t++){ 
	//   cout<<atomsDistance[g][t]<<"   ";
	//}cout<<endl;
	for(int t = 0; t < (molec.num_bonds+1); t++){ 
	cout<<bondTypes[g][t]<<"   ";
	}cout<<endl;
	*/
	//      }
	/*
	cout<<"list of Bonds = ";
	for(int hh = 0; hh < molec.num_atoms+1; hh++){
	cout<<nBonds[hh]<<" ";
	}cout<<endl;
	cout<<"list of Hydrogens = ";
	for(int hhh = 0; hhh < molec.num_atoms; hhh++){
	cout<<nHydrogens[hhh]<<" ";
	}cout<<endl; 
	*/

    molec.xlogp = XlogP;
} 

void xlogp::getBondsBetweenAllAtoms (int *adjacencyMatrix, int *bondsBetween, int baseAtom,
						 int currentAtom, int sepBonds, int numA)
{
	sepBonds++;

	for (int y = 0; y < numA; y++)
	{
		if (adjacencyMatrix[(currentAtom * numA) + baseAtom] == 1)
		{
			sepBonds = 2;
		}
		for (int x = 0; x < numA; x++)
		{
			if (adjacencyMatrix[(currentAtom * numA) + x] == 1)
			{
				if (adjacencyMatrix[(x * numA) + baseAtom] == 1)
				{
					if (sepBonds > 3)
					{
						sepBonds = 3;
					}
				}
			}
		}
		if (adjacencyMatrix[(currentAtom * numA) + y] == -9)
		{
		}
		else if (adjacencyMatrix[(currentAtom * numA) + y] == 0)
		{
			bondsBetween[(currentAtom * numA) + y] =
				bondsBetween[(numA * y) + currentAtom] = 0;
		}
		else
		{
			if (sepBonds < bondsBetween[(baseAtom * numA) + y])
			{
				bondsBetween[(baseAtom * numA) + y] =
					bondsBetween[(y * numA) + baseAtom] = sepBonds;
			}
			else if (bondsBetween[(baseAtom * numA) + y] == -1
				&& bondsBetween[(y * numA) + baseAtom] == -1)
			{
				bondsBetween[(baseAtom * numA) + y] =
					bondsBetween[(y * numA) + baseAtom] = sepBonds;
				getBondsBetweenAllAtoms (adjacencyMatrix, bondsBetween, baseAtom,
					y, sepBonds, numA);
			}
			else
			{
			}
		}
	}
}
