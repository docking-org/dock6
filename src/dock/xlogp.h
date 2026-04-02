#ifndef XLOGP_H
#define XLOGP_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "dockmol.h"
using namespace std;
class xlogp
{
public:
           xlogp ();
	   bool isAnOxygen (string s);
	   bool isANitrogen (string s);
	   bool isASulfur (string s);
	   bool isACarbon (string s);
	   bool isAromatic (string s);
	   bool isONSP_HAL (string s);
	   bool isONSP (string s);
	   bool isHalogen (string s);
	   bool isPiBonded (string s);
	   double getAtomMW (string s);
	   int nHdonnors (string, string, int);
	   void getXLogP (DOCKMol &);
	   void get_H_BOND_Don_Acc (DOCKMol &, string);
	   void getInterMolec_H_Bonds (DOCKMol &, DOCKMol & );
	   void getBondsBetweenAllAtoms (int *, int *, int, int, int, int);
};


#endif // XLOGP_H
