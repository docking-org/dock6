#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include "dockmol.h"


//std::string map_gasteiger_type(std::string);
const double converg  = 0.00001;
const int    max_iter = 11;

bool  rms_converg(double *, double *, int);
float oco2charge(DOCKMol &, int);
float o2charge(DOCKMol &, int);
float n2charge(DOCKMol &, int);
float sulfcharge(DOCKMol &, int);
float compute_gast_charges(DOCKMol &);

