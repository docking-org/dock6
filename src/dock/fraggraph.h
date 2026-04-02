#ifndef FRAGGRAPH_H
#define FRAGGRAPH_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "amber_typer.h"
#include "dockmol.h"
#include "master_score.h"
#include "utils.h"

// +++++++++++++++++++++++++++++++++++++++++
// Fragment Graph: A data structure that stores a multi-directional graph of similarities for
// fragment libraries based on Tanimoto (or something similar).
class           FragGraph {
   public:

       std::vector < std::pair <float, int> >   tanvec;
       std::vector <int>                        rankvec;
       bool                                     visited;

       FragGraph();
       ~FragGraph();
};


#endif //FRAGGRAPH_H

