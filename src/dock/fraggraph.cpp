#include <iostream>
#include <string.h>
#include "dockmol.h"
#include "fraggraph.h"

using namespace std;


// +++++++++++++++++++++++++++++++++++++++++
// Some constructors and destructors

FragGraph::FragGraph(){
    visited = false;
}
FragGraph::~FragGraph(){
    tanvec.clear();
    rankvec.clear();
}

