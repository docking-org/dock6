// General headers
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <utility>
#include <fstream>
#include <tuple>
#include <memory>

#ifdef BUILD_DOCK_WITH_RDKIT
 
// RDKit stuff
#include <boost/range/combine.hpp>
#include <RDGeneral/export.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/FilterCatalog/FilterCatalog.h>
#include <DataStructs/SparseIntVect.h>
#include <GraphMol/Fingerprints/MACCS.h>
#include <DataStructs/ExplicitBitVect.h>

// DOCK6 headers
#include "dockmol.h"
#include "rdtyper.h"


/////////////////////////////////
//     Accessory functions     //
/////////////////////////////////

//This is a tuple that stores coefficients of Asymmetric Double Sigmoidal functions for EACH QED parameter. 
//Each column represent the coffcients corresponding in float "ADS( float x , string parameter_type )" function 
std::tuple<double, double, double, double, double, double, double> get_adsParameters(std::string parameter_type) { 
    if (parameter_type == "MW") return std::make_tuple(2.817065973, 392.5754953, 290.7489764, 2.419764353, 49.22325677, 65.37051707, 104.9805561);
    if (parameter_type == "ALOGP") return std::make_tuple(3.172690585, 137.8624751, 2.534937431, 4.581497897, 0.822739154, 0.576295591, 131.3186604);
    if (parameter_type == "HBA") return std::make_tuple(2.948620388, 160.4605972, 3.615294657, 4.435986202, 0.290141953, 1.300669958, 148.7763046);
    if (parameter_type == "HBD") return std::make_tuple(1.618662227, 1010.051101, 0.985094388, 0.000000001, 0.713820843, 0.920922555, 258.1632616);
    if (parameter_type == "PSA") return std::make_tuple(1.876861559, 125.2232657, 62.90773554, 87.83366614, 12.01999824, 28.51324732, 104.5686167);
    if (parameter_type == "ROTB") return std::make_tuple(0.010000000, 272.4121427, 2.558379970, 1.565547684, 1.271567166, 2.758063707, 105.4420403);
    if (parameter_type == "AROM") return std::make_tuple(3.217788970, 957.7374108, 2.274627939, 0.000000001, 1.317690384, 0.375760881, 312.3372610);
    if (parameter_type == "ALERTS") return std::make_tuple(0.010000000, 1199.094025, -0.09002883, 0.000000001, 0.185904477, 0.875193782, 417.7253140);
    throw std::invalid_argument("id");
}

//vectors that store QED properties and different parameter weights in each parameter 
const std::vector <std::string> QED_properties = {"MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS"};
//const std::vector <double> WEIGHT_MAX = {0.50, 0.25, 0.00, 0.50, 0.00, 0.50, 0.25, 1.00};
const std::vector <double> WEIGHT_MEAN = {0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95};
//const std::vector <double> WEIGHT_NONE = {1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00};

// Function that inputs the molecule in question and outputs the number of H-bond acceptors 
// via the definition of "vector <string> AcceptorSmarts". 
int getAcceptorSmarts (RDKit::ROMol &mol ) {
    //initializing variables and vectors
    std::vector <std::string> AcceptorSmarts = {
      "[oH0;X2]",
      "[OH1;X2;v2]",
      "[OH0;X2;v2]",
      "[OH0;X1;v2]",
      "[O-;X1]",
      "[SH0;X2;v2]",
      "[SH0;X1;v2]",
      "[S-;X1]",
      "[nH0;X2]",
      "[NH0;X1;v3]",
      "[$([N;+0;X3;v3]);!$(N[C,S]=O)]"
    };
    int sum_results = 0;
    RDKit::ROMol arlt;
    int num_accvec = AcceptorSmarts.size(); 
    std::vector<std::vector <std::pair<int,int>>> matchVect = {};
    //testing if the mol has any "Acceptors"
    for(int i=0; i < num_accvec; i++){
        //initiaizing a pointer because SmartsToMol returns a pointer
        RDKit::RWMol * rdmolListOfacc = RDKit::SmartsToMol(AcceptorSmarts[i]); 
        sum_results = sum_results + RDKit::SubstructMatch( mol, *rdmolListOfacc, matchVect, true, true, false, false, 1000, 1);
        //delete pointer
        delete rdmolListOfacc;
    }
    matchVect.clear();

    return sum_results;
}

// Function that inputs the molecule in question and outputs the number of StructuralAlertSmarts
// via the definition of "vector <string> StructuralAlertSmarts ".
int getStructuralALERTS (RDKit::ROMol &mol) {
    //initializing variables and vectors
    std::vector <std::string>  StructuralAlertSmarts {
        "*1[O,S,N]*1",
        "[S,C](=[O,S])[F,Br,Cl,I]",
        "[CX4][Cl,Br,I]",
        "[#6]S(=O)(=O)O[#6]",
        "[$([CH]),$(CC)]#CC(=O)[#6]",
        "[$([CH]),$(CC)]#CC(=O)O[#6]",
        "n[OH]",
        "[$([CH]),$(CC)]#CS(=O)(=O)[#6]",
        "C=C(C=O)C=O",
        "n1c([F,Cl,Br,I])cccc1",
        "[CH1](=O)",
        "[#8][#8]",
        "[C;!R]=[N;!R]",
        "[N!R]=[N!R]",
        "[#6](=O)[#6](=O)",
        "[#16][#16]",
        "[#7][NH2]",
        "C(=O)N[NH2]",
        "[#6]=S",
        "[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]=[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]",
        "C1(=[O,N])C=CC(=[O,N])C=C1",
        "C1(=[O,N])C(=[O,N])C=CC=C1",
        "a21aa3a(aa1aaaa2)aaaa3",
        //"a1aa2a(aa3c(S2)aaaa3)aa1",
        "a31a(a2a(aa1)aaaa2)aaaa3",
        "a1aa2a3a(a1)A=AA=A3=AA=A2",
        "c1cc([NH2])ccc1",
        "[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si,Na,Ca,Ge,Ag,Mg,K,Ba,Sr,Be,Ti,Mo,Mn,Ru,Pd,Ni,Cu,Au,Cd]",
        "[Al,Ga,Sn,Rh,Tl,Bi,Nb,Li,Pb,Hf,Ho]",
        "I",
        "OS(=O)(=O)[O-]",
        "[N+](=O)[O-]",
        "C(=O)N[OH]",
        "C1NC(=O)NC(=O)1",
        "[SH]",
        "[S-]",
        "c1ccc([Cl,Br,I,F])c([Cl,Br,I,F])c1[Cl,Br,I,F]",
        "c1cc([Cl,Br,I,F])cc([Cl,Br,I,F])c1[Cl,Br,I,F]",
        "[CR1]1[CR1][CR1][CR1][CR1][CR1][CR1]1",
        "[CR1]1[CR1][CR1]cc[CR1][CR1]1",
        "[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2][CR2]1",
        "[CR2]1[CR2][CR2]cc[CR2][CR2][CR2]1",
        "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
        "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
        "C#C",
        "[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]",
        "[$([N+R]),$([n+R]),$([N+]=C)][O-]",
        "[#6]=N[OH]",
        "[#6]=NOC=O",
        "[#6](=O)[CX4,CR0X3,O][#6](=O)",
        "c1ccc2c(c1)ccc(=O)o2",
        "[O+,o+,S+,s+]",
        "N=C=O",
        "[NX3,NX4][F,Cl,Br,I]",
        "c1ccccc1OC(=O)[#6]",
        "[CR0]=[CR0][CR0]=[CR0]",
        "[C+,c+,C-,c-]",
        "N=[N+]=[N-]",
        "C12C(NC(N1)=O)CSC2",
        "c1c([OH])c([OH,NH2,NH])ccc1",
        "P",
        "[N,O,S]C#N",
        "C=C=O",
        "[Si][F,Cl,Br,I]",
        "[SX2]O",
        "[SiR0,CR0](c1ccccc1)(c2ccccc2)(c3ccccc3)",
        "O1CCCCC1OC2CCC3CCCCC3C2",
        "N=[CR0][N,n,O,S]",
        "[cR2]1[cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2][cR2]1[cR2]2[cR2][cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2]2",
        "C=[C!r]C#N",
        "[cR2]1[cR2]c([N+0X3R0,nX3R0])c([N+0X3R0,nX3R0])[cR2][cR2]1",
        "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2]c([N+0X3R0,nX3R0])[cR2]1",
        "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2][cR2]c1([N+0X3R0,nX3R0])",
        "[OH]c1ccc([OH,NH2,NH])cc1",
        "c1ccccc1OC(=O)O",
        "[SX2H0][N]",
        "c12ccccc1(SC(S)=N2)",
        "c12ccccc1(SC(=S)N2)",
        "c1nnnn1C=O",
        "s1c(S)nnc1NC=O",
        "S1C=CSC1=S",
        "C(=O)Onnn",
        "OS(=O)(=O)C(F)(F)F",
        "N#CC[OH]",
        "N#CC(=O)",
        "S(=O)(=O)C#N",
        "N[CH2]C#N",
        "C1(=O)NCC1",
        "S(=O)(=O)[O-,OH]",
        "NC[F,Cl,Br,I]",
        "C=[C!r]O",
        "[NX2+0]=[O+0]",
        "[OR0,NR0][OR0,NR0]",
        "C(=O)O[C,H1].C(=O)O[C,H1].C(=O)O[C,H1]",
        "[CX2R0][NX3R0]",
        "c1ccccc1[C;!R]=[C;!R]c2ccccc2",
        "[NX3R0,NX4R0,OR0,SX2R0][CX4][NX3R0,NX4R0,OR0,SX2R0]",
        "[s,S,c,C,n,N,o,O]~[n+,N+](~[s,S,c,C,n,N,o,O])(~[s,S,c,C,n,N,o,O])~[s,S,c,C,n,N,o,O]",
        "[s,S,c,C,n,N,o,O]~[nX3+,NX3+](~[s,S,c,C,n,N])~[s,S,c,C,n,N]",
        "[*]=[N+]=[*]",
        "[SX3](=O)[O-,OH]",
        "N#N",
        "F.F.F.F",
        "[R0;D2][R0;D2][R0;D2][R0;D2]",
        "[cR,CR]~C(=O)NC(=O)~[cR,CR]",
        "C=!@CC=[O,S]",
        "[#6,#8,#16][#6](=O)O[#6]",
        "c[C;R0](=[O,S])[#6]",
        "c[SX2][C;!R]",
        "C=C=C",
        "c1nc([F,Cl,Br,I,S])ncc1",
        "c1ncnc([F,Cl,Br,I,S])c1",
        "c1nc(c2c(n1)nc(n2)[F,Cl,Br,I])",
        "[#6]S(=O)(=O)c1ccc(cc1)F",
        "[15N]",
        "[13C]",
        "[18O]",
        "[34S]"
    };
    int sum_results = 0;
    RDKit::ROMol arlt;
    int num_altvec = StructuralAlertSmarts.size();
    std::vector <std::pair<int,int>> matchVect = {}; 
    //testing if the mol has any "StructuralAlerts"
    for(int i=0; i < num_altvec; i++){  
        //initializing pointer because SmartsToMol returns pointers
        RDKit::RWMol *rdmolListOfStructAlerts = RDKit::SmartsToMol( StructuralAlertSmarts[i]);
        if (RDKit::SubstructMatch(mol, *rdmolListOfStructAlerts,matchVect)){
             //This line of code returns the results with the number of matches, but in QED.py in RDKIT
             //shows that it is only if it there or not. So it should be a boolean.
             //sum_results = sum_results + RDKit::SubstructMatch( mol, *rdmolListOfStructAlerts, matchVect, true, true, false, false, 1000, 1);
            sum_results +=1;
        }
        //delete pointer
        delete rdmolListOfStructAlerts;
    }
    matchVect.clear();
 
return sum_results;
}

// This is the core of QED scoring. Inputs: individual QED parameters calculated by 
// "double calculate_qed_score(string smi, vector <float> w)". 
// Outputs: results stemming from the ADS function.  
double ADS(double x, std::string parameter_type) {
    //initializing parameters and variables
    double        exp1{};
    double        exp2{};
    double          dx{};
    double results_ads{};

    //gets ADS parameters from get_adsParameters
    auto ADS_parameters = get_adsParameters(parameter_type);

    //calculating the ADS
    exp1 = 1. + exp(-1.*(x - std::get<2>(ADS_parameters) + std::get<3>(ADS_parameters)/2.)/std::get<4>(ADS_parameters));
    exp2 = 1. + exp(-1.*(x - std::get<2>(ADS_parameters) - std::get<3>(ADS_parameters)/2)/std::get<5>(ADS_parameters));
    dx = std::get<0>(ADS_parameters) + std::get<1>(ADS_parameters)/exp1*(1 - 1/exp2);
    results_ads = dx/std::get<6>(ADS_parameters);

    return results_ads;
}


//calculating number of PAINS and the names of the PAINS associated with the PAINS hits
//this is changed because RDKit has a way to calculate the pains using RDKIT functions only
int calculate_pns (DOCKMol &mol) {
    int total_results {0};
    RDKit::FilterCatalogParams params;
    std::vector<RDKit::FilterCatalog::CONST_SENTRY> matches{};

    params.addCatalog(RDKit::FilterCatalogParams::PAINS_A);
    params.addCatalog(RDKit::FilterCatalogParams::PAINS_B);
    params.addCatalog(RDKit::FilterCatalogParams::PAINS_C);
    RDKit::FilterCatalog catalog(params);

    RDKit::ROMol *romol = RDKit::SmilesToMol(mol.smiles, 0, true, nullptr);

    if(catalog.hasMatch(*romol)){
        matches = catalog.getMatches(*romol);
         
        total_results = matches.size();
        for(int i = 0; i<matches.size(); i++){
            std::string pain_name = matches[i]->getDescription();
            mol.pns_name.push_back(pain_name+" ");
        }
    }
    
    delete romol;

    return total_results;
}

int getSSSR(RDKit::ROMol &mol) {
    RDKit::VECT_INT_VECT rings;
    int nr = RDKit::MolOps::findSSSR(mol, rings);
    return nr;
}


/////////////////////////////////////
////     RDTyper class methods     //
/////////////////////////////////////

int RDTYPER::get_pns (DOCKMol &mol, RDKit::ROMol &rdmol, std::map<std::string,std::string> &PAINSmap) {
    //initializing a variable = 0
    int score {0};
  
    //calculate pains  
    score = calculate_pns(mol);

    return score;
}

double RDTYPER::calculate_sa_score(RDKit::ROMol &rdmol, int nChiralCenters, int rtn_nSpi, std::map<unsigned int, double> &fragMap) {

    //variables initialize
    double                         score1            {0};
    double                         score2            {0};
    double                         score3            {0};
    int                            nf                {0};
    unsigned int                   sfp               {0};
    double                         tmp_val          {-4}; 
    double                         nAtoms            {0};
    RDKit::RingInfo*               ri;
    double                         rtn_nBdghd;
    double                         nMacrocycles      {0};
    std::vector<std::vector<int>>  atom_rngs          {};
    double                         sizePenalty       {0};
    double                         stereoPenalty     {0};
    double                         spiroPenalty      {0};
    double                         bridgePenalty     {0};
    double                         macrocyclePenalty {0};
    double                         min               {0};
    double                         max               {0};
    double                         sascore           {0};

    //SASCORE = SCORE 1 - SCORE 2 - SCORE 3
    //SCORE 1
    //calculates MorganFingerprints and sift the molecular fingerprints through the 
    //list of fragments (loaded up above) to count the number of recognized fragments
 
    //initializing pointers 
    RDKit::SparseIntVect<unsigned int> *fp = RDKit::MorganFingerprints::getFingerprint(rdmol, 2); 
    
    const std::map<unsigned int, int, std::less<unsigned int>, std::allocator<std::pair<const unsigned int, int>>> fps = fp ->getNonzeroElements();
    //for loops through the fragment file to recognize fragments
    for (auto &temp: fps) {
        nf = nf + temp.second;
        //temp.second is v
        sfp = temp.first;
        //temp.first is bitId 
        auto search = fragMap.find(sfp);
        if (search != fragMap.end()){
            score1 = score1 + (fragMap[sfp] * temp.second);
        } else {
            score1 = score1 + ( tmp_val * temp.second);
        }
        
    }
    //DELETE POINTERS
    delete fp;

    score1 = score1 / nf;
    nAtoms = rdmol.getNumAtoms();
    nChiralCenters = RDKit::Descriptors::numAtomStereoCenters( rdmol );
    ri = rdmol.getRingInfo();
    rtn_nSpi = RDKit::Descriptors::calcNumSpiroAtoms( rdmol );
    rtn_nBdghd = RDKit::Descriptors::calcNumBridgeheadAtoms( rdmol );

    atom_rngs = ri -> atomRings();
    for(auto tmp_atomRings : atom_rngs){
        if (tmp_atomRings.size() > 8) {
            nMacrocycles++;
        }
    }
   
    //SCORE 2
    //the score that penalizes SCORE 1
    sizePenalty = pow(nAtoms,1.005) - nAtoms;
    stereoPenalty = log10(nChiralCenters + 1.);
    spiroPenalty = log10(rtn_nSpi + 1.);
    bridgePenalty = log10(rtn_nBdghd  + 1.);

    if (nMacrocycles > 0) {
        macrocyclePenalty = log10(2.);
    }

    //calculates SCORE 2
    score2 = 0.0 - sizePenalty - stereoPenalty - spiroPenalty - bridgePenalty - macrocyclePenalty;

    //SCORE 3
    //the score penalizes SCORE 1
    //if the number of ATOMS in the molecule of question is bigger than the 
    //number of recognizes fragments, then PENALIZED
    if ( nAtoms > fps.size() ) {
        score3 = log(nAtoms / fps.size()) * 0.5;
    }
     
    //calculates sascore in total
    sascore = score1 + score2 + score3;
    //please uncomment these if you want to know the breakdown of all scores
    //std::cout << "score1: " << score1 << std::endl;
    //std::cout << "score2: " << score2 << std::endl;
    //std::cout << "score3: " << score3 << std::endl; 
    //normalizes SAscore
    min = -4.0;
    max = 2.5;
    sascore = 11.0 - (sascore - min + 1) / (max - min) * 9.0;
    if (sascore > 8.0) {
        sascore = 8.0 + log(sascore + 1.0 - 9.0);
    }
    if (sascore > 10.0) {
        sascore = 10.0;
    } else if (sascore < 1.0) {
        sascore = 1.0;

    }
    
    return sascore;
}
// END SA SCORE

double RDTYPER::calculate_qed_score(RDKit::ROMol &rdmol, double CLOGP, double TPSA, std::vector <double> w) {

    //variables initialize
    double     MW       {};
    double     HBA      {};
    double     HBD      {};
    double     ROTB     {};
    double     AROM     {};
    double     ALERTS   {};

    //initializing pointers
    RDKit::ROMol *AliphaticRings = RDKit::SmartsToMol("[$([A;R][!a])]");   
    RDKit::ROMol *removeHrdmol = {RDKit::MolOps::removeHs( rdmol )};  
    RDKit::ROMol *deleteSubs = RDKit::deleteSubstructs(rdmol, *AliphaticRings);
   
    //initial QED parameters calculations   
    MW = RDKit::Descriptors::calcAMW( *removeHrdmol );
    HBA = getAcceptorSmarts ( *removeHrdmol );
    HBD = RDKit::Descriptors::calcNumHBD( *removeHrdmol );
    ROTB = RDKit::Descriptors::calcNumRotatableBonds( *removeHrdmol , true );
    //AROM = RDKit::MolOps::findSSSR( *deleteSubs );
    AROM = getSSSR(*deleteSubs); 
    ALERTS = getStructuralALERTS( *removeHrdmol );
    std::vector <double> QED_properties_numbers {MW, CLOGP, HBA, HBD, TPSA, ROTB, AROM, ALERTS};
 
    //uncomment if you want to print each QED property
    //for(auto itr : QED_properties_numbers)
    //    cout << itr << endl;
    //    

    //DELETE POINTERS
    delete deleteSubs;
    delete removeHrdmol; 

    // In case hydrogen atoms were not removed before
    //RDKit::MolOps::removeHs( rdmol );
    std::vector <double> d {};
    double        t {0};
    double    w_sum {0};

    //calculates QED score by calling the "float  ADS( float x , string parameter_type )" function
    for (auto tup : boost::combine( QED_properties_numbers , QED_properties) ) {
        std::string prop;
        double prop_numbers;
        boost::tie(prop_numbers, prop) = tup;

        d.push_back(ADS(prop_numbers,prop.c_str()));

    }  
   
    //after the ADS function is executed, weights for the coefficients will be applied to QED score 
    for (auto tup1 : boost::combine( w , d) ) {
        double wi {0};
        double di {0};
        boost::tie(wi, di) = tup1;

        t = t + ( wi  * log(di));
    }
    
    //sums up all of the components of the ADS function for each QED parameter and spits out results
    for (auto itr : w)
        w_sum = w_sum + itr;
    double results = exp(t / w_sum);
   
    //deleting pointers 
    delete AliphaticRings; 
   
    return results;
} 
// END QED



double RDTYPER::calculate_esol(RDKit::ROMol &rdmol, double CLOGP, bool delaney) {
    
    //initializes coefficients for ESOL.
    //These coefficients are RDKit optimized
    double esol_score        {0};
    double esol_intercept    {0.26121066137801696};
    double esol_coef_logp    {-0.7416739523408995};
    double esol_coef_mw      {-0.0066138847738667125};
    double esol_coef_rotors  {0.003451545565957996};
    double esol_coef_ap      {-0.42624840441316975}; 
    double esol_ap_counter   {0};
     
    //These are the original coeffice    nts
    double esol_orig_intercept       {0.16};
    double esol_coef_orig_logp       {-0.63};
    double esol_coef_orig_mw         {-0.0062};
    double esol_coef_orig_rotors     {0.066};
    double esol_coef_orig_ap         {-0.74};
    double esol_ap_orig_counter      {0};
    
    // variables initialized 
    double MW     {};               
    double ROTB   {};      
    double AP     {};       
    std::vector <std::vector <std::pair<int, int> > > matchVect1;
    
    //pointers initialized
    RDKit::RWMol *mol_ptr = RDKit::SmartsToMol("a" , 0 , false , 0);

    //ESOL parameters calcualted 
    MW = RDKit::Descriptors::calcAMW( rdmol );
    ROTB = RDKit::Descriptors::calcNumRotatableBonds( rdmol, true );
 
    //logic to calculate the RDKIT optimized ESOL or the original version
    if (delaney == true) {
    esol_ap_counter = esol_ap_counter + 
    RDKit::SubstructMatch( rdmol, *mol_ptr, matchVect1, true, true, false, false, 1000, 1);
    AP = esol_ap_counter / rdmol.getNumAtoms();
    esol_score = esol_intercept           + 
                 esol_coef_logp   * CLOGP + 
                 esol_coef_mw     * MW    + 
                 esol_coef_rotors * ROTB  + 
                 esol_coef_ap     * AP; 
    } else {
    esol_ap_orig_counter = esol_ap_orig_counter +  RDKit::SubstructMatch( rdmol, *mol_ptr, matchVect1, true, true, false, false, 1000, 1);
    AP = esol_ap_orig_counter / rdmol.getNumAtoms();
    esol_score = esol_orig_intercept           + 
                 esol_coef_orig_logp   * CLOGP + 
                 esol_coef_orig_mw     * MW    + 
                 esol_coef_orig_rotors * ROTB  + 
                 esol_coef_orig_ap     * AP; 
    }  

    //delete pointers
    delete mol_ptr; 

    return esol_score; 
}

void RDTYPER::calculate_descriptors(DOCKMol &mol, std::map<unsigned int, double> &fragMap, bool create_smiles, std::map<std::string,std::string> &PAINSmap){    
    if (create_smiles == true){
        // Start calculation
        try {
            // Create RDKit mol object
            RDKit::ROMol tmp_rdmol = mol.DOCKMol_to_ROMol( create_smiles );
            if (mol.smiles == ""){
                throw mol.title + "_has_no_SMILES_string"; 
            }
            // Create RDKit mol object that needs to be used (no H)
            // intializing pointers
            RDKit::RWMol *rwmol = RDKit::SmilesToMol(mol.smiles, 0, true, nullptr);
            RDKit::ROMol noH_rdmol {*rwmol};

            //initialize datastructs that is unique to MACCSfingerprinting 
            //ExplicitBitVect *MACCSfp = RDKit::MACCSFingerprints::getFingerprintAsBitVect(tmp_rdmol);
            ExplicitBitVect *MACCSfp = RDKit::MACCSFingerprints::getFingerprintAsBitVect(noH_rdmol);
            boost::dynamic_bitset<> tmp_bitset = *(MACCSfp->dp_bits); 

            mol.num_arom_rings = RDKit::Descriptors::calcNumAromaticRings(noH_rdmol);
            mol.num_alip_rings = RDKit::Descriptors::calcNumAliphaticRings(noH_rdmol);
            mol.num_sat_rings = RDKit::Descriptors::calcNumSaturatedRings(noH_rdmol);
            mol.num_stereocenters = RDKit::Descriptors::numAtomStereoCenters(noH_rdmol);
            mol.num_spiro_atoms = RDKit::Descriptors::calcNumSpiroAtoms(noH_rdmol);
            mol.tpsa = RDKit::Descriptors::calcTPSA(noH_rdmol);
            mol.clogp = RDKit::Descriptors::calcClogP(noH_rdmol);



            //This uses a rdkit mol object that comes from SMILES. 
            mol.qed_score = calculate_qed_score(noH_rdmol, mol.clogp, mol.tpsa, WEIGHT_MEAN);
            mol.sa_score = calculate_sa_score(noH_rdmol, mol.num_stereocenters, mol.num_spiro_atoms, fragMap);
            mol.esol = calculate_esol(noH_rdmol, mol.clogp, true);
            mol.pns = get_pns(mol, tmp_rdmol, PAINSmap);



            //to assign MACCS fingerprints
            mol.MACCS = tmp_bitset;
            mol.MACCS_size = tmp_bitset.size();

            //delete MACCS pointer
            delete MACCSfp;

            //delete pointer
            delete rwmol;

        } catch (const char* message){
            std::cout << message << std::endl;
            mol.smiles = mol.title + " RDKIT ERROR";
            mol.bad_molecule = true;
          
        } catch (...) {
            // RDKit could not read this molecule
            mol.smiles = mol.title + " RDKIT ERROR";
            mol.bad_molecule = true;
        }
    } else {
        std::cout << "Warning: SMILES are needed to calculate the descriptors" << std::endl;
        mol.smiles = mol.title + " NO SMILES";
    }
         
}

void RDTYPER::calculate_descriptors_drive(DOCKMol &mol, std::map<unsigned int, double> fragMap, bool create_smiles, std::map<std::string,std::string> PAINSmap){    
    if (create_smiles == true){
        // Start calculation
        try {
            // Create RDKit mol object
            RDKit::ROMol tmp_rdmol = mol.DOCKMol_to_ROMol( create_smiles );
            if (mol.smiles == ""){
                throw mol.title + "_has_no_SMILES_string"; 
            }
            // Create RDKit mol object that needs to be used (no H)
            // intializing pointers
            RDKit::RWMol *rwmol = RDKit::SmilesToMol(mol.smiles, 0, true, nullptr);
            RDKit::ROMol noH_rdmol {*rwmol};

            //initialize datastructs that is unique to MACCSfingerprinting 
            //ExplicitBitVect *MACCSfp = RDKit::MACCSFingerprints::getFingerprintAsBitVect(tmp_rdmol);
            ExplicitBitVect *MACCSfp = RDKit::MACCSFingerprints::getFingerprintAsBitVect(noH_rdmol);
            boost::dynamic_bitset<> tmp_bitset = *(MACCSfp->dp_bits); 

            mol.num_arom_rings = RDKit::Descriptors::calcNumAromaticRings(noH_rdmol);
            mol.num_alip_rings = RDKit::Descriptors::calcNumAliphaticRings(noH_rdmol);
            mol.num_sat_rings = RDKit::Descriptors::calcNumSaturatedRings(noH_rdmol);
            mol.num_stereocenters = RDKit::Descriptors::numAtomStereoCenters(noH_rdmol);
            mol.num_spiro_atoms = RDKit::Descriptors::calcNumSpiroAtoms(noH_rdmol);
            mol.tpsa = RDKit::Descriptors::calcTPSA(noH_rdmol);
            mol.clogp = RDKit::Descriptors::calcClogP(noH_rdmol);



            //This uses a rdkit mol object that comes from SMILES. 
            mol.qed_score = calculate_qed_score(noH_rdmol, mol.clogp, mol.tpsa, WEIGHT_MEAN);
            mol.sa_score = calculate_sa_score(noH_rdmol, mol.num_stereocenters, mol.num_spiro_atoms, fragMap);
            mol.esol = calculate_esol(noH_rdmol, mol.clogp, true);
            mol.pns = get_pns(mol, tmp_rdmol, PAINSmap);



            //to assign MACCS fingerprints
            mol.MACCS = tmp_bitset;
            mol.MACCS_size = tmp_bitset.size();

            //delete MACCS pointer
            delete MACCSfp;

            //delete pointer
            delete rwmol;

        } catch (const char* message){
            std::cout << message << std::endl;
            mol.smiles = mol.title + " RDKIT ERROR";
            mol.bad_molecule = true;
          
        } catch (...) {
            // RDKit could not read this molecule
            mol.smiles = mol.title + " RDKIT ERROR";
            mol.bad_molecule = true;
        }
    } else {
        std::cout << "Warning: SMILES are needed to calculate the descriptors" << std::endl;
        mol.smiles = mol.title + " NO SMILES";
    }
         
}

#else
  // Configured to Not BUILD_DOCK_WITH_RDKIT
#endif
