#include "dockmol.h"

#include <fstream>

#ifdef BUILD_DOCK_WITH_RDKIT
#include "rdtyper.h"
#include "rddrive.h"



/// RD_Calc

RD_Calc::Calc::Calc(){

    defCreate_smiles = true;

    char result[ 256]; 
    ssize_t count = readlink( "/proc/self/exe", result, 256); 

    std::string EXE_PATH( result, (count > 0) ? count : 0 ) ; 
    std::size_t botDirPos = EXE_PATH.find_last_of("/");
    std::string BINDIR    = EXE_PATH.substr(0, botDirPos); 
    botDirPos             = BINDIR.find_last_of("/");
    std::string BASEDIR   = BINDIR.substr(0, botDirPos); 

    std::ostringstream oss_pains;
    oss_pains << BASEDIR << "/parameters/pains_table_2019_09_01.dat";
    std::string paramater_path_pains = oss_pains.str();

    // Add PAINSMap to memory once
    std::vector<std::string> PAINStmp;
    if (PAINStmp.empty() == true) {
        std::ifstream fin( paramater_path_pains );
        std::string tmp_string;
        while (fin >> tmp_string) {
            PAINStmp.push_back(tmp_string);
        }
        fin.close();
    }
    //if the PAINStmp is still empty, quit 
    if (PAINStmp.empty() == true){
        std::string tmp_string;
        std::cout << "pains_table.dat is empty. Please give the correct PAINS table...";
        exit(0);
    }

    //processing the pains maps.
    int numofvectors = PAINStmp.size() - 1;
    for (int i{0}; i < numofvectors; ) {
        defPAINSmap[PAINStmp[i]] = " " + PAINStmp[i + 1];
        i += 2;
    }
    PAINStmp.clear();




    std::ostringstream oss_sa_fraglib;
    oss_sa_fraglib << BASEDIR << "/parameters/sa_fraglib.dat";
    std::string paramater_path_sa_fraglib = oss_sa_fraglib.str();

    // Add fragMap to memory once
    if ( defFragMap.empty() == true ) {
        std::ifstream fin( paramater_path_sa_fraglib );
        double key;
        double val;
        while (fin >> key >> val) {
            defFragMap[key] = val;
        }
        //the fragMap is empty so quit 
        if ( defFragMap.empty() == true){
            std::cout << "sa_fraglib.dat is empty. Please give the correct fraglib data table...";
            exit(0);
        }
        fin.close();
    }

}

RD_Calc::Calc::~Calc(){

    defCreate_smiles = false;
    defPAINSmap      = {};
    defFragMap       = {};
}

void 
RD_Calc::Calc::calc_des( DOCKMol & mol, RD_Parm::RD_Drive_Parm  RD_parm = RD_Parm::RD_Drive_Parm()) {

    this->rd_typer.calculate_descriptors_drive( mol , 
                                          RD_parm.get_fragMap(),
                                          RD_parm.get_createsmiles(),
                                          RD_parm.get_PAINSmap());

    return;
}


// bool drive_clogp, bool drive_esol, bool drive_qed, bool _drive_sa, 
// and bool drive_stereocenters are attributes of the DN_Build object
bool RD_Calc::Drive::drive_growth( DOCKMol temp_molecule, 
                                   RD_Parm::RD_Drive_Parm  RD_parm = RD_Parm::RD_Drive_Parm() ){
    //Initialize LOCAL variables    
    bool drive_clogp;
    bool drive_esol;
    bool drive_qed;
    bool drive_sa;
    bool drive_stereocenters;
    bool drive_pains;
    bool drive_tpsa;

    // Does temp_molecule passes filters?
    //    
    // Remember that only one false is enough to make this 
    // method return `false`. That's why the `else` clause
    // should be always true. Remember of TRUTH TABLES.
    //    
    // STEREOCENTERS
    drive_stereocenters = stereocenter_cutoff( temp_molecule, RD_parm.get_upper_stereocenter() );
    temp_molecule.fail_stereo = !drive_stereocenters;

    // TPSA
    drive_tpsa = tpsa_cutoff( temp_molecule, RD_parm.get_lower_tpsa(),
                                             RD_parm.get_upper_tpsa(), 
                                             RD_parm.get_tpsa_std_dev());
    temp_molecule.fail_tpsa = !drive_tpsa;

    // CLOGP
    drive_clogp = clogp_cutoff( temp_molecule, RD_parm.get_lower_clogp(),
                                               RD_parm.get_upper_clogp(),
                                               RD_parm.get_clogp_std_dev());
    temp_molecule.fail_clogp = !drive_clogp;

    // ESOL
    drive_esol = esol_cutoff( temp_molecule, RD_parm.get_lower_esol(),
                                             RD_parm.get_upper_esol(),
                                             RD_parm.get_esol_std_dev());
    temp_molecule.fail_esol = !drive_esol;

    // QED
    drive_qed = qed_cutoff( temp_molecule, RD_parm.get_lower_qed(), 
                                           RD_parm.get_qed_std_dev());
    temp_molecule.fail_qed = !drive_qed;

    // SA 
    drive_sa = sa_cutoff( temp_molecule, RD_parm.get_upper_sa(),
                                         RD_parm.get_sa_std_dev());
    temp_molecule.fail_sa = !drive_sa;

    // PAINS
    drive_pains = pains_cutoff( temp_molecule, RD_parm.get_upper_pains() );
    temp_molecule.fail_pains = !drive_pains;

    // Analyze truth value and return result
    if (drive_clogp && drive_esol && drive_tpsa && drive_qed &&
        drive_sa && drive_stereocenters && drive_pains) {
        // These are AND operations. If a single one of these booleans is
        // `false,` this method returns `false.`
        return true; 
    } else { return false; }
} // End  RD_Calc::Drive::drive_growth()

bool RD_Calc::Drive::stereocenter_cutoff( DOCKMol temp_molecule , int upper_stereocenter ){ 
    // Initialize variables
    float  rand_num1{};
    float  rand_num2{};
    float  rand_num1_dec{};
    float  rand_num2_dec{};

    // smoothing #stereocenters
    if (temp_molecule.num_stereocenters > upper_stereocenter ) {
        rand_num1 = (rand() % 100 + 1); 
        rand_num1_dec = rand_num1 / 100;
        rand_num2 = (rand() % 100 + 1) / (temp_molecule.num_stereocenters - 1); 
        rand_num2_dec = rand_num2 / 100;
        if (rand_num1_dec < rand_num2_dec) {
            return true; 
        } else {
            return false;
        }
    } else {
        return true; 
    }
} // End RD_Calc::Drive::stereocenter_cutoff()

bool RD_Calc::Drive::tpsa_cutoff( DOCKMol temp_molecule, float lower_tpsa, float upper_tpsa, float tpsa_std_dev ){
    // Initialize variables
    float  rand_num{};
    float  rand_num_dec{};
    double excessTPSA{};
    double Z_scoreExcess{};
    double acceptRate{};

    // soft cutoff
    if (temp_molecule.tpsa > lower_tpsa) {

        // Check if value is higher than upper boundary
        if (temp_molecule.tpsa > upper_tpsa) {

            // Do soft cutoff procedure
            rand_num = (rand() % 100 + 1); 
            rand_num_dec = rand_num / 100;
            excessTPSA = temp_molecule.tpsa - upper_tpsa;
            Z_scoreExcess = excessTPSA / tpsa_std_dev;
            acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
            if (rand_num_dec < acceptRate) {
                return true; // molecule accepted and given to push_back
            } else {
                return false;
            }
        } else { // mol.tpsa > lower_tpsa && mol.tpsa < upper_esol
            return true; 
        }
    } else { // mol.tpsa < lower_tpsa. Needs soft cutoff

        // Do soft cutoff procedure
        rand_num = (rand() % 100 + 1); 
        rand_num_dec = rand_num / 100;
        excessTPSA = lower_tpsa - temp_molecule.tpsa;
        Z_scoreExcess = excessTPSA / tpsa_std_dev;
        acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
        if (rand_num_dec < acceptRate) {
            return true; // molecule accepted and given to push_back
        } else {
            return false;
        }
    }

} // End RD_Calc::Drive::tpsa_cutoff()

bool RD_Calc::Drive::clogp_cutoff( DOCKMol temp_molecule, float lower_clogp, float upper_clogp, float clogp_std_dev ){
    // Define variables
    float  rand_num{};
    float  rand_num_dec{};
    double excessCLOGP{};
    double Z_scoreExcess{};
    double acceptRate{};
    
    // directing clogp

    // soft cutoff
    if (temp_molecule.clogp > lower_clogp) {

        // Check if value is higher than upper boundary
        if (temp_molecule.clogp > upper_clogp) {

            // Do soft cutoff procedure
            rand_num = (rand() % 100 + 1); 
            rand_num_dec = rand_num / 100;
            excessCLOGP = temp_molecule.clogp - upper_clogp;
            Z_scoreExcess = excessCLOGP / clogp_std_dev;
            acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
            if (rand_num_dec < acceptRate) {
                return true; // molecule accepted and given to push_back
            } else {
                return false;
            }
        } else { // mol.clogp > lower_clogp && mol.clogp < upper_clogp
            return true; 
        }
    } else { // mol.clogp < lower_clogp. Needs soft cutoff

        // Do soft cutoff procedure
        rand_num = (rand() % 100 + 1); 
        rand_num_dec = rand_num / 100;
        excessCLOGP = lower_clogp - temp_molecule.clogp;
        Z_scoreExcess = excessCLOGP / clogp_std_dev;
        acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
        if (rand_num_dec < acceptRate) {
            return true; // molecule accepted and given to push_back
        } else {
            return false;
        }
    }
      
} // End RD_Calc::Drive::clogp_cutoff()

bool RD_Calc::Drive::esol_cutoff( DOCKMol temp_molecule, float lower_esol, float upper_esol, float esol_std_dev ){
    // Initialize variables
    float  rand_num{};
    float  rand_num_dec{};
    double excessESOL{};
    double Z_scoreExcess{};
    double acceptRate{};

    // soft cutoff
    if (temp_molecule.esol > lower_esol) {

        // Check if value is higher than upper boundary
        if (temp_molecule.esol > upper_esol) {

            // Do soft cutoff procedure
            rand_num = (rand() % 100 + 1);
            rand_num_dec = rand_num / 100;
            excessESOL = temp_molecule.esol - upper_esol;
            Z_scoreExcess = excessESOL / esol_std_dev;
            acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
            if (rand_num_dec < acceptRate) {
                return true; // molecule accepted and given to push_back
            } else {
                return false;
            }
        } else { // mol.esol > lower_esol && mol.esol < upper_esol
            return true;
        }
    } else { // mol.esol < lower_esol. Needs soft cutoff

        // Do soft cutoff procedure
        rand_num = (rand() % 100 + 1);
        rand_num_dec = rand_num / 100;
        excessESOL = lower_esol - temp_molecule.esol;
        Z_scoreExcess = excessESOL / esol_std_dev;
        acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
        if (rand_num_dec < acceptRate) {
            return true; // molecule accepted and given to push_back
        } else {
            return false;
        }
    }

} // End RD_Calc::Drive::esol_cutoff()

bool RD_Calc::Drive::qed_cutoff( DOCKMol temp_molecule, float lower_qed, float qed_std_dev){ 
    // Initialize variables
    float  rand_num{};
    float  rand_num_dec{};
    double excessQED{};
    double Z_scoreExcess{};
    double acceptRate{};

    // directing qed
    if (temp_molecule.qed_score < lower_qed) {
        rand_num = (rand() % 100 + 1);
        rand_num_dec = rand_num / 100;
        excessQED = lower_qed - temp_molecule.qed_score;
        Z_scoreExcess = excessQED  / qed_std_dev;
        acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
        if (rand_num_dec < acceptRate) {
            return true; // molecule accepted and given to push_back
        } else {
            return false;
        }
    } else {
        return true; // molecule accepted and given to push_back
    }
} // End RD_Calc::Drive::qed_cutoff()

bool RD_Calc::Drive::sa_cutoff( DOCKMol temp_molecule, float upper_sa, float sa_std_dev){
    // Initialize variables
    float  rand_num{};
    float  rand_num_dec{};
    double excessSA{};
    double Z_scoreExcess{};
    double acceptRate{};

    // directing sa
    if (temp_molecule.sa_score > upper_sa) {

        // proper method
        rand_num = (rand() % 100 + 1);
        rand_num_dec = rand_num / 100;
        excessSA = temp_molecule.sa_score - upper_sa;
        Z_scoreExcess = excessSA / sa_std_dev;
        acceptRate = exp(-1 * Z_scoreExcess * Z_scoreExcess);
        if (rand_num_dec < acceptRate) {
            return true;
        } else {
            return false;
        }
    } else {
        return true;
    }
} // End RD_Calc::Drive::sa_cutoff() 

bool RD_Calc::Drive::pains_cutoff( DOCKMol temp_molecule, float upper_pains ){ 
    // Initialize variables
    float  rand_num1{};
    float  rand_num2{};
    float  rand_num1_dec{};
    float  rand_num2_dec{};

    // smoothing #pains
    if (temp_molecule.pns > upper_pains) {
        rand_num1 = (rand() % 100 + 1);
        rand_num1_dec = rand_num1 / 100;
        rand_num2 = (rand() % 100 + 1) / (temp_molecule.pns - 1);
        rand_num2_dec = rand_num2 / 100;
        if (rand_num1_dec < rand_num2_dec) {
            return true;
        } else {
            return false;
        }
    } else {
        return true;
    }
} // End RD_Calc::Drive::pains_cutoff()


//// End RD_Calc


//// RD_Parm::RD_Drive_Parm
RD_Parm::RD_Drive_Parm::RD_Drive_Parm(){

    //RD_Parm::RD_Drive_Parm::set_createsmiles( true );
    create_smiles = true;

    char result[ 256];
    ssize_t count = readlink( "/proc/self/exe", result, 256);

    std::string EXE_PATH( result, (count > 0) ? count : 0 ) ;
    std::size_t botDirPos = EXE_PATH.find_last_of("/");
    std::string BINDIR    = EXE_PATH.substr(0, botDirPos);
    botDirPos             = BINDIR.find_last_of("/");
    std::string BASEDIR   = BINDIR.substr(0, botDirPos);

    std::ostringstream oss_pains;
    oss_pains << BASEDIR << "/parameters/pains_table_2019_09_01.dat";
    std::string paramater_path_pains = oss_pains.str();

    // Add PAINSMap to memory once
    std::vector<std::string> PAINStmp;
    if (PAINStmp.empty() == true) {
        std::ifstream fin( paramater_path_pains );
        std::string tmp_string;
        while (fin >> tmp_string) {
            PAINStmp.push_back(tmp_string);
        }
        fin.close();
    }   
    //if the PAINStmp is still empty, quit 
    if (PAINStmp.empty() == true){
        std::string tmp_string;
        std::cout << "pains_table.dat is empty. Please give the correct PAINS table...";
        exit(0);
    }   

    //processing the pains maps.
    int numofvectors = PAINStmp.size() - 1;
    for (int i{0}; i < numofvectors; ) { 
        PAINSmap[PAINStmp[i]] = " " + PAINStmp[i + 1]; 
        i += 2;
    }   
    PAINStmp.clear();




    std::ostringstream oss_sa_fraglib;
    oss_sa_fraglib << BASEDIR << "/parameters/sa_fraglib.dat";
    std::string paramater_path_sa_fraglib = oss_sa_fraglib.str();

    // Add fragMap to memory once
    if ( fragMap.empty() == true ) { 
        std::ifstream fin( paramater_path_sa_fraglib );
        double key;
        double val;
        while (fin >> key >> val) {
            fragMap[key] = val;
        }
        //the fragMap is empty so quit 
        if ( fragMap.empty() == true){
            std::cout << "sa_fraglib.dat is empty. Please give the correct fraglib data table...";
            exit(0);
        }
        fin.close();
    }



    
    this->lower_tpsa = 28.53;
    this->upper_tpsa = 113.20;
    this->tpsa_std_dev = 42.33;

    this->lower_clogp = -0.30;
    this->upper_clogp = 3.75;
    this->clogp_std_dev = 2.02;

    this->lower_esol = -5.23;
    this->upper_esol = -1.35;
    this->esol_std_dev = 1.94;

    this->lower_qed = 0.61;
    this->qed_std_dev = 0.19;

    this->upper_sa = 3.34;
    this->sa_std_dev = 0.9;

    this->upper_pains = 1;
}


RD_Parm::RD_Drive_Parm::~RD_Drive_Parm(){

    this->create_smiles = true;

    this->fragMap = {};
    this->PAINSmap ={};

    this->lower_tpsa    = 0; 
    this->upper_tpsa    = 0; 
    this->tpsa_std_dev  = 0; 

    this->lower_clogp   = 0; 
    this->upper_clogp   = 0; 
    this->clogp_std_dev = 0; 

    this->lower_esol    = 0; 
    this->upper_esol    = 0; 
    this->esol_std_dev  = 0; 

    this->lower_qed     = 0; 
    this->qed_std_dev   = 0; 

    this->upper_sa      = 0; 
    this->sa_std_dev    = 0; 

    this->upper_pains   = 0; 
}
void
RD_Parm::RD_Drive_Parm::set_fragMap( std::map<unsigned int, double> input_fragMap ){
   this->fragMap = input_fragMap; 
}

std::map<unsigned int, double>
RD_Parm::RD_Drive_Parm::get_fragMap( ){
    return this->fragMap;
}

void
RD_Parm::RD_Drive_Parm::set_createsmiles( bool input_createsmiles ){
   this->create_smiles = input_createsmiles; 
}

bool
RD_Parm::RD_Drive_Parm::get_createsmiles( ){
    return this->create_smiles;
}

void
RD_Parm::RD_Drive_Parm::set_PAINSmap( std::map<std::string,std::string>  input_PAINSmap ){
   this->PAINSmap = input_PAINSmap; 
}

std::map<std::string,std::string>
RD_Parm::RD_Drive_Parm::get_PAINSmap( ){
    return this->PAINSmap;
}

// STEREOCENTER
void
RD_Parm::RD_Drive_Parm::set_upper_stereocenter( int input_upper_stereocenter ) {
   this->upper_stereocenter = input_upper_stereocenter ;
}

int
RD_Parm::RD_Drive_Parm::get_upper_stereocenter(){
    return this->upper_stereocenter;
}


// TPSA
void
RD_Parm::RD_Drive_Parm::set_lower_tpsa( float input_lower_tpsa ) {
   this->lower_tpsa = input_lower_tpsa ;
}

float
RD_Parm::RD_Drive_Parm::get_lower_tpsa(){
    return this->lower_tpsa;
}

void
RD_Parm::RD_Drive_Parm::set_upper_tpsa( float input_upper_tpsa ) {
   this->upper_tpsa = input_upper_tpsa ;
}

float
RD_Parm::RD_Drive_Parm::get_upper_tpsa(){
    return this->upper_tpsa;
}

void
RD_Parm::RD_Drive_Parm::set_tpsa_std_dev( float input_tpsa_std_dev) {
   this->tpsa_std_dev = input_tpsa_std_dev;
}

float
RD_Parm::RD_Drive_Parm::get_tpsa_std_dev(){
    return this->tpsa_std_dev;
}

// CLOGP
void
RD_Parm::RD_Drive_Parm::set_lower_clogp( float input_lower_clogp ) {
   this->lower_clogp = input_lower_clogp ;
}

float
RD_Parm::RD_Drive_Parm::get_lower_clogp(){
    return this->lower_clogp;
}

void
RD_Parm::RD_Drive_Parm::set_upper_clogp( float input_upper_clogp ) {
   this->upper_clogp = input_upper_clogp ;
}

float
RD_Parm::RD_Drive_Parm::get_upper_clogp(){
    return this->upper_clogp;
}

void
RD_Parm::RD_Drive_Parm::set_clogp_std_dev( float input_clogp_std_dev) {
   this->clogp_std_dev = input_clogp_std_dev;
}

float
RD_Parm::RD_Drive_Parm::get_clogp_std_dev(){
    return this->clogp_std_dev;
}

// ESOL
void
RD_Parm::RD_Drive_Parm::set_lower_esol( float input_lower_esol ) {
   this->lower_esol = input_lower_esol ;
}

float
RD_Parm::RD_Drive_Parm::get_lower_esol(){
    return this->lower_esol;
}

void
RD_Parm::RD_Drive_Parm::set_upper_esol( float input_upper_esol ) {
   this->upper_esol = input_upper_esol ;
}

float
RD_Parm::RD_Drive_Parm::get_upper_esol(){
    return this->upper_esol;
}
void
RD_Parm::RD_Drive_Parm::set_esol_std_dev( float input_esol_std_dev) {
   this->esol_std_dev = input_esol_std_dev;
}

float
RD_Parm::RD_Drive_Parm::get_esol_std_dev(){
    return this->esol_std_dev;
}

// QED
void
RD_Parm::RD_Drive_Parm::set_lower_qed( float input_lower_qed ) {
   this->lower_qed = input_lower_qed ;
}

float
RD_Parm::RD_Drive_Parm::get_lower_qed(){
    return this->lower_qed;
}

void
RD_Parm::RD_Drive_Parm::set_qed_std_dev( float input_qed_std_dev) {
   this->qed_std_dev = input_qed_std_dev;
}

float
RD_Parm::RD_Drive_Parm::get_qed_std_dev(){
    return this->qed_std_dev;
}

// SYNTHA
void
RD_Parm::RD_Drive_Parm::set_upper_sa( float input_upper_sa ) {
   this->upper_sa = input_upper_sa ;
}

float
RD_Parm::RD_Drive_Parm::get_upper_sa(){
    return this->upper_sa;
}

void
RD_Parm::RD_Drive_Parm::set_sa_std_dev( float input_sa_std_dev) {
   this->sa_std_dev = input_sa_std_dev;
}

float
RD_Parm::RD_Drive_Parm::get_sa_std_dev(){
    return this->sa_std_dev;
}

// PAINS

void
RD_Parm::RD_Drive_Parm::set_upper_pains( int input_upper_pains ) {
   this->upper_pains = input_upper_pains ;
}

int
RD_Parm::RD_Drive_Parm::get_upper_pains(){
    return this->upper_pains;
}

//// End RD_Drive_Parm



bool
RD_Drive::drive( DOCKMol test_mol ){

    DOCKMol copy_mol;

    copy_molecule( copy_mol, test_mol );

    this->calc_des.calc_des( copy_mol );

    return this->drive_des.drive_growth( copy_mol );
}

bool
RD_Drive::drive( DOCKMol test_mol, RD_Parm::RD_Drive_Parm& RD_parm ){

    DOCKMol copy_mol;

    copy_molecule( copy_mol, test_mol );

    this->calc_des.calc_des( copy_mol );

    return this->drive_des.drive_growth( copy_mol, RD_parm );
}
#else
  // Configured to Not BUILD_DOCK_WITH_RDKIT
#endif
