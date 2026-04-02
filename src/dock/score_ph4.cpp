#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include "score_ph4.h"
#include "amber_typer.h"

using namespace std;

// +++++++++++++++++++++++++++++++++++++++++
// This function is used to locate ring members of one ring (up to 7 member ring)
bool
Ph4_Score::find_ring_member(DOCKMol & mol, int & id, vector <int> & visited,vector <int> & ring_member)
{
 ring_member.push_back(id);//include the current atom in the current ring

 if ((ring_member.size() > 7)){
     return false;//if the total number of atoms in the current ring, then stop and return false (not a ring) 
 }

 vector <int> nbrs;
 nbrs.clear();

 for (int i=0; i < mol.num_bonds; i++) {
      if (mol.bonds_origin_atom[i] == id){
          nbrs.push_back(mol.bonds_target_atom[i]);
      }
      else if (mol.bonds_target_atom[i] == id){
          nbrs.push_back(mol.bonds_origin_atom[i]);
      }
 }

 bool old_atom, through;
 
 for (int i=0; i < nbrs.size(); i++) {
      if (nbrs[i] == ring_member[0]){
          if(ring_member.size()>=3 && ring_member.size()<=7) {
             for (int j=0; j<ring_member.size(); j++) {
                  visited.push_back(ring_member[j]);//save atoms in all identified rings to "visited"
             }
             return true;//find a ring!
          }
      }
      old_atom = false;
      for (int j=0; j<ring_member.size(); j++) {
           if (nbrs[i] == ring_member[j]) {
               old_atom = true;//don't go to an atom that's already in the current ring.
               break;
           }
      }
      if (!old_atom){ 
          //if(mol.ph4_types[nbrs[i]].compare("aromatic")==0 || mol.ph4_types[nbrs[i]].compare("aroAcc")==0) {
          //used the above if statement to only account for atoms that are labeled as aromatic.
             through = find_ring_member(mol, nbrs[i], visited, ring_member);//iterativly searching for ring member
             if(through){
                return true;
             }
             else {//clear ring_member when find_ring_member doesn't end up finding a ring.
                int j= ring_member.size();
                ring_member.pop_back();
             }
          //}
      }
 }        
 return false;
} // end Ph4_Score::find_ring_member()


// +++++++++++++++++++++++++++++++++++++++++
//read in Ph4_ref from a ph4.txt file. to be modified by lingling
//PH4_ELEMENT temp_ph4; //need to define Ph4_ELEMENT in dockmol
bool read_ph4_txt(Ph4_Struct & Ph4_ref, istream &ifs)
{     
     string       line; //read in individual lines from ref_file.

     std::vector <int>          id;
     std::vector <float>        x_crd, y_crd, z_crd, v_x, v_y, v_z, radius;
     std::vector <std::string>  ph4_types;

     int   i;
     float k;
     std::string label;

     Ph4_ref.num_features=0;
     while (getline(ifs,line,'\n')) {
            vector < string > tokens;

            if (!line.empty()) {
                Tokenizer(line, tokens, ' ');
                if (tokens.size() > 0){
                    if (!(tokens[0].compare("ph4name") == 0) && !(tokens[0].compare("#") == 0)){
                        //ignore the column names and comments lines
                        if (tokens.size() != 9){//terminate if wrong number of coloms in a line.
                            cout << "A line in the pharmamcophore reference txt file has an incorrect"
                                 << endl
                                 << "number of entries. Lines must have 9 entries separated by spaces as follows:"
                                 << endl
                                 << "ph4name ph4id coox cooy cooz vecx vecy vecz radius"
                                 << endl
                                 << "Note that lines beginning with a comment \"#\" will be ignored"
                                 << endl
                                 << "Program will terminate."
                                 << endl;
                            cout << "line which fails on:" <<endl
                                 << line << endl;
                            exit(1);
                        }
                        stringstream ss; ss.clear();
                        label = tokens[0];       ph4_types.push_back(label);
                        ss <<tokens[1]; ss >> i; id.push_back(i);         ss.clear();
                        ss <<tokens[2]; ss >> k; x_crd.push_back(k);      ss.clear();
                        ss <<tokens[3]; ss >> k; y_crd.push_back(k);      ss.clear();
                        ss <<tokens[4]; ss >> k; z_crd.push_back(k);      ss.clear();
                        ss <<tokens[5]; ss >> k; v_x.push_back(k);        ss.clear();
                        ss <<tokens[6]; ss >> k; v_y.push_back(k);        ss.clear();
                        ss <<tokens[7]; ss >> k; v_z.push_back(k);        ss.clear();
                        ss <<tokens[8]; ss >> k; radius.push_back(k);     ss.clear();
                        Ph4_ref.num_features++;
                    }
                }
            }
            tokens.clear();
     }

     if (Ph4_ref.num_features > 0){
         Ph4_ref.initialize();
         for(i=0; i<Ph4_ref.num_features; i++){
             Ph4_ref.id[i]=id[i];
             Ph4_ref.x[i]=x_crd[i];
             Ph4_ref.y[i]=y_crd[i];
             Ph4_ref.z[i]=z_crd[i];
             Ph4_ref.v_x[i]=v_x[i];
             Ph4_ref.v_y[i]=v_y[i];
             Ph4_ref.v_z[i]=v_z[i];
             Ph4_ref.radius[i]=radius[i];
             if (ph4_types[i].compare("PHO")==0){
                 Ph4_ref.ph4_types[i]="hydrophobic";
             }
             if (ph4_types[i].compare("HBD")==0){
                 Ph4_ref.ph4_types[i]="donor";
             }
             if (ph4_types[i].compare("HBA")==0){
                 Ph4_ref.ph4_types[i]="acceptor";
             }
             if (ph4_types[i].compare("ARO")==0){
                 Ph4_ref.ph4_types[i]="aromatic";
             }
             if (ph4_types[i].compare("RNG")==0){
                 Ph4_ref.ph4_types[i]="ring";
             }
             if (ph4_types[i].compare("POS")==0){
                 Ph4_ref.ph4_types[i]="positive";
             }
             if (ph4_types[i].compare("NEG")==0){
                 Ph4_ref.ph4_types[i]="negative";
             }
         }
         id.clear();
         x_crd.clear(); y_crd.clear(); z_crd.clear();
         v_x.clear(); v_y.clear(); v_z.clear();
         radius.clear(); 
         ph4_types.clear();
           
         return true;
     }

     return false;
} //end read_ph4_txt()

// +++++++++++++++++++++++++++++++++++++++++
// determine which ph4 feature the atom ph4_type corresponds to.
int component_number( string s, vector <string> list)
{
    for (int i=0; i < list.size(); i++){
         if (s == list[i]){
             return i;
         }
    }

    cout << "Error: feature " << s << " not defined in the ph4.defn" << endl;
    exit(0);

} //end component_number() 

// +++++++++++++++++++++++++++++++++++++++++
Ph4_Struct::Ph4_Struct()
{
   x               = NULL;
   y               = NULL;
   z               = NULL;
   v_x             = NULL;
   v_y             = NULL;
   v_z             = NULL;
   radius          = NULL;
   id              = NULL;
   ph4_types       = NULL;
   contri_to_score = NULL;
}

// +++++++++++++++++++++++++++++++++++++++++
Ph4_Struct::~Ph4_Struct()
{
   this->clear();
}

// +++++++++++++++++++++++++++++++++++++++++
//initialize the labeled atom set or the ph4 model.
//Need to write another function to initialize the ph4 model
//max num_features? dummy atom v.s vector
void
Ph4_Struct::initialize()
{  
   

   x          = new float  [num_features];
   y          = new float  [num_features];
   z          = new float  [num_features];
   v_x        = new float  [num_features];
   v_y        = new float  [num_features];
   v_z        = new float  [num_features];
   radius     = new float  [num_features];
   id         = new int    [num_features];
   ph4_types  = new string [num_features];
   contri_to_score = new bool [num_features];
} // end Ph4_Struct::initialize()

// +++++++++++++++++++++++++++++++++++++++++
//clear contents of c++ strings
void
Ph4_Struct::clear()
{
   delete [] x;
   delete [] y;
   delete [] z;
   delete [] v_x;
   delete [] v_y;
   delete [] v_z;
   delete [] radius;
   delete [] id;
   delete [] ph4_types;
   delete [] contri_to_score;
} //end Ph4_Struct::clear()

// +++++++++++++++++++++++++++++++++++++++++
Ph4_Score::Ph4_Score()
{
   match_comp = NULL;
   components = NULL;
   match_num  = NULL;
   match_term = NULL;
}

// +++++++++++++++++++++++++++++++++++++++++
Ph4_Score::~Ph4_Score()
{
   this ->close();
}

// +++++++++++++++++++++++++++++++++++++++++
// When this class is closed, delete some of the C++ strings.
void
Ph4_Score::close()
{
   delete [] match_comp;
   delete [] components;
   delete [] match_num;
   delete [] match_term;
   list_comp.clear();
   constant_ph4_ref_file.clear();
   ref_ph4_out_mol2_filename.clear();
   ref_ph4_out_txt_filename.clear();
   cad_ph4_out_filename.clear();
   mat_ph4_out_filename.clear();
   ph4_compare_type.clear();
   //output_score_summary.clear();
} //Ph4_Score::close()

// +++++++++++++++++++++++++++++++++++++++++
// Input parameters for pharmacophore similarity scoring function
// called in master score
// standard dock scoring function component, will show up as entries in the dock.in files
void
Ph4_Score::input_parameters(Parameter_Reader & parm,
                                        bool & primary_score,
                                        bool & secondary_score)
{
    string          tmp;

    use_primary_score = false;
    use_secondary_score = false;
    bool_ph4_ref_mol2 = false;
    bool_ph4_ref_txt = false;

    cout << "\nPharmacophore Similarity Score Parameters" << endl;
    cout <<
         "------------------------------------------------------------------------------------------"
         << endl;

    if (!primary_score) {
        tmp = parm.query_param("pharmacophore_score_primary", "no", "yes no");
        if (tmp == "yes")
            use_primary_score = true;
        else
            use_primary_score = false;

        primary_score = use_primary_score;
    }

    if (!secondary_score) {
        //tmp = parm.query_param("pharmacophore_score_secondary", "no", "yes no");
        tmp = "no";
	if (tmp == "yes")
            use_secondary_score = true;
        else
            use_secondary_score = false;

        secondary_score = use_secondary_score;
    }
    if (use_primary_score || use_secondary_score)
        use_score = true;
    else
        use_score = false;
    // input_parameters_main is used here and in descriptor_score. 
    input_parameters_main(parm,"fms_score");
} //end Ph4_Score::input_parameters()

// this is also called in descriptor score // Trent Balius, 2016
// It is beter to used the same input_parameter function for both
// so that there that it is easer to make modifications and so that 
// there are not errors (doing this fix some issues with valgrind). 
void
Ph4_Score::input_parameters_main(Parameter_Reader & parm, string parm_head)
{
    string          tmp;

    if (use_score) {
        // if Pharmacophore is precomputed, read in fms.txt;
        // else ask for the reference mol2 file.
        tmp = parm.query_param(parm_head+"_use_ref_mol2", "no", "yes no");

        if (tmp == "yes")
            bool_ph4_ref_mol2 = true;
        else
            bool_ph4_ref_mol2 = false;

        if (bool_ph4_ref_mol2)
            constant_ph4_ref_file =
                parm.query_param(parm_head+"_ref_mol2_filename", "Ph4.mol2");
        else{
            tmp = parm.query_param(parm_head+"_use_ref_txt", "no", "yes no");

            if (tmp == "yes")
                bool_ph4_ref_txt = true;
            else
                bool_ph4_ref_txt = false;

            if (bool_ph4_ref_txt)
                constant_ph4_ref_file =
                    parm.query_param(parm_head+"_ref_txt_filename", "Ph4.txt");
            else{
                cout << "Error: No pharmacophore reference filename is specified, Program will terminate" << endl;
                exit(0);
            }
        }

        tmp = parm.query_param(parm_head+"_write_reference_pharmacophore_mol2", "no", "yes no");

        if (tmp == "yes")
            bool_out_ref_ph4_mol2 = true;
        else
            bool_out_ref_ph4_mol2 = false;

        tmp = parm.query_param(parm_head+"_write_reference_pharmacophore_txt", "no", "yes no");

        if (tmp == "yes")
            bool_out_ref_ph4_txt = true;
        else
            bool_out_ref_ph4_txt = false;
        
        if (bool_out_ref_ph4_mol2)
            ref_ph4_out_mol2_filename = parm.query_param(parm_head+"_reference_output_mol2_filename", "ref_ph4.mol2");

        if (bool_out_ref_ph4_txt)
            ref_ph4_out_txt_filename = parm.query_param(parm_head+"_reference_output_txt_filename", "ref_ph4.txt");

        tmp = parm.query_param(parm_head+"_write_candidate_pharmacophore", "no", "yes no");

        if (tmp == "yes")
            bool_out_cad_ph4 = true;
        else
            bool_out_cad_ph4 = false;

        if (bool_out_cad_ph4)
            cad_ph4_out_filename = parm.query_param(parm_head+"_candidate_output_filename", "cad_pharmacophore.mol2");

        tmp = parm.query_param(parm_head+"_write_matched_pharmacophore", "no", "yes no");

        if (tmp == "yes")
            bool_out_mat_ph4 = true;
        else
            bool_out_mat_ph4 = false;

        if (bool_out_mat_ph4)
            mat_ph4_out_filename = parm.query_param(parm_head+"_matched_output_filename", "mat_ph4.mol2");
            //mat_ph4 is a subset of cad_ph4


        //to be done: 
        //1. receptor/ligand based matching 
        //2. full match or partial match 
        //3. whether or not to use excluded volume: input receptor mol2/ how to model it?
     
        //ask what matrix to use to compare the ph4.
        //if ref is a ligand based reference, use "o" for "overlap"
        //if ref is a receptor based reference, use "c" for "compatible", to be modified
        ph4_compare_type = parm.query_param(parm_head+"_compare_type", "overlap", "overlap compatible");

        if (ph4_compare_type.compare("overlap") == 0)
            cout << "-------------------------------------------------------------"   << endl
                 << "You are using a ligand based reference for computing the FMS."   << endl
                 << "When the value is 0 then there is a perfect overlap."            << endl
                 << "When the value is negative then you have multi-matched ph4."     << endl
                 << "When the value is positive then you have matches with residual." << endl
                 << "-------------------------------------------------------------"   << endl;
        else if (ph4_compare_type.compare("compatible") == 0) 
            cout << "-------------------------------------------------------------"   << endl
                 << "You are using a receptor based reference for computing the FMS." << endl
                 << "When the value is ? then there is a perfect overlap."            << endl
                 << "When the value is ? then you have multi-matched ph4."            << endl
                 << "When the value is ? then you have matches with residual."        << endl
                 << "This is under development and not currently available."          << endl
                 << "-------------------------------------------------------------"   << endl;

        //use full match or partial match
        tmp = parm.query_param(parm_head+"_full_match", "yes", "yes no");
        if (tmp == "yes")
            ph4_full_match = true;
        if (tmp == "no")
            ph4_full_match = false;

        if (!ph4_full_match){
            tmp = parm.query_param(parm_head+"_partial_match_type", "num_mismatch", "num_mismatch percentage");
            if (tmp == "num_mismatch"){
                ph4_specify_mismatch_num = true;
                cout << "-------------------------------------------------------------" << endl
                     << "You chose to specify the number of mismatch allowed. If this " << endl
                     << "number is no smaller than total number of feature, then all  " << endl
                     << "will be considered success. Recommended num_mismatch: 1-3.   " << endl
                     << "-------------------------------------------------------------" <<endl;
            }
            if (tmp == "percentage") {
                ph4_specify_mismatch_percent = true;
                cout << "-------------------------------------------------------------" << endl
                     << "You chose to specify the percentage of mismatch allowed. If  " << endl
                     << "Recommended percentage of mismatch allowed: 1-20. "            << endl
                     << "-------------------------------------------------------------" << endl;
            }
            if (ph4_specify_mismatch_num) {
                ph4_mismatch_num     = atoi (parm.query_param(parm_head+"_mismatch_num", "1-3").c_str());
                if ((ph4_mismatch_num < 1) || (ph4_mismatch_num > 3)) {
                    cout << "ERROR: parameter must be an integer between 1 to 3. Program will terminate." << endl;
                    exit (0);
                }
            }
 
            if (ph4_specify_mismatch_percent) {
                ph4_mismatch_percent = atoi (parm.query_param(parm_head+"_mismatch_percent", "1-20").c_str());
                if ((ph4_mismatch_percent < 1) || (ph4_mismatch_percent > 20)) {
                    cout << "ERROR: parameter must be an integer between 1 to 20. Program will terminate." << endl;
                    exit (0);
                }
            }
        }

        //input the constant parameter k, which is the weight on the match rate term in FMS
        ph4_rate_k = atoi (parm.query_param(parm_head+"_match_rate_weight", "5.0").c_str());
        if (ph4_rate_k <= 0.0) {
            cout << "ERROR: Parameter must be a double greater than zero. Program will terminate."
                 << endl;
            exit(0);
        }

        //input the default distance cutoff r in the FMS matching protocol
        ph4_dist_r = atoi (parm.query_param(parm_head+"_match_dist_cutoff", "1.0").c_str());
        if (ph4_dist_r <= 0) {
            cout << "ERROR: Parameter must be a float greater than zero. Program will terminate."
                 << endl;
            exit(0);
        }

        //input the default scalar projection cutoff cos(sigma) in the FMS matching protocol
        //v_proj cutoff, 60 deg:0.5, 45 deg:0.7071, 30 deg:0.8660.
        ph4_proj_cos = atoi (parm.query_param(parm_head+"_match_proj_cutoff", "0.7071").c_str());
        if (ph4_proj_cos <0 or ph4_proj_cos >1) {
            cout << "ERROR: Parameter must be a double between zero and one. Program will terminate."
                 << endl;
            exit(0);
        }

        //input the default FMS score for pharmacophore model pairs with no matches
        ph4_max_x = atoi (parm.query_param(parm_head+"_max_score", "20").c_str());
        if (ph4_max_x <= 0) {
            cout << "ERROR: Parameter must be an integer greater than zero. Program will terminate."
                 << endl;
            exit(0);
        }

/*
        //use excluded volume from a receptor file or not
        tmp = parm.query_param("ph4_score_use_excluded_volume", "no", "yes no");
        if (tmp == "yes")
            bool_ph4_use_excluded_volume = true;
        if (tmp == "no")
            bool_ph4_use_excluded_volume = false;

        if (bool_ph4_use_excluded_volume)
            receptor_filename = parm.query_param("ph4_score_rec_filename", "receptor.mol2");
*/  
    }
} //end Ph4_Score::input_parameters_main()


// +++++++++++++++++++++++++++++++++++++++++
// This functions reads in the reference based on the DOCKMol ref file provided by input
void 
Ph4_Score::initialize_read_Mol2_helper(AMBER_TYPER & typer){
    ifstream ref_file;
    ref_file.open(constant_ph4_ref_file.c_str());

    if (ref_file.fail()) {
        cout << "Error Opening Reference File!" << endl;
        exit(0);
    }

    if (!Read_Mol2(reference, ref_file, false, false, false)) {
        cout << "Error Reading Reference Molecule!" << endl;
        exit(0);
    }

    ref_file.close();
}

// +++++++++++++++++++++++++++++++++++++++++
// This function reads in the reference based on a Ph4 file provided by input
void
Ph4_Score::initialize_read_ph4_helper(AMBER_TYPER & typer){
    ifstream ref_file;
    ref_file.open(constant_ph4_ref_file.c_str());
    if (ref_file.fail()) {
        cout << "Error Opening Reference File!" << endl;
        exit(0);
    }

    if (!read_ph4_txt(Ph4_ref, ref_file)) {
                     //make ph4 ref from ref txt
                     cout << "Error reading in Reference ph4!" << endl;
                     exit(0);
    }

    ref_file.close();
    // read in ph4 ref from txt and set ph4
}


// +++++++++++++++++++++++++++++++++++++++++
// Write out a ph4 struct or mol2 from the other
void
Ph4_Score::initialize_output_helper(){
    if(bool_out_ref_ph4_mol2){
       refresh_file(ref_ph4_out_mol2_filename);
       submit_ph4_mol2(Ph4_ref, ref_ph4_out_mol2_filename, reference.title, false);
    }
    if(bool_out_ref_ph4_txt){
       refresh_file(ref_ph4_out_txt_filename);
       submit_ph4_txt(Ph4_ref, ref_ph4_out_txt_filename);
    }

    if(bool_out_cad_ph4)
       refresh_file(cad_ph4_out_filename);

    if(bool_out_mat_ph4)
       refresh_file(mat_ph4_out_filename);
}

// +++++++++++++++++++++++++++++++++++++++++
//
void
Ph4_Score::initialize_prepare_Ph4_components_helper(AMBER_TYPER & typer){
    int att_exp = 6, rep_exp = 9;

    // Hardcoded for now.  Read from the ph4.defn file in the future;
    // also change the number for the definition of component in Ph4_Score.
    //char*  list_types[] = {"hydrophobic", "donor", "acceptor", "aromatic",  "aroAcc", "ring"};
    //char*  list_types[] = {"aromatic","positive", "negative","ring"};
    char const *const list_types[] = {"hydrophobic", "donor", "acceptor",
                          "aromatic", "aroAcc", "positive", "negative", "ring"};
    int const num_types = sizeof list_types / sizeof list_types[0];

    list_comp.clear();// initialize the list of components.
    for (int i=0;i<num_types;i++) {
         list_comp.push_back(list_types[i]);
    }

    match_comp = new int [list_comp.size()];
    components = new double [list_comp.size()];
    match_num  = new int [3];
    match_term = new double [2];
    //we also initialize the vdw energy here


    if (use_score || true){
        init_vdw_energy(typer, att_exp, rep_exp);
    }

}

// +++++++++++++++++++++++++++++++++++++++++
// Takes a reference DOCKMol and prepares it for use in pharmacohphore scoring
void 
Ph4_Score::initialize_prepare_Mol2_helper(AMBER_TYPER & typer, DOCKMol & reference){
    bool read_vdw = true;
    bool use_chem = false;
    bool use_ph4  = true;
    bool use_volume = false;
    typer.prepare_molecule(reference, read_vdw, use_chem, use_ph4, use_volume);

    compute_ph4(reference,Ph4_ref);
    // make ph4 ref from ref mol2

}

// +++++++++++++++++++++++++++++++++++++++++
// subroutine designed to initialize important values for score calculations when the parameters are not read in from input (i.e. in dmax)
void
Ph4_Score::initialize_define_score_coeffecients_helper(){
    this->ph4_rate_k = 5.0; //hard coded default values
    this->ph4_dist_r = 1.0;
    this->ph4_proj_cos = 0.7071;
    this->ph4_max_x = 20;
}

// +++++++++++++++++++++++++++++++++++++++++
// An overloaded version of initialize meant to be more lightweight and use any DOCKMol object as a reference without having to worry about IO (i.e. for dmax)
void
Ph4_Score::initialize(AMBER_TYPER & typer, DOCKMol & ref_mol){
	Ph4_Score::use_score = true;
	initialize_define_score_coeffecients_helper();
	initialize_prepare_Ph4_components_helper(typer);
	this->reference = ref_mol;
	initialize_prepare_Mol2_helper(typer, ref_mol);
}

// +++++++++++++++++++++++++++++++++++++++++
// The primary initialize function used in routine DOCK6 FMS function, dependent on IO from input
void
Ph4_Score::initialize(AMBER_TYPER & typer)
{
    ifstream        ref_file;
//    ifstream        rec_file; // in case we are doing receptor based scoring
     
    initialize_prepare_Ph4_components_helper(typer); //prepare the definitions for different pharmacophore component, init vdw energy term 
    if (use_score) {
        if (bool_ph4_ref_mol2){
            initialize_read_Mol2_helper(typer); 
	    initialize_prepare_Mol2_helper(typer, reference);
	    initialize_output_helper();
        }
        else if (bool_ph4_ref_txt){
	    initialize_read_ph4_helper(typer);
	    initialize_output_helper();
	}
        else {
                cout << "Initialize Error: No reference specified!" << endl;
                exit(0);
        }
/*
        if (bool_ph4_use_excluded_volume){
            rec_file.open(ph4_score_rec_filename.c_str());
            
            if (rec_file.fail() {
                cout << "Error Opening Receptor File!" << endl;
                exit(0);
            }

            if (!Read_Mol2(receptor, rec_file, false, false, false)) {
                cout << "Error Reading Receptor Molecule!" << endl;
                exit(0);
            }

            rec_file.close();

            read_vdw = true;
            use_chem = false;
            use_ph4  = true;
            use_volume = false;
            typer.prepare_molecule(receptor, read_vdw, use_chem, use_ph4, use_volume);                  
        }
*/
    }
    return;
}

// +++++++++++++++++++++++++++++++++++++++++
void
Ph4_Score::refresh_file(string filename)
{
    ostringstream file;
    file << filename;
    fstream fout;
    fout.open (file.str().c_str(), fstream::out);

    fout.close();
    return;
}
// +++++++++++++++++++++++++++++++++++++++++
void
Ph4_Score::submit_ph4_mol2(Ph4_Struct & ph4_out, string filename, string molname, bool matched_ph4_only)
{
//write out a ph4.mol2 file for visualization.
//now it will write out ph4.mol2 for both the reference and the candidate molecules.
    ostringstream file;
    file << filename;
    fstream fout;
    fout.open (file.str().c_str(), fstream::out | fstream::app);

    int num_ph4=1, num_ph4_atom=0;
    char line[200];
    bool print_all_ph4;
    bool print_current_ph4;

    if(matched_ph4_only)
       print_all_ph4 = false;
    else
       print_all_ph4 = true;

    for (int i=0; i< ph4_out.num_features; i++) {

         if (print_all_ph4)
            print_current_ph4=true;
         else if (ph4_out.contri_to_score[i])
            print_current_ph4=true;
         else
            print_current_ph4=false;

         if (print_current_ph4){
             if(ph4_out.ph4_types[i].compare("aromatic")==0||ph4_out.ph4_types[i].compare("ring")==0){
                num_ph4_atom+=3;//two dummy atom per ring
             }
             else if((ph4_out.ph4_types[i].compare("donor")==0)||(ph4_out.ph4_types[i].compare("acceptor")==0)){
                num_ph4_atom+=2;//one dummy atom per HBD/HBA
             }
             else{
                num_ph4_atom++;
             }
         }
    }
          
    fout << "@<TRIPOS>MOLECULE" << endl;
    fout << molname << endl;
    
    sprintf(line, " %d 0 1 0 0", num_ph4_atom);
    fout << line << endl;

    fout << "Pharmcophore" << endl;
    fout << "NO_CHARGES"   << endl << endl << endl;
    
    fout << "@<TRIPOS>ATOM" << endl;

    for (int i=0; i< ph4_out.num_features; i++) {

         if (print_all_ph4)
            print_current_ph4=true;
         else if (ph4_out.contri_to_score[i])
            print_current_ph4=true;
         else
            print_current_ph4=false;

         if (print_current_ph4){
             if(ph4_out.ph4_types[i].compare("aromatic")==0){
                sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", num_ph4, " S", ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], " S", "1", "SPH", "0.0000");
                fout << line  << endl;
                num_ph4++;
                sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", num_ph4, " H1", ph4_out.x[i]+ph4_out.v_x[i], ph4_out.y[i]+ph4_out.v_y[i], ph4_out.z[i]+ph4_out.v_z[i], " H", "1", "SPH", "0.0000");
                fout << line  << endl;
                num_ph4++;
                sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", num_ph4, " H2", ph4_out.x[i]-ph4_out.v_x[i], ph4_out.y[i]-ph4_out.v_y[i], ph4_out.z[i]-ph4_out.v_z[i], " H", "1", "SPH", "0.0000");
                fout << line  << endl;
                num_ph4++;
             }

             else if(ph4_out.ph4_types[i].compare("ring")==0){
                sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", num_ph4, " P", ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], " S", "1", "SPH", "0.0000");
                fout << line  << endl;
                num_ph4++;
                sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", num_ph4, " H3", ph4_out.x[i]+ph4_out.v_x[i], ph4_out.y[i]+ph4_out.v_y[i], ph4_out.z[i]+ph4_out.v_z[i], " H", "1", "SPH", "0.0000");
                fout << line  << endl;
                num_ph4++;
                sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", num_ph4, " H4", ph4_out.x[i]-ph4_out.v_x[i], ph4_out.y[i]-ph4_out.v_y[i], ph4_out.z[i]-ph4_out.v_z[i], " H", "1", "SPH", "0.0000");
                fout << line  << endl;
                num_ph4++;
             }

             else if(ph4_out.ph4_types[i].compare("donor")==0){
                sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", num_ph4, " HD", ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], " H", "1", "SPH", "0.0000");
                fout << line  << endl;
                num_ph4++;
                sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", num_ph4, " N", ph4_out.x[i]-ph4_out.v_x[i], ph4_out.y[i]-ph4_out.v_y[i], ph4_out.z[i]-ph4_out.v_z[i], " N", "1", "SPH", "0.0000");
                fout << line  << endl;
                num_ph4++;
             }
             else if(ph4_out.ph4_types[i].compare("acceptor")==0){
                sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", num_ph4, " O", ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], " O", "1", "SPH", "0.0000");
                fout << line  << endl;
                num_ph4++;
                sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", num_ph4, " HA", ph4_out.x[i]+ph4_out.v_x[i], ph4_out.y[i]+ph4_out.v_y[i], ph4_out.z[i]+ph4_out.v_z[i], " H", "1", "SPH", "0.0000");
                fout << line  << endl;
                num_ph4++;
             }
             else if(ph4_out.ph4_types[i].compare("hydrophobic")==0){
                sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", num_ph4, " C", ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], " C", "1", "SPH", "0.0000");
                fout << line  << endl;
                num_ph4++;
             }
             else if(ph4_out.ph4_types[i].compare("positive")==0){
                sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%3s%7s%5s%15s", num_ph4, " Na", ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], " Na", "1", "SPH", "1.0000");
                fout << line  << endl;
                num_ph4++;
             }
             else if(ph4_out.ph4_types[i].compare("negative")==0){
                sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%3s%7s%5s%15s", num_ph4, " Cl", ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], " Cl", "1", "SPH", "-1.0000");
                fout << line  << endl;
                num_ph4++;
             }
             else if(ph4_out.ph4_types[i].compare("polar")==0){// shouldn't happen as we changed the set up of polar for ph4_mol
                sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", num_ph4, " OP", ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], " O", "1", "SPH", "0.0000");
                fout << line  << endl;
                num_ph4++;
             }
         }
    }
    

    fout << "@<TRIPOS>BOND" << endl;
    fout << "@<TRIPOS>SUBSTRUCTURE" << endl;
    fout << " 1   ****   1" << endl;
    if (num_ph4-1 != num_ph4_atom){
    fout << "ERROR: num_ph4 and num_ph4_atom not consistent!" << endl;
    }
    fout << endl;
    fout.close();

    return;
}
// +++++++++++++++++++++++++++++++++++++++++
void
Ph4_Score::submit_ph4_txt(Ph4_Struct & ph4_out, string filename)
{
//write out a ph4_reference.txt for further refinement
    ostringstream file;
    file << filename;
    fstream fout;
    fout.open (file.str().c_str(), fstream::out | fstream::app);

    int num_ph4=1;
    char line[200];

    fout << "# ################################ #" << endl;
    fout << "# ### Reference Pharmacophore  ### #" << endl;
    fout << "# ################################ #" << endl;

    fout << "ph4name ph4id   coox       cooy     cooz       vecx      vecy      vecz      radius" << endl;

    for (int i=0; i< ph4_out.num_features; i++) {
        if(ph4_out.ph4_types[i].compare("aromatic")==0){
           sprintf(line, "%-7s%6d%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f", "ARO", num_ph4, ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], ph4_out.v_x[i], ph4_out.v_y[i], ph4_out.v_z[i], ph4_out.radius[i] );
                fout << line  << endl;
           num_ph4++;
        }
        else if(ph4_out.ph4_types[i].compare("ring")==0){
           sprintf(line, "%-7s%6d%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f", "RNG", num_ph4, ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], ph4_out.v_x[i], ph4_out.v_y[i], ph4_out.v_z[i], ph4_out.radius[i] );
                fout << line  << endl;
           num_ph4++;
        }
        else if(ph4_out.ph4_types[i].compare("donor")==0){
           sprintf(line, "%-7s%6d%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f", "HBD", num_ph4, ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], ph4_out.v_x[i], ph4_out.v_y[i], ph4_out.v_z[i], ph4_out.radius[i] );
                fout << line  << endl;
           num_ph4++;
        }
        else if(ph4_out.ph4_types[i].compare("acceptor")==0){
           sprintf(line, "%-7s%6d%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f", "HBA", num_ph4, ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], ph4_out.v_x[i], ph4_out.v_y[i], ph4_out.v_z[i], ph4_out.radius[i] );
                fout << line  << endl;
           num_ph4++;
        }
        else if(ph4_out.ph4_types[i].compare("hydrophobic")==0){
           sprintf(line, "%-7s%6d%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f", "PHO", num_ph4, ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], ph4_out.v_x[i], ph4_out.v_y[i], ph4_out.v_z[i], ph4_out.radius[i] );
                fout << line  << endl;
           num_ph4++;
        }
        else if(ph4_out.ph4_types[i].compare("positive")==0){
           sprintf(line, "%-7s%6d%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f", "POS", num_ph4, ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], ph4_out.v_x[i], ph4_out.v_y[i], ph4_out.v_z[i], ph4_out.radius[i] );
                fout << line  << endl;
           num_ph4++;
        }
        else if(ph4_out.ph4_types[i].compare("negative")==0){
           sprintf(line, "%-7s%6d%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f", "NEG", num_ph4, ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], ph4_out.v_x[i], ph4_out.v_y[i], ph4_out.v_z[i], ph4_out.radius[i] );
                fout << line  << endl;
           num_ph4++;
        }
        else if(ph4_out.ph4_types[i].compare("polar")==0){// shouldn't happen as we changed the set up of polar for ph4_mol
           sprintf(line, "%-7s%6d%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f", "POL", num_ph4, ph4_out.x[i], ph4_out.y[i], ph4_out.z[i], ph4_out.v_x[i], ph4_out.v_y[i], ph4_out.v_z[i], ph4_out.radius[i] );
                fout << line  << endl;
           num_ph4++;
        }
    }

    if (num_ph4-1 != ph4_out.num_features){
    fout << "ERROR: num_ph4 and num_ph4_atom not consistent!" << endl;
    }
    fout << endl;
    fout.close();

    return;
}
// +++++++++++++++++++++++++++++++++++++++++
bool
Ph4_Score::compute_ph4(DOCKMol & mol, Ph4_Struct & ph4_mol)
{//calculate ph4 of mol with compute_ph4
        Ph4_Struct atom_mol;//declare the colored atom set first, as a type Ph4_Struct variable

        atom_mol.num_features = 0;
        ph4_mol.num_features  = 0;  
      
        for ( int i=0; i<mol.num_atoms; i++) {
             if (mol.ph4_types_assigned){
                 if (mol.ph4_types[i] != "null")
                     atom_mol.num_features++;
             }
        }
        
        atom_mol.initialize();
        //compute colored atom set atom_mol from ligand pose mol.
        // cout xyz and color, and the id in atom_mol record the atom number in the mol file.
        //*cout << "\n\n";
        //*cout << "colored atom set: " << endl;

        int count = 0;
        for ( int i=0; i<mol.num_atoms; i++) {
             //*cout << i + 1 << ",  "
                  //*<< mol.atom_names[i] << ",  "
                  //*<< mol.x[i] << ",  "
                  //*<< mol.y[i] << ",  "
                  //*<< mol.z[i] << ",  "
                  //*<< mol.ph4_types[i]
                  //*<< endl;
             //output the list of colored atom set.
             if   (mol.ph4_types_assigned){
                   if (mol.ph4_types[i] != "null"){
                        atom_mol.id[count] = i;
                        atom_mol.x[count] = mol.x[i];
                        atom_mol.y[count] = mol.y[i];
                        atom_mol.z[count] = mol.z[i];
                        atom_mol.ph4_types[count] = mol.ph4_types[i];
                        count++;
                   }
             }//define colored ph4 atom with atoms in mol with ph4_type!="null"
        }
        
        if (count != atom_mol.num_features) {
            //*cout << "warning: number of features in atom_mol don't match size of atom_mol." << endl;
            //*cout << "reset atom_mol.number_features" << endl;
            atom_mol.num_features = count;
        }

        ph4_mol.num_features = 4*atom_mol.num_features;
        //hardcoded, depending on how ph4 feature are defined vectors/dummy atoms
        //for the definition now, an atom will correspond to 3 dummy atoms at most
        ph4_mol.initialize();

 
        //*cout << "\n\n";
        //*cout << "ph4 point set: " << endl;
        count = 0;
        //float r_ring, r_vec, r_temp, r=1.0; //default radius of ph4 point r, to be modified as a parameter.
        float r_ring, r_vec, r_temp, r=ph4_dist_r; //default radius of ph4 point r        

        float * crd;
        float * v_crd, * va_crd, * vb_crd, * vc_crd;
 
        vector <int> nbrs; //neighor atom list, for donot/acceptor features
        
        crd   = new float [3];
        v_crd = new float [3];
        va_crd= new float [3];
        vb_crd= new float [3];
        vc_crd= new float [3];

        vector <int> visited, ring_member; 
        // a list of atoms in any ring that is pre-identified the mol file; and list of current ring members

        visited.clear();
        bool  find_ring, aro_ring;
        int   num_ring = 0;

        for (int i=0; i<atom_mol.num_features; i++) {
        //find all rings and then determine if they are aromatic (planar)
             find_ring = true;
             for(int j=0; j < visited.size(); j++) {
                 if (atom_mol.id[i] == visited[j]) {
                     find_ring = false;
                 }
             }
             if (find_ring){
                 ring_member.clear();
                 if(find_ring_member(mol, atom_mol.id[i], visited, ring_member)){
                    num_ring++;
                    crd[0]= 0.0;
                    crd[1]= 0.0;
                    crd[2]= 0.0;
                    //*cout << "size of the " << num_ring << "th ring: " << ring_member.size() << endl;
                    //*cout << "atom list in the ring: ";
                    for (int j=0; j< ring_member.size(); j++){
                         //*cout << ring_member[j]+1 << " ";
                         crd[0]=crd[0]+mol.x[ring_member[j]];
                         crd[1]=crd[1]+mol.y[ring_member[j]];
                         crd[2]=crd[2]+mol.z[ring_member[j]];
                    }
                    //*cout << endl;
                    //*cout << "center of the ring";
                    crd[0] = crd [0]/(float)ring_member.size(); //*cout << crd[0] << ", " ;
                    crd[1] = crd [1]/(float)ring_member.size(); //*cout << crd[1] << ", " ;
                    crd[2] = crd [2]/(float)ring_member.size(); //*cout << crd[2] << ", " << endl;
                    //compute the center of ring  coordinate crd[]
                    v_crd[0]  = 0;
                    v_crd[1]  = 0;
                    v_crd[2]  = 0;
                    //v_crd is the cross product (norm vector) of adjacent center-vertex vectors
                    for (int j=0; j < ring_member.size()-1; j++) {
                          va_crd[0] = mol.x[ring_member[j]]-crd[0];
                          vb_crd[0] = mol.x[ring_member[j+1]]-crd[0];
                          va_crd[1] = mol.y[ring_member[j]]-crd[1];
                          vb_crd[1] = mol.y[ring_member[j+1]]-crd[1];
                          va_crd[2] = mol.z[ring_member[j]]-crd[2];
                          vb_crd[2] = mol.z[ring_member[j+1]]-crd[2];
                          //va_crd : vector from center of ring to current ring_member
                          //vb_crd : vector from center of ring to the next ring_member
                          vc_crd[0] =(va_crd[1]*vb_crd[2]) - (va_crd[2]*vb_crd[1]);
                          vc_crd[1] =(va_crd[2]*vb_crd[0]) - (va_crd[0]*vb_crd[2]);
                          vc_crd[2] =(va_crd[0]*vb_crd[1]) - (va_crd[1]*vb_crd[0]);
                          //vc_crd : norm vector of adjacent center-vertex vectors;
                          v_crd[0] += vc_crd[0];
                          v_crd[1] += vc_crd[1];
                          v_crd[2] += vc_crd[2];
                    }
                    va_crd[0] = mol.x[ring_member[ring_member.size()-1]]-crd[0];
                    vb_crd[0] = mol.x[ring_member[0]]-crd[0];
                    va_crd[1] = mol.y[ring_member[ring_member.size()-1]]-crd[1];
                    vb_crd[1] = mol.y[ring_member[0]]-crd[1];
                    va_crd[2] = mol.z[ring_member[ring_member.size()-1]]-crd[2];
                    vb_crd[2] = mol.z[ring_member[0]]-crd[2];
                    vc_crd[0] =(va_crd[1]*vb_crd[2]) - (va_crd[2]*vb_crd[1]);
                    vc_crd[1] =(va_crd[2]*vb_crd[0]) - (va_crd[0]*vb_crd[2]);
                    vc_crd[2] =(va_crd[0]*vb_crd[1]) - (va_crd[1]*vb_crd[0]);
                    //compute vi_crd of  the last ring_member and the first ring_member
                    v_crd[0] += vc_crd[0];
                    v_crd[1] += vc_crd[1];
                    v_crd[2] += vc_crd[2];
                    v_crd[0]  = v_crd[0]/(float)ring_member.size();
                    v_crd[1]  = v_crd[1]/(float)ring_member.size();
                    v_crd[2]  = v_crd[2]/(float)ring_member.size();
                    //average norm vector. to be modified to combine with the aromatic identification
                    //define v_crd = [[]] as two dimensional array.
                    r_vec     = (v_crd[0]*v_crd[0]) + (v_crd[1]*v_crd[1]) + (v_crd[2]*v_crd[2]);
                    r_vec     = sqrt(r_vec);
                    v_crd[0]  = v_crd[0]/r_vec;
                    v_crd[1]  = v_crd[1]/r_vec;
                    v_crd[2]  = v_crd[2]/r_vec;
                    //normalize v_crd[]
                    r_ring = 0;
                    aro_ring = false;
                    if (ring_member.size() != 3){
                        aro_ring = true;//3 membered ring is not aromatic but planar
                    }
                    for (int j=0; j < ring_member.size()-1; j++) {
                         if (ring_member.size() != 3){
                             va_crd[0] = mol.x[ring_member[j]]-crd[0];
                             vb_crd[0] = mol.x[ring_member[j+1]]-crd[0];
                             va_crd[1] = mol.y[ring_member[j]]-crd[1];
                             vb_crd[1] = mol.y[ring_member[j+1]]-crd[1];
                             va_crd[2] = mol.z[ring_member[j]]-crd[2];
                             vb_crd[2] = mol.z[ring_member[j+1]]-crd[2];
                             vc_crd[0] =(va_crd[1]*vb_crd[2]) - (va_crd[2]*vb_crd[1]);
                             vc_crd[1] =(va_crd[2]*vb_crd[0]) - (va_crd[0]*vb_crd[2]);
                             vc_crd[2] =(va_crd[0]*vb_crd[1]) - (va_crd[1]*vb_crd[0]);
                             r_temp    = (vc_crd[0]*vc_crd[0]) + (vc_crd[1]*vc_crd[1]) + (vc_crd[2]*vc_crd[2]);
                             r_temp    = sqrt(r_temp);
                             if(fabs((v_crd[0]*vc_crd[0]+v_crd[1]*vc_crd[1]+v_crd[2]*vc_crd[2])/r_temp) < 0.99){//inner product of the avg norm vector to each individual one, need to be close to 1 for aromatic
                                aro_ring = false;
                             }
                         }
                         r_ring += sqrt((mol.x[ring_member[j]]-crd[0]*mol.x[ring_member[j]]-crd[0]) + (mol.y[ring_member[j]]-crd[1]*mol.y[ring_member[j]]-crd[1]) + (mol.z[ring_member[j]]-crd[2]*mol.z[ring_member[j]]-crd[2]));
                         //cout << "radius of the ring now " << r_ring << endl;
                    }
                    if (ring_member.size() != 3){
                        va_crd[0] = mol.x[ring_member[ring_member.size()-1]]-crd[0];
                        vb_crd[0] = mol.x[ring_member[0]]-crd[0];
                        //adajent atom is the first atom, instead of the next one
                        va_crd[1] = mol.y[ring_member[ring_member.size()-1]]-crd[1];
                        vb_crd[1] = mol.y[ring_member[0]]-crd[1];
                        va_crd[2] = mol.z[ring_member[ring_member.size()-1]]-crd[2];
                        vb_crd[2] = mol.z[ring_member[0]]-crd[2];
                        vc_crd[0] =(va_crd[1]*vb_crd[2]) - (va_crd[2]*vb_crd[1]);
                        vc_crd[1] =(va_crd[2]*vb_crd[0]) - (va_crd[0]*vb_crd[2]);
                        vc_crd[2] =(va_crd[0]*vb_crd[1]) - (va_crd[1]*vb_crd[0]);
                        r_temp    = (vc_crd[0]*vc_crd[0]) + (vc_crd[1]*vc_crd[1]) + (vc_crd[2]*vc_crd[2]);
                        r_temp    = sqrt(r_temp);
                        if(fabs((v_crd[0]*vc_crd[0]+v_crd[1]*vc_crd[1]+v_crd[2]*vc_crd[2])/r_temp) < 0.99){
                           aro_ring = false;
                        }
                    }
                    r_ring += sqrt(((mol.x[ring_member[ring_member.size()-1]]-crd[0])*(mol.x[ring_member[ring_member.size()-1]]-crd[0])) + ((mol.y[ring_member[ring_member.size()-1]]-crd[1])*(mol.y[ring_member[ring_member.size()-1]]-crd[1])) + ((mol.z[ring_member[ring_member.size()-1]]-crd[2])*(mol.z[ring_member[ring_member.size()-1]]-crd[2])));
                    //cout << "radius of the ring now " << r_ring << endl;
                    r_ring = r_ring/(float)ring_member.size();
                    //cout << "radius of the ring now " << r_ring << endl;
                    //*cout << "first atom in the ring: " << atom_mol.id[i]+1 << endl << endl;
                    ph4_mol.id[count]= count + 1;        //*cout << ph4_mol.id[count] << ", " << atom_mol.id[i] +1 << ", ";
                    ph4_mol.x[count] = crd[0];           //*cout << ph4_mol.x[count]  << ", ";
                    ph4_mol.y[count] = crd[1];           //*cout << ph4_mol.y[count]  << ", ";
                    ph4_mol.z[count] = crd[2];           //*cout << ph4_mol.z[count]  << ", ";
                    ph4_mol.v_x[count] = v_crd[0];       //*cout << ph4_mol.v_x[count]<< ", ";
                    ph4_mol.v_y[count] = v_crd[1];       //*cout << ph4_mol.v_y[count]<< ", ";
                    ph4_mol.v_z[count] = v_crd[2];       //*cout << ph4_mol.v_z[count]<< ", ";
                    ph4_mol.radius[count] = r_ring;      //*cout << ph4_mol.radius[count]  << ", ";
                    if(aro_ring){
                       ph4_mol.ph4_types[count] = "aromatic";//hardcoded!
                    }
                    else{
                       ph4_mol.ph4_types[count] = "ring";
                    }
                    //*cout << ph4_mol.ph4_types[count]  << endl;
                    //*cout << endl;
                    count++;
                 }
             }

        }

        for (int i=0; i<atom_mol.num_features; i++) {
        //compute ph4 set (features other than rings)  from the colored atom set.
             if (atom_mol.ph4_types[i].compare("hydrophobic") == 0){
                 //to be modified so that hydrophobic atoms can be grouped/only ones that is not in the ring/not used;
                 
                 bool use_hydro = true;
                 //only consider hydrophobic points that is not in a ring. or turn hydrophobic on/off
                 
                 for(int j=0; j < visited.size(); j++) {
                     if (atom_mol.id[i] == visited[j]) {
                         use_hydro = false;
                    }
                 }
                 
                 if (use_hydro == true){
                     ph4_mol.id[count]= count + 1;        //*cout << ph4_mol.id[count] << ", " << atom_mol.id[i]+1 << ", ";
                     ph4_mol.x[count] = atom_mol.x[i];    //*cout << ph4_mol.x[count]  << ", ";
                     ph4_mol.y[count] = atom_mol.y[i];    //*cout << ph4_mol.y[count]  << ", ";
                     ph4_mol.z[count] = atom_mol.z[i];    //*cout << ph4_mol.z[count]  << ", ";
                     ph4_mol.v_x[count] = 1;              //*cout << ph4_mol.v_x[count]<< ", ";
                     ph4_mol.v_y[count] = 0;              //*cout << ph4_mol.v_y[count]<< ", ";
                     ph4_mol.v_z[count] = 0;              //*cout << ph4_mol.v_z[count]<< ", ";
                     //default directory of ph4 is (1,0,0) 
                     ph4_mol.radius[count] = r;           //*cout << ph4_mol.radius[count]  << ", ";
                     ph4_mol.ph4_types[count] = atom_mol.ph4_types[i];
                     //*cout << ph4_mol.ph4_types[count]  << endl;
                     count++;
                 }
             }

             else if (atom_mol.ph4_types[i].compare("donor") == 0){
                 nbrs.clear();
                 for (int j=0; j < mol.num_bonds; j++){
                      if (mol.bonds_origin_atom[j] == atom_mol.id[i]) {
                              nbrs.push_back(mol.bonds_target_atom[j]);
                      }
                      else if (mol.bonds_target_atom[j] == atom_mol.id[i]) {
                          nbrs.push_back(mol.bonds_origin_atom[j]);
                     }
                 }
                 for(int j=0; j <nbrs.size(); j++) {
                     //create a ph4 feature for each atom that is connected to a donor (H)
                     ph4_mol.id[count]= count + 1;        //*cout << ph4_mol.id[count] << ", " << atom_mol.id[i]+1 << ", ";
                     ph4_mol.x[count] = atom_mol.x[i];    //*cout << ph4_mol.x[count]  << ", ";
                     ph4_mol.y[count] = atom_mol.y[i];    //*cout << ph4_mol.y[count]  << ", ";
                     ph4_mol.z[count] = atom_mol.z[i];    //*cout << ph4_mol.z[count]  << ", ";
                     v_crd[0] = atom_mol.x[i]-mol.x[nbrs[j]];
                     v_crd[1] = atom_mol.y[i]-mol.y[nbrs[j]];
                     v_crd[2] = atom_mol.z[i]-mol.z[nbrs[j]];
                     //direction of the donor feature : from heavy atom to hydrogen
                     r_vec  = (v_crd[0]*v_crd[0]) + (v_crd[1]*v_crd[1]) + (v_crd[2]*v_crd[2]);
                     r_vec  = sqrt(r_vec);
                     ph4_mol.v_x[count] = v_crd[0]/r_vec; //*cout << ph4_mol.v_x[count]<< ", ";
                     ph4_mol.v_y[count] = v_crd[1]/r_vec; //*cout << ph4_mol.v_y[count]<< ", ";
                     ph4_mol.v_z[count] = v_crd[2]/r_vec; //*cout << ph4_mol.v_z[count]<< ", ";
                     ph4_mol.radius[count] = r;           //*cout << ph4_mol.radius[count]  << ", ";
                     ph4_mol.ph4_types[count] = atom_mol.ph4_types[i];
                     //*cout << ph4_mol.ph4_types[count]  << endl;
                     count++;
                 }
             }

             else if (atom_mol.ph4_types[i].compare("acceptor") == 0){
                 nbrs.clear();
                 for (int j=0; j < mol.num_bonds; j++){
                      if (mol.bonds_origin_atom[j] == atom_mol.id[i]) {
                          nbrs.push_back(mol.bonds_target_atom[j]);
                      }
                      else if (mol.bonds_target_atom[j] == atom_mol.id[i]) {
                          nbrs.push_back(mol.bonds_origin_atom[j]);
                      }
                 }
                 ph4_mol.id[count]= count + 1;        //*cout << ph4_mol.id[count] << ", " << atom_mol.id[i]+1 << ", ";
                 ph4_mol.x[count] = atom_mol.x[i];    //*cout << ph4_mol.x[count]  << ", ";
                 ph4_mol.y[count] = atom_mol.y[i];    //*cout << ph4_mol.y[count]  << ", ";
                 ph4_mol.z[count] = atom_mol.z[i];    //*cout << ph4_mol.z[count]  << ", ";

                 bool use_acc_dir = true; //turn the direction for acceptor on/off.
    
                 if (use_acc_dir) {
                     v_crd[0] = 0;
                     v_crd[1] = 0;
                     v_crd[2] = 0;
                     for(int j=0; j <nbrs.size(); j++) {
                        //the HB vector is the sum of all vectors from the neighbor atoms to the acceptor
                         v_crd[0] += (atom_mol.x[i]-mol.x[nbrs[j]])/(float(nbrs.size()));
                         v_crd[1] += (atom_mol.y[i]-mol.y[nbrs[j]])/(float(nbrs.size()));
                         v_crd[2] += (atom_mol.z[i]-mol.z[nbrs[j]])/(float(nbrs.size()));
                     }
                     r_vec  = (v_crd[0]*v_crd[0]) + (v_crd[1]*v_crd[1]) + (v_crd[2]*v_crd[2]);
                     r_vec  = sqrt(r_vec);
                     ph4_mol.v_x[count] = v_crd[0]/r_vec; //*cout << ph4_mol.v_x[count]<< ", ";
                     ph4_mol.v_y[count] = v_crd[1]/r_vec; //*cout << ph4_mol.v_y[count]<< ", ";
                     ph4_mol.v_z[count] = v_crd[2]/r_vec; //*cout << ph4_mol.v_z[count]<< ", ";
                 }
                 else {
                     ph4_mol.v_x[count] = 1;              //*cout << ph4_mol.v_x[count]<< ", ";
                     ph4_mol.v_y[count] = 0;              //*cout << ph4_mol.v_y[count]<< ", ";
                     ph4_mol.v_z[count] = 0;              //*cout << ph4_mol.v_z[count]<< ", ";
                     //default directory of ph4 is (1,0,0)
                 }
                 
                 ph4_mol.radius[count] = r;           //*cout << ph4_mol.radius[count]  << ", ";
                 ph4_mol.ph4_types[count] = atom_mol.ph4_types[i];
                 //*cout << ph4_mol.ph4_types[count]  << endl;
                 count++;
             }
/*
             else if (atom_mol.ph4_types[i].compare("aromatic") == 0){
                 find_ring = true;
                 for(int j=0; j < visited.size(); j++) {
                     if (atom_mol.id[i] == visited[j]) {
                         find_ring = false;
                     }
                 }
                 if (find_ring){
                     ring_member.clear();
                     find_ring = find_ring_member(mol, atom_mol.id[i], visited, ring_member);
                     crd[0]= 0.0;
                     crd[1]= 0.0;
                     crd[2]= 0.0;
                     cout << "size of the ring: " << ring_member.size() << endl;
                     cout << "atom list in the ring: "; 
                     for (int j=0; j< ring_member.size(); j++){
                          cout << ring_member[j] << " ";
                          crd[0]=crd[0]+mol.x[ring_member[j]];
                          crd[1]=crd[1]+mol.y[ring_member[j]];
                          crd[2]=crd[2]+mol.z[ring_member[j]];
                     }
                     cout << endl;
                     crd[0]    = crd [0]/(float)ring_member.size(); cout << crd[0] << ", " ;
                     crd[1]    = crd [1]/(float)ring_member.size(); cout << crd[1] << ", " ;
                     crd[2]    = crd [2]/(float)ring_member.size(); cout << crd[2] << ", " << endl;
                     v_crd[0]  = 0;
                     v_crd[1]  = 0;
                     v_crd[2]  = 0;
                     for (int j=0; j < ring_member.size()-1; j++) {
                           va_crd[0] = mol.x[ring_member[j]]-crd[0];
                           vb_crd[0] = mol.x[ring_member[j+1]]-crd[0];
                           va_crd[1] = mol.y[ring_member[j]]-crd[1];
                           vb_crd[1] = mol.y[ring_member[j+1]]-crd[1];
                           va_crd[2] = mol.z[ring_member[j]]-crd[2];
                           vb_crd[2] = mol.z[ring_member[j+1]]-crd[2];
                           v_crd[0]  += (va_crd[1]*vb_crd[2]) - (va_crd[2]*vb_crd[1]);
                           v_crd[1]  += (va_crd[2]*vb_crd[0]) - (va_crd[0]*vb_crd[2]);
                           v_crd[2]  += (va_crd[0]*vb_crd[1]) - (va_crd[1]*vb_crd[0]);
                     }
                     va_crd[0] = mol.x[ring_member[ring_member.size()-1]]-crd[0];
                     vb_crd[0] = mol.x[ring_member[0]]-crd[0];
                     va_crd[1] = mol.y[ring_member[ring_member.size()-1]]-crd[1];
                     vb_crd[1] = mol.y[ring_member[0]]-crd[1];
                     va_crd[2] = mol.z[ring_member[ring_member.size()-1]]-crd[2];
                     vb_crd[2] = mol.z[ring_member[0]]-crd[2];
                     v_crd[0]  += (va_crd[1]*vb_crd[2]) - (va_crd[2]*vb_crd[1]);
                     v_crd[1]  += (va_crd[2]*vb_crd[0]) - (va_crd[0]*vb_crd[2]);
                     v_crd[2]  += (va_crd[0]*vb_crd[1]) - (va_crd[1]*vb_crd[0]);
                     v_crd[0]   = v_crd[0]/(float)ring_member.size();
                     v_crd[1]   = v_crd[1]/(float)ring_member.size();
                     v_crd[2]   = v_crd[2]/(float)ring_member.size();
                     r_vec  = pow(v_crd[0], 2) + pow(v_crd[1], 2) + pow(v_crd[2], 2);
                     r_ring = pow(atom_mol.x[i]-crd[0], 2) + pow(atom_mol.y[i]-crd[1], 2) + pow(atom_mol.z[i]-crd[2], 2);
                     cout << "first atom in the ring: " << atom_mol.id[i] << endl;
                     r_vec  = sqrt(r_vec);
                     r_ring = sqrt(r_ring);
                     ph4_mol.id[count]= count + 1;        cout << ph4_mol.id[count] << ", " << atom_mol.id[i] << ", ";
                     ph4_mol.x[count] = crd[0];           cout << ph4_mol.x[count]  << ", ";
                     ph4_mol.y[count] = crd[1];           cout << ph4_mol.y[count]  << ", ";
                     ph4_mol.z[count] = crd[2];           cout << ph4_mol.z[count]  << ", ";
                     ph4_mol.v_x[count] = v_crd[0]/r_vec; cout << ph4_mol.v_x[count]<< ", ";
                     ph4_mol.v_y[count] = v_crd[1]/r_vec; cout << ph4_mol.v_y[count]<< ", ";
                     ph4_mol.v_z[count] = v_crd[2]/r_vec; cout << ph4_mol.v_z[count]<< ", ";
                     ph4_mol.radius[count] = r_ring;      cout << ph4_mol.radius[count]  << ", ";
                     ph4_mol.ph4_types[count] = atom_mol.ph4_types[i];
                     cout << ph4_mol.ph4_types[count]  << endl;
                     count++;
                 }
             }
*/
             else if (atom_mol.ph4_types[i].compare("aroAcc") == 0){
/*                 find_ring = true;
                 for(int j=0; j < visited.size(); j++) {
                     if (atom_mol.id[i] == visited[j]) {
                         find_ring = false;
                     }
                 }
                 if (find_ring){
                     ring_member.clear();
                     find_ring = find_ring_member(mol, atom_mol.id[i], visited, ring_member);
                     crd[0]= 0.0;
                     crd[1]= 0.0;
                     crd[2]= 0.0;
                     cout << "size of the ring: " << ring_member.size() << endl;
                     cout << "atom list in the ring: ";
                     for (int j=0; j< ring_member.size(); j++){
                          cout << ring_member[j] << " ";
                          crd[0]=crd[0]+mol.x[ring_member[j]];
                          crd[1]=crd[1]+mol.y[ring_member[j]];
                          crd[2]=crd[2]+mol.z[ring_member[j]];
                     }
                     cout << endl;
                     crd[0] = crd [0]/(float)ring_member.size(); cout << crd[0] << ", " ;
                     crd[1] = crd [1]/(float)ring_member.size(); cout << crd[1] << ", " ;
                     crd[2] = crd [2]/(float)ring_member.size(); cout << crd[2] << ", " << endl;
                     v_crd[0]  = 0;
                     v_crd[1]  = 0;
                     v_crd[2]  = 0;
                     for (int j=0; j < ring_member.size()-1; j++) {
                           va_crd[0] = mol.x[ring_member[j]]-crd[0];
                           vb_crd[0] = mol.x[ring_member[j+1]]-crd[0];
                           va_crd[1] = mol.y[ring_member[j]]-crd[1];
                           vb_crd[1] = mol.y[ring_member[j+1]]-crd[1];
                           va_crd[2] = mol.z[ring_member[j]]-crd[2];
                           vb_crd[2] = mol.z[ring_member[j+1]]-crd[2];
                           v_crd[0]  += (va_crd[1]*vb_crd[2]) - (va_crd[2]*vb_crd[1]);
                           v_crd[1]  += (va_crd[2]*vb_crd[0]) - (va_crd[0]*vb_crd[2]);
                           v_crd[2]  += (va_crd[0]*vb_crd[1]) - (va_crd[1]*vb_crd[0]);
                     }
                     va_crd[0] = mol.x[ring_member[ring_member.size()-1]]-crd[0];
                     vb_crd[0] = mol.x[ring_member[0]]-crd[0];
                     va_crd[1] = mol.y[ring_member[ring_member.size()-1]]-crd[1];
                     vb_crd[1] = mol.y[ring_member[0]]-crd[1];
                     va_crd[2] = mol.z[ring_member[ring_member.size()-1]]-crd[2];
                     vb_crd[2] = mol.z[ring_member[0]]-crd[2];
                     v_crd[0]  += (va_crd[1]*vb_crd[2]) - (va_crd[2]*vb_crd[1]);
                     v_crd[1]  += (va_crd[2]*vb_crd[0]) - (va_crd[0]*vb_crd[2]);
                     v_crd[2]  += (va_crd[0]*vb_crd[1]) - (va_crd[1]*vb_crd[0]);
                     v_crd[0]  = v_crd[0]/(float)ring_member.size();
                     v_crd[1]  = v_crd[1]/(float)ring_member.size();
                     v_crd[2]  = v_crd[2]/(float)ring_member.size();
                     v_crd[0]  = (va_crd[1]*vb_crd[2]) - (va_crd[2]*vb_crd[1]);
                     v_crd[1]  = (va_crd[2]*vb_crd[0]) - (va_crd[0]*vb_crd[2]);
                     v_crd[2]  = (va_crd[0]*vb_crd[1]) - (va_crd[1]*vb_crd[0]);
                     r_vec  = pow(v_crd[0], 2) + pow(v_crd[1], 2) + pow(v_crd[2], 2);
                     r_ring = pow(atom_mol.x[i]-crd[0], 2) + pow(atom_mol.y[i]-crd[1], 2) + pow(atom_mol.z[i]-crd[2], 2);
                     cout << "first atom in the ring: " << atom_mol.id[i] << endl;
                     r_vec  = sqrt(r_vec);
                     r_ring = sqrt(r_ring);
                     ph4_mol.id[count]= count + 1;        cout << ph4_mol.id[count] << ", " << atom_mol.id[i] << ", ";
                     ph4_mol.x[count] = crd[0];           cout << ph4_mol.x[count]  << ", ";
                     ph4_mol.y[count] = crd[1];           cout << ph4_mol.y[count]  << ", ";
                     ph4_mol.z[count] = crd[2];           cout << ph4_mol.z[count]  << ", ";
                     ph4_mol.v_x[count] = v_crd[0]/r_vec; cout << ph4_mol.v_x[count]<< ", ";
                     ph4_mol.v_y[count] = v_crd[1]/r_vec; cout << ph4_mol.v_y[count]<< ", ";
                     ph4_mol.v_z[count] = v_crd[2]/r_vec; cout << ph4_mol.v_z[count]<< ", ";
                     ph4_mol.radius[count] = r_ring;      cout << ph4_mol.radius[count]  << ", ";
                     ph4_mol.ph4_types[count] = "aromatic";//hardcoded!
                     cout << ph4_mol.ph4_types[count]  << endl;
                     count++;
                 }
*/
                 nbrs.clear();
                 for (int j=0; j < mol.num_bonds; j++){
                      if (mol.bonds_origin_atom[j] == atom_mol.id[i]) {
                          nbrs.push_back(mol.bonds_target_atom[j]);
                      }
                      else if (mol.bonds_target_atom[j] == atom_mol.id[i]) {
                          nbrs.push_back(mol.bonds_origin_atom[j]);
                      }
                 }

                 ph4_mol.id[count]= count + 1;        //*cout << ph4_mol.id[count] << ", " << atom_mol.id[i]+1 << ", ";
                 ph4_mol.x[count] = atom_mol.x[i];    //*cout << ph4_mol.x[count]  << ", ";
                 ph4_mol.y[count] = atom_mol.y[i];    //*cout << ph4_mol.y[count]  << ", ";
                 ph4_mol.z[count] = atom_mol.z[i];    //*cout << ph4_mol.z[count]  << ", ";

                 v_crd[0] = 0;
                 v_crd[1] = 0;
                 v_crd[2] = 0;
                 for(int j=0; j <nbrs.size(); j++) {
                     //the HB vector is the sum of all vectors of the neighbor atoms to the acceptor
                     v_crd[0] += (atom_mol.x[i]-mol.x[nbrs[j]])/(float(nbrs.size()));
                     v_crd[1] += (atom_mol.y[i]-mol.y[nbrs[j]])/(float(nbrs.size()));
                     v_crd[2] += (atom_mol.z[i]-mol.z[nbrs[j]])/(float(nbrs.size()));
                 }
                 r_vec  = (v_crd[0]*v_crd[0]) + (v_crd[1]*v_crd[1]) + (v_crd[2]*v_crd[2]);
                 r_vec  = sqrt(r_vec);
                 ph4_mol.v_x[count] = v_crd[0]/r_vec; //*cout << ph4_mol.v_x[count]<< ", ";
                 ph4_mol.v_y[count] = v_crd[1]/r_vec; //*cout << ph4_mol.v_y[count]<< ", ";
                 ph4_mol.v_z[count] = v_crd[2]/r_vec; //*cout << ph4_mol.v_z[count]<< ", ";
                 ph4_mol.radius[count] = r;           //*cout << ph4_mol.radius[count]  << ", ";                     
                 ph4_mol.ph4_types[count] = "acceptor";
                 //*cout << ph4_mol.ph4_types[count]  << endl;
                 count++;
             }
 
             else if((atom_mol.ph4_types[i].compare("positive") == 0) || (atom_mol.ph4_types[i].compare("negative") == 0)){
                     //to be modified by lingling
                     ph4_mol.id[count]= count + 1;        //*cout << ph4_mol.id[count] << ", " << atom_mol.id[i]+1 << ", ";
                     ph4_mol.x[count] = atom_mol.x[i];    //*cout << ph4_mol.x[count]  << ", ";
                     ph4_mol.y[count] = atom_mol.y[i];    //*cout << ph4_mol.y[count]  << ", ";
                     ph4_mol.z[count] = atom_mol.z[i];    //*cout << ph4_mol.z[count]  << ", ";
                     ph4_mol.v_x[count] = 1;              //*cout << ph4_mol.v_x[count]<< ", ";
                     ph4_mol.v_y[count] = 0;              //*cout << ph4_mol.v_y[count]<< ", ";
                     ph4_mol.v_z[count] = 0;              //*cout << ph4_mol.v_z[count]<< ", ";
                     ph4_mol.radius[count] = r;           //*cout << ph4_mol.radius[count]  << ", ";
                     ph4_mol.ph4_types[count] = atom_mol.ph4_types[i];
                     //*cout << ph4_mol.ph4_types[count]  << endl;
                     count++;
             }

        }
        
        //*cout << "count=" << count << ", num_features=" << ph4_mol.num_features <<endl; 
        if (count != ph4_mol.num_features) {
            //*cout << "warning: number of features in ph4_mol don't match size of ph4_mol." << endl;
            //*cout << "Reset ph4_mol.number_features to " << count  << endl;
            ph4_mol.num_features = count;
        }

        for(int j=0; j<ph4_mol.num_features; j++)
            ph4_mol.contri_to_score[j] = false;
        //initialize the bool vector contri_to_score. in the later stage of the code, only the ph4 points in cad
        //are labeled as true so that we can write out matched_ph4_only in the .submit_ph4_mol2 function.

        delete [] crd;
        delete [] v_crd;
        delete [] va_crd;
        delete [] vb_crd;
        delete [] vc_crd;
        visited.clear();
        ring_member.clear();
        return true;
}

// +++++++++++++++++++++++++++++++++++++++++
int
Ph4_Score::calc_num_match(Ph4_Struct & Ph4_ref, Ph4_Struct & ph4_mol){//compute the number of matches of any ph4 model to itself
   int  num_ph4_match=0;
   //double proj_tol = 0.5; // v_proj cutoff, 60 deg:0.5, 45 deg:0.7071, 30 deg:0.8660.
   //double proj_tol = 0.7071;
   //double proj_tol = 0.8660;
   double tmp_residual2;
   int    temp_m_e;

   for (int i=0;i<Ph4_ref.num_features;i++){
        temp_m_e=-1;
        tmp_residual2=100.0;//hardcoded, have to be larger than any allowable dist2/v_proj.
        for (int j=0;j<ph4_mol.num_features;j++){
            if (Ph4_ref.ph4_types[i] == ph4_mol.ph4_types[j]){
                double dist2 =  ((Ph4_ref.x[i]-ph4_mol.x[j])*(Ph4_ref.x[i]-ph4_mol.x[j]))
                              + ((Ph4_ref.y[i]-ph4_mol.y[j])*(Ph4_ref.y[i]-ph4_mol.y[j]))
                              + ((Ph4_ref.z[i]-ph4_mol.z[j])*(Ph4_ref.z[i]-ph4_mol.z[j]));
                double v_proj=  Ph4_ref.v_x[i]*ph4_mol.v_x[j]
                              + Ph4_ref.v_y[i]*ph4_mol.v_y[j]
                              + Ph4_ref.v_z[i]*ph4_mol.v_z[j];
                if (Ph4_ref.ph4_types[i].compare("aromatic")==0||Ph4_ref.ph4_types[i].compare("ring")==0){
                       v_proj= fabs(v_proj);
                }
                if (dist2 <= (Ph4_ref.radius[i]*Ph4_ref.radius[i])){
                   if (v_proj >= ph4_proj_cos) {
                       if (dist2/v_proj <= tmp_residual2){
                           temp_m_e=j;
                           tmp_residual2 = dist2/v_proj;
                       }
                  }
                }
            }
        }
        if (temp_m_e>=0){
           num_ph4_match++;
        }
    }

    cout << endl;

return num_ph4_match;
}

// +++++++++++++++++++++++++++++++++++++++++
double
Ph4_Score::calc_ph4_similarity(Ph4_Struct & ph4_mol, char comparison_type){
    // to be modified by lingling
    // this function quantifies the similarity between a reference ph4 and a pose ph4
    // comparison_type depends on if the reference is the receptor ("c") or the ligand ("o")
    double         value = 0.0;
    int        num_match = 0;
    double     residual2 = 0.0;
    double tmp_residual2 = 0.0;
    int    temp_m_e;

    int num_types = list_comp.size(); 
    //double proj_tol = 0.5; // v_proj cutoff, 60 deg:0.5, 45 deg:0.7071, 30 deg:0.8660.
    //double proj_tol = 0.7071;
    //double proj_tol = 0.8660;    
           
    //double ris_comp[num_types];
    //double tolris_comp[num_types];
    //match_comp = new int [num_types];
    //components = new double [num_types];
    //match_num  = new int [3]; //match_num = [num_match_tot, max_match_ref, max_match_mol]:
    //match_term = new double [2]; // match_term = [match_rate, match_resid]; 
    // initialize the components for decomposition.
    for (int i=0;i<num_types;i++) {
         //ris_comp[i]=0.0;
         //tolris_comp[i]=0.0;
         match_comp[i]=0;
         components[i]=0.0;
    }//residual and matches per component, to be modified. not used for FMS score now

    //*cout << "\n\n";
    
    for (int i=0;i<Ph4_ref.num_features;i++){
        int k = component_number(Ph4_ref.ph4_types[i],list_comp);
        //*cout << "Ph4 point " <<i+1 << "in Ph4_ref" << "...";
        temp_m_e=-1;
        tmp_residual2=100.0;//hardcoded, have to be larger than any allowable dist2/v_proj.
        for (int j=0;j<ph4_mol.num_features;j++){
            if (Ph4_ref.ph4_types[i] == ph4_mol.ph4_types[j]){
                double dist2 =  ((Ph4_ref.x[i]-ph4_mol.x[j])*(Ph4_ref.x[i]-ph4_mol.x[j]))
                              + ((Ph4_ref.y[i]-ph4_mol.y[j])*(Ph4_ref.y[i]-ph4_mol.y[j]))
                              + ((Ph4_ref.z[i]-ph4_mol.z[j])*(Ph4_ref.z[i]-ph4_mol.z[j]));
                double v_proj=  Ph4_ref.v_x[i]*ph4_mol.v_x[j]
                              + Ph4_ref.v_y[i]*ph4_mol.v_y[j]
                              + Ph4_ref.v_z[i]*ph4_mol.v_z[j];
                if (Ph4_ref.ph4_types[i].compare("aromatic")==0||Ph4_ref.ph4_types[i].compare("ring")==0){
                       v_proj= fabs(v_proj);
                }
                if (dist2 <= (Ph4_ref.radius[i]*Ph4_ref.radius[i])){
                   if (v_proj >= ph4_proj_cos) {
                       if (dist2/v_proj <= tmp_residual2){
                           temp_m_e=j;
                           tmp_residual2 = dist2/v_proj;
                       }
                  }
                }
            }
        }
        if (temp_m_e>=0){
            num_match++;
            match_comp[k]++;
            residual2 += tmp_residual2;
            ph4_mol.contri_to_score[temp_m_e]=true;
            //temp_m_e stores the final index of the ph4 point in ph4_mol contributing to score.
            //*cout << "match to Ph4 point" << temp_m_e+1 << "in ph4_mol" << endl;
        }
        else {
            //*cout << "has no match" << endl;
        }
    }
    
    //num_match = calc_num_match(Ph4_ref, ph4_mol);

    //*cout << "v_proj cutoff is:" << proj_tol << " angstrom." << endl;
    //*cout << "\n\n" << endl;
    
    //double wm = 5.0; //weight of match_rate, wm=1, 5, 10; 
 
    //int num_match_tot = num_match;
    match_num[0] = num_match;
    //*cout << "total number of match:" << " " << num_match << " " << endl;

     if (comparison_type == 'o'){
        int max_match_ref = Ph4_ref.num_features;
        int max_match_mol = ph4_mol.num_features;
        match_num[1] = max_match_ref;
        match_num[2] = max_match_mol;
        if (num_match > 0 ){
           ///*
           value = (float(num_match)/max_match_ref);//hardcoded for now, to be modified.
           //*cout << "max_match_ref num = " << max_match_ref << endl;
           //*cout << "Ph4_ref.num_features=" << Ph4_ref.num_features <<endl;
           //use max_match_ref all the time for virtual screening. and the max() for horizontal pruning; to be modified
           //*/
           /*
           if (max_match_ref > max_match_mol) {
              value = (float(num_match)/max_match_ref);//hardcoded for now, to be modified.
              cout << "max_match_ref num = " << max_match_ref << endl;
            }
           else {
              value = (float(num_match)/max_match_mol);
              cout << "max_match_mol num = " << max_match_mol << endl;
           }
           */
         
           value = 1-value; //equation for computing FMS. to be modified.
           match_term[0] = value;                            //*cout << "match rate term: " << value << endl;
           match_term[1] = sqrt(residual2/float(num_match)); //*cout << "residual   term: " << sqrt(residual2/num_match) << endl;
           //*cout << "ph4_rate_k*value " << ph4_rate_k*match_term[0] << "and value " << value  << endl;
           //*cout << "resid2/num_match " << residual2/float(num_match)  << " and sqrt" << match_term[1]<< endl;
           value = ph4_rate_k*match_term[0] + match_term[1];
        }
        else {
             value = ph4_max_x; // maximum FMS score for ph4 model pairs with no matches
             match_term[0] = value/ph4_rate_k;  //*cout << "match rate term: " << value << endl;
             match_term[1] = 0;         //*cout << "residual   term: 0"  << endl;
        }
     }

     if  (comparison_type == 'c'){
        //change the way FMS score is computed when compare to a receptor reference. to be modified.
        cout << "receptor based ph4 similarity functionality under construction" <<endl;
/*
        int max_match_ref = Ph4_ref.num_features;
        int max_match_mol = ph4_mol.num_features;

        if (num_match > 0 ){
        }
        else {
        }
*/
     }

/*    for(int i = 0;i < num_types; i++) {
         if (match_comp[i] > 0) {
            components[i]= -sqrt(tolris_comp[i]/match_comp[i]) + sqrt(ris_comp[i]/match_comp[i]) + (para/match_comp[i]);//FMS decomposition on the types of feature.
         }
         else {
             components[i] = 1;
         }
      }
*/

    return value;
}

// +++++++++++++++++++++++++++++++++++++++++
// A light weight version of compute score which uses the ligand based pharmacophore
// Reads in a mol to score, a mol to use as the reference, and a compute type. Will modify the input mol's score and featues. Initializes the ref_mol as a reference and computes the score, should not overwrite the score of mol
// ph4_compare_type can either be 'overlap' or 'compatible'
bool Ph4_Score::compute_score(DOCKMol & mol, DOCKMol & ref_mol, std::string ph4_compare_type, AMBER_TYPER & c_typer){
      Ph4_Struct Ph4_pose;
      compute_ph4(mol,Ph4_pose); // compare Ph4_pose with the Ph4_ref
      this->ph4_compare_type = ph4_compare_type;
      char *compare_type;
      compare_type = new char [ph4_compare_type.size()+1];
      strcpy (compare_type, ph4_compare_type.c_str());

      initialize(c_typer, ref_mol); //initilize the ref mol as Ph4 reference
      //float Ph4_sim      = calc_ph4_similarity(Ph4_pose);
      mol.current_score  = calc_ph4_similarity(Ph4_pose,compare_type[0]);
      mol.current_data   = output_score_summary(mol);
      delete [] compare_type;
      return true;

}

// +++++++++++++++++++++++++++++++++++++++++
bool
Ph4_Score::compute_score(DOCKMol & mol)
{

      Ph4_Struct Ph4_pose;
      compute_ph4(mol,Ph4_pose); // compare Ph4_pose with the Ph4_ref

      char *compare_type;
      compare_type = new char [ph4_compare_type.size()+1];
      strcpy (compare_type, ph4_compare_type.c_str());

      //float Ph4_sim      = calc_ph4_similarity(Ph4_pose);
      mol.current_score  = calc_ph4_similarity(Ph4_pose,compare_type[0]);
      if(bool_out_cad_ph4){
         submit_ph4_mol2(Ph4_pose, cad_ph4_out_filename, mol.title, false);
      }
      if(bool_out_mat_ph4){
         submit_ph4_mol2(Ph4_pose, mat_ph4_out_filename, mol.title, true);
      }

 
      mol.current_data   = output_score_summary(mol);
      delete [] compare_type;
      return true;
}

// +++++++++++++++++++++++++++++++++++++++++
string
Ph4_Score::output_score_summary(DOCKMol & mol)
{

    ostringstream text;

    if (use_score) {
        text << DELIMITER << setw(STRING_WIDTH) << "Pharmacophore_Score:"
             << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl;
        text << DELIMITER << setw(STRING_WIDTH) << "FMS_num_match_tot:"
             << setw(FLOAT_WIDTH) << fixed << match_num[0] << endl;
        text << DELIMITER << setw(STRING_WIDTH) << "FMS_max_match_ref:"
             << setw(FLOAT_WIDTH) << fixed << match_num[1] << endl;
        text << DELIMITER << setw(STRING_WIDTH) << "FMS_max_match_mol:"
             << setw(FLOAT_WIDTH) << fixed << match_num[2] << endl;
        text << DELIMITER << setw(STRING_WIDTH) << "FMS_match_rate:"
             << setw(FLOAT_WIDTH) << fixed << match_term[0] << endl;
        text << DELIMITER << setw(STRING_WIDTH) << "FMS_match_resid:"
             << setw(FLOAT_WIDTH) << fixed << match_term[1] << endl;
        text << DELIMITER << setw(STRING_WIDTH) << "FMS_num_hydrophobic_matched:"
             << setw(FLOAT_WIDTH) << fixed << match_comp[0] << endl;
        text << DELIMITER << setw(STRING_WIDTH) << "FMS_num_donor_matched:"
             << setw(FLOAT_WIDTH) << fixed << match_comp[1] << endl;
        text << DELIMITER << setw(STRING_WIDTH) << "FMS_num_acceptor_matched:"
             << setw(FLOAT_WIDTH) << fixed << match_comp[2] << endl;
        text << DELIMITER << setw(STRING_WIDTH) << "FMS_num_aromatic_matched:"
             << setw(FLOAT_WIDTH) << fixed << match_comp[3] << endl;
        //text << DELIMITER << setw(STRING_WIDTH) << "FMS_num_aroAcc_matched:"
        //     << setw(FLOAT_WIDTH) << fixed << match_comp[4] << endl;
        text << DELIMITER << setw(STRING_WIDTH) << "FMS_num_positive_matched:"
             << setw(FLOAT_WIDTH) << fixed << match_comp[5] << endl;
        text << DELIMITER << setw(STRING_WIDTH) << "FMS_num_negative_matched:"
             << setw(FLOAT_WIDTH) << fixed << match_comp[6] << endl;
        text << DELIMITER << setw(STRING_WIDTH) << "FMS_num_ring_matched:"
             << setw(FLOAT_WIDTH) << fixed << match_comp[7] << endl;

        if (use_internal_energy){
            text << DELIMITER << setw(STRING_WIDTH) << "Internal_energy_repulsive:"
                 << setw(FLOAT_WIDTH) << fixed << mol.internal_energy << endl;
            }
        //Numer of matched features changed from a loop to hardcoded for formatting purpose Yuchen 10/24/2016
        //Each matched term has a component score that is not printed here in this version Yuchen 10/24/2016
    }
    return text.str();
}


// +++++++++++++++++++++++++++++++++++++++++
