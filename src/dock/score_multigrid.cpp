#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <math.h>
#include <sstream>
#include "amber_typer.h"
#include "base_score.h"
#include "dockmol.h"
#include "grid.h" 
#include "score.h"
#include "score_multigrid.h" 
#include "utils.h"
#include "trace.h"

using namespace std;



// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++
Multigrid_Energy_Score::Multigrid_Energy_Score()
{
    grid_file_names = NULL;
    energy_grids = NULL;
    vdw_ref_array = NULL;
    es_ref_array = NULL;
    mgweights_array = NULL;
    bltzmn_weight_array = NULL;
}

// +++++++++++++++++++++++++++++++++++++++++
Multigrid_Energy_Score::~Multigrid_Energy_Score()
{
   close();
}

// +++++++++++++++++++++++++++++++++++++++++
// This function could be called in the function in Master_score called close_all(),
// but is curently not.  the function is called in the destructure.
// This is important to clear the memory of scoring function at the end of
// the docking calculations.
void 
Multigrid_Energy_Score::close()
{
/*  
    WE NEED TO DELETE VECTORS 
*/
//    cout << "Multigrid_Energy_Score::close()[1]" << endl;
    footprint_ref_file.clear();
    footprint_txt_file.clear();
    gridweight_txt_file.clear();

//    for (int i = 0; i < numgrids;i++) {
//        energy_grids[i].clear_grid();
//        grid_file_names[i].clear(); // string
//    }
    delete[] grid_file_names;
    delete[] energy_grids;
//    cout << "Multigrid_Energy_Score::close()[2]" << endl;
//    delete energy_grids;
//    delete[] grid_file_names;
// the pose vectors are deleted after assigned.
//    delete[] vdw_pose_array;
//    delete[] es_pose_array;
    delete[] vdw_ref_array;
    delete[] es_ref_array;
    delete[] mgweights_array;
    delete[] bltzmn_weight_array; //BCF
}

// +++++++++++++++++++++++++++++++++++++++++
void
Multigrid_Energy_Score::input_parameters(Parameter_Reader & parm, bool & primary_score,
                               bool & secondary_score)
{
    use_primary_score = false;
    use_secondary_score = false;

    cout << "\nMultiGrid Score Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    string          tmp;

    if (!primary_score) {
        tmp = parm.query_param("multigrid_score_primary", "no", "yes no");
        use_primary_score = tmp == "yes";
        primary_score = use_primary_score;
    }

    if (!secondary_score) {
        //tmp = parm.query_param("multigrid_score_secondary", "no", "yes no");
        tmp = "no";
	    use_secondary_score = tmp == "yes";
        secondary_score = use_secondary_score;
    }

    if (use_primary_score || use_secondary_score)
        use_score = true;
    else
        use_score = false;

   input_parameters_main(parm, "multigrid_score");

}//end Multigrid_Energy_Score::input_parameters()



//Begin Multigrid_Energy_Score::input_parameters_main()
//This function reads the input parameters needed for necessary.
//In addition this same function is called in "score_descriptor.cpp"
//whenever multigrid is used when descriptor score is called.
//Employing the same function ensures that the question tree no matter
//whether it is called through descriptor score or multigrid remains consistent

void
Multigrid_Energy_Score::input_parameters_main(Parameter_Reader & parm, string parm_head)
{


    if (use_score) {
        rep_radius_scale = atof(parm.query_param(parm_head+"_rep_rad_scale", "1.0").c_str());
        if (rep_radius_scale <= 0.0) {
                cout << "ERROR:  Parameter must be a float greater than zero."
                        " Program will terminate."
                     << endl;
                exit(0);
        }
        vdw_scale = atof(parm.query_param(parm_head+"_vdw_scale", "1.0").c_str());
        es_scale = atof(parm.query_param(parm_head+"_es_scale", "1.0").c_str());

        numgrids       = atoi(parm.query_param( parm_head+"_number_of_grids", "20").c_str());
        if (numgrids > 50) 
            cout << "Warning: due to limited memory one, should not use more than 50." << endl; 
        grid_file_names = new string[numgrids];
        for (int i = 0; i < numgrids;i++) {
             stringstream ss1,ss2;
             ss1 << parm_head << "_grid_prefix" << i; 
             ss2 << "multigrid" << i; 
             grid_file_names[i] = parm.query_param( ss1.str(), ss2.str() );
             //grid_file_names.push_back(temp);
        }

        // note: modify code to make sure we can select only one of these options
        ir_ensemble = false;
        wt_txt = false;
        bltzmn = false; //BCF
        fp_txt = false;
        fp_mol = false;

	// Added MIR score - LEP 
        ir_ensemble = parm.query_param( "multigrid_score_individual_rec_ensemble", "no",
                                         "yes no") == "yes";
        if (!ir_ensemble){
		// Remove weights.txt for 6.6 release - sudipto
                wt_txt       = parm.query_param( "multigrid_score_weights_text", "no",
                                         "yes no") == "yes";
                if (!wt_txt){
                       // bltzmn =  parm.query_param( "dynamic_boltzmann_weighting", "no", //BCF
                       //                          "yes no") == "yes";
         
                       // if (!bltzmn){
                            fp_mol       = parm.query_param( parm_head+"_fp_ref_mol", "no",
                                                     "yes no") == "yes";
         
                            if (!fp_mol){
                                fp_txt       = parm.query_param( parm_head+"_fp_ref_text", "no",
                                                             "yes no") == "yes";
                            }
                }
        }

        //if (!ir_ensemble && !wt_txt && !fp_mol && !fp_txt && !bltzmn)  {  //BCF
        //    cout << "Choice one of the following: fp_mol,fp_txt,wt_txt, ir_ensemble, or bltzmn(in development)" << endl;

        //if (!ir_ensemble && !wt_txt && !fp_mol && !fp_txt)  {
        //   cout << "Choose one of the following: fp_mol,fp_txt,wt_txt or ir_ensemble" << endl;
        //   exit(0);
        //}


        if ( fp_mol )
        {
           footprint_ref_file = parm.query_param( parm_head+"_footprint_ref", "reference.mol2" );
        }
        else if (fp_txt )
        {
           footprint_txt_file = parm.query_param(parm_head+"_footprint_text", "reference.txt" );
        }
        else if ( wt_txt)
        {
           gridweight_txt_file = parm.query_param( parm_head+"_gridweight_text", "weights.txt" );
        }
        else if ( bltzmn ) {   //BCF
            float k = 0.0019872041; // Boltzmann Constant in kcal/mol/K
            //float bltzmn_temp;
            bltzmn_temp = atof(parm.query_param("scale_factor", "1.0").c_str());
            bltzmn_kt = k * bltzmn_temp;
         }

        else if ( ir_ensemble )
        {
           cout << "using ir_ensemble" << endl;
        }

        vdw_cor_scale = 0;
        vdw_euc_scale = 0;
        vdw_norm_scale = 0;
        es_cor_scale = 0;
        es_euc_scale = 0;
        es_norm_scale = 0;

        // select which footprint measure to use
        if (fp_mol || fp_txt) {
           
           //change input parameter logic to be exactly the same as footprint score 
           //for clarity and consistency.  Yuchen 11/3/2016
           string                     tmp;
           string      mg_fp_compare_type;

           mg_fp_compare_type = 
               parm.query_param(parm_head+"_foot_compare_type", "Euclidean", "Pearson Euclidean");

           if (mg_fp_compare_type == "Pearson") {
               cout << "-------------------------------------------------------------\n"
                       "You chose to use the correlation coefficient as the metric \n"
                       "to compare the footprints.  When the value is 1 then \n"
                       "there is perfect agreement between the two footprints.\n"
                       "When the value is 0 then there is poor agreement\n"
                       "between the two footprints.\n"
                       "-------------------------------------------------------------"
                    << endl;
               use_cor  = true;
               use_euc  = false;
               use_norm = false;
           }
           else if (mg_fp_compare_type == "Euclidean") {
               cout << "-------------------------------------------------------------\n"
                       "You chose to use the Euclidean distance as the metric \n"
                       "to compare the footprints.  When the value is 0 then \n"
                       "there is perfect agreement between the two footprints.\n"
                       "As the agreement gets worse between the two\n"
                       "footprints the value increases.\n"
                       "-------------------------------------------------------------"
                    << endl;
               use_euc  = true;
               use_cor  = false;
               use_norm = false;
               tmp = parm.query_param(parm_head+"_normalize_foot", "no", "yes no");
               if (tmp == "yes"){ 
                   use_norm = true;
                   use_euc = false;
                   use_cor = false;
               }
           }

           // note that only one of these should be chosen
           if (use_euc) {
               cout << "-------------------------------------------------------------\n"
                       "You chose to not normalize footprints which penalizes\n"
                       "the score of the pose if the footprint is bad.\n"
                       "-------------------------------------------------------------"
                    << endl;
               vdw_euc_scale =  atof(parm.query_param(parm_head+"_vdw_euc_scale", "1.0").c_str());
               es_euc_scale =  atof(parm.query_param(parm_head+"_es_euc_scale", "1.0").c_str());
           }
           else if (use_norm) {
               cout << "-------------------------------------------------------------\n"
                       "You chose to normalize footprints which penalizes\n"
                       "the score of the pose if the footprint is bad.\n"
                       "-------------------------------------------------------------"
                    << endl;
               vdw_norm_scale =  atof(parm.query_param(parm_head+"_vdw_norm_scale", "10.0").c_str());
               es_norm_scale =  atof(parm.query_param(parm_head+"_es_norm_scale", "10.0").c_str());
           }
           else if (use_cor) {
               cout << "-------------------------------------------------------------\n"
                       "You chose the Pearson metric which makes the score more\n"
                       "favorable if the footprint is good and penalizes\n"
                       "the score if the footprint is bad.\n"
                       "-------------------------------------------------------------"
                    << endl;
               vdw_cor_scale =  atof(parm.query_param(parm_head+"_vdw_cor_scale", "-10.0").c_str());
               es_cor_scale =  atof(parm.query_param(parm_head+"_es_cor_scale", "-10.0").c_str());
           }
        }
    }
} //end Multigrid_Energy_Score::input_parameters_main()


// +++++++++++++++++++++++++++++++++++++++++
// modified from score_descriptor
// 
bool 
Multigrid_Energy_Score::read_mgfootprint_txt(istream & ifs){

   string line;
   // get every string or number from file.

   vdw_ref_array = new float[numgrids];
   es_ref_array  = new float[numgrids];
   int count = 0;
 
   while (getline(ifs,line,'\n')) {
       vector < string >  tokens;
       
       if (!line.empty()){
           Tokenizer(line, tokens,' ');
           if (tokens.size() > 0){
               // ignore the column names and comments lines
               if (!(tokens[0].compare("resname") == 0) && !(tokens[0].compare(0,1,"#") == 0)){
               //if ((tokens[0].compare("gridname") == 0) ){
                   // terminate if wrong number of coloms in a line.
                   if (tokens.size() != 8){
                       cout << "A line in the  multigrid fp reference txt file has an incorrect" << endl
                            << "number of entries. Lines must have 6 entries separated by"
                            << " spaces as follows:" << endl
                            << "\t\tresname  resid  vdw_ref  es_ref hb_ref  vdw_pose  es_pose hb_pose" << endl
                            << "Note that lines beginning with a comment \"#\" will be ignored" << endl
                            << "Program will terminate." << endl;
                       cout << "line which fails on:" <<endl
                            << line << endl;
                       exit(1);

                   }
                   stringstream ss; ss.clear();
                   //temp_foot.resname = tokens[0];
                   if ( count > numgrids ) {
                       cout << "Error: num of line in multigrid fp input txt file"
                            << endl;
                       exit(1);
                   }
                   string gridname;
                   ss <<tokens[1]; ss >> gridname; ss.clear();
                   ss <<tokens[2]; ss >> vdw_ref_array[count]; ss.clear();
                   ss <<tokens[3]; ss >> es_ref_array[count]; ss.clear();
                   count = count + 1;
               }
           }
       }
   }
// else cout << "line is empty" << line << endl;
   if (count != numgrids){
       cout << "Error: multigrid fp input txt." << endl
            << "Num of lines:" << count << endl
            << "Num of grids:" << numgrids << endl;
       return false;
   }
   return true;
}

// 
bool 
Multigrid_Energy_Score::read_mgweights_txt(istream & ifs){

   string line;
   // get every string or number from file.

   mgweights_array = new float[numgrids];
   int count = 0;
 
   while (getline(ifs,line,'\n')) {
       vector < string >  tokens;
       
       if (!line.empty()){
           Tokenizer(line, tokens,' ');
           if (tokens.size() > 0){
               //cout << "I AM HERE" << endl;
               // ignore the column names and comments lines
               //if (!(tokens[0].compare("gridname") == 0) && !(tokens[0].compare(0,1,"#") == 0)){
               //if ((tokens[0].compare("gridname") == 0) ){
                   // terminate if wrong number of coloms in a line.
                   if (tokens.size() != 2){
                       cout << "A line in the multigrid weights txt file has an incorrect"
                            << endl 
                            << "number of entries. Lines must have 2 entries separated by"
                            << " spaces as follows:"
                            << endl
                            << "gridname grid_weight" 
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
                   //temp_foot.resname = tokens[0];
                   if ( count > numgrids ) {
                       cout << "Error: num of line in multigrid fp input txt file"
                            << endl;
                       exit(1);
                   }
                   string gridname;
                   ss <<tokens[0]; ss >> gridname; ss.clear();
                   //cout << gridname;
                   ss <<tokens[1]; ss >> mgweights_array[count]; ss.clear();
                   //cout << " " << mgweights_array[count] << endl;
                   count = count + 1;
               //}
           }
       }
   }
//   else cout << "line is empty" << line << endl;
   if (count != numgrids){
       cout << "Error: multigrid fp input txt." << endl
            << "Num of lines:" << count << endl
            << "Num of grids:" << numgrids << endl;
       return false;
   }
   return true;
}

// +++++++++++++++++++++++++++++++++++++++++
// modified from score_descriptor
// get footprint reference
void
Multigrid_Energy_Score::submit_footprint_reference(AMBER_TYPER & typer)
{ Trace trace("Multigrid_Energy_Score::submit_footprint_reference()");
    ifstream        footprint_ref;
    ifstream        footprint_txt;

    // if footprint is cal from a ref in mol2 file
    if (fp_mol){
       footprint_ref.open(footprint_ref_file.c_str());
       bool temp_mol2_OK = Read_Mol2(footprint_reference, footprint_ref, false, false, false);
       if(!(temp_mol2_OK)){
           cout << "there is a problem with footprint_reference mol2 file" <<endl;
           exit(0);
       }
       footprint_ref.close();
       bool read_vdw = true;
       bool use_chem = false;
       bool use_ph4  = false;
       bool use_volume = false;
       typer.prepare_molecule(footprint_reference, read_vdw, use_chem, use_ph4, use_volume);

       vdw_ref_array = new float[numgrids];
       es_ref_array  = new float[numgrids];
       bool reference_prob = compute_multigrid(footprint_reference,vdw_ref_array, es_ref_array);
      
       //storing footprint information in the footprints field of DOCKMol object, 
       //to enable printing of footprint as a text file.  Yuchen 11/1/2016
       for (int i=0; i<numgrids; i++) {
            FOOTPRINT_ELEMENT temp_footprint_ele;
            temp_footprint_ele.resname = "grid";
            temp_footprint_ele.resid   = i;
            temp_footprint_ele.vdw     = vdw_ref_array[i];
            temp_footprint_ele.es      = es_ref_array[i];
            temp_footprint_ele.hb      = 0;

            footprint_reference.footprints.push_back(temp_footprint_ele);

        }

       if(!reference_prob) {
           cout << "ref not in grid" << endl;
           exit(0);
       }
    }
    // if footprint is in a text file
    trace.boolean("footpint_txt_bool=", fp_txt);
    if(fp_txt){
       footprint_txt.open(footprint_txt_file.c_str());
       bool temp_txt_OK = read_mgfootprint_txt(footprint_txt);
       if(!(temp_txt_OK)){
           cout << "The footprint text reference file, " << footprint_txt_file << ", does not exist." <<endl;
           cout << "there is a problem with fp txt file" <<endl;
           exit(0);
       }

       for (int i=0; i<numgrids; i++) {
            FOOTPRINT_ELEMENT temp_footprint_ele;
            temp_footprint_ele.resname = "grid";
            temp_footprint_ele.resid   = i;
            temp_footprint_ele.vdw     = vdw_ref_array[i];
            temp_footprint_ele.es      = es_ref_array[i];
            temp_footprint_ele.hb      = 0;

            footprint_reference.footprints.push_back(temp_footprint_ele);
        }
    }
//    // if tresholds are spesified then calculate ranges for ref.
//    if(desc_foot_specify_a_threshold){
//         ref_ranges = range_satisfying_threshold_all(footprint_reference.footprints);
//    }

    return;
}

// +++++++++++++++++++++++++++++++++++++++++
// get grid weights
void
Multigrid_Energy_Score::submit_wts()
{
    ifstream        gridweights_txt;
    gridweights_txt.open(gridweight_txt_file.c_str());
    bool temp_txt_OK = read_mgweights_txt(gridweights_txt);
    if(!(temp_txt_OK)){
        cout << "The mgweight file, " << gridweight_txt_file << ", does not exit" << endl;
        cout << "there is a problem with fp txt file" <<endl;
        exit(0);
    }

}

// +++++++++++++++++++++++++++++++++++++++++
void
Multigrid_Energy_Score::initialize(AMBER_TYPER & typer)
{
    cout << "To cite Multigrid Energy Score use: \n Balius, T. E.; Mukherjee, S.; Rizzo, R. C. Implementation and Evaluation of a Docking-rescoring Method using Molecular Footprint Comparisons. J. Comput. Chem., 2011, 32, 2273-2289\n" << endl;

    if (use_score) {
        cout << "Initializing Multi Grid Score Routines..." << endl;
        // array of pointers.
        energy_grids = new Energy_Grid[numgrids];
        cout << "Reading in " << numgrids << " grids" << endl;
        for (int i = 0; i < numgrids;i++) {
          //// 2011-01-07 -- trent e balius 
          //// reset got_the_grid so that more than one grid can be read in.
          //energy_grids[i] = new Energy_Grid();
          energy_grids[i].got_the_grid = false;
          energy_grids[i].get_instance(grid_file_names[i]); 
          init_vdw_energy(typer, energy_grids[i].att_exp, energy_grids[i].rep_exp);
        }
        //compute_multigrid(footprint_reference,vdw_ref_array, es_ref_array);
        if (fp_mol || fp_txt )
           submit_footprint_reference(typer);
        if (wt_txt)
           submit_wts();
        if (bltzmn)     //BCF
           bltzmn_weight_array = new float[numgrids];
    } 
}

// +++++++++++++++++++++++++++++++++++++++++
bool
Multigrid_Energy_Score::compute_multigrid(DOCKMol & mol,float * vdw_array, float * es_array)
{
    int atom;

    if (use_score) {

        float total = 0.0;


        for (int i = 0; i < numgrids;i++) { 
            // check to see if molecule is inside grid box
            for (atom = 0; atom < mol.num_atoms; atom++) {
                if (!energy_grids[i].is_inside_grid_box(mol.x[atom], mol.y[atom], mol.z[atom])) {
                    mol.current_score = -MIN_FLOAT;  // arbitrarily large score
                    mol.current_data = "ERROR:  Conformation could not be scored."
                        "\nConformation not completely within grid box.\n";
                    return false;
                }
            }
            float es_val = 0.0;
            float vdw_val = 0.0;
            for (atom = 0; atom < mol.num_atoms; atom++) {
   
                if (mol.atom_active_flags[atom]) {
   
                    energy_grids[i].find_grid_neighbors(mol.x[atom], mol.y[atom], mol.z[atom]);
   
                    // if we account for weighting here we can't report the true grid values.
                    // move scaling to compute_score
                    //vdw_val +=
                    //    ((vdwA[mol.amber_at_id[atom]] * energy_grids[i].interpolate(energy_grids[i].avdw)) -
                    //     (vdwB[mol.amber_at_id[atom]] * energy_grids[i].interpolate(energy_grids[i].bvdw))) *
                    //    vdw_scale;
                    //es_val += mol.charges[atom] * energy_grids[i].interpolate(energy_grids[i].es) * es_scale;

                    vdw_val +=
                        ((vdwA[mol.amber_at_id[atom]] * energy_grids[i].interpolate(energy_grids[i].avdw)) -
                         (vdwB[mol.amber_at_id[atom]] * energy_grids[i].interpolate(energy_grids[i].bvdw))) ;
                    es_val += mol.charges[atom] * energy_grids[i].interpolate(energy_grids[i].es);
                }
            }
   
            // if we account for weighting here we can't report the true grid values.
            // move scaling to compute_score
            //if (wt_txt){
            //    vdw_val = mgweights_array[i]*vdw_val ;
            //    es_val  = mgweights_array[i]*es_val;
            //}

            total = total + vdw_val + es_val;
            vdw_array[i] = vdw_val;
            es_array[i]  = es_val;
        }
    return true;
    }
}

// modified from code in descriptor score;
float correlation(float * x, float * y,int size){
   float sumx=0,sumy=0,sumxy=0,sumsqx=0,sumsqy=0;
   float r=0,r_1=0;
   for (int i=0;i<size;i++)
   {
       sumx = sumx + x[i];
       sumy = sumy + y[i];
       sumxy = sumxy + x[i]*y[i];
       sumsqx = sumsqx + x[i]*x[i];
       sumsqy = sumsqy + y[i]*y[i];
   }
   r_1 = sqrt((size*sumsqx-sumx*sumx)*(size*sumsqy-sumy*sumy));
   if(r_1!=0.0)
   {
       //r = (n*sumxy - sumx*sumy)/ r_1;
       r = (size*sumxy - sumx*sumy)/ r_1;
   // the r value formular is here;
   // note: the r value here can be negative
   }else {
       cout << "r is not meaningful" << endl;
       r = -0.0;
   }
   return r;
}
float norm(float * x,int size){
   float sumsq = 0;
   for (int i=0;i<size;i++){
       sumsq = sumsq + x[i]*x[i];
   }
   return sqrt(sumsq);
}
float euclidean(float * x, float * y,int size){
   //float *xny; // x-y
   //xny = new float [size];
   float xny[size]; // x-y
   for (int i=0;i<size;i++)
       xny[i] = x[i] - y[i];
   //delete[] xny;
   return norm(xny,size); 
}

void normalize(float * array, int size, float * array_new){
// caculate 2nd norm
   float n = norm(array,size);
// devied all entrees by 2nd norm
   for (int i=0;i<size;i++)
        array_new[i] = array[i]/n;
   return;
}

float  sum(float *x,int size){
   float sumx = 0;
   for (int i=0;i<size;i++){
       sumx = sumx + x[i];
   }
   return sumx;
}

float  wt_sum(float *x, float *wt,int size){

   float sumx = 0;

   for (int i=0;i<size;i++){
       sumx = sumx + wt[i] * x[i];
   }
   return sumx;
}


void printarrays(float *x,float *y, int size, ostringstream &sout){
sout << "#  array length = " << size << endl; 
sout << "#  printing grid arrays:" << endl;
sout << "## name: VDW, ES"<< endl;
for (int i=0;i<size;i++) 
   sout << "# grid"<< i << ": "<< x[i] << ", " << y[i] << endl; 
return;
}

void printarrays(float *x,float *y,float *w, int size, ostringstream &sout){
sout << "#  array length = " << size << endl;
sout << "#  printing grid arrays:" << endl;
sout << "## name: VDW, ES, weight"<< endl;
for (int i=0;i<size;i++)
   sout <<  "# grid"<< i << ": "<< x[i] << ", " << y[i] << ", " << w[i] << endl;
return;
}


// +++++++++++++++++++++++++++++++++++++++++
// compute_ir_ensemble_score used in individual ensemble rececptor docking
// Returns the score of mol wrt grid# grid_num
// grid_num is now stored inside the dockmol object
bool
Multigrid_Energy_Score::compute_ir_ensemble_score(DOCKMol & mol)
{
    int atom;

    if (!ir_ensemble) {
        cout << "Cannot use compute_score(DOCKMol&, int) unless using ir_emsemble";
        return false;
    }

    if (use_score) {

        // check to see if molecule is inside grid box
        for (atom = 0; atom < mol.num_atoms; atom++) {
            if (!energy_grids[mol.grid_num].is_inside_grid_box(mol.x[atom], mol.y[atom], mol.z[atom])) {
                mol.current_score = -MIN_FLOAT;  // arbitrarily large score
                mol.current_data = "ERROR:  Conformation could not be scored."
                    "\nConformation not completely within grid box.\n";
                return false;
            }
        }

        // compute es and vdw energies for mol
        es_sum  = 0.0;
        vdw_sum = 0.0;
        for (atom = 0; atom < mol.num_atoms; atom++) {

            if (mol.atom_active_flags[atom]) {

                energy_grids[mol.grid_num].find_grid_neighbors(mol.x[atom], mol.y[atom], mol.z[atom]);

                vdw_sum +=
                    ((vdwA[mol.amber_at_id[atom]] * energy_grids[mol.grid_num].interpolate(energy_grids[mol.grid_num].avdw)) -
                     (vdwB[mol.amber_at_id[atom]] * energy_grids[mol.grid_num].interpolate(energy_grids[mol.grid_num].bvdw))) ;
                es_sum += mol.charges[atom] * energy_grids[mol.grid_num].interpolate(energy_grids[mol.grid_num].es);
            }
        }

        return true;
    }
    return false; // score_ligands = no
}


// +++++++++++++++++++++++++++++++++++++++++
bool
Multigrid_Energy_Score::compute_score(DOCKMol & mol)
{
        // see whether pose is on the grid.
        vdw_pose_array = new float[numgrids];
        es_pose_array  = new float[numgrids];
        bool ongrid = compute_multigrid(mol,vdw_pose_array, es_pose_array);

        //storing footprint information in the footprints field of DOCKMol object,
        //to enable printing of footprint as a text file.  Yuchen 11/1/2016
        for (int i=0; i<numgrids; i++) {
            FOOTPRINT_ELEMENT temp_footprint_ele;
            temp_footprint_ele.resname = "grid";
            temp_footprint_ele.resid   = i;
            temp_footprint_ele.vdw     = vdw_pose_array[i];
            temp_footprint_ele.es      = es_pose_array[i];
            temp_footprint_ele.hb      = 0;

            mol.footprints.push_back(temp_footprint_ele);

        }
 
        if (!ongrid) {
           //cout << "pose out of grid"<< endl;
           return false;
        }

        //Compute score for ir_emsemble
        if (ir_ensemble) compute_ir_ensemble_score(mol);
        else {
                // weight the grids         
                if (wt_txt) {
                  vdw_sum = wt_sum(vdw_pose_array,mgweights_array,numgrids);
                  es_sum = wt_sum(es_pose_array,mgweights_array,numgrids);
                }
		else {
                         if (bltzmn) {   //BCF
                            partitionFunc(bltzmn_Z,vdw_pose_array,es_pose_array,bltzmn_weight_array,numgrids,bltzmn_temp);
                            vdw_sum = bltzmn_sum(vdw_pose_array,bltzmn_weight_array,numgrids);
                            es_sum  = bltzmn_sum(es_pose_array,bltzmn_weight_array,numgrids);
                         }
                         else {
                            vdw_sum = sum(vdw_pose_array,numgrids);
                            es_sum = sum(es_pose_array,numgrids);
                         }           
                }
        }

        mol.current_score = vdw_sum * vdw_scale + es_sum * es_scale;
        mol.score_nrg = mol.current_score; // CS 06-06-16, save grid score to DOCKMol Object:774


        if ( fp_mol || fp_txt ) {

           //change from calculating all three metrics to calculating only what user specifies.  Yuchen 11/3/2016
           //CS 06-06-16, Save the unweighted FPS score as a separate DOCKMol Object
           if (use_cor){

               vdw_cor = correlation(vdw_pose_array,vdw_ref_array,numgrids);
               es_cor = correlation(es_pose_array,es_ref_array,numgrids);
               mol.score_fps = vdw_cor + es_cor;//CS
               mol.current_score = mol.current_score + vdw_cor_scale * vdw_cor + es_cor_scale * es_cor;
       
           }

           else if (use_euc){
 
               vdw_euc = euclidean(vdw_pose_array,vdw_ref_array,numgrids);
               es_euc = euclidean(es_pose_array,es_ref_array,numgrids);
               mol.score_fps = vdw_euc + es_euc;//CS
               mol.current_score = mol.current_score + vdw_euc_scale * vdw_euc + es_euc_scale * es_euc;

           }

           else if (use_norm){

               float norm_vdw_pose[numgrids];
               float norm_vdw_ref[numgrids];
               normalize(vdw_pose_array,numgrids,norm_vdw_pose);
               normalize(vdw_ref_array,numgrids,norm_vdw_ref);
               vdw_norm = euclidean(norm_vdw_pose,norm_vdw_ref,numgrids);
           
               float norm_es_pose[numgrids];
               float norm_es_ref[numgrids];
               normalize(es_pose_array,numgrids,norm_es_pose);
               normalize(es_ref_array,numgrids,norm_es_ref);
               es_norm = euclidean(norm_es_pose,norm_es_ref,numgrids);
               mol.score_fps =  vdw_norm + es_norm;//CS
               mol.current_score = mol.current_score + vdw_norm_scale * vdw_norm + es_norm_scale * es_norm;
            
            }

        }
         
        mol.current_data = output_score_summary(mol);
        delete[] vdw_pose_array; delete[] es_pose_array;
    return true;
}

// +++++++++++++++++++++++++++++++++++++++++
string
Multigrid_Energy_Score::output_score_summary(DOCKMol & mol)
{
    ostringstream text;

    if (use_score) {

        text << DELIMITER << setw(STRING_WIDTH) << "MultiGrid_Score:"
             << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl;
        text << DELIMITER << setw(STRING_WIDTH) << "MGS_vdw_energy:" 
             << setw(FLOAT_WIDTH) << fixed << vdw_sum << endl
             << DELIMITER << setw(STRING_WIDTH) << "MGS_es_energy:"  
             << setw(FLOAT_WIDTH) << fixed << es_sum  << endl
             << DELIMITER << setw(STRING_WIDTH) << "MGS_vdw+es_energy:" 
             << setw(FLOAT_WIDTH) << fixed << vdw_sum+es_sum  << endl;

        if ( fp_mol || fp_txt ) {

           //change from outputting all three metrics to outputting what user specifies.  Yuchen 11/3/2016
           if (use_cor){
               text << DELIMITER << setw(STRING_WIDTH) << "MGS_vdw_fps:"
                    << setw(FLOAT_WIDTH) << fixed << vdw_cor << endl;
               text << DELIMITER << setw(STRING_WIDTH) << "MGS_es_fps:"
                    << setw(FLOAT_WIDTH) << fixed << es_cor << endl;
               text << DELIMITER << setw(STRING_WIDTH) << "MGS_vdw+es_fps:"
                    << setw(FLOAT_WIDTH) << fixed << vdw_cor+es_cor << endl;
           }
           else if (use_euc){
               text << DELIMITER << setw(STRING_WIDTH) << "MGS_vdw_fps:"
                    << setw(FLOAT_WIDTH) << fixed << vdw_euc << endl;
               text << DELIMITER << setw(STRING_WIDTH) << "MGS_es_fps:"
                    << setw(FLOAT_WIDTH) << fixed << es_euc << endl;
               text << DELIMITER << setw(STRING_WIDTH) << "MGS_vdw+es_fps:"
                    << setw(FLOAT_WIDTH) << fixed << vdw_euc+es_euc << endl;
           }
           else if (use_norm){
               text << DELIMITER << setw(STRING_WIDTH) << "MGS_vdw_fps:"
                    << setw(FLOAT_WIDTH) << fixed << vdw_norm << endl;
               text << DELIMITER << setw(STRING_WIDTH) << "MGS_es_fps:"
                    << setw(FLOAT_WIDTH) << fixed << es_norm << endl;
               text << DELIMITER << setw(STRING_WIDTH) << "MGS_vdw+es_fps:"
                    << setw(FLOAT_WIDTH) << fixed << vdw_norm+es_norm << endl;
           }
        }
        // Compute lig internal energy with segments (sudipto 12-12-08)
        //float int_vdw_att, int_vdw_rep, int_es;
        //compute_ligand_internal_energy(mol, int_vdw_att, int_vdw_rep, int_es);

        if (use_internal_energy) 
            text << DELIMITER << setw(STRING_WIDTH) << "Internal_energy_repulsive:"
                 << setw(FLOAT_WIDTH) << fixed << mol.internal_energy << endl;
        //if ( fp_mol || fp_txt )
            //printarrays(vdw_pose_array,es_pose_array,numgrids,text);
        //if ( wt_txt )
            //printarrays(vdw_pose_array,es_pose_array,mgweights_array,numgrids,text);
        //if ( bltzmn )
            //printarrays(vdw_pose_array,es_pose_array,bltzmn_weight_array,numgrids,text);

        if ( ir_ensemble )
                text << DELIMITER << setw(STRING_WIDTH) << "Grid_num:"
                   << setw(FLOAT_WIDTH) << fixed << mol.grid_num << endl;

        // To get the correct internal energy, it MUST be calculated before calling the primary score

    }
    return text.str();
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++
//  THis is a boltman like weighting function:
//  E_i/(psudoT*abs(E_best)) is unitless like E_i/(k*T);
//  k is a small number. to get greater sepperation it is nessesary for psudoT to be small as well.
//  smaller k*T (or psudoT*abs(E_best)) results is greater sepperation will larger, less seperation.
//
void                     //BCF calculates Z and weighting, assigns both by reference
Multigrid_Energy_Score::partitionFunc(float &Z, float *vdw, float *es, float *wt, int size, float  pseudoT) {

  //cout << "Temp ==" << pseudoT << endl;
  //
  //exit(0);

  float E_best;
  //float pseudoT = t;
  float EMAX = 100;
  float E[size];
  float P[size];
  Z = 0;
  for(int h = 0; h < size; h++) {    // find E_best
      E[h] = vdw[h] + es[h];
      if (h == 0){
          E_best = E[h];
      }
      if (E[h] < E_best) {
          E_best = E[h];
      }
      if (E[h] > EMAX){
          E[h] = EMAX;
          //vdw[h] = 50.;
          //es[h]  = 50.;
      }
  }
  if (E_best < 1 && E_best > -1){
     E_best = 1; //no scaling if E_best is close to zero
  }

  for(int i = 0; i < size; i++){
    //Z += exp( -1 * ((vdw[i] + es[i]) / ( pseudoT * fabs(E_best)) ));  //0.593 = kT at room temp in kcal
    //Z += exp( -1 * ((E[i]) / ( pseudoT * fabs(E_best)) ));  //0.593 = kT at room temp in kcal
    if ( (-1.0 * E[i]) / (pseudoT * fabs(E_best)) > 10.0  ){
       P[i]  =  exp(10.0);
    }else{
       P[i] = exp( -1 * ((E[i]) / ( pseudoT * fabs(E_best)) ));  //0.593 = kT at room temp in kcal
    }
    Z += P[i];
  }
  if (Z > 100000000){
//    for (int i = 0; i < size; i++) {
//        cout << vdw[i] + es[i] << endl;
 //   }
    cout << "best energy: " << E_best << endl;
    cout << "partition function: " << Z << endl;

 }else if (Z < 0.0000001){
    cout << "partition function is to small: " << Z << endl;
    Z = 1;
 }


/*  if (Z < 0.000001) {    //if energy is high for all states Z will approach zero - don't want instability
     for(int j = 0; j < size; j++) {
       wt[j] = (1.0 / size);  //equal weighting, this pose will be pruned anyway
     }
  }   */
//  else {
     for(int j = 0; j < size; j++){
       //wt[j] = (1 / Z) * exp( -1 * ((vdw[j] + es[j]) / (pseudoT * fabs(E_best)) ));    // here N = 1 because we want fractional occupancy
       //wt[j] = (1 / Z) * exp( -1 * ((E[j]) / (pseudoT * fabs(E_best)) ));    // here N = 1 because we want fractional occupancy
       wt[j] = (1 / Z) * P[j];
      //  wt[j] = 0.5 * ( wt[j] + (1 / Z) * exp( -1 * ((vdw[j] + es[j]) / (kt) )));    // inertia, but previos wt is not for same system
     }
//  }
  return;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++
float
Multigrid_Energy_Score::bltzmn_sum(float *x, float *wt,int size) {

   float sumx = 0;

   for (int i=0;i<size;i++){
       sumx = sumx + wt[i] * x[i];
   }
   return sumx;
}

