#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <sstream>
#include "amber_typer.h"
#include "base_score.h"
#include "dockmol.h"
#include "grid.h" 
//#include "score_descriptor.h"
#include "score_hbond.h"
#include "utils.h"

#include <math.h>

using namespace std;

// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++

// +++++++++++++++++++++++++++++++++++++++++

Hbond_Energy_Score::Hbond_Energy_Score()
{
}

// +++++++++++++++++++++++++++++++++++++++++
Hbond_Energy_Score::~Hbond_Energy_Score()
{
}

void
Hbond_Energy_Score::close()
{
    // clear contance of c++ strings 
    receptor_filename.clear();
}

// +++++++++++++++++++++++++++++++++++++++++
void
Hbond_Energy_Score::input_parameters(Parameter_Reader & parm,
                                          bool & primary_score,
                                          bool & secondary_score)
{
    string          tmp;
   // string          tmp2;

    use_primary_score = false;
    use_secondary_score = false;


    // intailize parameters.

    cout << "\nHbond Energy Score Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    if (!primary_score) {
        tmp = parm.query_param("hbond_score_primary", "no", "yes no");
        if (tmp == "yes")
            use_primary_score = true;
        else
            use_primary_score = false;

        primary_score = use_primary_score;
    }

    if (!secondary_score) {
        //tmp = parm.query_param("hbond_score_secondary", "no", "yes no");
        //if (tmp == "yes")
            //use_secondary_score = true;
        //else
            use_secondary_score = false;

        secondary_score = use_secondary_score;
    }

    if (use_primary_score || use_secondary_score)
        use_score = true;
    else
        use_score = false;



    if (use_score) {
        receptor_filename =
            parm.query_param("hbond_score_rec_filename", "receptor.mol2");

/*
        att_exp = atoi(parm.query_param("hbond_score_att_exp", "6").c_str());
        if (att_exp <= 0) {
            cout <<
                "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }
        rep_exp = atoi(parm.query_param("hbond_score_rep_exp", "12").c_str());
        if (rep_exp <= 0) {
            cout <<
                "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }
        rep_radius_scale = atof(parm.query_param("hbond_score_rep_rad_scale", "1").c_str());
        if (rep_radius_scale <= 0.0) {
                cout <<
                    "ERROR:  Parameter must be a float greater than zero. Program will terminate."
                    << endl;
                exit(0);
        }

        //Ask whether we should use a distance dependent dieleectric

        tmp = parm.query_param("desc_use_distance_dependent_dielectric", "yes", "yes no");
        if (tmp == "yes")
            use_ddd = true;
        else
            use_ddd = false;

        //Use default dielectric = 4.0 for distance dependent dielectric
        //else use d=1.0 assuming gas phase

        char default_diel[4];
        if (use_ddd) strcpy(default_diel, "4.0");
        else strcpy(default_diel,"1.0");

        diel_screen =
            atof(parm.query_param("hbond_score_dielectric", default_diel).c_str());
        if (diel_screen <= 0.0) {
            cout <<
                "ERROR: Parameter must be a float greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }

        vdw_scale = atoi(parm.query_param("hbond_score_vdw_scale", "1").c_str());
        if (vdw_scale <= 0) {
            cout << "WARNING: Parameter should be an integer greater than zero."<< endl;
        }
        es_scale = atoi(parm.query_param("hbond_score_es_scale", "1").c_str());
        if (es_scale <= 0) {
           cout <<"WARNING: Parameter should be an integer greater than zero." << endl;
        }
        hb_scale = atoi(parm.query_param("hbond_score_hb_scale", "0").c_str());
        if (hb_scale>0){
           cout <<"WARNING: Parameter should be a number less than zero."<< endl;
        }
*/
        acc_don_rec_distance  = atoi(parm.query_param("hbond_acc_don_rec_shell_distance", "3.0").c_str()); // this is a parameter for ditermining the h-bond acceptors/doners on the receptor. (distance from ligand)
        threshold = atoi(parm.query_param("hbond_distance_threshold", "2.5").c_str()); // this is the distance used for defining an h-bond
        min_angle = atoi(parm.query_param("hbond_min_angle", "120").c_str()); // this is the minimum angle permited for defining an b-bond

        internal_scale = atoi(parm.query_param("hbond_score_internal_scale", "0").c_str());
    }
}

// +++++++++++++++++++++++++++++++++++++++++
void
Hbond_Energy_Score::initialize(AMBER_TYPER & typer)
{
    ifstream        rec_file;
    bool            read_vdw,
                    use_chem,
                    use_ph4,
                    use_volume;

    if (use_score) {

        init_vdw_energy(typer, att_exp, rep_exp);

        rec_file.open(receptor_filename.c_str());

        if (rec_file.fail()) {
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
        use_ph4 = false;
        use_volume = false;
        typer.prepare_molecule(receptor, read_vdw, use_chem, use_ph4, use_volume);

        hbond = rec_hb_don = rec_hb_acc = 
        lig_hb_don = lig_hb_acc =
        rec_lig_hb_don = rec_lig_hb_acc =
        lig_rec_hb_don = lig_rec_hb_acc = 0;  

    }
}

/*
*/
// +++++++++++++++++++++++++++++++++++++++++
// hbond_cal is defined in footprint code
// +++++++++++++++++++++++++++++++++++++++++
bool Hbond_Energy_Score::hbond_cal_func1(DOCKMol & mol_a,DOCKMol & mol_d,int i,int j,double & hbond_dist, double & hbond_angle)
{
     // distances between Hd --- A.  non-covalent
     //double threshold = 2.5;
     //double min_angle = 120;
     hbond_dist  = 0;
     hbond_angle = 0;
     double x_diff = mol_d.x[j] - mol_a.x[i];
     double y_diff = mol_d.y[j] - mol_a.y[i];
     double z_diff = mol_d.z[j] - mol_a.z[i];
     // eliminate any thing outside the bounding box
     if (   x_diff <= threshold && x_diff >= -threshold
         && y_diff <= threshold && y_diff >= -threshold
         && z_diff <= threshold && z_diff >= -threshold)
     { 
          hbond_dist = sqrt(pow(x_diff,2)+pow(y_diff,2) + pow(z_diff,2));
          // eliminate any thing outside the spheir with radius of 3.
          if (hbond_dist <= threshold) 
          {
               // distances beteewn Ad-Hd. covalent bond.
               double x_diff_covalent = mol_d.x[j]-mol_d.x[mol_d.acc_heavy_atomid[j]];
               double y_diff_covalent = mol_d.y[j]-mol_d.y[mol_d.acc_heavy_atomid[j]];
               double z_diff_covalent = mol_d.z[j]-mol_d.z[mol_d.acc_heavy_atomid[j]];
               double covalent_dist   = sqrt(pow(x_diff_covalent,2)+pow(y_diff_covalent,2) + pow(z_diff_covalent,2));
               //calculate the dot product
               double dot_coval_hbond =  x_diff*x_diff_covalent+y_diff*y_diff_covalent+z_diff*z_diff_covalent;
               hbond_angle = acos(dot_coval_hbond / (covalent_dist*hbond_dist))*180/PI; 
               if (hbond_angle <= 180 && hbond_angle >= min_angle){
                     return true;
               }
          }
     }
     return false;
}

// +++++++++++++++++++++++++++++++++++++++++
    // I want to do the following:
    // (0) Get the receptor residues close to the ligand.  
    //     (we may want to define the pocket like for fooprint score?).  
    // (0) ID which hbond aceptor or donars are satisfied 
    //     Generate an array of hydrogen bond aceptor 
    //     Generate an array of hydrogen bond aceptor
    //
    // 
    // (1) caluclate intramolecuar hbonds for receptor.
    //     (a)only look at residues close to the ligand (frist shell) 
    //     and those close to the frist shell (second shell)
    // 
    // (2) calculate intramolecuar hbonds for the ligand.
    // (3) calculate intermolecular hbonds for ligand-receptor complex. 
    //
// +++++++++++++++++++++++++++++++++++++++++
bool
Hbond_Energy_Score::compute_hbonds(DOCKMol & mol)
{
    cout << "In Hbond_Energy_Score::compute_hbonds(DOCKMol & mol)" << endl;
    double          vdw_val,
                    es_val,
                    total;
    int             i,
                    j;
    double          dist;
    stringstream    hbtxt; 
    //ostringstream    hbtxt; 
    //ostream   hbtxt = cout; 
    //these three variables are used to store energies for atom pairs;

    total = vdw_val = es_val = 0.0;
    vector <int> rec_acceptors, rec_donators, rec_acceptors_close, rec_donators_close, lig_acceptors, lig_donators;
    vector <bool> rec_sat_acc, rec_sat_don, lig_sat_acc, lig_sat_don, lig_rec_sat_don, lig_rec_sat_acc, rec_lig_sat_don, rec_lig_sat_acc; 
    
        hbond = rec_hb_don = rec_hb_acc = 
        lig_hb_don = lig_hb_acc =
        rec_lig_hb_don = rec_lig_hb_acc =
        lig_rec_hb_don = lig_rec_hb_acc = 0;  
    // vectors of boolians 
    // len should be the number of 
    // if acceptor or donar is satified then true if not then false. 

    //mol.num_residues_in_rec = receptor.num_residues;

    // start debuging
    // receptor
    int tempcountdonators  =0;
    int tempcountacceptors =0;
    for (j = 0; j < receptor.num_atoms; j++) {
             if (receptor.flag_donator[j]){
                  tempcountdonators++;
             }
             if (receptor.flag_acceptor[j] ){
                  tempcountacceptors++;
             }
    }
    cout << "receptor: number of donators  = " << tempcountdonators<<endl;
    cout << "receptor: number of acceptors = " << tempcountacceptors<<endl;
    // ligands
    tempcountdonators  =0;
    tempcountacceptors =0;
    for (j = 0; j < mol.num_atoms; j++) {
             if (mol.flag_donator[j]){
                  tempcountdonators++;
             }
             if (mol.flag_acceptor[j] ){
                  tempcountacceptors++;
             }
    }
    cout << "ligand: number of donators  = " << tempcountdonators<<endl;
    cout << "ligand: number of acceptors = " << tempcountacceptors<<endl;
    // end debuging
    
    //cout << "I AM HERE" << endl;
    if (use_score == 1) {
         // cout << "I AM HERE" << endl;
         // get list of receptor doners or acceptors close to ligand. 
         for (j = 0; j < receptor.num_atoms; j++) {
             if (receptor.flag_donator[j] || receptor.flag_acceptor[j] ){
             
                  for (i = 0; i < mol.num_atoms; i++) {
                       dist =
                            sqrt(pow((mol.x[i] - receptor.x[j]), 2) +
                                 pow((mol.y[i] - receptor.y[j]), 2) + 
                                 pow((mol.z[i] - receptor.z[j]), 2));
                       if (dist < acc_don_rec_distance){
                           if (receptor.flag_donator[j])  rec_donators_close.push_back(j);
                           if (receptor.flag_acceptor[j]) rec_acceptors_close.push_back(j);
                           break; // move do next donator or acceptor on receptor
                       } // check the dist sqared
                  } // loop over ligand atoms         

             } // if acceptor of doner
             if (receptor.flag_donator[j])  rec_donators.push_back(j);
             if (receptor.flag_acceptor[j]) rec_acceptors.push_back(j);
         } // loop over receptor atoms

         // cout << "I AM HERE" << endl;
         // get list of ligands doners or acceptors
         for (i = 0; i < mol.num_atoms; i++) {
                  if (mol.flag_donator[i]) lig_donators.push_back(i);
                  if (mol.flag_acceptor[i]) lig_acceptors.push_back(i);
         } // loop over ligand atoms

         cout << "receptor: number of donators  = " << rec_donators.size() << endl;
         cout << "receptor: number of acceptors = " << rec_acceptors.size() << endl;
         cout << "receptor shell: number of donators  = " << rec_donators_close.size() << endl;
         cout << "receptor shell: number of acceptors = " << rec_acceptors_close.size() << endl;
         cout << "ligand: number of donators  = " << lig_donators.size() << endl;
         cout << "ligand: number of acceptors = " << lig_acceptors.size() << endl;
         // cout << "I AM HERE" << endl;
         // cacluate rec - rec interactions.
         rec_sat_acc.resize(rec_acceptors_close.size(),false);
         rec_sat_don.resize(rec_donators_close.size(),false);
         for (j = 0; j < rec_donators_close.size(); j++) {
              for (i = 0; i < rec_acceptors.size(); i++) {
                   double hbond_dist,hbond_angle;
                   rec_sat_don[j] = (rec_sat_don[j] || hbond_cal_func1(receptor,receptor,rec_acceptors[i],rec_donators_close[j],hbond_dist, hbond_angle));
              }         
         }
         for (j = 0; j < rec_acceptors_close.size(); j++) {
              for (i = 0; i < rec_donators.size(); i++) {
                   double hbond_dist,hbond_angle;
                   rec_sat_acc[j] = (rec_sat_acc[j] || hbond_cal_func1(receptor,receptor,rec_acceptors_close[j],rec_donators[i],hbond_dist, hbond_angle));
              }
         } 
         // cout << "I AM HERE" << endl;
         // cout << "cacluate lig - lig interactions." << endl;
         // cacluate lig - lig interactions.
         lig_sat_acc.resize(lig_acceptors.size(),false);
         lig_sat_don.resize(lig_donators.size(),false);
         for (j = 0; j < lig_donators.size(); j++) {
              for (i = 0; i < lig_acceptors.size(); i++) {
                   double hbond_dist,hbond_angle;
                   lig_sat_don[j] = (lig_sat_don[j] || hbond_cal_func1(mol,mol,lig_acceptors[i],lig_donators[j],hbond_dist, hbond_angle));
              }
         }
         for (j = 0; j < lig_acceptors.size(); j++) {
              for (i = 0; i < lig_donators.size(); i++) {
                   double hbond_dist,hbond_angle;
                   lig_sat_acc[j] = (lig_sat_acc[j] || hbond_cal_func1(mol,mol,lig_acceptors[j],lig_donators[i],hbond_dist, hbond_angle));
              }
         }

         // cout << "I AM HERE" << endl;
         // cacluate rec - lig interactions.

         rec_lig_sat_acc.resize(rec_acceptors_close.size(),false);
         rec_lig_sat_don.resize(rec_donators_close.size(),false);
         for (j = 0; j < rec_donators_close.size(); j++) {
              for (i = 0; i < lig_acceptors.size(); i++) {
                   double hbond_dist,hbond_angle;
                   rec_lig_sat_don[j] = (rec_lig_sat_don[j] || hbond_cal_func1(mol,receptor,lig_acceptors[i],rec_donators_close[j],hbond_dist, hbond_angle));
              }
         }
         for (j = 0; j < rec_acceptors_close.size(); j++) {
              for (i = 0; i < lig_donators.size(); i++) {
                   double hbond_dist,hbond_angle;
                   rec_lig_sat_acc[j] = (rec_lig_sat_acc[j] || hbond_cal_func1(receptor,mol,rec_acceptors_close[j],lig_donators[i],hbond_dist, hbond_angle));
              }
         }

         // cout << "I AM HERE" << endl;
         // cacluate lig - rec interactions.
         lig_rec_sat_acc.resize(lig_acceptors.size(),false);
         lig_rec_sat_don.resize(lig_donators.size(),false);
         for (j = 0; j < lig_donators.size(); j++) {
              for (i = 0; i < rec_acceptors_close.size(); i++) {
                   double hbond_dist,hbond_angle;
                   //lig_rec_sat_don[j] = (lig_rec_sat_don[j] || hbond_cal(receptor,mol,rec_acceptors_close[i],lig_donators[j],hbond_dist, hbond_angle))
                   bool hflag = hbond_cal_func1(receptor,mol,rec_acceptors_close[i],lig_donators[j],hbond_dist, hbond_angle);
                   lig_rec_sat_don[j] = (lig_rec_sat_don[j] || hflag);
                   //cout << "hflag = " << hflag << endl;
                   if (hflag) {
                              int ind_don = lig_donators[j];
                              int ind_acc = rec_acceptors_close[i];
                              hbtxt << "ligand donator and receptor acceptor" << endl;
                              hbtxt << "hbond detected: " << mol.subst_names[ind_don]<<" ";
                              hbtxt << "(" << mol.atom_residue_numbers[ind_don] << ")" << " ";
                              hbtxt << mol.atom_types[mol.acc_heavy_atomid[ind_don]]<<" ";
                              hbtxt << mol.atom_types[ind_don]<<"---"<<receptor.subst_names[ind_acc]<<" ";
                              hbtxt << "(" << receptor.atom_residue_numbers[ind_acc] << ")" << " ";
                              hbtxt << receptor.atom_types[ind_acc] << " " << hbond_dist <<" " << hbond_angle  << endl;
                              //cout << hbtxt.str();
                              hbond ++;
                  }
              }
         }
         // cout << "I AM HERE" << endl;
         for (j = 0; j < lig_acceptors.size(); j++) {
              for (i = 0; i < rec_donators_close.size(); i++) {
                   double hbond_dist,hbond_angle;
                   //lig_rec_sat_acc[j] = (lig_rec_sat_acc[j] || hbond_cal(mol,receptor,lig_acceptors[j],rec_donators_close[i],hbond_dist, hbond_angle))
                   bool hflag =  hbond_cal_func1(mol,receptor,lig_acceptors[j],rec_donators_close[i],hbond_dist, hbond_angle);
                   lig_rec_sat_acc[j] = (lig_rec_sat_acc[j] || hflag);

                   if (hflag) {
                              int ind_don = rec_donators_close[i];
                              int ind_acc = lig_acceptors[j]; 
                              hbtxt << "recptor donator and ligand acceptor" << endl;
                              hbtxt << "hbond detected: " << receptor.subst_names[ind_don]<<" ";
                              hbtxt << "(" << receptor.atom_residue_numbers[ind_don] << ")" << " ";
                              hbtxt << receptor.atom_types[receptor.acc_heavy_atomid[ind_don]]<<" ";
                              hbtxt << receptor.atom_types[ind_don]<<"---"<< mol.subst_names[ind_acc]<<" ";
                              hbtxt << "(" << mol.atom_residue_numbers[ind_acc] << ")" << " ";
                              hbtxt << mol.atom_types[ind_acc] << " " << hbond_dist <<" " << hbond_angle  << endl;
                              //cout << hbtxt.str();
                              hbond++;
                  }


              }
         }
        
        
         //cout << "I AM HERE" << endl;

         hbtxt << "hbond receptor acceptor: \n";

         rec_rec_hb_acc = rec_lig_hb_acc = rec_hb_acc = rec_acceptors_close.size();

         for (j = 0; j < rec_acceptors_close.size(); j++) {
              int ind = rec_acceptors_close[j];
              hbtxt << ind << " " <<receptor.atom_residue_numbers[ind] << " " << receptor.atom_types[ind] << " " <<rec_sat_acc[j] << " " << rec_lig_sat_acc[j] << endl;
              if (rec_sat_acc[j])  rec_rec_hb_acc -- ;
              if (rec_lig_sat_acc[j])  rec_lig_hb_acc -- ;
              if (rec_lig_sat_acc[j] || rec_sat_acc[j]) rec_hb_acc -- ;
         }
         hbtxt << "hbond receptor donator: \n";
         rec_rec_hb_don = rec_lig_hb_don = rec_hb_don = rec_donators_close.size();
         for (j = 0; j < rec_donators_close.size(); j++) {
              int ind = rec_donators_close[j];
              hbtxt << ind << " " << receptor.atom_types[receptor.acc_heavy_atomid[ind]] <<" " << receptor.atom_residue_numbers[ind] << " " << receptor.atom_types[ind] << " " <<rec_sat_don[j] << " " << rec_lig_sat_don[j] << endl;
              if (rec_sat_don[j])  rec_rec_hb_don --;
              if (rec_lig_sat_don[j])  rec_lig_hb_don -- ;
              if (rec_lig_sat_don[j] || rec_sat_don[j])  rec_hb_don -- ;

         }
         hbtxt << "hbond ligand acceptor: \n";
         lig_lig_hb_acc = lig_rec_hb_acc = lig_hb_acc = lig_acceptors.size();
         for (j = 0; j < lig_acceptors.size(); j++) {
              int ind = lig_acceptors[j];
              hbtxt << ind << " " <<mol.atom_residue_numbers[ind] << " " << mol.atom_types[ind] << " " <<lig_sat_acc[j] << " " << lig_rec_sat_acc[j] << endl;
              if (lig_sat_acc[j])  lig_lig_hb_acc --;
              if (lig_rec_sat_acc[j])  lig_rec_hb_acc -- ;
              if (lig_rec_sat_acc[j] || lig_sat_acc[j])  lig_hb_acc -- ;

         }
         hbtxt << "hbond ligand donator: \n";
         lig_lig_hb_don = lig_rec_hb_don = lig_hb_don = lig_donators.size();
         for (j = 0; j < lig_donators.size(); j++) {
              int ind = lig_donators[j];
              hbtxt << ind << " "<< mol.atom_types[mol.acc_heavy_atomid[ind]] <<" " << mol.atom_residue_numbers[ind] << " " << mol.atom_types[ind] << " " <<lig_sat_don[j] << " " << lig_rec_sat_don[j] << endl;
              if (lig_sat_don[j])  lig_lig_hb_don --;
              if (lig_rec_sat_don[j])  lig_rec_hb_don -- ;
              if (lig_rec_sat_don[j] || lig_sat_don[j])  lig_hb_don -- ;
             
         }      
        mol.hbond_text_data = hbtxt.str();
        rec_acceptors_close.clear();
        rec_donators_close.clear();
        rec_acceptors.clear();
        rec_donators.clear();
        lig_acceptors.clear();
        lig_donators.clear();
        //cout << mol.hbond_text_data;
        mol.hbond_text_data = hbtxt.str();
        hbtxt.clear();
        cout << "Out Hbond_Energy_Score::compute_hbonds(DOCKMol & mol)" << endl;
        return true;
     }
     return false;
}


// +++++++++++++++++++++++++++++++++++++++++
bool
Hbond_Energy_Score::compute_score(DOCKMol & mol)
{
      compute_hbonds(mol);
      mol.current_score = hbond;
      mol.current_data = output_score_summary(mol);
      return true;
}

// +++++++++++++++++++++++++++++++++++++++++
string
Hbond_Energy_Score::output_score_summary(DOCKMol & mol)
{
    ostringstream text;
//    cout << "inside output_score_summary" << endl;
    if (use_score) {
        text << DELIMITER << setw(STRING_WIDTH) << "Hbond_Score:"
             << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl
             << DELIMITER << setw(STRING_WIDTH) << "unsatified_rec_donator:"
             << setw(FLOAT_WIDTH) << fixed << rec_hb_don << endl
             << DELIMITER << setw(STRING_WIDTH) << "unsatified_rec_acceptor:"
             << setw(FLOAT_WIDTH) << fixed << rec_hb_acc << endl
             << DELIMITER << setw(STRING_WIDTH) << "unsatified_lig_donator:"
             << setw(FLOAT_WIDTH) << fixed << lig_hb_don << endl
             << DELIMITER << setw(STRING_WIDTH) << "unsatified_lig_acceptor:"
             << setw(FLOAT_WIDTH) << fixed << lig_hb_acc << endl
             << DELIMITER << setw(STRING_WIDTH) << "hbond:"
             << setw(FLOAT_WIDTH) << fixed << hbond  << endl ;
        if (use_internal_energy)
            text << DELIMITER << setw(STRING_WIDTH) << "Int_energy:"
                 << setw(FLOAT_WIDTH) << fixed << mol.internal_energy << endl;
        text << mol.hbond_text_data << endl;
    }
    return text.str();
}


// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++
