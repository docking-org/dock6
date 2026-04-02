#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include "amber_typer.h"
#include "base_score.h"
#include "dockmol.h"
#include "grid.h" 
#include "score_volume.h"
#include "utils.h"

#include <math.h>

using namespace std;



// +++++++++++++++++++++++++++++++++++++++++
Volume_Score::Volume_Score()
{

}



// +++++++++++++++++++++++++++++++++++++++++
Volume_Score::~Volume_Score()
{

}



// +++++++++++++++++++++++++++++++++++++++++
void
Volume_Score::input_parameters(Parameter_Reader & parm, bool & primary_score,
                               bool & secondary_score)
{
    use_primary_score = false;
    use_secondary_score = false;

    cout << "\nVolume Score Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    string          tmp;

    if (!primary_score) {
        tmp = parm.query_param("volume_score_primary", "no", "yes no");
        use_primary_score = tmp == "yes";
        primary_score = use_primary_score;
    }

    if (!secondary_score) {
        //tmp = parm.query_param("volume_score_secondary", "no", "yes no");
        tmp = "no";
	use_secondary_score = tmp == "yes";
        secondary_score = use_secondary_score;
    }

    if (use_primary_score || use_secondary_score)
        use_score = true;
    else
        use_score = false;

    if (use_score) {
        input_parameters_main(parm, "volume_score");
    }
}



// +++++++++++++++++++++++++++++++++++++++++
void
Volume_Score::input_parameters_main(Parameter_Reader & parm, string parm_head)
{
    volume_ref_file = parm.query_param(parm_head + "_reference_mol2_filename",
                                       parm_head + "_reference.mol2");
    volume_compare_type = parm.query_param(parm_head + "_overlap_compute_method",
                                           "analytical", "analytical grid");
    if (volume_compare_type == "grid") {
        steps = atoi(parm.query_param(parm_head + "_number_of_layers_in_one_dimension",
                                      "30").c_str()); 
        separation = atof(parm.query_param(parm_head + "_maximum_grid_point_separation",
                                           "0.3").c_str());
        tmp = parm.query_param(parm_head + "_write_points_cloud", "no", "yes no");
        output_cloud = (tmp == "yes");
    }
    //define the method

    /*
    rep_radius_scale = atof(parm.query_param("grid_score_rep_rad_scale", "1").c_str());
    if (rep_radius_scale <= 0.0) {
            cout << "ERROR:  Parameter must be a float greater than zero."
                    " Program will terminate."
                 << endl;
            exit(0);
    }
    vdw_scale = atof(parm.query_param("grid_score_vdw_scale", "1").c_str());
    if (vdw_scale <= 0.0) {
        bool off = parm.query_param( "grid_score_turn_off_vdw", "yes",
                                         "yes no") == "yes";
        if (!off) {
            cout << "ERROR:  Parameter must be a float greater than zero."
                    " Program will terminate."
                 << endl;
            exit(0);
        }
    }
    es_scale = atof(parm.query_param("grid_score_es_scale", "1").c_str());
    if (es_scale <= 0.0) {
        bool off = parm.query_param( "grid_score_turn_off_es", "yes",
                                     "yes no") == "yes";
        if (!off) {
            cout << "ERROR:  Parameter must be a float greater than zero."
                    " Program will terminate."
                 << endl;
            exit(0);
        }
    }
    grid_file_name = parm.query_param( "grid_score_grid_prefix", "grid" );
    */

} // end Volume_Score::input_parameters()



// +++++++++++++++++++++++++++++++++++++++++
void
Volume_Score::initialize(AMBER_TYPER & typer)
{

    ifstream        volume_ref;

    if (use_score) {

        volume_ref.open(volume_ref_file.c_str());

        if (volume_ref.fail()) {
            cout << "Error Opening Reference File!" << endl;
            exit(0);
        }

        Read_Mol2(volume_ref_mol, volume_ref, false, false, false);
        cout <<" Volume reference mol2 file was read in successfully" <<endl;

        volume_ref.close();

        typer.prepare_molecule(volume_ref_mol, true, false, false, true);

    }

    int         i,
                j;
    float       pi = 3.14,
                d_sqr,
                overlap_temp;
    overlap_ref_ref = 0;
    overlap_ref_ref_hvy = 0;
    overlap_ref_ref_pho = 0;
    overlap_ref_ref_phi = 0;
    overlap_ref_ref_neg = 0;
    overlap_ref_ref_pos = 0;


    for (i = 0; i < volume_ref_mol.num_atoms; i++){
        for (j = 0; j < volume_ref_mol.num_atoms; j++){
            d_sqr = (volume_ref_mol.x[i] - volume_ref_mol.x[j]) * (volume_ref_mol.x[i] - volume_ref_mol.x[j]) + (volume_ref_mol.y[i] - volume_ref_mol.y[j]) * (volume_ref_mol.y[i] - volume_ref_mol.y[j]) + (volume_ref_mol.z[i] - volume_ref_mol.z[j]) * (volume_ref_mol.z[i] - volume_ref_mol.z[j]);
            if (d_sqr < (volume_ref_mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j]) * (volume_ref_mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j])) {
                if (d_sqr != 0)
                    overlap_temp = pi / 12 * (volume_ref_mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j] - sqrt(d_sqr)) * (volume_ref_mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j] - sqrt(d_sqr)) * (sqrt(d_sqr) + 2 * (volume_ref_mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j]) - 3 / sqrt(d_sqr) * (volume_ref_mol.amber_at_radius[i] - volume_ref_mol.amber_at_radius[j]) * (volume_ref_mol.amber_at_radius[i] - volume_ref_mol.amber_at_radius[j]));
                else
                    overlap_temp = pi / 12 * 2 * (volume_ref_mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j]) * (volume_ref_mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j]) * (volume_ref_mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j]);
                overlap_ref_ref = overlap_ref_ref + overlap_temp;
                if (volume_ref_mol.atom_types[i] != "H" && volume_ref_mol.atom_types[j] != "H")
                    overlap_ref_ref_hvy = overlap_ref_ref_hvy + overlap_temp;
                if (volume_ref_mol.chem_types[i] == "hydrophobic" && volume_ref_mol.chem_types[j] == "hydrophobic")
                    overlap_ref_ref_pho = overlap_ref_ref_pho + overlap_temp;
                if ((volume_ref_mol.chem_types[i] == "donar" || volume_ref_mol.chem_types[i] == "acceptor" || volume_ref_mol.chem_types[i] == "polar") && (volume_ref_mol.chem_types[j] == "donar" || volume_ref_mol.chem_types[j] == "acceptor" || volume_ref_mol.chem_types[j] == "polar"))
                    overlap_ref_ref_phi = overlap_ref_ref_phi + overlap_temp;
                if (volume_ref_mol.charges[i] < 0 && volume_ref_mol.charges[j] < 0)
                    overlap_ref_ref_neg = overlap_ref_ref_neg + overlap_temp;
                if (volume_ref_mol.charges[i] > 0 && volume_ref_mol.charges[j] > 0)
                    overlap_ref_ref_pos = overlap_ref_ref_pos + overlap_temp;

            }
        }
    }



    return;

    /*
    if (use_score) {
        cout << "Initializing Grid Score Routines..." << endl;
        energy_grid = Energy_Grid :: get_instance(grid_file_name);
        init_vdw_energy(typer, energy_grid->att_exp, energy_grid->rep_exp);
    } 
    */

} // end Volume_Score::initialize()

/*
// +++++++++++++++++++++++++++++++++++++++++
Cloud_Struct::Cloud_Struct()
{
}

// +++++++++++++++++++++++++++++++++++++++++
Cloud_Struct::~Cloud_Struct()
{
   this->clear();
}

// +++++++++++++++++++++++++++++++++++++++++

// +++++++++++++++++++++++++++++++++++++++++
void
Cloud_Struct::initialize(){
//initialize the points cloud for the output mol2 file
    new vector <float>      x;
    new vector <float>      y;
    new vector <float>      z;
}

*/

// +++++++++++++++++++++++++++++++++++++++++    

bool
Volume_Score::Grid_method(DOCKMol & mol)
{
    int             i,j,k;
    int             count = 0;
    float           temp = 0;
    float           d1,d2;
    float           maxX,maxY,maxZ,minX,minY,minZ,length,depth,height,xStep,yStep,zStep,xQ,yQ,zQ;
    float           positive_d1,negative_d1,positive_d2,negative_d2,pho_d1,pho_d2,phi_d1,phi_d2;
    float           in_mol = 0,
                    in_mol_hvy = 0,
	            in_mol_pos = 0,
	            in_mol_neg = 0,
	            in_mol_pho = 0,
	            in_mol_phi = 0,
	            in_ref = 0,
                    in_ref_hvy = 0,
	            in_ref_pos = 0,
	            in_ref_neg = 0,
	            in_ref_pho = 0,
	            in_ref_phi = 0,
	            in_both = 0,
                    in_both_hvy = 0,
	            in_both_pos = 0,
	            in_both_neg = 0,
	            in_both_pho = 0,
	            in_both_phi = 0;
    bool            in_mol_tag,
                    in_mol_hvy_tag,
	            in_mol_pos_tag,
                    in_mol_neg_tag,
                    in_mol_pho_tag,
                    in_mol_phi_tag,
                    in_ref_tag,
                    in_ref_hvy_tag,
                    in_ref_pos_tag,
                    in_ref_neg_tag,
                    in_ref_pho_tag,
                    in_ref_phi_tag;
    float           number,
                    total;
    Cloud_Struct    cloud;
    Cloud_Struct    charge_cloud;
    Cloud_Struct    hydrophobicity_cloud;
    maxX = 0;
    maxY = 0;
    maxZ = 0;
    minX = 100;
    minY = 100;
    minZ = 100;

    cout <<"Making the box..." <<endl;

    for (i = 0; i < mol.num_atoms; i++) {
        if (mol.atom_active_flags[i]) {
            if (maxX < mol.x[i] + mol.amber_at_radius[i]) {
                maxX = mol.x[i] + mol.amber_at_radius[i];
            }
            if (maxY < mol.y[i] + mol.amber_at_radius[i]) {
                maxY = mol.y[i] + mol.amber_at_radius[i];
            }
            if (maxZ < mol.z[i] + mol.amber_at_radius[i]) {
                maxZ = mol.z[i] + mol.amber_at_radius[i];
            }
            if (minX > mol.x[i] - mol.amber_at_radius[i] - 1) {
                minX = mol.x[i] - mol.amber_at_radius[i] - 1;
            }
            if (minY > mol.y[i] - mol.amber_at_radius[i] - 1) {
                minY = mol.y[i] - mol.amber_at_radius[i] - 1;
            }
            if (minZ > mol.z[i] - mol.amber_at_radius[i] - 1) {
                minZ = mol.z[i] - mol.amber_at_radius[i] - 1;
            }
        }
    }

    for (i = 0; i < volume_ref_mol.num_atoms; i++) {
        if (maxX < volume_ref_mol.x[i] + volume_ref_mol.amber_at_radius[i]) {
            maxX = volume_ref_mol.x[i] + volume_ref_mol.amber_at_radius[i];
        }
        if (maxY < volume_ref_mol.y[i] + volume_ref_mol.amber_at_radius[i]) {
            maxY = volume_ref_mol.y[i] + volume_ref_mol.amber_at_radius[i];
        }
        if (maxZ < volume_ref_mol.z[i] + volume_ref_mol.amber_at_radius[i]) {
            maxZ = volume_ref_mol.z[i] + volume_ref_mol.amber_at_radius[i];
        }
        if (minX > volume_ref_mol.x[i] - volume_ref_mol.amber_at_radius[i] - 1) {
            minX = volume_ref_mol.x[i] - volume_ref_mol.amber_at_radius[i] - 1;
        }
        if (minY > volume_ref_mol.y[i] - volume_ref_mol.amber_at_radius[i] - 1) {
            minY = volume_ref_mol.y[i] - volume_ref_mol.amber_at_radius[i] - 1;
        }
        if (minZ > volume_ref_mol.z[i] - volume_ref_mol.amber_at_radius[i] - 1) {
            minZ = volume_ref_mol.z[i] - volume_ref_mol.amber_at_radius[i] - 1;
        }
    }

    length = maxX - minX;
    depth = maxY - minY;
    height = maxZ - minZ;

    //cout << length << endl;


    xStep = min((length/steps),separation);
    yStep = min((depth/steps),separation);
    zStep = min((height/steps),separation); 

    //cout << volume_ref_mol.chem_types[3] << endl;

    xQ = minX;
    while (xQ <= maxX){
        yQ = minY;
	while (yQ <= maxY){
	    zQ = minZ;
            while (zQ <= maxZ){
                in_mol_tag = false;
                in_mol_hvy_tag = false;
                in_mol_pos_tag = false;
                in_mol_neg_tag = false;
                in_mol_pho_tag = false;
                in_mol_phi_tag = false;
                in_ref_tag = false;
                in_ref_hvy_tag = false;
                in_ref_pos_tag = false;
                in_ref_neg_tag = false;
                in_ref_pho_tag = false;
                in_ref_phi_tag = false;
                positive_d1 = 100;
                negative_d1 = 100;
                positive_d2 = 100;
                negative_d2 = 100;
                pho_d1 = 100;
                pho_d2 = 100;
                phi_d1 = 100;
                phi_d2 = 100;
                
                for (j = 0; j < mol.num_atoms; j++){
                    if (mol.atom_active_flags[j]) {
                        d2 = (mol.x[j] - xQ) * (mol.x[j] - xQ) + (mol.y[j] - yQ) * (mol.y[j] - yQ) + (mol.z[j] - zQ) * (mol.z[j] - zQ);
                        //cout<<d2<<endl;
                        //cout<<mol.amber_at_radius[j]<<endl;
                        if (d2 <= mol.amber_at_radius[j] * mol.amber_at_radius[j]){
                            in_mol_tag = true;
                            if (mol.atom_types[j] != "H"){
                                in_mol_hvy_tag = true;
                                if (mol.chem_types[j] == "hydrophobic"){
                                    if (d2 < pho_d2) pho_d2 = d2;
                                    if (pho_d2 < phi_d2) in_mol_pho_tag = true;
                                    else in_mol_pho_tag = false;
                                }
                                else if (mol.chem_types[j] == "donor" || mol.chem_types[j] == "acceptor" || mol.chem_types[j] == "polar"){
                                    if (d2 < phi_d2) phi_d2 = d2;
                                    if (phi_d2 < pho_d2) in_mol_phi_tag = true;  
                                    else in_mol_phi_tag = false;
                                }
                            }
                            if (mol.charges[j] > 0) {
                                if (d2 < positive_d2) positive_d2 = d2;
                                if (positive_d2 < negative_d2) in_mol_pos_tag = true;
                                else in_mol_pos_tag = false;
                            }
                            else {
                                if (d2 < negative_d2) negative_d2 = d2;
                                if (negative_d2 < positive_d2) in_mol_neg_tag = true; 
                                else in_mol_neg_tag = false;
                            }
                        }
                    }
                }
                for (i = 0; i < volume_ref_mol.num_atoms; i++){
                    d1 = (volume_ref_mol.x[i] - xQ) * (volume_ref_mol.x[i] - xQ) + (volume_ref_mol.y[i] - yQ) * (volume_ref_mol.y[i] - yQ) + (volume_ref_mol.z[i] - zQ) * (volume_ref_mol.z[i] - zQ);
                    //cout<<d1<<endl;
                    //cout<<volume_ref_mol.amber_at_radius[i]<<endl;
                    if (d1 <= volume_ref_mol.amber_at_radius[i] * volume_ref_mol.amber_at_radius[i]){
                        in_ref_tag = true;
                        if (volume_ref_mol.atom_types[i] != "H"){
                            in_ref_hvy_tag = true;
                            if (volume_ref_mol.chem_types[i] == "hydrophobic") {
                                if (d1 < pho_d1) pho_d1 = d1;
                                if (pho_d1 < phi_d1) in_ref_pho_tag = true;
                                else in_ref_pho_tag = false;
                            }
                            else if (volume_ref_mol.chem_types[i] == "donot" || volume_ref_mol.chem_types[i] == "acceptor" || volume_ref_mol.chem_types[i] == "polar") { 
                                if (d1 < phi_d1) phi_d1 = d1;
                                if (phi_d1 < pho_d1) in_ref_phi_tag = true;
                                else in_ref_phi_tag = false;
                            }
                        }
                        if (volume_ref_mol.charges[i] > 0) {
                            if (d1 < positive_d1) positive_d1 = d1;
                            if (positive_d1 < negative_d1) in_ref_pos_tag = true;
                            else in_ref_pos_tag = false;
                        }
                        else {
                            if (d1 < negative_d1) negative_d1 = d1;
                            if (negative_d1 < positive_d1) in_ref_neg_tag = true;
                            else in_ref_neg_tag = false;
                        }
                    }
                }
                if (in_mol_tag) {
                    in_mol++;
                }
                if (in_mol_hvy_tag) {
                    in_mol_hvy += 1;
                }
                if (in_mol_pho_tag) {
                    in_mol_pho += 1;
                }
                if (in_mol_phi_tag) {
                    in_mol_phi += 1;
                }
                if (in_mol_pos_tag) {
                    in_mol_pos += 1;
                }
                if (in_mol_neg_tag) {
                    in_mol_neg += 1;
                }
                if (in_ref_tag) {
                    in_ref++;
                }
                if (in_ref_hvy_tag) {
                    in_ref_hvy += 1;
                }
                if (in_ref_pho_tag) {
                    in_ref_pho += 1;
                }
                if (in_ref_phi_tag) {
                    in_ref_phi += 1;
                }
                if (in_ref_pos_tag) {
                    in_ref_pos += 1;
                }
                if (in_ref_neg_tag) {
                    in_ref_neg += 1;
                }
                if (in_mol_tag && in_ref_tag) {
                    in_both += 1;
                    if (output_cloud == true){
                        cloud.x.push_back(xQ);
                        cloud.y.push_back(yQ);
                        cloud.z.push_back(zQ);
                        cloud.type.push_back("N");
                    }
                }
                if (in_mol_hvy_tag && in_ref_hvy_tag) {
                    in_both_hvy += 1;
                }
                if (in_mol_pho_tag && in_ref_pho_tag) {
                    in_both_pho += 1;
                    if(output_cloud == true){
                        hydrophobicity_cloud.x.push_back(xQ);
                        hydrophobicity_cloud.y.push_back(yQ);
                        hydrophobicity_cloud.z.push_back(zQ);
                        hydrophobicity_cloud.type.push_back("C");
                    }
                }
                if (in_mol_phi_tag && in_ref_phi_tag) {
                    in_both_phi += 1;
                    if(output_cloud == true){
                        hydrophobicity_cloud.x.push_back(xQ);
                        hydrophobicity_cloud.y.push_back(yQ);
                        hydrophobicity_cloud.z.push_back(zQ);
                        hydrophobicity_cloud.type.push_back("S");
                    }
                }
                if (in_mol_pos_tag && in_ref_pos_tag) {
                    in_both_pos += 1;
                    if (output_cloud == true){
                        charge_cloud.x.push_back(xQ);
                        charge_cloud.y.push_back(yQ);
                        charge_cloud.z.push_back(zQ);
                        charge_cloud.type.push_back("H");
                    }
                }
                if (in_mol_neg_tag && in_ref_neg_tag) {
                    in_both_neg += 1;
                    if (output_cloud == true){
                        charge_cloud.x.push_back(xQ);
                        charge_cloud.y.push_back(yQ);
                        charge_cloud.z.push_back(zQ);
                        charge_cloud.type.push_back("O");
                    }
                }
                zQ += zStep;
            }
            yQ += yStep;
        }
        xQ += xStep;                
    }
    
    total = 0;
    number = 0;


    if (in_mol + in_ref - in_both != 0) {
        total += in_both / (in_mol + in_ref - in_both);
        total_component = in_both / (in_mol + in_ref - in_both);
        number++;
    }
    else    total_component = 0;
    if (in_mol_hvy + in_ref_hvy - in_both_hvy != 0) {
        total += in_both_hvy / (in_mol_hvy + in_ref_hvy - in_both_hvy);
        heavy_atom_component = in_both_hvy / (in_mol_hvy + in_ref_hvy - in_both_hvy);
        number++;
    }
    else    heavy_atom_component = 0;
    if (in_mol_pos + in_ref_pos - in_both_pos != 0) {
        total += in_both_pos / (in_mol_pos + in_ref_pos - in_both_pos);
        positive_component = in_both_pos / (in_mol_pos + in_ref_pos - in_both_pos);
        number++;
    }
    else    positive_component = 0;
    if (in_mol_neg + in_ref_neg - in_both_neg != 0) {
        total += in_both_neg / (in_mol_neg + in_ref_neg - in_both_neg);
        negative_component = in_both_neg / (in_mol_neg + in_ref_neg - in_both_neg);
        number++;
    }
    else    negative_component = 0;
    if (in_mol_pho + in_ref_pho - in_both_pho != 0) {
        total += in_both_pho / (in_mol_pho + in_ref_pho - in_both_pho);
        hydrophobic_component = in_both_pho / (in_mol_pho + in_ref_pho - in_both_pho);
        number++;
    }
    else    hydrophobic_component = 0;
    if (in_mol_phi + in_ref_phi - in_both_phi != 0) {
        total += in_both_phi / (in_mol_phi + in_ref_phi - in_both_phi);
        hydrophilic_component = in_both_phi / (in_mol_phi + in_ref_phi - in_both_phi);
        number++;
    }
    else    hydrophilic_component = 0;

 
    mol.current_score = total / number;
    mol.current_data = output_score_summary(mol);

    if (output_cloud == true)    write_cloud(cloud);
    if (output_cloud == true)    write_cloud(charge_cloud);
    if (output_cloud == true)    write_cloud(hydrophobicity_cloud);


    return true;
}

bool
Volume_Score::Analytical_method(DOCKMol & mol,float & overlap_ref_ref,float & overlap_ref_ref_hvy,float & overlap_ref_ref_pho,float & overlap_ref_ref_phi,float & overlap_ref_ref_neg,float & overlap_ref_ref_pos)
{
    int         i,
                j;
    float       pi = 3.14,
                d_sqr,
                overlap_temp,
                overlap_mol_ref = 0,
                overlap_mol_ref_hvy = 0,
                overlap_mol_ref_pho = 0,
                overlap_mol_ref_phi = 0,
                overlap_mol_ref_neg = 0,
                overlap_mol_ref_pos = 0,
                overlap_mol_mol = 0,
                overlap_mol_mol_hvy = 0,
                overlap_mol_mol_pho = 0,
                overlap_mol_mol_phi = 0,
                overlap_mol_mol_neg = 0,
                overlap_mol_mol_pos = 0;


    for (i = 0; i < mol.num_atoms; i++){
        if (mol.atom_active_flags[i]) {
            for (j = 0; j < volume_ref_mol.num_atoms; j++){
                d_sqr = (mol.x[i] - volume_ref_mol.x[j]) * (mol.x[i] - volume_ref_mol.x[j]) + (mol.y[i] - volume_ref_mol.y[j]) * (mol.y[i] - volume_ref_mol.y[j]) + (mol.z[i] - volume_ref_mol.z[j]) * (mol.z[i] - volume_ref_mol.z[j]);
                if (d_sqr < (mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j]) * (mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j])) {
                    if (d_sqr != 0)
                        overlap_temp = pi / 12 * (mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j] - sqrt(d_sqr)) * (mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j] - sqrt(d_sqr)) * (sqrt(d_sqr) + 2 * (mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j]) - 3 / sqrt(d_sqr) * (mol.amber_at_radius[i] - volume_ref_mol.amber_at_radius[j]) * (mol.amber_at_radius[i] - volume_ref_mol.amber_at_radius[j]));
                    else
                        overlap_temp = pi / 12 * 2 * (mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j]) * (mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j]) * (mol.amber_at_radius[i] + volume_ref_mol.amber_at_radius[j]);
                    overlap_mol_ref = overlap_mol_ref + overlap_temp;
                    if (mol.atom_types[i] != "H" && volume_ref_mol.atom_types[j] != "H") 
                        overlap_mol_ref_hvy = overlap_mol_ref_hvy + overlap_temp;
                    if (mol.chem_types[i] == "hydrophobic" && volume_ref_mol.chem_types[j] == "hydrophobic")
                        overlap_mol_ref_pho = overlap_mol_ref_pho + overlap_temp;
                    if ((mol.chem_types[i] == "donar" || mol.chem_types[i] == "acceptor" || mol.chem_types[i] == "polar") && (volume_ref_mol.chem_types[j] == "donar" || volume_ref_mol.chem_types[j] == "acceptor" || volume_ref_mol.chem_types[j] == "polar"))
                        overlap_mol_ref_phi = overlap_mol_ref_phi + overlap_temp;
                    if (mol.charges[i] < 0 && volume_ref_mol.charges[j] < 0)
                        overlap_mol_ref_neg = overlap_mol_ref_neg + overlap_temp;
                    if (mol.charges[i] > 0 && volume_ref_mol.charges[j] > 0)
                        overlap_mol_ref_pos = overlap_mol_ref_pos + overlap_temp;
    
                }
            }
        }
    }

    for (i = 0; i < mol.num_atoms; i++){
        if (mol.atom_active_flags[i]) {
            for (j = 0; j < mol.num_atoms; j++){
                if (mol.atom_active_flags[j]) {
                    d_sqr = (mol.x[i] - mol.x[j]) * (mol.x[i] - mol.x[j]) + (mol.y[i] - mol.y[j]) * (mol.y[i] - mol.y[j]) + (mol.z[i] - mol.z[j]) * (mol.z[i] - mol.z[j]);
                    if (d_sqr < (mol.amber_at_radius[i] + mol.amber_at_radius[j]) * (mol.amber_at_radius[i] + mol.amber_at_radius[j])) {
                        if (d_sqr != 0)
                            overlap_temp = pi / 12 * (mol.amber_at_radius[i] + mol.amber_at_radius[j] - sqrt(d_sqr)) * (mol.amber_at_radius[i] + mol.amber_at_radius[j] - sqrt(d_sqr)) * (sqrt(d_sqr) + 2 * (mol.amber_at_radius[i] + mol.amber_at_radius[j]) - 3 / sqrt(d_sqr) * (mol.amber_at_radius[i] - mol.amber_at_radius[j]) * (mol.amber_at_radius[i] - mol.amber_at_radius[j]));
                        else
                            overlap_temp = pi / 12 * 2 * (mol.amber_at_radius[i] + mol.amber_at_radius[j]) * (mol.amber_at_radius[i] + mol.amber_at_radius[j]) * (mol.amber_at_radius[i] + mol.amber_at_radius[j]);
                        overlap_mol_mol = overlap_mol_mol + overlap_temp;
                        if (mol.atom_types[i] != "H" && mol.atom_types[j] != "H")
                            overlap_mol_mol_hvy = overlap_mol_mol_hvy + overlap_temp;
                        if (mol.chem_types[i] == "hydrophobic" && mol.chem_types[j] == "hydrophobic")
                            overlap_mol_mol_pho = overlap_mol_mol_pho + overlap_temp;
                        if ((mol.chem_types[i] == "donar" || mol.chem_types[i] == "acceptor" || mol.chem_types[i] == "polar") && (mol.chem_types[j] == "donar" || mol.chem_types[j] == "acceptor" || mol.chem_types[j] == "polar"))
                            overlap_mol_mol_phi = overlap_mol_mol_phi + overlap_temp;
                        if (mol.charges[i] < 0 && mol.charges[j] < 0)
                            overlap_mol_mol_neg = overlap_mol_mol_neg + overlap_temp;
                        if (mol.charges[i] > 0 && mol.charges[j] > 0)
                            overlap_mol_mol_pos = overlap_mol_mol_pos + overlap_temp;

                    }
                }
            }
        }
    }


    //cout << overlap_mol_mol << endl;

    if (overlap_mol_mol == 0 && overlap_ref_ref == 0)
        total_component = 0;
    else
        total_component = overlap_mol_ref / max(overlap_mol_mol, overlap_ref_ref);
    if (overlap_mol_mol_hvy == 0 && overlap_ref_ref_hvy == 0)
        heavy_atom_component = 0;
    else
        heavy_atom_component = overlap_mol_ref_hvy / max(overlap_mol_mol_hvy, overlap_ref_ref_hvy);
    if (overlap_mol_mol_pho == 0 && overlap_ref_ref_pho == 0)
        hydrophobic_component = 0;
    else
        hydrophobic_component = overlap_mol_ref_pho / max(overlap_mol_mol_pho, overlap_ref_ref_pho);
    if (overlap_mol_mol_phi == 0 && overlap_ref_ref_phi == 0)
        hydrophilic_component = 0;
    else
        hydrophilic_component = overlap_mol_ref_phi / max(overlap_mol_mol_phi, overlap_ref_ref_phi);
    if (overlap_mol_mol_neg == 0 && overlap_ref_ref_neg == 0)
        negative_component = 0;
    else
        negative_component = overlap_mol_ref_neg / max(overlap_mol_mol_neg, overlap_ref_ref_neg);
    if (overlap_mol_mol_pos == 0 && overlap_ref_ref_pos == 0)
        positive_component = 0;
    else
        positive_component = overlap_mol_ref_pos / max(overlap_mol_mol_pos, overlap_ref_ref_pos);


    mol.current_score = (total_component + heavy_atom_component + hydrophobic_component + hydrophilic_component + negative_component + positive_component) / 6;
    mol.current_data = output_score_summary(mol);

    return true;

}


bool
Volume_Score::compute_score(DOCKMol & mol)
{

    
    if (volume_compare_type == "grid")
        Grid_method(mol);

    if (volume_compare_type == "analytical")
        Analytical_method(mol,overlap_ref_ref,overlap_ref_ref_hvy,overlap_ref_ref_pho,overlap_ref_ref_phi,overlap_ref_ref_neg,overlap_ref_ref_pos);

    //cout << mol.current_score << endl;




/*
    if (volume method == meth1){
        function here to compute score, number


          mol.current_score = score;
     }
     elsif ( == methd2){
          function to compute score in a different way

          mol.current_score = score;
     }

     return true;

*/
//end here

/*
    int atom;

    if (use_score) {

        // check to see if molecule is inside grid box
        for (atom = 0; atom < mol.num_atoms; atom++) {
            if (!energy_grid->is_inside_grid_box(mol.x[atom], mol.y[atom], mol.z[atom])) {
                mol.current_score = -MIN_FLOAT;  // arbitrarily large score
                mol.current_data = "ERROR:  Conformation could not be scored."
                    "\nConformation not completely within grid box.\n";
                return false;
            }
        }

        float es_val = 0.0;
        float vdw_val = 0.0;
        float total = 0.0;
        for (atom = 0; atom < mol.num_atoms; atom++) {

            if (mol.atom_active_flags[atom]) {

                energy_grid->find_grid_neighbors(mol.x[atom], mol.y[atom], mol.z[atom]);

                vdw_val +=
                    ((vdwA[mol.amber_at_id[atom]] * energy_grid->interpolate(energy_grid->avdw)) -
                     (vdwB[mol.amber_at_id[atom]] * energy_grid->interpolate(energy_grid->bvdw))) *
                    vdw_scale;
                es_val += mol.charges[atom] * energy_grid->interpolate(energy_grid->es) * es_scale;
            }
        }

        total = vdw_val + es_val;

        vdw_component = vdw_val;
        es_component = es_val;

        mol.current_score = total;
        mol.current_data = output_score_summary(mol);

    }
*/

    return true;

} //end Volume_Score::compute_score()


// +++++++++++++++++++++++++++++++++++++++++
void
Volume_Score::write_cloud(Cloud_Struct & cloud)
{
//write out a mol2 file of the cloud of points of volume
    ostringstream file;
    if (cloud.type[0] == "N") file << "cloud.mol2";
    if (cloud.type[0] == "H" || cloud.type[0] == "O") file << "charge_cloud.mol2";
    if (cloud.type[0] == "C" || cloud.type[0] == "S") file << "hydrophobicity_cloud.mol2";
    fstream fout;
    fout.open (file.str().c_str(), fstream::out);

    int num_points = cloud.x.size();
    char line[200];

    fout << "@<TRIPOS>MOLECULE" << endl;
    //file << "cloud.mol2";
    if (cloud.type[0] == "C") fout << "cloud.mol2" << endl;
    if (cloud.type[0] == "H" || cloud.type[0] == "O") fout << "charge_cloud.mol2" << endl;
    if (cloud.type[0] == "C" || cloud.type[0] == "S") fout << "hydrophobicity_cloud.mol2" << endl;

    sprintf(line, " %d 0 1 0 0", num_points);
    fout << line << endl;

    fout << "Volume" << endl;
    fout << "NO_CHARGES"   << endl << endl << endl;

    fout << "@<TRIPOS>ATOM" << endl;

    for (int i=0; i < num_points; i++){
        if (cloud.type[i] == "N") sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", i+1, " N", cloud.x[i], cloud.y[i], cloud.z[i], " N", "1", "VOL", "0.0000");
        if (cloud.type[i] == "H") sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", i+1, " H", cloud.x[i], cloud.y[i], cloud.z[i], " H", "1", "VOL", "0.0000");
        if (cloud.type[i] == "O") sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", i+1, " O", cloud.x[i], cloud.y[i], cloud.z[i], " O", "1", "VOL", "0.0000");
        if (cloud.type[i] == "C") sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", i+1, " C", cloud.x[i], cloud.y[i], cloud.z[i], " C", "1", "VOL", "0.0000");
        if (cloud.type[i] == "S") sprintf(line, "%7d%-7s%12.4f%10.4f%10.4f%2s%8s%5s%15s", i+1, " S", cloud.x[i], cloud.y[i], cloud.z[i], " S", "1", "VOL", "0.0000");
        fout << line << endl;
    }
    
    fout << "@<TRIPOS>BOND" << endl;
    fout << "@<TRIPOS>SUBSTRUCTURE" << endl;
    fout << " 1   ****   1" << endl;
    fout << endl;
    fout.close();

    return;
}




// +++++++++++++++++++++++++++++++++++++++++
string
Volume_Score::output_score_summary(DOCKMol & mol)
{
    ostringstream text;

    //changed output score component names to manage length and be consistent.  Yuchen 10/24/2016
    text << DELIMITER << setw(STRING_WIDTH) << "Property_Volume_Score:"
         << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl;
    text << DELIMITER << setw(STRING_WIDTH) << "VOS_geometric_vo:"
         << setw(FLOAT_WIDTH) << fixed << total_component << endl;
    text << DELIMITER << setw(STRING_WIDTH) << "VOS_hvy_atom_vo:"
         << setw(FLOAT_WIDTH) << fixed << heavy_atom_component << endl;
    text << DELIMITER << setw(STRING_WIDTH) << "VOS_pos_chrg_atm_vo:"
         << setw(FLOAT_WIDTH) << fixed << positive_component << endl;
    text << DELIMITER << setw(STRING_WIDTH) << "VOS_neg_chrg_atm_vo:"
         << setw(FLOAT_WIDTH) << fixed << negative_component << endl;
    text << DELIMITER << setw(STRING_WIDTH) << "VOS_hydrophobic_atm_vo:"
         << setw(FLOAT_WIDTH) << fixed << hydrophobic_component << endl;
    text << DELIMITER << setw(STRING_WIDTH) << "VOS_hydrophilic_atm_vo:"
         << setw(FLOAT_WIDTH) << fixed << hydrophilic_component << endl;
    if (use_internal_energy){
         text << DELIMITER << setw(STRING_WIDTH) << "Internal_energy_repulsive:"
              << setw(FLOAT_WIDTH) << fixed << mol.internal_energy << endl;
         }


/*
    if (use_score) {

        text << DELIMITER << setw(STRING_WIDTH) << "Grid Score:"
             << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl;
        text << DELIMITER << setw(STRING_WIDTH) << "Grid_vdw:" 
             << setw(FLOAT_WIDTH) << fixed << vdw_component << endl
             << DELIMITER << setw(STRING_WIDTH) << "Grid_es:"  
             << setw(FLOAT_WIDTH) << fixed << es_component  << endl;

        // Compute lig internal energy with segments (sudipto 12-12-08)
        //float int_vdw_att, int_vdw_rep, int_es;
        //compute_ligand_internal_energy(mol, int_vdw_att, int_vdw_rep, int_es);

        if (use_internal_energy) 
            text << DELIMITER << setw(STRING_WIDTH) << "Int_energy:"
                 << setw(FLOAT_WIDTH) << fixed << mol.internal_energy << endl;

        // To get the correct internal energy, it MUST be calculated before calling the primary score

    }
*/

    return text.str();

} //end Volume_Score::output_score_summary()



