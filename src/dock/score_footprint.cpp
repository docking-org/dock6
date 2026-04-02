 #include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <sstream>
#include "amber_typer.h"
#include "base_score.h"
#include "dockmol.h"
#include "grid.h" 
#include "score_footprint.h"
#include "utils.h"
#include <math.h>

using namespace std;


// +++++++++++++++++++++++++++++++++++++++++
RANGE::RANGE()
{
    fps_range_vdw = "";
    fps_range_es = "";
    fps_range_hb ="";
}


// +++++++++++++++++++++++++++++++++++++++++
// This function is used to parse the footprint txt file, when applicable
bool
read_footprint_txt(DOCKMol & ref,istream & ifs)
{

    FOOTPRINT_ELEMENT temp_foot;
    string line;

    // get every string or number from file. 
    while (getline(ifs,line,'\n')) {
        vector < string > tokens;
       
        if (!line.empty()){
            Tokenizer(line, tokens,' ');
            if (tokens.size() > 0){
                // ignore the column names and comments lines
                if (!(tokens[0].compare("resname") == 0) && !(tokens[0].compare(0,1,"#") == 0)){
                    // terminate if wrong number of columns in a line
                    if (tokens.size() != 5){
                        cout << "A line in the footprint reference txt file has an incorrect" << endl 
                             << "number of entries. Lines must have 5 entries separated by spaces" << endl
                             << "as follows:" << endl
                             << "\t\tresname  resid  vdw_fp  es_fp  hb_fp" << endl
                             << "Note that lines beginning with a comment \"#\" will be ignored" << endl
                             << "Program will terminate." << endl;
                        cout << "line which fails on: " << line << endl;
                        exit(1);

                    }
                    stringstream ss; ss.clear();
                    temp_foot.resname = tokens[0];
                    ss <<tokens[1]; ss >> temp_foot.resid; ss.clear();
                    ss <<tokens[2]; ss >> temp_foot.vdw; ss.clear();
                    ss <<tokens[3]; ss >> temp_foot.es; ss.clear();
                    ss <<tokens[4]; ss >> temp_foot.hb; ss.clear();
                    // cout << temp_foot.resname <<" "<< temp_foot.resid<<" "
                    //     << temp_foot.vdw <<" "<< temp_foot.es <<" "
                    //     << temp_foot.hb << endl;
                    ref.footprints.push_back(temp_foot); 
                }
            }
        }
        tokens.clear();
        // else cout << "line is empty" << line << endl;
    }

    if (ref.footprints.size() > 0)
        return true;

    return false;

} // end Footprint_Similarity_Score::read_footprint_txt()


// +++++++++++++++++++++++++++++++++++++++++
// Constructor and desctructor - nothing is initialized at the moment
Footprint_Similarity_Score::Footprint_Similarity_Score()
{
}

Footprint_Similarity_Score::~Footprint_Similarity_Score()
{
}


// +++++++++++++++++++++++++++++++++++++++++
// When this class is closed, delete some of the C++ strings. Is there a reason
// we don't put this in the destructor?
void
Footprint_Similarity_Score::close()
{
    // clear contents of c++ strings 
    // fps_range_vdw.clear();
    // fps_range_es.clear();
    // fps_range_hb.clear();
    receptor_filename.clear();
    fps_foot_compare_type.clear();
    fps_fp_info.clear();
    constant_footprint_ref_file.clear();

} // end Footprint_Similarity_Score::close()


// +++++++++++++++++++++++++++++++++++++++++
// Input parameters for footprint similarity scoring function - called in
// master score
void
Footprint_Similarity_Score::input_parameters(Parameter_Reader & parm,
                                          bool & primary_score,
                                          bool & secondary_score)
{
    string          tmp;

    use_primary_score = false;
    use_secondary_score = false;

    fps_foot_specify_a_range = false;
    fps_foot_specify_a_threshold = false;
    fps_normalize_foot = false;
    // remainder is an element of the footprint vector of the remaining residues
    fps_use_remainder = false;

    // initialize parameters
    fps_vdw_num_resid = 0;
    fps_es_num_resid = 0;
    fps_hb_num_resid = 0;
    fps_fp_info = "";

    cout << "\nFootprint Similarity Score Parameters\n"
         << "--------------------------------------------------------"
            "----------------------------------"
         << endl;

    if (!primary_score) {
        tmp = parm.query_param("footprint_similarity_score_primary", "no", "yes no");
        use_primary_score = (tmp == "yes");
        primary_score = use_primary_score;
    }

    if (!secondary_score) {
        //tmp = parm.query_param("footprint_similarity_score_secondary", "no", "yes no");
        tmp = "no";
	use_secondary_score = (tmp == "yes");
        secondary_score = use_secondary_score;
    }

    use_score = (use_primary_score || use_secondary_score);

    if (use_score) {
        input_parameters_main(parm, "fps_score");
    }
    return;
} // end input_parameters


// +++++++++++++++++++++++++++++++++++++++++
// Input parameters for footprint similarity scoring function - called above and in
// description score
void
Footprint_Similarity_Score::input_parameters_main(Parameter_Reader & parm,
                                          string parm_head)
{
    string          tmp;

    // if footprints are precomputed, read in txt file with footprints
    // else ask for the reference mol2 file
    tmp = parm.query_param(parm_head+"_use_footprint_reference_mol2", "no", "yes no");

    bool_footref_mol2 = (tmp == "yes");
    bool_footref_txt = false;

    if (bool_footref_mol2)
        constant_footprint_ref_file = parm.query_param(
            parm_head+"_footprint_reference_mol2_filename", "ligand_footprint.mol2");
    else {
        tmp = parm.query_param(parm_head+"_use_footprint_reference_txt", "no", "yes no");

        bool_footref_txt = (tmp == "yes");

        if (bool_footref_txt)
            constant_footprint_ref_file = parm.query_param(
                parm_head+"_footprint_reference_txt_filename", "ligand_footprint.txt");
        else {
            cout << "Error: No footprint reference filename is specified."
                    "Program will terminate."
                 << endl;
            exit(0);
        }
    }

    // ask what metric to use to compare the footprints
    fps_foot_compare_type = 
        parm.query_param(parm_head+"_foot_compare_type", "Euclidean", "Pearson, Euclidean");
 
    if (fps_foot_compare_type.compare("Pearson") == 0)
        cout <<"-------------------------------------------------------------"<<endl
             << "You chose to use the correlation coefficient as the metric " << endl
             << "to compare the footprints.  When the value is 1 then " << endl
             << "there is perfect agreement between the two footprints." << endl
             << "When the value is 0 then there is poor agreement" << endl
             << "between the two footprints." << endl
             <<"-------------------------------------------------------------"<<endl;
    else if (fps_foot_compare_type.compare("Euclidean") == 0){
        cout <<"-------------------------------------------------------------"<<endl
             << "You chose to use the Euclidean distance as the metric " << endl
             << "to compare the footprints.  When the value is 0 then " << endl
             << "there is perfect agreement between the two footprints." << endl
             << "As the agreement gets worse between the two" << endl
             << "footprints the value increases." << endl
             <<"-------------------------------------------------------------"<<endl;
        tmp = parm.query_param(parm_head+"_normalize_foot", "no", "yes no");
        if (tmp == "yes")
            fps_normalize_foot = true;
    }

    // use all residues; yes or no
    tmp = parm.query_param(parm_head+"_foot_comp_all_residue", "yes", "yes no");

    if (tmp == "yes")
        fps_foot_comp_all_residue = true;
    if (tmp == "no")
        fps_foot_comp_all_residue = false;

    if (!fps_foot_comp_all_residue){
        tmp =
            parm.query_param(parm_head+"_choose_foot_range_type", "specify_range", "specify_range threshold");
        if (tmp == "specify_range"){
            fps_foot_specify_a_range = true;
            cout <<"-------------------------------------------------------------"<<endl
                 << "You chose to use a residue range that you specify.  All footprints" << endl
                 << "will be evaluated (compared to ref.) only on this residue range." << endl
                 << "First residues id = 1, not zero." << endl
                 <<"-------------------------------------------------------------"<<endl;
        }
        if (tmp == "threshold"){
            fps_foot_specify_a_threshold = true;
            cout <<"-------------------------------------------------------------"<<endl
                 << "You chose to use a residue range that is defined by only residues that" << endl
                 << "have magnitudes that excede the specified thresholds (the 3 following" << endl
                 << "parameters).  Each footprint type will have a distinct residue range."<< endl
                 << "The union of both the ref. and pose ranges will define the"<<endl
                 << "evaluated range." << endl
                 <<"-------------------------------------------------------------"<<endl;
        }
        if (fps_foot_specify_a_range) {
            string fps_range = parm.query_param(parm_head+"_range", "20-40,60-80");
            ref_ranges.fps_range_vdw = fps_range;
            ref_ranges.fps_range_es  = fps_range;
            ref_ranges.fps_range_hb  = fps_range;
            pose_ranges.fps_range_vdw = fps_range;
            pose_ranges.fps_range_es  = fps_range;
            pose_ranges.fps_range_hb  = fps_range;
        }
        if (fps_foot_specify_a_threshold) {
            fps_vdw_threshold = atof(parm.query_param(parm_head+"_vdw_threshold", "1").c_str());
            fps_es_threshold  = atof(parm.query_param(parm_head+"_es_threshold", "1").c_str());
            fps_hb_threshold  = atof(parm.query_param(parm_head+"_hb_threshold", "0.5").c_str());
        }

    // if yes no nothing.
    // if no ask for range to evaluate footprint comparison. 
    // examples:
    //     1-50,57-63
    //     1,3,4,6,10
    //    

        tmp = parm.query_param(parm_head+"_use_remainder", "yes", "yes no");
        if (tmp == "yes")
            fps_use_remainder = true;
    }

    receptor_filename =
        parm.query_param(parm_head+"_receptor_filename", "receptor.mol2");

    att_exp = atoi(parm.query_param(parm_head+"_vdw_att_exp", "6").c_str());
    if (att_exp <= 0) {
        cout << "ERROR: Parameter must be an integer greater than zero. Program will terminate."
             << endl;
        exit(0);
    }
    rep_exp = atoi(parm.query_param(parm_head+"_vdw_rep_exp", "12").c_str());
    if (rep_exp <= 0) {
        cout << "ERROR: Parameter must be an integer greater than zero. Program will terminate."
             << endl;
        exit(0);
    }
    rep_radius_scale = atof(parm.query_param(parm_head+"_vdw_rep_rad_scale", "1").c_str());
    if (rep_radius_scale <= 0.0) {
        cout << "ERROR: Parameter must be a float greater than zero. Program will terminate."
             << endl;
        exit(0);
    }

    //Ask whether we should use a distance dependent dieleectric

    tmp = parm.query_param(parm_head+"_use_distance_dependent_dielectric", "yes", "yes no");
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
        atof(parm.query_param(parm_head+"_dielectric", default_diel).c_str());
    if (diel_screen <= 0.0) {
        cout << "ERROR: Parameter must be a float greater than zero. Program will terminate."
             << endl;
        exit(0);
    }

    /* These will now only appear in descriptor score
    vdw_scale = atoi(parm.query_param(parm_head+"_score_vdw_scale", "1").c_str());
    if (vdw_scale <= 0) {
        cout << "WARNING: Parameter should be an integer greater than zero."<< endl;
    }
    es_scale = atoi(parm.query_param(parm_head+"_score_es_scale", "1").c_str());
    if (es_scale <= 0) {
        cout <<"WARNING: Parameter should be an integer greater than zero." << endl;
    }
    hb_scale = atoi(parm.query_param(parm_head+"_score_hb_scale", "0").c_str());
    if (hb_scale>0){
        cout <<"WARNING: Parameter should be a number less than zero."<< endl;
    }
    internal_scale = atoi(parm.query_param(parm_head+"_score_internal_scale", "0").c_str());
    */

    fp_vdw_scale = atoi(parm.query_param(parm_head+"_vdw_fp_scale", "1").c_str());
    //if (fp_vdw_scale<0){
    //    cout <<"WARNING: Parameter should be a number less than zero."<< endl;
    //}

    fp_es_scale = atoi(parm.query_param(parm_head+"_es_fp_scale", "1").c_str());
    //if (fp_es_scale<0){
    //    cout <<"WARNING: Parameter should be a number less than zero."<< endl;
    //}

    fp_hb_scale = atoi(parm.query_param(parm_head+"_hb_fp_scale", "0").c_str());
    //if (fp_hb_scale<0){
    //    cout <<"WARNING: Parameter should be a number less than zero."<< endl;
    //}

    /*
    rot_bonds_scale  = atoi(parm.query_param(parm_head+"_score_rot_bonds_scale", "0").c_str());
    heavyatoms_scale = atoi(parm.query_param(parm_head+"_score_Heavy_Atoms_scale", "0").c_str());
    hb_don_scale = atoi(parm.query_param(parm_head+"_score_HB_Donors_scale", "0").c_str());
    hb_acc_scale = atoi(parm.query_param(parm_head+"_score_HB_Acceptors_scale", "0").c_str());
    molecular_wt_scale = atoi(parm.query_param(parm_head+"_score_Molecular_Wt_scale", "0").c_str());
    formal_charge_scale = atoi(parm.query_param(parm_head+"_score_Formal_Charge_scale", "0").c_str());
    XlogP_scale = atoi(parm.query_param(parm_head+"_score_XLogP_scale", "0").c_str());
    */
    return;

} // end Footprint_Similarity_Score::input_parameters()


// +++++++++++++++++++++++++++++++++++++++++
// Initialize some values for the Amber typer
void
Footprint_Similarity_Score::initialize(AMBER_TYPER & typer)
{   
    //Footprint citation
    cout << "To cite Footprint Similarity Score use: \n Balius, T. E.; Mukherjee, S.; Rizzo, R. C. Implementation and Evaluation of a Docking-rescoring Method using Molecular Footprint Comparisons. J. Comput. Chem., 2011, 32, 2273-2289.\n" << endl;

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
        use_ph4  = false;
        use_volume = false;
        typer.prepare_molecule(receptor, read_vdw, use_chem, use_ph4, use_volume);
    }

    return;

} // end Footprint_Similarity_Score::initialize()


// +++++++++++++++++++++++++++++++++++++++++
// Get footprint reference from disk
void
Footprint_Similarity_Score::submit_footprint_reference(AMBER_TYPER & typer)
{ 
    ifstream        footprint_ref;
    footprint_ref.open(constant_footprint_ref_file.c_str());

    if (bool_footref_mol2){
        bool temp_mol2_OK = Read_Mol2(footprint_reference, footprint_ref, false, false, false);
        if(!(temp_mol2_OK)){
            cout << "there is a problem reading the footprint_reference mol2 file" <<endl;
            exit(0);
        }
        footprint_ref.close();
        bool read_vdw = true;
        bool use_chem = false;
        bool use_ph4  = false;
        bool use_volume = false;
        typer.prepare_molecule(footprint_reference, read_vdw, use_chem, use_ph4, use_volume);
        bool temp = compute_footprint(footprint_reference);
    }

    if (bool_footref_txt){
        bool temp_txt_OK = read_footprint_txt(footprint_reference,footprint_ref);
        if (!(temp_txt_OK)){
            cout << "there is a problem reading the footprint_reference txt file" <<endl;
            exit(0);
        }
    }

    // if tresholds are specified then calculate ranges for ref.
    if (fps_foot_specify_a_threshold){
        ref_ranges = range_satisfying_threshold_all(footprint_reference.footprints);
    }

    return;

} // end Footprint_Similarity_Score::submit_footprint_reference()

// Do the floating point variables in this file really need to be double ?
// This is inconsistent with most of dock which is float
// and has caused minor compilation issues.
// srb Aug 2013.


// +++++++++++++++++++++++++++++++++++++++++
// Compute h-bonds
bool
hbond_cal(DOCKMol & mol_a,DOCKMol & mol_d,int i,int j,double & hbond_dist, double & hbond_angle)
{
    // distances between Hd --- A.  non-covalent
    double threshold = 2.5;
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
        hbond_dist = sqrt((x_diff*x_diff)+(y_diff*y_diff) + (z_diff*z_diff));
        // eliminate any thing outside the spheir with radius of 3.
        if (hbond_dist <= threshold) 
        {
            // distances beteewn Ad-Hd. covalent bond.
            double x_diff_covalent = mol_d.x[j]-mol_d.x[mol_d.acc_heavy_atomid[j]];
            double y_diff_covalent = mol_d.y[j]-mol_d.y[mol_d.acc_heavy_atomid[j]];
            double z_diff_covalent = mol_d.z[j]-mol_d.z[mol_d.acc_heavy_atomid[j]];
            double covalent_dist   = sqrt((x_diff_covalent*x_diff_covalent)+(y_diff_covalent*y_diff_covalent) + (z_diff_covalent*z_diff_covalent));
            //calculate the dot product
            double dot_coval_hbond =  x_diff*x_diff_covalent+y_diff*y_diff_covalent+z_diff*z_diff_covalent;
            hbond_angle = acos(dot_coval_hbond / (covalent_dist*hbond_dist))*180/PI; 
            if (hbond_angle <= 180 && hbond_angle >= 120){
                return true;
            }
        }
    }

    return false;

} // end hbond_cal()


// +++++++++++++++++++++++++++++++++++++++++
//    calculate footprint of mol with compute_footprint
//    if reference footprint is calculated then
//        calculate distance between reference and mol footprints
//    else 
//        calculate reference footprint with compute_footprint 
//        and compute distance
//    output score.
bool
Footprint_Similarity_Score::compute_footprint(DOCKMol & mol)
{
    double          vdw_val,
                    es_val,
                    total;
    int             i,
                    j;
    double          dist;
    double          single_vdw_val,
                    single_es_val,
                    single_total_val;
    stringstream    hbtxt; 
    //these three variables are used to store energies for atom pairs;

    total = vdw_val = es_val = 0.0;
    // footprints must be empty
    mol.footprints.clear();

    //mol.num_residues_in_rec = receptor.num_residues;
    int num_residues = receptor.num_residues;
    double current_resnum_vdw = 0;
    double current_resnum_es  = 0;
    int    current_resnum_hb  = 0;
    int    current_resnum     = 0;// this is a counter
    string  resname_old = receptor.subst_names[0];
    int  resid_old = atoi(receptor.atom_residue_numbers[0].c_str());

    if (use_score == 1) {
        int hbond_total = 0; 
        int prev_resid  = 0; 
        for (j = 0; j < receptor.num_atoms; j++) {
            // This is a check that the residue atoms are all together.
            // Some programs ie chimera will place hydrogens at the end of the
            // atom list rather than with the proper residue.
            // This was an issue with the current implamentation.
            // Could sort by resnum internally?
            if ( atoi(receptor.atom_residue_numbers[j].c_str()) < prev_resid) {
                cout << "ERROR: Can not process footprints, resids in receptor mol2 are not monotonic:"
                     << receptor.atom_residue_numbers[j] << "<" << prev_resid
                     << endl;
                exit(0);
                // We may wish to change the way we are generating the footprints.
                // for example we could initialize an array with all 0.0 elements 
                // then sum the energy to that residue hashing (accessing the
                // element) based on resid
            }
            single_vdw_val           = 0.0;
            single_es_val            = 0.0;
            int hbond_count              = 0; 
            int hbond_acc_from_pro_count = 0; 
            int hbond_don_to_pro_count   = 0; 
 
            // in this loop, calculate energy for each atom pair
            for (i = 0; i < mol.num_atoms; i++) {

                dist = sqrt(((mol.x[i] - receptor.x[j])*(mol.x[i] - receptor.x[j])) +
                            ((mol.y[i] - receptor.y[j])*(mol.y[i] - receptor.y[j])) + 
                            ((mol.z[i] - receptor.z[j])*(mol.z[i] - receptor.z[j])));
                //added by jwu for test
                double temp;
                //jwu_temp is just a temp variable to store the increase
                temp = (((vdwA[mol.amber_at_id[i]] *
                          vdwA[receptor.amber_at_id[j]]) / pow(dist,
                                                               double(rep_exp))) -
                        ((vdwB[mol.amber_at_id[i]] *
                          vdwB[receptor.amber_at_id[j]]) / pow(dist,
                                                               double(att_exp))));
 
                single_vdw_val += temp;
                vdw_val += temp;
 
                // IF LIG atom is donator and REC atom is acceptor then
                // see if h-bond and return angle and distance
                // IF LIG atom is acceptor and REC atom is donator then
                // see if h-bond and return angle and distance
                if (mol.flag_acceptor[i] && receptor.flag_donator[j])
                {
                    double hbond_dist,hbond_angle;
                    if (hbond_cal(mol,receptor,i,j,hbond_dist, hbond_angle)){
                    //cout << "hbond detected: "<< hbond_dist << " " << hbond_angle << endl;
                    hbtxt << "hbond detected: " << receptor.subst_names[j] << " ";
                    hbtxt << "(" << receptor.atom_residue_numbers[j] << ")" << " ";
                    hbtxt << receptor.atom_types[receptor.acc_heavy_atomid[j]]<<" ";
                    hbtxt << receptor.atom_types[j]<<"---"<<mol.subst_names[i]<<" ";
                    hbtxt << "(" << mol.atom_residue_numbers[i] << ")" << " ";
                    hbtxt << mol.atom_types[i] <<" "<< hbond_dist <<" " << hbond_angle  <<endl;
                    hbond_total++;
                    hbond_count++;
                    hbond_acc_from_pro_count++;
                    }
                } 

                if (mol.flag_donator[i] && receptor.flag_acceptor[j])
                {
                    double hbond_dist,hbond_angle;
                    if (hbond_cal(receptor,mol,j,i,hbond_dist, hbond_angle)){
                        hbtxt << "hbond detected: " << mol.subst_names[i]<<" ";
                        hbtxt << "(" << mol.atom_residue_numbers[i] << ")" << " ";
                        hbtxt << mol.atom_types[mol.acc_heavy_atomid[i]]<<" ";
                        hbtxt << mol.atom_types[i]<<"---"<<receptor.subst_names[j]<<" ";
                        hbtxt << "(" << receptor.atom_residue_numbers[j] << ")" << " ";
                        hbtxt << receptor.atom_types[j] << " " << hbond_dist <<" " << hbond_angle  << endl;
                        hbond_total++;
                        hbond_count++;
                        hbond_don_to_pro_count++;
                    }
                } 

                if (use_ddd){
                   temp = ((332 * mol.charges[i] * receptor.charges[j]) /
                          ((dist*dist) * diel_screen));
                     es_val += temp;
                     single_es_val += temp;
                }
                else {
                    temp = ((332 * mol.charges[i] * receptor.charges[j]) /
                           (dist * diel_screen));
                    es_val += temp;
                    single_es_val += temp;
                } 
            }

            // The perresidue footprint 
            if ( atoi(receptor.atom_residue_numbers[j].c_str()) == resid_old)
            {
                current_resnum_vdw = current_resnum_vdw + single_vdw_val;
                current_resnum_es = current_resnum_es + single_es_val;
                current_resnum_hb = current_resnum_hb + hbond_count;
            }
            else 
            {
                FOOTPRINT_ELEMENT temp_footprint_ele;
                temp_footprint_ele.resname = resname_old;
                temp_footprint_ele.resid   = current_resnum + 1;
                temp_footprint_ele.vdw     = current_resnum_vdw;
                temp_footprint_ele.es      = current_resnum_es;
                temp_footprint_ele.hb      = current_resnum_hb;

                mol.footprints.push_back(temp_footprint_ele);
 
                current_resnum_vdw = 0 + single_vdw_val;
                current_resnum_es = 0  + single_es_val;
                current_resnum_hb = 0 + hbond_count;

                resname_old = receptor.subst_names[j];
                resid_old = atoi(receptor.atom_residue_numbers[j].c_str());
                current_resnum++;
            }
            prev_resid = atoi(receptor.atom_residue_numbers[j].c_str());
        }
        //hbtxt <<"#### h-bond total: "<<hbond_total<< endl; 
        total = vdw_val + es_val;

        vdw_component = vdw_val;
        es_component = es_val;
        hbond        = hbond_total;

        FOOTPRINT_ELEMENT temp_footprint_ele;
        temp_footprint_ele.resname = resname_old;
        temp_footprint_ele.resid   = current_resnum + 1;
        temp_footprint_ele.vdw     = current_resnum_vdw;
        temp_footprint_ele.es      = current_resnum_es;
        temp_footprint_ele.hb      = current_resnum_hb;

        mol.footprints.push_back(temp_footprint_ele);
        
    }
    mol.hbond_text_data = hbtxt.str();
    //cout << mol.hbond_text_data;
    if (num_residues != mol.footprints.size()) {
        cout << "Error: inconsistance in number of residues: "<< mol.footprints.size() <<" "<< num_residues <<endl;
        exit(1);
    }

    return true;

} // end Footprint_Similarity_Score::compute_footprint()


// +++++++++++++++++++++++++++++++++++++++++
// This function takes a resid as input and determines if it is within the
// target subset. Example format: 30-47,56,60,88-112
bool 
Footprint_Similarity_Score::on_range(int resid,string fps_range)
{
    vector <string> list,temp;
    int val, rangestart, rangestop;
 
    // fps_foot_comp_all_residue is true then all residues are included in the 
    // footprint comparison calculation. 
    if (fps_foot_comp_all_residue){
        return true;
    } else {
        // continuous chunks of the subsets that are to be included in the comparison
        // are separated by commas
        Tokenizer(fps_range, list,',');
        for (int i = 0;i<list.size();i++){
            //cout << list[i] << endl;
            Tokenizer(list[i], temp,'-');
 
            // if there are no "-" then Tokenizer returns 1 entry in vector temp.
            // this is a single residue chunk.
            if (temp.size() == 1){
                val = atoi(list[i].c_str());
                if (val == resid){
                    //cout << resid << " = " << val << endl;
                    list.clear(); temp.clear();
                    return true;
                }
            }
 
            // this specifies a sub-range between 2 residues. 
            else if (temp.size() == 2){
                rangestart = atoi(temp[0].c_str());
                rangestop = atoi(temp[1].c_str());
                if ((resid >= rangestart) && (resid <= rangestop)){
                    //cout << resid << "is on the range" << rangestart<<"-" << rangestop
                    //     << endl;
                    list.clear(); temp.clear();
                    return true;
                }
            }
 
            else {
                cout << "range specified is not in correct format. program will terminate"
                     << endl;
                exit(1);
            }
        }
    }

    //cout << "return false" << endl;
    list.clear(); temp.clear();
 
    return false;

} // end Footprint_Similarity_Score::on_range()


// +++++++++++++++++++++++++++++++++++++++++
// Range is modified in this function
void
range_satisfying_threshold(string & range, double threshold, vector <double> footprint)
{
    stringstream    range_stream;
    bool first_flage = true;
    for (int i = 0; i < footprint.size(); i++){
        if (first_flage && (fabs(footprint[i]) >= threshold)){
            range_stream << (i+1);
            first_flage = false;
        }
        else if (fabs(footprint[i]) >= threshold){
            range_stream <<"," <<(i+1);
        }
        // else do nothing
    }
    range = range_stream.str();

    return;

} //end range_satisfying_threshold()


// +++++++++++++++++++++++++++++++++++++++++
// Find Range. this function should find a range given a footprint
RANGE
Footprint_Similarity_Score::range_satisfying_threshold_all(vector <FOOTPRINT_ELEMENT> footprints)
{
    RANGE temp;
    temp.fps_range_vdw = "";
    temp.fps_range_es  = "";
    temp.fps_range_hb = "";
    //string range = "";
   
    vector <double> footprint_vdw;
    vector <double> footprint_es;
    vector <double> footprint_hb;

    // separate component (vdw, es, hb) for footprints into 3 vectors
    for (int i = 0; i < footprints.size(); i++){
        footprint_vdw.push_back(footprints[i].vdw);
        footprint_es.push_back(footprints[i].es);
        footprint_hb.push_back(footprints[i].hb); // int footprints[i].hb
    }

    range_satisfying_threshold(temp.fps_range_vdw ,fps_vdw_threshold,footprint_vdw);
    range_satisfying_threshold(temp.fps_range_es  ,fps_es_threshold ,footprint_es);
    range_satisfying_threshold(temp.fps_range_hb  ,fps_hb_threshold ,footprint_hb);
    footprint_vdw.clear(); footprint_es.clear(); footprint_hb.clear();

    return temp;

} // end Footpring_Similarity_Score::range_satisfying_threshold_all()


// +++++++++++++++++++++++++++++++++++++++++
// Find the union of vdw / es / hb ranges. range1 is always the pose, range2
// is always the reference.
string 
Footprint_Similarity_Score::union_of_ranges(string range1,string range2)
{
    if(fps_foot_specify_a_threshold){
        stringstream    range_union;
        // range_union.clear();

        bool first_flage = true;
       
        // loop over footprint residues
        for (int i=0;i<footprint_reference.footprints.size();i++){
            if (on_range((i+1),range1) or on_range((i+1),range2)){
                if (first_flage){
                    range_union << (i+1);
                    first_flage = false;
                }
                else{
                    range_union <<","<< (i+1);
                }
            }
        }
 
        return range_union.str();

    }

    return range1;

} // end Footprint_Similarity_Score::union_of_ranges()


// +++++++++++++++++++++++++++++++++++++++++ 
// Computes the normalization factor (norm)
double 
Norm2(const vector <double> *footprint)
{
    //if (footprint->size() == 0) return 1;
    double sum2 = 0;

    for (int i=0;i<footprint->size();i++){
        sum2 = sum2 + ((footprint->at(i))*(footprint->at(i)));
    }

    if (fabs(sum2) <= 0.0000001){
        cout << "WARNING: sum of squares close to zero in norm calculation."
             << endl;
        return 1;
    }

    return sqrt(sum2);

} // end Norm2()


// +++++++++++++++++++++++++++++++++++++++++ 
// Normalizes a (Euclidean) footprint vector
void
Normalization(vector <double> *footprint)
{
    double  norm = Norm2(footprint);
    for (int i=0;i<footprint->size();i++) {
        //(*footprint)[i] = (*footprint)[i]/norm;
        footprint->at(i) = footprint->at(i)/norm;
    }

    return;

} // end Normalization()


// +++++++++++++++++++++++++++++++++++++++++
// Print the residue range
void
print_range(string range1, string range2, string range_union, string type, string * print_info)
{
    stringstream ss;
    ss << print_info->c_str();
    ss << "#" << " footprint type:"<< type<<endl;   
    ss << "#" << "   range_pose:" <<endl;
    ss << "#" << "      " << range1 <<endl;
    ss << "#" << "   range_ref:" <<endl;
    ss << "#" << "      " << range2 <<endl;
    ss << "#" << "   range_union:" <<endl;
    ss << "#" << "     " << range_union <<endl;
    //print_info->at() = ss.str();
    print_info->clear();
    print_info->append(ss.str());

    return;

} // end print_range()


// +++++++++++++++++++++++++++++++++++++++++
// Print a footprint
void
print_foot(vector <string> *names,vector <double> *fp1,vector <double> *fp2, string * print_info)
{
    stringstream ss;
    ss << print_info->c_str();
    if (fp1->size() != fp2->size() || fp1->size() != names->size()) {
        ss << "fp1->size() != fp2->size()" << endl;
        return;
    }
    ss << "#" << "compare footprints " << endl;
    for (int i = 0; i < fp1->size(); i++){
        ss << "#" << i << "," << names->at(i) <<"," << fp1->at(i) <<"," << fp2->at(i) << endl;
    }
    print_info->clear();
    print_info->append(ss.str());

    return;

} //end print_foot()


// +++++++++++++++++++++++++++++++++++++++++
// Return Euclidean distance between two footprints
double
Footprint_Similarity_Score::Euclidean_distance( vector <double> *footprint_ref,vector <double> *footprint)
{
    double distance2 = 0;
    if (footprint_ref->size() != footprint->size()){
        cout << "Error: reference and pose footprints do not have same size." << endl;
        exit(1);
    }
    if (footprint->size() == 0){
        cout << "Note: footprint range is empty" << endl;
        return 0;
    }

    // the simple Euclidean distances
    distance2 = 0;

    for (int i=0;i<footprint->size();i++){
        double delfp = footprint->at(i) - footprint_ref->at(i);
        distance2=distance2 + delfp*delfp;
    }

    return sqrt(distance2);

} // end Footprint_Similarity_Score::Euclidean_distance


// +++++++++++++++++++++++++++++++++++++++++
// Return Pearson correlation between two footprints
double
Footprint_Similarity_Score::Pearson_correlation( vector <double> *footprint_ref,vector <double> *footprint)
{
    double r,sumx,sumy,sumxy,sumsqx,sumsqy; //x for footprint and y for footprint_ref
                                           //to store types of summations
    r = 0;
    sumx = 0;
    sumy = 0;
    sumxy = 0;
    sumsqx = 0;
    sumsqy = 0;

    if (footprint_ref->size() != footprint->size()){
        cout << "Error: reference and pose footprints are not the same size." << endl;
        exit(1);
    }

    if (footprint->size() == 0){
        cout << "Note: footprint range is empty" << endl;
        return 0;
    }
    else if (footprint->size() < 3){
        cout << "Warning: r-value is not meaningful: footprint range" << endl
             << "consists of only one or two residues." << endl;
    }

    for (int i=0;i<footprint->size();i++){
        sumx = sumx + footprint->at(i);
        sumy = sumy + footprint_ref->at(i);
        sumxy = sumxy + footprint->at(i)*footprint_ref->at(i);
        sumsqx = sumsqx + footprint->at(i)*footprint->at(i);
        sumsqy = sumsqy + footprint_ref->at(i)*footprint_ref->at(i);
    }

    // this value is calculated to avoid 'nan' in scores
    double r_1=sqrt((footprint->size()*sumsqx-sumx*sumx)*(footprint->size()*sumsqy-sumy*sumy));

    if (r_1!=0.0){
        //r = (n*sumxy - sumx*sumy)/ r_1;
        r = (footprint->size()*sumxy - sumx*sumy)/ r_1;
        // the r-value formula is here. note: the r value here can be negative
    }
    else {
        cout << "Note: For one footprint, the r_1=sqrt((N*sumsqx-sumx*sumx)*(N*sumsqy-sumy*sumy)) is 0." << endl
             << "This is most likely because all points are the same point eg (0.0,0.0)." << endl;
        r = 0.0;
    }

    return r;

} // end Footprint_Similarity_Score::Pearson_correlation()


// because of the datastructure footprints this function is not going to work
//void get_footprints( vector <double> * fp_ref_orig, vector <double> * fp_orig,
//                     vector <double> *  fp_ref_new,  vector <double> * fp_new,
//                     string range, string text){
//   fp_new->clear(); footprint_ref_tmp.clear();
//   for (int i=0;i<n;i++)
//   {
//      if (on_range((i+1),range)){
//         fp_new.push_back();
//         footprint_ref_tmp.push_back(footprint_reference.footprints[i].es);
//      }
//      else {
//         remainder = remainder + mol.footprints[i].vdw;
//         remainder_ref = remainder_ref + footprint_reference.footprints[i].vdw;
//      }
//   }
//   if (fps_use_remainder){
//       fp_new.push_back(remainder);
//       footprint_ref_tmp.push_back(remainder_ref);
//   }
//   fps_es_num_resid = fp_new.size();
//}


// +++++++++++++++++++++++++++++++++++++++++
// This function calculates the FP similarity between two footprints.
// char E_type specifies which footprint is used (h = h-bond; v = vdw; e = es)
// char comparison_type specifies the comparison metric (E = Euclidean; P = Pearson)
double 
Footprint_Similarity_Score::calc_fp_similarity(DOCKMol & mol,char E_type,char comparison_type)
{
    int n = receptor.num_residues;
    string range;

    if (mol.footprints.size() != n ) {
        cout << "Error: There is an inconsistancy between"
             <<" number of resiudes in footprint and mol2." << endl;
        exit(1);
    }

    if (mol.footprints.size() != footprint_reference.footprints.size()){
        cout << "Error: The reference footprint does not have the same number" << endl
             << "of entries as the pose being evaluated." << endl
             << "Verify that the reference footprint text file has a line for" << endl
             << "every residue (including ions and water) in your target." << endl
             << "Program will terminate" << endl;
        exit(1);
    }

    // if thresholds are specified then calculate ranges for ref.
    if (fps_foot_specify_a_threshold)
        pose_ranges = range_satisfying_threshold_all(mol.footprints);

    // take union 
    vector <double> footprint_tmp;
    vector <double> footprint_ref_tmp;
    vector <string> footprint_str_tmp;
    double remainder,remainder_ref; // only used if fps_use_remainder; 

    //fps_fp_info = "";

    // if h, copy hbond footprints to tmp
    if (E_type == 'h'){

        if (fps_foot_specify_a_threshold || fps_foot_specify_a_range){
            range = union_of_ranges(pose_ranges.fps_range_hb,ref_ranges.fps_range_hb);
            print_range(pose_ranges.fps_range_hb,ref_ranges.fps_range_hb,range,"hb",&fps_fp_info);
        }
         
        footprint_tmp.clear(); footprint_ref_tmp.clear();
        footprint_str_tmp.clear();
        remainder = 0;
        remainder_ref = 0;
  
        for (int i=0;i<n;i++){
            if (on_range((i+1),range)){
                footprint_tmp.push_back(mol.footprints[i].hb);
                footprint_ref_tmp.push_back(footprint_reference.footprints[i].hb);
                footprint_str_tmp.push_back(mol.footprints[i].resname);
            }
            else {
                remainder = remainder + mol.footprints[i].hb;
                remainder_ref = remainder_ref + footprint_reference.footprints[i].hb;
            }
        }

        if (fps_use_remainder){
            footprint_tmp.push_back(remainder);
            footprint_ref_tmp.push_back(remainder_ref);
            footprint_str_tmp.push_back("remainder");
        }
        fps_hb_num_resid = footprint_tmp.size();
    }

    //if v, copy vdw footprints to tmp
    if (E_type == 'v'){

        if (fps_foot_specify_a_threshold || fps_foot_specify_a_range){
            range = union_of_ranges(pose_ranges.fps_range_vdw,ref_ranges.fps_range_vdw);
            print_range(pose_ranges.fps_range_vdw,ref_ranges.fps_range_vdw,range,"vdw",&fps_fp_info);
        }

        footprint_tmp.clear(); footprint_ref_tmp.clear();
        remainder = 0;
        remainder_ref = 0;
  
        for (int i=0;i<n;i++){
            if (on_range((i+1),range)){
                footprint_tmp.push_back(mol.footprints[i].vdw);
                footprint_ref_tmp.push_back(footprint_reference.footprints[i].vdw);
                footprint_str_tmp.push_back(mol.footprints[i].resname);
            }
            else {
                remainder = remainder + mol.footprints[i].vdw;
                remainder_ref = remainder_ref + footprint_reference.footprints[i].vdw;
            }
        }

        if (fps_use_remainder){
            footprint_tmp.push_back(remainder);
            footprint_ref_tmp.push_back(remainder_ref);
            footprint_str_tmp.push_back("remainder");
        }
        fps_vdw_num_resid = footprint_tmp.size();
    }

    //if e, copy es footprints to tmp
    if (E_type == 'e'){

        if (fps_foot_specify_a_threshold || fps_foot_specify_a_range){
            range = union_of_ranges(pose_ranges.fps_range_es,ref_ranges.fps_range_es);
            print_range(pose_ranges.fps_range_es,ref_ranges.fps_range_es,range,"es",&fps_fp_info);
        }

        footprint_tmp.clear(); footprint_ref_tmp.clear();
        remainder = 0;
        remainder_ref = 0;
  
        for (int i=0;i<n;i++){
            if (on_range((i+1),range)){
                footprint_tmp.push_back(mol.footprints[i].es);
                footprint_ref_tmp.push_back(footprint_reference.footprints[i].es);
                footprint_str_tmp.push_back(mol.footprints[i].resname);
            }
            else {
                remainder = remainder + mol.footprints[i].es;
                remainder_ref = remainder_ref + footprint_reference.footprints[i].es;
            }
        }

        if (fps_use_remainder){
            footprint_tmp.push_back(remainder);
            footprint_ref_tmp.push_back(remainder_ref);
            footprint_str_tmp.push_back("remainder");
        }
        fps_es_num_resid = footprint_tmp.size();
    }

    double value = -10000; // value could be distance or r 

    //begin to compute distances
    if (comparison_type == 'E'){
        if (fps_normalize_foot){
            Normalization(&footprint_ref_tmp);
            Normalization(&footprint_tmp);
        }
        double distance = Euclidean_distance(&footprint_ref_tmp, &footprint_tmp);
        value = distance;
    }

    if (comparison_type == 'P'){
        double r = Pearson_correlation(&footprint_ref_tmp, &footprint_tmp);
        value = r;
    }

    if (fps_foot_specify_a_threshold || fps_foot_specify_a_range)
        print_foot(&footprint_str_tmp,&footprint_ref_tmp,&footprint_tmp,&fps_fp_info);

    // set string to lenth 0 before leaving the function.
    footprint_tmp.clear();
    footprint_ref_tmp.clear();
    footprint_str_tmp.clear();

    return value;

} // end Footprint_Similarity_Score::calc_fp_similarity()


// +++++++++++++++++++++++++++++++++++++++++
// Assigns the FPS score to mol.current_score
bool
Footprint_Similarity_Score::compute_score(DOCKMol & mol)
{
    compute_footprint(mol);
    char *compare_type;
    compare_type = new char [fps_foot_compare_type.size()+1];
    strcpy (compare_type, fps_foot_compare_type.c_str());

    // pass only the frist char of compare_type to the function calc_fp_similarity
    fps_fp_info = "";
    vdw_foot_dist     = calc_fp_similarity(mol,'v',compare_type[0]);
    es_foot_dist      = calc_fp_similarity(mol,'e',compare_type[0]);
    hbond_foot_dist   = calc_fp_similarity(mol,'h',compare_type[0]);

    delete [] compare_type;

    mol.current_score = fp_vdw_scale*vdw_foot_dist + 
                        fp_es_scale*es_foot_dist + 
                        fp_hb_scale*hbond_foot_dist;

    mol.current_data = output_score_summary(mol);

    return true;

} // end Footprint_Similarity_Score::compute_score()


// +++++++++++++++++++++++++++++++++++++++++
// Write the footprint info to mol.current_data
string
Footprint_Similarity_Score::output_score_summary(DOCKMol & mol)
{
    ostringstream text;

    if (use_score) {

        //changed output score component names to manage length and be consistent.  Yuchen 10/24/2016
        text << DELIMITER << setw(STRING_WIDTH) << "Footprint_Similarity_Score:"
             << setw(FLOAT_WIDTH) << fixed << mol.current_score << endl
             << DELIMITER << setw(STRING_WIDTH) << "FPS_vdw_energy:"
             << setw(FLOAT_WIDTH) << fixed << vdw_component << endl
             << DELIMITER << setw(STRING_WIDTH) << "FPS_es_energy:"
             << setw(FLOAT_WIDTH) << fixed << es_component << endl
             << DELIMITER << setw(STRING_WIDTH) << "FPS_vdw+es_energy:"
             << setw(FLOAT_WIDTH) << fixed << vdw_component+es_component << endl
             << DELIMITER << setw(STRING_WIDTH) << "FPS_num_hbond:"
             << setw(FLOAT_WIDTH) << fixed << hbond  << endl;

        if (use_internal_energy)
            text << DELIMITER << setw(STRING_WIDTH) << "Internal_energy_repulsive:"
                 << setw(FLOAT_WIDTH) << fixed << mol.internal_energy << endl;

        text << DELIMITER << setw(STRING_WIDTH) << "FPS_vdw_fps:"
             << setw(FLOAT_WIDTH) << fixed << vdw_foot_dist << endl
             << DELIMITER << setw(STRING_WIDTH) << "FPS_es_fps:"
             << setw(FLOAT_WIDTH) << fixed << es_foot_dist << endl
             << DELIMITER << setw(STRING_WIDTH) << "FPS_hb_fps:"
             << setw(FLOAT_WIDTH) << fixed << hbond_foot_dist << endl
             << DELIMITER << setw(STRING_WIDTH) << "FPS_vdw_fp_numres:"
             << setw(FLOAT_WIDTH) << fixed << fps_vdw_num_resid << endl
             << DELIMITER << setw(STRING_WIDTH) << "FPS_es_fp_numres:"
             << setw(FLOAT_WIDTH) << fixed << fps_es_num_resid << endl
             << DELIMITER << setw(STRING_WIDTH) << "FPS_hb_fp_numres:"
             << setw(FLOAT_WIDTH) << fixed << fps_hb_num_resid  << endl;
         if (fps_fp_info != "") 
             text << fps_fp_info << endl;
    }

    return text.str();

} // end Footprint_Similarity_Score::output_score_summary()


