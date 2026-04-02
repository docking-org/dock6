#include <algorithm>
#include <iostream>
#include <iomanip>
#include <list>
#include <sstream>
#include <assert.h>
#include <string.h>
//#include <string>
#include "library_file.h"
#include "master_score.h"
#include "simplex.h"
#include "hungarian.h"
#include "trace.h"
#include "filter.h"
//#include <gzstream.h>
#include "gzstream/gzstream.h"
#include "xlogp.h"

// Pak and GDRM, October 2020 
#ifdef BUILD_DOCK_WITH_RDKIT
#include "rdtyper.h"
//#include <GraphMol/Descriptors/MolDescriptors.h>
//#include <GraphMol/SmilesParse/SmilesWrite.h>
#endif

using namespace std;


// +++++++++++++++++++++++++++++++++++++++++
// static member initializers

// These are the same as those in Base_Score.
const string Library_File::DELIMITER    = "########## ";
const int    Library_File::FLOAT_WIDTH  = 20;
const int    Library_File::STRING_WIDTH = 17 + 19;


/************************************************/
Library_File::Library_File()
    :   total_mols( 0 ),
        initial_skip( 0 )
{
}

/************************************************/
string
Library_File::hdb_atom_type_converter(int numeric_type)
{
    string          type;

    switch (numeric_type) {

    case 1:
        type = "C.";
        break;

    case 5:
        type = "C.3";
        break;

    case 6:
        type = "H";
        break;

    case 7:
        type = "H";
        break;

    case 8:
        type = "N.";
        break;

    case 9:
        type = "N.4";
        break;

    case 10:
        type = "N.3";
        break;

    case 11:
        type = "O.";
        break;

    case 12:
        type = "O.3";
        break;

    case 13:
        type = "P.";            // possible degeneracy with P. or P2
        break;

    case 14:
        type = "S.";
        break;

    case 15:
        type = "F";
        break;

    case 16:
        type = "Cl";
        break;

    case 17:
        type = "Br";
        break;

    case 18:
        type = "I";
        break;

    case 19:
        type = "Na";            // possible degeneracy with Na or K
        break;

    case 20:
        type = "Li";            // possible degeneracy with Li, Al or B
        break;

    case 21:
        type = "Ca";
        break;

    case 24:
        type = "Si";
        break;

    case 25:
        type = "Du";
        break;

    default:
        type = "Du";
        cout << "Atom type conversion not found for atom type " << numeric_type
            << ".  Du type assumed." << endl;
        break;
    }

    return type;

}

/************************************************/
/* 
 * Reading in db2 formate hierarchys . . .
 * info from : https://sites.google.com/site/dock37wiki/home/mol2db2-format
 * Mol2db2 Format 2
 * File Format
 *
 * T type information (implicitly assumed)
 * M molecule (4 lines req'd, after that they are optional, 24 lines max)
 * A atoms
 * B bond
 * X xyz
 * R rigid xyz for matching (can actually be any xyzs)
 * C conformation
 * S sets
 * D clusters
 * E end of molecule
 * T ## namexxxx (implicitly assumed to be the standard 7)
 * M zincname protname #atoms #bonds #xyz #confs #sets #rigid #Mlines #clusters
 * M charge polar_solv apolar_solv total_solv surface_area
 * M smiles
 * M longname
 * [M arbitrary information preserved for writing out]
 * A stuff about each atom, 1 per line 
 * B stuff about each bond, 1 per line
 * X coordnum atomnum confnum x y z 
 * R rigidnum color x y z
 * C confnum coordstart coordend
 * S setnum #lines #confs_total broken hydrogens omega_energy
 * S setnum linenum #confs confs [until full column]
 * D clusternum setstart setend matchstart matchend #additionalmatching
 * D matchnum color x y z
 * E 
 * */
bool
//Library_File::read_hierarchy_db2(DOCKMol & mol, ifstream & ifs)
//Library_File::read_hierarchy_db2(string filename, HDB_Mol & db2_data)
//Library_File::read_hierarchy_db2(ifstream & db2, HDB_Mol & db2_data)
Library_File::read_hierarchy_db2(igzstream & db2, HDB_Mol & db2_data)
{
    cout << "Library_File::read_hierarchy_db2 entering..." << endl;
    char line[1000];
    char zname[100];
    char pname[100];
    int  num_atoms,num_bonds,num_xyz,num_seg,num_branch,num_rigid,num_M,num_cluster;
    int count = 1;
    int count_atom   = 0;
    int count_bond   = 0;
    int count_coord  = 0;
    int count_rigid  = 0;
    int count_seg    = 0;
    int count_branch = 0;
    int setline = 1; // set also a branch. 
    int numofsetline; // read in this from the frist line
    //bool flag_header = false;
    bool flag_atom   = false;
    bool flag_bond   = false;
    bool flag_coord  = false;
    bool flag_rigid  = false;
    bool flag_seg    = false;
    bool flag_branch = false;
    int count_segs_in_conf=0; // we need to fill the list of segments in each conformation with will help us keep track. 

    char s_char = 'S';
    char d_char = 'D';

    ranked_poses.clear();
    num_confs = 0;

    while (db2.getline(line, 1000)) {
        //cout << line << endl;
        char first_char = line[0];
        if (count == 1) {
            //cout << line << endl;
            sscanf(line, "M %s %s %d %d %d %d %d %d %d %d", &zname, &pname, &num_atoms, &num_bonds, &num_xyz, &num_seg, &num_branch, &num_rigid, &num_M, &num_cluster);
            // flag_header = true;
            //db2_data.name          = zname ;
            strncpy(db2_data.name, zname, 100);
            db2_data.num_of_atoms  = num_atoms;
            db2_data.num_of_bonds  = num_bonds;
            db2_data.num_of_coords = num_xyz;
            db2_data.num_of_rigid  = num_rigid;
            db2_data.num_of_segs   = num_seg;  
            db2_data.num_of_confs  = num_branch; //branches
            db2_data.atoms = new HDB_atom [num_atoms];
            db2_data.bonds = new HDB_bond [num_bonds];
            db2_data.coords = new HDB_coordinate [num_xyz];
            db2_data.rigid = new HDB_coordinate [num_rigid];
            db2_data.segs = new HDB_segment [num_seg];
            db2_data.confs = new HDB_conformer [num_branch]; // num_branch is the number of lines and this is an upper bond on the accual number. 
            
            //break;
        }
        if (count == num_M+1){//start atom
            //cout << num_M+1 << endl;
            //cout << "Read in atoms" << endl;
            //cout << line << endl;
            flag_atom = true; flag_bond   = false; flag_coord  = false; flag_rigid  = false; flag_seg    = false; flag_branch = false;
            //break;
        }
        if (count == num_M+num_atoms+1){//start bonds
            //cout << num_M+num_atoms+1 << endl;
            //cout << "Read in bonds.  " << endl;
            //cout << line << endl;
            flag_atom = false; flag_bond   = true; flag_coord  = false; flag_rigid  = false; flag_seg    = false; flag_branch = false;
        }
        if (count == num_M+num_atoms+num_bonds+1){//start coord
            //cout << num_M+num_atoms+num_bonds+1 << endl;
            //cout << "Read in coord.  " << endl;
            //cout << line << endl;
            flag_atom = false; flag_bond   = false; flag_coord  = true; flag_rigid  = false; flag_seg    = false; flag_branch = false;
        }
        if (count == num_M+num_atoms+num_bonds+num_xyz+1){//start rigid.  
            //cout << "Read in rigid.  " << endl;
            //cout << line << endl;
            flag_atom = false; flag_bond   = false; flag_coord  = false; flag_rigid  = true; flag_seg    = false; flag_branch = false;
        }
        if (count == num_M+num_atoms+num_bonds+num_xyz+num_rigid+1){//start segments
            //cout << num_M+num_atoms+num_bonds+num_xyz+1 << endl;
            //cout << "Read in segments.  " << endl;
            //cout << line << endl;
            flag_atom = false; flag_bond   = false; flag_coord  = false; flag_rigid  = false; flag_seg    = true; flag_branch = false;
        }
        //if (count == num_M+num_atoms+num_bonds+num_xyz+num_rigid+num_seg+1){//start branches.  
        //if (strcmp(line[0],'S')==0) {
        if (first_char==s_char) {
            //cout << num_M+num_atoms+num_bonds+num_xyz+num_seg+1 << endl;
            //cout << "Read in branches.  " << endl;
            //cout << line << endl;
            flag_atom = false; flag_bond   = false; flag_coord  = false; flag_rigid  = false; flag_seg    = false; flag_branch = true;
        }
        //if (count == num_M+num_atoms+num_bonds+num_xyz+num_rigid+num_seg+num_branch+1) {
        //if (strcmp(line[0],'D')==0) {
        if (first_char==d_char) {
            flag_branch = false;
            //here is where we would read in clusters. but we aren't for the time being. 
        }
        if (first_char=='E') { // at end of database.
            return true;
        }
        // not reading in clusters for now. 
        if (flag_atom){
        // example: 
        //A   1 C1   C.1    1  7   -0.1700     -0.300     +0.570     +0.270    19.250
        //A   2 C2   C.1    1  7   -0.1500     -1.050     +0.710     -0.340    13.240
        //
            int tmp_atom_num, tmp_docktype, tmp_seg_num;
            char tmp_atom_name[10], tmp_atom_sybyl[10];
            string tmp_s_atom_name, tmp_s_atom_sybyl;
            float tmp_crg, tmp_pdesolv, tmp_adesolv, tmp_tot_desolv, tmp_sasa;
            sscanf(line, "A %d %s %s %d %d %f %f %f %f %f", &tmp_atom_num, &tmp_atom_name, &tmp_atom_sybyl, &tmp_docktype, &tmp_seg_num, &tmp_crg, &tmp_pdesolv, &tmp_adesolv, &tmp_tot_desolv, &tmp_sasa);
            //cout <<  tmp_atom_num <<" "<< tmp_atom_name <<" "<< tmp_atom_sybyl <<" "<< tmp_docktype <<" "<< tmp_seg_num <<" "<< tmp_crg <<" "<< tmp_pdesolv <<" "<< tmp_adesolv <<" "<< tmp_tot_desolv <<" "<< tmp_sasa << endl;
            tmp_s_atom_name = string(tmp_atom_name); tmp_s_atom_sybyl = string(tmp_atom_sybyl);
            
            //db2_data.atoms[count_atom].initialize(tmp_atom_num,tmp_s_atom_name,tmp_s_atom_sybyl,tmp_docktype,tmp_seg_num,tmp_crg,tmp_pdesolv,tmp_adesolv,tmp_tot_desolv,tmp_sasa);
            db2_data.atoms[count_atom].initialize(tmp_atom_num,tmp_s_atom_name,tmp_s_atom_sybyl,tmp_docktype,tmp_seg_num,tmp_crg,tmp_pdesolv,tmp_adesolv,tmp_sasa);
            
            count_atom++;
        }else if (flag_bond){
        //B   1   1   2 3
        //B   2   1  30 1
        //
           int bond_num, atom_num1, atom_num2;
           char bond_type[3]; 
           sscanf(line, "B %d %d %d %s", &bond_num, &atom_num1, &atom_num2, &bond_type);
           //cout << bond_num <<" "<< atom_num1 <<" "<< atom_num2 <<" "<< bond_type << endl;
           db2_data.bonds[count_bond].initialize(bond_num,atom_num1,atom_num2,bond_type);
           count_bond++;
        }else if (flag_coord){
        // X         1   1      1   -3.1875   +4.1115   -2.0685
        // X         2   2      1   -2.2603   +3.7291   -1.4583
           int coord_num, atom_num, seg_num; 
           float x,y,z;
           sscanf(line, "X %d %d %d %f %f %f", &coord_num, &atom_num, &seg_num, &x, &y, &z);
           //cout << coord_num << " " << atom_num << " " << seg_num << " " << x << " " << y << " " << z << endl;
           db2_data.coords[count_coord].initialize(coord_num,atom_num,seg_num,x,y,z);
           count_coord++;
        }else if (flag_rigid){
        //R      1  7   -3.1875   +4.1115   -2.0685
        //R      2  7   -2.2603   +3.7291   -1.4583
           int rig_num, t_num;
           float x, y, z;
           sscanf(line, "R %d %d %f %f %f", &rig_num, &t_num, &x, &y, &z);
           //cout << rig_num << " " << t_num << " " << x << " " << y << " " << z << endl;
           db2_data.rigid[count_rigid].initialize(rig_num,t_num,1,x,y,z);
           count_rigid++;
        }else if (flag_seg){
        //C      1         1        14
        //C      2        15        18
        //C      3        19        32
           int seg_num, coord_num_start, coord_num_stop;
           sscanf(line, "C %d %d %d", &seg_num, &coord_num_start, &coord_num_stop);  
           //cout << seg_num << " " << coord_num_start << " " << coord_num_stop << endl;
           db2_data.segs[count_seg].initialize(seg_num, coord_num_start, coord_num_stop);
           count_seg++;
        }else if (flag_branch){
        //S      1      3  17 0 0      +0.000
        //S      1      1 8      1    193    194    195    271    311    312    313
        //S      1      2 8    387    477    478    719    720    721    854   1417
        //S      1      3 1   1500
        //S      2      3  17 0 0      +0.020
        //S      2      1 8      1     13    145    146    147    457    501    502
        //S      2      2 8    503    682    683    684    927   1169   1178   1179
        //S      2      3 1   1253
           int set_num, num_of_seg, temp1, temp2;
           float internal_energy;
           int t_set_num, t_line_count, num_of_seg_line, t_seg_c_num;
           string stemp1, stemp2, stemp3, stemp4;
           //cout << "I AM HERE in seg" << endl;
           if (setline == 1) {
             sscanf(line, "S %d %d %d %d %d %f", &set_num, &numofsetline, &num_of_seg, &temp1, &temp2, &internal_energy);
             //cout << "Line 1. " << set_num << " " << numofsetline << " " << num_of_seg << " " << temp1 << " " << temp2 << " " << internal_energy << endl;
             db2_data.confs[count_branch].initialize(set_num, num_of_seg, internal_energy);
             //cout << "I AM HERE in seg" << endl;
           } else {
             //cout << "I AM HERE in seg" << endl;
             sscanf(line, "S %d %d %d", &t_set_num, &t_line_count, &num_of_seg_line);
             //cout << "Line 2+. " << t_set_num << " " << t_line_count << " " << num_of_seg_line << endl;
             stringstream ss, ss2;
             ss << line; // put string in stringstream
             // skip over the first 4 entries (the S, set_num, line count, snf num of segments in the line).   
             ss >> stemp1 >> stemp2 >> stemp3 >> stemp4 ;  
             //cout <<" "<< stemp1 <<" "<< stemp2 <<" "<< stemp3 <<" "<< stemp4 << endl;
             //cout << num_of_seg_line << endl;
             for ( int lc=0; lc<num_of_seg_line;lc++){
                  ss >> t_seg_c_num;
                  //cout << "conf segs" << lc << " "<<count_segs_in_conf<<" "<< t_seg_c_num << endl;
                  db2_data.confs[count_branch].list_of_seg[count_segs_in_conf] = t_seg_c_num; 
                  count_segs_in_conf++;
             }
             if (setline == numofsetline+1) {//go to new set (ie banch)
                 setline = 0;
                 count_segs_in_conf=0;
                 count_branch++;
                 //exit(0);
             }
           } 
           setline++; 
        } 
        //else {
        //    cout << "Something is wrong in read_hierarchy_db2 " << count << endl;
        //}

        count++;
    }
    cout << "Library_File::read_hierarchy_db2 exiting..." << endl;
    return false;
}
 
/************************************************/
// bool Library_File::read_hierarchy_db(vector<DOCKMol> &mol_branches, ifstream 
// &ifs) {
bool
Library_File::read_hierarchy_db(DOCKMol & mol, ifstream & ifs)
{
    int             i,
                    j,
                    k,
                    l;
    char            line[1000];
    string          tmp;
    int             total_atoms;
    vector < DOCKMol > mol_branches;

    // family level data
    int             clunum,
                    nmol,
                    nbr,
                    fbr;

    // molecule level data
    string          name,
                    refcode;

    // branch level data
    int             beatm,
                    bmatm,
                    bnhvy,
                    bnhyd,
                    confnum,
                    bi,
                    iconf;
    float           solvat,
                    apol;
    int             int_color;  // kxr

    // new molecule data
    DOCKMol         new_mol;
    FLOATVec        charges;
    vector < string > types;
    vector < string > color;    // kxr
    vector < string > summary_lines;
    string          family_line,
                    name_line,
                    branch_line;
    INTVec          level_counts,
                    level_values,
                    level_offsets;
    int             current_level,
                    current_count,
                    index;
    INTVec          atom_levels;
    int             next_level,
                    last_level;

    mol_branches.clear();

    // if file has no lines, return false
    if (!ifs.getline(line, 1000))
        return false;

    do {

        if (!strncmp(line, "Family", 6)) {

            // read vars from Family line
            sscanf(line, "%*s %d %d %d %d", &clunum, &nmol, &nbr, &fbr);
            family_line = line;

            // hack to bypass nbr and fbr = 0 kxr
            if (fbr < 5) {
                // cout << "Error: records are zero. Skipping to next entry" << 
                // endl;
                // ;
                // }

                nbr++;
                fbr++;

                // loop over nmols
                for (i = 0; i < nmol; i++) {

                    // ////////////////////
                    // skip first MBR loop
                    // ////////////////////

                    // read in mol name & refcode
                    if (!ifs.getline(line, 1000)) {
                        cout << "Error: Incomplete hierarchy molecule record."
                            << endl;
                        return false;
                    }

                    name_line = line;
                    name = name_line.substr(0, 47);
                    refcode = name_line.substr(47, 9);

                    // /////////////////////
                    // skip the other MBR variables
                    // /////////////////////

                    // loop over branches
                    for (j = 0; j < fbr; j++) {

                        // read branch header line
                        if (!ifs.getline(line, 1000)) {
                            cout <<
                                "Error: Incomplete hierarchy molecule record."
                                << endl;
                            return false;
                        }

                        sscanf(line, "%d %d %d %d %f %d %f %d %d", &beatm,
                               &bmatm, &bnhvy, &bnhyd, &solvat, &confnum, &apol,
                               &bi, &iconf);
                        branch_line = line;

                        // create new DOCKMol to store each branch
                        new_mol.clear_molecule();
                        new_mol.allocate_arrays(beatm, 0);
                        // new_mol.title = name;
                        new_mol.title = refcode;

                        // assign heirarchy data to mol_data string
                        if (j == 0) {
                            new_mol.mol_data = family_line + "\n";
                            new_mol.mol_data += name_line + "\n";
                        }
                        new_mol.mol_data += branch_line;

                        // allocate heirarchy data structures
                        charges.resize(bmatm);
                        types.resize(bmatm);
                        color.resize(bmatm);
                        summary_lines.resize(bmatm);

                        level_counts.clear();
                        level_counts.push_back(0);
                        level_values.clear();
                        level_values.push_back(0);

                        atom_levels.resize(beatm, 0);

                        current_level = 0;
                        next_level = 0;
                        last_level = 0;

                        // read the info into infovec struct
                        for (k = 0; k < bmatm; k++) {
                            ifs.getline(line, 1000);
                            tmp = line;

                            charges[k] = atof(tmp.substr(5, 5).c_str()) * 0.001;
                            int_color = atoi(tmp.substr(3, 2).c_str());
                            types[k] = hdb_atom_type_converter(int_color);

                            // hard code ligand coloring for now kxr
                            switch (int_color) {

                            case 1:
                                color[k] = "positive";
                                break;

                            case 2:
                                color[k] = "negative";
                                break;

                            case 3:
                                color[k] = "acceptor";
                                break;

                            case 4:
                                color[k] = "donor";
                                break;

                            case 5:
                                color[k] = "ester_o";
                                break;

                            case 6:
                                color[k] = "amide_o";
                                break;

                            case 7:
                                color[k] = "neutral";
                                break;

                            default:
                                color[k] = "null";
                                break;
                            }



                            // cout << " color " << color[k] << endl;
                            summary_lines[k] = "";
                            summary_lines[k] += tmp.substr(0, 3) + "\n";        // level
                            summary_lines[k] += tmp.substr(3, 2) + "\n";        // vdwtype
                            sprintf(line, "%f\n", charges[k]);
                            summary_lines[k] += line;   // charge
                            summary_lines[k] += tmp.substr(10, 1) + "\n";       // flagat
                            summary_lines[k] += tmp.substr(11, 2) + "\n";       // lcolor
                            summary_lines[k] += tmp.substr(13, 9) + "\n";       // polsolv
                            summary_lines[k] += tmp.substr(22, 9) + "\n";       // apolsolv
                            summary_lines[k] += types[k] + "\n";        // converted 
                                                                        // type

                            next_level = atoi(tmp.substr(0, 3).c_str());

                            // count the number of atoms in each level
                            if (k > 0) {
                                if (next_level != last_level) { // !=, >
                                    level_counts.push_back(0);
                                    level_values.push_back(0);
                                    current_level++;
                                    last_level = next_level;
                                }
                            } else
                                last_level = next_level;

                            level_values[current_level] = next_level;
                            level_counts[current_level]++;
                        }

                        // calculate the level offsets
                        level_offsets.clear();
                        level_offsets.resize(level_counts.size(), 0);

                        for (k = 1; k < level_offsets.size(); k++) {
                            level_offsets[k] =
                                level_offsets[k - 1] + level_counts[k - 1];
                        }

                        // add the level count information to mol_data
                        new_mol.mol_data += "\n";

                        for (k = 0; k < level_counts.size(); k++) {
                            sprintf(line, "%d %d,", level_counts[k],
                                    level_values[k]);
                            new_mol.mol_data += line;
                        }

                        current_level = 0;
                        current_count = 0;

                        // read coords for all confs
                        for (k = 0; k < beatm; k++) {

                            // read one line
                            ifs.getline(line, 1000);
                            tmp = line;

                            // assign coordinates
                            new_mol.x[k] =
                                atoi(tmp.substr(3, 6).c_str()) * 0.001;
                            new_mol.y[k] =
                                atoi(tmp.substr(9, 6).c_str()) * 0.001;
                            new_mol.z[k] =
                                atoi(tmp.substr(15, 6).c_str()) * 0.001;

                            atom_levels[k] = atoi(tmp.substr(0, 3).c_str());

                            // scan for level transitions
                            if (k > 0) {

                                // if level has increased
                                if (atom_levels[k] > atom_levels[k - 1]) {

                                    // cout <<
                                    // ".";////////////////////////////////////////

                                    current_level = 0;

                                    for (l = 0; l < level_values.size(); l++) {

                                        // cout << level_values[l] << "
                                        // ";//////////////

                                        if (atom_levels[k] == level_values[l]) {
                                            current_level = l;
                                            break;
                                        }
                                    }

                                    current_count = 0;
                                }
                                // if level has decreased
                                if (atom_levels[k] < atom_levels[k - 1]) {
                                    current_level = 0;

                                    for (l = 0; l < level_values.size(); l++) {
                                        if (atom_levels[k] == level_values[l]) {
                                            current_level = l;
                                            break;
                                        }
                                    }

                                    current_count = 0;
                                }

                                if (current_count >=
                                    level_counts[current_level])
                                    current_count = 0;
                            }
                            // cout << current_level <<
                            // endl;/////////////////////////////////

                            // apply vdw & charge values from info section to
                            // atom records
                            index =
                                level_offsets[current_level] + current_count;

                            new_mol.charges[k] = charges[index];
                            new_mol.atom_data[k] = summary_lines[index];
                            new_mol.atom_types[k] = types[index];
                            new_mol.atom_names[k] = types[index];
                            new_mol.subst_names[k] = types[index];
                            new_mol.atom_color[k] = color[index];
                            current_count++;

                        }

                        // assign branches of molecule to vector of mols
                        mol_branches.push_back(new_mol);

                    }

                }

                // merge molecules into one total mol

                // count all the atoms
                total_atoms = 0;
                for (i = 0; i < mol_branches.size(); i++)
                    total_atoms += mol_branches[i].num_atoms;

                // create new hybrid molecule with all branches
                mol.clear_molecule();
                mol.allocate_arrays(total_atoms, 0);

                index = 0;

                for (i = 0; i < mol_branches.size(); i++) {

                    for (j = 0; j < mol_branches[i].num_atoms; j++) {

                        // copy molecule atom data
                        mol.x[index] = mol_branches[i].x[j];
                        mol.y[index] = mol_branches[i].y[j];
                        mol.z[index] = mol_branches[i].z[j];

                        mol.charges[index] = mol_branches[i].charges[j];
                        mol.atom_types[index] = mol_branches[i].atom_types[j];
                        mol.atom_names[index] = mol_branches[i].atom_names[j];
                        mol.subst_names[index] = mol_branches[i].subst_names[j];
                        mol.atom_color[index] = mol_branches[i].atom_color[j];
                        // cout << " mol_atom_color " << mol.atom_color[index]
                        // << endl;

                        mol.atom_data[index] = mol_branches[i].atom_data[j];

                        // no bonds, hence no rings
                        mol.atom_ring_flags[index] = false;

                        // activate atoms
                        mol.atom_active_flags[index] = true;    // /////////////

                        index++;

                    }

                    mol.mol_data += mol_branches[i].mol_data;
                    mol.mol_data += "\n";

                }

                mol.num_active_atoms = total_atoms;

                mol_branches.clear();
                // return true if mol is found
                return true;
            }

        }

    } while (ifs.getline(line, 1000));

    // if no mol found, return false
    return false;

}

/************************************************/
void
Library_File::input_parameters_input(Parameter_Reader & parm)
{
    string          tmp;

    cout << "\nMolecule Library Input Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    input_file_name = parm.query_param("ligand_atom_file", "database.mol2");
   
    max_mol_limit =
        (parm.query_param("limit_max_ligands", "no", "yes no") ==
         "yes") ? true : false;
    if (max_mol_limit) {
        max_mols = atoi(parm.query_param("max_ligands", "1000").c_str());
        if (max_mols <= 0) {
            cout <<
                "ERROR:  Parameter must be an integer greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }
    } else
        max_mols = -1;
    
    bool            skip_molecule;
    skip_molecule =
        (parm.query_param("skip_molecule", "no", "yes no") ==
         "yes") ? true : false;
    if (skip_molecule) {
        initial_skip = atoi(parm.query_param("initial_skip", "0").c_str());
        if (initial_skip <= 0) {
            cout <<
                "ERROR: Parameter must be an integer greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }
    } else {
        initial_skip = 0;
    }

    read_solvation =
        parm.query_param("read_mol_solvation", "no", "yes no") == "yes";
    // read_color = parm.query_param("read_mol_color","no","yes no") == "yes";
    read_color = false;         // silenced keyword for now

    calc_rmsd = parm.query_param("calculate_rmsd", "no", "yes no") == "yes";
    if (calc_rmsd) {
        constant_rmsd_ref =
            parm.query_param("use_rmsd_reference_mol", "no", "yes no") == "yes";
        if (constant_rmsd_ref) {
            constant_rmsd_ref_file =
                parm.query_param("rmsd_reference_filename", "ligand_rmsd.mol2");
        }
    } else {
        constant_rmsd_ref = false;
    }
}

/************************************************/
void
Library_File::input_parameters_output(Parameter_Reader & parm, Master_Score & score, bool USE_MPI)
{
    Trace trace( "Library_File::input_parameters_output" );
    max_ranked = 0;
    rank_secondary_ligands = false;
    max_secondary_ranked = 0;
    use_secondary_score = false;
    cluster_ranked_poses = false;
    write_conformers = false;
    write_secondary_conformers = false;
    num_scored_poses = 1;
    num_clusterheads_rescore = 1;

    db2flag = false;

    if (score.use_secondary_score)
        use_secondary_score = true;

    cout << "\nMolecule Library Output Parameters" << endl;
    cout <<
        "------------------------------------------------------------------------------------------"
        << endl;

    output_file_prefix = parm.query_param("ligand_outfile_prefix", "output");

    trace.boolean( "score.use_score", score.use_score );
    if (score.use_score &&
        (score.c_fps.use_score ||
         score.c_mg_nrg.use_score ||
         (score.c_desc.use_score && score.c_desc.desc_use_fps && score.c_desc.desc_c_fps.use_score) ||
         (score.c_desc.use_score && score.c_desc.desc_use_mg_nrg && score.c_desc.desc_c_mg_nrg.use_score))) {
        trace.boolean( "score.c_fps.use_score", score.c_fps.use_score );
        trace.boolean( "score.c_mg_nrg.use_score", score.c_mg_nrg.use_score );
        trace.boolean( "score.c_desc.use_score", score.c_desc.use_score );
        write_footprints = (parm.query_param("write_footprints",
                                            "no", "yes no") =="yes");
        if (score.use_score && (score.c_fps.use_score ||(score.c_desc.use_score
            && score.c_desc.desc_use_fps && score.c_desc.desc_c_fps.use_score)))
            write_hbonds = (parm.query_param("write_hbonds", "no", "yes no") =="yes");
        else write_hbonds = false;
    } else {
        write_footprints = false;
        write_hbonds = false;
    }
    write_orients = (parm.query_param("write_orientations", "no", "yes no") == "yes");
    if (write_orients){
        if (USE_MPI){
            cout << "ERROR:  To protect against filling up the disk, \n"
                 << "        DOCK cannot write orientations while running in parallel."
                 << endl;
            finalize_mpi();
            exit(0);
        }
    } 

    if(score.use_secondary_score){
        num_scored_poses =
            atoi(parm.query_param("num_primary_scored_conformers_rescored", "1").
                 c_str());
        if (num_scored_poses <= 0) {
            cout <<
                "ERROR:  Parameter must be an integer greater than zero.  Program will terminate."
                << endl;
            exit(0);
        }
        if (num_scored_poses > 1) {
                
            write_conformers =  (parm.query_param("write_primary_conformations", "yes", "yes no") ==
                 "yes") ? true : false;
                if (write_conformers){
                     if (USE_MPI){
                          cout << "ERROR:  To protect against filling up the disk, DOCK cannot write ";
                          cout << "conformations while running in parallel." << endl;
                          finalize_mpi();
                          exit(0);
                     }
                }

            cluster_ranked_poses =
                (parm.query_param("cluster_primary_conformations", "yes", "yes no") ==
                 "yes") ? true : false;

            if (cluster_ranked_poses) {
                cluster_rmsd_threshold = atof(parm.query_param(
                        "cluster_rmsd_threshold", "2.0").c_str());
                if (cluster_rmsd_threshold <= 0.0) {
                    cout << "ERROR:  Parameter must be a float greater than zero.  Program will terminate."
                        << endl;
                    exit(0);
                }
                num_clusterheads_rescore = atoi(parm.query_param("num_clusterheads_for_rescore", "5").c_str());
            }

        } 
        if(!rank_secondary_ligands) {
            num_secondary_scored_poses =
                atoi(parm.query_param("num_secondary_scored_conformers", "1").c_str());
            if(num_secondary_scored_poses > 1){
                write_secondary_conformers = (parm.query_param("write_secondary_conformations", "yes", "yes no") == "yes") ? true : false;
                if (write_secondary_conformers){
                    if (USE_MPI){
                         cout << "ERROR:  To protect against filling up the disk, DOCK cannot write ";
                         cout << "conformations while running in parallel." << endl;
                         finalize_mpi();
                         exit(0);
                    }
                }

                if(num_secondary_scored_poses > num_scored_poses){
                    cout << "ERROR:  Number of secondary poses written cannot exceed number of primary ";
                    cout << "poses rescored.  Program will terminate." << endl;
                    exit(0);
                }
                if (num_secondary_scored_poses <= 0) {
                    cout << "ERROR:  Parameter must be an integer greater than zero.  Program will terminate."
                    << endl;
                    exit(0);
                }
            }
        }
    } else {
        num_secondary_scored_poses = 0;
        num_scored_poses = atoi(parm.query_param("num_scored_conformers", "1")
                                .c_str());
        if (num_scored_poses <= 0) {
            cout << "ERROR:  Parameter must be an integer greater than zero."
                    "  Program will terminate."
                 << endl;
            exit(0);
        }
        if (num_scored_poses > 1) {
            write_conformers = (parm.query_param("write_conformations", "yes",
                    "yes no") == "yes");
            if (write_conformers){ 
                 if (USE_MPI){
                      cout << "ERROR:  To protect against filling up the disk, DOCK cannot write ";
                      cout << "conformations while running in parallel." << endl;
                      finalize_mpi();
                      exit(0);
                 }
            }

            cluster_ranked_poses = (parm.query_param("cluster_conformations",
                                    "yes", "yes no") == "yes");

            if (cluster_ranked_poses) {
                cluster_rmsd_threshold = atof(parm.query_param(
                        "cluster_rmsd_threshold", "2.0").c_str());
                if (cluster_rmsd_threshold <= 0.0) {
                    cout << "ERROR:  Parameter must be a float greater than zero."
                            "  Program will terminate."
                         << endl;
                    exit(0);
                }
            }
        } 
    }
    if(score.use_secondary_score){
        rank_ligands = (parm.query_param("rank_primary_ligands", "no", "yes no") == "yes");
        if(rank_ligands){
            max_ranked = atoi(parm.query_param("max_primary_ranked", "500").c_str());
            if (max_ranked <= 0) {
                cout <<
                    "ERROR:  Parameter must be integer greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }
        }
        rank_secondary_ligands = (parm.query_param("rank_secondary_ligands", "no", "yes no") == "yes");
        if(rank_secondary_ligands){
            max_secondary_ranked = atoi(parm.query_param("max_secondary_ranked", "500").c_str());
            if (max_secondary_ranked <= 0) {
                cout <<
                    "ERROR:  Parameter must be integer greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }
            // This restriction on numbers of poses should be removed;
            // note that it also implies that primary ranking is a prerequisite
            // for secondary ranking; srb.
            if(max_ranked <  max_secondary_ranked){
                cout << "ERROR:  Maximum number of poses rescored and reranked "
                                 "by secondary score must be\n"
                        "        no more than maximum number of poses scored "
                                 "and ranked by primary score."
                     << endl;
                exit(0);
            }
        }
    } else {
        rank_ligands = (parm.query_param("rank_ligands", "no", "yes no") == "yes");
        if (rank_ligands) {
            max_ranked = atoi(parm.query_param("max_ranked_ligands", "500").c_str());
            if (max_ranked <= 0) {
                cout <<
                    "ERROR:  Parameter must be integer greater than zero.  Program will terminate."
                    << endl;
                exit(0);
            }
        }
    }
}

/************************************************/
void
Library_File::initialize(int argc, char **argv, bool USE_MPI)
{
    cout << "Initializing Library File Routines...\n";
    total_mols = 0;
    completed = 0;

    output_file_orient = output_file_prefix + "_orients.mol2";
    if(use_secondary_score){
        output_file_confs = output_file_prefix + "_primary_conformers.mol2";
        output_file_scored = output_file_prefix + "_primary_scored.mol2";
        output_file_ranked = output_file_prefix + "_primary_ranked.mol2";
        if(write_footprints){
            output_file_footprint_confs = output_file_prefix + "_primary_footprint_conformers.txt";
            output_file_footprint_scored = output_file_prefix + "_primary_footprint_scored.txt";
            output_file_footprint_ranked = output_file_prefix + "_primary_footprint_ranked.txt";
        }
        if(write_hbonds){
            output_file_hbond_confs = output_file_prefix + "_primary_hbond_conformers.txt";
            output_file_hbond_scored = output_file_prefix + "_primary_hbond_scored.txt";
            output_file_hbond_ranked = output_file_prefix + "_primary_hbond_ranked.txt";
        }
        secondary_output_file_confs = output_file_prefix + "_secondary_conformers.mol2";
        secondary_output_file_scored = output_file_prefix + "_secondary_scored.mol2";
        secondary_output_file_ranked = output_file_prefix + "_secondary_ranked.mol2";
    } else {
        output_file_confs = output_file_prefix + "_conformers.mol2";
        output_file_scored = output_file_prefix + "_scored.mol2";
        output_file_ranked = output_file_prefix + "_ranked.mol2";
        #ifdef BUILD_DOCK_WITH_RDKIT 
        output_file_rejected = output_file_prefix + "_rejected.mol2";
        #endif
        if(write_footprints){
            output_file_footprint_confs = output_file_prefix + "_footprint_conformers.txt";
            output_file_footprint_scored = output_file_prefix + "_footprint_scored.txt";
            output_file_footprint_ranked = output_file_prefix + "_footprint_ranked.txt";
        }
        if(write_hbonds){
            output_file_hbond_confs = output_file_prefix + "_hbond_conformers.txt";
            output_file_hbond_scored = output_file_prefix + "_hbond_scored.txt";
            output_file_hbond_ranked = output_file_prefix + "_hbond_ranked.txt";
        }
    }
    //cout << "I AM HERE AAAAA" <<endl;
    //bestmol.clear_molecule();
    //bestscore = 0.0;

    if (rank_ligands)
        ranked_list.reserve(max_ranked);


    if ((!USE_MPI) || (is_master_node()))
        open_files();

    rank_list_highest_energy = 0.0;
    rank_list_highest_pos = -1;

    read_success = false;
}

/************************************************/
void
Library_File::open_files()
{
    DOCKMol         tmp_mol;
    int             i;
    string          extension;
    //cout << "I AM HERE in Library_File::open_files." << endl;
    // identify the input file type: MOL2 or HDB (if not clear- assume MOL2)
    if (!db2flag){
       extension =
           input_file_name.substr(input_file_name.rfind("."),
                                  input_file_name.size());
 
       if (extension == ".db") {
           db_type = 2;
       } else if (extension == ".mol2") {
           db_type = 1;
       } else {
           cout << "Warning: Filetype not recognized for " << input_file_name <<
               ".  Assuming MOL2 format" << endl;
           db_type = 1;
       }
 
       ligand_in.open(input_file_name.c_str());
    }
    if (ligand_in.fail()) {
        cout << "\n\nCould not open " << input_file_name <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    if (write_orients) {
        ligand_out_orients.open(output_file_orient.c_str());
        if (ligand_out_orients.fail()) {
            cout << "\n\nCould not open " << output_file_orient <<
                " for writing.  Program will terminate." << endl << endl;
            exit(0);
        }
    }

    if (write_conformers) {
        ligand_out_confs.open(output_file_confs.c_str());
        if (ligand_out_confs.fail()) {
            cout << "\n\nCould not open " << output_file_confs <<
                " for writing.  Program will terminate." << endl << endl;
            exit(0);
        }
        if (write_footprints) {
            ligand_out_footprint_confs.open(output_file_footprint_confs.c_str());
            if (ligand_out_footprint_confs.fail()) {
                cout << "\n\nCould not open " << output_file_footprint_confs <<
                    " for writing.  Program will terminate." << endl << endl;
                exit(0);
            }
        }
        if (write_hbonds) {
            ligand_out_hbond_confs.open(output_file_hbond_confs.c_str());
            if (ligand_out_hbond_confs.fail()) {
                cout << "\n\nCould not open " << output_file_hbond_confs <<
                    " for writing.  Program will terminate." << endl << endl;
                exit(0);
            }
        }
    }

    if(write_secondary_conformers) {
        secondary_ligand_out_confs.open(secondary_output_file_confs.c_str());
        if (secondary_ligand_out_confs.fail()) {
            cout << "\n\nCould not open " << secondary_output_file_confs <<
                " for writing.  Program will terminate." << endl << endl;
            exit(0);
        }
    }

    if(!rank_ligands){
        ligand_out_scored.open(output_file_scored.c_str());
        if (ligand_out_scored.fail()) {
            cout << "\n\nCould not open " << output_file_scored <<
                " for writing.  Program will terminate." << endl << endl;
            exit(0);
        }
        #ifdef BUILD_DOCK_WITH_RDKIT
        ligand_out_rejected.open(output_file_rejected.c_str());
        if (ligand_out_rejected.fail()) {
            cout << "\n\nCould not open " << output_file_rejected <<
                " for writing.  Program will terminate." << endl << endl;
            exit(0);
        } 
        #endif
        if (write_footprints) {
            ligand_out_footprint_scored.open(output_file_footprint_scored.c_str());
            if (ligand_out_footprint_scored.fail()) {
                cout << "\n\nCould not open " << output_file_footprint_scored <<
                    " for writing.  Program will terminate." << endl << endl;
                exit(0);
            }
        }
        if (write_hbonds) {
            ligand_out_hbond_scored.open(output_file_hbond_scored.c_str());
            if (ligand_out_hbond_scored.fail()) {
                cout << "\n\nCould not open " << output_file_hbond_scored <<
                    " for writing.  Program will terminate." << endl << endl;
                exit(0);
            }
        }
    }
    if(use_secondary_score){
        if(!rank_secondary_ligands){
            secondary_ligand_out_scored.open(secondary_output_file_scored.c_str());
            if (secondary_ligand_out_scored.fail()) {
                cout << "\n\nCould not open " << secondary_output_file_scored <<
                    " for writing.  Program will terminate." << endl << endl;
                exit(0);
            }
        }
    }

    if (rank_ligands) {
        ligand_out_ranked.open(output_file_ranked.c_str());
        if (ligand_out_ranked.fail()) {
            cout << "\n\nCould not open " << output_file_ranked <<
                " for writing.  Program will terminate." << endl << endl;
            exit(0);
        }
        if (write_footprints) {
            ligand_out_footprint_ranked.open(output_file_footprint_ranked.c_str());
            if (ligand_out_footprint_ranked.fail()) {
                cout << "\n\nCould not open " << output_file_footprint_ranked <<
                    " for writing.  Program will terminate." << endl << endl;
                exit(0);
            }
        }
        if (write_hbonds) {
            ligand_out_hbond_ranked.open(output_file_hbond_ranked.c_str());
            if (ligand_out_hbond_ranked.fail()) {
                cout << "\n\nCould not open " << output_file_hbond_ranked <<
                    " for writing.  Program will terminate." << endl << endl;
                exit(0);
            }
        }
    }

    if (rank_secondary_ligands) {
        secondary_ligand_out_ranked.open(secondary_output_file_ranked.c_str());
        if (secondary_ligand_out_ranked.fail()) {
            cout << "\n\nCould not open " << secondary_output_file_ranked <<
                " for writing.  Program will terminate." << endl << endl;
            exit(0);
        }
    }

    // read through the initial skip ligands
    int             successful_skip = 0;
    for (i = 0; i < initial_skip; i++) {
        if (read_mol(tmp_mol, false))
            successful_skip++;
        else {
            cout <<
                "ERROR:  initial_skip value exceeds number of ligands in database.  Program will terminate."
                << endl << endl;
            exit(0);
        }

    }

    if (initial_skip > 0)
        cout << " Skipped " << successful_skip << " ligands" << endl << endl;

}

/************************************************/
void
Library_File::close_files()
{

    ligand_in.close();

    ligand_out_orients.close();
    ligand_out_confs.close();
    ligand_out_scored.close();
    #ifdef BUILD_DOCK_WITH_RDKIT
    ligand_out_rejected.close();
    #endif
    
    ligand_out_ranked.close(); 
    if(write_footprints) {
        ligand_out_footprint_confs.close();
        ligand_out_footprint_ranked.close();
        ligand_out_footprint_scored.close();
    }
    if(write_hbonds) {
        ligand_out_hbond_confs.close();
        ligand_out_hbond_ranked.close();
        ligand_out_hbond_scored.close();
    }
    if(use_secondary_score){
        secondary_ligand_out_confs.close();
        secondary_ligand_out_ranked.close();
        secondary_ligand_out_scored.close();
    }

}

/************************************************/
bool
Library_File::read_mol(DOCKMol & mol, bool read_amber)
{

    mol.clear_molecule();

    if (db_type == 1) {         // if file is a MOL2 file

        if (max_mol_limit && (total_mols >= max_mols))
            return false;

        if (Read_Mol2(mol, ligand_in, read_color, read_solvation, read_amber)) {
            total_mols++;
            return true;
        } else {
            return false;
        }

    } else if (db_type == 2) {  // if file is a hierarchy db
/**
        // read a new mol if branch list is empty
        if(mol_branches.size() == 0) {
            if(!read_hierarchy_db(mol_branches, ligand_in))
                return false;
        }

        // take one branch and copy it to mol
        copy_molecule(mol, mol_branches[mol_branches.size()-1]);
        mol_branches.pop_back();

        // if last branch, count the molecule as read
        if(mol_branches.size() == 0)
            total_mols++;
**/

        if (read_hierarchy_db(mol, ligand_in)) {
            total_mols++;

            // fail if too many mols have been read
            if (max_mol_limit && (total_mols > max_mols))
                return false;
            else
                return true;

        } else {
            return false;
        }

    } else {

        return false;

    }

}

/************************************************/
void
Library_File::write_mol(DOCKMol & mol, ofstream & ofs)
{

    Write_Mol2(mol, ofs);

}



/************************************************/
void
Library_File::write_footprint_file(DOCKMol & pose, Master_Score & score, std::ofstream & ofs)
{

ofs.setf(ios::fixed,ios::floatfield); 

// Reads in a dockmol and output the footprint of ref and pose in the following way:
// resid  resnum  VDW_ref  es_ref  hb_ref  VDW_pose1  es_pose1  hb_pose1  

////////////////////////////////////////////
//  output information                    //
//  assosiated with each pose             //
////////////////////////////////////////////

 ofs << "####################################" << endl;
 ofs << "###  Molecule: " << pose.title   << endl;
 ofs << "####################################" << endl;

 ofs << pose.score_text_data
     << endl;

 if (calc_rmsd)
     ofs << calc_rmsd_string(pose);
 if(score.use_score)
     ofs << pose.current_data << endl;

const int num_residues = pose.footprints.size();

int width1 =  8;
int width2 = 12;
int width3 =  8;

////////////////////////////////////////////
//          output footprints             //
////////////////////////////////////////////
//TODO
//fix this so footprints are written to file when using FPS component of new Descriptor Score
/*
if (score.c_desc_nrg.use_primary_score){// || score.c_desc_nrg.use_secondary_score){
      for (int i = 0; i < score.c_desc_nrg.receptor.num_residues; i++){
        if ( pose.footprints[i].resname != " "){
            // output the reference footprint
            ofs << setw (width1)
                << pose.footprints[i].resname                             << setw (width1)
                << pose.footprints[i].resid                               << setw (width2)
                << score.c_desc_nrg.footprint_reference.footprints[i].vdw << setw (width2)
                << score.c_desc_nrg.footprint_reference.footprints[i].es  << setw (width3)
                << score.c_desc_nrg.footprint_reference.footprints[i].hb  << setw (width2);
            // output the pose footprint    
            ofs << pose.footprints[i].vdw << setw (width2)
                << pose.footprints[i].es << setw (width1)
                << pose.footprints[i].hb;
            ofs << endl;
        }
    }
}
*/
if (score.c_fps.use_primary_score){// || score.c_fps.use_secondary_score){

      ////////////////////////////////////////////
      //            write header                //
      ////////////////////////////////////////////
      ofs << setw (width1) << "resname" << setw (width1)<<"resid"<< setw (width2) <<"vdw_ref"<< setw (width2)
          << "es_ref" << setw (width3) <<"hb_ref" << setw (width2) << "vdw_pose"<< setw (width2)
          << "es_pose" << setw (width3) <<"hb_pose"<< endl;

      for (int i = 0; i < score.c_fps.receptor.num_residues; i++){
        if ( pose.footprints[i].resname != " "){
            // output the reference footprint
            ofs << setw (width1)
                << pose.footprints[i].resname                             << setw (width1)
                << pose.footprints[i].resid                               << setw (width2)
                << score.c_fps.footprint_reference.footprints[i].vdw << setw (width2)
                << score.c_fps.footprint_reference.footprints[i].es  << setw (width3)
                << score.c_fps.footprint_reference.footprints[i].hb  << setw (width2);
            // output the pose footprint    
            ofs << pose.footprints[i].vdw << setw (width2)
                << pose.footprints[i].es << setw (width1)
                << pose.footprints[i].hb;
            ofs << endl;
        }
    }
}

if (score.c_desc.desc_c_fps.use_score){// || score.c_fps.use_secondary_score){

      ////////////////////////////////////////////
      //            write header                //
      ////////////////////////////////////////////
      ofs << setw (width1) << "resname" << setw (width1)<<"resid"<< setw (width2) <<"vdw_ref"<< setw (width2)
          << "es_ref" << setw (width3) <<"hb_ref" << setw (width2) << "vdw_pose"<< setw (width2)
          << "es_pose" << setw (width3) <<"hb_pose"<< endl;

      for (int i = 0; i < score.c_desc.desc_c_fps.receptor.num_residues; i++){
        if ( pose.footprints[i].resname != " "){
            // output the reference footprint
            ofs << setw (width1)
                << pose.footprints[i].resname                             << setw (width1)
                << pose.footprints[i].resid                               << setw (width2)
                << score.c_desc.desc_c_fps.footprint_reference.footprints[i].vdw << setw (width2)
                << score.c_desc.desc_c_fps.footprint_reference.footprints[i].es  << setw (width3)
                << score.c_desc.desc_c_fps.footprint_reference.footprints[i].hb  << setw (width2);
            // output the pose footprint    
            ofs << pose.footprints[i].vdw << setw (width2)
                << pose.footprints[i].es << setw (width1)
                << pose.footprints[i].hb;
            ofs << endl;
        }
    }
}

if (score.c_mg_nrg.use_primary_score){

      ////////////////////////////////////////////
      //            write header                //
      ////////////////////////////////////////////
      ofs << setw (width1) << "resname" << setw (width1)<<"resid"<< setw (width2) <<"vdw_ref"<< setw (width2)
          << "es_ref" << setw (width2) << "vdw_pose"<< setw (width2) << "es_pose" << endl;

      for (int i = 0; i < score.c_mg_nrg.numgrids; i++){
        if ( pose.footprints[i].resname != " "){
            // output the reference footprint
            ofs << setw (width1)
                << pose.footprints[i].resname                             << setw (width1)
                << pose.footprints[i].resid                               << setw (width2)
                << score.c_mg_nrg.footprint_reference.footprints[i].vdw << setw (width2)
                << score.c_mg_nrg.footprint_reference.footprints[i].es  << setw (width2);
            // output the pose footprint
            ofs << pose.footprints[i].vdw << setw (width2)
                << pose.footprints[i].es;
            ofs << endl;
        }
    }
}

if (score.c_desc.desc_c_mg_nrg.use_score){

      ////////////////////////////////////////////
      //            write header                //
      ////////////////////////////////////////////
      ofs << setw (width1) << "resname" << setw (width1)<<"resid"<< setw (width2) <<"vdw_ref"<< setw (width2)
          << "es_ref" << setw (width2) << "vdw_pose"<< setw (width2) << "es_pose" << endl;

      for (int i = 0; i < score.c_desc.desc_c_mg_nrg.numgrids; i++){
        if ( pose.footprints[i].resname != " "){
            // output the reference footprint
            ofs << setw (width1)
                << pose.footprints[i].resname                             << setw (width1)
                << pose.footprints[i].resid                               << setw (width2)
                << score.c_desc.desc_c_mg_nrg.footprint_reference.footprints[i].vdw << setw (width2)
                << score.c_desc.desc_c_mg_nrg.footprint_reference.footprints[i].es  << setw (width2);
            // output the pose footprint
            ofs << pose.footprints[i].vdw << setw (width2)
                << pose.footprints[i].es;
            ofs << endl;
        }
    }
}

ofs << endl;

}

/************************************************/
void
Library_File::write_hbond_file(DOCKMol & pose, Master_Score & score, std::ofstream & ofs)
{
 ofs << "####################################" << endl;
 ofs << "###  Molecule: " << pose.title   << endl;
 ofs << "####################################" << endl;

 ofs << pose.score_text_data
     << endl;
 if (calc_rmsd)
     ofs << calc_rmsd_string(pose);
 if(score.use_score)
     ofs << pose.current_data << endl;

 ofs << pose.hbond_text_data
     << endl;
}

// TEB added this. 
//
void Library_File::sort_write(bool USE_FILT, bool USE_MPI,Master_Score & score, Simplex_Minimizer & min){
    if (ranked_poses.size() > 0) {
        // sort ranked_poses
        sort(ranked_poses.begin(), ranked_poses.end());

        // cluster ranked poses if requested
        if (cluster_ranked_poses) {
            cluster_list();
        }
        
        //write out data for conformers
        write_scored_poses(USE_FILT, USE_MPI, score );
        
        //analyze and write data for secondary scored conformers
        if(num_scored_poses > 1)
            submit_secondary_conformation(score, min);
    }

}

/************************************************/
bool
Library_File::get_mol(DOCKMol & mol,bool USE_FILT, bool USE_MPI, bool amber, AMBER_TYPER & typer, Master_Score & score, Simplex_Minimizer & min)
{
    RANKMol         tmp_mol;
    int             i;

    if (!USE_MPI) {             // single proc code
        if (read_success) {
            sort_write(USE_FILT, USE_MPI, score, min);
        }

        // pose ranking code
        ranked_poses.clear();
        high_position = 0;
        high_score = 0.0;

        if (read_mol(mol, amber)) {

            num_orients = 0;
            num_confs = 0;
            num_anchors = 0;

            if (calc_rmsd) {
                submit_rmsd_reference(mol, typer);
            }
            read_success = true;

            return true;
        } else {
            if (rank_ligands) {
                write_ranked_ligands(USE_FILT,score);
            }
            return false;
        }

    } else {
/***/
        if (is_master_node()) { // master node code
            while (listen_for_request()) {
                if (is_send_request()) {
                    while (continue_reading_mols()) {
                        if (read_mol(mol, amber)) {
                            add_to_send_queue(mol);
                            // if HDB mols in send queue, then set continue
                            // mols true to allow all to be sent
                            if (mol_branches.size() > 0)
                                continue_mols = true;

                            read_success = true;
                        } else {
                            continue_mols = false;
                        }

                    }           // End loop over work unit of mols
                    transmit_send_queue();

                }               // End send mol code

                if (is_receive_request()) {
                    receive_flag = true;

                    ranked_poses.clear();
                    while (get_from_recv_queue(tmp_mol)) {
                        ranked_poses.push_back(tmp_mol);

                    }           // End loop over received mols

                    // copy docking stats from MPI stats struct
                    num_anchors = stats[0];
                    num_orients = stats[1];
                    num_confs = stats[2];
                    if (read_success) {
                        sort_write(USE_FILT, USE_MPI, score, min);
                    }

                }               // End receive mol code

            }                   // End looping over MPI listener

            if (rank_ligands)
                write_ranked_ligands(USE_FILT,score);

            close_files();
            return false;

        } else if (is_client_node()) {  // client node code
            // sort ranked poses
            if (ranked_poses.size() > 0) {
                completed++;
                sort(ranked_poses.begin(), ranked_poses.end());
            }
            // add ranked poses to return list
            recv_queue.clear();
            for (i = 0; i < ranked_poses.size(); i++)
                add_to_recv_queue(ranked_poses[i]);

            // copy docking stats to MPI data struct
            stats[0] = num_anchors;
            stats[1] = num_orients;
            stats[2] = num_confs;
            // init pose ranking
            ranked_poses.clear();
            high_position = 0;
            high_score = 0.0;
            // get molecule from server
            if (get_from_send_queue(mol)) {

                num_anchors = 0;
                num_orients = 0;
                num_confs = 0;
                if (calc_rmsd)
                    submit_rmsd_reference(mol, typer);
                read_success = true;

                return true;
            } else
                return false;

        }
/***/
        return false;           // to squelch warnings
    }


}

/************************************************/
bool
Library_File::submit_orientation(DOCKMol & mol, Master_Score & score, bool orient_ligand)
{
    //Trace("Library_File::submit_orientation enter"); 
    //cout << "Library_File::submit_orientation enter" << endl; 
    bool valid_orientation = false;
    num_orients++;
    if (orient_ligand) {
        // compute the score for the molecule
        if (score.use_score) {
            if (score.use_primary_score) 
                valid_orientation = score.compute_primary_score(mol);
                // in the case for grid, returns false if outside grid.
                if (!valid_orientation) 
                    return false;
        } else { // no scoring function spesified then
                 // all orientations are valid.
            valid_orientation = true;
        }
    }

    if (write_orients) {

        ligand_out_orients << "\n" << DELIMITER << setw(STRING_WIDTH) << "Name:" 
            << setw(FLOAT_WIDTH) << mol.title << "_" << num_orients << endl;

        if (calc_rmsd) 
            ligand_out_orients << calc_rmsd_string(mol);

        if(score.use_score)
            ligand_out_orients << mol.current_data << endl;
        
        write_mol(mol, ligand_out_orients);

    }
    //Trace ("Library_File::submit_orientation exit"); 
    //cout << "Library_File::submit_orientation exit" << endl; 
    if (!valid_orientation) {
        return false;
    }
    return true;

}

/************************************************/
void
Library_File::submit_conformation(Master_Score & score)
{

    if (write_conformers) {
    
        if(score.use_score)
            sort(ranked_poses.begin(), ranked_poses.end());

        for(int i=0;i<ranked_poses.size();i++){
            ligand_out_confs << "\n" << DELIMITER << setw(STRING_WIDTH) << "Name:" 
                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.title << "_" << i + 1 << endl;

            if (calc_rmsd)
                ligand_out_confs << calc_rmsd_string(ranked_poses[i].mol);
            if(score.use_score)
                ligand_out_confs << ranked_poses[i].mol.current_data << endl;

            write_mol(ranked_poses[i].mol, ligand_out_confs);
            if (write_footprints)
               write_footprint_file(ranked_poses[i].mol, score, ligand_out_footprint_confs);
            if (write_hbonds)
               write_hbond_file(ranked_poses[i].mol, score, ligand_out_hbond_confs);

        }

    }
    
}
/************************************************/
void
Library_File::submit_secondary_conformation(Master_Score & score, Simplex_Minimizer & min)
{
    Trace trace( "Library_File::submit_secondary_conformation" );
    int i;

    //reduce list to maximum clusterheads
    if(cluster_ranked_poses){
        for( i = 0; i < ranked_poses.size(); ) {
            if( cluster_assignments[i] == 0 ) {
                ranked_poses.erase( ranked_poses.begin() + i );
                cluster_assignments.erase( cluster_assignments.begin() + i );
            }
            else {
                ++i;
            }
        }
        while(ranked_poses.size() > num_clusterheads_rescore){
            ranked_poses.pop_back();
            cluster_assignments.pop_back();
        }
    }
    // assert( ranked_poses.size() == cluster_assignments.size() );
  
    ostringstream   text;    
    for(i=0;i<ranked_poses.size();i++){
       
        //collect data about primary score
        if (calc_rmsd) 
            text <<  calc_rmsd_string(ranked_poses[i].mol);

        if(cluster_ranked_poses){
            text << DELIMITER << setw(STRING_WIDTH) << "Cluster_Size:" 
                 << setw(FLOAT_WIDTH) << fixed << cluster_assignments[i] << endl;
        }
        text << ranked_poses[i].mol.current_data << endl;

        ranked_poses[i].mol.primary_data = text.str();

        text.str("");            

        //minimize or score using secondary function
        if(min.secondary_min_pose){
            min.secondary_minimize_pose(ranked_poses[i].mol, score);
        } else {
            score.compute_secondary_score(ranked_poses[i].mol);
        }
        ranked_poses[i].score = ranked_poses[i].mol.current_score;
        
    }

    //sort by secondary scores
    sort(ranked_poses.begin(), ranked_poses.end());

    trace.integer( "num_secondary_scored_poses", num_secondary_scored_poses );
    trace.integer( "ranked_poses.size", ranked_poses.size() );
    while(ranked_poses.size() > num_secondary_scored_poses){
        ranked_poses.pop_back();
    }

    //loop until user parameter reached
    for(i=0;i<ranked_poses.size();i++){
        secondary_ligand_out_confs << "\n" << DELIMITER << setw(STRING_WIDTH) << "Name:" 
            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.title << "_" << i + 1 << endl;
        
        secondary_ligand_out_confs << DELIMITER << setw(STRING_WIDTH) << "SECONDARY SCORE" << endl;
        
        if (calc_rmsd)
            secondary_ligand_out_confs << calc_rmsd_string(ranked_poses[i].mol);

        secondary_ligand_out_confs << ranked_poses[i].mol.current_data;
        
        secondary_ligand_out_confs << DELIMITER << setw(STRING_WIDTH) << "PRIMARY SCORE" << endl;
        secondary_ligand_out_confs << ranked_poses[i].mol.primary_data << endl;

        write_mol(ranked_poses[i].mol, secondary_ligand_out_confs);
    }


    if(ranked_poses[0].mol.num_atoms > 0) {
        DOCKMol tmp_dockmol;
        copy_molecule(tmp_dockmol, ranked_poses[0].mol);
        poses_for_rescore.push_back(tmp_dockmol);
    }

    
}
/************************************************/
void
Library_File::submit_scored_pose(DOCKMol & mol, Master_Score & score,
                                 Simplex_Minimizer & min)
{
    char            score_convert[100];
    RANKMol         tmp_rankmol;
    DOCKMol         tmp_dockmol;
    int             i;
    float           max_score;
    int             max_position;
    bool            success;

    num_confs++;

    if (score.use_score) {
        if(!min.minimize_ligand){
                score.compute_primary_score(mol);
        }
    } else {
        mol.current_score = 0.0;
        mol.internal_energy = 0.0; // trent 2009-02-13
        mol.current_data = "";
    }

    // Multi-pose output
    if (ranked_poses.size() < num_scored_poses) {

        // Add the molecule
        ranked_poses.push_back(tmp_rankmol);
        ranked_poses[ranked_poses.size() - 1].score = mol.current_score;
        ranked_poses[ranked_poses.size() - 1].data = mol.current_data;
        copy_molecule(ranked_poses[ranked_poses.size() - 1].mol, mol);
        
        // find position with the highest score
        max_score = ranked_poses[0].score;
        max_position = 0;

        for (i = 0; i < ranked_poses.size(); i++) {
            if (ranked_poses[i].score > max_score)
                max_position = i;
        }

        high_score = max_score;
        high_position = max_position;


    } else if (mol.current_score < high_score) {

        // add mol to the high position
        ranked_poses[high_position].data = mol.current_data;
        ranked_poses[high_position].score = mol.current_score;
        copy_molecule(ranked_poses[high_position].mol, mol);
        
        // find position with the highest score
        max_score = ranked_poses[0].score;
        max_position = 0;

        for (i = 0; i < ranked_poses.size(); i++) {
            if (ranked_poses[i].score > max_score)
                max_position = i;
        }

        high_score = max_score;
        high_position = max_position;


    }


}

/************************************************/
void
Library_File::write_scored_poses(bool USE_FILT, bool USE_MPI, Master_Score & score)
{
    Trace trace( "Library_File::write_scored_poses" );
    int             i,
                    j;
    // int insert_point;
    ofstream        ofs;
    SCOREMol        tmp_mol;
    DOCKMol         tmp_dockmol;
    //PAK

    /**
    change if statement- only skip the writing if it is a client
    for both single proc. and master mode, writing should occur
    **/
 
    if (ranked_poses.size() > 0) {

        if (USE_MPI) {
            cout << "\n" << "-----------------------------------" << "\n";
            cout << "Molecule: " << ranked_poses[0].mol.title << "\n" << "\n";
        }
        // print out the ligand info
        cout << " Anchors:\t\t" << num_anchors << "\n";
        cout << " Orientations:\t\t" << num_orients << "\n";
        cout << " Conformations:\t\t" << num_confs << "\n";

        cout << endl;
          
        if (ranked_poses[0].mol.num_atoms > 0)
            if(use_secondary_score)
                cout << " Primary Score" << "\n";
            for (i = 0; i < ranked_poses[0].data.size(); i++){
                if (ranked_poses[0].data[i] != '#'){
                    cout << ranked_poses[0].data[i];
                    cout.flush();
                }
            }

        if (cluster_ranked_poses) {

            for (i = 0; i < ranked_poses.size(); i++) {
                //PAK
                calc_num_HBA_HBD(ranked_poses[i].mol);

                // GDRM, October 2020
                // If possible, the lines below should not be here. Variables should be 
                // calculated elsewhere and stored in DOCKMol. Check copying function in
                // dockmol.cpp

                // Calculate xlogp
                //xlogp  filterXlogP {};
                //filterXlogP.getXLogP(ranked_poses[i].mol);

                if (ranked_poses[i].mol.num_atoms > 0) {
                    // if clustering is used, and molecule is a
                    // clusterhead, write it out
                    if (cluster_assignments[i] > 0) {
                        ligand_out_scored << "\n" << DELIMITER << setw(STRING_WIDTH) << "Name:" 
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.title << endl;
                        ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "DOCK_Rotatable_Bonds:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.rot_bonds << endl;
                        if (ranked_poses[i].mol.heavy_atoms > 0)
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "Heavy_Atoms:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.heavy_atoms << endl;
                        ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "Molecular_Weight:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.mol_wt << endl;
                        ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "DOCK_Rotatable_Bonds:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.rot_bonds << endl;
                        ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "Formal_Charge:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.formal_charge << endl;
                        ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "HBond_Acceptors:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.hb_acceptors << endl;
                        ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "HBond_Donors:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.hb_donors << endl;
                        ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "Heavy_Atoms:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.heavy_atoms << endl;
                        #ifdef BUILD_DOCK_WITH_RDKIT
                        //if dbfilter is used...write theses scores
                        if (USE_FILT) {
                            std::vector<std::string> vecpnstmp{};
                            vecpnstmp = ranked_poses[i].mol.pns_name;
                            std::string molpns_name = std::accumulate(vecpnstmp.begin(), vecpnstmp.end(), std::string(""));
                            if (ranked_poses[i].mol.pns == 0 && molpns_name.empty() == true){
                                molpns_name = "NO_PAINS"; 
                            }else if (ranked_poses[i].mol.pns != 0 && molpns_name.empty() == false) {}
                            else{
                                molpns_name = "ERROR_IN_PAINS_MATCHING_PLEASE_INVESTIGATE";
                            }


                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_num_arom_rings:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.num_arom_rings << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_num_alip_rings:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.num_alip_rings << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_num_sat_rings:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.num_sat_rings << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_Stereocenters:" 
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.num_stereocenters << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_Spiro_atoms:" 
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.num_spiro_atoms << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_LogP:" 
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.clogp << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_TPSA:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.tpsa << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_SYNTHA:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.sa_score << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_QED:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.qed_score << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_LogS:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.esol << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_num_of_PAINS:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.pns << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_PAINS_names:"
                                << setw(FLOAT_WIDTH) << molpns_name << endl;   
                            }
                        #endif
                        //ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "xLogP:" 
                        //    << setw(FLOAT_WIDTH) << ranked_poses[i].mol.xlogp << endl;
                        
                        #ifdef BUILD_DOCK_WITH_RDKIT 
                        if (USE_FILT) {
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_SMILES:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.smiles << endl;
                        }
                        #endif
                     
                        if (calc_rmsd)
                            ligand_out_scored << calc_rmsd_string(ranked_poses[i].mol);

                        ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "Cluster_Size:" 
                            << setw(FLOAT_WIDTH) << fixed << cluster_assignments[i] << "\n";
                        ligand_out_scored << ranked_poses[i].data << "\n";
                        write_mol(ranked_poses[i].mol, ligand_out_scored);
                        if (write_footprints)
                           write_footprint_file(ranked_poses[i].mol, score, ligand_out_footprint_scored);
                        if (write_hbonds)
                           write_hbond_file(ranked_poses[i].mol, score, ligand_out_hbond_scored);
                    } 
                }
            } 
        } else {

            for (i = 0; i < ranked_poses.size(); i++) {
                 calc_num_HBA_HBD(ranked_poses[i].mol); 
              
                 // GDRM, October 2020
                 // If possible, the lines below should not be here. Variables should be
                 // calculated elsewhere and stored in DOCKMol. Check copying function in
                 // dockmol.cpp

                 //// Calculate xlogp
                 //xlogp  filterXlogP {};
                 //filterXlogP.getXLogP(ranked_poses[i].mol);
               
                if (ranked_poses[i].mol.num_atoms > 0) {
                    #ifdef BUILD_DOCK_WITH_RDKIT 
                    if ( USE_FILT && (ranked_poses[i].mol.fails_filt || ranked_poses[i].mol.bad_molecule) ) {
                        ligand_out_rejected << "\n" << DELIMITER << setw(STRING_WIDTH) << "Name:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.title << endl;
                        ligand_out_rejected << DELIMITER << setw(STRING_WIDTH) << "Molecular_Weight:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.mol_wt << endl;
                        ligand_out_rejected << DELIMITER << setw(STRING_WIDTH) << "DOCK_Rotatable_Bonds:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.rot_bonds << endl;
                        ligand_out_rejected << DELIMITER << setw(STRING_WIDTH) << "Formal_Charge:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.formal_charge << endl;
                        ligand_out_rejected << DELIMITER << setw(STRING_WIDTH) << "HBond_Acceptors:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.hb_acceptors << endl;
                        ligand_out_rejected << DELIMITER << setw(STRING_WIDTH) << "HBond_Donors:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.hb_donors << endl;
                        ligand_out_rejected << DELIMITER << setw(STRING_WIDTH) << "Heavy_Atoms:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.heavy_atoms << endl;
                        if (ranked_poses[i].mol.bad_molecule && ranked_poses[i].mol.fails_filt){
                            ligand_out_rejected << DELIMITER << setw(STRING_WIDTH) << "Cause_of_rej:"
                                << setw(FLOAT_WIDTH) << "bad_input_filt_fail" << endl;
                        }else if (ranked_poses[i].mol.bad_molecule){
                            ligand_out_rejected << DELIMITER << setw(STRING_WIDTH) << "Cause_of_rej:"
                                << setw(FLOAT_WIDTH) << "bad_input" << endl; 
                        }else if (ranked_poses[i].mol.fails_filt){
                            ligand_out_rejected << DELIMITER << setw(STRING_WIDTH) << "Cause_of_rej:"
                                << setw(FLOAT_WIDTH) << "filter_fails" << endl;
                        }
                        ligand_out_rejected << ranked_poses[i].data << "\n";
                        write_mol(ranked_poses[i].mol, ligand_out_rejected); 
                    }else{
                    #endif
                        ligand_out_scored << "\n" << DELIMITER << setw(STRING_WIDTH) << "Name:" 
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.title << "\n";
                        ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "Molecular_Weight:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.mol_wt << endl;
                        ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "DOCK_Rotatable_Bonds:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.rot_bonds << endl;
                        ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "Formal_Charge:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.formal_charge << endl;
                        ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "HBond_Acceptors:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.hb_acceptors << endl;
                        ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "HBond_Donors:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.hb_donors << endl;
                        ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "Heavy_Atoms:"
                            << setw(FLOAT_WIDTH) << ranked_poses[i].mol.heavy_atoms << endl;
                        #ifdef BUILD_DOCK_WITH_RDKIT
                        //if dbfilter is used...write theses scores
                        if (USE_FILT) {
                            std::vector<std::string> vecpnstmp{};
                            vecpnstmp = ranked_poses[i].mol.pns_name;
                            std::string molpns_name = std::accumulate(vecpnstmp.begin(), vecpnstmp.end(), std::string(""));
                            if (ranked_poses[i].mol.pns == 0 && molpns_name.empty() == true){
                                molpns_name = "NO_PAINS"; 
                            }else if(ranked_poses[i].mol.pns != 0 && molpns_name.empty() == false) {}
                            else{
                                molpns_name = "ERROR_IN_PAINS_MATCHING_PLEASE_INVESTIGATE";
                            }

                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_num_arom_rings:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.num_arom_rings << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_num_alip_rings:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.num_alip_rings << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_num_sat_rings:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.num_sat_rings << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_Stereocenters:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.num_stereocenters << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_Spiro_atoms:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.num_spiro_atoms << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_LogP:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.clogp << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_TPSA:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.tpsa << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_SYNTHA:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.sa_score << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_QED:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.qed_score << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_LogS:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.esol << endl;
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_num_of_PAINS:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.pns << endl; 
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_PAINS_names:"
                               << setw(FLOAT_WIDTH) << molpns_name << endl;                   
                        }
                        #endif
                        //ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "xLogP:"
                        //    << setw(FLOAT_WIDTH) << ranked_poses[i].mol.xlogp << endl;
                        #ifdef BUILD_DOCK_WITH_RDKIT
                        if (USE_FILT) {
                            ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "RD_SMILES:"
                                << setw(FLOAT_WIDTH) << ranked_poses[i].mol.smiles << endl;
                        }
                       
                        #endif

                        if (calc_rmsd)
                            ligand_out_scored << calc_rmsd_string(ranked_poses[i].mol); 
                        ligand_out_scored << ranked_poses[i].mol.current_data << "\n";
                        write_mol(ranked_poses[i].mol, ligand_out_scored);
                        if (write_footprints)
                           write_footprint_file(ranked_poses[i].mol, score, ligand_out_footprint_scored);
                        if (write_hbonds)
                           write_hbond_file(ranked_poses[i].mol, score, ligand_out_hbond_scored);
                    #ifdef BUILD_DOCK_WITH_RDKIT
                    }
                    #endif
                }

           }
        }

        if (rank_ligands) {

            if (ranked_poses[0].mol.num_atoms > 0) {

                tmp_mol.first = ranked_poses[0].score;
                copy_molecule(tmp_mol.second, ranked_poses[0].mol);
                tmp_mol.second.score_text_data = ranked_poses[0].data;

                ranked_list.push_back(tmp_mol);

                if (ranked_list.size() > 2 * max_ranked) {
                    sort(ranked_list.begin(), ranked_list.end(),
                         less_than_pair);

                    while (ranked_list.size() > max_ranked)
                        ranked_list.pop_back();
                }

            }
        }
        if(num_scored_poses < 2){
            if(use_secondary_score){
                if(ranked_poses[0].mol.num_atoms > 0) {
                    copy_molecule(tmp_dockmol, ranked_poses[0].mol);
                    poses_for_rescore.push_back(tmp_dockmol);
                }
            }
        }

    }

}

/************************************************/
void
Library_File::write_ranked_ligands(bool USE_FILT , Master_Score & score )
{
    int             i;
    sort(ranked_list.begin(), ranked_list.end(), less_than_pair);
    while (ranked_list.size() > max_ranked)
        ranked_list.pop_back();

    if (rank_ligands) {
        for (i = 0; i < ranked_list.size(); i++) {

            calc_num_HBA_HBD(ranked_list[i].second);
            // GDRM, October 2020
            // If possible, the lines below should not be here. Variables should be
            // calculated elsewhere and stored in DOCKMol. Check copying function in
            // dockmol.cpp

            //// Calculate xlogp
            //xlogp  filterXlogP {};
            ligand_out_ranked << "\n"<< DELIMITER << setw(STRING_WIDTH) << "Name:" 
                << setw(FLOAT_WIDTH) << ranked_list[i].second.title << endl;
            ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "DOCK_Rotatable_Bonds:"
                << setw(FLOAT_WIDTH) << ranked_list[i].second.rot_bonds << endl;
            if (ranked_list[i].second.heavy_atoms > 0)
                ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "Heavy_Atoms:"
                    << setw(FLOAT_WIDTH) << ranked_list[i].second.heavy_atoms << endl;
            ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "Molecular_Weight:"
                << setw(FLOAT_WIDTH) << ranked_list[i].second.mol_wt << endl;
            ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "DOCK_Rotatable_Bonds:"
                << setw(FLOAT_WIDTH) << ranked_list[i].second.rot_bonds << endl;
            ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "Formal_Charge:"
                << setw(FLOAT_WIDTH) << ranked_list[i].second.formal_charge << endl;
            ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "HBond_Acceptors:"
                << setw(FLOAT_WIDTH) << ranked_list[i].second.hb_acceptors << endl;
            ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "HBond_Donors:"
                << setw(FLOAT_WIDTH) << ranked_list[i].second.hb_donors << endl;
            ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "Heavy_Atoms:"
                << setw(FLOAT_WIDTH) << ranked_list[i].second.heavy_atoms << endl;
            #ifdef BUILD_DOCK_WITH_RDKIT
            //if dbfilter is used...write theses scores
            if (USE_FILT) {
                std::vector<std::string> vecpnstmp{};
                vecpnstmp = ranked_list[i].second.pns_name;
                std::string molpns_name = std::accumulate(vecpnstmp.begin(), vecpnstmp.end(), std::string(""));
                if (ranked_list[i].second.pns == 0 && molpns_name.empty() == true){
                    molpns_name = "NO_PAINS"; 
                }else if ( ranked_list[i].second.pns != 0 && molpns_name.empty() == false) {
               
                }else {
                    molpns_name = "ERROR_IN_PAINS_MATCHING_PLEASE_INVESTIGATE";
                }

                ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "RD_num_arom_rings:"
                    << setw(FLOAT_WIDTH) << ranked_list[i].second.num_arom_rings << endl;
                ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "RD_num_alip_rings:"
                    << setw(FLOAT_WIDTH) << ranked_list[i].second.num_alip_rings << endl;
                ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "RD_num_sat_rings:"
                    << setw(FLOAT_WIDTH) << ranked_list[i].second.num_sat_rings << endl;
                ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "RD_Stereocenters:"
                    << setw(FLOAT_WIDTH) << ranked_list[i].second.num_stereocenters << endl;
                ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "RD_Spiro_atoms:"
                    << setw(FLOAT_WIDTH) << ranked_list[i].second.num_spiro_atoms << endl;
                ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "RD_LogP:"
                    << setw(FLOAT_WIDTH) << ranked_list[i].second.clogp << endl;
                ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "RD_TPSA:"
                    << setw(FLOAT_WIDTH) << ranked_list[i].second.tpsa << endl;
                ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "RD_SYNTHA:"
                    << setw(FLOAT_WIDTH) << ranked_list[i].second.sa_score << endl;
                ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "RD_QED:"
                    << setw(FLOAT_WIDTH) << ranked_list[i].second.qed_score << endl;
                ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "RD_LogS:"
                    << setw(FLOAT_WIDTH) << ranked_list[i].second.esol << endl;
                ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "RD_num_of_PAINS:"
                    << setw(FLOAT_WIDTH) << ranked_list[i].second.pns << endl; 
                ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "RD_PAINS_names:"
                    << setw(FLOAT_WIDTH) << molpns_name << endl;
            }
            #endif
            //ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "xLogP:"
            //    << setw(FLOAT_WIDTH) << ranked_list[i].second.xlogp << endl;
            #ifdef BUILD_DOCK_WITH_RDKIT
            if (USE_FILT) {
                ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "RD_SMILES:"
                    << setw(FLOAT_WIDTH) << ranked_list[i].second.smiles << endl;
            }
            #endif


            if (calc_rmsd)
                ligand_out_ranked << calc_rmsd_string(ranked_list[i].second);

            ligand_out_ranked << ranked_list[i].second.score_text_data << endl;

            write_mol(ranked_list[i].second, ligand_out_ranked);
            if (write_footprints){
               /////////////////////////////////////////////////////////////////
               // RANKMol and SCOREMol are reduenent
               // Dockmol has a score_text_data, it is the same as data in RANKMol.
               // so RANKMol has the same information stored in two places
               /////////////////////////////////////////////////////////////////
               //  write footprint file
               write_footprint_file(ranked_list[i].second, score, ligand_out_footprint_ranked);
            }
            if (write_hbonds)
               write_hbond_file(ranked_list[i].second, score, ligand_out_hbond_ranked);
        }
    }

}

/************************************************/
void
Library_File::submit_rmsd_reference(DOCKMol & mol, AMBER_TYPER & typer)
{
    if (constant_rmsd_ref) {
        ifstream rmsd_ref(constant_rmsd_ref_file.c_str());
        if (!rmsd_ref) {
            cout << "\n\nError: Something is wrong with the RMSD reference."<< endl
                 <<     "       Could not open file " << constant_rmsd_ref_file
                 << " for reading.\nProgram will terminate.\n" << endl;
            exit(0);
        }

        Read_Mol2(rmsd_reference, rmsd_ref, false, false, false);

        rmsd_ref.close();
        // prepare_molecule(DOCKMol & mol, bool read_vdw, bool use_chem, bool use_ph4, bool use_volume)
        // the prepare_molecule function will assign heavy atom flag needed by the
        // min RMSD function.
        // typer.prepare_molecule(rmsd_reference,true,false,false);
        // Bugfix for DOCK6.7

    } else {
        // the rmsd reference molecule is the input molecule, mol
        copy_molecule(rmsd_reference, mol);
    }

    // Initialize RMSD reference in all cases:
    typer.prepare_molecule(rmsd_reference,true,false,false,false);
}


/************************************************/
void
Library_File::calculate_rmsd(DOCKMol & refmol, DOCKMol & mol,double *rmsds)
{
    rmsds[0] = calculate_std_rmsd(refmol,mol); // upper bound
    Hungarian_RMSD a;
    rmsds[1] = a.calc_Hungarian_RMSD(refmol,mol);
    rmsds[2] = calculate_min_rmsd(refmol,mol); // lower bound
}


/************************************************/
void
Library_File::calculate_rmsd(DOCKMol & mol, double *rmsds)
{
    calculate_rmsd(rmsd_reference, mol, rmsds);
}


/************************************************/
// Calculates standard heavy atom rmsd, does not account for symmetry
float
Library_File::calculate_std_rmsd(DOCKMol & refmol, DOCKMol & mol)
{

    if (refmol.num_atoms != mol.num_atoms)
        return -1000.0;

    float rmsd = 0.0;
    int hcount = 0;
    int i;
    for (i = 0; i < mol.num_atoms; i++) {
        if (mol.amber_at_heavy_flag[i]) {
            rmsd +=
                (((mol.x[i] - refmol.x[i])*(mol.x[i] - refmol.x[i])) +
                 ((mol.y[i] - refmol.y[i])*(mol.y[i] - refmol.y[i])) +
                 ((mol.z[i] - refmol.z[i])*(mol.z[i] - refmol.z[i])));
            hcount++;
        }
    }

    // rmsd = rmsd / mol.num_atoms;
    rmsd = rmsd/hcount;
    rmsd = sqrt(rmsd);
    return rmsd;
}

/************************************************/
// One-way min rmsd function, should only be called by calculate_min_rmsd(mol,mol)
float
Library_File::calculate_min_rmsd_one_way(DOCKMol & refmol, DOCKMol & mol)
{
    if (refmol.num_atoms != mol.num_atoms) return -1000.0;
    float rmsd = 0.0;
    float dist = 0.0;
    float mindist = 100000;
    int hcount = 0;
    for (int i = 0; i < mol.num_atoms; i++) {
        mindist = 100000;
        if (mol.amber_at_heavy_flag[i]) {
            for (int j = 0; j < refmol.num_atoms; j++) {
                //if ((mol.amber_at_heavy_flag[j]) && (refmol.atom_types[i]==mol.atom_types[j])) {
                if (mol.atom_types[i]==refmol.atom_types[j]) {
                    dist =  (((refmol.x[j] - mol.x[i])*(refmol.x[j] - mol.x[i])) +
                             ((refmol.y[j] - mol.y[i])*(refmol.y[j] - mol.y[i])) +
                             ((refmol.z[j] - mol.z[i])*(refmol.z[j] - mol.z[i])));  // compute the distance between ai and bj.
                    if (dist < mindist) {
                        mindist = dist;     // Update minimum distance of ith ref to the pose.
                    }
                    
                }
            }    
            rmsd += mindist;
            hcount++;
        }
    }

    rmsd = rmsd/hcount;
    rmsd = sqrt(rmsd);
    return rmsd;

}

/************************************************/
// Minimum rmsd: same as autodock symmetry function
// Computes min distance between atoms of the same type; not 1-to-1; keeps max(a->b,b->a)
// should be called directly by wrapper function
float
Library_File::calculate_min_rmsd(DOCKMol & mol1, DOCKMol & mol2)
{
    if (mol1.num_atoms != mol2.num_atoms)
        return 0.0;
    // mol2 is always the mol, so heavy atoms are known
    // mol1 is always the ref, heavy atoms are NOT known
    // Set all the values in dockmol to be the same as mol2 except the crds 
    // in mol1_new; the crd in mol1_new are those from mol1
/*    DOCKMol mol1_new;
    copy_molecule(mol1_new, mol2);
    copy_crds(mol1_new,mol1);


    float   reftopose_rmsd = calculate_min_rmsd_one_way(mol1_new, mol2);
    float   posetoref_rmsd = calculate_min_rmsd_one_way(mol2, mol1_new);
*/
    
    float   reftopose_rmsd = calculate_min_rmsd_one_way(mol1, mol2);
    float   posetoref_rmsd = calculate_min_rmsd_one_way(mol2, mol1);

    //cout << "RMSD_AB: " << reftopose_rmsd << "; RMSD_BA: " << posetoref_rmsd << endl;

    return max(reftopose_rmsd,posetoref_rmsd);
}
/************************************************/
//Keep looking for the pairs of same type with minimum distance, and one atom from ref should be uniquely matched to one from pose.
//This algorithm is actually greedy algorithm, which behaves well in most case, but cannot ensure finding the exact minimum of square distance EVERY TIME.
float
Library_File::calculate_min_cor_rmsd(DOCKMol & refmol, DOCKMol & mol)
{
    if (refmol.num_atoms != mol.num_atoms) return 0.0;
    const int N = refmol.num_atoms;
    float rmsd = 0.0;
    float dist[N][N];     // matrix of square distance between atoms
    float mindist ;
    int hcount = 0;            
    bool row_bool[N];         // rows that has been used
    bool col_bool[N];         //columns that has been used
    int cor_col_index = -1;
    int cor_row_index = -1;
    for (int i = 0; i < N; i++)
        row_bool[i]= false;            //First set the boolean of all rows be false.
    for (int j = 0; j < N; j++)
        col_bool[j]= false;           // First set the boolean of all columns be false.
    
    for (int i = 0; i < N; i++) {
        if (refmol.amber_at_heavy_flag[i]) {
            for (int j = 0; j < N; j++) {
                if ((mol.amber_at_heavy_flag[j]) && (refmol.atom_types[i]==mol.atom_types[j])) {
                    dist[i][j] =  (((mol.x[j] - refmol.x[i])*(mol.x[j] - refmol.x[i])) +
                                   ((mol.y[j] - refmol.y[i])*(mol.y[j] - refmol.y[i])) +
                                   ((mol.z[j] - refmol.z[i])*(mol.z[j] - refmol.z[i])));  // compute the square distance between ith refmol atom  and jth mol.
                                    }

                }
            }
            hcount++;                    //heavy atom count noumber
        }

   for (int k = 0; k < N; k++){  // we find one pair in each loop
        mindist = 1000000.0;     // initial value of minimum distance is infinity.
        cor_col_index = -1;
        cor_row_index = -1;

        for (int i = 0; i < N; i++){
              if ( row_bool[i] || !refmol.amber_at_heavy_flag[i] )   // if the ith ref atom has been used, skip it. Also check if atom is heavy.
                 continue;
              for (int j = 0; j < N; j++){                                         
                  if (col_bool[j] ||( !(mol.amber_at_heavy_flag[j])) ||(refmol.atom_types[i]==mol.atom_types[j]) )
                        //if the jth atom has been used,or the type of atoms from ref and pose does not match, skip it                  
                     continue;
                  if (dist[i][j] < mindist){
                      mindist = dist[i][j];       // Update minimum distance of ith ref to the pose.
                     cor_row_index = i;   // index of ref atom that has the correspondence with min distance. Keep update it.
                     cor_col_index = j;   // index of pose atom that has the correspondence with min distance. Keep update it.
                     }
                }
        row_bool[cor_row_index] = true;      // change ith row boolean to be true, i.e never use the ith ref atom in future.
        col_bool[cor_col_index] = true;      // change ith row boolean to be true, i.e never use the ith pose atom in future.
        rmsd = rmsd + mindist;                // sum up the square distance
      }
      //  print "#### the ref atom " ,ref.atom_list[cor_row_index].name, " corresponds to the pose atom" ,pose.atom_list[cor_col_index].name," and the square distance is ",mindist

    rmsd = rmsd/hcount;
    rmsd = sqrt(rmsd);
    return rmsd;
   }
}

/************************************************/
void
Library_File::cluster_list(){

    // cluster_assignments[i] is the number of ranked_poses in the cluster
    // with head ranked_poses[i].  Thus,
    // assert( ranked_poses[i] corresponds to cluster_assignments[i] );

    // initialize clusters to every ranked pose is in its own cluster.
    cluster_assignments.clear();
    cluster_assignments.resize(ranked_poses.size(), 1);

    // then recluster via rmsd.
    float           crmsd;
    int i, j;
    // loop over all ranked poses
    for (i = 0; i < ranked_poses.size() - 1; i++) {

        // if i-pose is not assigned to a cluster already
        if (cluster_assignments[i] > 0) {

            // loop over all subsequent poses
            for (j = i + 1; j < ranked_poses.size(); j++) {

                // if j-pose is not assigned to a cluster already
                if (cluster_assignments[j] > 0) {
                    double rmsds[3]; 
                    calculate_rmsd(ranked_poses[i].mol, ranked_poses[j].mol,rmsds);
                    crmsd = rmsds[0];
                    // if j-pose is close enough to i-pose, add it to i-cluster
                    if (crmsd < cluster_rmsd_threshold) {
                        ++ cluster_assignments[i];
                        cluster_assignments[j] = 0;
                    }
                }
            }
        }
    }
    // assert( ranked_poses.size() == cluster_assignments.size() );
}

/************************************************/
void
Library_File::secondary_rescore_poses(Master_Score & score, Simplex_Minimizer & min)
{
    Trace trace( "Library_File::secondary_rescore_poses" );
    int i, j;
    ostringstream text;
    SCOREMol tmp_mol;
    
    for(i=0; i < poses_for_rescore.size(); ++i){

        if (num_scored_poses < 2){
            //collect data from primary scoring function
            if (calc_rmsd) 
                text << calc_rmsd_string(poses_for_rescore[i]);
            text << poses_for_rescore[i].current_data;

            poses_for_rescore[i].primary_data = text.str();
        }

        //collect ligands for preliminary ranking by primary score before rescoring
        if(rank_ligands || rank_secondary_ligands){
            
            tmp_mol.first = poses_for_rescore[i].current_score;
            copy_molecule(tmp_mol.second, poses_for_rescore[i]);

            ranked_secondary_list.push_back(tmp_mol);

        } else {
            //rescore if simply rescoring molecules
            if(num_scored_poses < 2){
                if(min.secondary_min_pose){
                    min.secondary_minimize_pose(poses_for_rescore[i], score);
                } else {
                    score.compute_secondary_score(poses_for_rescore[i]);
                }
            }
            //write score to output file
            //cout << endl << "-----------------------------------" << endl;
            cout << "Molecule: " << poses_for_rescore[i].title << endl << endl;
            cout << " Secondary Score"  << endl;
            for (j = 0; j < poses_for_rescore[i].current_data.size(); j++){
                if (poses_for_rescore[i].current_data[j] != '#')
                    cout << poses_for_rescore[i].current_data[j];
            }
        }
        text.str("");
    }

    // the above loop is the only way to populate ranked_secondary_list; srb.
    assert( ranked_secondary_list.size() <= poses_for_rescore.size() );

    //sort ligands by primary scoring function
    sort(ranked_secondary_list.begin(), ranked_secondary_list.end(), less_than_pair);
    //prune poses by user defined maximum for primary scoring
    trace.integer( "max_ranked", max_ranked );
    trace.integer( "ranked_secondary_list.size", ranked_secondary_list.size() );
    while(ranked_secondary_list.size() > max_ranked)
        ranked_secondary_list.pop_back();

    //rescore pruned list
    for(i=0; i < ranked_secondary_list.size(); ++i){
        if(num_scored_poses < 2){ 
            if(min.secondary_min_pose){
                min.secondary_minimize_pose(ranked_secondary_list[i].second, score);
            } else {
                score.compute_secondary_score(ranked_secondary_list[i].second);
            }
        }
        //write score to output file
        cout << "\n" "-----------------------------------" "\n"
                "Molecule: " << ranked_secondary_list[i].second.title << "\n\n"
                " Secondary Score" << endl;
        for (j = 0; j < ranked_secondary_list[i].second.current_data.size(); j++){
            if (ranked_secondary_list[i].second.current_data[j] != '#')
                cout << ranked_secondary_list[i].second.current_data[j];
        }
    }
}

/************************************************/
void
Library_File::submit_secondary_pose()
{
    int i;
    
    if(poses_for_rescore.size() > 0){
        
        if((!rank_ligands) && (!rank_secondary_ligands)){
            
            for(i=0;i<poses_for_rescore.size();i++){
            
                secondary_ligand_out_scored << "\n" << DELIMITER << setw(STRING_WIDTH) << "Name:"
                    << setw(FLOAT_WIDTH) << poses_for_rescore[i].title << endl;
        
                secondary_ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "SECONDARY SCORE" << endl;
                if (calc_rmsd) 
                    secondary_ligand_out_scored << calc_rmsd_string(poses_for_rescore[i]); 
                secondary_ligand_out_scored << poses_for_rescore[i].current_data;
            
                secondary_ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "PRIMARY SCORE" << endl;
                secondary_ligand_out_scored << poses_for_rescore[i].primary_data << endl;

                write_mol(poses_for_rescore[i], secondary_ligand_out_scored);
            }
        }

        if((rank_ligands) && (!rank_secondary_ligands)) {
        
            for(i=0;i<ranked_secondary_list.size();i++){
           
                secondary_ligand_out_scored << "\n" << DELIMITER << setw(STRING_WIDTH) << "Name:" 
                    << setw(FLOAT_WIDTH) << ranked_secondary_list[i].second.title << endl;

                secondary_ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "SECONDARY SCORE" << endl;
                if (calc_rmsd)
                    secondary_ligand_out_scored << calc_rmsd_string(ranked_secondary_list[i].second);
                secondary_ligand_out_scored << ranked_secondary_list[i].second.current_data; 

                secondary_ligand_out_scored << DELIMITER << setw(STRING_WIDTH) << "PRIMARY SCORE" << endl;
                secondary_ligand_out_scored << ranked_secondary_list[i].second.primary_data << endl;

                write_mol(ranked_secondary_list[i].second, secondary_ligand_out_scored);
           }
        }

    
        if(rank_secondary_ligands){

            for(i=0;i<ranked_secondary_list.size();i++)
                ranked_secondary_list[i].first = ranked_secondary_list[i].second.current_score;
            sort(ranked_secondary_list.begin(), ranked_secondary_list.end(), less_than_pair);
        
            //prune poses by user defined maximum for secondary scoring
            //which is    <= user defined maximum for primary   scoring.
            while(ranked_secondary_list.size() > max_secondary_ranked)
                ranked_secondary_list.pop_back();

            for(i=0;i<ranked_secondary_list.size();i++) {

                secondary_ligand_out_ranked << "\n" << DELIMITER << setw(STRING_WIDTH) << "Name:" 
                    << setw(FLOAT_WIDTH) << ranked_secondary_list[i].second.title << endl;

                secondary_ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "SECONDARY SCORE" << endl;
                if (calc_rmsd)
                    secondary_ligand_out_ranked << calc_rmsd_string(ranked_secondary_list[i].second);
                secondary_ligand_out_ranked << ranked_secondary_list[i].second.current_data;

                secondary_ligand_out_ranked << DELIMITER << setw(STRING_WIDTH) << "PRIMARY SCORE" << endl;
                secondary_ligand_out_ranked << ranked_secondary_list[i].second.primary_data << endl;

                write_mol(ranked_secondary_list[i].second, secondary_ligand_out_ranked);
            }

        }
    }

    close_files();
}


/************************************************
 Formats the RMSD string to be accurate to 4 decimal places (fixed point)
 instead of writing in the scientific notation. This prevents minor RMSD
 changes of 10^-7 A etc, which are likely a result of rounding, from 
 showing up in the output files.
*/

string Library_File::calc_rmsd_string(DOCKMol & mol)  {
    Trace trace( "Library_File::calc_rmsd_string" );
    ostringstream   text;
    text.flags(ios::fixed);
    double rmsds[3];
    calculate_rmsd(mol,rmsds);
    //cout << rmsds[0]<<","<<rmsds[1]<<","<<rmsds[2]<< endl;
    text << DELIMITER << setw(STRING_WIDTH) << "HA_RMSDs:" 
         << setw(FLOAT_WIDTH) << fixed << setprecision (3) << rmsds[0] << endl;
    text << DELIMITER << setw(STRING_WIDTH) << "HA_RMSDh:" 
         << setw(FLOAT_WIDTH) << fixed << setprecision (3) << rmsds[1] << endl;
    text << DELIMITER << setw(STRING_WIDTH) << "HA_RMSDm:" 
         << setw(FLOAT_WIDTH) << fixed << setprecision (3) << rmsds[2] << endl;
    return text.str();
}

// +++++++++++++++++++++++++++++++++++++++++
// LEP: Given a mol, calculate the hydrogen acceptors and donor
void
Library_File::calc_num_HBA_HBD( DOCKMol & mol )
{
    // Populate HD fields
    int counter = 0;
    for (int i=0; i<mol.num_atoms; i++){
        if (mol.flag_acceptor[i] == true){
           counter++;
        }
    }
    mol.hb_acceptors = counter;

    // Populate HA fields
    counter = 0;
    for (int i=0; i<mol.num_atoms; i++){
        if (mol.flag_donator[i] == true){
           counter++;
        }
    }
    mol.hb_donors = counter;

   return;

} //end LIbrary_File::calc_num_HBA_HBD();

