#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include "amber_typer.h"
#include "dockmol.h"
#include "fingerprint.h"

using namespace std;

/*
 * =================================================================== 
 */
char           *
white_line(char *line)
{
    for (unsigned int i = 0; i < strlen(line); i++)
        if (isspace(line[i]))
            line[i] = ' '; 

    return line;
}




/*
 * ==================================================================== 
 */
void
print_node(ATOM_TYPE_NODE & node, int level)
{
    int             i;

    if (level) {
        if (node.include)
            cout << " (";
        else
            cout << " [";
    }

    if (node.multiplicity)
        cout << node.multiplicity << " ";

    cout << node.type;

    for (i = 0; i < node.next_total; i++)
        print_node(node.next[i], level + 1);

    if (level) {
        if (node.include)
            cout << ")";
        else
            cout << "]";
    }
}

/*
 * =================================================================== 
 */
int
assign_node(ATOM_TYPE_NODE & node, int include)
{
    char            temp[6];
    char           *branch;

    node.next_total = 0;
    strcpy(temp, strtok(NULL, " "));

    if (isdigit(temp[0])) {
        node.multiplicity = atoi(temp);

        if (node.multiplicity < 0) {
            cout << "Cannot specify negative multiplicity in definition.\n";
            return 0;
        }

        strcpy(node.type, strtok(NULL, " "));

    } else {
        node.multiplicity = 0;
        strcpy(node.type, temp);
    }

    node.include = include;

    while (branch = strtok(NULL, " ")) {

        if ((!strncmp(branch, "(", 1)) || (!strncmp(branch, "[", 1))) {

            if (node.next_total >= 6) {
                cout <<
                    "Cannot exceed 6 substituents for each atom in definition.\n";
                return 0;
            }

            ATOM_TYPE_NODE  tmp_node;
            node.next.push_back(tmp_node);

            if ((!strncmp(branch, "(", 1))
                && (!assign_node(node.next[node.next_total], 1)))
                return 0;

            if ((!strncmp(branch, "[", 1))
                && (!assign_node(node.next[node.next_total], 0)))
                return 0;

            node.next_total++;
        }

        else if ((!strncmp(branch, ")", 1)) || (!strncmp(branch, "]", 1)))
            return 1;

        else
            return 0;
    }

    return 1;
}

/*
 * =================================================================== 
 */
int
check_type(const char *candidate, const char *reference)
{
    if ((strstr(candidate, reference)) || (reference[0] == '*'))
        return 1;
    else
        return 0;
}

/*
 * =================================================================== 
 */
bool
count_bond_neighbors(DOCKMol & mol, int bond_num)
{
    int             c1,
                    c2;

    c1 = mol.get_atom_neighbors(mol.bonds_origin_atom[bond_num]).size();
    c2 = mol.get_atom_neighbors(mol.bonds_target_atom[bond_num]).size();


    //*****//
    // WJA - Add this part in if you want halogens chopped off during fragment library gen
    //string          oa_type,
    //                ta_type;
    //
    //oa_type = mol.atom_types[mol.bonds_origin_atom[bond_num]];
    //ta_type = mol.atom_types[mol.bonds_target_atom[bond_num]];
    //
    //if (oa_type == "F" || oa_type == "Cl" || oa_type == "Br" || oa_type == "I")
    //    c1++;
    //
    //if (ta_type == "F" || ta_type == "Cl" || ta_type == "Br" || ta_type == "I")
    //    c2++;
    //*****//


    if ((c1 <= 1) || (c2 <= 1))
        return false;
    else
        return true;
}

/*
 * =================================================================== 
 */
int
count_atom_neighbors(DOCKMol & mol, int atom_num)
{
    int             count;

    count = mol.get_atom_neighbors(atom_num).size();

    return count;
}

/*
 * =================================================================== 
 */
int
check_bonded_atoms(DOCKMol & mol, int current_atom, int previous_atom,
                   ATOM_TYPE_NODE & node)
{
    int             i,
                    j,
                    next_atom;
    int             match,
                    match_count;
    vector < int   >nbrs;

    match_count = 0;

    // loop over neighbors of current atom
    nbrs = mol.get_atom_neighbors(current_atom);
    for (i = 0; i < nbrs.size(); i++) {

        next_atom = nbrs[i];
        match = 0;

        if ((next_atom != previous_atom)
            && (match =
                check_type(mol.atom_types[next_atom].c_str(), node.type))) {

            for (j = 0; j < node.next_total; j++) {
                if (!check_bonded_atoms
                    (mol, next_atom, current_atom, node.next[j]))
                    match = 0;
            }

        }

        if (match)
            match_count++;
    }

    if (node.multiplicity) {
        if (node.multiplicity == match_count)
            match = 1;
        else
            match = 0;
    } else {
        if (match_count)
            match = 1;
        else
            match = 0;
    }

    if (match == node.include)
        return 1;
    else
        return 0;

}

/*
 * =================================================================== 
 */
int
check_atom(DOCKMol & mol, int current_atom, ATOM_TYPE_NODE & node)
{
    int             match = 0;
    int             i;

    if (match = check_type(mol.atom_types[current_atom].c_str(), node.type)) {

        for (i = 0; i < node.next_total; i++) {

            if (!check_bonded_atoms
                (mol, current_atom, current_atom, node.next[i]))
                match = 0;

        }

    }

    return match;
}
/***********************************/
ATOM_TYPE::ATOM_TYPE()
{

    name[0] = '\0';
    atom_model = '\0';
    radius = MIN_FLOAT;
    well_depth = MIN_FLOAT;
    heavy_flag = INT_MIN;
    valence = INT_MIN;
    bump_id = INT_MIN;
    gbradius = MIN_FLOAT;
    gbscale = MIN_FLOAT;
    definitions.reserve(0);

}
/***********************************/
ATOM_TYPE::~ATOM_TYPE()
{

}
/***********************************/
void
ATOM_TYPER::get_vdw_labels(string fname, bool read_gb_parm)
{
    FILE           *ifp;
    char            line[100],
                    model[100];
    ATOM_TYPE       tmp_type;
    ATOM_TYPE_NODE  tmp_node;
    int             i;

    ifp = fopen(fname.c_str(), "r");

    if (ifp == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    while (fgets(line, 100, ifp) != NULL) {

        if (!strncmp(line, "name", 4)) {        // read in name field
            types.push_back(tmp_type);

            if (sscanf(line, "%*s %s", types[types.size() - 1].name) < 1) {
                cout << "Incomplete vdw member declaration.\n";
                exit(0);
            }

        } else if (!strncmp(line, "atom_model", 10)) {  // read in atom model
                                                        // field
            if (sscanf(line, "%*s %s", model) != 1) {
                cout << "Incomplete atom_model specification.\n";
                exit(0);
            }

            types[types.size() - 1].atom_model = tolower(model[0]);

            if ((types[types.size() - 1].atom_model != 'a')
                && (types[types.size() - 1].atom_model != 'u')
                && (types[types.size() - 1].atom_model != 'e')) {
                cout <<
                    "Atom_model specification restricted to ALL, UNITED, or EITHER.\n";
                exit(0);
            }
        } else if (!strncmp(line, "heavy_flag", 8)) {   // read in heavy flag
                                                        // field
            if (sscanf(line, "%*s %d", &types[types.size() - 1].heavy_flag) !=
                1) {
                cout << "Incomplete heavy_flag specification.\n";
                exit(0);
            }
        } else if (!strncmp(line, "radius", 6)) {       // read in radius field
            if (sscanf(line, "%*s %f", &types[types.size() - 1].radius) != 1) {
                cout << "Incomplete radius specification.\n";
                exit(0);
            }
        } else if (!strncmp(line, "well_depth", 10)) {  // read in well depth
                                                        // field
            if (sscanf(line, "%*s %f", &types[types.size() - 1].well_depth) !=
                1) {
                cout << "Incomplete well_depth specification.\n";
                exit(0);
            }
        } else if (!strncmp(line, "valence", 7)) {      // read in valence
                                                        // field
            if (sscanf(line, "%*s %d", &types[types.size() - 1].valence) != 1) {
                cout << "Incomplete valence specification.\n";
                exit(0);
            }
        } else if (!strncmp(line, "gbradii", 7)) {
            if (sscanf(line, "%*s %f", &types[types.size() - 1].gbradius) != 1) {
                cout << "Incomplete GBRadius specification.\n";
                exit(0);
            }
        } else if (!strncmp(line, "gbscale", 7)) {
            if (sscanf(line, "%*s %f", &types[types.size() - 1].gbscale) != 1) {
                cout << "Incomplete GBScale specification.\n";
                exit(0);
            }
        } else if (!strncmp(line, "definition", 10)) {  // read in definition
                                                        // fields
            strtok(white_line(line), " ");
            types[types.size() - 1].definitions.push_back(tmp_node);

            if (!assign_node
                (types[types.size() - 1].
                 definitions[types[types.size() - 1].definitions.size() - 1],
                 true)) {
                cout << "Error assigning vdw member definitions.\n";
                exit(0);
            }
        }

    }                           // end While

    // check to see if any vdw values have been read in
    if (types.size() < 1) {
        cout << "ERROR:  VDW parameter file empty." << endl;
        exit(0);
    }
    // check to see if all values for vdw parameters have been read in
    for (i = 0; i < types.size() - 1; i++) {
        if (types[i].atom_model == '\0') {
            cout << "ERROR:  No atom_model assigned for " << types[i].
                name << " in VDW parameter file." << endl;
            exit(0);
        }
        if (types[i].radius == MIN_FLOAT) {
            cout << "ERROR:  No radius assigned for " << types[i].
                name << " in VDW parameter file." << endl;
            exit(0);
        }
        if (types[i].well_depth == MIN_FLOAT) {
            cout << "ERROR:  No well_depth assigned for " << types[i].
                name << " in VDW parameter file." << endl;
            exit(0);
        }
        if (types[i].heavy_flag == INT_MIN) {
            cout << "ERROR:  No heavy_flag assigned for " << types[i].
                name << " in VDW parameter file." << endl;
            exit(0);
        }
        if (types[i].valence == INT_MIN) {
            cout << "ERROR:  No valence assigned for " << types[i].
                name << " in VDW parameter file." << endl;
            exit(0);
        }
        if (read_gb_parm) {
            if (types[i].gbradius == MIN_FLOAT) {
                cout << "ERROR:  No gbradii assigned for " << types[i].
                    name << " in VDW parameter file." << endl;
                exit(0);
            }

            if (types[i].gbscale == MIN_FLOAT) {
                cout << "ERROR:  No gbscale assigned for " << types[i].
                    name << " in VDW parameter file." << endl;
                exit(0);
            }
        }
        if (types[i].definitions.size() < 1) {
            cout << "ERROR:  No definitions assigned for " << types[i].
                name << " in VDW parameter file." << endl;
            exit(0);
        }

    }

    for (i = 0; i < types.size() - 2; i++) {    // calculate bump_id values
        if (types[i].heavy_flag)
            types[i].bump_id = NINT(10.0 * types[i].radius);
        else
            types[i].bump_id = 0;
    }

    fclose(ifp);
}

/***********************************/
int
ATOM_TYPER::assign_vdw_labels(DOCKMol & mol, int atom_model)
{
    int             i,
                    j,
                    k;
    int             vdw_assigned = false;
    vector < int   >nbrs;

    atom_types.clear();

    // loop over all atoms in the mol
    for (i = 0; i < mol.num_atoms; i++) {

        // by default, assign dummy type
        atom_types.push_back(0);

        // loop over types read in from vdw defn file
        for (j = 0; j < types.size(); j++) {

            if ((atom_model == 'a') && (types[j].atom_model == 'u'))
                continue;

            if ((atom_model == 'u') && (types[j].atom_model == 'a'))
                continue;

            for (k = 0; k < types[j].definitions.size(); k++) {

                if (check_atom(mol, i, types[j].definitions[k])) {
                    atom_types[i] = j;
                    vdw_assigned = true;
                }               // end k loop

            }

        }                       // end j loop

        // check for vdw label assignment
        if (!vdw_assigned) {
            cout << "WARNING assign_vdw_labels: No vdw parameters for ";
            cout << i << " " << " " << mol.atom_types[i] << endl;
            // return false;
        }
        // check that no atom's valence is violated
        if (count_atom_neighbors(mol, i) > types[atom_types[i]].valence) {

            // loop over neighbor atoms
            nbrs = mol.get_atom_neighbors(i);
            k = 0;

            for (j = 0; j < nbrs.size(); j++)
                if ((strcmp(mol.atom_types[nbrs[j]].c_str(), "LP"))
                    && (strcmp(mol.atom_types[nbrs[j]].c_str(), "Du")))
                    k++;

            if (k > types[atom_types[i]].valence) {
                cout << "WARNING assign_vdw_labels: Atom valence violated for ";
                cout << mol.title << " ";
                cout << "atom number: " << i << endl;

                // return false;
            }

        }
        // transfer partial charges for united models
        if ((atom_model == 'u')
            && (fabs(types[atom_types[i]].well_depth) < 0.00001)) {
            if (count_atom_neighbors(mol, i) == 1) {
                nbrs = mol.get_atom_neighbors(i);
                mol.charges[nbrs[0]] = mol.charges[nbrs[0]] + mol.charges[i];
                mol.charges[i] = 0.0;
            } else {
                cout <<
                    "WARNING assign_vdw_labels: Unable to transfer partial charge away from ";
                cout << mol.title << " ";
                cout << "atom number: " << i << endl;

                mol.charges[i] = 0.0;
            }
        }

    }                           // end i loop

    return true;

}

/***********************************/
void
BOND_TYPER::get_flex_labels(string fname)
{
    FILE           *ifp;
    char            line[100];
    BOND_TYPE       tmp_flex;
    int             i;
    int             definition_count;

    ifp = fopen(fname.c_str(), "r");

    if (ifp == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    while (fgets(line, 100, ifp)) {

        if (!strncmp(line, "name", 4)) {
            definition_count = 0;

            types.push_back(tmp_flex);

            if (types.size() > 1) {
            }

            if (sscanf(line, "%*s %s", types[types.size() - 1].name) < 1) {
                cout << "Incomplete Flex Definition Failure.\n";
                exit(0);
            }

            for (i = 0; i < strlen(types[types.size() - 1].name); i++)
                types[types.size() - 1].name[i] =
                    (char) tolower(types[types.size() - 1].name[i]);

            types[types.size() - 1].drive_id = -1;
            types[types.size() - 1].minimize = -1;

        }                       // End "name" if

        else if (!strncmp(line, "drive_id", 6))
            sscanf(line, "%*s %d", &types[types.size() - 1].drive_id);

        else if (!strncmp(line, "minimize", 8))
            sscanf(line, "%*s %d", &types[types.size() - 1].minimize);

        else if (!strncmp(line, "definition", 10)) {

            strtok(white_line(line), " ");
            assign_node(types[types.size() - 1].definition[definition_count],
                        1);
            definition_count++;

        }                       // End "definition" if
    }

    fclose(ifp);
}

/***********************************/
void
BOND_TYPER::get_flex_search(string fname)
{
    FILE           *ifp;
    char            line[100];
    char           *token;
    int             i;

    int             drive_id;
    int             torsion_total;
    FLOATVec        torsions;

    ifp = fopen(fname.c_str(), "r");

    if (ifp == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    for (i = 0; i < types.size(); i++)
        types[i].torsions.clear();

    while (fgets(line, 100, ifp)) {

        token = strtok(white_line(line), " ");

        if (!strcmp(token, "drive_id")) {

            torsions.clear();

            if (token = strtok(NULL, " "))
                drive_id = atoi(token);

            else {
                cout <<
                    "ERROR get_flex_search: Search_id value not specified in ";
                cout << fname << endl;
                exit(0);
            }

            if (!fgets(line, 100, ifp)
                || !(token = strtok(white_line(line), " "))
                || strcmp(token, "positions")) {
                cout <<
                    "ERROR get_flex_search: Positions field doesn't follow Id in ";
                cout << fname << endl;
                exit(0);
            }

            if (token = strtok(NULL, " "))
                torsion_total = atoi(token);
            else {
                cout <<
                    "ERROR get_flex_search: Postions value not specified in ";
                cout << fname << endl;
                exit(0);
            }

            if (!fgets(line, 100, ifp)
                || !(token = strtok(white_line(line), " "))
                || strcmp(token, "torsions")) {
                cout <<
                    "ERROR get_flex_search: Torsions doesn't follow Positions in ";
                cout << fname << endl;
                exit(0);
            }

            for (i = 0; i < torsion_total; i++) {
                if (token = strtok(NULL, " ")) {
                    torsions.push_back(atof(token));
                } else {
                    cout <<
                        "ERROR get_flex_search: Insufficient number of torsions in ";
                    cout << fname << endl;
                    exit(0);
                }

            }

            for (i = 0; i < types.size(); i++) {
                if (types[i].drive_id == drive_id) {
                    types[i].torsion_total = torsion_total;
                    types[i].torsions = torsions;
                }
            }

        }                       // End if drive_id

    }                           // End while


    for (i = 0; i < types.size(); i++) {
        if (types[i].torsion_total < 1) {
            cout << "ERROR get_flex_search: Missing torsion parameters in ";
            cout << fname << endl;
            exit(0);
        }
    }

    fclose(ifp);


}

/***********************************/
void
BOND_TYPER::apply_flex_labels(DOCKMol & mol)
{
    int             i,
                    j;

    flex_ids.clear();
    total_torsions = 0;

    for (i = 0; i < mol.num_bonds; i++) {

        flex_ids.push_back(0);
        flex_ids[i] = -1;

        if (mol.bond_ring_flags[i])
            continue;

        if (!count_bond_neighbors(mol, i))
            continue;

        for (j = 0; j < types.size(); j++) {    // Loop over bond type (1 to
                                                // size)

            if (check_atom
                (mol, mol.bonds_origin_atom[i], types[j].definition[0])
                && check_atom(mol, mol.bonds_target_atom[i],
                              types[j].definition[1])) {
                flex_ids[i] = j;
            } else
                if (check_atom
                    (mol, mol.bonds_origin_atom[i], types[j].definition[1])
                    && check_atom(mol, mol.bonds_target_atom[i],
                                  types[j].definition[0])) {
                flex_ids[i] = j;
            }

        }                       // End loop over bond types

        if (flex_ids[i] != -1)
            total_torsions++;
    }

}

/***********************************/
bool
BOND_TYPER::is_rotor(int bond)
{

    if (flex_ids[bond] > -1)
        return true;
    else
        return false;
}

/***********************************/
void
CHEM_TYPER::get_chem_labels(string fname)
{
    FILE           *ifp;
    char            line[100];
    CHEM_TYPE       tmp_type;
    ATOM_TYPE_NODE  tmp_node;

    ifp = fopen(fname.c_str(), "r");

    if (ifp == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    while (fgets(line, 100, ifp) != NULL) {

        if (!strncmp(line, "name", 4)) {        // read in name field
            types.push_back(tmp_type);

            if (sscanf(line, "%*s %s", types[types.size() - 1].name) < 1) {
                //cout << "Incomplete vdw member declaration.\n";
                cout << "Incomplete chem type member declaration.\n";//LINGLING, typo fixed
                exit(0);
            }
        } else if (!strncmp(line, "definition", 10)) {  // read in definition
                                                        // fields
            strtok(white_line(line), " ");
            types[types.size() - 1].definitions.push_back(tmp_node);

            if (!assign_node
                (types[types.size() - 1].
                 definitions[types[types.size() - 1].definitions.size() - 1],
                 true)) {
                //cout << "Error assigning vdw member definitions.\n";
                cout << "Error assigning chem type member definition.\n";//LINGLING, typo fixed
                exit(0);
            }
        }
    }
    fclose(ifp);//LINGLING
}

/***********************************/
void
CHEM_TYPER::apply_chem_labels(DOCKMol & mol)
{
    int             i,
                    j,
                    k;

    chem_type_ids.clear();

    // loop over all atoms in the mol
    for (i = 0; i < mol.num_atoms; i++) {

        // by default, assign dummy type
        chem_type_ids.push_back(-1);

        // loop over types read in from chem defn file
        for (j = 0; j < types.size(); j++) {
            for (k = 0; k < types[j].definitions.size(); k++) {
                if (check_atom(mol, i, types[j].definitions[k])) {
                    chem_type_ids[i] = j;
                }
            }                   // end k loop
        }                       // end j loop
    }                           // end i loop

}

/***********************************/
void
PH4_TYPER::get_ph4_labels(string fname)
{
    FILE           *ifp;
    char            line[100];
    PH4_TYPE       tmp_type;
    ATOM_TYPE_NODE  tmp_node;

    ifp = fopen(fname.c_str(), "r");

    if (ifp == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    while (fgets(line, 100, ifp) != NULL) {

        if (!strncmp(line, "name", 4)) {        // read in name field
            types.push_back(tmp_type);

            if (sscanf(line, "%*s %s", types[types.size() - 1].name) < 1) {
                cout << "Incomplete ph4 type member declaration.\n";
                exit(0);
            }
        } else if (!strncmp(line, "definition", 10)) {  // read in definition
                                                        // fields
            strtok(white_line(line), " ");
            types[types.size() - 1].definitions.push_back(tmp_node);

            if (!assign_node
                (types[types.size() - 1].
                 definitions[types[types.size() - 1].definitions.size() - 1],
                 true)) {
                cout << "Error assigning ph4 type member definition.\n";
                exit(0);
            }
        }

   }
   fclose(ifp);
}

/***********************************/
void
PH4_TYPER::apply_ph4_labels(DOCKMol & mol)
{
    int             i,
                    j,
                    k;

    ph4_type_ids.clear();

    // loop over all atoms in the mol
    for (i = 0; i < mol.num_atoms; i++) {

        // by default, assign dummy type
        ph4_type_ids.push_back(-1);

        // loop over types read in from ph4 defn file
        for (j = 0; j < types.size(); j++) {
            for (k = 0; k < types[j].definitions.size(); k++) {
                if (check_atom(mol, i, types[j].definitions[k])) {
                    ph4_type_ids[i] = j;
                }
            }                   // end k loop
        }                       // end j loop
    }                           // end i loop

}

/***********************************/
void
AMBER_TYPER::input_parameters(Parameter_Reader & parm, bool read_vdw,
                              bool use_chem, bool use_ph4, bool use_volume)
{

    if (read_vdw || use_chem || use_ph4 || use_volume) {
        cout << "\nAtom Typing Parameters" << endl;
        cout <<
            "------------------------------------------------------------------------------------------"
            << endl;

        if (read_vdw) {
            atom_model =
                parm.query_param("atom_model", "all", "all | united").c_str()[0];
            vdw_defn_file = parm.query_param("vdw_defn_file", "vdw.defn");
            flex_defn_file = parm.query_param("flex_defn_file", "flex.defn");
            flex_drive_tbl =
                parm.query_param("flex_drive_file", "flex_drive.tbl");
        }
        if (use_chem || use_volume) {
            chem_defn_file = parm.query_param("chem_defn_file", "chem.defn");
        }
        if (use_ph4) {
            ph4_defn_file = parm.query_param("pharmacophore_defn_file", "ph4.defn");
        }
    }
    verbose = parm.verbosity_level();
    
}

/***********************************/
void
AMBER_TYPER::initialize(bool read_vdw, bool read_gb_parm, bool use_chem, bool use_ph4, bool use_volume)
{

    if (read_vdw) {
        atom_typer.get_vdw_labels(vdw_defn_file.c_str(), read_gb_parm);
        bond_typer.get_flex_labels(flex_defn_file.c_str());
        bond_typer.get_flex_search(flex_drive_tbl.c_str());
    }
    if (use_chem || use_volume) {
        chem_typer.get_chem_labels(chem_defn_file.c_str());
    }
    if (use_ph4)  {
        ph4_typer.get_ph4_labels(ph4_defn_file.c_str());
    }
    skip_verbose_flag = false;
}

/***********************************/
// H-bond add
void
AMBER_TYPER::assign_hbond_labels( DOCKMol & mol )
{
 // assign hbond labels
 // Courtney's hydrogen bond fix
 for (int i=0;i<mol.num_atoms;i++){
      stringstream ss;
      stringstream ss1;
      char atom_type_char[20];
      ss << mol.atom_types[i]; 
      ss >> atom_type_char;
      // see if the atom is a acceptor
      if (    atom_type_char[0] == 'O' 
          ||  atom_type_char[0] == 'N' 
          ||  atom_type_char[0] == 'S'
          ||  atom_type_char[0] == 'F'
          || (atom_type_char[0] == 'C' && atom_type_char[1] == 'l')
          || (atom_type_char[0] == 'C' && atom_type_char[1] == 'L'))
      {
         mol.flag_acceptor[i] = true;
      }
      // see if atom is a donor
      // atoms must be an h and must be conected to a acceptor through a bond
      if(   atom_type_char[0] == 'H')
      {
            for (int j=0;j<mol.num_bonds;j++)
            {
                 if( i == mol.bonds_origin_atom[j])// atom i (which is an h) the start a bond
                 {
                    ss1 << mol.atom_types[mol.bonds_target_atom[j]]; 
                    ss1 >> atom_type_char;
                    if (    atom_type_char[0] == 'O'
                        ||  atom_type_char[0] == 'N'
                        ||  atom_type_char[0] == 'S'
                        ||  atom_type_char[0] == 'F'
                        || (atom_type_char[0] == 'C' && atom_type_char[1] == 'l')
                        || (atom_type_char[0] == 'C' && atom_type_char[1] == 'L'))
                    {
                           mol.flag_donator[i] = true;
                           mol.acc_heavy_atomid[i] = mol.bonds_target_atom[j];
                    }
                 }

                 if( i == mol.bonds_target_atom[j])// atom i (which is an h) the terminates a bond
                 {
                    ss1 << mol.atom_types[mol.bonds_origin_atom[j]]; 
                    ss1 >> atom_type_char;
                    if (    atom_type_char[0] == 'O'
                        ||  atom_type_char[0] == 'N'
                        ||  atom_type_char[0] == 'S'
                        ||  atom_type_char[0] == 'F'
                        || (atom_type_char[0] == 'C' && atom_type_char[1] == 'l')
                        || (atom_type_char[0] == 'C' && atom_type_char[1] == 'L'))
                    {
                           mol.flag_donator[i] = true;
                           mol.acc_heavy_atomid[i] = mol.bonds_origin_atom[j];
                    }
                 }

            }
      }
 }

}


/***********************************/
void
AMBER_TYPER::prepare_molecule(DOCKMol & mol, bool read_vdw, bool use_chem, bool use_ph4, bool use_volume)
{
    int             i;

//cout <<"Entering amber_typer::prepare_molecule " <<endl;
//cout <<"num of atoms in dockmol = " <<mol.num_atoms <<endl;
//cout <<"read_vdw = " <<read_vdw <<";  use_chem = " <<use_chem <<endl;

    // assign hbond acc and donor labels
    assign_hbond_labels(mol);

    // use the amber atom typers to assign atom and bond types
    if (read_vdw) {
        atom_typer.assign_vdw_labels(mol, atom_model);
        bond_typer.apply_flex_labels(mol);
    }
    if (use_chem || use_volume)
        chem_typer.apply_chem_labels(mol);
    if (use_ph4)
        ph4_typer.apply_ph4_labels(mol);


    // copy atom types to molecule
    if (read_vdw) {
        for (i = 0; i < mol.num_atoms; i++) {
            mol.amber_at_id[i] = atom_typer.atom_types[i];
            mol.amber_at_radius[i] =
                atom_typer.types[atom_typer.atom_types[i]].radius;
            mol.amber_at_well_depth[i] =
                atom_typer.types[atom_typer.atom_types[i]].well_depth;
            mol.amber_at_heavy_flag[i] =
                atom_typer.types[atom_typer.atom_types[i]].heavy_flag;
            mol.amber_at_valence[i] =
                atom_typer.types[atom_typer.atom_types[i]].valence;
            mol.amber_at_bump_id[i] =
                atom_typer.types[atom_typer.atom_types[i]].bump_id;
            mol.gb_hawkins_radius[i] =
                atom_typer.types[atom_typer.atom_types[i]].gbradius;
            mol.gb_hawkins_scale[i] =
                atom_typer.types[atom_typer.atom_types[i]].gbscale;
        }
        mol.amber_at_assigned = true;

	if (!skip_verbose_flag){
            if ( 0 != verbose ){
            	int nheavy_atoms=0;
            	for (i = 0; i < mol.num_atoms; i++) {
                    if (mol.amber_at_heavy_flag[i])
                    	nheavy_atoms++;
            	}
            	cout << endl <<  "-----------------------------------" << endl;
            	cout << "VERBOSE MOLECULE STATS" << endl << endl;
            	cout << "  Number of heavy atoms = " << nheavy_atoms << endl;
            }
	}

        // copy bond types to molecule
        for (i = 0; i < mol.num_bonds; i++) {
            mol.amber_bt_id[i] = bond_typer.flex_ids[i];

            if (bond_typer.flex_ids[i] != -1) {

                mol.amber_bt_minimize[i] =
                    bond_typer.types[bond_typer.flex_ids[i]].minimize;
                mol.amber_bt_torsion_total[i] =
                    bond_typer.types[bond_typer.flex_ids[i]].torsion_total;
                mol.amber_bt_torsions[i] =
                    bond_typer.types[bond_typer.flex_ids[i]].torsions;

            } else {
                mol.amber_bt_minimize[i] = 0;
                mol.amber_bt_torsion_total[i] = 0;
                mol.amber_bt_torsions[i].clear();
            }

        }
        mol.amber_bt_assigned = true;

        mol.rot_bonds=0;  //defined in dockmol.cpp as part of dbfilter code
        for (i = 0; i < mol.num_bonds; i++)
            if (mol.bond_is_rotor(i)) mol.rot_bonds++;

        mol.formal_charge = 0.0;  //defined in dockmol.cpp as part of dbfilter code
        for (i = 0; i < mol.num_atoms; i++)
            mol.formal_charge += mol.charges[i];

        mol.heavy_atoms = 0;
        // Count the number of heavy atoms
        for (i = 0; i < mol.num_atoms; i++)
            if (mol.amber_at_heavy_flag[i]) mol.heavy_atoms++;

        mol.mol_wt = getMW(mol);  //defined in dockmol.cpp as part of dbfilter code

        //The following codes that generate atom_keys for tanimoto calculation is adopted
        //from compute_tanimoto function in fingerprint.cpp.  It is moved here to be 
        //calculated when preparing molecules to speed up tanimoto calculation.  YZ

        Fingerprint tmp_finger;

        vector <int> tmp_atom_vec;
        vector < pair< pair<int,int>, string > > tmp_bond_vec;

        tmp_finger.prepare_noH_vectors(mol, tmp_atom_vec, tmp_bond_vec);

        mol.atom_keys_0 = tmp_finger.generate_keys( mol, 0, tmp_atom_vec, tmp_bond_vec );
        mol.atom_keys_1 = tmp_finger.generate_keys( mol, 1, tmp_atom_vec, tmp_bond_vec );
        mol.atom_keys_2 = tmp_finger.generate_keys( mol, 2, tmp_atom_vec, tmp_bond_vec );

        tmp_atom_vec.clear();
        tmp_bond_vec.clear();


	if (!skip_verbose_flag){
            if ( 0 != verbose ) {
                cout << "  Number of rotatable bonds = " << mol.rot_bonds << endl;

                ostringstream charge_text;
                //LEP set precision
                charge_text << "  Formal Charge = " << fixed << setprecision(3)
                        << mol.formal_charge << endl
                        << "  Molecular Weight = " << fixed << setprecision(3) 
		        <<  mol.mol_wt << endl
                        << "  Heavy Atoms = " << mol.heavy_atoms << endl;

                // sudipto: Print warning for non-integral formal charge here

                cout << charge_text.str();
            }
	}
    }

    if (use_chem || use_volume) {
        for (i = 0; i < mol.num_atoms; i++) {
            mol.chem_types[i] =
                chem_typer.types[chem_typer.chem_type_ids[i]].name;
        }
        mol.chem_types_assigned = true;
    }
    if (use_ph4) {
        for (i = 0; i < mol.num_atoms; i++) {
            mol.ph4_types[i] =
                ph4_typer.types[ph4_typer.ph4_type_ids[i]].name;
        }
        mol.ph4_types_assigned = true;    
    }    

}

/***********************************/
float AMBER_TYPER::getMW(DOCKMol & mol)
{
     //atomic weights from General Chemistry 3rd Edition Darrell D. Ebbing
     float mw = 0.0;
     char sybyl_type[10];

     for (int i = 0; i < mol.num_atoms; i++) {

        strcpy(sybyl_type, mol.atom_types[i].c_str());
        string element = strtok(sybyl_type,".");

        if (element == "O") mw += 15.9994;
        else if (element == "N") mw += 14.00674;
        else if (element == "C") mw += 12.011;
        else if (element == "F") mw += 18.9984032;
        else if (element == "Cl") mw += 35.4527;
        else if (element == "Br") mw += 79.904;
        else if (element == "I") mw += 126.90447;
        else if (element == "H" ) mw += 1.00794;
        else if (element == "B")  mw += 10.811;
        else if (element == "S" ) mw += 32.066;
        else if (element == "P")  mw += 30.973762;
        else if (element == "Li") mw += 6.941;
        else if (element == "Na") mw += 22.98968;
        else if (element == "Mg") mw += 24.3050;
        else if (element == "Al") mw += 26.981539;
        else if (element == "Si") mw += 28.0855;
        else if (element == "K") mw += 39.0983;
        else if (element == "Ca") mw += 40.078;
        else if (element == "Cr") mw += 51.9961;
        else if (element == "Mn") mw += 54.93805;
        else if (element == "Fe") mw += 55.847;
        else if (element == "Co") mw += 58.93320;
        else if (element == "Cu") mw += 63.546;
        else if (element == "Zn") mw += 65.39;
        else if (element == "Se") mw += 78.96;
        else if (element == "Mo") mw += 95.94;
        else if (element == "Sn") mw += 118.710;
        else if (element == "LP")  mw += 0.0;
        else if (element == "Du") mw += 0.0;
        else cout << "Element " << element << " not found in MW code\n";
    }
    return mw;
}

//#ifdef BUILD_DOCK_WITH_RDKIT
//// RDTYPER methods
//
//void RDTYPER::initialize_parameters(Parameter_Reader &parm){
//}
//
//double RDTYPER::calculate_qed_score(std::string smi, std::vector<double> w){
//}
//
//double RDTYPER::calculate_sa_score(std::string smi){
//}
//
//double RDTYPER::calculate_esol(std::string smi, bool delaney){
//}
//
//void RDTYPER::calculate_descriptors(DOCKMol &mol, bool ){
//}
//
//void RDTYPER::assign_descriptors(DOCKMol &mol){
//}
//
//#endif


