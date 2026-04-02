#include "fingerprint.h"
using namespace std;

class DOCKMol;
class Parameter_Reader;


// +++++++++++++++++++++++++++++++++++++++++
// Constructors and destructors
Fingerprint::Fingerprint()
{
}
Fingerprint::~Fingerprint()
{
}
Node::Node()
{
    atom_index = -1;
    parent_index = -1;
    atom_type = "\0";
    bond_to_parent = 0;
    children_identities = "\0";
    children_bond_types = "\0";
    grandchildren_num = "\0";
    grandchildren_identities = "\0";
    grandchildren_bond_types = "\0";
}
Node::~Node()
{
}



// +++++++++++++++++++++++++++++++++++++++++
// Get parameters from dock.in file 
void
Fingerprint::input_parameters( Parameter_Reader & parm )
{
    // The only time a user would want to set 'compute_fingerprint = yes' is if they are performing a
    // virtual screen and they only want to dock molecules that are similiar to a reference molecule.


    //cout << "\nFingerprint Parameters" << endl;
    //cout <<
    //    "------------------------------------------------------------------------------------------"
    //    << endl;


    // Compute fingerprint?
    dbfilter_compute_fingerprint = (parm.query_param("dbfilter_compute_fingerprint", "no", "yes no") == "yes") ? true : false;

    if (dbfilter_compute_fingerprint){
        dbfilter_tanimoto_reference = parm.query_param("dbfilter_tanimoto_reference", "reference.mol2");

        ifstream fin;
        fin.open(dbfilter_tanimoto_reference.c_str());
        if (fin.fail()){
            cout <<"Warning: Could not open " <<dbfilter_tanimoto_reference <<endl;
        } else {
            Read_Mol2(ref_mol, fin, false, false, false);
        }

        dbfilter_tanimoto_cutoff = atof(parm.query_param("dbfilter_tanimoto_cutoff", "0.5").c_str());
    }


    return;

} // end Fingerprint::input_parameters()



// +++++++++++++++++++++++++++++++++++++++++
// Given a dockmol object, this function will compute 3-atom deep keys for all heavy atoms,
// canonicalize each of the trees, turn them into strings, alphabatize, and concatenate. The
// string it returns can be used to check for uniqueness. 
string
Fingerprint::compute_unique_string( DOCKMol & mol )
{

    // Declare some temporary vectors to describe the mol without Hydrogens
    vector <int> tmp_atom_vec;
    vector < pair< pair<int,int>, string > > tmp_bond_vec;

    // Fill the vectors with atom / bond information for all atoms except hydrogens
    prepare_noH_vectors(mol, tmp_atom_vec, tmp_bond_vec);

    // Generate sets of atom keys for each molecule to different depths
    vector <string> atom_keys_3 = generate_keys( mol, 3, tmp_atom_vec, tmp_bond_vec );

    // Print the keys to the screen
    // cout <<"\nKeys that are 3 bonds deep:\n";
    // for (int i=0; i<atom_keys_3.size(); i++){ cout <<atom_keys_3[i] <<"\n"; }


    // Alphabatize and concatenate keys
    sort(atom_keys_3.begin(), atom_keys_3.end());
    stringstream alpha;
    for (int i=0; i<atom_keys_3.size(); i++){ alpha <<atom_keys_3[i]; }
    string alpha_string = alpha.str();

    // Print the unique string for this molecule
    // cout <<"\nAlphabatized and concatenated:\n";
    // cout <<alpha.str() <<"\n";


    // Clear the atom and bond vecs for this molecule
    alpha.clear();
    atom_keys_3.clear();
    tmp_atom_vec.clear();
    tmp_bond_vec.clear();


    return alpha_string;

} // end Fingerprint::compute_unique_string()



// +++++++++++++++++++++++++++++++++++++++++
// Given a dockmol object, this function will compute 3-atom deep keys for all heavy atoms,
// canonicalize each of the trees, turn them into strings, alphabatize, and concatenate. The
// string it returns can be used to check for uniqueness.
string
Fingerprint::compute_unique_string_active( DOCKMol & mol )
{

    // Declare some temporary vectors to describe the mol without Hydrogens
    vector <int> tmp_atom_vec;
    vector < pair< pair<int,int>, string > > tmp_bond_vec;

    // Fill the vectors with atom / bond information for all atoms except hydrogens OR un-active atoms
    prepare_noH_vectors_for_active_atoms(mol, tmp_atom_vec, tmp_bond_vec);

    // Generate sets of atom keys for each molecule to different depths
    vector <string> atom_keys_3 = generate_keys( mol, 3, tmp_atom_vec, tmp_bond_vec );

    // Print the keys to the screen
    // cout <<"\nKeys that are 3 bonds deep:\n";
    // for (int i=0; i<atom_keys_3.size(); i++){ cout <<atom_keys_3[i] <<"\n"; }


    // Alphabatize and concatenate keys
    sort(atom_keys_3.begin(), atom_keys_3.end());
    stringstream alpha;
    for (int i=0; i<atom_keys_3.size(); i++){ alpha <<atom_keys_3[i]; }
    string alpha_string = alpha.str();

    // Print the unique string for this molecule
    // cout <<"\nAlphabatized and concatenated:\n";
    // cout <<alpha.str() <<"\n";


    // Clear the atom and bond vecs for this molecule
    alpha.clear();
    atom_keys_3.clear();
    tmp_atom_vec.clear();
    tmp_bond_vec.clear();


    return alpha_string;

} // end Fingerprint::compute_unique_string_active()



// +++++++++++++++++++++++++++++++++++++++++
// Given two DOCKMols, this function computes all of the 'atom key' bits 0-, 1-, and 2-bonds 
// deep for all heavy atoms. Then it computes and returns the Tanimoto coefficient between
// two DOCKMols using all atom keys as bits.
float
Fingerprint::compute_tanimoto( DOCKMol & mol1, DOCKMol & mol2 )
{
/*
    // Declare some temporary vectors to describe each mol without Hydrogens
    vector <int> tmp_atom_vec1;
    vector <int> tmp_atom_vec2;
    vector < pair< pair<int,int>, string > > tmp_bond_vec1;
    vector < pair< pair<int,int>, string > > tmp_bond_vec2;

    // Fill the vectors with atom / bond information for all atoms except hydrogens
    prepare_noH_vectors(mol1, tmp_atom_vec1, tmp_bond_vec1);
    prepare_noH_vectors(mol2, tmp_atom_vec2, tmp_bond_vec2);


    // Generate sets of atom keys for each molecule to different depths
    // mol1
    vector <string> atom_keys1_0 = generate_keys( mol1, 0, tmp_atom_vec1, tmp_bond_vec1 );
    vector <string> atom_keys1_1 = generate_keys( mol1, 1, tmp_atom_vec1, tmp_bond_vec1 );
    vector <string> atom_keys1_2 = generate_keys( mol1, 2, tmp_atom_vec1, tmp_bond_vec1 );
    // mol2
    vector <string> atom_keys2_0 = generate_keys( mol2, 0, tmp_atom_vec2, tmp_bond_vec2 );
    vector <string> atom_keys2_1 = generate_keys( mol2, 1, tmp_atom_vec2, tmp_bond_vec2 );
    vector <string> atom_keys2_2 = generate_keys( mol2, 2, tmp_atom_vec2, tmp_bond_vec2 );

    // Clear some of the vector data from memory
    tmp_atom_vec1.clear();
    tmp_atom_vec2.clear();
    tmp_bond_vec1.clear();
    tmp_bond_vec2.clear();
*/

    //The atom_keys of the molecules are now computed when amber_typer is preparing the 
    //molecule, and stored in the DOCKMol object to speed up tanimoto calculation.  The 
    //above commented out codes are moved to amber_typer.cpp.  YZ

    // Some variables for computing the tanimoto coefficient
    int tc_both = 0;
    int tc_mol1 = 0;
    int tc_mol2 = 0;
    float tc_final = 0.0;

    // These weigths are optional - only leave one pair uncommented. (These should probably be 
    // input parameters)
    float weight_alpha = 1.0;  float weight_beta = 1.0;    // Pure Tanimoto coefficient
    //float weight_alpha = 0.5;  float weight_beta = 0.5;    // Dice's coefficient
    //float weight_alpha = 0.75; float weight_beta = 0.25;   // Tversky index, weighted prototype
    //float weight_alpha = 0.25; float weight_beta = 0.75;   // Tversky index, weighted variant


    running_tanimoto( tc_both, tc_mol1, tc_mol2, mol1.atom_keys_0, mol2.atom_keys_0 );
    running_tanimoto( tc_both, tc_mol1, tc_mol2, mol1.atom_keys_1, mol2.atom_keys_1 );
    running_tanimoto( tc_both, tc_mol1, tc_mol2, mol1.atom_keys_2, mol2.atom_keys_2 );


    // Compute and return tanimoto coefficient
    tc_final = ( (float)tc_both / ((float)tc_both + weight_alpha*(float)tc_mol1 + weight_beta*(float)tc_mol2) );

    return tc_final;

} // end Fingerprint::compute_tanimoto()


//The compute_tanimoto_fragment function is the exact same function as the old 
//compute_tanimoto function, now it is only used for computing tanimoto between
//fragments to generate the fragment graph.  YZ
float
Fingerprint::compute_tanimoto_fragment( DOCKMol & mol1, DOCKMol & mol2 )
{

    // Declare some temporary vectors to describe each mol without Hydrogens
    vector <int> tmp_atom_vec1;
    vector <int> tmp_atom_vec2;
    vector < pair< pair<int,int>, string > > tmp_bond_vec1;
    vector < pair< pair<int,int>, string > > tmp_bond_vec2;

    // Fill the vectors with atom / bond information for all atoms except hydrogens
    prepare_noH_vectors(mol1, tmp_atom_vec1, tmp_bond_vec1);
    prepare_noH_vectors(mol2, tmp_atom_vec2, tmp_bond_vec2);


    // Generate sets of atom keys for each molecule to different depths
    // mol1
    vector <string> atom_keys1_0 = generate_keys( mol1, 0, tmp_atom_vec1, tmp_bond_vec1 );
    vector <string> atom_keys1_1 = generate_keys( mol1, 1, tmp_atom_vec1, tmp_bond_vec1 );
    vector <string> atom_keys1_2 = generate_keys( mol1, 2, tmp_atom_vec1, tmp_bond_vec1 );
    // mol2
    vector <string> atom_keys2_0 = generate_keys( mol2, 0, tmp_atom_vec2, tmp_bond_vec2 );
    vector <string> atom_keys2_1 = generate_keys( mol2, 1, tmp_atom_vec2, tmp_bond_vec2 );
    vector <string> atom_keys2_2 = generate_keys( mol2, 2, tmp_atom_vec2, tmp_bond_vec2 );

    // Clear some of the vector data from memory
    tmp_atom_vec1.clear();
    tmp_atom_vec2.clear();
    tmp_bond_vec1.clear();
    tmp_bond_vec2.clear();


    // Some variables for computing the tanimoto coefficient
    int tc_both = 0;
    int tc_mol1 = 0;
    int tc_mol2 = 0;
    float tc_final = 0.0;

    // These weigths are optional - only leave one pair uncommented. (These should probably be
    // input parameters)
    float weight_alpha = 1.0;  float weight_beta = 1.0;    // Pure Tanimoto coefficient

    running_tanimoto( tc_both, tc_mol1, tc_mol2, atom_keys1_0, atom_keys2_0 );
    running_tanimoto( tc_both, tc_mol1, tc_mol2, atom_keys1_1, atom_keys2_1 );
    running_tanimoto( tc_both, tc_mol1, tc_mol2, atom_keys1_2, atom_keys2_2 );


    // Compute and return tanimoto coefficient
    tc_final = ( (float)tc_both / ((float)tc_both + weight_alpha*(float)tc_mol1 + weight_beta*(float)tc_mol2) );

    return tc_final;

} // end Fingerprint::compute_tanimoto()

// +++++++++++++++++++++++++++++++++++++++++
// In de novo growth, torsion environments are used to determine whether the connection
// of two fragments is 'allowed'. Pass this function a DOCKMol with only two active atoms -
// the 1-bond-deep environments will be passed back via a referenced pair of strings, which
// is also passed to this function.
void
Fingerprint::return_torsion_environments( DOCKMol & mol, pair <string, string> & tmp_pair )
{
    // Declare some temporary vectors to describe the mol without Hydrogens
    vector <int> tmp_atom_vec;
    vector < pair< pair<int,int>, string > > tmp_bond_vec;

    // Fill the vectors with atom / bond information for all atoms except hydrogens
    prepare_noH_vectors_for_active_atoms(mol, tmp_atom_vec, tmp_bond_vec);

    // Generate sets of atom keys for each molecule to different depths
    vector <string> atom_keys_1 = generate_keys( mol, 1, tmp_atom_vec, tmp_bond_vec );

    // Take keys two at a time, sort, concatenate, and print

        sort(atom_keys_1.begin(), atom_keys_1.end());
        tmp_pair.first = atom_keys_1[0];
        tmp_pair.second = atom_keys_1[1];


        //stringstream alpha;
        //alpha << atom_keys_1[0];
        //alpha << "-";
        //alpha << atom_keys_1[1];
        //string alpha_string = alpha.str();
        //alpha.clear();

    // Print the unique string for this molecule
    // cout <<"\nAlphabatized and concatenated:\n";
    //cout <<alpha.str() <<"\n";


    // Clear the atom and bond vecs for this molecule
    // alpha.clear();
    atom_keys_1.clear();
    tmp_atom_vec.clear();
    tmp_bond_vec.clear();

    return;

} // end Fingerprint::return_torsion_environments();



// +++++++++++++++++++++++++++++++++++++++++
//
string
Fingerprint::return_environment( DOCKMol & mol, int atom_num )
{
    // Declare some temporary vectors to describe the mol without Hydrogens
    vector <int> tmp_atom_vec;
    vector < pair< pair<int,int>, string > > tmp_bond_vec;

    // Fill the vectors with atom / bond information for all atoms except hydrogens
    prepare_vectors(mol, tmp_atom_vec, tmp_bond_vec);

    // Generate sets of atom keys for each molecule to different depths
    vector <string> atom_keys_1 = generate_keys( mol, 1, tmp_atom_vec, tmp_bond_vec );

    string return_string = atom_keys_1[atom_num];

    // Clear the atom and bond vecs for this molecule
    // alpha.clear();
    atom_keys_1.clear();
    tmp_atom_vec.clear();
    tmp_bond_vec.clear();

    return return_string;

} // end Fingerprint::return_environment();

string
Fingerprint::return_noH_environment( DOCKMol & mol, int atom_num )
{
    // Declare some temporary vectors to describe the mol without Hydrogens
    vector <int> tmp_atom_vec;
    vector < pair< pair<int,int>, string > > tmp_bond_vec;
    // Fill the vectors with atom / bond information for all atoms except hydrogens
    prepare_vectors(mol, tmp_atom_vec, tmp_bond_vec);

/*
    cout << endl;
    for (int i=0; i<tmp_atom_vec.size(); ++i){
        cout << "ATOM_VEC: " << tmp_atom_vec[i] << endl;
        cout << "ATOM_NUM: " << mol.atom_number[i] << endl << endl;
    }
    cout << "PASSED_NUM: " << atom_num << endl << endl; 
    for (int i=0; i<tmp_bond_vec.size(); ++i){
        cout << "BOND_VEC: " << tmp_bond_vec[i].first.first << " " << tmp_bond_vec[i].first.second << " " 
            << tmp_bond_vec[i].second << 
            endl;
    }
    
*/
    // Generate sets of atom keys for each molecule to different depths
    vector <string> atom_keys_1 = generate_keys( mol, 1, tmp_atom_vec, tmp_bond_vec );

    /*
    for (int i=0; i<atom_keys_1.size(); ++i){
        cout << "atom_keys: " << atom_keys_1[i] << endl;
        cout << "stringy s: " << i << endl;
    }
    */

    string return_string = atom_keys_1[atom_num];


    // Clear the atom and bond vecs for this molecule
    // alpha.clear();

    if (return_string.find("Z") != std::string::npos) {
        //cout << "THERE IS A HYDROGEN IN HERE!!!!" << endl;
        while (return_string.find("Z") != std::string::npos) {
            size_t pos = return_string.find("Z");
            string new_string = return_string.substr(0,pos-2) + return_string.substr(pos+1);
            return_string = new_string;
        }

    }
    //cout << "RETURN STRING: " << return_string << endl; 
    
    atom_keys_1.clear();
    tmp_atom_vec.clear();
    tmp_bond_vec.clear();

    return return_string;

} // end Fingerprint::return_noH_environment();


// +++++++++++++++++++++++++++++++++++++++++
// During fragment library generation for de novo growth, this function can optionally be
// called to write a list of 'allowable' torsion environments. The environments originate
// from pairs of heavy atoms about a rotatable bond.
void
Fingerprint::write_torsion_environments( DOCKMol & mol, vector <string> & tmp_torsions )
{

    // Declare some temporary vectors to describe the mol without Hydrogens
    vector <int> tmp_atom_vec;
    vector < pair< pair<int,int>, string > > tmp_bond_vec;

    // Fill the vectors with atom / bond information for all atoms except hydrogens
    prepare_noH_vectors_for_table(mol, tmp_atom_vec, tmp_bond_vec);

    // Generate sets of atom keys for each molecule to different depths
    vector <string> atom_keys_1 = generate_keys( mol, 1, tmp_atom_vec, tmp_bond_vec );

    // Print the keys to the screen
    //cout <<"\nKeys that are 1 bond deep:\n";
    //for (int i=0; i<atom_keys_1.size(); i++){ cout <<atom_keys_1[i] <<"\n"; }

    // Open up the filestream for writing torsion environments to file
    //fstream fout_torenv;
    //fout_torenv.open (filename.c_str(), fstream::out|fstream::app);

    // Take keys two at a time, sort, concatenate, and print
    for (int i=0; i<tmp_atom_vec.size(); i+=2){

        vector <string> torsion_environments;
        torsion_environments.push_back(atom_keys_1[i]);
        torsion_environments.push_back(atom_keys_1[i+1]); 

        sort(torsion_environments.begin(), torsion_environments.end());

        stringstream alpha;
        alpha << torsion_environments[0];
        alpha << "-";
        alpha << torsion_environments[1];
        string alpha_string = alpha.str();
        alpha.clear();

        tmp_torsions.push_back(alpha_string);

        // This probably needs to be a write-to-file
        //cout <<alpha_string <<"\n";
        //fout_torenv <<alpha_string <<"\n";
        //torsion_environments.clear();
    }


    // Clear the atom and bond vecs for this molecule
    // alpha.clear();
    atom_keys_1.clear();
    tmp_atom_vec.clear();
    tmp_bond_vec.clear();
    //fout_torenv.close();

    return;

} // end Fingerprint::write_torsion_environments()



// +++++++++++++++++++++++++++++++++++++++++
// Given a DOCKMol, prepare some vectors that include atom information for all atoms
// except for hydrogen. (Iterating through these vectors instead of the DOCKMol should
// speed up things). 
void
Fingerprint::prepare_noH_vectors( DOCKMol & mol, vector <int> & tmp_atom_vec, 
                                  vector < pair < pair <int,int>, string > > & tmp_bond_vec )
{

    // Make a list of non-hydrogen atom indices
    for (int i=0; i<mol.num_atoms; i++){
        if (mol.atom_types[i] != "H"){
            tmp_atom_vec.push_back(i);
        }
    }

    // cout <<"These atoms are not H: ";
    // for (int i=0; i<tmp_atom_vec.size(); i++){
    //     cout <<tmp_atom_vec[i] <<" ";
    // }
    // cout <<"\n\n";


    // Make a list of non-hydrogen involved bonds and bond types. The data structure is a little weird - it is
    // a vector that contains three values in each position: 'first.first', 'first.second', and 'second'.
    // 'first.first' and 'first.second' are indices (ints) of bonded atoms, 'second' is the bond type (string)
    for (int i=0; i<mol.num_bonds; i++){
        if (mol.atom_types[mol.bonds_origin_atom[i]] != "H" && 
            mol.atom_types[mol.bonds_target_atom[i]] != "H")   {

            // Adding onto tmp_bond_vec (expanded for readability)
            tmp_bond_vec.push_back( 
                make_pair(
                    make_pair(
                        mol.bonds_origin_atom[i],
                        mol.bonds_target_atom[i]
                    ), 
                    mol.bond_types[i] 
                )
            );
            // End adding onto tmp_bond_vec
        }
    }

    // for (int i=0; i<tmp_bond_vec.size(); i++){
    //     cout <<"Atom " <<tmp_bond_vec[i].first.first <<" is bound to atom " <<tmp_bond_vec[i].first.second 
    //          <<" and it is type " <<tmp_bond_vec[i].second <<"\n";
    // } cout <<"\n";

    return;

} // end Fingerprint::prepare_noH_vectors()



// +++++++++++++++++++++++++++++++++++++++++
//
void
Fingerprint::prepare_vectors( DOCKMol & mol, vector <int> & tmp_atom_vec,
                                  vector < pair < pair <int,int>, string > > & tmp_bond_vec )
{

    // Make a list of all atoms
    for (int i=0; i<mol.num_atoms; i++){
            tmp_atom_vec.push_back(i);
    }


    // Make a list of non-hydrogen involved bonds and bond types. The data structure is a little weird - it is
    // a vector that contains three values in each position: 'first.first', 'first.second', and 'second'.
    // 'first.first' and 'first.second' are indices (ints) of bonded atoms, 'second' is the bond type (string)
    for (int i=0; i<mol.num_bonds; i++){

            // Adding onto tmp_bond_vec (expanded for readability)
            tmp_bond_vec.push_back(
                make_pair(
                    make_pair(
                        mol.bonds_origin_atom[i],
                        mol.bonds_target_atom[i]
                    ),
                    mol.bond_types[i]
                )
            );
            // End adding onto tmp_bond_vec
    }

    // for (int i=0; i<tmp_bond_vec.size(); i++){
    //     cout <<"Atom " <<tmp_bond_vec[i].first.first <<" is bound to atom " <<tmp_bond_vec[i].first.second
    //          <<" and it is type " <<tmp_bond_vec[i].second <<"\n";
    // } cout <<"\n";

    return;

} // end Fingerprint::prepare_vectors()



// +++++++++++++++++++++++++++++++++++++++++
// Given a DOCKMol, prepare vectors of atom and bond information that does not include
// hydrogen. In this special case, only return atoms that are part of a rotatable bond.
void
Fingerprint::prepare_noH_vectors_for_table( DOCKMol & mol, vector <int> & tmp_atom_vec, 
                                  vector < pair < pair <int,int>, string > > & tmp_bond_vec )
{

    // Make a list of non-hydrogen involved bonds and bond types. The data structure is a little weird - it is
    // a vector that contains three values in each position: 'first.first', 'first.second', and 'second'.
    // 'first.first' and 'first.second' are indices (ints) of bonded atoms, 'second' is the bond type (string)
    for (int i=0; i<mol.num_bonds; i++){
        
        // If that bond is rotatable, remember target and origin atom indices
        if (mol.bond_is_rotor(i)){
            //cout <<mol.bonds_origin_atom[i] <<"\t" <<mol.bonds_target_atom[i] <<"\n";
            tmp_atom_vec.push_back(mol.bonds_origin_atom[i]);
            tmp_atom_vec.push_back(mol.bonds_target_atom[i]);
        }


        if (mol.atom_types[mol.bonds_origin_atom[i]] != "H" && 
            mol.atom_types[mol.bonds_target_atom[i]] != "H")   {

            // Adding onto tmp_bond_vec (expanded for readability)
            tmp_bond_vec.push_back( 
                make_pair(
                    make_pair(
                        mol.bonds_origin_atom[i],
                        mol.bonds_target_atom[i]
                    ), 
                    mol.bond_types[i] 
                )
            );
            // End adding onto tmp_bond_vec
        }
    }

    return;

} // end Fingerprint::prepare_noH_vectors_for_table()



// +++++++++++++++++++++++++++++++++++++++++
// Given a DOCKMol, prepare vectors of atom and bond information that does not include
// hydrogen. In this special case, only return atoms that are active.
void
Fingerprint::prepare_noH_vectors_for_active_atoms( DOCKMol & mol, vector <int> & tmp_atom_vec, 
                                  vector < pair < pair <int,int>, string > > & tmp_bond_vec )
{
    // Add only active atoms to this vector (should only be two during de novo growth)
    for (int i=0; i<mol.num_atoms; i++){
        //if (mol.atom_active_flags[i]){
        if (mol.atom_types[i] != "H" && mol.atom_active_flags[i]){

            tmp_atom_vec.push_back(i);
        }
    }


    // Make a list of non-hydrogen involved bonds and bond types. The data structure is a little weird - it is
    // a vector that contains three values in each position: 'first.first', 'first.second', and 'second'.
    // 'first.first' and 'first.second' are indices (ints) of bonded atoms, 'second' is the bond type (string)
    for (int i=0; i<mol.num_bonds; i++){

        if (mol.atom_types[mol.bonds_origin_atom[i]] != "H" && 
            mol.atom_types[mol.bonds_target_atom[i]] != "H")   {

        if (mol.bond_active_flags[i]){


            // Adding onto tmp_bond_vec (expanded for readability)
            tmp_bond_vec.push_back( 
                make_pair(
                    make_pair(
                        mol.bonds_origin_atom[i],
                        mol.bonds_target_atom[i]
                    ), 
                    mol.bond_types[i] 
                )
            );
            // End adding onto tmp_bond_vec
        }
        }
    }

    return;

} // end Fingerprint::prepare_noH_vectors_for_active_atoms()




// +++++++++++++++++++++++++++++++++++++++++
// Generate keys to a specified depth for all of the non-H atoms in a mol
vector <string>
Fingerprint::generate_keys( DOCKMol & mol, int depth, vector <int> & tmp_atom_vec,
                            vector < pair< pair<int,int>, string > > & tmp_bond_vec )
{

    // Create a temporary vector of atom_keys
    vector <string> tmp_atom_keys;

    // For every atom that isn't a hydrogen, form a connection map and write a string
    for (int i=0; i<tmp_atom_vec.size(); i++){

        // Create the root node (atom i) and assign some values 
        Node root;
        root.atom_index = tmp_atom_vec[i];
        root.atom_type = return_atom_type( mol.atom_types[ tmp_atom_vec[i] ] );

        // Then make a connection map for that root node up to 'depth' bonds away
        populate_tree(root, depth, mol, tmp_atom_vec, tmp_bond_vec);

        // Sort the tree according to 8 sorting rules (see below) - makes the tree canonical
        sort_tree(root);

        // Function to print tree
        // cout <<"Sorted tree: \n";
        // print_tree(root, 1);

        // The string will be wrtten to a stringstream as the tree is traversed
        stringstream ss;
        write_string(root, ss);
        // End each string with a pound sign
        ss <<"#";

        // Push each atom key onto the vector of atom keys
        tmp_atom_keys.push_back(ss.str());
        ss.clear();
 

        // Note: I don't think we can use the class destructor to recursively deallocate
        // the memory of child vectors for a given root node. Something about as soon 
        // as the destructor is called for a given object (root, in this case), then the
        // lifetime of root has ended. Therefore, you can't call a destructor within
        // a destructor because the behavior is now undefined. Simple fix - write a 
        // function to destruct a recursive object rather than let that object pass out
        // of scope.
        destruct(root);

    }


    return tmp_atom_keys;

} // end Fingerprint::generate_keys()



// +++++++++++++++++++++++++++++++++++++++++
// A recursive function to create a tree given a starting node and a depth
void
Fingerprint::populate_tree( Node & node_ref, int depth, DOCKMol & mol, vector <int> & tmp_atom_vec, 
                            vector < pair < pair <int,int>, string > > & tmp_bond_vec )
{

    // Continue populating the tree so long as depth!=0 when this function is entered. Given a different
    // starting 'depth', trees of different complexity can be formed. As an example, given a depth of '3',
    // the tree will contain a root node and paths through all atoms up to 3 bonds away.
    if (depth <= 0){ return; }
    depth--;

    // Look through every bond and identify those that include node_ref.atom_index
    for (int i=0; i<tmp_bond_vec.size(); i++){

        // If a bond is found that contains the reference node atom index, and the partner in the bond
        // is not the parent's atom index, then create a new child node.
        if (tmp_bond_vec[i].first.first  == node_ref.atom_index  && 
            tmp_bond_vec[i].first.second != node_ref.parent_index  ){

            // Assign some values to the child node
            Node tmp_child;
            tmp_child.atom_index = tmp_bond_vec[i].first.second;
            tmp_child.parent_index = node_ref.atom_index;
            tmp_child.atom_type = return_atom_type( mol.atom_types[ tmp_bond_vec[i].first.second ] );
            tmp_child.bond_to_parent = return_bond_type( tmp_bond_vec[i].second ); 

            // Then recursively enter this function for that child node (until depth <= 0)
            populate_tree( tmp_child, depth, mol, tmp_atom_vec, tmp_bond_vec);

            // Once a node's children have all been assigned, push that node onto its parent's vector
            node_ref.children.push_back(tmp_child);
        }

        // If a bond is found that contains the reference node atom index, and the partner in the bond
        // is not the parent's atom index, then create a new child node.
        if (tmp_bond_vec[i].first.second == node_ref.atom_index  && 
            tmp_bond_vec[i].first.first  != node_ref.parent_index  ){

            // Assign some values to the child node
            Node tmp_child;
            tmp_child.atom_index = tmp_bond_vec[i].first.first;
            tmp_child.parent_index = node_ref.atom_index;
            tmp_child.atom_type = return_atom_type( mol.atom_types[ tmp_bond_vec[i].first.first ] );
            tmp_child.bond_to_parent = return_bond_type( tmp_bond_vec[i].second );

            // Then recursively enter this function for that child node (until depth <= 0)
            populate_tree( tmp_child, depth, mol, tmp_atom_vec, tmp_bond_vec);

            // Once a node's children have all been assigned, push that node onto its parent's vector
            node_ref.children.push_back(tmp_child);
        }
    }

    // When all children for a given node have been defined and added to the child_vector, go back
    // up a level in the tree
    return;

} // end Fingerprint::populate_tree()



// +++++++++++++++++++++++++++++++++++++++++
// Given a DOCKMol atom_type, returns a string corresponding to a sybyl atom_type
string
Fingerprint::return_atom_type( string atom_type )
{

    // The form of this function makes it very easy to add new atom types. After all of
    // the capital letters run out, lowercase letters can be used. Lowercase letters
    // should follow the form:
    //       A < B < C ... < Z < a < b < c ... < z
    // allowing for 52 atom types.

    if      ( atom_type == "C.3"   ){ return "A"; }
    else if ( atom_type == "C.2"   ){ return "B"; }
    else if ( atom_type == "C.1"   ){ return "C"; }
    else if ( atom_type == "C.ar"  ){ return "D"; }
    else if ( atom_type == "C.cat" ){ return "E"; }
    else if ( atom_type == "N.3"   ){ return "F"; }
    else if ( atom_type == "N.2"   ){ return "G"; }
    else if ( atom_type == "N.1"   ){ return "H"; }
    else if ( atom_type == "N.ar"  ){ return "I"; }
    else if ( atom_type == "N.am"  ){ return "J"; }
    else if ( atom_type == "N.pl3" ){ return "K"; }
    else if ( atom_type == "N.4"   ){ return "L"; }
    else if ( atom_type == "O.3"   ){ return "M"; }
    else if ( atom_type == "O.t3p" ){ return "M"; }
    else if ( atom_type == "O.2"   ){ return "N"; }
    else if ( atom_type == "O.co2" ){ return "O"; }
    else if ( atom_type == "S.3"   ){ return "P"; }
    else if ( atom_type == "S.2"   ){ return "Q"; }
    else if ( atom_type == "S.O"   ){ return "R"; } // Both of these seem to
    else if ( atom_type == "S.o"   ){ return "R"; } // occur in Zinc mols
    else if ( atom_type == "S.O2"  ){ return "S"; }    // These ones
    else if ( atom_type == "S.o2"  ){ return "S"; }    // too
    else if ( atom_type == "P.3"   ){ return "T"; }
    else if ( atom_type == "F"     ){ return "U"; }
    else if ( atom_type == "Cl"    ){ return "V"; }
    else if ( atom_type == "Br"    ){ return "W"; }
    else if ( atom_type == "I"     ){ return "X"; }
    else if ( atom_type == "Du"    ){ return "Y"; }
    else if ( atom_type == "H"     ){ return "Z"; }
    else if ( atom_type == "S"     ){ return "a"; } // ZINC15 dropped S.o and S.o2
    else if ( atom_type == "P"     ){ return "b"; } // ZINC15 dropped P.3
    else if ( atom_type == "Si"    ){ return "c"; }

    else {
        cout <<"Warning: atom type " <<atom_type <<" in Fingerprint::return_atom_type not recognized\n";
    }


    // If you return z, something bad happened
    return "z";

} // end Fingerprint::return_atom_type();



// +++++++++++++++++++++++++++++++++++++++++
// Given a DOCKMol bond_type, returns an int corresponding to a sybyl bond_type
int
Fingerprint::return_bond_type( string bond_type )
{

    // The form of this function makes it very easy to add new bond types. The way the sort function is
    // structured right now, however, we should only use 1-digit ints in this function. In other words,
    // if there happens to be a bunch of new bond types, don't define them using '10', '11', etc..
    // Instead, it would probably probably be easier to change the type to strings and start returning
    // some punctuation marks that follow these digits in the ASCII table.

    if      ( bond_type == "1"  ){ return 1; }
    else if ( bond_type == "2"  ){ return 2; }
    else if ( bond_type == "3"  ){ return 3; }
    else if ( bond_type == "am" ){ return 4; }
    else if ( bond_type == "ar" ){ return 5; }
    else if ( bond_type == "du" ){ return 6; }

    else { 
        cout <<"Warning: bond type " <<bond_type <<" in Fingerprint::return_bond_type not recognized\n";
    }

    // If you return 9, something bad happened
    return 9;

} // end Fingerprint::return_bond_type();



// +++++++++++++++++++++++++++++++++++++++++
// Given the root node of a tree, sort it in a canonical way
// (see 8 sorting rules in tree_sort_function())
void
Fingerprint::sort_tree( Node & tmp_node )
{

    // If this function has reached a leaf - or a node without children, there is nothing
    // below it to sort, so go back up one step.
    if (tmp_node.children.size() == 0) { return; }

    // If there are children, enter into this function for each child recursively until you 
    // reach a node that only contains leaves (or, later, sorted branches)
    for (int i=0; i<tmp_node.children.size(); i++){
        sort_tree(tmp_node.children[i]);
    }

    // If this branch has more than one child, sort according to tree_sort_function()
    if (tmp_node.children.size() > 1){
        sort( tmp_node.children.begin(), tmp_node.children.end(), tree_sort_function);
    }


    // Now that this branch's leaves are sorted, assign some weights to it. This way it
    // can also be treated like a leaf for sorting purposes

    stringstream str_bond_types;
    stringstream str_grandchildren_num;
    stringstream str_grandchildren_bond_types;

    // For each child of the current node, 
    for (int i=0; i<tmp_node.children.size(); i++){
        // Append their atom_types together and add to tmp_node.children_identities
        tmp_node.children_identities.append( tmp_node.children[i].atom_type );

        // Append their bond types together and add to tmp_node.children_bond_types
        str_bond_types <<tmp_node.children[i].bond_to_parent;

        // Append the number of children each has together
        str_grandchildren_num <<tmp_node.children[i].children.size();

        // Append the identities of their grandchildren together
        tmp_node.grandchildren_identities.append( tmp_node.children[i].children_identities );

        // Append the bond types between children and grandchildren together
        str_grandchildren_bond_types <<tmp_node.children[i].children_bond_types;
    }

    tmp_node.children_bond_types = str_bond_types.str();
    tmp_node.grandchildren_num = str_grandchildren_num.str();
    tmp_node.grandchildren_bond_types = str_grandchildren_bond_types.str();


    // This should now be a 'sorted branch', and for the purposes of sorting, can be treated
    // as a leaf. Go back up the tree one step.
    return;

} // end Fingerprint::sort_tree()



// +++++++++++++++++++++++++++++++++++++++++
// Sort a tree based on 8 sorting rules
bool 
tree_sort_function(const Node & node1, const Node & node2)
{
    /////////////////////////////////////////////////////////////
    // Sorting Rules:                                          //
    // 1.) Identity of node              | low  |      |       //
    // 2.) Bond type to node             |      |      |       //
    // 3.) # of children                        | mid  |       //
    // 4.) Identity of children                 |      |       //
    // 5.) Bond types to children               |      | high  //
    // 6.) # of grand-children                         |       //
    // 7.) Identity of grand-children                  |       //
    // 8.) Bond types to grand-children                |       //
    /////////////////////////////////////////////////////////////


    if (node1.atom_type < node2.atom_type)
        return true;
    else if (node2.atom_type < node1.atom_type)
        return false;
    else

      if (node1.bond_to_parent < node2.bond_to_parent)
          return true;
      else if (node2.bond_to_parent < node1.bond_to_parent)
          return false;
      else

        if (node1.children.size() > node2.children.size())
            return true;
        else if (node2.children.size() > node1.children.size())
            return false;
        else

          if (node1.children_identities < node2.children_identities)
              return true;
          else if (node2.children_identities < node1.children_identities)
              return false;
          else

            if (node1.children_bond_types < node2.children_bond_types)
                return true;
            else if (node2.children_bond_types < node1.children_bond_types)
                return false;
            else

              if (node1.grandchildren_num > node2.grandchildren_num)
                  return true;
              else if (node2.grandchildren_num > node1.grandchildren_num)
                  return false;
              else

                if (node1.grandchildren_identities < node2.grandchildren_identities)
                    return true;
                else if (node2.grandchildren_identities < node1.grandchildren_identities)
                    return false;
                else

                  if (node1.grandchildren_bond_types < node2.grandchildren_bond_types)
                      return true;
                  else if (node2.grandchildren_bond_types < node1.grandchildren_bond_types)
                      return false;
    
    return false;

} // end tree_sort_function();



// +++++++++++++++++++++++++++++++++++++++++
// Given the root node of a sorted tree, print the canonical string 
void
Fingerprint::write_string( Node & tmp_node, stringstream & ss )
{

    // Write to the stringstream the bond type between the current node and the parent (which will
    // be 0 for the root node only), the the atom type of the current node
    ss <<tmp_node.bond_to_parent <<tmp_node.atom_type;

    // So long as there are children, enter them in order (already sorted) and continue printing
    for (int i=0; i<tmp_node.children.size(); i++){

         write_string(tmp_node.children[i], ss);
    }

    // Write a caret every time you go back up a level
    ss <<"^";

    return;

} // end Fingerprint::write_string()



// +++++++++++++++++++++++++++++++++++++++++
// Recursively deallocate the memory of child vectors for a given tree 
void
Fingerprint::destruct( Node & tmp_node )
{

    // If there are no children, there is nothing to destruct
    if (tmp_node.children.size() == 0){ return; }

    // Recursively enter this function until you get to the lowest nodes in the tree
    for (int i=0; i<tmp_node.children.size(); i++){
        destruct(tmp_node.children[i]);
    }

    // Destruct those ones first, then go back up a step
    tmp_node.children.clear();

    return;

} // end Fingerprint::destruct()



// +++++++++++++++++++++++++++++++++++++++++
// Computes part of the terms needed for tanimoto for a specific subset of atom keys.
// Breaking the problem up into subsets like this should reduce the complexity.
void
Fingerprint::running_tanimoto( int & tc_both, int & tc_mol1, int & tc_mol2, 
                               vector <string> & str1, vector <string> & str2 )
{

    // Create a new string, str3, and copy both str1 and str2 onto it. This string
    // represents all of the keys from the two molecules.
    vector <string> str3 = str1;
    vector <string> str1_origin = str1;
    vector <string> str2_origin = str2;
    str3.insert( str3.end(), str2.begin(), str2.end() );

    //cout <<"mol1 (str1):\n";
    //for (int i=0; i<str1.size(); i++){ cout <<str1[i] <<"\n"; }
    //cout <<"mol2 (str2):\n";
    //for (int i=0; i<str2.size(); i++){ cout <<str2[i] <<"\n"; }
    //cout <<"mol2 (str3):\n";
    //for (int i=0; i<str3.size(); i++){ cout <<str3[i] <<"\n"; }


    // Iterate through every atom key in str3
    for (int i=0; i<str3.size(); i++){

        // Begin by setting both booleans to false (these represent whether the atom key
        // for this iteration is found in either of the reference molecules)
        bool temp1 = false;
        bool temp2 = false;

        // Iterate through every atom key in the first molecule
        for (int j=0; j<str1.size(); j++){

            // If the key (i) is found in this molecule...
            if (str3[i].compare(str1[j]) == 0){

                // ...set the correct boolean to true
                temp1 = true;

                // and delete the key from the molecule so it isn't counted twice
                str1.erase(str1.begin() + j);

                // Break the inner for loop
                break;
            }
        }

        // Iterate through every atom key in the second molecule
        for (int j=0; j<str2.size(); j++){

            // If the key (i) is found in this molecule...
            if (str3[i].compare(str2[j]) == 0){

                // ...set the correct boolean to true
                temp2 = true;

                // and delete the key from the molecule so it isn't counted twice
                str2.erase(str2.begin() + j);

                // Break the inner for loop
                break;
            }
        }


        // Evaluate temp1 and temp2, and increment the appropriate (reference to) tanimoto
        // coefficient variable
        if (temp1 && temp2){
            // Both molecules have the key (i)
            tc_both++;
        }

        else if (temp1 && !temp2){
            // Only mol1 has the key (i)
            tc_mol1++;
        }

        else if (!temp1 && temp2){
            // Only mol2 has the key (i)
            tc_mol2++;
        }
    }

    // Clear the vector data that is no longer needed. Note: str1 and str2 were passed here as
    // references, so the original vectors are being cleared here.
    str1 = str1_origin;
    str2 = str2_origin;
    str3.clear();


    return;

} // end Fingerprint::running_tanimoto()



////////////////////////////////////////////////////////////////////////////////
// The following functions are for de-bugging purposes                        //
////////////////////////////////////////////////////////////////////////////////

// +++++++++++++++++++++++++++++++++++++++++
// Given the root node, print a graphical depiction of the whole tree
void
Fingerprint::print_tree( Node & tmp_node, int depth )
{

    // Assuming a depth of '3', print a number of leading spaces so the data structure
    // appears more tree-like
    for (int j=0; j<(3-depth); j++){ cout <<"    "; }

    // Print the atom index followed by bond type and atom type
    cout <<"|-- " <<tmp_node.atom_index <<".(" <<tmp_node.bond_to_parent <<"." <<tmp_node.atom_type <<")\n";

    // If the node that was just printed has no children, go back up a level
    if (tmp_node.children.size() == 0){ return; }

    // Remember the depth during all stages in order to print an accurate-looking tree
    if (depth <= 0){ return; }
    depth--;

    // Recursively enter this function for every child of the current node
    for (int i=0; i<tmp_node.children.size(); i++){
        print_tree(tmp_node.children[i], depth);
    }

    // Once all children have been explored at a given level, go back up a level
    return;

} // end Fingerprint::print_tree()



// +++++++++++++++++++++++++++++++++++++++++
// Print all of the elements of a single 'Node' data structure 
void
Fingerprint::print_node_info( Node & tmp_node )
{

    cout <<"Atom index: " <<tmp_node.atom_index <<"\n";
    cout <<"Atom type: " <<tmp_node.atom_type <<"\n";
    cout <<"Parent index: " <<tmp_node.parent_index <<"\n";
    cout <<"Bond to parent: " <<tmp_node.bond_to_parent <<"\n";
    cout <<"Num of children: " <<tmp_node.children.size() <<"\n";
    cout <<"Children identities: " <<tmp_node.children_identities <<"\n";
    cout <<"Children bond types: " <<tmp_node.children_bond_types <<"\n";
    cout <<"Num of grand children: " <<tmp_node.grandchildren_num <<"\n";
    cout <<"Grand children identities: " <<tmp_node.grandchildren_identities <<"\n";
    cout <<"Grand children bond types: " <<tmp_node.grandchildren_bond_types <<"\n\n";


    return;

} // end Fingerprint::print_node_info()



// +++++++++++++++++++++++++++++++++++++++++
// Print all of the elements of the 'Node' data structure for each node in a tree 
void
Fingerprint::print_tree_info( Node & tmp_node )
{

    cout <<"Atom index: " <<tmp_node.atom_index <<"\n";
    cout <<"Atom type: " <<tmp_node.atom_type <<"\n";
    cout <<"Parent index: " <<tmp_node.parent_index <<"\n";
    cout <<"Bond to parent: " <<tmp_node.bond_to_parent <<"\n";
    cout <<"Num of children: " <<tmp_node.children.size() <<"\n";
    cout <<"Children identities: " <<tmp_node.children_identities <<"\n";
    cout <<"Children bond types: " <<tmp_node.children_bond_types <<"\n";
    cout <<"Num of grand children: " <<tmp_node.grandchildren_num <<"\n";
    cout <<"Grand children identities: " <<tmp_node.grandchildren_identities <<"\n";
    cout <<"Grand children bond types: " <<tmp_node.grandchildren_bond_types <<"\n\n";


    for (int i=0; i<tmp_node.children.size(); i++){

         print_node_info(tmp_node.children[i]);
    }

    return;

} // end Fingerprint::print_tree_info()



