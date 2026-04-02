//
#ifndef FINGERPRINT_H
#define FINGERPRINT_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include "dockmol.h"

class Parameter_Reader;


// +++++++++++++++++++++++++++++++++++++++++
// Node class - for building molecular connection maps 
class Node {

  public:

    Node();        // Constructor - set all things to a null character or impossible value
    ~Node();       // Destructor - empty the children vector

    int                 atom_index;                // Index of atom according to DOCKMol
    std::string         atom_type;                 // Atom type stored as a single letter
    int                 parent_index;              // Index of parent atom according to DOCKMol
    int                 bond_to_parent;            // Bond type between this atom and parent atom

    std::vector <Node>  children;                  // Vector of children (nodes) that this node is connected to
    std::string         children_identities;       // Concatenated string of (sorted) children atom types
    std::string         children_bond_types;       // Concatenated string of (sorted) children bond types

    std::string         grandchildren_num;         // Concatenated string of number of children each child has
    std::string         grandchildren_identities;  // Concatenated string of grandchildren atom types
    std::string         grandchildren_bond_types;  // Concatenated string of grandchildren bond types
};


// +++++++++++++++++++++++++++++++++++++++++
// Fingerprint class - for writing molecular bitstrings
class Fingerprint {

  public:

    Fingerprint();      // Constructor - does nothing right now
    ~Fingerprint();     // Destructor - does nothing right now

    // Input parameters
    bool            dbfilter_compute_fingerprint;
    float           dbfilter_tanimoto_cutoff;
    std::string     dbfilter_tanimoto_reference;
    DOCKMol         ref_mol;

    // Functions
    void            input_parameters(Parameter_Reader &);

    // There are the different ways to call the fingerprint functions depending on the goal
    std::string   compute_unique_string( DOCKMol & );
    std::string   compute_unique_string_active( DOCKMol & );
    float         compute_tanimoto( DOCKMol &, DOCKMol & );
    float         compute_tanimoto_fragment( DOCKMol &, DOCKMol & );
    void          running_tanimoto( int &, int &, int &, std::vector <std::string> &, std::vector <std::string> & );
    void          write_torsion_environments( DOCKMol &, std::vector <std::string> & );
    void          return_torsion_environments( DOCKMol &, std::pair <std::string, std::string> & );
    std::string   return_environment( DOCKMol &, int );
    std::string   return_noH_environment( DOCKMol&, int); //JDB


    // Goes through the DOCKMol and writes two vectors of necessary atom / bond information
    void   prepare_vectors(DOCKMol &, std::vector <int> &,
                               std::vector < std::pair < std::pair <int,int>, std::string > > &);
    void   prepare_noH_vectors(DOCKMol &, std::vector <int> &, 
                               std::vector < std::pair < std::pair <int,int>, std::string > > &);
    void   prepare_noH_vectors_for_table(DOCKMol &, std::vector <int> &, 
                               std::vector < std::pair < std::pair <int,int>, std::string > > &);
    void   prepare_noH_vectors_for_active_atoms(DOCKMol &, std::vector <int> &, 
                               std::vector < std::pair < std::pair <int,int>, std::string > > &);

    // Generate a set of atomic keys for a molecule given a certain depth
    std::vector <std::string> generate_keys( DOCKMol &, int, std::vector <int> &,
                                             std::vector < std::pair< std::pair<int,int>, std::string > > & );

    // Given the atom/bond vectors, a root atom, and a depth, creates a connection map
    void   populate_tree(Node &, int, DOCKMol &, std::vector <int> &, 
                         std::vector < std::pair <std::pair <int,int>, std::string> > &);

    // These functions help create and sort the tree
    std::string  return_atom_type(std::string);        // Converts SYBYL atomtype to a letter
    int          return_bond_type(std::string);        // Converts SYBYL bondtype to an int

    void   sort_tree(Node &);                          // Make the tree canonical
    void   write_string(Node &, std::stringstream &);  // Write a string that represents connection map
    void   destruct(Node &);                           // Deallocate memory of child vectors


    // These functions are not necessary, but useful for debugging
    void   print_tree( Node &, int );                  // Print a graphical depiction of a tree
    void   print_node_info( Node & );                  // Print all values stored in a node
    void   print_tree_info( Node & );                  // Print all values stored in a tree

};


// +++++++++++++++++++++++++++++++++++++++++
// Auxilliary functions 
bool tree_sort_function(const Node &, const Node &);   // Checks hierarchy of eight tree-sorting rules

#endif // FINGERPRINT_H
