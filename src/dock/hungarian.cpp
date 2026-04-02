#include "hungarian.h"


// +++++++++++++++++++++++++++++++++++++++++
// Function to calculate hungarian (symmetry corrected) RMSD
//
double Hungarian_RMSD::calc_Hungarian_RMSD(DOCKMol & refmol, DOCKMol & mol){

    if (refmol.num_atoms != mol.num_atoms){		// check to make sure reference molecule
         rmsd =  -1000.0;  		  		// and pose molecule have same # of atoms
         return rmsd;					// if not - return an invalid rmsd
    }

    NUM = mol.num_atoms;		// number of atoms in the molecule
    MAX = 2147483647;			// A large double used to find minimum values
    total_assignment = 0;		// the 'cost' of the solved matrix
    rmsd = 0;				// the rmsd computed from cost

    unique_atom_types = new string [NUM];	// A list of unique atom types in mol
    num_unique = 0;				// the number of unique atom types

    // Compile a list of unique atom types that occur in mol - ignoring hydrogen
    bool unique_flag = true;
    for (int i=0; i<NUM; i++){						// for each atom in mol
      if (mol.atom_types[i].compare("H") != 0){				// (ignore hydrogens)
        for (int j=0; j<num_unique+1; j++){				// for every entry in unique_atom_types
            if (mol.atom_types[i].compare(unique_atom_types[j]) == 0){	// if that atom type has already been
                unique_flag = false;					// added to the list, then break the 
                break;							// inner for-loop
            }
	}
        if (unique_flag){						// if that atom type is not found in
            unique_atom_types[num_unique] = mol.atom_types[i];		// the list, then add it to the list
            num_unique++;						// in the next position
        }
        unique_flag = true;						// reset this flag and start over with
      }									// next atom in mol
    }

    // Go through unique_atom_types one by one, forming a distance matrix for each set of
    // atoms of the same atom type. Use hungarian to solve each matrix, and add the assignment
    // to total_assignment
    for (int z=0; z<num_unique; z++){							

        num_type = 0;		// the number of times unique_atom_types[z] occurs in mol

        for (int i=0; i<NUM; i++){					// for every atom in mol, check the atom type
            if (mol.atom_types[i].compare(unique_atom_types[z]) == 0){  // and increment num_type if it matches the
                num_type++;						// atom type for this iteration of the loop (z)
            }
        }

        //cout <<endl <<"Iteration " <<z+1 <<" : atomtype " <<unique_atom_types[z] <<" : num_type " <<num_type;

        int num_type2 = 0;    // temporary counter for the next check

        for (int i=0; i<NUM; i++){                                         // for every atom in refmol,
            if (refmol.atom_types[i].compare(unique_atom_types[z]) == 0){  // if the atom type is equal to 'z',
                num_type2++;                                               // then increment the counter
            }
        }

        if (num_type != num_type2){          // if the refmol and mol do not have the same number
            cout <<"Warning: atom type " <<unique_atom_types[z] <<" has " <<num_type <<" in the reference and "
                 <<num_type2 <<" in the mol object " <<endl;
            delete[] unique_atom_types;      // of a certain atomtype, return an invalid RMSD
            return -1000.0;
        }

        initialize();		// initialize the matrices to the appropriate size for this iteration of the loop (z)

        // This loop goes through each element in mol and refmol - if they both match this rounds' atom type (z),
        // then the distance between them is computed, and it is added to next cell in matrix[][]
        int count_i = -1;	// counter for forming matrix[][]
        int count_j = -1;	// counter for forming matrix[][]

        for (int i=0; i<NUM; i++){							// for each atom in mol,
            if (mol.atom_types[i].compare(unique_atom_types[z]) == 0){			// if the atom_type matches,
                count_i++;								// move to the appropriate row
                count_j = -1;								// in matrix[][]
                for (int j=0; j<NUM; j++){						// then for each atom in refmol,
                    if (refmol.atom_types[j].compare(unique_atom_types[z]) == 0){	// if the atom type matches,
                        count_j++;							// move to the appropriate column
                        matrix[count_i][count_j] = ((mol.x[i] - refmol.x[j])*(mol.x[i] - refmol.x[j])) +	// in matrix[][] and enter the 
                                                   ((mol.y[i] - refmol.y[j])*(mol.y[i] - refmol.y[j])) + 	// the distance between those atoms
                                                   ((mol.z[i] - refmol.z[j])*(mol.z[i] - refmol.z[j]));
			// There is some debate about whether to use a pure distance matrix, or a matrix of squared
			// distances. As it is now, this matrix is a matrix of squared distances. Using this matrix,
			// the algorithm gives a weighted preference to matches of shorter distances and a weighted 
			// cost to matches of longer distances. In many test cases, the matrix of squared distances
			// always returned if not the same, a slighty smaller total assignment than did the matrix of
			// pure distances due to slightly different final matching.
                    } 
                }

		//cout <<mol.atom_names[i] <<"=" <<count_i <<",  ";

            } 
        } 

        hungarian();					// hungarian algorithm assigns worker-job pairs for matrix[][],
							// the results are saved in matrix_match[]
        double assignment = sum_assignment();		// then the assignmment (cost) of matrix_match[] is computed
        total_assignment += assignment;			// and added to the total cost

	//cout <<endl <<"(refmol : mol)\n";
	//for (int i=0; i<num_type; i++){
	//	cout <<i <<" : " <<matrix_match[i] <<endl;
	//}
	//cout <<"assignment = " <<assignment <<endl;

        clear(); 				// de-allocate memory of the arrays

    } // this for-loop repeats for every atom type in mol (excluding H), the end result is total_assignment

    int heavy_count = 0;			// number of heavy atoms in mol

    for (int i=0; i<NUM; i++){ 			// for each atom in mol,
        if (mol.amber_at_heavy_flag[i]){ 	// if the atom is a 'heavy atom' (not hydrogen), then
            heavy_count++;			// increment the counter 'heavy_count'
        }
    }

    rmsd = sqrt( (total_assignment/heavy_count) );	// compute symmetry-corrected rmsd

    //cout <<endl <<endl <<"total_assignment = " <<total_assignment <<endl;
    //cout <<"number of heavy atoms = " <<heavy_count <<endl;
    //cout <<"rmsd = " <<rmsd <<endl;

    delete[] unique_atom_types;				// de-allocate memory of unique_atom_types
    return rmsd;					// return the final number

} // end Hungarian_RMSD::calc_Hungarian_RMSD() 



// +++++++++++++++++++++++++++++++++++++++++
// Function to calculate hungarian RMSD for two dissimilar molecules. The return values are the
// symmetry-corrected RMSD for the matched atoms as a double and the number of unmatched-atoms 
// as an int.
// Modified to calculate hungarian RMSD for active atoms only - CS: 2016-06-23
pair <double, int> Hungarian_RMSD::calc_Hungarian_RMSD_dissimilar(DOCKMol & refmol, DOCKMol & mol){

    //cout <<endl <<endl <<"########## Entering calc_Hungarian_RMSD_dissimilar " <<endl;

    int num_unmatched = 0;
    int num_matched = 0;
    int num_matched_total = 0;
    MAX = 2147483647;
    total_assignment = 0;
    rmsd = 0;

    int num_unique_refmol = 0;
    int num_unique_mol = 0;
    int num_unique = 0;
    string * unique_atom_types_refmol = new string [refmol.num_atoms];
    string * unique_atom_types_mol = new string [mol.num_atoms];
    vector <string> unique_atom_types_vector;


    // Compile a list of unique atom types that occur in refmol - ignoring hydrogen
    bool unique_flag = true;
    for (int i=0; i<refmol.num_atoms; i++){
        // CS:2016-06-23, Adjusted to count active atoms only
        if (refmol.atom_active_flags[i] == true){
            if (refmol.atom_types[i].compare("H") != 0){
              for (int j=0; j<num_unique_refmol+1; j++){
                  if (refmol.atom_types[i].compare(unique_atom_types_refmol[j]) == 0){
                     unique_flag = false; 
                     break;
                  }
              }
              if (unique_flag){
                 unique_atom_types_refmol[num_unique_refmol] = refmol.atom_types[i];
                 num_unique_refmol++;
              }
              unique_flag = true;
            }
        }
    }

    //cout <<num_unique_refmol <<" unique atom types that occur in refmol: ";

    // Compile a list of unique atom types that occur in mol - ignoring hydrogen
    unique_flag = true;
    for (int i=0; i<mol.num_atoms; i++){
        // CS:2016-06-23, Adjusted to count active atoms only
        if (mol.atom_active_flags[i] == true){
            if (mol.atom_types[i].compare("H") != 0){
               for (int j=0; j<num_unique_mol+1; j++){
                   if (mol.atom_types[i].compare(unique_atom_types_mol[j]) == 0){
                       unique_flag = false; 
                       break;
                   }
               }
               if (unique_flag){
                   unique_atom_types_mol[num_unique_mol] = mol.atom_types[i];
                   num_unique_mol++;
               }
               unique_flag = true;
           }
        }
    }

    //cout <<num_unique_mol <<" unique atom types that occur in mol: ";
    /*cout << "IN HUNGARIAN MOL unique ";
    for (int i=0; i<num_unique_mol; i++){ cout <<unique_atom_types_mol[i] <<" "; }
    cout <<endl;*/


    // Make a list of atom types that occur in both, and enumerate unmatched atoms
    for (int i=0; i<num_unique_refmol; i++){
        for (int j=0; j<num_unique_mol; j++){
            if (unique_atom_types_refmol[i].compare(unique_atom_types_mol[j])==0){
                unique_atom_types_vector.push_back(unique_atom_types_refmol[i]);
                num_unique++;
            }
        }
    }

    //cout <<num_unique <<" types occur in both: ";
    //for (int i=0; i<num_unique; i++){ cout <<unique_atom_types_vector[i] <<" "; }
    //cout <<endl;


    // Find unmatched atoms in refmol
    bool temp_flag = false;
    for (int i=0; i<num_unique_refmol; i++){
        for (int j=0; j<unique_atom_types_vector.size(); j++){
            if (unique_atom_types_refmol[i].compare(unique_atom_types_vector[j])==0){
                temp_flag = true;
            }
        }
        if (!temp_flag){
            // CS:2016-06-23, Adjusted to count active atoms only
            for (int j=0; j<refmol.num_atoms; j++){
                if (refmol.atom_active_flags[j] == true){
                   if (refmol.atom_types[j].compare(unique_atom_types_refmol[i])==0){
                      //cout <<refmol.atom_types[j] <<" is not found in both!" <<endl;
                      num_unmatched++;
                   }
                }
            }
        }
        temp_flag = false;
    }

    //cout <<"After searching through refmol, there are " <<num_unmatched <<" total unmatched atoms so far" <<endl;


    // Find unmatched atoms in mol
    temp_flag = false;
    for (int i=0; i<num_unique_mol; i++){
        for (int j=0; j<unique_atom_types_vector.size(); j++){
            if (unique_atom_types_mol[i].compare(unique_atom_types_vector[j])==0){
                temp_flag = true;
            }
        }
        if (!temp_flag){
            // CS:2016-06-23, Adjusted to count active atoms only
            for (int j=0; j<mol.num_atoms; j++){
                if (mol.atom_active_flags[j] == true){
                   if (mol.atom_types[j].compare(unique_atom_types_mol[i])==0){
                      //cout <<mol.atom_types[j] <<" is not found in both!" <<endl;
                      num_unmatched++;
                   }
                }
            }
        }
        temp_flag = false;
    }

    //cout <<"After searching through mol and refmol, there are " <<num_unmatched <<" total unmatched atoms so far" <<endl;
    // Delete some stuff
    delete [] unique_atom_types_refmol;
    delete [] unique_atom_types_mol;


    // Go through unique_atom_types one by one, forming a distance matrix for each set of
    // atoms of the same atom type. Use hungarian to solve each matrix, and add the assignment
    // to total_assignment
    for (int z=0; z<num_unique; z++){

        num_type = 0;
        int num_type1 = 0;
        int num_type2 = 0;
        int difference = 0;

        // For every atom in refmol, count how many match the current atom type
        // CS:2016-06-23, Adjusted to count active atoms only
        for (int i=0; i<refmol.num_atoms; i++){
            if (refmol.atom_active_flags[i] == true){ 
                if (refmol.atom_types[i].compare(unique_atom_types_vector[z]) == 0){
                   num_type1++;
                }
            }
        }

        // For every atom in mol, count how many match the current atom type
        // CS:2016-06-23, Adjusted to count active atoms only
        for (int i=0; i<mol.num_atoms; i++){
            if (mol.atom_active_flags[i] == true){ 
               if (mol.atom_types[i].compare(unique_atom_types_vector[z]) == 0){
                  num_type2++;
               }
            }
        }

        //cout <<"Iteration " <<z+1 <<", atomtype " <<unique_atom_types_vector[z] <<", num_type1 = " <<num_type1 <<", num_type2 = " <<num_type2 <<endl;


        // Store the larger of the two as 'num_type' and increment num_unmatched
        if (num_type1 > num_type2){
            num_type = num_type1;
            num_matched = num_type2;
            difference = num_type1-num_type2;
            num_unmatched += difference;
            //cout <<"found " <<difference <<" more unmatched, num_unmatched = " <<num_unmatched <<endl;

        } else if (num_type2 > num_type1){
            num_type = num_type2;
            num_matched = num_type1;
            difference = num_type2-num_type1;
            num_unmatched += difference;
            //cout <<"found " <<difference <<" more unmatched, num_unmatched = " <<num_unmatched <<endl;

        } else {
            num_type = num_type2;
            num_matched = num_type1;
            //cout <<"here is a case where there is no difference" <<endl;
        }


        // Initialize the matrices based on num_type
        //cout <<"num_type when initializing = " <<num_type <<endl;
        initialize();
        NUM = 1000;


        // Populate 'matrix[][]' with distances (needs to be done slightly different depending on
        // whether num_type1 or num_type2 is larger, or if they are the same size)

        // refmol has more of the current atom type than mol
        if (num_type1 > num_type2){
            //cout <<"num_type1 is greater than num_type2" <<endl;
            int count_i = -1;
            int count_j = -1;

            // for each row in the matrix
            // CS:2016-06-23, Adjusted to count active atoms only
            for (int i=0; i<refmol.num_atoms; i++){
                if (refmol.atom_active_flags[i] == true){
                   if (refmol.atom_types[i].compare(unique_atom_types_vector[z]) == 0){
                      count_i++;
                      count_j = -1;

                      // fill in as many distances as possible
                      // CS:2016-06-23, Adjusted to count active atoms only
                      for (int j=0; j<mol.num_atoms; j++){
                          if (mol.atom_active_flags[j] == true){
                             if (mol.atom_types[j].compare(unique_atom_types_vector[z]) == 0){
                                count_j++;
                                matrix[count_i][count_j] = ((refmol.x[i] - mol.x[j])*(refmol.x[i] - mol.x[j])) +
                                                           ((refmol.y[i] - mol.y[j])*(refmol.y[i] - mol.y[j])) +
                                                           ((refmol.z[i] - mol.z[j])*(refmol.z[i] - mol.z[j]));
                             }
                         } 
                      }

                      // then put some dummy values at the end to make it square
                      for (int j=0; j<difference; j++){
                          count_j++;
                          matrix[count_i][count_j] = 100000.0;
                      }
                  }
               }
            }

            //cout <<"Here is the resulting matrix: " <<endl;
            //for (int i=0; i<num_type1; i++){
            //    for (int j=0; j<num_type1; j++){
            //        cout <<matrix[i][j] <<" ";
            //    }
            //    cout <<endl;
            //}

        // mol has more of the current atom type than refmol
        } else if (num_type2 > num_type1){
            //cout <<"num_type2 is greater than num_type1" <<endl;
            int count_i = -1;
            int count_j = -1;

            // for each row in the matrix
            // CS:2016-06-23, Adjusted to count active atoms only
            for (int i=0; i<mol.num_atoms; i++){
                if (mol.atom_active_flags[i] == true){
                   if (mol.atom_types[i].compare(unique_atom_types_vector[z]) == 0){
                      count_i++;
                      count_j = -1;

                      // fill in as many distances as possible
                      // CS:2016-06-23, Adjusted to count active atoms only
                      for (int j=0; j<refmol.num_atoms; j++){
                          if (refmol.atom_active_flags[j] == true){
                             if (refmol.atom_types[j].compare(unique_atom_types_vector[z]) == 0){
                                count_j++;
                                matrix[count_i][count_j] = ((mol.x[i] - refmol.x[j])*(mol.x[i] - refmol.x[j])) +
                                                           ((mol.y[i] - refmol.y[j])*(mol.y[i] - refmol.y[j])) +
                                                           ((mol.z[i] - refmol.z[j])*(mol.z[i] - refmol.z[j]));
                             }
                          }
                      }

                      // then put some dummy vales at the end to make it square
                      for (int j=0; j<difference; j++){
                          count_j++;
                          matrix[count_i][count_j] = 100000.0;
                      }                    
                  }
               }
            }


            //cout <<"Here is the resulting matrix: " <<endl;
            //for (int i=0; i<num_type2; i++){
            //    for (int j=0; j<num_type2; j++){
            //        cout <<matrix[i][j] <<" ";
            //    }
            //    cout <<endl;
            //}

        // mol and refmol have the same number of the current atom type
        } else {
            //cout <<"num_type1 and num_type2 are equal" <<endl;
            int count_i = -1;
            int count_j = -1;

            // CS:2016-06-23, Adjusted to count active atoms only
            for (int i=0; i<mol.num_atoms; i++){
                if (mol.atom_active_flags[i] == true){
                   if (mol.atom_types[i].compare(unique_atom_types_vector[z]) == 0){
                      count_i++;
                      count_j = -1;
                      for (int j=0; j<refmol.num_atoms; j++){
                          if (refmol.atom_active_flags[j] == true){
                             if (refmol.atom_types[j].compare(unique_atom_types_vector[z]) == 0){
                                count_j++;
                                matrix[count_i][count_j] = ((mol.x[i] - refmol.x[j])*(mol.x[i] - refmol.x[j])) +
                                                           ((mol.y[i] - refmol.y[j])*(mol.y[i] - refmol.y[j])) +
                                                           ((mol.z[i] - refmol.z[j])*(mol.z[i] - refmol.z[j]));
                             }
                          } 
                      }
                  }
               } 
            } 

            //cout <<"Here is the resulting matrix: " <<endl;
            //for (int i=0; i<num_type; i++){
            //    for (int j=0; j<num_type; j++){
            //        cout <<matrix[i][j] <<" ";
            //    }
            //    cout <<endl;
            //}

        }
        //cout <<endl;


        // Solve matrix[][] with hungarian algorithm, store results in matrix_match[]
        hungarian();

        //cout <<"(refmol: mol)" <<endl;
        //for (int i=0; i<num_type; i++){
        //    cout <<i <<" : " <<matrix_match[i] <<" - " <<matrix_original[i][matrix_match[i]] <<endl;
        //}
        //cout <<endl;


        // Compute the cost of matrix_match[] for the smallest num_matched matches
        double assignment = sum_assignment_dissimilar(num_matched);
        num_matched_total+=num_matched;
        total_assignment += assignment;


        // De-allocate memory of all matrices given the size of num_type
        clear();

    } // end for every atom_type


    // Compute the rmsd for the matched atoms
    // If there are matched atoms, else return very large number - CS: 09/19/2016
    if ( num_matched_total > 0 ){rmsd = sqrt( (total_assignment/num_matched_total) );}
    else { rmsd = 999999999; }
    
    //cout <<"RMSD for matched atoms = " <<rmsd <<endl;
    //cout <<"Total number of matched atoms = " <<num_matched_total <<endl;
    //cout <<"Total number of unmatched atoms = " <<num_unmatched <<endl <<endl;


    // Combine the results as a pair and return them
    unique_atom_types_vector.clear();
    pair <double, int> result;
    result.first = rmsd;
    result.second = num_unmatched;
    return result;

} // end Hungarian_RMSD::calc_Hungarian_RMSD_dissimilar() 



// +++++++++++++++++++++++++++++++++++++++++
// This function allocates the memory required for arrays
//
void Hungarian_RMSD::initialize(){
   
    // Citation for Hungarian Score
    //cout << "To cite Hungarian RMSD use: \n Allen, W. J.; Rizzo, R. C. Implementation of the Hungarian Algorithm to Account for Ligand Symmetry and Similarity in Structure-Based Design, J. Chem. Inf. Model., 2014, 54, 518-529 \n" << endl;

    // Allocate the memory for the 2D arrays 
    matrix           =  new double* [num_type];		// cost matrix
    matrix_original  =  new double* [num_type];		// backup copy of cost matrix
    matrix_case      =  new int* [num_type];		// matrix of flags (0 || 1 || 2)

    // Allocate the rest of the memory for the 2D arrays
    for (int i=0; i<num_type; i++){
        matrix[i]           =  new double [num_type];
        matrix_original[i]  =  new double [num_type];
        matrix_case[i]      =  new int [num_type];
    }

    // Allocate the memory for the 1D arrays
    row_count        =  new int [num_type];		// number of 0s in each row
    column_count     =  new int [num_type];		// number of 0s in each column
    row_assigned     =  new int [num_type];		// flag for row state (0 || 1)
    column_assigned  =  new int [num_type];		// flag for column state (0 || 1)
    matrix_match     =  new int [num_type];		// worker / job assignments
}



// +++++++++++++++++++++++++++++++++++++++++
// This is a step-wise implementation of the O(n^4) solution of the Hungarian algorithm to solve
// the minimum assignment problem
//
void Hungarian_RMSD::hungarian(){

    // Step 1: Prepare the matrix by making a copy and reducing the cells.
    //
    reduce_matrix();

    // This while loop will end when every worker is assigned a job
    int count_cycles = 0;
    while (count_cycles < (NUM*NUM)){

        // Step 2: Attempt to assign one job to each worker - if number of assignments is less
        // than number of workers, go to Step 3, else skip to end.
        //
        int number_assigned = assign_jobs();
        if (number_assigned == num_type){ break; }

        // Step 3: Draw the minimum number of lines required to cross out all 0s from matrix[][].
        //
        draw_line();

        // Step 4: Update the cells of matrix[][] according to how lines were drawn.
        //
        update_matrix();
	count_cycles++;

    } // end loop when number of assignments is equal to number of workers/jobs

    //cout <<"  cycles = " <<count_cycles+1 <<endl;

    return;

} // end Hungarian_RMSD::hungarian()



// +++++++++++++++++++++++++++++++++++++++++
// Each of the next 6 functions ( reset_xxx() ) resets all of the values of a certain matrix
//
void Hungarian_RMSD::reset_match(){
    for (int i=0; i<num_type; i++){
        matrix_match[i] = -100;
    }
}
void Hungarian_RMSD::reset_row_assigned(){
    for (int i=0; i<num_type; i++){
        row_assigned[i] = 0;
    }
}
void Hungarian_RMSD::reset_column_assigned(){
    for (int i=0; i<num_type; i++){
        column_assigned[i] = 0;
    }
}
void Hungarian_RMSD::reset_row_count(){
    for (int i=0; i<num_type; i++){
        row_count[i] = 0;
    }
}
void Hungarian_RMSD::reset_column_count(){
    for (int i=0; i<num_type; i++){
        column_count[i] = 0;
    }
}
void Hungarian_RMSD::reset_case(){
    for (int i=0; i<num_type; i++){
        for (int j=0; j<num_type; j++){
            matrix_case[i][j]=0;
        }
    }
}



// +++++++++++++++++++++++++++++++++++++++++
// This function has three tasks. First, it backs up matrix[] as matrix_original[][]. Then it 
// finds the minimum value in each row and subtracts that from every value in the row. Finally, 
// it does that same thing, but for columns
//
void Hungarian_RMSD::reduce_matrix(){

    // 1. Read through the original matrix and make a backup copy
    //
    for (int i=0; i<num_type; i++){
        for (int j=0; j<num_type; j++){
            matrix_original[i][j] = matrix[i][j];
        }
    }

    // 2. find the minimum value in each row and subtract that from every other value in the row
    //
    for (int i=0; i<num_type; i++){
        double min_value = MAX;
        for (int j=0; j<num_type; j++){			// look through each column, j, in row i, and 
            if (matrix[i][j] < min_value){		// remember which is the minimum value
                min_value = matrix[i][j];
            }
        }
        for (int j=0; j<num_type; j++){			// then go back through the same row and subtract
            matrix[i][j] -= min_value;			// the minimum value from each cell
        }
    }

    // 3. find the minimum value in each column and subtract that from every other value in the column
    //
    for (int j=0; j<num_type; j++){
        double min_value = MAX;
        for (int i=0; i<num_type; i++){			// look through each row, i, in column j, and 
            if (matrix[i][j] < min_value){		// remember which is the minimum value
                min_value = matrix[i][j];
            }
        }
        for (int i=0; i<num_type; i++){			// then go back through the same column and 
            matrix[i][j] -= min_value;			// subtract the minimum value from each cell
        } 
    }

    return;

} // end Hungarian_RMSD::reduce_matrix()



// +++++++++++++++++++++++++++++++++++++++++
// In this function assignments are first made in any row or column that only contains one
// 0. All other 0s (or any number, for that matter) in the row and in the column containing
// assignment are 'crossed out', and the row and column are remembered in row_assigned[] and
// column_assigned[]. In the case that each remaining row and column has two (or more) 0s,
// then an assignment is made arbitrarily in any row or column with the minimum number of 0s.
//
int Hungarian_RMSD::assign_jobs(){

    // Reset the matrices that will be used in this loop and set number_assigned to 0
    reset_match();
    reset_row_assigned();
    reset_column_assigned();
    int number_assigned = 0;	// the number of worker/job pairs that have been assigned

    // This while loop will proceed until all possible worker/job assignments are made
    while (number_assigned  < num_type+1){

        // 0. Reset the row_count[] and column_count[] counters at the beginning of each round
        //
        reset_row_count();
        reset_column_count();

        // 1. Count the number of 0s in each row and column
        //
        for (int i=0; i<num_type; i++){ 		// The result will be two 1-dimensional arrays that contain
            for (int j=0; j<num_type; j++){		// the number of 0s in each row or column (these counts are
                if (matrix[i][j] == 0){			// reset at each new iteration of assign_jobs())
                    if (row_assigned[i] == 0 && column_assigned[j] == 0){	// and 0s are only counted if that row
                        row_count[i]++;						// or column has not been crossed out
                        column_count[j]++;
                    }
                }
            }
        }

        // 2. Identify the row or column that has the least number of zeroes in it
        //
        double min_value = MAX;		// a baseline value for the upcoming comparisons
        int row_or_column = 0;		// can be 0 (no more 0s left in the matrix), 1 (the least number
					// of 0s is found in a row), or 2 (the least number of 0s is
					// found in a column)
        int row = -1;			// the row with the least number of 0s
        int column = -1;		// the column with the least number of 0s

        // The next loop looks through row_count[] and column_count[] to find the row or column that contains the
        // least number of 0s. However, it will only consider that row or column so long as row_assigned[] or
        // column_assigned[] is 0 in that same position. (if, for example, row_assigned[x] is 1, then an assignment
        // has already been made in row x, so row x should be ignored). The flag row_or_column is set to remember
        // whether a row or column contains the least number of 0s, and the variables row OR column are assigned to
        // remember in which row or column the least number of zeros are located.
        for (int i=0; i<num_type; i++){
            if (row_count[i] != 0 && row_count[i] < min_value){
                min_value = row_count[i];
                row = i;
                row_or_column = 1;
                column = -1;
            }
            if (column_count[i] != 0 && column_count[i] < min_value){
                min_value = column_count[i];
                column = i;
                row_or_column = 2;
                row = -1;
            }
        }

        // 3. Make an assignment in the row or column with the least number of 0s	
        //
        if (row_or_column == 0){		// if no row or column contains any more 0s, assignment is
            break;				// complete and assign_jobs() is exited

        } else if (row_or_column == 1){					// if a row has the least number of 0s,
            for (int i=0; i<num_type; i++){				// then check through each column [i] of that row and
                if (matrix[row][i] == 0 && column_assigned[i] == 0){	// make an assignment, so long as now other assignment
                    matrix_match[row] = i;				// has already been made in that column
                    row_assigned[row] = 1;
                    column_assigned[i] = 1;
                    break;			// this breaks the for loop
                }
            }
        } else if (row_or_column == 2){					// do the same thing, but if a column has the least
            for (int i=0; i<num_type; i++){				// number of 0s
                if (matrix[i][column] == 0 && row_assigned[i] == 0){
                    matrix_match[i] = column;
                    column_assigned[column] = 1;
                    row_assigned[i] = 1;
                    break;			// this breaks the for loop
                }
            }
        }

        // 4. Increment the counter and repeat
        //
        number_assigned++;

    } // end while: this will loop until all possible worker / job pairs are formed

    return(number_assigned);

} // end Hungarian_RMSD::assign_jobs()



// +++++++++++++++++++++++++++++++++++++++++
// This function figures out the least number of lines required to cross out all of the
// 0s in matrix[][]. First, all rows with no assignment are marked. Second, all columns
// with 0s in marked rows are marked. Third, all rows with assignments in marked columns
// are marked. That is repeated from step 2 until loop is closed. Lines are drawn over 
// marked columns and unmarked rows.
//
void Hungarian_RMSD::draw_line(){

    // 0. Reset the matrices that will be used in this function
    //
    reset_row_assigned();
    reset_column_assigned();
    reset_case();
    int lines = 0;	// the number of lines required to cross out all 0s in matrix[][]

    // 1. Mark all rows that have no assignment
    //
    for (int i=0; i<num_type; i++){
        if (matrix_match[i] == -100){		// if matrix_match[i] is still -100, then no assignment was ever made in row i
            row_assigned[i] = 1;		// therefore, row_assigned[i] is changed from 0 to 1
        }
    }

    bool flag_update = true;		// a flag that remembers whether an update was made during the last iteration
					// if no update was made, then the loop is 'closed'
    int count_update = 0;		// just to make sure it doesn't loop indefinitely

    // This while loop will continue until it figures out which lines to draw
    while (flag_update && count_update < num_type*2){

        flag_update = false;		// set the flag to false initially, it will be set to true if an update is made

        // 2. Mark all columns that have 0s in those rows
        //
        for (int i=0; i<num_type; i++){			// for all of the rows in matrix[][],
            if (row_assigned[i] == 1){			// if an assignement was made in that row,
                for (int j=0; j<num_type; j++){		// then go through each column of matrix[row][]
                    if (matrix[i][j] == 0){		// and check to see if there is 0
                        column_assigned[j] = 1;		// and remember that position in column_assigned
                    }
                }
            }
        }

        // 3. Mark all rows having assignments in those columns
        //
        for (int i=0; i<num_type; i++){					// for all of the columns in matrix[][]
            if (column_assigned[i] == 1){				// if the column has been 'marked',
                for (int j=0; j<num_type; j++){				// then check every row in matrix[][column]
                    if (matrix_match[j] == i && row_assigned[j] != 1){	// and see if an assignment has been made - if so,
                        row_assigned[j] = 1;				// then mark that row
                        flag_update = true;				// an update has been made, so repeat
                    }
                }
            }
        }

    } // end while loop when no more updates are made

    // 4. Draw lines through marked columns and unmarked rows
    //
    for (int i=0; i<num_type; i++){			// for all of the cells in matrix[][]...
        if (row_assigned[i] == 0){			// if a row has not been marked,
            for (int j=0; j<num_type; j++){		// draw a line through that row by incrementing every entry in
                matrix_case[i][j]++;			// matrix_case[row][] by 1
            }
            lines++;					// a line was just drawn, increment lines by 1
        }
        if (column_assigned[i] == 1){			// if a column has been marked,
            for (int j=0; j<num_type; j++){		// draw a line through that column by incrementing every entry in
                matrix_case[j][i]++;			// matrix_case[][column] by 1
            }
            lines++;					// a line was just drawn, increment lines by 1
        }
    }
    // every cell in matrix_case[][] will now either be a 0 (no lines covering), 1 (one line covering),
    // or 2 (two lines covering)

    return;

} // end Hungarian_RMSD::draw_line()



// +++++++++++++++++++++++++++++++++++++++++
// This function updates the cells in matrix[][] according to the state of matrix_case[][]. First it
// finds the minimum value of all cells in matrix[i][j] where matrix_case[i][j] == 0. Then, it goes
// back through matrix[][], subtracting the minimum value from each of those cells. It then adds that
// value to each of the cells in matrix[i][j] where matrix_case[i][j] == 2. Any cells in matrix[i][j]
// where matrix_case[i][j] == 1 are ignored.
//
void Hungarian_RMSD::update_matrix(){

    // 1. Find the minimum value in all cells of matrix[][] where the corresponding cell in
    // matrix_case[][] is 0.
    //
    double min_value = MAX;				// a baseline for finding minimum values

    for (int i=0; i<num_type; i++){			// for all cells in matrix_case[][],
        for (int j=0; j<num_type; j++){
            if (matrix_case[i][j] == 0){		// if it is '0', then consider the value in the same 
                if (matrix[i][j] < min_value){		// position in matrix[][] for the minimum value of all '0s'
                    min_value = matrix[i][j];
                }
            }
        }
    }

    // 2. Update all the cells of matrix[][] according to the rules below
    //
    for (int i=0; i<num_type; i++){			// for every cell in matrix[][], 
        for (int j=0; j<num_type; j++){
            if (matrix_case[i][j] == 0){		// if the same position in matrix_case[][] is 0, then subtract
                matrix[i][j] -= min_value;		// the min_value from step 1
            } else if (matrix_case[i][j] == 1){		// if the same position in matrix_case[][] is 1, then do nothing
                // do nothing
            } else if (matrix_case[i][j] == 2){		// if the same position in matrix_casep[][] is 2, then add the
                matrix[i][j] += min_value;		// min_value from step 1
            }
        }
    }

    return;

} // end Hungarian_RMSD::update_matrix()



// +++++++++++++++++++++++++++++++++++++++++
// Once all of the worker / job assignments are made, this function looks at which worker is 
// assigned to which job, grabs the original distance from matrix_original[][], and sums across 
// all pairs.
//
double Hungarian_RMSD::sum_assignment(){
    double running_total = 0;		// a running total for finding the sum of individual matrices

    // The sum of the original 'costs' in matrix_original[][] based on the final matching in matrix_match[]
    for (int i=0; i<num_type; i++){
            running_total += matrix_original[i][matrix_match[i]];
    }

    return(running_total);

} // end Hungarian_RMSD::sum_assignment()



// +++++++++++++++++++++++++++++++++++++++++
// Once all of the worker / job assignments are made, this function looks at which worker is 
// assigned to which job, grabs the original distance from matrix_original[][], and sums across 
// all pairs.
//
double Hungarian_RMSD::sum_assignment_dissimilar(int num_matched){
    double running_total = 0;
    vector <double> temp_vector;

    // First push all of the values onto this vector
    for (int i=0; i<num_type; i++){
        temp_vector.push_back(matrix_original[i][matrix_match[i]]);
    }

    // Sort the vector
    sort(temp_vector.begin(), temp_vector.end());

    // Sum the first 'num_matched' values on the vector
    for (int i=0; i<num_matched; i++){
        running_total += temp_vector[i];
    }

    temp_vector.clear();
    return(running_total);

} // end Hungarian_RMSD::sum_assignment_dissimilar()



// +++++++++++++++++++++++++++++++++++++++++
// De-allocate the memory of each array at the end of each round of atom types
//
void Hungarian_RMSD::clear(){

    for (int i=0; i<num_type; i++){
        delete[] matrix[i];
        delete[] matrix_original[i];
        delete[] matrix_case[i];
    }

    delete[] matrix;
    delete[] matrix_original;
    delete[] matrix_case;

    delete[] row_count;
    delete[] column_count;	
    delete[] row_assigned;
    delete[] column_assigned;	
    delete[] matrix_match;
}



