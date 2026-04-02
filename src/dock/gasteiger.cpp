#include "gasteiger.h"
#include <math.h>
using namespace std;

class DOCKMol;


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

bool
rms_converg(double* tot_charges, double* prev_charges, int num_atoms) {


//
    double converg = 0.00001;
//

    double rms = 0.0;
    for (int i = 0; i < num_atoms; i++) {
        rms += ( (tot_charges[i] - prev_charges[i]) * (tot_charges[i] - prev_charges[i]) );
    }

    rms = rms / num_atoms;
    //rms = pow(rms, 0.5);
    rms = sqrt(rms);


    if (rms <= converg) {return true;}
    return false;
} // rmsConverg



///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


float
oco2charge(DOCKMol & mol, int a) {
    //printf("entering oco2charge() function, a = %d \n", a);
    int i = 0;
    while (mol.atom_envs[a][i] != '#') { //read fingerprint
        //printf("loop interation = %d \n", i); 
        // phosphate groups
        if (mol.atom_envs[a][i] == 'T') { //connected to P.3
            for (int j = 0; j < 6; j++) { //foreach atom connected to O.co2
                if (mol.atom_envs[mol.neighbor_list[a][j]][1] == 'T') {  //check for index of phosphate
                    int k = 0;
                    int p = mol.neighbor_list[a][j];  // index of phos
                    int n = 0;
                    while (mol.atom_envs[p][k] != '#') { //count up how many o.co2s are bound
                        if (mol.atom_envs[p][k] == 'O') {n++;}
                        k++;                
                    } //count up O.c02
   
                  /*if (n == 3) return -0.6666667;
                    if (n == 2) return -0.5;
                    return -0.0; */
                    return ( (1.0 - n) / n);

                } //index of phos
            } //get index of phos
        } // if bound to phosphorous   

        // add (-) charge to oxygens on sulfoxides and sulfones
        if (mol.atom_envs[a][i] == 'P' || //S.3
            mol.atom_envs[a][i] == 'R' || //S.o
            mol.atom_envs[a][i] == 'S' ) { //S.o2     
                return -1.0;
            } //sulfoxide or sulfone



    // check if protonated acid
        if (mol.atom_envs[a][i] == 'B') { //connected to C.2, check for protonation
            //printf("gasteiger O.co2 connected to C.2 \n");
            int c = -1;
            //printf("O.co2 neighbor_list size %d \n", mol.neighbor_list[a].size());
            for (int j = 0; j < mol.neighbor_list[a].size(); j++) { //foreach atom connected to O.co2

                //printf("neighbor_list iter %d \n", j);
                int index = mol.neighbor_list[a][j];
                if (mol.atom_envs[index][1] == 'Z') {  // protonated? (Z is hydrogen)
                    return 0.0; // if protonated, no charge 
                }
//   BCF June 26, 2013 also need to change parameters to C.3
//
//
                if (mol.atom_envs[mol.neighbor_list[a][j]][1] == 'B') {  // get index of C.2
                    c = mol.neighbor_list[a][j];  // index of C.2
                }

                if (mol.atom_envs[mol.neighbor_list[a][j]][1] == '#') {  // end of fingerprint, should have C.2 index
                    break;
                }
            } //foreach atom connected to O.co2  
            for (int j = 0; j < mol.neighbor_list[c].size(); j++) { //get index of attached O.co2
                if (mol.neighbor_list[c][j] != a) { //don't consider O.co2 we already looked at
                    int cur_atom = mol.neighbor_list[c][j];
                    if (mol.atom_envs[cur_atom][1] == 'O') { //if other O.co2
                        int k = 0;
                        while (mol.atom_envs[cur_atom][k] != '#') { // read env of other O.co2 to check protonation
                            if (mol.atom_envs[cur_atom][k] == 'Z') { // if other O.co2 is protonated
                                return 0.0;
                            } // if protonated
                            k++;
                        } // read env of other O.co2
                    } //if other O.co2                
                } // if not O.co2 we already considered
            } // foreach atom connected to central C.2


        } //if connected to C.2



        i++;  //move onto next connected atom
    } // read original O.co2 fingerprint
    return -0.5;
} //oco2charge



//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////




float
o2charge(DOCKMol & mol, int a) {
    int i = 0;
    while (mol.atom_envs[a][i] != '#') { //read fingerprint
    

    // check if nitro
        if (mol.atom_envs[a][i] == 'G') { //connected to N.2, check nitro 
            int n;
            for (int j = 0; j < mol.neighbor_list[a].size(); j++) { //foreach atom connected to N.2

                int index = mol.neighbor_list[a][j];
                if (mol.atom_envs[index][1] == 'G') {  // get index of N.2 if connected
                    n = mol.neighbor_list[a][j];  
                    break; 
                }

                if (mol.atom_envs[mol.neighbor_list[a][j]][1] == '#') {  // end of fingerprint, not nitro, is neutral
                    return 0.0;
                }
            } //foreach atom connected to O.co2  
            for (int j = 0; j < mol.neighbor_list[n].size(); j++) { //check N.2 connection list for O.co2
                if (mol.neighbor_list[n][j] != a) { //don't consider atom we already looked at
                    int cur_atom = mol.neighbor_list[n][j];
                    if (mol.atom_envs[cur_atom][1] == 'O') { //O.co2
                        int k = 0;
                        while (mol.atom_envs[cur_atom][k] != '#') { // read env of other O.co2 to check protonation
                            if (mol.atom_envs[cur_atom][k] == 'Z') { // if other O.co2 is protonated
                                return 0.0;
                            } // if protonated
                            k++;
                        } // read env of other O.co2
                        return -0.5;
                    } //if other O.co2                
                } // if not atom we already considered
            } // foreach atom connected to central N.2


        } //if connected to N.2


    // add (-) charge to oxygens on sulfoxides and sulfones, same code as for O.co2 (sometimes these are mixed)
        if (mol.atom_envs[a][i] == 'P' || //S.3
            mol.atom_envs[a][i] == 'R' || //S.o
            mol.atom_envs[a][i] == 'S' ) { //S.o2
                return -1.0;
            } //sulfoxide or sulfone


        i++;  //move onto next connected atom
    } // read original O.2 fingerprint
    return 0.0;
} //o2charge


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////



float n2charge(DOCKMol & mol, int a) {
    int i = 0;
    int oxy_bound = 0;  //check if nitro
    int total_bound = -1; //in joeprint there will be 1 more ^ than atom bonded if atom is bonded at all
    while (mol.atom_envs[a][i] != '#') { //read fingerprint
        if (mol.atom_envs[a][i] == 'O' || mol.atom_envs[a][i] == 'N') { //O.2 or O.co2
            oxy_bound++;
        } // if bound to O.2 or ).co2
        if (mol.atom_envs[a][i] == '^') { //in joeprint ^ separates elements, there are 2 at end
        // 1 ^ -> not possible
        // 2 ^ -> nothing bonded, 1 atom bonded
        // 3 ^ -> 2 atoms bonded
        // 4 ^ -> 3 atoms bonded
            total_bound++;
        } 
        i++;
    } // read N.2 joeprint
    if (oxy_bound == 2 && total_bound == 3) { //nitro group
        return 1.0;
    }
    return 0.0;
} //n2charge


/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////


// need to adjust charge on sulfers (and oxygens) in sulfoxides and sulfones to get proper polarity
// modified to account for sulfonium ions
float sulfcharge(DOCKMol & mol, int a) {
    int i = 0;
    int oxy_bound = 0;  //number of O.2 and O.co2 bound to sulfer
    int num_singlebond = 0;
    bool start_counting = false;
    while (mol.atom_envs[a][i] != '#') { //read fingerprint
        if (mol.atom_envs[a][i] == 'O' || 
            mol.atom_envs[a][i] == 'N') { //O.2 or O.co2
            
            oxy_bound++;
        } // if bound to O.2 or ).co2
        if (mol.atom_envs[a][i] == '1') num_singlebond++;

        i++;
    } // read joeprint
 /*   switch (oxy_bound) {
        case 1:  //sulfoxide initial charge +1
            return 1.0;
        case 2:  //sulfone initial charge +2
        case 3:  //sulfonic acid initial charge +2
            return 2.0;
        default:
            return 0.0;
    } //switch on number of oxygen bound */
    if      (oxy_bound == 1) return 1.0;
    else if (oxy_bound == 2) return 2.0;
    else if (oxy_bound == 3) return 2.0; //sulfonic acid -SO3
    else if (oxy_bound == 4) return 2.0; //sulfuric acid SO4
    else if ( (mol.atom_envs[a][1] == 'P') && (num_singlebond == 3) ) return 1.0; //S.3 with 3 single bonds - sulfonium
    else                     return 0.0;
} //sulfcharge


/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////


float
compute_gast_charges(DOCKMol & mol) {
    //double gast_params[2][4] = { {7.98, 9.18, 1.88, 0.0}, {7.17, 6.24, -0.56, 0.0} };
    double q = 0.0; // charge transfered

    //printf("length of arrays %d \n", mol.num_atoms);
    double tot_charges[mol.num_atoms];
    double prev_charges[mol.num_atoms]; //stores last charge value for convergence
    double electroneg[mol.num_atoms];
    double param_a[mol.num_atoms];
    double param_b[mol.num_atoms];
    double param_c[mol.num_atoms];
    double param_d[mol.num_atoms];


    //Fochmod
    //printf("gasteiger mark #1 \n");

    for (int a = 0; a < mol.num_atoms; a++) {
        //printf("switchloop %d  \n", a);
        switch (mol.atom_envs[a][1]) {
            case 'A':   //C.3  tetetete
                tot_charges[a] = 0.0;
                electroneg[a]  = 7.98;
                param_a[a]     = 7.98;
                param_b[a]     = 9.18;
                param_c[a]     = 1.88;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'B':   //C.2  trtrtr pi
                tot_charges[a] = 0.0;
                electroneg[a]  = 8.79;
                param_a[a]     = 8.79;
                param_b[a]     = 9.32;
                param_c[a]     = 1.51;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'C':   //C.1  didi pi pi
                tot_charges[a] = 0.0;
                electroneg[a]  = 10.39;
                param_a[a]     = 10.39;
                param_b[a]     = 9.45;
                param_c[a]     = 0.73;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'D':   //C.ar, C.2 params 
                tot_charges[a] = 0.0;
                electroneg[a]  = 8.79;
                param_a[a]     = 8.79;
                param_b[a]     = 9.32;
                param_c[a]     = 1.51;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'E':   //C.cat  guanadinium, C.3 params with +1 formal charge
                tot_charges[a] = 1.0;
                electroneg[a]  = 7.98;
                param_a[a]     = 7.98;
                param_b[a]     = 9.18;
                param_c[a]     = 1.88;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'F':   //N.3  te^2 tetete
                tot_charges[a] = 0.0;
                electroneg[a]  = 11.54;
                param_a[a]     = 11.54;
                param_b[a]     = 10.82;
                param_c[a]     = 1.36;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'G':   //N.2  tr^2 trtr pi
                tot_charges[a] = n2charge(mol, a);;
                electroneg[a]  = 12.87;
                param_a[a]     = 12.87;
                param_b[a]     = 11.15;
                param_c[a]     = 0.85;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'H':   //N.1  di^2 di pi pi
                tot_charges[a] = 0.0;
                electroneg[a]  = 15.68;
                param_a[a]     = 15.68;
                param_b[a]     = 11.70;
                param_c[a]     = -0.27;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'I':   //N.ar N.2 params
                tot_charges[a] = 0.0;
                electroneg[a]  = 12.87;
                param_a[a]     = 12.87;
                param_b[a]     = 11.15;
                param_c[a]     = 0.85;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'J':   //N.am N.2 params
                tot_charges[a] = 0.0;
                electroneg[a]  = 12.87;
                param_a[a]     = 12.87;
                param_b[a]     = 11.15;
                param_c[a]     = 0.85;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'K':   //N.pl3  trigonal planar, used sp2 params
                tot_charges[a] = 0.0;
                electroneg[a]  = 12.87;
                param_a[a]     = 12.87;
                param_b[a]     = 11.15;
                param_c[a]     = 0.85;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'L':   //N.4, using na+ from GASPARM.DAT
//GASPARM na+      12.32   11.20    1.34   24.86    1.00
                tot_charges[a] = 1.0;
                electroneg[a]  = 12.32;
                param_a[a]     = 12.32;
                param_b[a]     = 11.20;
                param_c[a]     = 1.34;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'M':   //O.3  te^2 te^2 tete
                tot_charges[a] = 0.0;
                electroneg[a]  = 14.18;
                param_a[a]     = 14.18;
                param_b[a]     = 12.92;
                param_c[a]     = 1.39;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'N':   //O.2  tr^2 trtr pi
                tot_charges[a] = o2charge(mol, a);
                electroneg[a]  = 17.07;
                param_a[a]     = 17.07;
                param_b[a]     = 13.79;
                param_c[a]     = 0.47;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'O':   //O.co2 O.2 with -0.5 formal charge
                tot_charges[a] = oco2charge(mol, a);
                electroneg[a]  = 17.07;
                param_a[a]     = 17.07;
                param_b[a]     = 13.79;
                param_c[a]     = 0.47;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'P':   //S.3  s3 from GASPARM.DAT
                tot_charges[a] = sulfcharge(mol, a);
                electroneg[a]  = 10.14;
                param_a[a]     = 10.14;
                param_b[a]     = 9.13;  
                param_c[a]     = 1.38;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'Q':   //S.2  s2 from GASPARM.DAT
                tot_charges[a] = 0.0;
                electroneg[a]  = 10.88;
                param_a[a]     = 10.88;
                param_b[a]     = 9.485; 
                param_c[a]     = 1.325;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'R':   //S.O, S.o  copied s2 params
                tot_charges[a] = sulfcharge(mol, a);
                electroneg[a]  = 10.88;
                param_a[a]     = 10.88;
                param_b[a]     = 9.485; 
                param_c[a]     = 1.325;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'S':   //S.O2, S.o2 copied s2
                tot_charges[a] = sulfcharge(mol, a);
                electroneg[a]  = 10.88;
                param_a[a]     = 10.88;
                param_b[a]     = 9.485; 
                param_c[a]     = 1.325;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'T':   //P.3  p from GASPARM.DAT
                tot_charges[a] = 0.0;
                electroneg[a]  = 8.90;
                param_a[a]     = 8.90;
                param_b[a]     = 8.24;
                param_c[a]     = 0.96;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'U':   //F  s^2 p^2 p^2 p^2
                tot_charges[a] = 0.0;
                electroneg[a]  = 14.66;
                param_a[a]     = 14.66;
                param_b[a]     = 13.85;
                param_c[a]     = 2.31;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'V':   //Cl  s^2 p^2 p^2 p^2
                tot_charges[a] = 0.0;
                electroneg[a]  = 11.00;
                param_a[a]     = 11.00;
                param_b[a]     = 9.69;
                param_c[a]     = 1.35;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'W':   //Br  s^2 p^2 p^2 p^2
                tot_charges[a] = 0.0;
                electroneg[a]  = 10.08;
                param_a[a]     = 10.08;
                param_b[a]     = 8.47;
                param_c[a]     = 1.16;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'X':   //I  s^2 p^2 p^2 p^2
                tot_charges[a] = 0.0;
                electroneg[a]  = 9.90;
                param_a[a]     = 9.90;
                param_b[a]     = 7.96;
                param_c[a]     = 0.96;
                param_d[a]     = param_a[a] + param_b[a] + param_c[a];
                break;
            case 'Y':   //Du  will replicate Hydrogen parameters
                tot_charges[a] = 0.0;
                electroneg[a]  = 7.17;
                param_a[a]     = 7.17;
                param_b[a]     = 6.24;
                param_c[a]     = -0.56;
                param_d[a]     = 20.02;
                break;
            case 'Z':   //H
                tot_charges[a] = 0.0;
                electroneg[a]  = 7.17;
                param_a[a]     = 7.17;
                param_b[a]     = 6.24;
                param_c[a]     = -0.56;
                param_d[a]     = 20.02;
                break;
            //Fochmod to find segfault, default to H - this is not correct
         /*  default:
                printf("had unknown fingerprint");
                tot_charges[a] = 0.0;
                electroneg[a]  = 7.17;
                param_a[a]     = 7.17;
                param_b[a]     = 6.24;
                param_c[a]     = -0.56;
                param_d[a]     = 20.02;
                break;
         */
        }
    }

    // Count the total number of charges for return
    float total_charges = 0.0;
    for (int a = 0; a < mol.num_atoms; a++) {
        if (tot_charges[a] < 0.0){
            total_charges += tot_charges[a]*(-1);
        } else if (tot_charges[a] > 0.0){
            total_charges += tot_charges[a];
        }
    }

    //Fochmod
    //printf("gasteiger mark #2 \n");

int iteration = 0;
# define GASMAXITER 10
# define DAMPFACTOR 0.5
 for (int i = 0; i < mol.num_atoms; i++) {
     prev_charges[i] = tot_charges[i];
 }


 do {
                for (int i = 0; i < mol.num_atoms; i++) {
                        electroneg[i] =  param_a[i] + 
                                        (param_b[i] * tot_charges[i]) + 
                                        (param_c[i] * tot_charges[i] * tot_charges[i]);
                        if (electroneg[i] == 0.0)
                                electroneg[i] = 0.0000000001;
                }


                for (int i = 0; i < mol.num_atoms; i++) {
                        for (int j = i + 1; j < mol.num_atoms; j++) {  //force j > i
                                    for (int a = 0; a < mol.neighbor_list[i].size(); a++) {
                                        int index = mol.neighbor_list[i][a];
                                        if ( j == index) {
	                                        if (electroneg[i] <= electroneg[j]) {
        	                                        q = (electroneg[j] - electroneg[i]) / param_d[i] * pow(DAMPFACTOR,
                                                                                                         iteration + 1);
                	                                tot_charges[i] += q;
                        	                        tot_charges[j] -= q;
                                	        }
                                        	if (electroneg[i] > electroneg[j]) {
                                               	 	q = (electroneg[i] - electroneg[j]) / param_d[j] * pow(DAMPFACTOR,
                                                                                                         iteration + 1);
                                                	tot_charges[i] -= q;
                                                	tot_charges[j] += q;
                                        	}
                                        } //if j is on that list
                                    } //for each atom connected to i
                        } // for each j
                } //for each i
                iteration++;
        } while ( ( ! rms_converg(tot_charges, prev_charges, mol.num_atoms) ) && iteration < GASMAXITER);








/*

    int iter = 1;
    bool converged = false;
    //int num_converg = 0;
    //while (iter < max_iter && num_converg < mol.num_atoms) {
    while (iter < max_iter && ! converged) {

        //cout << "iter: " << iter << endl;
      /*  for (int a = 0; a < mol.num_atoms; a++) { //foreach atom calc electronegativity
            electroneg[a] = param_a[a] + 
                            (param_b[a] * tot_charges[a]) + 
                            (param_c[a] * tot_charges[a] * tot_charges[a]);
                    
        //cout << mol.atom_types[a] << " has electroneg " << electroneg[a] << endl;
        }
                    */



/*
        //num_converg = 0;
        for (int a = 0; a < mol.num_atoms; a++) { //foreach atom compute charge transfer
            q = 0.0;
        
            for (int b = 0; b < mol.neighbor_list[a].size(); b++) { //foreach connected atom
                int index = mol.neighbor_list[a][b];
                if (electroneg[a] > electroneg[index]) {
                    q +=  (1 / param_d[index]) * 
                          (electroneg[index] - electroneg[a]); 
                }
                else if (electroneg[a] < electroneg[index]) {
                    q +=  (1 / param_d[a]) * 
                          (electroneg[index] - electroneg[a]); 
                }
            } //foreach connected atom

            //cout << q << " charge transfered" << endl;
            prev_charges[a] = tot_charges[a];
            tot_charges[a] += (q * (pow(0.5, iter)));
          //  if (q < converg) { num_converg++;}
        } //foreach atom

        for (int a = 0; a < mol.num_atoms; a++) { //foreach atom calc electronegativity
            electroneg[a] = param_a[a] + 
                            (param_b[a] * tot_charges[a]) + 
                            (param_c[a] * tot_charges[a] * tot_charges[a]);
        }

        converged = rms_converg( tot_charges, prev_charges, mol.num_atoms );


    iter++;
        
    } // while


*/

    //Fochmod
    //printf("gasteiger mark #3 \n");

    for (int a = 0; a < mol.num_atoms; a++) { //foreach atom update mol object
        mol.charges[a] = tot_charges[a];
    }
 
    //Fochmod
    //printf("gasteiger mark #4 \n");

    // return an estimation of the number of charged groups in the molecule
    return total_charges; 

}//computegastcharges

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

