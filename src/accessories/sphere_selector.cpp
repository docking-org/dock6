#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

#include "sphere_selector.h"


int main(int argc, char **argv){

	Spheres		c_spheres;

	//check command line arguments
	if(argc != 4){
		cout << "Usage:  sphere_selector <sphere_file.sph> <ligand.mol2> <distance>." << endl;
		exit(0);
	}

	//import a non-redundant list of spheres and parse it
	FILE		*sphere_file;
	sphere_file = fopen(argv[1], "r");
	if(sphere_file == NULL){
		cout << "Could not open " << argv[1] << " for reading.  Program will terminate." << endl << endl;
		exit(0);
	}
	c_spheres.read_spheres(sphere_file);

	//import target atom coordinates
	FILE		*atom_file;
	atom_file = fopen(argv[2], "r");
	if(atom_file == NULL){
		cout << "Could not open " << argv[2] << " for reading.  Program will terminate." << endl << endl;
		exit(0);
	}
	c_spheres.read_atoms(atom_file);

	//loop over spheres and flag if within cutoff
	float	dist, cutoff;
	bool	selected;
	c_spheres.num_selected = 0;
	cutoff = atof(argv[3]);
	for(int i=0;i<c_spheres.num_spheres;i++){
		selected = false;
		for(int j=0;j<c_spheres.atoms;j++){
			dist = c_spheres.distance(i, j); 
			if(dist < cutoff){
				selected = true;
			}
		}
		if(selected){
			c_spheres.tot_spheres[i].selected = true;
			c_spheres.num_selected++;
		}
		else{c_spheres.tot_spheres[i].selected = false;}
	}

	//print out selected spheres
	c_spheres.print_output(argv[3]);	
		
	return 0;
}

//////////////////////////////////////////////////////
Spheres::~Spheres(){

	delete[] atom_coord;
}

//////////////////////////////////////////////////////
void Spheres::read_spheres(FILE *sphere_file){

	char		line[500];
	bool		end;
	int			tmp, num, j;
	
	//flip over header
	fgets(line, 500, sphere_file);
	if( NULL == fgets(line, 500, sphere_file) ){
		cout << "Error:  premature end of sphere file!" << endl;
		exit(0);
	}

	end = false;
	num_spheres = j = 0;
	const int INITIAL_SIZE_OF_SPHERE_VECTOR = 19999;  // a guess by SRB
	tot_spheres.resize( INITIAL_SIZE_OF_SPHERE_VECTOR );

	//read in all spheres but cluster 0
	while(!end){
		if(sscanf(line, "cluster     0   number of spheres in cluster %i", &tmp)){
			end = true;
		}
		else if(sscanf(line, "cluster     %i    number of spheres in cluster   %i"
		               , &tmp, &num)){
			num_spheres += num;
			if ( tot_spheres.capacity() < num_spheres ) {
				tot_spheres.resize( num_spheres );
			}
			int		idx, idx2;
			char	x[10], y[10], z[10], radius[8];
			for(int i=0;i<num;i++){
				fgets(line, 500, sphere_file);
				sscanf(line, "%5i%s%s%s%s%i %*s", &idx, x, y, z, radius, &idx2);
				tot_spheres[j].idx = idx;
				tot_spheres[j].x = atof(x);
				tot_spheres[j].y = atof(y);
				tot_spheres[j].z = atof(z);
				tot_spheres[j].radius = atof(radius);
				tot_spheres[j].idx2 = idx2;
				j++;
			}
		}
		if( NULL == fgets(line, 500, sphere_file) ){
			cout << "Error:  premature end of sphere file!" << endl;
			exit(0);
		}
	}

	//debug
	/*for(int i=0;i<num_spheres;i++){
		cout << tot_spheres[i].x << endl;
	}*/

}

//////////////////////////////////////////////////////
void Spheres::read_atoms(FILE *atom_file){

	bool		head;
	char		line[500];

	//flip over header
	head = false;
	while(!head){
		fgets(line, 500, atom_file);
		if(sscanf(line, "@<TRIPOS>MOLECULE%*s")){
			fgets(line, 500, atom_file);
			fgets(line, 500, atom_file);
			head = true;
		}
	}
	sscanf(line, "%i %*s", &atoms);

	//flip over more header
	head = false;
	while(!head){
		fgets(line, 500, atom_file);
		if(sscanf(line, "@<TRIPOS>ATOM%*s")){
			fgets(line, 500, atom_file);
			head = true;
		}
	}

	//create array of atoms
	atom_coord = new XYZ[ atoms ];
	float	x, y, z;
	for(int i=0;i<atoms;i++){
		sscanf(line, "%*i %*s %f %f %f %*s", &x, &y, &z);
		atom_coord[i].coord[0] = x;
		atom_coord[i].coord[1] = y;
		atom_coord[i].coord[2] = z;
		fgets(line, 500, atom_file);
	}

	//debug
	/*for(int i=0;i<atoms;i++){
		cout << atom_coord[i].coord[0] << endl;
	}*/
}

/////////////////////////////////////////////////////
float Spheres::distance(int i, int j){

	return sqrt(pow((atom_coord[j].coord[0] - tot_spheres[i].x),2) + pow((atom_coord[j].coord[1] - tot_spheres[i].y),2) + pow((atom_coord[j].coord[2] - tot_spheres[i].z),2));

}

/////////////////////////////////////////////////////
void Spheres::print_output(char *dist){

	FILE * outfile;
	char	line[100];

	outfile = fopen("selected_spheres.sph", "w");
	sprintf(line, "DOCK spheres within %s ang of ligands\n", dist);
	fputs(line, outfile);
	sprintf(line, "cluster     1   number of spheres in cluster%6i\n", num_selected);
	fputs(line, outfile);
	for(int i=0;i<num_spheres;i++){
		if(tot_spheres[i].selected){
			sprintf(line, "%5i%10.5f%10.5f%10.5f%8.3f%5i 0  0\n", tot_spheres[i].idx, tot_spheres[i].x, tot_spheres[i].y, tot_spheres[i].z, tot_spheres[i].radius, tot_spheres[i].idx2); 
			fputs(line, outfile);
		}
	}
	fclose(outfile);

}
