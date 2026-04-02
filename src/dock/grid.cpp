// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// grid.cpp
//

#include <iostream>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "grid.h"
#include "trace.h"
#include "utils.h"

using namespace std;


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class Bump_Grid

// static member initializers
bool Bump_Grid :: got_the_grid = false;

// +++++++++++++++++++++++++++++++++++++++++
// Default constructor is a no op.
Bump_Grid :: Bump_Grid()
{
}

// +++++++++++++++++++++++++++++++++++++++++
//creates instance of bump grid and reads in data if not yet read in
void
Bump_Grid :: get_instance(string filename)
{
    if ( ! got_the_grid ) {
        assert( 0 != filename.size() );
        read_bump_grid( filename );
        // if no errors then
        got_the_grid = true;
    }
    return;
}

// +++++++++++++++++++++++++++++++++++++++++
//read in dimensions of Bump Grid and bumps calculated
//for receptor by GRID accessory program
void
Bump_Grid::read_bump_grid(string filename)
{
    string fname = filename + ".bmp";
    FILE* grid_in = fopen( fname.c_str(), "rb");
    if (grid_in == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    cout << " Reading the bump grid from " << fname << endl;

    fread(&size, sizeof(int), 1, grid_in);
    fread(&spacing, sizeof(float), 1, grid_in);
    fread(origin, sizeof(float), 3, grid_in);
    fread(span, sizeof(int), 3, grid_in);

    bump = new unsigned char[size];
    fread(bump, sizeof(unsigned char), size, grid_in);

    cout << " Done reading the bump grid." << endl;
    fclose(grid_in);

    calc_corner_coords();
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class Contact_Grid

// static member initializers
bool Contact_Grid :: got_the_grid = false;

// +++++++++++++++++++++++++++++++++++++++++
// Default constructor.
Contact_Grid :: Contact_Grid()
{
    cnt = NULL;
}

// +++++++++++++++++++++++++++++++++++++++++
// Destructor.
Contact_Grid :: ~Contact_Grid()
{
    delete [] cnt;
}

// +++++++++++++++++++++++++++++++++++++++++
//creates instance of contact grid and reads in data if not yet read in
void
Contact_Grid :: get_instance(string filename)
{
    if ( ! got_the_grid ) {
        assert( 0 != filename.size() );
        read_contact_grid( filename );
        // if no errors then
        got_the_grid = true;
    }
    return;
}

// +++++++++++++++++++++++++++++++++++++++++
//read in dimensions of Contact Grid and contacts calculated
//for receptor by GRID accessory program
void
Contact_Grid::read_contact_grid(string filename)
{

    string fname = filename + ".cnt";
    FILE* grid_in = fopen(fname.c_str(), "rb");
    if (grid_in == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    cout << " Reading the contact grid from " << fname << endl;

    fread(&size, sizeof(int), 1, grid_in);
    cnt = new short int[size];
    fread(cnt, sizeof(short int), size, grid_in);

    cout << " Done reading the contact grid." << endl;
    fclose(grid_in);

    read_header(filename);
    
    calc_corner_coords();
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class Energy_Grid

// static member initializers
bool Energy_Grid :: got_the_grid = false;

// +++++++++++++++++++++++++++++++++++++++++
// Default constructor is a no op.
Energy_Grid :: Energy_Grid()
{
  avdw = NULL;
  bvdw = NULL;
  es   = NULL;
}

// +++++++++++++++++++++++++++++++++++++++++
Energy_Grid::~Energy_Grid()
{
    Trace trace( "Energy_Grid::~Energy_Grid()" );
    this->clear_grid();
}


// +++++++++++++++++++++++++++++++++++++++++
void Energy_Grid::clear_grid()
{
    delete[] avdw;
    delete[] bvdw;
    delete[] es;
}


// +++++++++++++++++++++++++++++++++++++++++
//creates instance of energy grid and reads in data if not yet read in
void
Energy_Grid :: get_instance(string filename)
{
    if ( ! got_the_grid ) {
        assert( 0 != filename.size() );
        read_energy_grid( filename );
        // if no errors then
        got_the_grid = true;
    }
    return;
}


// +++++++++++++++++++++++++++++++++++++++++
void
Energy_Grid::read_energy_grid(string filename)
{

    string fname = filename + ".nrg";
    FILE* grid_in;
    grid_in = fopen(fname.c_str(), "rb");
    if (grid_in == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    cout << " Reading the energy grid from " << fname << endl;

    fread(&size, sizeof(int), 1, grid_in);
    fread(&atom_model, sizeof(int), 1, grid_in);
    fread(&att_exp, sizeof(int), 1, grid_in);
    fread(&rep_exp, sizeof(int), 1, grid_in);

    avdw = new float[size];
    bvdw = new float[size];
    es = new float[size];

    fread( bvdw, sizeof(float), size, grid_in);
    fread( avdw, sizeof(float), size, grid_in);
    fread( es, sizeof(float), size, grid_in);

    cout << " Done reading the energy grid." << endl;
    fclose(grid_in);
    
    read_header(filename);

    calc_corner_coords();
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class GIST_Grid

// static member initializers
bool GIST_Grid :: got_the_grid = false;

// +++++++++++++++++++++++++++++++++++++++++
// Default constructor is a no op.
GIST_Grid :: GIST_Grid()
{
  gist   = NULL;
  //got_the_grid = false;
}

// +++++++++++++++++++++++++++++++++++++++++
GIST_Grid::~GIST_Grid()
{
    Trace trace( "GIST_Grid::~GIST_Grid()" );
    this->clear_grid();
}


// +++++++++++++++++++++++++++++++++++++++++
void GIST_Grid::clear_grid()
{
    delete[] gist;
}


// +++++++++++++++++++++++++++++++++++++++++
//creates instance of gist grid and reads in data if not yet read in
void
GIST_Grid :: get_instance(string filename)
{
    if ( ! got_the_grid ) {
        assert( 0 != filename.size() );
        read_gist_grid( filename );
        // if no errors then
        got_the_grid = true;
    }
    return;
}


// +++++++++++++++++++++++++++++++++++++++++
void
GIST_Grid::read_gist_grid(string filename)
{

    string fname = filename;
    float  *read_gist; // array to read in.  
    char temp[20];
    int  xn,yn,zn;   // grid demintions
    float  xc,yc,zc; // origin of the grid box
    float  xs,ys,zs; // step size
    float  xl,yl,zl; // length of grid box in all dimintions
    int gridlen, numlines, remainder,count,i;   // grid length
    FILE* grid_in;
    grid_in = fopen(fname.c_str(), "r");
    if (grid_in == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    cout << " Reading the gist grid from " << fname << endl;
//    fread(&temp, sizeof(char[100]), 1, grid_in);
//    cout << temp << endl;
//    fread(&temp, sizeof(char[100]), 1, grid_in);
//    cout << temp << endl;
//    fread(&temp, sizeof(char[100]), 1, grid_in);
//    cout << temp << endl;

    char line[500];
    int  inittemp;
    //char header[20];
    // read in frist line and get grid dimintions
    fgets(line, 500, grid_in);
    sscanf(line, "%s %d %s %s %s %d %d %d\n", &temp, &inittemp,&temp,&temp,&temp,&xn,&yn,&zn);
    cout << xn <<" "<<yn<<" "<<zn << endl;

    // read in second line and get the center of grid
    fgets(line, 500, grid_in);
    sscanf(line, "%s %f %f %f\n", &temp,&xc,&yc,&zc);
    cout << xc <<" "<<yc<<" "<<zc << endl;
    // read in 3rd to 5th lines and get the stepsize
    fgets(line, 500, grid_in);
    sscanf(line, "%s %f %f %f\n", &temp,&xs,&inittemp,&inittemp);
    fgets(line, 500, grid_in);
    sscanf(line, "%s %f %f %f\n", &temp,&inittemp,&ys,&inittemp);
    fgets(line, 500, grid_in);
    sscanf(line, "%s %f %f %f\n", &temp,&inittemp,&inittemp,&zs);
    cout << xs <<" "<<ys<<" "<<zs << endl;

    if (xs != ys || xs != zs || ys != zs){ 
       cout << "gist grid spacing ("<< xs <<", "<<ys<<", "<<zs << ") must be equal" << endl;
       exit(0);
    }
    fgets(line, 500, grid_in); // grid connections
    fgets(line, 500, grid_in); // number of gridpoints
    sscanf(line, "%s %d %s %s %s %s %s %d %s %d %s %s\n", &temp,&inittemp,&temp,&temp,&temp,&temp,&temp,&inittemp,&temp,&gridlen,&temp,&temp);
    cout << gridlen << endl;
    read_gist = new float[gridlen];// read in array to this one. 
    gist = new float[gridlen];
    numlines = gridlen / 3;
    remainder = gridlen % 3;
    cout << gridlen << " " << numlines << " " << remainder << endl;

    // calculate length 
    xl = float(xn)*xs;
    yl = float(yn)*ys;
    zl = float(zn)*zs;
    cout << xl <<" "<<yl<<" "<<zl << endl;

    size = gridlen;
    origin[0] = xc; origin[1] = yc; origin[2] = zc; // populate origin base_grid, calc_corners needs this.  
    span[0] = xn; span[1] = yn; span[2] = zn; // populate span.  
    //span[0] = xn+1; span[1] = yn+1; span[2] = zn+1; // populate span.  
    spacing = xs; // xs == ys == zs 
    vol = xs*ys*zs;

    // loop over all the lines, they will contain 3 floats each
    count=0;  
    for (i = 0; i < numlines; i++) {
        fgets(line, 500, grid_in); // get next line
        sscanf(line, "%f %f %f\n",&read_gist[count],&read_gist[count+1],&read_gist[count+2]);
        count=count+3;
    }
    // accept the last one might have 2 or 1 float(s). 
    if (remainder == 2) { 
       fgets(line, 500, grid_in); // get last line
       sscanf(line, "%f %f\n",&read_gist[count],&read_gist[count+1]);
    } 
    else if(remainder == 1) {
       fgets(line, 500, grid_in); // get last line
       sscanf(line, "%f\n",&read_gist[count]);
    }
    else if(remainder == 0) {
       // do nothing
    }
    else {
       cout << "remainder = " << remainder << endl;
       cout << "something is wrong . . ." << endl;
       exit(0);
    }
    cout << " Done reading the gist grid." << endl;
    fclose(grid_in);

    // flip the order of x,y,z -> z,y,x
    int ind = 0;
    int ind2 = 0;
    for(int i =0; i < span[0]; i++){//x
       for(int j =0; j < span[1]; j++){//y
          for(int k =0; k < span[2]; k++){//z
             ind = find_grid_index(i, j, k);
             gist[ind] = read_gist[ind2];
             //cout<<"Index: "<<ind<<" "<<ind2<<endl;
             ind2++;
          }
       }
    }
     
/*
    // make sure things are right
    float  cx,cy,cz; // cords. 
    count = 0;
    for (int xi = 0; xi < xn; xi++){
       cx = float(xi)*spacing + origin[0];
       for (int yi = 0; yi < yn; yi++){
           cy = float(yi)*spacing + origin[1];
           for (int zi = 0; zi < zn; zi++){
               cz = float(zi)*spacing + origin[2];
               if (gist[count] > 2.0 ) {
                   cout << cx <<" " << cy << " " << cz << " " << gist[count] << endl;
               } 
               count=count+1;
           }
       }
    } 
*/    

    calc_corner_coords();

    //cout << spacing << " " << x_min << " " << x_max << " " << y_min << " " << y_max << " " << z_min << " " << z_max << endl;
    // pad
    /*
    float pad = 1.0;
    x_min = x_min + pad;
    x_max = x_max - pad;
    y_min = y_min + pad;
    y_max = y_max - pad;
    z_min = z_min + pad;
    z_max = z_max - pad;
    */
    delete[] read_gist; 
    //exit(0);
}

/******************************************/
float
//GIST_Grid :: atomic_displacement(float x, float y, float z, float radius, bool* boolgrid)
GIST_Grid :: atomic_displacement(float x, float y, float z, float radius, std::map<unsigned int, bool> &displacemap)
{   
    //float radius2 = pow(radius,2); // radius squard
    float radius2 = radius*radius; // radius squard
    float dist2;
    float sum = 0;
    int ind;
    //int box_about_side = ceil(radius / spacing * 2.0 ); // this is the size of one side of  box that contains the atom.
    //int hbox_about_side = ceil(radius / spacing ); // this is half the size of one side of  box that contains the atom.
    int box_about_side  = (ceil(radius / spacing)+1.0) * 2.0 ; // this is the size of one side of  box that contains the atom.
    int hbox_about_side =  ceil(radius / spacing)+1.0;         // this is half the size of one side of  box that contains the atom.
    //cout << radius << " " << spacing << " " << box_about_side <<endl;

    float x_int = x - origin[0];
    float y_int = y - origin[1];
    float z_int = z - origin[2];

    int x_nearest = NINT(x_int / spacing);
    int y_nearest = NINT(y_int / spacing);
    int z_nearest = NINT(z_int / spacing);

    nearest_neighbor = find_grid_index(x_nearest, y_nearest, z_nearest);

    // we will start at a corner of the bounding box. 
    int start_corner_x = x_nearest - hbox_about_side;
    int start_corner_y = y_nearest - hbox_about_side;
    int start_corner_z = z_nearest - hbox_about_side;

    for (int xa = start_corner_x; xa < start_corner_x + box_about_side; xa++){
        for (int ya = start_corner_y; ya < start_corner_y + box_about_side; ya++){
            for (int za = start_corner_z; za < start_corner_z + box_about_side; za++){
                float xcor = float(xa) * spacing + origin[0];
                float ycor = float(ya) * spacing + origin[1];
                float zcor = float(za) * spacing + origin[2];
                dist2 = pow((xcor-x),2.0)+pow((ycor-y),2.0) + pow((zcor-z),2.0);
                //cout << dist2<< " " << radius2 << endl;
                if (dist2<=radius2) {
                   //cout << "I AM HERE" << endl;
                   ind = find_grid_index(xa, ya, za);
                   if (ind > size) { 
                      cout << "Warning index too big for grid ..." << endl;
                   } 
                   else {
                      bool dflag = displacemap[ind];
                      if (!(dflag)){ // check if grid point is already displaced by another atom. 
                         //cout << "I AM HERE" << endl;
                         //cout <<"ATOM      1  N   LEU A   4    "<<xcor<<" "<<ycor<<" "<<zcor<<" "<< gist[ind] <<endl;
                         displacemap[ind] = true;
                         sum = sum+gist[ind];
                      }
                  }
                }
            }//za  
        }//ya
    }//xa
    return sum;
}

/******************************************/
float
GIST_Grid :: atomic_blurry_displacement(float x, float y, float z, float radius, float sigma2)
{   
    // float sigma2 = pow(radius/2.0,2.0);
    //float radius2 = pow(radius,2); // radius squard
    float radius2 = radius*radius; // radius squard
    float dist2;
    float sum = 0;
    int ind;
    //int box_about_side = ceil(radius / spacing * 2.0 ); // this is the size of one side of  box that contains the atom.
    int box_about_side  = (ceil(radius / spacing)+1.0) * 2.0 ; // this is the size of one side of  box that contains the atom.
    int hbox_about_side =  ceil(radius / spacing)+1.0;         // this is half the size of one side of  box that contains the atom.
    //cout << radius << " " << spacing << " " << box_about_side <<endl;

    float x_int = x - origin[0];
    float y_int = y - origin[1];
    float z_int = z - origin[2];

    int x_nearest = NINT(x_int / spacing);
    int y_nearest = NINT(y_int / spacing);
    int z_nearest = NINT(z_int / spacing);

    nearest_neighbor = find_grid_index(x_nearest, y_nearest, z_nearest);

    // we will start at a corner of the bounding box. 
    int start_corner_x = x_nearest - hbox_about_side;
    int start_corner_y = y_nearest - hbox_about_side;
    int start_corner_z = z_nearest - hbox_about_side;

    for (int xa = start_corner_x; xa < start_corner_x + box_about_side; xa++){
        for (int ya = start_corner_y; ya < start_corner_y + box_about_side; ya++){
            for (int za = start_corner_z; za < start_corner_z + box_about_side; za++){
                float xcor = float(xa) * spacing + origin[0];
                float ycor = float(ya) * spacing + origin[1];
                float zcor = float(za) * spacing + origin[2];
                //cout << xa << " " << ya << " " << za << endl;
                //cout << xcor << " " << ycor << " " << zcor << endl;
                dist2 = pow((xcor-x),2.0)+pow((ycor-y),2.0) + pow((zcor-z),2.0);
                if (dist2<=radius2) {
                   //float gausian_weight = (1.0/sqrt(2.0 * PI * sigma2)) * pow(EXP,(-1.0*dist2/(2.0 * sigma2)));  
                   float gausian_weight = (1.0/sqrt(2.0 * PI * sigma2)) * exp(-1.0*dist2/(2.0 * sigma2));  
                   ind = find_grid_index(xa, ya, za);
                   if (ind > size) { 
                      cout << "Warning index too big for grid ..." << endl;
                   } 
                   else {
                      //cout << "gist = " << gist[ind] << endl;
                      //cout << "g_weight = " << gausian_weight << endl;
                      sum = sum+gausian_weight*gist[ind];
                   }
                }
            }//za  
        }//ya
    }//xa
    //cout << "sum = " << sum << endl;
    return sum;
}


// This will write out the grid for the displaced region 
void
GIST_Grid::write_gist_grid(string fname, bool* boolgrid)
{
    FILE* grid_out = fopen(fname.c_str(), "w");
    int count = 0;
    int ind = 0;
    int ind2 = 0;
//    for (int ind = 0; ind < size;ind++){
    for(int i =0; i < span[0]; i++){//x
       for(int j =0; j < span[1]; j++){//y
          for(int k =0; k < span[2]; k++){//z
               //ind = find_grid_index(k, j, i);
               ind = find_grid_index(i, j, k);
               if (boolgrid[ind]){
                  fprintf(grid_out, "%f",gist[ind]);
                  //fprintf(grid_out, "1.0");
               } else{
                  fprintf(grid_out, "0.0");   
               }
               if (count < 2){
                   fprintf(grid_out, " ");
                   count = count + 1;
               }else{
                   //cout << "I AM HERE" << endl;
                   fprintf(grid_out, "\n");
                   count = 0;
                   //exit(0);
               }
               cout<<"Index: "<<ind<<" "<<ind2<<endl;
               ind2++;
          }   
       }   
    }   
//    }   
    fclose(grid_out);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class GB_Grid for gbsa_zou score

// static member initializers
bool GB_Grid :: got_the_grid = false;

// +++++++++++++++++++++++++++++++++++++++++
// Default constructor is a no op.
GB_Grid :: GB_Grid()
{
}

// +++++++++++++++++++++++++++++++++++++++++
//creates instance of gb grid and reads in data if not yet read in
void
GB_Grid :: get_instance(string filename)
{
    if ( ! got_the_grid ) {
        assert( 0 != filename.size() );
        read_gb_grid( filename );
        // if no errors then
        got_the_grid = true;
    }
    return;
}

// +++++++++++++++++++++++++++++++++++++++++
void
GB_Grid::read_gb_grid(string file_prefix)
{
    int             i;

    // read from .bmp file
    string fname = file_prefix + ".bmp";
    FILE* bmpfp = fopen(fname.c_str(), "r");
    if (bmpfp == NULL) {
        cout << "\n\nCould not open " << fname << " for reading."
            << "  Program will terminate." << endl << endl;
        exit(0);
    }

    cout << " Reading from GB bump grid file " << fname << endl;
    char line[500];
    fgets(line, 500, bmpfp);
    fgets(line, 500, bmpfp);
    sscanf(line, "%f %f %f %f %d %d %d", &gb_grid_spacing,
           &gb_grid_origin[0], &gb_grid_origin[1], &gb_grid_origin[2],
           &gb_grid_span[0], &gb_grid_span[1], &gb_grid_span[2]);
    fgets(line, 500, bmpfp);
    sscanf(line, "%f %f %d %f %f %f", &gb_grid_grdcut, &gb_grid_grdcuto,
           &gb_rec_total_atoms, &gb_grid_cutsq, &gb_grid_f_scale,
           &gb_rec_solv_rec);

    cout << " Done reading from GB bump grid file." << endl;
    fclose(bmpfp);

    gb_grid_size = gb_grid_span[0] * gb_grid_span[1] * gb_grid_span[2];

    // resize receptor vectors
    gb_rec_charge.clear();
    gb_rec_coords.clear();
    gb_rec_inv_a.clear();
    gb_rec_inv_a_rec.clear();
    gb_rec_radius_eff.clear();
    gb_rec_receptor_screen.clear();

    gb_rec_charge.resize(gb_rec_total_atoms);
    gb_rec_coords.resize(gb_rec_total_atoms);
    gb_rec_inv_a.resize(gb_rec_total_atoms);
    gb_rec_inv_a_rec.resize(gb_rec_total_atoms);
    gb_rec_radius_eff.resize(gb_rec_total_atoms);
    gb_rec_receptor_screen.resize(gb_rec_total_atoms);

    // read from .rec file
    fname = file_prefix + ".rec";
    FILE* recfp = fopen(fname.c_str(), "r");
    if (recfp == NULL) {
        cout << "\n\nCould not open " << fname << " for reading."
            << "  Program will terminate." << endl << endl;
        exit(0);
    }

    cout << " Reading from GB receptor grid file " << fname << endl;
    for (i = 0; i < gb_rec_total_atoms; i++) {
        fscanf(recfp, "%f %f %f %f %f %f %f", &gb_rec_coords[i].crd[0],
                &gb_rec_coords[i].crd[1], &gb_rec_coords[i].crd[2],
                &gb_rec_charge[i], &gb_rec_radius_eff[i],
                &gb_rec_inv_a_rec[i], &gb_rec_receptor_screen[i] );
    }

    cout << " Done reading from GB receptor grid file." << endl;
    fclose(recfp);
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class SA_Grid for gbsa_zou score

// static member initializers
bool SA_Grid :: got_the_grid = false;

// +++++++++++++++++++++++++++++++++++++++++
// Default constructor is a no op.
SA_Grid :: SA_Grid()
{
}

// +++++++++++++++++++++++++++++++++++++++++
//creates instance of sa grid and reads in data if not yet read in
void
SA_Grid :: get_instance(string filename)
{
    if ( ! got_the_grid ) {
        assert( 0 != filename.size() );
        read_sa_grid( filename );
        // if no errors then
        got_the_grid = true;
    }
    return;
}

// +++++++++++++++++++++++++++++++++++++++++
void
SA_Grid::read_sa_grid(string file_prefix)
{
    int             i;
    int             j;

    // read the .bmp file
    string fname = file_prefix + ".bmp";
    FILE* ifp = fopen(fname.c_str(), "r");
    if (ifp == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    cout << " Reading from SA bump grid file " << fname << endl;
    char line[500];
    fgets(line, 500, ifp);
    char header[20];
    sscanf(line, "%19s", header);
    fgets(line, 500, ifp);
    sscanf(line, "%f %f %f %f %d %d %d", &sa_grid_spacing,
           &sa_grid_origin[0], &sa_grid_origin[1], &sa_grid_origin[2],
           &sa_grid_span[0], &sa_grid_span[1], &sa_grid_span[2]);
    fgets(line, 500, ifp);
    sscanf(line, "%f %f %d %f %f %f %d", &sa_grid_grdcut, &sa_grid_grdcuto,
           &sa_rec_total_atoms, &sa_grid_cutsq, &sa_grid_f_scale,
           &sa_rec_solv_rec, &sa_grid_nvtyp);

    cout << " Done reading from SA bump grid file." << endl;
    fclose(ifp);

    sa_grid_size = sa_grid_span[0] * sa_grid_span[1] * sa_grid_span[2];

    // resize receptor vectors
    sa_rec_coords.clear();
    sa_rec_radius_eff.clear();
    sa_rec_vdwn_rec.clear();
    sa_rec_nsphgrid.clear();
    sa_rec_nhp.clear();

    sa_rec_coords.resize(sa_rec_total_atoms);
    sa_rec_radius_eff.resize(sa_rec_total_atoms);
    sa_rec_vdwn_rec.resize(sa_rec_total_atoms);
    sa_rec_nsphgrid.resize(sa_grid_nvtyp);
    sa_rec_nhp.resize(sa_grid_nvtyp);

    // read the .sas file
    fname = file_prefix + ".sas";
    ifp = fopen(fname.c_str(), "r");
    if (ifp == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    cout << " Reading from SA receptor grid file " << fname << endl;

    for (i = 0; i < sa_rec_total_atoms; i++) {
        fscanf(ifp, "%f %f %f %f %d", &sa_rec_coords[i].crd[0],
                &sa_rec_coords[i].crd[1], &sa_rec_coords[i].crd[2],
                &sa_rec_radius_eff[i], &sa_rec_vdwn_rec[i] );
    }

    for (i = 0; i < sa_grid_nvtyp; i++) {
        fscanf(ifp, "%d %d", &sa_rec_nsphgrid[i], &sa_rec_nhp[i] );
    }

    cout << " Done reading from SA receptor grid file." << endl;
    fclose(ifp);

    sa_grid_sphgrid_crd.resize(sa_grid_nvtyp);
    for (i = 0; i < sa_grid_nvtyp; i++) {
        sa_grid_sphgrid_crd[i].resize(sa_rec_nsphgrid[i]);
    }

    sa_grid_mark_sas.resize(sa_rec_total_atoms);
    for (i = 0; i < sa_rec_total_atoms; i++) {
        sa_grid_mark_sas[i].resize(sa_rec_nsphgrid[sa_rec_vdwn_rec[i] - 1]);
    }

    // read the .sasmark file
    fname = file_prefix + ".sasmark";
    ifp = fopen(fname.c_str(), "r");
    if (ifp == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    cout << " Reading from SA receptor info grid file " << fname << endl;
    fgets(line, 500, ifp);
    sscanf(line, "%d %d %d", &sa_rec_nsas, &sa_rec_nsas_hp,
           &sa_rec_nsas_pol);
    fgets(line, 500, ifp);
    sscanf(line, "%f %f %f", &sa_grid_r_probe, &sa_grid_sspacing,
           &sa_grid_r2_cutoff);

    for (i = 0; i < sa_grid_nvtyp; i++) {
        for (j = 0; j < sa_rec_nsphgrid[i]; j++) {
            fscanf(ifp, "%f %f %f", &sa_grid_sphgrid_crd[i][j].crd[0],
                   &sa_grid_sphgrid_crd[i][j].crd[1],
                   &sa_grid_sphgrid_crd[i][j].crd[2]);
        }
    }

    for (i = 0; i < sa_rec_total_atoms; i++) {
        for (j = 0; j < sa_rec_nsphgrid[sa_rec_vdwn_rec[i] - 1]; j++) {
            fscanf(ifp, "%d", &sa_grid_mark_sas[i][j]);
        }
    }

    cout << " Done reading from SA receptor info grid file." << endl;
    fclose(ifp);
}

