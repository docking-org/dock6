//definition of class Base_Grid
//
#ifndef BASE_GRID_H
#define BASE_GRID_H 

#include <string>
#include "utils.h"  // XYZCRD


class           Base_Grid {

  public:

    Base_Grid();
    virtual ~Base_Grid() ;

    int             span[3];          //number of bins for each dimension
    int             size;             //total number of items in array
    float           spacing,          //grid space size (in Ang)
                    origin[3];        //starting corner of grid
    float           x_min, x_max,     //min and max coordinates
                    y_min, y_max,     //in xyz dimensions
                    z_min, z_max;  
    XYZCRD          corners[8];       //xyz coordinates of corners of box
    int             neighbors[8];     //array index of all cubes surrounding
                                      //atom of interest
    float           cube_coords[3];   //xyz coordinates of cube of interest
                                      //translated to origin of 0,0,0
    int             nearest_neighbor; //closest corner of cube to atom
                                      //of interest
    unsigned char  *bump;             //array to hold bumps
                                      //needed as base variable because
                                      //header information about all grids
                                      //is contained in bump grid and thus
                                      //must be read in by all grids

    void            read_header(std::string filename);
    void            calc_corner_coords();
    bool            is_inside_grid_box(float x, float y, float z);
    void            find_grid_neighbors(float x, float y, float z);
    float           interpolate(float *grid);

  //private: //inorder to compile properly after GIST addition LEP

    int             find_grid_index(int, int, int);

};

#endif  // BASE_GRID_H

