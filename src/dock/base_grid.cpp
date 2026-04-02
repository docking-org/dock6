// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// base_grid.cpp
//
// This class contains functions that are used to identify a) where the grid is 
// in space and b) what grid box the atom of interest is in.
//

#include <iostream>
#include <assert.h>
#include "base_grid.h"
#include "trace.h"

using namespace std;

// +++++++++++++++++++++++++++++++++++++++++
// delete header information from .bmp file
Base_Grid::Base_Grid()
{
     bump = NULL;
}

// +++++++++++++++++++++++++++++++++++++++++
// delete header information from .bmp file
Base_Grid::~Base_Grid()
{
    Trace trace( "Base_Grid::~Base_Grid()" );
    delete [] bump;
}

// +++++++++++++++++++++++++++++++++++++++++
//reads header information from .bmp file
void
Base_Grid::read_header(string filename)
{
    Trace trace( "Base_Grid::read_header(string filename)" );
    string fname = filename + ".bmp";
    FILE* grid_in;
    grid_in = fopen( fname.c_str(), "rb");
    if (grid_in == NULL) {
        cout << "\n\nCould not open " << fname <<
            " for reading.  Program will terminate." << endl << endl;
        exit(0);
    }

    cout << " Reading the grid box quantifiers from " << fname << endl;
    fread(&size, sizeof(int), 1, grid_in);
    fread(&spacing, sizeof(float), 1, grid_in);
    fread(origin, sizeof(float), 3, grid_in);
    fread(span, sizeof(int), 3, grid_in);

    if (Parameter_Reader::verbosity_level() > 0) {
        cout <<   '\t' << "size = " << size
             << "\n\t" << "spacing = " << spacing
             << "\n\t" << "origin = " << origin[0]<<","<<origin[1]<<","<<origin[2]
             << "\n\t" << "span = " << span[0]<<","<<span[1]<<","<<span[2]
             << '\n';
    }

    cout << " Done reading the grid box quantifiers." << endl;
    fclose(grid_in);
}

// +++++++++++++++++++++++++++++++++++++++++
//computes the minimum and maximum xyz coordinates for the dimensions
//of the grid, and collects them as the corners of the grid box
void
Base_Grid::calc_corner_coords()
{
    Trace trace( "Base_Grid::calc_corner_coords()" );
    float dx = (span[0] - 1) * spacing;
    float dy = (span[1] - 1) * spacing;
    float dz = (span[2] - 1) * spacing;

    x_min = origin[0];
    x_max = origin[0] + dx;
    y_min = origin[1];
    y_max = origin[1] + dy;
    z_min = origin[2];
    z_max = origin[2] + dz;

    corners[0].assign_vals(x_min, y_min, z_min);
    corners[1].assign_vals(x_max, y_min, z_min);
    corners[2].assign_vals(x_max, y_min, z_max);
    corners[3].assign_vals(x_min, y_min, z_max);
    corners[4].assign_vals(x_min, y_max, z_min);
    corners[5].assign_vals(x_max, y_max, z_min);
    corners[6].assign_vals(x_max, y_max, z_max);
    corners[7].assign_vals(x_min, y_max, z_max);
}

// +++++++++++++++++++++++++++++++++++++++++
//determines whether atom is inside the boundaries of the grid box
//if not, return false
bool
Base_Grid::is_inside_grid_box(float x, float y, float z)
{
    return   (x > x_min+spacing && x < x_max-spacing &&
              y > y_min+spacing && y < y_max-spacing &&
              z > z_min+spacing && z < z_max-spacing) ;
    // Since we are using bicubic interpolation in the grid energy compuation,
    // any (x,y,z) point needs to be between x_min+spacing and x_max-spacing
    // in order to return valid scores for all 8 points being interpolated
    // For comparing floats, we need to account for round off error
}

// +++++++++++++++++++++++++++++++++++++++++
//computes index of array where value is stored
int 
Base_Grid::find_grid_index(int x, int y, int z)
{
    return span[0] * span[1] * z + span[0] * y + x;
}

// +++++++++++++++++++++++++++++++++++++++++
//finds the indices of the cubes that are neighboring
//the atom of interest
void
Base_Grid::find_grid_neighbors(float x, float y, float z)
{
    float x_int = x - origin[0];
    float y_int = y - origin[1];
    float z_int = z - origin[2];

    int x_nearest = NINT(x_int / spacing);
    int x_below = INTFLOOR(x_int / spacing);
    int x_above = x_below + 1;
    if (x_nearest >= span[0]) {
        if (x_below >= span[0])
            x_nearest = span[0] - 1;
        else
            x_nearest = x_below;
    }
    if (x_nearest < 0) {
        if (x_above < 0)
            x_nearest = 0;
        else
            x_nearest = x_above;
    }

    int y_nearest = NINT(y_int / spacing);
    int y_below = INTFLOOR(y_int / spacing);
    int y_above = y_below + 1;
    if (y_nearest >= span[1]) {
        if (y_below >= span[1])
            y_nearest = span[1] - 1;
        else
            y_nearest = y_below;
    }
    if (y_nearest < 0) {
        if (y_above < 0)
            y_nearest = 0;
        else
            y_nearest = y_above;
    }

    int z_nearest = NINT(z_int / spacing);
    int z_below = INTFLOOR(z_int / spacing);
    int z_above = z_below + 1;
    if (z_nearest >= span[2]) {
        if (z_below >= span[2])
            z_nearest = span[2] - 1;
        else
            z_nearest = z_below;
    }
    if (z_nearest < 0) {
        if (z_above < 0)
            z_nearest = 0;
        else
            z_nearest = z_above;
    }

    nearest_neighbor = find_grid_index(x_nearest, y_nearest, z_nearest);

    neighbors[0] = find_grid_index(x_above, y_above, z_above);
    neighbors[1] = find_grid_index(x_above, y_above, z_below);
    neighbors[2] = find_grid_index(x_above, y_below, z_above);
    neighbors[3] = find_grid_index(x_below, y_above, z_above);
    neighbors[4] = find_grid_index(x_above, y_below, z_below);
    neighbors[5] = find_grid_index(x_below, y_above, z_below);
    neighbors[6] = find_grid_index(x_below, y_below, z_above);
    neighbors[7] = find_grid_index(x_below, y_below, z_below);

    float corrected_coords[3];
    corrected_coords[0] = x - origin[0];
    corrected_coords[1] = y - origin[1];
    corrected_coords[2] = z - origin[2];

    cube_coords[0] = corrected_coords[0] / spacing - (float) (INTFLOOR(corrected_coords[0] / spacing)); 
    cube_coords[1] = corrected_coords[1] / spacing - (float) (INTFLOOR(corrected_coords[1] / spacing));
    cube_coords[2] = corrected_coords[2] / spacing - (float) (INTFLOOR(corrected_coords[2] / spacing)); 
}

// +++++++++++++++++++++++++++++++++++++++++
//function for interpolating values from a set of 8 grid points
float
Base_Grid::interpolate(float *grid)
{
    float           a1,
                    a2,
                    a3,
                    a4,
                    a5,
                    a6,
                    a7,
                    a8;
    float           value;
    int             out_of_bounds,
                    i;

    out_of_bounds = 0;
    for (i = 0; i < 8; i++)
        if ((neighbors[i] > size) || (neighbors[i] < 0))
            out_of_bounds = 1;

    if (out_of_bounds == 0) {
        a8 = grid[neighbors[7]];
        a7 = grid[neighbors[6]] - a8;
        a6 = grid[neighbors[5]] - a8;
        a5 = grid[neighbors[4]] - a8;
        a4 = grid[neighbors[3]] - a8 - a7 - a6;
        a3 = grid[neighbors[2]] - a8 - a7 - a5;
        a2 = grid[neighbors[1]] - a8 - a6 - a5;
        a1 = grid[neighbors[0]] - a8 - a7 - a6 - a5 - a4 - a3 - a2;

        value =
            a1 * cube_coords[0] * cube_coords[1] * cube_coords[2] +
            a2 * cube_coords[0] * cube_coords[1] +
            a3 * cube_coords[0] * cube_coords[2] +
            a4 * cube_coords[1] * cube_coords[2] + a5 * cube_coords[0] +
            a6 * cube_coords[1] + a7 * cube_coords[2] + a8;

        return value;
    } else {
        return grid[nearest_neighbor];
    }
}

