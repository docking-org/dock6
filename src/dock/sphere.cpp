#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "sphere.h"

using namespace std;


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class Sphere


// +++++++++++++++++++++++++++++++++++++++++
void
Sphere::clear()
{
    crds.assign_vals(0.0, 0.0, 0.0);
    radius = 0.0;
    surface_point_i = 0;
    surface_point_j = 0;
    critical_cluster = 0;
    color.clear();
}


// +++++++++++++++++++++++++++++++++++++++++
float
Sphere::distance(Sphere & sph)
{
    return crds.distance(sph.crds);
}


// +++++++++++++++++++++++++++++++++++++++++
//
// input the spheres from sphere_file_name into vector spheres and
// return the number of spheres read.
//
int
read_spheres( string sphere_file_name, SphereVec & spheres )
{
    char            junk[50],
                    line[500],
                    label[50];
    FILE           *sphere_file;
    string          color_label;
    int             color_int;
    vector < string > site_color_labels;  // sphere file chem type labels
    INTVec          site_color_ints;      // sphere file chem type numbers

    int num_spheres = 0;

    // open the sphere file
    sphere_file = fopen(sphere_file_name.c_str(), "r");

    if (sphere_file == NULL) {
        cout << "\n\nERROR:  Could not open " << sphere_file_name
             << " for reading.  Program will terminate.\n" << endl;
        exit(0);
    }

    // read in the color table (if any)
    while (fgets(line, 500, sphere_file)) {
        sscanf(line, "%s", junk);

        // if line is a color table entry
        if (strcmp(junk, "color") == 0) {
            sscanf(line, "%s %s %d", junk, label, &color_int);
            color_label = label;

            site_color_labels.push_back(color_label);
            site_color_ints.push_back(color_int);
        }
        // if color table is finished
        if (strcmp(junk, "cluster") == 0)
            break;
    }

    // read in the spheres
    color_int = 0;
    Sphere          tmp;
    while (fgets(line, 500, sphere_file)) {

        sscanf(line, "%s", junk);
        if (strcmp(junk, "cluster") == 0)
            break;
        else
            sscanf(line, "%d %f %f %f %f %d %d %d", &tmp.surface_point_i,
                   &tmp.crds.x, &tmp.crds.y, &tmp.crds.z, &tmp.radius,
                   &tmp.surface_point_j, &tmp.critical_cluster, &color_int);

        // assign the proper color label
        tmp.color = "null";

        //cout <<  "color =" << color_int << endl;
        if (color_int > 0) {
            for (size_t i = 0; i < site_color_ints.size(); ++i) {
                if (color_int == site_color_ints[i]) {
                    //cout << " color = " <<  color_int << " " << site_color_ints[i] << endl;
                    tmp.color = site_color_labels[i];
                    break;
                }
            }
        }

        spheres.push_back(tmp);
        num_spheres++;
    }

    fclose(sphere_file);

    return num_spheres;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class Active_Site_Spheres

// static member initializers

string Active_Site_Spheres :: sphere_file_name = string();
SphereVec Active_Site_Spheres :: the_instance = vector< Sphere >();


// +++++++++++++++++++++++++++++++++++++++++
SphereVec &
Active_Site_Spheres :: get_instance()
{
    if ( 0 == the_instance.size() ) {
        assert( 0 != sphere_file_name.size() );
        size_t count = read_spheres( sphere_file_name, the_instance );
        assert( the_instance.size() == count );
    }
    return the_instance;
}


// +++++++++++++++++++++++++++++++++++++++++
void
Active_Site_Spheres :: set_sphere_file_name( Parameter_Reader & parm )
{
    if ( 0 == sphere_file_name.size() ) {
        sphere_file_name = parm.query_param("receptor_site_file",
                                            "receptor.sph");
    }
}

