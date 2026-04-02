#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "sphere.h"

using namespace std;

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// class Sphere

void
write_spheres( string sphere_file_name, SphereVec & spheres, std::string context)
{
    // open the sphere file
    FILE           *sphere_file;
    sphere_file = fopen(sphere_file_name.c_str(), "w");

    if (context.find('\n') == std::string::npos) {
        // Append a newline character to the string
        context += '\n';
    }

    //convert c-string to char buffer
    size_t length = context.length();

    // Allocate memory for a buffer to hold the cstring
    char buffer[length + 1]; // Add 1 for the null terminator

    // Copy the cstring into the buffer
    strcpy(buffer, context.c_str());

    if (sphere_file == NULL) {
        cout << "\n\nERROR:  Could not open " << sphere_file_name
             << " for reading.  Program will terminate.\n" << endl;
        exit(0);
    }
    
    //write the first line to include context for the spheres
    fprintf(sphere_file, buffer);
    //print cluster info
    int cluster_index = 1;
    fprintf(sphere_file, "cluster     %d   number of spheres in cluster    %d\n",
            cluster_index, spheres.size());

    // read in the color table (if any)
    for (size_t i = 0; i < spheres.size(); ++i) {
        char sphere_line[BUFFER_SIZE];
        strcpy(sphere_line, spheres[i].write_to_buffer().c_str());
        fprintf(sphere_file, "%s \n", sphere_line);
    }
    fclose(sphere_file);
}

// +++++++++++++++++++++++++++++++++++++++++
std::string
Sphere::write_to_buffer()
{
    //string buffer;
    char buffer[BUFFER_SIZE];
    int length; // Store the length of the formatted string
    double temp_radius = 1.0; // temporary radius for sphere generated from atom,
                              // this should probably be vdw radius.

    length = snprintf(buffer, BUFFER_SIZE, "%4d %9.5f %9.5f %9.5f %7.3f %4d %2d %2s",
                this->surface_point_i, this->crds.x, this->crds.y, this->crds.z,
                temp_radius, this->surface_point_j, this->critical_cluster,
                this->color.c_str());

    // Check if snprintf encountered an error
    if (length < 0 || length >= BUFFER_SIZE) {
        // Handle error
        printf("Error: snprintf encountered an error or buffer overflow.\n");
    } 
    std::string sphere_line(buffer);
    return sphere_line;
}



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

        if (color_int > 0) {
            for (size_t i = 0; i < site_color_ints.size(); ++i) {
                if (color_int == site_color_ints[i]) {
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

