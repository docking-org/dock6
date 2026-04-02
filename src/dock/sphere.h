// definition of classes Sphere and Active_Site_Spheres

#ifndef SPHERE_H
#define SPHERE_H

#include <string>
#include <vector>
#include "utils.h"  // XYZCRD
class Parameter_Reader;



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class Sphere is a representation of a sphere file produced by sphgen.
//
class           Sphere {

  public:
    XYZCRD          crds;
    float           radius;
    int             surface_point_i;
    int             surface_point_j;
    int             critical_cluster;
    std::string     color;

    Sphere() {
        clear();
    };
    ~Sphere() {
        clear();
    };

    void            clear();
    float           distance(Sphere &);
    Sphere & operator=(const Sphere & s) {
        crds = s.crds;
        radius = s.radius;
        surface_point_i = s.surface_point_i;
        surface_point_j = s.surface_point_j;
        critical_cluster = s.critical_cluster;
        color = s.color;
        return (*this);
    };

};

typedef         std::vector< Sphere > SphereVec;

// input the spheres from sphere_file_name into vector spheres and
// return the number of spheres read.
int read_spheres( std::string sphere_file_name, SphereVec & spheres );


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class Active_Site_Spheres is the sphere list from the sphgen sphere file.
// It is a community resource.
// This class uses the Singleton pattern: 
// use Active_Site_Spheres :: get_instance() to access the sphere list.
//
class Active_Site_Spheres {
public:
    typedef SphereVec :: const_iterator const_iterator ;
    static SphereVec & get_instance();
    static void set_sphere_file_name( Parameter_Reader & parm );
private:
    Active_Site_Spheres();
    static std::string sphere_file_name;
    static SphereVec the_instance;
};

#endif  // SPHERE_H

