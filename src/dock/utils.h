//
#ifndef UTILS_H
#define UTILS_H 

#include <map>
#include <string>
#include <utility>  // pair
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
class DOCKMol;


#define PI  3.1415926535897932384626433f
#define EXP 2.7182818284590452353602875f
#define NINT(x) (int) ((x) > 0 ? ((x) + 0.5) : ((x) - 0.5))
#define INTFLOOR(x) (int) (floor(x + 0.00001))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

const float     MIN_FLOAT = -999999986991104.00;
// Large negative number for initializing float values.
// Defined to overcome compiler specific NAN definitions.
// Note that this float has an exact IEEE-754 representation.


// +++++++++++++++++++++++++++++++++++++++++
class           TORSION {

  public:

    int             atom1;      // first atom bound to atom 2 (to define angle)
    int             atom2;      // first atom of rotatable bond
    int             atom3;      // second atom of rotatable bond
    int             atom4;      // first atom bound to atom 3 (to define angle)
    int             bond_num;   // bond number in Mol

};


// +++++++++++++++++++++++++++++++++++++++++
class           XYZCRD {
  public:
    float           x,
                    y,
                    z;

    XYZCRD(float _x = 0.0, float _y = 0.0, float _z = 0.0)
        : x(_x), y(_y), z(_z) {
    };
    ~XYZCRD() {
        x = y = z = 0.0;
    };

    void assign_vals(float a, float b, float c) {
        x = a;
        y = b;
        z = c;
    };
    float distance_squared( const XYZCRD & p ) const {
        return (p.x - x)*(p.x - x) + (p.y - y)*(p.y - y) + (p.z - z)*(p.z - z);
    };
    float distance( const XYZCRD & p ) const {
        return sqrt( distance_squared( p ) );
    };

    XYZCRD & operator=(const XYZCRD & xyz) {
        x = xyz.x;
        y = xyz.y;
        z = xyz.z;
        return (*this);
    };

};


// +++++++++++++++++++++++++++++++++++++++++
class           DOCKVector {

  public:
    double           x;
    double           y;
    double           z;

    DOCKVector      operator=(const DOCKVector &);
                    DOCKVector & operator+=(const DOCKVector &);
                    DOCKVector & operator+=(const float &);
                    DOCKVector & operator-=(const DOCKVector &);
                    DOCKVector & operator-=(const float &);
                    DOCKVector & operator*=(const float &);
                    DOCKVector & operator/=(const float &);
                    DOCKVector & normalize_vector();
    float           length() const;
    float           squared_length() const;
    float           squared_dist(const DOCKVector &) const;

    friend int      operator==(const DOCKVector &, const DOCKVector &);
    friend int      operator!=(const DOCKVector &, const DOCKVector &);

    friend DOCKVector operator+(const DOCKVector &, const DOCKVector &);
    friend DOCKVector operator+(const float &, const DOCKVector &);
    friend DOCKVector operator+(const DOCKVector &, const float &);
    friend DOCKVector operator-(const DOCKVector &, const DOCKVector &);
    friend DOCKVector operator-(const DOCKVector &);
    friend DOCKVector operator-(const DOCKVector &, const float &);
    friend DOCKVector operator*(const DOCKVector &, const DOCKVector &);
    friend DOCKVector operator*(const float &, const DOCKVector &);
    friend DOCKVector operator*(const DOCKVector &, const float &);
    friend DOCKVector operator/(const DOCKVector &, const float &);

    friend float    dot_prod(const DOCKVector &, const DOCKVector &);
    friend DOCKVector cross_prod(const DOCKVector &, const DOCKVector &);
    friend float    get_vector_angle(const DOCKVector &, const DOCKVector &);
    friend float    get_torsion_angle(DOCKVector &, DOCKVector &, DOCKVector &,
                                      DOCKVector &);

};


// +++++++++++++++++++++++++++++++++++++++++
typedef         std::vector < int >INTVec;
typedef         std::pair < int, int > INTPair;
typedef         std::vector < float >FLOATVec;
typedef         std::vector < XYZCRD > XYZVec;
typedef         std::vector < INTVec > INTVecVec;
typedef         std::vector < FLOATVec > FLOATVecVec;
typedef         std::vector < INTVecVec > INTVecVecVec;

typedef char    STRING20[21];

int             get_matrix_from_quaternion(float m[3][3], float qin[3]);


// +++++++++++++++++++++++++++++++++++++++++
class           Parameter_Reader {

  public:

    Parameter_Reader();
    // workaround until this class is made a singleton.
    static int      verbosity_level();

    void            initialize(int argc, char **argv);
    bool            parameter_input_successful();
    std::string     query_param(std::string name, std::string default_value,
                                std::string legal_values = "");
    void            read_params();
    void            write_params();

  private:

    // for consistent I/O formatting of name/value input pairs
    static const int         NAME_WIDTH;
    // workaround until this class is made a singleton.
    static const int         UNDEFINED_VERBOSITY;
    static int      verbosity;

    bool            no_stdin;
    int             num_missing_params;
    std::string     param_file_name;

    typedef std::map< std::string, std::string > InputPairsList;
    InputPairsList params_in;

    typedef std::vector < std::pair< std::string, std::string > >
        ListOfUsedInputPairs;
    ListOfUsedInputPairs params_out;

};

int             check_commandline_flag(char **, int, const char *);
int             check_commandline_argument(char **, int, const char *);
std::string     parse_commandline_argument(char **, int, const char *);

void            Tokenizer(std::string , std::vector < std::string > & , char );
float           check_neg_angle(float v1[3],float v2[3],float M[3][3]);
void            rotate_on_x_axis(DOCKMol &, int, int, float);


#endif  // UTILS_H

