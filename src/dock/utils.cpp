#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <sstream>
#include "utils.h"

using namespace std;


// +++++++++++++++++++++++++++++++++++++++++
int
check_commandline_flag(char **argv, int argc, const char * flag)
{
    int             i;

    for (i = 0; i < argc; i++) {
        if (argv[i] != NULL) {
            if (strcmp(argv[i], flag) == 0) {
                if (i < (argc - 1)) {
                    if (argv[i + 1][0] == '-')
                        return i;
                } else
                    return i;
            }
        }
    }

    return -1;
}

// +++++++++++++++++++++++++++++++++++++++++
int
check_commandline_argument(char **argv, int argc, const char * flag)
{
    int             i;

    for (i = 0; i < argc; i++) {
        if (argv[i] != NULL) {
            if (strcmp(argv[i], flag) == 0) {
                if ((i < (argc - 1)) && (argv[i + 1][0] != '-'))
                    return i;
            }
        }
    }

    return -1;
}

// +++++++++++++++++++++++++++++++++++++++++
string
parse_commandline_argument(char **argv, int argc, const char * flag)
{
    int             i;
    string          value;

    value = "";
    for (i = 0; i < argc; i++) {
        if (argv[i] != NULL) {
            if (strcmp(argv[i], flag) == 0) {
                value.append(argv[i + 1]);
            }
        }
    }
    return value;
}


// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
DOCKVector::operator=(const DOCKVector & v)
{

    x = v.x;
    y = v.y;
    z = v.z;

    return *this;

}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector & DOCKVector::operator+=(const DOCKVector & v)
{

    x += v.x;
    y += v.y;
    z += v.z;

    return *this;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector & DOCKVector::operator+=(const float &f)
{

    x += f;
    y += f;
    z += f;

    return *this;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector & DOCKVector::operator-=(const DOCKVector & v)
{

    x -= v.x;
    y -= v.y;
    z -= v.z;

    return *this;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector & DOCKVector::operator-=(const float &f)
{

    x -= f;
    y -= f;
    z -= f;

    return *this;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector & DOCKVector::operator*=(const float &f)
{

    x *= f;
    y *= f;
    z *= f;

    return *this;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector & DOCKVector::operator/=(const float &f)
{

    x /= f;
    y /= f;
    z /= f;

    return *this;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector & DOCKVector::normalize_vector()
{
    float
        len;

    len = length();

    if (len != 0) {
        x = x / len;
        y = y / len;
        z = z / len;
    }

    return *this;
}

// +++++++++++++++++++++++++++++++++++++++++
float
DOCKVector::length() const
{
    float           len;

    len = sqrt(x * x + y * y + z * z);

    return len;
}

// +++++++++++++++++++++++++++++++++++++++++
float
DOCKVector::squared_length() const
{
    float           len;

    len = x * x + y * y + z * z;

    return len;
}

// +++++++++++++++++++++++++++++++++++++++++
float
DOCKVector::squared_dist(const DOCKVector & v) const
{
    float           dist;

    dist =
        (x - v.x) * (x - v.x) + (y - v.y) * (y - v.y) + (z - v.z) * (z - v.z);

    return dist;
}

// +++++++++++++++++++++++++++++++++++++++++
int
operator==(const DOCKVector & v1, const DOCKVector & v2)
{

    if ((v1.x == v2.x) && (v1.y == v2.y) && (v1.z == v2.z))
        return (true);
    else
        return (false);

}

// +++++++++++++++++++++++++++++++++++++++++
int
operator!=(const DOCKVector & v1, const DOCKVector & v2)
{

    if ((v1.x != v2.x) || (v1.y != v2.y) || (v1.z != v2.z))
        return (true);
    else
        return (false);

}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator+(const DOCKVector & v1, const DOCKVector & v2)
{
    DOCKVector      vec;

    vec.x = v1.x + v2.x;
    vec.y = v1.y + v2.y;
    vec.z = v1.z + v2.z;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator+(const float &f, const DOCKVector & v2)
{
    DOCKVector      vec;

    vec.x = f + v2.x;
    vec.y = f + v2.y;
    vec.z = f + v2.z;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator+(const DOCKVector & v1, const float &f)
{
    DOCKVector      vec;

    vec.x = v1.x + f;
    vec.y = v1.y + f;
    vec.z = v1.z + f;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator-(const DOCKVector & v1, const DOCKVector & v2)
{
    DOCKVector      vec;

    vec.x = v1.x - v2.x;
    vec.y = v1.y - v2.y;
    vec.z = v1.z - v2.z;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator-(const DOCKVector & v)
{
    DOCKVector      vec;

    vec.x = -v.x;
    vec.y = -v.y;
    vec.z = -v.z;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator-(const DOCKVector & v1, const float &f)
{
    DOCKVector      vec;

    vec.x = v1.x - f;
    vec.y = v1.y - f;
    vec.z = v1.z - f;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator*(const DOCKVector & v1, const DOCKVector & v2)
{
    DOCKVector      vec;

    vec.x = v1.x * v2.x;
    vec.y = v1.y * v2.y;
    vec.z = v1.z * v2.z;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator*(const float &f, const DOCKVector & v2)
{
    DOCKVector      vec;

    vec.x = f * v2.x;
    vec.y = f * v2.y;
    vec.z = f * v2.z;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator*(const DOCKVector & v1, const float &f)
{
    DOCKVector      vec;

    vec.x = v1.x * f;
    vec.y = v1.y * f;
    vec.z = v1.z * f;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
operator/(const DOCKVector & v1, const float &f)
{
    DOCKVector      vec;

    vec.x = v1.x / f;
    vec.y = v1.y / f;
    vec.z = v1.z / f;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
float
dot_prod(const DOCKVector & v1, const DOCKVector & v2)
{
    float           dp;

    dp = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;

    return dp;
}

// +++++++++++++++++++++++++++++++++++++++++
DOCKVector
cross_prod(const DOCKVector & v1, const DOCKVector & v2)
{
    DOCKVector      vec;

    vec.x = v1.y * v2.z - v1.z * v2.y;
    vec.y = -v1.x * v2.z + v1.z * v2.x;
    vec.z = v1.x * v2.y - v1.y * v2.x;

    return vec;
}

// +++++++++++++++++++++++++++++++++++++++++
float
get_vector_angle(const DOCKVector & v1, const DOCKVector & v2)
{
    float           mag,
                    prod;
    float           result;

    mag = v1.length() * v2.length();
    prod = dot_prod(v1, v2) / mag;

    if (prod < -0.999999)
        prod = -0.9999999f;

    if (prod > 0.9999999)
        prod = 0.9999999f;

    if (prod > 1.0)
        prod = 1.0f;

    result = (acos(prod) / PI) * 180;

    return result;
}

// +++++++++++++++++++++++++++++++++++++++++
float
get_torsion_angle(DOCKVector & v1, DOCKVector & v2, DOCKVector & v3,
                  DOCKVector & v4)
{
    float           torsion;
    DOCKVector      b1,
                    b2,
                    b3,
                    c1,
                    c2,
                    c3;

    b1 = v1 - v2;
    b2 = v2 - v3;
    b3 = v3 - v4;

    c1 = cross_prod(b1, b2);
    c2 = cross_prod(b2, b3);
    c3 = cross_prod(c1, c2);

    if (c1.length() * c2.length() < 0.001) {
        torsion = 0.0;
    } else {
        torsion = get_vector_angle(c1, c2);
        if (dot_prod(b2, c3) > 0.0)
            torsion *= -1.0;
    }

    return (torsion);
}


// +++++++++++++++++++++++++++++++++++++++++
// static member initializers

// a magic number
const int    Parameter_Reader::NAME_WIDTH = 61;

const int    Parameter_Reader::UNDEFINED_VERBOSITY = -1;
int          Parameter_Reader::verbosity = UNDEFINED_VERBOSITY;


// +++++++++++++++++++++++++++++++++++++++++
Parameter_Reader::Parameter_Reader()
{
    //
}

// +++++++++++++++++++++++++++++++++++++++++
void
Parameter_Reader::initialize(int argc, char **argv)
{

    if (check_commandline_argument(argv, argc, "-i") != -1) {
        param_file_name = parse_commandline_argument(argv, argc, "-i");
        read_params();
    } else {
        cout << endl <<
            "Usage:\n\tdock6 -i filename.in [-o filename.out] [-v]" <<
            endl;
        exit(0);
    }

    verbosity = 0;
    // the verbose flag is used for extra scoring information
    if (check_commandline_flag(argv, argc, "-v") != -1)
        verbosity = 1;

    // LEP - included a debug verbose flag
    if (check_commandline_flag(argv, argc, "-V") != -1)
        verbosity = 2;

    if (check_commandline_argument(argv, argc, "-o") != -1)
        no_stdin = true;
    else
        no_stdin = false;

}

// +++++++++++++++++++++++++++++++++++++++++
// Currently there are only 2 levels: 
// 0 no verbose
// 1 verbose
int
Parameter_Reader::verbosity_level()
{
    if ( verbosity == UNDEFINED_VERBOSITY ) {
        cout << "\nError:  internal program issue processing verbosity! \n"
                "Please report this to dock-fans@docking.org\n"
                "Docking results will not be affected.\n"
             << endl;
    }
    return verbosity;
}

// +++++++++++++++++++++++++++++++++++++++++
// Read all the parameters as name/value pairs from the input file.
// Warn on duplicates and ignore all but the final value.
void
Parameter_Reader::read_params()
{
    // constructor opens the file; destructor closes it.
    ifstream param_file_in( param_file_name.c_str() );
    if ( ! param_file_in ) {
        cout << "ERROR:  Failure opening the input file: "
             << param_file_name.c_str()
             << endl;
        exit(0); 
    }

    string name, value;
    while (param_file_in >> name >> value) {
        InputPairsList::const_iterator p = params_in.find( name );
        if ( p != params_in.end() ) {  // found - duplicate input name
            cout << "Warning: parameter " << name
                 << " is repeated in the input file.\n"
                 << "    Ignoring the previous value of \"" << p->second << "\""
                 << endl; 
        }
        params_in[ name ] = value;
    }

    num_missing_params = 0;
}

// +++++++++++++++++++++++++++++++++++++++++
// Search the list of inputed name/value pairs and validate the value.
// If a legal value is not found then query the user or use the
// default value. 
// Save valid name/value pairs as well as unsed pairs.
string
Parameter_Reader::query_param(string name, string default_value,
                              string legal_values)
{
    bool            found_legal_value = false;
    string          value;

    InputPairsList::iterator p = params_in.find( name );
    if ( p != params_in.end() ) {  // found
        if ( legal_values.empty() ||
                legal_values.find( p->second ) != string::npos ) {  // legal
            found_legal_value = true;
            value = p->second;
            // save this valid name/value pair
            params_out.push_back( make_pair( name, value ) );
            // remove this name/value pair
            params_in.erase( p );
            cout << setw( NAME_WIDTH ) << left << name << value << endl;
        } else {
            cout << "Error: parameter " << name
                 << " has the illegal value:  \"" << p->second << "\"" << endl; 
            // remove this name/value pair
            params_in.erase( p );
        }
    }

    if ( ! found_legal_value ) {
        if ( no_stdin ) {
            value = default_value;
            cout << "Warning:  No legal value found for parameter " << name
                 << ".\n  The default value of \"" << value
                 << "\" will be used." << endl;
            params_out.push_back( make_pair( name, value ) );
        } else {
            cout << setw( NAME_WIDTH ) << name << " [" << default_value
                 << "] (" << legal_values << "):" << flush;

            if (cin.peek() != '\n') {
                cin >> value;
                cin.get();
            } else {
                value = default_value;
                cin.get();
            }

            int prompt_count = 0;
            while ((!legal_values.empty()) && 
                   (legal_values.find(value) == string::npos) &&
                   (prompt_count < 5)) {
                cout << "Illegal value for " << name << ", please re-enter: ";
                cin >> value;
                cin.get();
                cout << endl;
                prompt_count++;
            }

            if ( legal_values.empty() ||
                    legal_values.find(value) != string::npos ) {
                params_out.push_back( make_pair( name, value ) );
            } else {
                ++num_missing_params;
                cout << "Error:  Missing parameter/value pair.  Too many illegal values."
                    << endl;
                cout << name << " [" << default_value << "] (" << legal_values <<
                    ")" << endl;
            }
        }
    }

    return value;
}

// +++++++++++++++++++++++++++++++++++++++++
// Emit the unused input parameters to the dock output file.
// Rewrite the dock input file with the used input parameters.
void
Parameter_Reader::write_params()
{
    cout <<
        "------------------------------------------------------------------------------------------\n";

    // Any remaining name/value pairs were unused.
    InputPairsList::const_iterator p;
    for ( p = params_in.begin(); p != params_in.end(); ++p ) {
        cout << "\nWarning: parameter " << p->first
             << " has not been used." << endl; 
    }

    // constructor opens the file; destructor closes it.
    ofstream param_file_out( param_file_name.c_str() );

    ListOfUsedInputPairs::const_iterator q;
    for ( q = params_out.begin(); q != params_out.end(); ++q ) {
        param_file_out << setw( NAME_WIDTH ) << left << q->first
                       << q->second << endl;
    }

}

// +++++++++++++++++++++++++++++++++++++++++
bool
Parameter_Reader::parameter_input_successful()
{

    if (num_missing_params > 0) {
        cout << endl << endl << "Parameterization Failed." << endl;
        cout << num_missing_params << " missing parameters." << endl;
        return false;
    }

    return true;
}


// +++++++++++++++++++++++++++++++++++++++++
int             get_matrix_from_quaternion(float m[3][3],       /* rotation
                                                                 * matrix */
                                           float qin[3] /* input independent
                                                         * quaternion elements */
    ) {
    int             i;          /* iterator */
    float           qn;         /* dependent quaternion element, q-naught */
    float           qn2;        /* square of q-naught */
    float           q[3];       /* independent quaternion elements */
    float           q2[3];      /* square of independent quaternion elements */
    float           sum2;       /* sum of squares of independent q values */
    float           sum;        /* sum of independent q values */

    float           temp1,
                    temp2;

    /*
     * Check that each q-value is between -1.0 and 1.0.
     * If not, remap into this range using wrap-around.
     * Compute q-squared values and their sum
     * 10/96 te
     */
    for (i = 0, sum2 = 0.0; i < 3; i++) {
        q[i] = qin[i];

        if (q[i] > 1.0)
            q[i] = fmod(q[i] + 1.0, 2.0) - 1.0;

        else if (q[i] < -1.0)
            q[i] = fmod(q[i] - 1.0, 2.0) + 1.0;

        sum2 += q2[i] = (q[i] * q[i]);
    }

    /*
     * If the sum-of-squares is less than 1.0, compute q-naught
     * 10/96 te
     */
    if (sum2 < 1.0) {
        qn2 = 1.0 - sum2;
        qn = sqrt(qn2);
    }

    /*
     * If the sum is 1.0, set q-naught to zero
     * 10/96 te
     */
    else if (sum2 == 1.0)
        qn = qn2 = 0.0;

    /*
     * Otherwise, renormalize q-values and set q-naught to zero
     * 10/96 te
     */
    else {
        for (i = 0; i < 3; i++)
            q2[i] /= sum2;

        for (i = 0, sum = sqrt(sum2); i < 3; i++)
            q[i] /= sum;

        qn = qn2 = 0.0;
    }

    /*
     * Compute rotation matrix elements
     * 10/96 te
     */

    m[0][0] = qn2 + q2[0] - q2[1] - q2[2];

    temp1 = q[0] * q[1];
    temp2 = qn * q[2];
    m[0][1] = 2.0 * (temp1 + temp2);

    temp1 = q[0] * q[2];
    temp2 = qn * q[1];
    m[0][2] = 2.0 * (temp1 - temp2);

    temp1 = q[0] * q[1];
    temp2 = qn * q[2];
    m[1][0] = 2.0 * (temp1 - temp2);

    m[1][1] = qn2 - q2[0] + q2[1] - q2[2];

    temp1 = q[1] * q[2];
    temp2 = qn * q[0];
    m[1][2] = 2.0 * (temp1 + temp2);

    temp1 = q[0] * q[2];
    temp2 = qn * q[1];
    m[2][0] = 2.0 * (temp1 + temp2);

    temp1 = q[1] * q[2];
    temp2 = qn * q[0];
    m[2][1] = 2.0 * (temp1 - temp2);

    m[2][2] = qn2 - q2[0] - q2[1] + q2[2];

    return true;
}

// +++++++++++++++++++++++++++++++++++++++++
// Trent E. Balius
// THE FOLLOWING FUNCTIONS ARE USED FOR CALCULATING DESCRIPTOR SCORE.
// In c there is a strtok() for character arrays, but no equal function exists for strings in c++.
// This function takes a string, will split on white space, and returns a vector of strings.
void Tokenizer(string line, vector < string > & tokens, char splitOn) {
    stringstream ss;
    ss.clear();
    tokens.clear();
    ss << line;
    string temp_token;
    // split on the charicture splitOn
    while(getline(ss,temp_token,splitOn)){
        //check that the string is not empty then append on the list
        if (!temp_token.empty())
            tokens.push_back(temp_token);
    }
    return;
}

//mat[3][3]
// added by Trent Balius used in covalent and de novo. 
float check_neg_angle(float v1[3],float v2[3],float M[3][3]){
// # use a matrix
// # https://en.wikipedia.org/wiki/Cross_product#Computational_geometry
// # let there be 3 points.
// # let v1 be p1->p2, and v2 be p1-p3
// # let M be a rotation matrix that places v1 and v2 in the x,y plan

     cout << "In float Orient::check_neg_angle(float v1[3],float v2[3],float M[3][3])" << endl;


     float P = ( (v1[0]*M[0][0] + v1[1]*M[0][1] + v1[2]*M[0][2]) *
                 (v2[0]*M[1][0] + v2[1]*M[1][1] + v2[2]*M[1][2]) -
                 (v1[0]*M[1][0] + v1[1]*M[1][1] + v1[2]*M[1][2]) *
                 (v2[0]*M[0][0] + v2[1]*M[0][1] + v2[2]*M[0][2]) );
     cout << "Transformed vector1 onto X-Y plan: " << endl;
     cout << "v1_0: " << (v1[0]*M[0][0] + v1[1]*M[0][1] + v1[2]*M[0][2]) << endl;
     cout << "v1_1: " << (v1[0]*M[1][0] + v1[1]*M[1][1] + v1[2]*M[1][2]) << endl;
     cout << "v1_2: " << (v1[0]*M[2][0] + v1[1]*M[2][1] + v1[2]*M[2][2]) << endl;
     cout << "Transformed vector2 onto X-Y plan: " << endl;
     cout << "v2_0: " << (v2[0]*M[0][0] + v2[1]*M[0][1] + v2[2]*M[0][2]) << endl;
     cout << "v2_1: " << (v2[0]*M[1][0] + v2[1]*M[1][1] + v2[2]*M[1][2]) << endl;
     cout << "v2_2: " << (v2[0]*M[2][0] + v2[1]*M[2][1] + v2[2]*M[2][2]) << endl;
     if (P >= 0){
         cout << "angle is positive" << endl;
         return +1.0;
     }
     else{ // P< 0
         cout << "angle is negative" << endl;
         return -1.0;
     }
}

