#include "iso_align.h"
#include "dockmol.h"
#include "hungarian.h"

#include <ostream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>

class A_Pair;
class Scored_Triangle;
class DOCKMol;
class Hungarian_RMSD;
class Domain;
class Score;
class Iso_Parm;


//Iso_Parm methods

void  Iso_Parm::set_bond_angle_tol_sid(float angle){
    this->bond_angle_tol_sid = angle;
}
void  Iso_Parm::set_bond_angle_tol_lnk(float angle){
    this->bond_angle_tol_lnk = angle;
}
void  Iso_Parm::set_bond_angle_tol_scf(float angle){
    this->bond_angle_tol_scf = angle;
}

float Iso_Parm::get_bond_angle_tol_sid(){
    return this->bond_angle_tol_sid;
}
float Iso_Parm::get_bond_angle_tol_lnk(){
    return this->bond_angle_tol_lnk;
}
float Iso_Parm::get_bond_angle_tol_scf(){
    return this->bond_angle_tol_scf;
}

void  Iso_Parm::set_dist_tol_sid(float dist){
    this->dist_tol_sid=dist; 
};
void  Iso_Parm::set_dist_tol_lnk(float dist){
    this->dist_tol_lnk=dist;
}
void  Iso_Parm::set_dist_tol_scf(float dist){
    this->dist_tol_scf=dist;
}

float Iso_Parm::get_dist_tol_sid(){
    return this->dist_tol_sid;
};
float Iso_Parm::get_dist_tol_lnk(){
    return this->dist_tol_lnk;
}
float Iso_Parm::get_dist_tol_scf(){
    return this->dist_tol_scf;
}

//Iso_Parm END

//Iso_Score::Domain methods

Iso_Score::Domain::Domain(){

    this->overlap     = 0.0;
    this->overlap_hvy = 0.0;
    this->overlap_pho = 0.0;
    this->overlap_phi = 0.0;
    this->overlap_neg = 0.0;
    this->overlap_pos = 0.0;

}

float Iso_Score::Domain::get_overlap(){

    return this->overlap;
};

float Iso_Score::Domain::get_overlap_hvy(){
    
    return this->overlap_hvy;
};

float Iso_Score::Domain::get_overlap_pho(){

    return this->overlap_pho;

};
float Iso_Score::Domain::get_overlap_phi(){

    return this->overlap_phi;

};
float Iso_Score::Domain::get_overlap_neg(){

    return this->overlap_neg;
};

float Iso_Score::Domain::get_overlap_pos(){

    return this->overlap_pos;
};

void  Iso_Score::Domain::set_overlap(float overlap){

    this->overlap=overlap;
};
void  Iso_Score::Domain::set_overlap_hvy(float overlap){

    this->overlap_hvy=overlap;
};
void  Iso_Score::Domain::set_overlap_pho(float overlap){

    this->overlap_pho=overlap;

};
void  Iso_Score::Domain::set_overlap_phi(float overlap){

    this->overlap_phi=overlap;
};
void  Iso_Score::Domain::set_overlap_neg(float overlap){

    this->overlap_neg=overlap;
};
void  Iso_Score::Domain::set_overlap_pos(float overlap){

    this->overlap_pos=overlap;
};

//Domain method end


//Iso_Score::Domain methods end


float Iso_Score::Analytical_method(Fragment & reffrag, Fragment & testfrag)
{
    int         i,
                j;
    float       pi = 3.14,
                d_sqr,
                overlap_temp;

    //create three different domains for volume overlap
    Iso_Score::Domain      ref_o;
    Iso_Score::Domain     test_o;
    Iso_Score::Domain ref_test_o;

    //these are the different types of scores within 
    //volume overlap
    float                 total_component=0.0;
    float            heavy_atom_component=0.0;
    float              positive_component=0.0;
    float              negative_component=0.0;
    float           hydrophobic_component=0.0;
    float           hydrophilic_component=0.0;


    //pair wise calculation where you caculate 
    //atomic level volume overlaps
    for (i = 0; i < testfrag.mol.num_atoms; i++){
        if (testfrag.mol.atom_active_flags[i]) {
            for (j = 0; j < reffrag.mol.num_atoms; j++){
                d_sqr = (testfrag.mol.x[i] - reffrag.mol.x[j]) * 
                        (testfrag.mol.x[i] - reffrag.mol.x[j]) + 
                        (testfrag.mol.y[i] - reffrag.mol.y[j]) * 
                        (testfrag.mol.y[i] - reffrag.mol.y[j]) + 
                        (testfrag.mol.z[i] - reffrag.mol.z[j]) * 
                        (testfrag.mol.z[i] - reffrag.mol.z[j]);
                if (d_sqr < (testfrag.mol.amber_at_radius[i] + 
                    reffrag.mol.amber_at_radius[j]) * 
                    (testfrag.mol.amber_at_radius[i]            + 
                    reffrag.mol.amber_at_radius[j])) {
                    if (d_sqr != 0)
                        overlap_temp = pi / 12 * 
                           (testfrag.mol.amber_at_radius[i] + reffrag.mol.amber_at_radius[j] - sqrt(d_sqr)) * 
                           (testfrag.mol.amber_at_radius[i] + reffrag.mol.amber_at_radius[j] - sqrt(d_sqr)) * 
                           (sqrt(d_sqr) + 2 * (testfrag.mol.amber_at_radius[i] + reffrag.mol.amber_at_radius[j]) - 3 / 
                                sqrt(d_sqr) * (testfrag.mol.amber_at_radius[i] - reffrag.mol.amber_at_radius[j]) *
                                (testfrag.mol.amber_at_radius[i] - reffrag.mol.amber_at_radius[j]));
                    else
                        overlap_temp = pi / 12 * 2 * 
                            (testfrag.mol.amber_at_radius[i] + reffrag.mol.amber_at_radius[j]) * 
                            (testfrag.mol.amber_at_radius[i] + reffrag.mol.amber_at_radius[j]) * 
                            (testfrag.mol.amber_at_radius[i] + reffrag.mol.amber_at_radius[j]);

                    ref_test_o.set_overlap(ref_test_o.get_overlap() + overlap_temp);

                    if (testfrag.mol.atom_types[i] != "H" && reffrag.mol.atom_types[j] != "H")


                        ref_test_o.set_overlap_hvy(ref_test_o.get_overlap_hvy() + overlap_temp);

                    if (testfrag.mol.chem_types[i] == "hydrophobic" && reffrag.mol.chem_types[j] == "hydrophobic")

                        ref_test_o.set_overlap_pho(ref_test_o.get_overlap_pho() + overlap_temp);

                    if ((testfrag.mol.chem_types[i] == "donar"    || 
                         testfrag.mol.chem_types[i] == "acceptor" || 
                         testfrag.mol.chem_types[i] == "polar")   && 
                         (reffrag.mol.chem_types[j]   == "donar"    || 
                            reffrag.mol.chem_types[j] == "acceptor" || 
                            reffrag.mol.chem_types[j] == "polar"))

                        ref_test_o.set_overlap_phi(ref_test_o.get_overlap_phi() + overlap_temp);

                    if (testfrag.mol.charges[i] < 0 && reffrag.mol.charges[j] < 0)

                        ref_test_o.set_overlap_neg(ref_test_o.get_overlap_neg() + overlap_temp);

                    if (testfrag.mol.charges[i] > 0 && reffrag.mol.charges[j] > 0)

                        ref_test_o.set_overlap_pos(ref_test_o.get_overlap_pos() + overlap_temp);

                }
            }
        }
    }

    for (i = 0; i < testfrag.mol.num_atoms; i++){
        if (testfrag.mol.atom_active_flags[i]) {
            for (j = 0; j < testfrag.mol.num_atoms; j++){
                if (testfrag.mol.atom_active_flags[j]) {
                    d_sqr = (testfrag.mol.x[i] - testfrag.mol.x[j]) * 
                            (testfrag.mol.x[i] - testfrag.mol.x[j]) + 
                            (testfrag.mol.y[i] - testfrag.mol.y[j]) * 
                            (testfrag.mol.y[i] - testfrag.mol.y[j]) + 
                            (testfrag.mol.z[i] - testfrag.mol.z[j]) * 
                            (testfrag.mol.z[i] - testfrag.mol.z[j]);
                    if (d_sqr < (testfrag.mol.amber_at_radius[i] + testfrag.mol.amber_at_radius[j]) *
                        (testfrag.mol.amber_at_radius[i] + testfrag.mol.amber_at_radius[j])) {
                        if (d_sqr != 0)
                            overlap_temp = pi / 12 * 
                                (testfrag.mol.amber_at_radius[i] + testfrag.mol.amber_at_radius[j] - sqrt(d_sqr)) * 
                                (testfrag.mol.amber_at_radius[i] + testfrag.mol.amber_at_radius[j] - sqrt(d_sqr)) * 
                                (sqrt(d_sqr) + 2 * (testfrag.mol.amber_at_radius[i] + testfrag.mol.amber_at_radius[j]) - 3 / 
                                 sqrt(d_sqr) * (testfrag.mol.amber_at_radius[i] - testfrag.mol.amber_at_radius[j])       * 
                                (testfrag.mol.amber_at_radius[i] - testfrag.mol.amber_at_radius[j]));
                        else
                            overlap_temp = pi / 12 * 2 * 
                                (testfrag.mol.amber_at_radius[i] + testfrag.mol.amber_at_radius[j]) * 
                                (testfrag.mol.amber_at_radius[i] + testfrag.mol.amber_at_radius[j]) *  
                                (testfrag.mol.amber_at_radius[i] + testfrag.mol.amber_at_radius[j]);

                        test_o.set_overlap(test_o.get_overlap() + overlap_temp);
                        if (testfrag.mol.atom_types[i] != "H" && testfrag.mol.atom_types[j] != "H")

                            test_o.set_overlap_hvy(test_o.get_overlap_hvy() + overlap_temp);

                        if (testfrag.mol.chem_types[i] == "hydrophobic" && testfrag.mol.chem_types[j] == "hydrophobic")
                             test_o.set_overlap_pho(test_o.get_overlap_pho() + overlap_temp); 
                        if ((testfrag.mol.chem_types[i] == "donar"    || 
                             testfrag.mol.chem_types[i] == "acceptor" || 
                             testfrag.mol.chem_types[i] == "polar")   && 
                            (testfrag.mol.chem_types[j] == "donar"    ||  
                             testfrag.mol.chem_types[j] == "acceptor" || 
                             testfrag.mol.chem_types[j] == "polar"))

                            test_o.set_overlap_phi(test_o.get_overlap_phi() + overlap_temp);
                        if (testfrag.mol.charges[i] < 0 && testfrag.mol.charges[j] < 0)

                            test_o.set_overlap_neg(test_o.get_overlap_neg() + overlap_temp);  

                        if (testfrag.mol.charges[i] > 0 && testfrag.mol.charges[j] > 0)

                            test_o.set_overlap_pos(test_o.get_overlap_pos() + overlap_temp);
                    }
                }
            }
        }
    }



    if (test_o.get_overlap() == 0 && ref_o.get_overlap() == 0)
        total_component = 0;

    else
        total_component = ref_test_o.get_overlap() / max(test_o.get_overlap(), ref_o.get_overlap() );

    if (test_o.get_overlap_hvy() == 0 && ref_o.get_overlap_hvy() == 0)
        heavy_atom_component = 0; 
    else
        heavy_atom_component = ref_test_o.get_overlap_hvy() / max(test_o.get_overlap_hvy(), ref_o.get_overlap_hvy() );

    if (test_o.get_overlap_pho() == 0 && ref_o.get_overlap_pho() == 0)
        hydrophobic_component = 0; 
    else
        hydrophobic_component = ref_test_o.get_overlap_pho() / max(test_o.get_overlap_pho(), ref_o.get_overlap_pho() );

    if (test_o.get_overlap_phi() == 0 && ref_o.get_overlap_phi() == 0)
        hydrophilic_component = 0;
    else
        hydrophilic_component = ref_test_o.get_overlap_phi() / max(test_o.get_overlap_phi(), ref_o.get_overlap_phi() );

    if (test_o.get_overlap_neg() == 0 && ref_o.get_overlap_neg() == 0)
        negative_component = 0;
    else
        negative_component = ref_test_o.get_overlap_neg() / max(test_o.get_overlap_neg(), ref_o.get_overlap_neg() );

    if (test_o.get_overlap_pos() == 0 && ref_o.get_overlap_pos() == 0)
        positive_component = 0;
    else
        positive_component = ref_test_o.get_overlap_pos() / max(test_o.get_overlap_pos(), ref_o.get_overlap_pos() );

    return total_component;

}





//Iso_Score class Methods


void Iso_Score::Score::set_vo_score(float score){

    this->vo_score = score;
}

void Iso_Score::Score::set_score(float score){

    this->score = score;
}

float Iso_Score::Score::get_score(){

    return this->score;
}

float Iso_Score::Score::get_vo_score(){

    return this->vo_score;
}
void Iso_Score::Score::calc_score(Fragment& reffrag,Fragment& testfrag){
 
    this->score=Iso_Score::Analytical_method(reffrag,testfrag);
}


//Iso_Score class Methods END


//A_Pair class methods

int Iso_Acessory::A_Pair::get_ref(){
    int ref_atomnum = (*this).first;
    return ref_atomnum;
}

int Iso_Acessory::A_Pair::get_test(){
    int test_atomnum = (*this).second;
    return test_atomnum;
};

void Iso_Acessory::A_Pair::set_ref(int atom_num){
    (*this).first = atom_num;
};

void Iso_Acessory::A_Pair::set_test(int atom_num){
    (*this).second = atom_num;
};

std::string Iso_Acessory::A_Pair::get_title(){
    std::string name;
    name = this->title;
    return name;
};
void Iso_Acessory::A_Pair::set_title(std::string name){
    this->title = name;
};

void Iso_Acessory::A_Pair::print_all(){

    float score = (*this).get_cos_sim(); 
    int ref =  (*this).get_ref();
    int test = (*this).get_test();
    std::cout << " " << score << " "<< ref << " " 
                    << test << std::endl;
};

void Iso_Acessory::A_Pair::print_all(Fragment& reffrag,Fragment& testfrag){

    float score = (*this).get_cos_sim(); 
    int ref =  (*this).get_ref();
    int test = (*this).get_test();
    std::string ref_atom = reffrag.mol.atom_names[ref];
    std::string test_atom =testfrag.mol.atom_names[test];

    std::cout << " " << score << " "<< ref_atom << " " 
                << test_atom << std::endl;

  
};
void Iso_Acessory::A_Pair::print_ref(){
    int ref = get_ref();
    std::cout << ref <<std::endl;
  
};

void Iso_Acessory::A_Pair::print_test(){
    int test = get_test();
    std::cout << test <<std::endl;
  
};
float Iso_Acessory::A_Pair::get_cos_sim(){
    float score = this->cos_sim;
    return score;
};
void Iso_Acessory::A_Pair::set_cos_sim(float score){
    this->cos_sim = score;
}

void Iso_Acessory::A_Pair::clear(){
    (*this).first   =  0;
    (*this).second  =  0;
    this->cos_sim   = 0.0;
}

//A_Pair class methods END





//Scored_Triangle class methods


Iso_Acessory::Scored_Triangle::Scored_Triangle(){
  
    this->alignable = false;

}

void Iso_Acessory::Scored_Triangle::push_back(A_Pair& one_pair){

    (*this).three_pairs.push_back(one_pair);

    if ((*this).three_pairs.size()>3){
        std::cout << "A Scored_Triangle is reached its limit:"
                  << " has more than 3 elements in its stack"
                  << ", please float check your code so it"
                  << " doesn't beyond 3." << std::endl;
        exit(0);
    }

};

void Iso_Acessory::Scored_Triangle::pop_back(){

    if ((*this).three_pairs.size()>3){
        std::cout << "A Scored_Triangle is reached its limit:"
                  << " has more than 3 elements in its stack"
                  << ", please float check your code so it"
                  << " doesn't beyond 3." << std::endl;
        exit(0);
    }

    (*this).three_pairs.pop_back();
};

void Iso_Acessory::Scored_Triangle::print(){
    for (unsigned int i =0; i<(*this).three_pairs.size();i++){
        float cos = (*this).three_pairs[i].get_ref();
        int   ref = (*this).three_pairs[i].get_ref();
        int   test= (*this).three_pairs[i].get_test();
        std::cout << cos << " "<< ref << " " << test << std::endl;
    }
};


void Iso_Acessory::Scored_Triangle::print(Fragment& reffrag,Fragment& testfrag){

    for (unsigned int i =0; i<(*this).three_pairs.size();i++){
        float cos = (*this).three_pairs[i].get_cos_sim();
        int   ref = (*this).three_pairs[i].get_ref();
        int   test= (*this).three_pairs[i].get_test();
        std::string ref_atom = reffrag.mol.atom_names[ref];
        std::string test_atom = testfrag.mol.atom_names[test];
        std::cout << cos << " "<< ref_atom << " " << test_atom << std::endl;
    }   
};



void Iso_Acessory::Scored_Triangle::clear(){
    (*this).three_pairs.clear();

}

void Iso_Acessory::Scored_Triangle::is_not_alignable(){

    this->alignable = false;

}

void Iso_Acessory::Scored_Triangle::is_alignable(){

    this->alignable = true;

}


//drive function to sort the vector elements by 
//the first element of pairs in descending order
bool Iso_Acessory::Scored_Triangle::check_if_redundant (){
//ENTER HERE 
//  start assuming its false. only need one redudnat atom type to be true. it will return as TRUE 
//  to say "HEY THERE IS REDUNDANCY"
    bool tmp_tf = false;

    for (unsigned int i = 0; i < (*this).three_pairs.size(); i++){
        for (unsigned int z = i+1; z < (*this).three_pairs.size(); z++){    
    
            //This checks if there are redundancies
            if ((*this).three_pairs[i].get_ref()  == (*this).three_pairs[z].get_ref())   {tmp_tf=true;} 
            if ((*this).three_pairs[i].get_test() == (*this).three_pairs[z].get_test()) {tmp_tf=true;}
        }
    }   
    return tmp_tf;
}


//Scored_Triangle class methods END

bool check_H_vs_atom (Fragment & reffrag, Fragment & testfrag){
    int ref_num_H      =  0; 
    int ref_num_heavy  =  0;
    int test_num_H     =  0;
    int test_num_heavy =  0; 
 
    for (unsigned int i = 0; i < reffrag.mol.num_atoms; i++){
        if (reffrag.mol.atom_types[i] == "H"){ref_num_H++;}
        else {ref_num_heavy++;}
    }
   
    for (unsigned int i = 0; i < testfrag.mol.num_atoms; i++) {
        if (testfrag.mol.atom_types[i] == "H") {test_num_H++;}
        else {test_num_heavy++;} 
    }
    
    if (ref_num_H>=ref_num_heavy || test_num_H>=test_num_heavy || reffrag.mol.num_atoms==3 || testfrag.mol.num_atoms==3) {return true;} 
    return false; 
}

float dot_product(std::vector<float> vec1,std::vector<float>vec2) {
    float dot_result = 0.0;
    int vec1_2_size  = 0;
    if (vec1.size() == vec2.size()){
        vec1_2_size = vec1.size();
    }else{
        std::cout<<"dot_product vector sizes are not the same"<<std::endl;
        
    }
    for (int i = 0; i<vec1_2_size; i++){
       dot_result += (vec1[i] * vec2[i]);
    }
    return dot_result;
}

//calculates cross_prod
std::vector<float> cross_prod(const std::vector<float> & v1, const std::vector<float> & v2) 
{
    std::vector <float> vec(3,0.0);

    vec[0] = v1[1] * v2[2] - v1[2] * v2[1];
    vec[1] = -v1[0] * v2[2] + v1[2] * v2[0];
    vec[2] = v1[0] * v2[1] - v1[1] * v2[0];

    return vec;
}

//calcualtes the length of a vector
float length_vec(std::vector<float> len_vec)
{
    float           len;

    len = sqrt( ( len_vec[0] * len_vec[0] )   
             +  ( len_vec[1] * len_vec[1] )
             +  ( len_vec[2] * len_vec[2] ));  
    return len;
}

//substracts a vector
std::vector<float> subtract_vec(const std::vector<float> &v1, const std::vector<float> &v2)
{
    std::vector<float> diff {};

    diff.push_back(v1[0] - v2[0]);
    diff.push_back(v1[1] - v2[1]);
    diff.push_back(v1[2] - v2[2]);

    return diff;
}

//you have to pass through by reference or the molecules will be distorted
void normalize_vec(std::vector <float> &vec_to_norm)
{
    float len;
    len = length_vec(vec_to_norm);
     
    if (len != 0) {
        vec_to_norm[0] = vec_to_norm[0] / len;
        vec_to_norm[1] = vec_to_norm[1] / len;
        vec_to_norm[2] = vec_to_norm[2] / len;
    }
}

//gets magnitude in a vector
float magnitude(std::vector<float> mag_vec){
    float mag_result = 0.0;
    for (unsigned int i =0; i<mag_vec.size(); i++){
        mag_result += (mag_vec[i] * mag_vec[i]);
    }

    return mag_result;
}

float get_torsion_angle(const std::vector<float> v1
                        ,const std::vector<float> v2
                        ,const std::vector<float> v3
                        ,const std::vector<float> v4){

    float           torsion;
    std::vector<float>      b1,
                            b2,
                            b3,
                            c1,
                            c2,
                            c3;
    b1 = subtract_vec(v1,v2);
    b2 = subtract_vec(v2,v3);
    b3 = subtract_vec(v3,v4);
   
    c1 = cross_prod(b1, b2);
    c2 = cross_prod(b2, b3);
    c3 = cross_prod(c1, c2);
    
    if (length_vec(c1) * length_vec(c2) < 0.001) {
        torsion = 0.0;
    } else {
        torsion = get_vector_angle(c1, c2);
        if (dot_product(b2, c3) > 0.0)
            torsion *= -1.0;
    }
    return torsion;

}

float get_vector_angle(const std::vector<float> & v1, const std::vector<float> & v2)
{
    float           mag,
                    prod;
    float           result;
    float           pi = 3.1415926536;

    mag = length_vec(v1) * length_vec(v2);
    prod = dot_product(v1, v2) / mag;

    if (prod < -0.999999)
        prod = -0.9999999f;

    if (prod > 0.9999999)
        prod = 0.9999999f;

    if (prod > 1.0)
        prod = 1.0f;

    result = (acos(prod) / pi) * 180;

    return result;
}


float calc_cos_similarity( Fragment &reffrag, Fragment & testfrag,int atomref, int atomtest){

    float results = 0;
    std::vector <float> ref_fd_vector  = {};
    std::vector <float> test_fd_vector = {};
    ref_fd_vector = reffrag.get_radial_dist_distri(atomref);
    test_fd_vector = testfrag.get_radial_dist_distri(atomtest); 

    //portions of the cosine similirity equation.
    float dotproduct = 0.00;
    float magref = 0.00;
    float magtest = 0.00;
    float magreftest = 0.00; 

    if (ref_fd_vector.size() == test_fd_vector.size()){
        dotproduct = dot_product(ref_fd_vector,test_fd_vector); 
        magref  = magnitude(ref_fd_vector);
        magtest = magnitude(test_fd_vector);
     
        magreftest = sqrt(magref) * sqrt(magtest); 

    }else{
        std::cout << "error in dp calculations, both cosine similarity vectors are not the same size" << std::endl;
    }
    
    //calculate cosine similarity
    results = dotproduct / magreftest; 

    return results;

};


bool 
get_diff_dist_two_mol (Fragment & reffrag, Fragment & testfrag,Iso_Acessory::Scored_Triangle triangle,float cutoff,bool sidechains_or_scaffold)
{

    bool diff         = false;
    float ref_0_1_dist  = 0.0; 
    float ref_0_2_dist  = 0.0;
    float ref_1_2_dist  = 0.0;

    float test_0_1_dist = 0.0;
    float test_0_2_dist = 0.0;
    float test_1_2_dist = 0.0;



    ref_0_1_dist = calc_2atoms_length(reffrag
                                     ,triangle.three_pairs[0].get_ref()
                                     ,triangle.three_pairs[1].get_ref());

    test_0_1_dist = calc_2atoms_length(testfrag
                                     ,triangle.three_pairs[0].get_test()
                                     ,triangle.three_pairs[1].get_test());

    //if sidechain or scaffold
    //calculate distances between atoms between two molecules (ref vs test)
    if (sidechains_or_scaffold == false){

        ref_0_2_dist = calc_2atoms_length(reffrag
                                         ,triangle.three_pairs[0].get_ref()
                                         ,triangle.three_pairs[2].get_ref());

        ref_1_2_dist = calc_2atoms_length(reffrag
                                         ,triangle.three_pairs[1].get_ref()
                                         ,triangle.three_pairs[2].get_ref());

        test_0_2_dist = calc_2atoms_length(testfrag
                                         ,triangle.three_pairs[0].get_test()
                                         ,triangle.three_pairs[2].get_test());
        test_1_2_dist = calc_2atoms_length(testfrag
                                         ,triangle.three_pairs[1].get_test()
                                         ,triangle.three_pairs[2].get_test());
    }

    //std::cout << "ref_0_1_dist: " << ref_0_1_dist << " ref_0_2_dist: " << ref_0_2_dist << " ref_1_2_dist: " << ref_1_2_dist << std::endl;
    //std::cout << "test_0_1_dist: " << test_0_1_dist << " test_0_2_dist: " << test_0_2_dist << " test_1_2_dist: " << test_1_2_dist << std::endl;

   
    //at anypoint, if the distances between two atoms between two molecules are more than
    //cutoff, return true.
    if (sidechains_or_scaffold ==false){ 
        if (abs(ref_0_1_dist - test_0_1_dist) > cutoff){
            diff = true;
        }

        if (abs(ref_0_2_dist - test_0_2_dist) > cutoff){
            diff = true;
        } 
        if (abs(ref_1_2_dist - test_1_2_dist) > cutoff){
            diff = true;
        }
    //if linker, same concept, but apply to du du distances 
    }else{
        if (abs(ref_0_1_dist - test_0_1_dist) > cutoff){
            diff = true;
        }
    }

    return diff; 
}

//local function that gets centroids: 𝑥=Σ𝑥𝑖/𝑛,𝑦=Σ𝑦𝑖/𝑛,z=Σz𝑖/𝑛
std::vector<std::vector<float>> get_centroid(DOCKMol & refmol
                                            ,DOCKMol & testmol
                                            ,Iso_Acessory::Scored_Triangle tri_atom_pairs
                                            ,bool getref){

    std::vector<std::vector<float>> return_vec {};
    
    float x_val_ref  = 0.0;
    float y_val_ref  = 0.0;
    float z_val_ref  = 0.0; 

    float x_val_test = 0.0;
    float y_val_test = 0.0;
    float z_val_test = 0.0;

    if(getref){
        //get ref molecule centroid
        x_val_ref = (refmol.x[tri_atom_pairs.three_pairs[0].get_ref()] + 
                     refmol.x[tri_atom_pairs.three_pairs[1].get_ref()] +
                     refmol.x[tri_atom_pairs.three_pairs[2].get_ref()])/3; 

        y_val_ref = (refmol.y[tri_atom_pairs.three_pairs[0].get_ref()] + 
                     refmol.y[tri_atom_pairs.three_pairs[1].get_ref()] +
                     refmol.y[tri_atom_pairs.three_pairs[2].get_ref()])/3;

        z_val_ref = (refmol.z[tri_atom_pairs.three_pairs[0].get_ref()] +
                     refmol.z[tri_atom_pairs.three_pairs[1].get_ref()] +
                     refmol.z[tri_atom_pairs.three_pairs[2].get_ref()])/3; 
        
        return_vec.push_back({x_val_ref,y_val_ref,z_val_ref});
    }

    //get testing molecule centroid 
    x_val_test = (testmol.x[tri_atom_pairs.three_pairs[0].get_test()] + 
		  testmol.x[tri_atom_pairs.three_pairs[1].get_test()] +
		  testmol.x[tri_atom_pairs.three_pairs[2].get_test()])/3; 

    y_val_test = (testmol.y[tri_atom_pairs.three_pairs[0].get_test()] + 
                  testmol.y[tri_atom_pairs.three_pairs[1].get_test()] +
                  testmol.y[tri_atom_pairs.three_pairs[2].get_test()])/3; 

    z_val_test = (testmol.z[tri_atom_pairs.three_pairs[0].get_test()] + 
                  testmol.z[tri_atom_pairs.three_pairs[1].get_test()] +
                  testmol.z[tri_atom_pairs.three_pairs[2].get_test()])/3; 

    return_vec.push_back({x_val_test,y_val_test,z_val_test});

    return return_vec;
}

//translates the molecules by the target vector
void target_translation (DOCKMol & trans_mol,std::vector<float> trans_vec) {

    for (unsigned int i = 0; i<trans_mol.num_atoms; i++){
        trans_mol.x[i]+=trans_vec[0];
        trans_mol.y[i]+=trans_vec[1];
        trans_mol.z[i]+=trans_vec[2];
    }
}

//calculates and returns rotation matrix
//that swivels the moving vec to the ref vec
std::vector<std::vector<double>> get_rotation_mat( std::vector<float>& vec1_test, std::vector<float>& vec2) {
    std::vector<std::vector<double>> final_mat_vec;

    double dot_test        = 0.0; 
    double vec1_magsq_test = 0.0; 
    double vec2_magsq      = 0.0; 

    double cos_theta_test  = 0.0; 
    double sin_theta_test  = 0.0; 
   
    dot_test        =dot_product(vec1_test,vec2);

    vec1_magsq_test =magnitude(vec1_test);
    vec2_magsq      =magnitude(vec2); 

    cos_theta_test = dot_test/ (sqrt(vec1_magsq_test* vec2_magsq));
    if (cos_theta_test >= 1.0){
        cos_theta_test = 1.0;
        sin_theta_test = 0.0;
    }
    else if (cos_theta_test <= -1.0){
        cos_theta_test = -1.0;
        sin_theta_test = 0.0;
    }
    else{
        sin_theta_test = sqrt (1 - (cos_theta_test * cos_theta_test));
    }
    //// if cos_theta is 1, vec 1 and vec 2 are perfect,
    //// otherwise enter this loop and calculate lots of stuff
    //else if (cos_theta_ref !=1 ||  cos_theta_test !=1) {
        std::vector<float> normalU_test = cross_prod(vec1_test,vec2);
    //    
        std::vector<float> normalW_test = cross_prod(vec2,normalU_test);
    //    //normalize vectors
        normalize_vec(vec2);

        normalize_vec(normalU_test);

        normalize_vec(normalW_test);

    //    // Make coorinate rotation matrix, which rotates coordinates to normalW, vec2, normalU coordinates    
        float coorRot_test[3][3];


        coorRot_test[0][0]    = normalW_test[0];
        coorRot_test[0][1]    =         vec2[0];
        coorRot_test[0][2]    = normalU_test[0];
        coorRot_test[1][0]    = normalW_test[1];
        coorRot_test[1][1]    =         vec2[1];
        coorRot_test[1][2]    = normalU_test[1];
        coorRot_test[2][0]    = normalW_test[2];
        coorRot_test[2][1]    =         vec2[2];
        coorRot_test[2][2]    = normalU_test[2];

    //    // make inverse matrix of coordRot matrix - since coorRot is an orthonal matrix 
    //    // the inverse if its transpose, coorRot^T
        float invcoorRot_test[3][3];

        invcoorRot_test[0][0] = coorRot_test[0][0];
        invcoorRot_test[0][1] = coorRot_test[1][0];
        invcoorRot_test[0][2] = coorRot_test[2][0];
        invcoorRot_test[1][0] = coorRot_test[0][1];
        invcoorRot_test[1][1] = coorRot_test[1][1];
        invcoorRot_test[1][2] = coorRot_test[2][1];
        invcoorRot_test[2][0] = coorRot_test[0][2];
        invcoorRot_test[2][1] = coorRot_test[1][2];
        invcoorRot_test[2][2] = coorRot_test[2][2];

        float planeRot_test[3][3];

        planeRot_test[0][0] =  cos_theta_test; planeRot_test[0][1] = -sin_theta_test; planeRot_test[0][2] = 0;
        planeRot_test[1][0] =  sin_theta_test; planeRot_test[1][1] =  cos_theta_test; planeRot_test[1][2] = 0;
        planeRot_test[2][0] =               0; planeRot_test[2][1] =               0; planeRot_test[2][2] = 1;


        float temp_test[3][3];

        double finalmat_test[3][3];

        for (int w=0; w<3; w++){
            for (int q=0; q<3; q++){
                temp_test[w][q] = 0.0;
                for (int k=0; k<3; k++){
                    temp_test[w][q] += coorRot_test[w][k]*planeRot_test[k][q];
                }
            }
        }

        for (int w=0; w<3; w++){
            for (int q=0; q<3; q++){
                finalmat_test[w][q] = 0.0;
                for (int k=0; k<3; k++){
                    finalmat_test[w][q] += temp_test[w][k]*invcoorRot_test[k][q];
                }
            }
        }
        std::vector <double> t_vec {};
        for (int w=0; w<3; w++){
            for (int q=0; q<3; q++){
                t_vec.push_back(finalmat_test[w][q]);
            }
            final_mat_vec.push_back(t_vec);
            t_vec.clear();
        }

    return final_mat_vec;
};

void align_molcentroid (DOCKMol& refmol, DOCKMol& testmol,std::vector<std::vector<float>> centroid_vec,Iso_Acessory::Scored_Triangle three_atom_pairs){


    //instantiate vecotrs that are used to translate mols
    //to the origin
    std::vector<float> test_trans_vec{};
    std::vector<float> ref_trans_vec{};

    //translate centroid coordinates for coresponding
    //refmol and testmol to origin
    ref_trans_vec.push_back(-centroid_vec[0][0]);
    ref_trans_vec.push_back(-centroid_vec[0][1]);
    ref_trans_vec.push_back(-centroid_vec[0][2]);  

    test_trans_vec.push_back(-centroid_vec[1][0]);
    test_trans_vec.push_back(-centroid_vec[1][1]);
    test_trans_vec.push_back(-centroid_vec[1][2]); 
 
    target_translation(refmol,ref_trans_vec); 
    target_translation(testmol,test_trans_vec);
       

    //this block of code tries to move one of the axis
    //of testmol to an axis of the refmol.

    ////get the single point of the triangle that we want to transform 
    std::vector<float> vec1_test{testmol.x[three_atom_pairs.three_pairs[0].get_test()],
                                 testmol.y[three_atom_pairs.three_pairs[0].get_test()],
                                 testmol.z[three_atom_pairs.three_pairs[0].get_test()]};

    ////Where I want to transform to from the test frag
    std::vector<float> vec2     {refmol.x[three_atom_pairs.three_pairs[0].get_ref()],
                                 refmol.y[three_atom_pairs.three_pairs[0].get_ref()],
                                 refmol.z[three_atom_pairs.three_pairs[0].get_ref()]};
    double finalmat_test[3][3];
    std::vector<std::vector <double>> tmp_rot_mat{};
    tmp_rot_mat = get_rotation_mat(vec1_test,vec2);

    //converts the tmp_rot_mat into [3][3] matrix
    for (int w=0; w<3; w++) {
        for (int q=0; q<3; q++) {
            finalmat_test[w][q] = tmp_rot_mat[w][q];
        }
    }

    ////This is to move test frag to the vec2 
    testmol.rotate_mol(finalmat_test);

    tmp_rot_mat.clear();
    vec1_test.clear();
    vec2.clear();

    //Now we have to move both fragment to the X axis to find out how much rotation to align both 
    //triangles of test and ref frag
    vec1_test = {refmol.x[three_atom_pairs.three_pairs[0].get_ref()],
                 refmol.y[three_atom_pairs.three_pairs[0].get_ref()],
                 refmol.z[three_atom_pairs.three_pairs[0].get_ref()]};   
   
    vec2 = {1,0,0};
    double finalmat_test_1[3][3];
    tmp_rot_mat = get_rotation_mat(vec1_test,vec2);

    //converts the tmp_rot_mat into [3][3] matrix
    for (int w=0; w<3; w++) {
        for (int q=0; q<3; q++) {
            finalmat_test_1[w][q] = tmp_rot_mat[w][q];
        }    
    }  
    testmol.rotate_mol(finalmat_test_1);
    refmol.rotate_mol(finalmat_test_1);

    tmp_rot_mat.clear();

    //Here we are going to calculate the torsion angle between the test frag and test frag 
    //Assuming one of median of the triangble is aligned with the X axis.

    //the first point
    std::vector <float>tr1 {refmol.x[three_atom_pairs.three_pairs[1].get_ref()],
                            refmol.y[three_atom_pairs.three_pairs[1].get_ref()],
                            refmol.z[three_atom_pairs.three_pairs[1].get_ref()]};
    //the second point  
    std::vector <float>tr2 {0,0,0};
    //the centroid point 
    std::vector <float>tr3 {1,0,0};
    //the ref test point
    std::vector <float>tr4 {testmol.x[three_atom_pairs.three_pairs[1].get_test()],
                            testmol.y[three_atom_pairs.three_pairs[1].get_test()],
                            testmol.z[three_atom_pairs.three_pairs[1].get_test()]};


    //this is an illustration of what a torsion environment looks like 
    /*
    ILLUSTRATION OF WHAT A TORSION ENVIRONMENT LOOKS LIKE
    //tr1(ref frag) - tr2 (0,0,0) -- tr3(1,0,0) - tr4(test frag)
    */


    //get tor_angle between the fragment.
    float degree_ang; 
    degree_ang = get_torsion_angle(tr1,tr2,tr3,tr4);
    double rad_ang = degree_ang*3.14159/180;
   
    //calculate cos_theta and sin_theta 
    double cos_theta = cos(rad_ang);
    double sin_theta = sin(rad_ang);

    double x_axis_rot[3][3]; 
  
    //X_axis rotation clockwise 
    x_axis_rot[0][0] =  1; x_axis_rot[0][1] =           0; x_axis_rot[0][2] =          0; 
    x_axis_rot[1][0] =  0; x_axis_rot[1][1] =   cos_theta; x_axis_rot[1][2] =  sin_theta; 
    x_axis_rot[2][0] =  0; x_axis_rot[2][1] =  -sin_theta; x_axis_rot[2][2] =  cos_theta; 
  
    //rotation the mol x axis 
    testmol.rotate_mol(x_axis_rot);


    //transform back to where it came from
    double finalmat_test_2[3][3];
    tmp_rot_mat = get_rotation_mat(vec2,vec1_test);
    
    //converts the tmp_rot_mat into [3][3] matrix
    for (int w=0; w<3; w++) {
        for (int q=0; q<3; q++) {
            finalmat_test_2[w][q] = tmp_rot_mat[w][q];
        }
    }   


    //rotate based on the angle difference two molecules
    testmol.rotate_mol(finalmat_test_2);
    refmol.rotate_mol(finalmat_test_2);

    tmp_rot_mat.clear();
    vec1_test.clear();
    vec2.clear();
    ref_trans_vec.clear();

    //translate back into place via centroid
    ref_trans_vec.push_back(centroid_vec[0][0]);
    ref_trans_vec.push_back(centroid_vec[0][1]);
    ref_trans_vec.push_back(centroid_vec[0][2]);  
   
    target_translation(refmol,ref_trans_vec);
    target_translation(testmol,ref_trans_vec);

}


//gets hrmsd scores in a pair (first: if Du-du-angles are tolerable, second: another PAIR where 
//FIRST: hrmsd, SECOND: num of dissimilar atoms)
//the score is calcualted between refmol and testmol based on the
//the testing three_atom_pairs
std::pair<bool,std::pair <double, int>>
get_hrmsd_and_du_angles(DOCKMol &refmol, DOCKMol & testmol, Iso_Acessory::Scored_Triangle atom_pairs){

    float angle_tolerance = 5.0;

    //instantiate template DOCKMol objects
    DOCKMol tmprefmol;
    DOCKMol tmptestmol;
  
    //copy DOCKMol objects to make templates 
    copy_molecule_shallow(tmprefmol,refmol);
    copy_molecule_shallow(tmptestmol,testmol);

    //initialize variables to get hrmsd 
    std::vector<std::vector<float>> tmp_centroid_coor {};
    Hungarian_RMSD h;
    std::pair<bool,std::pair <double, int>> result;
    result.second.first  = 9999.0;
    result.second.second = 0; 

    //get coordinates of the centroid for both tmprefmol and tmptestmol
    tmp_centroid_coor = get_centroid(tmprefmol,tmptestmol,atom_pairs,true);

    //alignment the molecules based on centroid coordinates
    align_molcentroid(tmprefmol,tmptestmol,tmp_centroid_coor,atom_pairs);
    
    ///get the results for hungarian 
    result.second  = h.calc_Hungarian_RMSD_dissimilar(tmprefmol,tmptestmol);
  
    //get result for du_axi_overlap
    result.first = are_Du_axis_overlapped(tmprefmol,tmptestmol,atom_pairs,angle_tolerance);

    tmprefmol.clear_molecule();
    tmptestmol.clear_molecule();

    return result;
}

bool 
are_Du_axis_overlapped(DOCKMol &refmol, DOCKMol &testmol
                       ,Iso_Acessory::Scored_Triangle atom_pairs
                       ,float cutoff)
{
    DOCKMol tmprefmol;
    DOCKMol tmptestmol;

    std::vector<float> test_trans_vec {};
    std::vector<float> ref_trans_vec  {};


    std::vector<int> refnbrs {};
    std::vector<int> testnbrs {};

    std::vector <float> all_angles {};
    
    for (unsigned int i = 0; i<atom_pairs.three_pairs.size(); i++){
        if (refmol.atom_types[atom_pairs.three_pairs[i].get_ref()] == "Du" &&
            testmol.atom_types[atom_pairs.three_pairs[i].get_test()] == "Du"){

            std::vector <float> vec_ref{};
            std::vector <float> vec_test{};
            
            refnbrs = refmol.get_atom_neighbors(atom_pairs.three_pairs[i].get_ref());
            testnbrs = testmol.get_atom_neighbors(atom_pairs.three_pairs[i].get_test());

            copy_molecule_shallow(tmprefmol, refmol);
            copy_molecule_shallow(tmptestmol, testmol);

            ref_trans_vec.push_back(-tmprefmol.x[atom_pairs.three_pairs[i].get_ref()]);
            ref_trans_vec.push_back(-tmprefmol.y[atom_pairs.three_pairs[i].get_ref()]);
            ref_trans_vec.push_back(-tmprefmol.z[atom_pairs.three_pairs[i].get_ref()]);

            test_trans_vec.push_back(-tmptestmol.x[atom_pairs.three_pairs[i].get_test()]);
            test_trans_vec.push_back(-tmptestmol.y[atom_pairs.three_pairs[i].get_test()]);
            test_trans_vec.push_back(-tmptestmol.z[atom_pairs.three_pairs[i].get_test()]);

            target_translation(tmprefmol,ref_trans_vec);
            target_translation(tmptestmol,test_trans_vec);

            vec_ref.push_back(tmprefmol.x[refnbrs[0]]);
            vec_ref.push_back(tmprefmol.y[refnbrs[0]]);
            vec_ref.push_back(tmprefmol.z[refnbrs[0]]); 

            vec_test.push_back(tmptestmol.x[testnbrs[0]]);
            vec_test.push_back(tmptestmol.y[testnbrs[0]]);
            vec_test.push_back(tmptestmol.z[testnbrs[0]]);
           
             
            if (i==0){all_angles.push_back(get_vector_angle(vec_ref,vec_test));}
            if (i==1){all_angles.push_back(get_vector_angle(vec_ref,vec_test));} 
            if (i==2){all_angles.push_back(get_vector_angle(vec_ref,vec_test));}

            vec_test.clear();
            vec_ref.clear();
            refnbrs.clear();
            testnbrs.clear();
            ref_trans_vec.clear();
            test_trans_vec.clear();
            tmprefmol.clear_molecule();
            tmptestmol.clear_molecule();
        } 
    }

 
    unsigned int true_counter = 0;
    for (unsigned int i=0;i<all_angles.size();i++){
        if (all_angles[i] < cutoff){
            true_counter++;
        }
    }
    if (true_counter == all_angles.size()){return true;}

    return false;
}


//gets hrmsd scores in a pair (first: if Du-du-angles are tolerable, second: another PAIR where 
//FIRST: hrmsd, SECOND: num of dissimilar atoms)
//the score is calcualted between refmol and testmol based on the
//the testing three_atom_pairs
std::pair<bool,std::pair <double, int>>
get_hrmsd_and_du_angles(Fragment &reffrag, Fragment & testfrag
                       ,Iso_Acessory::Scored_Triangle atom_pairs){


    //std::fstream fout_molecules;
    //fout_molecules.open("test_translated.mol2",std::ios_base::app);

    //std::fstream ffout_molecules;
    //ffout_molecules.open("test_translated_true.mol2",std::ios_base::app);

    //instantiate template DOCKMol objects
    DOCKMol tmprefmol;
    DOCKMol tmptestmol;
  
    //copy DOCKMol objects to make templates 
    copy_molecule_shallow(tmprefmol,reffrag.mol);
    copy_molecule_shallow(tmptestmol,testfrag.mol);

    //initialize variables to get hrmsd 
    std::vector<std::vector<float>> tmp_centroid_coor {};
    Hungarian_RMSD h;
    std::pair<bool,std::pair <double, int>> result;
    result.second.first  = 9999.0;
    result.second.second = 0; 

    //get coordinates of the centroid for both tmprefmol and tmptestmol
    tmp_centroid_coor = get_centroid(tmprefmol,tmptestmol,atom_pairs,true);

    //alignment the molecules based on centroid coordinates
    align_molcentroid(tmprefmol,tmptestmol,tmp_centroid_coor,atom_pairs);
    
    ///get the results for hungarian 
    result.second  = h.calc_Hungarian_RMSD_dissimilar(tmprefmol,tmptestmol);
  

     
    //Write_Mol2(tmprefmol,fout_molecules);
    //Write_Mol2(tmptestmol,fout_molecules);

    //get result for du_axi_overlap
    result.first = are_Du_axis_overlapped(tmprefmol,tmptestmol,atom_pairs,5.0);     
   
    //if (result.first == true){Write_Mol2(tmptestmol,ffout_molecules);}

    tmprefmol.clear_molecule();
    tmptestmol.clear_molecule();

    return result;
}







Iso_Acessory::Scored_Triangle get_three_atoms_pairs(Fragment& reffrag
                                                   ,Fragment& testfrag){

    std::cout << testfrag.mol.title << std::endl;
    bool more_than = false; 
    bool du_more_than = false;
    float more_than_para = 1.0;

    //tolerances for linkers
    float more_than_tol_lnk = 3;
    float du_more_than_para_lnk = 1.5;

    float cos_score = 0;

    std::vector<Iso_Acessory::A_Pair>                       heavy_and_du_pairs;
    std::vector<Iso_Acessory::A_Pair>                       du_du_pairs;
    std::vector<Iso_Acessory::A_Pair>                       return_vec;
    std::vector<std::pair<double
                          ,Iso_Acessory::Scored_Triangle>>  sorting_stack;


    Iso_Acessory::Scored_Triangle notalignable;

    //create a pair that can accept a bool and  hrmsd scores
    //from the hungarian function.
    std::pair<bool,std::pair <double, int>> should_align {};
   
    //pairwise calculations to calcualte cos_sim_score 
    //and get the du_du filled vectors.
    for (unsigned int i = 0; i < reffrag.mol.num_atoms; i++){
        for (unsigned int m = 0; m < testfrag.mol.num_atoms; m++){
            cos_score = calc_cos_similarity(reffrag,testfrag, i, m);
                Iso_Acessory::A_Pair one_pair;
                one_pair.set_cos_sim(cos_score);
                one_pair.set_ref(i);
                one_pair.set_test(m);

            if ((reffrag.mol.num_atoms <=5 || testfrag.mol.num_atoms<=5) 
                 && check_H_vs_atom(reffrag,testfrag)){
                heavy_and_du_pairs.push_back(one_pair);

            //Here is where you get all Heavy_Atom-Heavy_atom pairs
            //including du-du
            } else if ((reffrag.mol.atom_types[i] != "H" && testfrag.mol.atom_types[m]!= "H")
                      && (reffrag.mol.atom_types[i] != "Du" && testfrag.mol.atom_types[m]!= "Du")){
                heavy_and_du_pairs.push_back(one_pair);
            }
            //Here is where you get all du-du pairs
            if (reffrag.mol.atom_types[i] == "Du" && testfrag.mol.atom_types[m]== "Du"){
                du_du_pairs.push_back(one_pair);
            }
            one_pair.clear();
        }

    }

    //for (int i = 0;i<heavy_and_du_pairs.size(); i++){
    //    heavy_and_du_pairs[i].print_all(reffrag,testfrag);

    //}
    //for (int i = 0;i<du_du_pairs.size(); i++){
    //    du_du_pairs[i].print_all(reffrag,testfrag);
  
    //}

    //IDONT KNOW IF I HAVE TO USE
    //activate each heavy atom 
    for (unsigned int i =0; i<reffrag.mol.num_atoms;i++){
        reffrag.mol.atom_active_flags[i] = true;
    }
    for (unsigned int i =0; i<testfrag.mol.num_atoms;i++){
        testfrag.mol.atom_active_flags[i] = true;
    }

    //sort the pairs based on cos_similarity
    std::sort(heavy_and_du_pairs.begin()
             ,heavy_and_du_pairs.end()
             ,[](Iso_Acessory::A_Pair a, Iso_Acessory::A_Pair b){
                 return a.get_cos_sim() > b.get_cos_sim();}); 

    std::sort(du_du_pairs.begin()
             ,du_du_pairs.end() 
             ,[](Iso_Acessory::A_Pair a, Iso_Acessory::A_Pair b){ 
                 return a.get_cos_sim() > b.get_cos_sim();}); 

    //if ref and tet fragments are sidechains...
    if (reffrag.get_num_du() == 1 && testfrag.get_num_du() == 1){

        Iso_Acessory::Scored_Triangle temp_triangle;
 
        //There should be just one du-du pair combination
        //for sidechains
        if (du_du_pairs.size()==1){
            temp_triangle.push_back(du_du_pairs[0]);
        }else{
            std::cout << "ERROR: du_du_pairs do not have a size of 1,"
                      << " check your mol2 files...exiting"
                      << std::endl;
            exit(0);
        }

         for (unsigned int j = 0; j < heavy_and_du_pairs.size(); j++){


            if (reffrag.mol.atom_types[heavy_and_du_pairs[j].get_ref()] == "Du" ||
                testfrag.mol.atom_types[heavy_and_du_pairs[j].get_test()] == "Du"){continue;}

            temp_triangle.push_back(heavy_and_du_pairs[j]);

            for (unsigned int z = 0; z < heavy_and_du_pairs.size(); z++) {

                if (reffrag.mol.atom_types[heavy_and_du_pairs[z].get_ref()] == "Du" ||
                    testfrag.mol.atom_types[heavy_and_du_pairs[z].get_test()] == "Du"){continue;}

                temp_triangle.push_back(heavy_and_du_pairs[z]);

                if (temp_triangle.three_pairs.size() == 3 
                    && temp_triangle.check_if_redundant()){
            	    temp_triangle.pop_back(); 
                    continue;
                }

                if (temp_triangle.three_pairs.size() == 2 && temp_triangle.check_if_redundant()){
                    std::cout<< "ERROR: du-du first pair and second non du pair"
                             << "have atom redundancies, when they should not. check your mol2 file. exiting.."<< endl;
                    exit(0);
                }

                should_align.second.first  = 9999.0;
                should_align.second.second = 0; 
                more_than    = false;  

                //get if the distance differences of corresponding triangle sides
                //between refmol and testmol are too much, return true. 
                more_than = get_diff_dist_two_mol(reffrag,testfrag,temp_triangle,more_than_para,false);

                //if true, must skip because you will waste time aligning
                if (more_than == true){temp_triangle.pop_back(); continue;}

                //If false, do the alignment
                if (more_than == false && temp_triangle.three_pairs.size() == 3 && 
                    !temp_triangle.check_if_redundant()){

                    //should these molecules be aligned?
                    should_align  = get_hrmsd_and_du_angles(reffrag,testfrag,temp_triangle);
                    
                    //if du_du angles are too far away, SKIP;      
                    if (!should_align.first){
                        temp_triangle.pop_back();
                        continue;
                    }                        
 
                    //once you get the hrmsd score, push back the three atom pairs
                    sorting_stack.push_back(std::make_pair(should_align.second.first,temp_triangle));

                    //get rid of the top pair on top of LIFO stack
                    temp_triangle.pop_back();
                    continue;
                }
            }

            temp_triangle.pop_back();

        }    

        //sort tmp stack to get the best hrmsd score
        try {

            //deactivate each heavy atom 
            for (unsigned int i =0; i<reffrag.mol.num_atoms;i++){
                reffrag.mol.atom_active_flags[i] = false;
            }
            for (unsigned int i =0; i<testfrag.mol.num_atoms;i++){
                testfrag.mol.atom_active_flags[i] = false;
            }

            //sometimes the tmp_stack will have nothing
            //if so, throw a 0 value, so it can be caught
            if (sorting_stack.empty()) { throw 0;}

            //sort and get the best hrmsd scoring three_atom_pair
            //and return

            std::sort(sorting_stack.begin(),sorting_stack.end()
                      ,[](std::pair<double,Iso_Acessory::Scored_Triangle>  a
                      ,std::pair<double,Iso_Acessory::Scored_Triangle> b)
                      {return a.first < b.first;});
         
            return sorting_stack[0].second;

        } catch(int size_of_stack) {
            std::cout << testfrag.mol.title  
                      << " is not alignable to reference. Size of stack is: " 
                      << size_of_stack 
                      << std::endl;
            return notalignable;
        }

    }else if (reffrag.get_num_du() == 2 && testfrag.get_num_du() == 2){

        //instialize a special data struct just for linkers
        //it will contain all combinations of Du-Du.
        std::pair<Iso_Acessory::Scored_Triangle,Iso_Acessory::Scored_Triangle> pair_of_Du_Dus {}; 

        ////this is where we get all of those combinations of Du-Du
        for (unsigned int i = 0; i < du_du_pairs.size(); i++){
            pair_of_Du_Dus.first.push_back(du_du_pairs[i]); 
            if (pair_of_Du_Dus.first.three_pairs.size() == 2 && pair_of_Du_Dus.first.check_if_redundant()) {
                pair_of_Du_Dus.first.pop_back();
                pair_of_Du_Dus.second.push_back(du_du_pairs[i]);
            }
            else if (pair_of_Du_Dus.first.three_pairs.size() == 3 && pair_of_Du_Dus.first.check_if_redundant()) {
                pair_of_Du_Dus.first.pop_back(); 
                pair_of_Du_Dus.second.push_back(du_du_pairs[i]);
            }
        }

        //start the for loop to assess the du-du distances, alignment,
        //and return the du-dus 
        for (int i = 0; i<2;i++){
         
            //initialize the du_du_more_than variable
            du_more_than      = false;
           
            //template stack for gaining the three pairs 
            Iso_Acessory::Scored_Triangle temp_triangle;

            //since there are only two possible ways to two du-du combos...
            if (i==0){
                temp_triangle = pair_of_Du_Dus.first;
            }else if (i==1){
                temp_triangle = pair_of_Du_Dus.second;
            
            }

            //check if the du-du distances are at least in tolerable distances
            du_more_than = get_diff_dist_two_mol(reffrag,testfrag,temp_triangle,du_more_than_para_lnk,true);

            //if there is a big difference, skip for loop.
            if(du_more_than == true){ continue;}

            //this is where we loop thru to get the THIRD
            //atom-atom pair that excludes Du atoms
            for (unsigned int j = 0; j < heavy_and_du_pairs.size(); j++){

                //assign values to set up a fresh alignment
                should_align.second.first  = 9999.0;
                should_align.second.second = 0; 
                more_than         = false;
           
                //if any of the atomtypes are Du, skip.
                if (reffrag.mol.atom_types[heavy_and_du_pairs[j].get_ref()] == "Du" ||
                    testfrag.mol.atom_types[heavy_and_du_pairs[j].get_test()] == "Du"){ continue;}
               
                //input the next pair into the stack 
                temp_triangle.push_back(heavy_and_du_pairs[j]);

                 //If there are redundancies,  skip 
                if (temp_triangle.three_pairs.size() == 3 && temp_triangle.check_if_redundant()){
                    std::cout<< "ERROR: the du-du pairs and the non_du third pair"
                             << "have atom redundancies, even though it shouldn't be possible"
                             << "please check your mol2 file. exiting..." << endl;
                    exit(0);
                }

                //check if the two Du-Du pair distances are more than a tolerence 
                more_than = get_diff_dist_two_mol(reffrag,testfrag,temp_triangle,more_than_tol_lnk,false);
    
                //if tolerance is broken, popback and skip 
                if (more_than == true){temp_triangle.pop_back();  continue;}

                //if tolerance is within boundaries, align and calc hrmsd.
                if (more_than == false && temp_triangle.three_pairs.size() == 3){

                    should_align  = get_hrmsd_and_du_angles(reffrag,testfrag,temp_triangle);
    
                    //if du_du angles are too far away, SKIP; 
                    if (!should_align.first){ temp_triangle.pop_back();  continue;}
		    sorting_stack.push_back(std::make_pair(should_align.second.first,temp_triangle));
                    temp_triangle.pop_back();
                }
            }

        } 

        //sort tmp stack to get the best hrmsd score
        try {
            //deactivate each heavy atom 
            for (unsigned int i =0; i<reffrag.mol.num_atoms;i++){
                reffrag.mol.atom_active_flags[i] = false;
            }
            for (unsigned int i =0; i<testfrag.mol.num_atoms;i++){
                testfrag.mol.atom_active_flags[i] = false;
            }

            //sometimes the tmp_stack will have nothing
            //if so, throw a 0 value, so it can be caught 
            if (sorting_stack.empty()) { 
                throw 0;
            }

            //sort and get the best hrmsd scoring three_atom_pair
            //and return
            std::sort(sorting_stack.begin(),sorting_stack.end()
                      ,[](std::pair<double,Iso_Acessory::Scored_Triangle>  a
                      ,std::pair<double,Iso_Acessory::Scored_Triangle> b)
                      {return a.first < b.first;});

            return sorting_stack[0].second;
            
        } catch(int error_code) {
            if (error_code == 0) {
                std::cout <<  "ERROR_"<< error_code <<": " 
                          << testfrag.mol.title
                          <<" is not alignable to referece."
                          <<" Size of the linker return stack is: 0." 
                          << std::endl;
                }
            return notalignable; 
        }

    //if scaffolds for 3 du atoms and above
    }else if (reffrag.get_num_du() == testfrag.get_num_du() ){

        //a special data struct: a VECTOR that contains a VECTOR,
        //which contains a PAIR, where first: cos sim and second: an atom pair
        std::vector<Iso_Acessory::Scored_Triangle> vec_three_dus_dus {};
        Iso_Acessory::Scored_Triangle single_du_du {};

        //to collect a vector of three du-du pairs to test out the differences of distances
        for (unsigned int i = 0; i < du_du_pairs.size(); i++){
            single_du_du.push_back(du_du_pairs[i]); 
            for (unsigned int j = 0; j < du_du_pairs.size(); j++){
                single_du_du.push_back(du_du_pairs[j]);

                if (single_du_du.check_if_redundant()) { single_du_du.pop_back(); continue;} 

                for  (unsigned int z = 0; z < du_du_pairs.size(); z++){
                    single_du_du.push_back(du_du_pairs[z]);
                    if (single_du_du.check_if_redundant()) {single_du_du.pop_back(); continue;
                    }else{
                        vec_three_dus_dus.push_back(single_du_du);
                        single_du_du.pop_back();
                    }
                }
                single_du_du.pop_back();
            }
            single_du_du.pop_back();
        }
        single_du_du.clear();
    
        //loop through vec_three_du_du where we can assess if the 
        //distances of corresponding du-du distances between two mol
        //are tolerable 
        std::vector <bool> vec_of_bool {};
        for (unsigned int i = 0;i < vec_three_dus_dus.size(); i++){
            //too check if the du-du distances are at least matching
            du_more_than = get_diff_dist_two_mol(reffrag,testfrag,vec_three_dus_dus[i],1.0,false);
            vec_of_bool.push_back(du_more_than);
        }
 
        //to check of the vec_of_bool and vec_of_dudus have the same size
        try {
            if (vec_of_bool.size() != vec_three_dus_dus.size()){ throw 0;}
            

            unsigned int num_of_false =0;
            unsigned int num_of_true  =0;
            //this line checks of ANY of the elements 
            //are true. If all false, keep going
            for (unsigned int j = 0; j < vec_of_bool.size(); j++){
                if (!vec_of_bool[j]){num_of_false++;}
                if (vec_of_bool[j]){num_of_true++;}
            }
            if (num_of_true == 0) {throw 1;}
            if (num_of_false > 0 && num_of_true > 0) {throw 2;}
            if (num_of_false == 0) {throw 3;} 
 
        }catch (int error_code){

            if (error_code == 0){
                std::cout << "ERROR "<< error_code << ": "
                          <<" Aligning " <<testfrag.mol.title << " and" << reffrag.mol.title
                          <<" are giving an error...REASON:\n" 
                          <<" std::vector<std::pair<float,std::pair<int,int>>> Iso_Align::select_three_atoms_pairs():"  
                          <<" error in matching vector size between 'vec_of_bool' and 'vec_three_dus_dus'."
                          <<" returning empty scaffold vector..." 
                          << std::endl;
                return notalignable;
            }else if (error_code == 1){
                std::cout << "ALIGNABLE " << error_code <<": " 
                          << testfrag.mol.title 
                          <<" All Du-Du distances between ref and testmol"
                          <<" are tolerable and are alignable."
                          << std::endl;
            }else if (error_code == 2){
                std::cout << "WARNING " << error_code <<": " 
                          << testfrag.mol.title << " is could be differet from " 
                          << reffrag.mol.title  << "."
                          <<" Some Du-Du distances could be unalignable."
                          <<" Double check on the alignments."
                          << std::endl; 
            }else if (error_code == 3){
                std::cout << "ERROR " << error_code <<": " 
                          << testfrag.mol.title << " is too differet from " 
                          << reffrag.mol.title  << "."
                          <<" Du-Du distances are too different. Can't be aligned."
                          <<" Returning empty vector."
                          << std::endl;  
                return notalignable;
            }
        }


        //this is where we align and check.
        for (unsigned int i = 0; i < vec_three_dus_dus.size(); i++){

            //if the element is false, align!
            if (!vec_of_bool[i]){

                should_align.first         = false;
                should_align.second.first  = 9999.0;
                should_align.second.second = 0;

                should_align   = get_hrmsd_and_du_angles(reffrag,testfrag,vec_three_dus_dus[i]);
            
                //if du_du angles are too faraway SKIP!
                if (!should_align.first){ 
                    vec_three_dus_dus[i].is_not_alignable();
                    continue;
                }
                vec_three_dus_dus[i].is_alignable();
                sorting_stack.push_back(std::make_pair(should_align.second.first,vec_three_dus_dus[i]));
            }
            vec_three_dus_dus[i].is_not_alignable();
        }


        //try check, if the sorting_stack is 0, throw 0
        try {

            //deactivate each heavy atom 
            for (unsigned int i =0; i<reffrag.mol.num_atoms;i++){
                reffrag.mol.atom_active_flags[i] = false;
            }
            for (unsigned int i =0; i<testfrag.mol.num_atoms;i++){
                testfrag.mol.atom_active_flags[i] = false;
            }

            if (sorting_stack.empty()) { throw 0;}

            //sort tmp stack to get the best hrmsd score
            std::sort(sorting_stack.begin(),sorting_stack.end()
                      ,[](std::pair<double,Iso_Acessory::Scored_Triangle>  a
                      ,std::pair<double,Iso_Acessory::Scored_Triangle> b)
                      {return a.first < b.first;});

            return sorting_stack[0].second;
        }catch (int size_of_stack){
            std::cout << testfrag.mol.title 
                      << " is not alignable to reference. Size of stack is: " 
                      << size_of_stack 
                      << std::endl;
            return notalignable;  
        }

    }else if (reffrag.get_num_du() != testfrag.get_num_du()){
        std::cout << testfrag.mol.title 
                  << " is not alignable to reference due to different num of attachment points. " 
                  << "refmol num of Du: " << reffrag.get_num_du() << " testmol num of Du: " << testfrag.get_num_du()  
                  << " Returning an empty vector."
                  << std::endl; 
        return notalignable;
    } 


    return notalignable;
}

float twoD_length(std::pair<float,float> first, std::pair<float, float> second)
{
    float           len;

    len = sqrt( ( (first.first - second.first) * (first.first - second.first) )   
             +  ( (first.second - second.second) * (first.second - second.second) ) );  
    return len;
}

// calculates the area of a pentagon where the top most plot triangle point is in between two numbered bins
float area_pentagon(std::vector<std::pair <float,float>> half_sq, float bin )
{
    float b1 = 0.0;
    float b2 = 0.0;
    float b3 = 0.0;
    float h1 = 1.0;
    float h2 = 0.0;
    float area1 = 0.0;
    float area2 = 0.0;
    float y1 = 0;  
    float y2 = 0;
   
    if (half_sq[3].second != 0) {
        y1 = half_sq[3].second;
    }else{
        y1 = -bin + half_sq[3].first;
    } 

    if (half_sq[0].second != 0) {
        y2 = half_sq[0].second; 
    }else{
        y2 = (bin-1) + -half_sq[0].first;
    }
    //calcualte the area of the trapezoid
    //b1 is the length of the left of the trapezoid
    b1 =  twoD_length({half_sq[0].first,y2}, {half_sq[0].first, 0});
    
    //b2 is the lenght of the right of the trapzoid
    b2 = twoD_length({half_sq[3].first,y1},{half_sq[3].first, 0}); 
    area1 = (b1 + b2)/2 * h1;

    //std::cout << "area1 " << area1 <<std::endl;
    //calculate the area of the trianlge that sits on top of the trapezoid
    b3 = twoD_length(half_sq[1],{bin,y1});
    h2 = twoD_length(half_sq[1],{bin-1,y2});  
    area2 = (b3 * h2)/2;
    //std::cout << "area2 " << area2 << std::endl;
    //std::cout << b3 << " b3 and h2 " << h2 << std::endl;
    //std::cout << half_sq[3].first << " half_sq[3].first" << std::endl;
    //std::cout << b1 << " b1 and b2 " << b2 << std::endl;
    //std::cout << y1 << " y1 and y2 " << y2 << std::endl;
    //return the total area
    return area1 + area2;
}

//calcualtes the area of the unit triangles needed for the creation of 2D atom plots
std::pair<float,float> area_triangle(std::vector<std::pair <float,float>> half_sq, float bin )
{
    float b1 = 0.0;
    float h1 = 0.0;
    std::pair<float,float> two_area {};
    if (half_sq[0].first < bin){
        h1 = (bin-1) + -half_sq[0].first;   
        b1 = abs(half_sq[0].first - (bin-1)); 
        two_area.first = (b1*h1)/2; 
    } 
    if (half_sq[3].first > bin){
        h1 = -bin + half_sq[3].first;; 
        b1 = abs(half_sq[3].first - (bin));
        two_area.second = (b1*h1)/2;
    } 
    return two_area;
}

//this creates the coordinates of the a single triangle for the
//point. That is calcualted by taking a square placing the 
//diagonal distance and laying on top the of the x-axis and only gettin the positive
//values. it will return four coordinates 
//the structure of single_half_sq figure shown below:
//{left most point <, bottom most point *, top most point ^, right most point >}
//          ^
//        .   .
//      .       .
//    .           .
//  .               .
//<         *         >
std::vector<std::pair <float,float>> get_plot_triangle_coord(float center, float min_bin, float max_bin){
    //initializing the four points to represent them as coordinates of a halfsquare
    std::vector<std::pair <float,float>> half_sq_coordinates {};
    std::pair <float,float> point1 {center-1.0,0.0};
    std::pair <float,float> point2 {center,1.0};
    std::pair <float,float> point3 {center,0.0};
    std::pair <float,float> point4 {center+1.0,0.0};
    //assigning the four points in a vector pair space. 
    half_sq_coordinates = {point1,point2,point3,point4};

    for (unsigned int i=0; i<half_sq_coordinates.size(); i++){
         //Left most triangles are trimmed because some x coordinates could be less than min bin
         if (half_sq_coordinates[i].first < min_bin){
             float b = 0;
             b = half_sq_coordinates[i].second - half_sq_coordinates[i].first;
             half_sq_coordinates[i].first = min_bin;
             half_sq_coordinates[i].second = b;
         }
         //Right most triangles are trimmed beceause some x coordinates could be more than max bin
         if (half_sq_coordinates[i].first > max_bin){
             float b = half_sq_coordinates[i].first;
             float y = -max_bin + abs(b);
             half_sq_coordinates[i].second = y;
             half_sq_coordinates[i].first = max_bin;
         } 
         
         //std::cout<<half_sq_coordinates[i].first << " " << half_sq_coordinates[i].second << std::endl;
    } 
    return half_sq_coordinates;
}

//Where the intial step of iso_align protocol is called.
//1st step -> get radial distances
//2nd step -> get the 2D atom plots for each atom within the internal structure
//3rd step -> calculate the similarity between each atom freq_dist 
//4th step -> align
std::vector <float> get_twoD_atom_plots(std::vector<float> vec_of_centers){

    //these are set min and max bins
    float max_bin = 6.0;
    float min_bin = 0.0;

    //structures to hold triangle areas
    std::pair <float,float> tmp_pair {};

    //The vector that contians the freq distributions 
    std::vector<float> atom_plot(int(max_bin),0.0);

    std::vector<std::pair <float,float>> single_tri_coords {};
   
    for (unsigned int i = 0; i<vec_of_centers.size(); i++){
        single_tri_coords = get_plot_triangle_coord(vec_of_centers[i],min_bin,max_bin);  
        for (int z=0; z<int(max_bin); z++){
            if (single_tri_coords[1].first == z) {
                atom_plot[z-1] = atom_plot[z-1] + 0.5;         
                atom_plot[z] = atom_plot[z] + 0.5;
                continue;
            } 

            if (single_tri_coords[1].first > z-1 && single_tri_coords[1].first < z ) {
                //std::cout << "lol " <<  single_half_sq[1].first << " " << z<< std::endl;
                float  left_area=0;
                float   mid_area=0;
                float right_area=0;
                //std::cout << "midarea start ===========" << std::endl;
                ////check if left side of triangle is hitting a bin
                ////if (single_half_sq[0].second != 0) {
                mid_area =  area_pentagon(single_tri_coords,z);  
                //    std::cout<< "mid_area left" << mid_area << std::endl;
                ////} 
                ////or check if the right side of the triangle is hitting a bin
                ////if (single_half_sq[3].second != 0) {
                ////    mid_area = area_pentagon(single_half_sq,z);  
                ////    std::cout<< "mid_area right" << mid_area << std::endl;
                ////}
                //std::cout << "midarea end   ===========" << std::endl;
                //check if there is a triangle in the bin before or after the current bin number
                tmp_pair   = area_triangle(single_tri_coords,z);
                left_area  = tmp_pair.first;
                right_area = tmp_pair.second;
                     
                atom_plot[z-1] = atom_plot[z-1] + mid_area;
                atom_plot[z-2] = atom_plot[z-2] + left_area;
                atom_plot[z]   = atom_plot[z]   + right_area;
 
                //right_tmp_length =  twoD_length(single_half_sq[3], {single_half_sq[3].first, 0});

                //std::cout << "left LENGTH " << left_tmp_length << " |right LENGTH" << right_tmp_length << std::endl;

            }
        }
    }

    return atom_plot; 
};

float calc_2atoms_length(Fragment  frag, int origin, int target)
{
    float           len;

    len = sqrt( ( (frag.mol.x[origin] - frag.mol.x[target]) * (frag.mol.x[origin] - frag.mol.x[target]) )   
             +  ( (frag.mol.y[origin] - frag.mol.y[target]) * (frag.mol.y[origin] - frag.mol.y[target]) )
             +  ( (frag.mol.z[origin] - frag.mol.z[target]) * (frag.mol.z[origin] - frag.mol.z[target]) ) ); 
    return len;
};

std::vector<float> get_a_radial_dist(Fragment & frag, unsigned int atom_num){

    std::vector<float> list_of_radial_dist {};
    float tmplength = 0;  
    std::vector<float> return_radial_areas_vec {};    

    for (unsigned int i=0; i<frag.mol.num_atoms; i++){
        if (i == atom_num){
            continue; 
        }
           
        tmplength = calc_2atoms_length(frag,atom_num,i);
        //vector filled with radial distances amongst each atom in the internal structure that comes
        //from the "atom_num" atom position of the internal structure. 
        list_of_radial_dist.push_back(tmplength);
    } 

    sort(list_of_radial_dist.begin(), list_of_radial_dist.end());

    //get twoD_atom_plots
    return_radial_areas_vec = get_twoD_atom_plots(list_of_radial_dist);
    

    return return_radial_areas_vec;
}

void
Iso_Align::get_all_radial_dist(Fragment& frag)
{

    std::vector<float> list_of_radial_dist {};
    float tmplength = 0;
    std::vector<float> tmp_radial{};


    for (int origin; origin<frag.mol.num_atoms;origin++){
        for ( int target=0; target<frag.mol.num_atoms; target++){
            if (frag.mol.atom_names[origin] == frag.mol.atom_names[target]){
                continue;
            }

            tmplength = calc_2atoms_length(frag,origin,target);

            //vector filled with radial distances amongst each atom in the internal structure that comes
            //from the "atom_num" atom position of the internal structure. 
            list_of_radial_dist.push_back(tmplength);
        }

        sort(list_of_radial_dist.begin(), list_of_radial_dist.end());

        tmp_radial = get_twoD_atom_plots(list_of_radial_dist);

        frag.set_radial_dist_distri(origin,tmp_radial);

        list_of_radial_dist.clear();
        tmp_radial.clear();
    }
};

//The command where you call the align method
void
Iso_Align::align(Fragment& reffrag, std::vector<Fragment>& vec_testfrags, Iso_Parm parameters){

    int num_frags =                 vec_testfrags.size();
    std::vector<Iso_Acessory::Scored_Triangle>    three_atoms_pairs {};
    std::vector<std::vector<float>> tmp_centroid_coor {}; 

    reffrag.allocate_radial_dist_distri();
    reffrag.calc_num_du();
    reffrag.calc_radial_dist_distri();

    for (unsigned int i = 0; i<vec_testfrags.size(); i++){
        vec_testfrags[i].allocate_radial_dist_distri();
        vec_testfrags[i].calc_num_du();
        vec_testfrags[i].calc_radial_dist_distri();
    }

    for (unsigned int i =0; i<num_frags; i++){
        Iso_Acessory::Scored_Triangle best_triangle;
        best_triangle = get_three_atoms_pairs(reffrag,vec_testfrags[i]);


        if(best_triangle.three_pairs.empty()) {
            vec_testfrags[i].is_not_iso_aligned();
            continue;
        }    
        vec_testfrags[i].is_iso_aligned();
        tmp_centroid_coor = get_centroid(reffrag.mol
                                        ,vec_testfrags[i].mol
                                        ,best_triangle
                                        ,true);
        align_molcentroid(reffrag.mol,vec_testfrags[i].mol,tmp_centroid_coor,best_triangle);

        best_triangle.clear();
    }

    reffrag.clear_radial_dist_distri();

    for (unsigned int i = 0; i<vec_testfrags.size(); i++){
	vec_testfrags[i].clear_radial_dist_distri();
    }

};
