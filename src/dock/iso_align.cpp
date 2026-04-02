#include "fragment.h"
#include "iso_align.h"
#include "dockmol.h"
#include "hungarian.h"
#include "fingerprint.h"
#include "utils.h"

#include <ostream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <cstdio>
#include <sys/stat.h>
#include <functional>
#include <list>
#include <iomanip>
#include <dirent.h>
#include <numeric>  

class A_Pair;
class Scored_Triangle;
class DOCKMol;
class Hungarian_RMSD;
class Domain;
class Score;
class Iso_Parm;

void Iso_Align::frag_sort(std::vector<Fragment> & vec_frag, std::function<bool( Fragment&,  Fragment&)> func){

    std::list<Fragment> list_growing(vec_frag.begin(), vec_frag.end());
    list_growing.sort(func);

    vec_frag.clear();
    for (Fragment tmp_frag : list_growing){
        vec_frag.push_back(tmp_frag); 
    }    

}
bool
Iso_Align::fragment_sort(  Fragment & a, Fragment & b )
{   
    return (a.get_iso_score() > b.get_iso_score());
}
bool
Iso_Align::fragment_sort_reverse(  Fragment & a, Fragment & b )
{   
    return (a.get_iso_score() < b.get_iso_score());
}
Iso_Table::Iso_Tab::Iso_Tab(){
    (*this).best_triangles = {};
}

Iso_Table::Iso_Tab::~Iso_Tab(){

    (*this).clear_table();

}

Iso_Table::Iso_Tab::Iso_Tab(const Iso_Tab& original){

    (*this).clear_table();

    //This is for copying best triangles
    for (const std::vector<Iso_Acessory::Scored_Triangle>  &tail : original.best_triangles){
        std::vector<Iso_Acessory::Scored_Triangle> tmp_vec {};
        (*this).best_triangles.push_back(tmp_vec);
        int i = (*this).best_triangles.size()-1;

        for (const Iso_Acessory::Scored_Triangle &tri : tail){
            (*this).best_triangles[i].push_back(tri);
        }
    }

    //copies the the head and tails
    for (const std::pair<std::string,std::vector<Fragment>> &tmp_pair : original){
        std::string tmp_str = "";
        std::vector<Fragment> tmp_vec {};
        (*this).push_back(std::make_pair(tmp_str,tmp_vec)); 
        int i = (*this).size()-1;

        //assign head and tails
        (*this)[i].first = tmp_pair.first;
        for(const Fragment & frag : tmp_pair.second){
            (*this)[i].second.push_back(frag);
        }

    }

}

void Iso_Table::Iso_Tab::operator=(const Iso_Tab& original){

    (*this).clear_table();

    //This is for copying best triangles
    for (const std::vector<Iso_Acessory::Scored_Triangle>  &tail : original.best_triangles){
        std::vector<Iso_Acessory::Scored_Triangle> tmp_vec {};
        (*this).best_triangles.push_back(tmp_vec);
        int i = (*this).best_triangles.size()-1;

        for (const Iso_Acessory::Scored_Triangle &tri : tail){
            (*this).best_triangles[i].push_back(tri);
        }
    }

    //copies the the head and tails
    for (const std::pair<std::string,std::vector<Fragment>> &tmp_pair : original){
        std::string tmp_str = "";
        std::vector<Fragment> tmp_vec {};
        (*this).push_back(std::make_pair(tmp_str,tmp_vec)); 
        int i = (*this).size()-1;

        //assign head and tails
        (*this)[i].first = tmp_pair.first;
        for(const Fragment & frag : tmp_pair.second){
            (*this)[i].second.push_back(frag);
        }

    }
}

void Iso_Table::Iso_Tab::set(Fragment head, std::vector<Fragment> &tail){

    Fingerprint fing;
    std::vector<Fragment> tmp_vec {};

    std::string head_st = fing.compute_unique_string(head.mol);

    //push back the names of the fragment
    (*this).head_names.push_back(head.mol.energy);

    std::string head_frag_name = head.mol.energy;

    size_t last_index = head_frag_name.find_last_not_of("0123456789");
    int head_frag_num = atoi(head_frag_name.substr(last_index + 1).c_str());
 

    //Where we label frag names.
    for (unsigned int i =0; i < tail.size(); i++){
        Fragment tmp_frag;
        tmp_frag = tail[i];

        std::string tail_frag_name = tmp_frag.mol.energy;
        size_t last_index = tail_frag_name.find_last_not_of("0123456789");
        int tail_frag_num = atoi(tail_frag_name.substr(last_index + 1).c_str());

        //if "iso" is not in the frag_string add it
        std::string tmp_str = tmp_frag.mol.energy;
        if (tmp_str.find("iso_") == string::npos && 
            head_frag_num != tail_frag_num){
            //append the name of frag with an iso
            ostringstream new_name;
            new_name << "iso_" <<tmp_frag.mol.energy << 
                "_" << head_frag_num;
            tmp_frag.mol.energy = new_name.str();
        } else if (tmp_str.find("iso_") == string::npos && 
             head_frag_num == tail_frag_num){
            //append the name of frag with an iso
            ostringstream new_name;
            new_name << "iso_" <<tmp_frag.mol.energy;
            tmp_frag.mol.energy = new_name.str();
        }
        
        //This is here to prepend the subst name with iso_
        for (unsigned int j =0; j< tail[i].mol.num_atoms; j++){

            std::string subst_name = tmp_frag.mol.subst_names[j]; 

            if (subst_name.find("iso_") == string::npos && 
                 head_frag_num != tail_frag_num){
                //append the name of frag with an iso
                ostringstream new_subs_name;
                new_subs_name << "iso_" <<tmp_frag.mol.subst_names[j] <<
                    "_" << head_frag_num;
                tmp_frag.mol.subst_names[j] = new_subs_name.str();

            } else if (subst_name.find("iso_") == string::npos && 
                 head_frag_num == tail_frag_num){
                //append the name of frag with an iso
                ostringstream new_subs_name;
                new_subs_name << "iso_" <<tmp_frag.mol.subst_names[j];;
                tmp_frag.mol.subst_names[j] = new_subs_name.str();
            }else if (subst_name.empty()){
                //append the name of frag with an iso
                tmp_frag.mol.subst_names[j] = "iso";
            }
        }

        //activate all atoms for rotation and translation
        for (unsigned int z =0; z<tmp_frag.mol.num_atoms; z++){ 
            tmp_frag.mol.atom_active_flags[z] == true;
        }

        // Iterate through all bonds and set bond_active_flag to true
        for (int z=0; z<tmp_frag.mol.num_bonds; z++){
                  tmp_frag.mol.bond_active_flags[z] = true;
        }

        tmp_vec.push_back(tmp_frag);
    }

    for(unsigned int i=0; i<tmp_vec.size(); i++){
        if (head_st == fing.compute_unique_string(tmp_vec[i].mol)){
            (*this).head_frags.push_back(tmp_vec[i]);
            break;
        }
        if (i==(tmp_vec.size()-1)) {
            std::cout << "You don't have a head in your snake vector." << endl;
            std::cout << "Are you sure your head fragments are also included in your" << endl;
            std::cout << "`dn_fraglib_iso_[scaffold_linker_sidechain]_file` mol2 files?" << endl;
            std::cout << "Exiting..." << endl;
            exit(0);
        }
    }   
    (*this).push_back(std::make_pair(head_st,tmp_vec));
    tmp_vec.clear();
}

void Iso_Table::Iso_Tab::set_tri(Fragment head, 
std::vector<Iso_Acessory::Scored_Triangle> best_tri){

    if (!(*this).check_if(head)){
        std::cout<< "you must set the heads of the iso table. Exiting..."
                 << std::endl;
        exit(0);
    }

    Fingerprint fing;
    std::string head_st = fing.compute_unique_string(head.mol);

    for (unsigned int i =0; i < (*this).get_size(); i++){
       if (head_st == (*this).get_string(i)){
           (*this).best_triangles.push_back(best_tri);
           break;
       }
    }
}

int Iso_Table::Iso_Tab::get_size(){
    return (*this).size();
}

int Iso_Table::Iso_Tab::get_size_total(){


    int counter = 0;
    for (unsigned int i =0; i < (*this).get_size(); i++){
        counter = counter + (*this)[i].second.size();
    }


    return counter;
}
void Iso_Table::Iso_Tab::print_table(){

}

std::string Iso_Table::Iso_Tab::get_string(int i){

    return (*this)[i].first;
}

bool Iso_Table::Iso_Tab::check_if(Fragment maybe_head){

    Fingerprint fing;
    std::string maybe_head_st = fing.compute_unique_string(maybe_head.mol);
    bool ret_bool = false;
    for (unsigned int i =0; i < (*this).get_size(); i++){
       if (maybe_head_st == (*this).get_string(i)){return true;}
    }
    return ret_bool;

}

std::vector<Fragment> Iso_Table::Iso_Tab::get(Fragment head, int number=100){
    Fingerprint fing;
    std::string head_st = fing.compute_unique_string(head.mol);
    for (unsigned int i =0 ; i < (*this).get_size(); i++){
        if (head_st == (*this).get_string(i)){
            std::vector<Fragment> return_vec;
            for (unsigned int j=0; j<(*this)[i].second.size(); j++){
                if (j<number){
                    Fragment tmp_frag;
                    tmp_frag= (*this)[i].second[j];
                    return_vec.push_back(tmp_frag);
                }
            }
            return return_vec;
        }
    } 
    return {};
}

Fragment Iso_Table::Iso_Tab::get_head(Fragment head){
    Fingerprint fing;
    std::string head_st = fing.compute_unique_string(head.mol);
    Fragment iso_head;

    for (unsigned int i =0 ; i < (*this).get_size(); i++){
        //three_pair handling

        std::string fragment_name = head.mol.energy;
        size_t last_index = fragment_name.find_last_not_of("0123456789");
        int head_frag_num = atoi(fragment_name.substr(last_index + 1).c_str());

        fragment_name = head_frags[i].mol.energy;
        last_index = fragment_name.find_last_not_of("0123456789");
        int tail_frag_num = atoi(fragment_name.substr(last_index + 1).c_str());
        if (tail_frag_num==head_frag_num){
            iso_head = (*this).head_frags[i];
            return iso_head;
        }
    } 
    iso_head.mol.title = "NOT_FOUND";
    return iso_head;
}


//Writes the entire table
void Iso_Table::Iso_Tab::write_table(std::string BASEPATH,int num_top){
    const string DELIMITER    = "########## ";
    const int    STRING_WIDTH = 17 + 19;
    const int    FLOAT_WIDTH  = 20;

    std::ofstream     header;
    std::stringstream base;

    base << BASEPATH << "/isofrags";
    string base_str = base.str();
    mkdir(base_str.c_str(), 0744);

    for (unsigned int i =0 ; i < (*this).get_size(); i++){
        std::ofstream    fout_isofrag;

        std::stringstream name_file; 
        name_file << "isofrag_" << (*this).head_names[i] << ".mol2";

        std::string name_file_str = name_file.str();
     
        fout_isofrag.open(name_file_str,std::ios_base::app);

        fout_isofrag << "#HEAD " << (*this).head_names[i] << endl;


        //Writing head
        std::stringstream three_pairs_str; 
        three_pairs_str << head_frags[i].best_tri.three_pairs[0].get_ref() << "/" <<
            head_frags[i].best_tri.three_pairs[0].get_test()<< "/" <<
            head_frags[i].best_tri.three_pairs[1].get_ref() << "/" <<
            head_frags[i].best_tri.three_pairs[1].get_test()<< "/" <<
            head_frags[i].best_tri.three_pairs[2].get_ref() << "/" <<
            head_frags[i].best_tri.three_pairs[2].get_test()<< std::endl;

        fout_isofrag << "\n" <<

            DELIMITER << std::setw(STRING_WIDTH) <<
            "Head_Name:" << std::setw(FLOAT_WIDTH) <<head_frags[i].mol.energy <<
            endl <<

            DELIMITER << std::setw(STRING_WIDTH) <<
            "Iso_Score:" << std::setw(FLOAT_WIDTH) << head_frags[i].get_iso_score() <<
            endl <<

            DELIMITER << std::setw(STRING_WIDTH) <<
            "Molecular_Weight:" << std::setw(FLOAT_WIDTH) << head_frags[i].mol.mol_wt <<
            endl <<


            DELIMITER << std::setw(STRING_WIDTH) <<
            "FREQ:" << std::setw(FLOAT_WIDTH) << head_frags[i].freq_num <<
            endl <<


            DELIMITER << std::setw(STRING_WIDTH) <<
            "Three_pairs:" << std::setw(21) << three_pairs_str.str() << 
            endl;

        for (unsigned int z =0;z<head_frags[i].mol.num_atoms;z++){
            head_frags[i].mol.atom_active_flags[z] =true;
        }

        Write_Mol2(head_frags[i].mol,fout_isofrag);

        std::stringstream iso_head_name;
        iso_head_name << "iso_" << (*this).head_names[i];


        //Writing tail
        for (unsigned int j =0; j< (*this)[i].second.size(); j++){

            if (iso_head_name.str() == (*this)[i].second[j].mol.energy) {
                continue;
            }

            //when fragment writing counter reaches num_top, break
            if (j >= num_top){
                break;
            }
           
            std::stringstream three_pairs_str; 
            three_pairs_str << (*this)[i].second[j].best_tri.three_pairs[0].get_ref() << "/" <<
                (*this)[i].second[j].best_tri.three_pairs[0].get_test()<< "/" <<
                (*this)[i].second[j].best_tri.three_pairs[1].get_ref() << "/" <<
                (*this)[i].second[j].best_tri.three_pairs[1].get_test()<< "/" <<
                (*this)[i].second[j].best_tri.three_pairs[2].get_ref() << "/" <<
                (*this)[i].second[j].best_tri.three_pairs[2].get_test()<< std::endl;

            fout_isofrag << "\n" <<

                DELIMITER << std::setw(STRING_WIDTH) <<
                "Name:" << std::setw(FLOAT_WIDTH) <<(*this)[i].second[j].mol.energy <<
                endl <<

                DELIMITER << std::setw(STRING_WIDTH) <<
                "Iso_Score:" << std::setw(FLOAT_WIDTH) << (*this)[i].second[j].get_iso_score() <<
                endl <<

                DELIMITER << std::setw(STRING_WIDTH) <<
                "Molecular_Weight:" << std::setw(FLOAT_WIDTH) << (*this)[i].second[j].mol.mol_wt <<
                endl <<


                DELIMITER << std::setw(STRING_WIDTH) <<
                "FREQ:" << std::setw(FLOAT_WIDTH) << (*this)[i].second[j].freq_num <<
                endl <<


                DELIMITER << std::setw(STRING_WIDTH) <<
                "Three_pairs:" << std::setw(21) << three_pairs_str.str() << 
                endl;

            for (unsigned int z =0;z<(*this)[i].second[j].mol.num_atoms;z++){
                (*this)[i].second[j].mol.atom_active_flags[z] =true;
            }

            Write_Mol2((*this)[i].second[j].mol,fout_isofrag);
        }

        fout_isofrag.flush();
        fout_isofrag.close();
        std::stringstream name_file_to;  
        name_file_to << BASEPATH << "/isofrags/"<<name_file_str;
        std::string name_file_to_str = name_file_to.str();
   
        rename(name_file_str.c_str(),name_file_to_str.c_str());
        
    } 
}

std::vector<Iso_Acessory::Scored_Triangle> 
Iso_Table::Iso_Tab::get_tri(Fragment head){
    Fingerprint fing;
    std::string head_st = fing.compute_unique_string(head.mol);
    for (unsigned int i=0; i < (*this).get_size(); i++){
        if (head_st == (*this).get_string(i)){return (*this).best_triangles[i];}
    } 
    return {};
}


void Iso_Table::Iso_Tab::clear_table(){

    (*this).clear();
    (*this).best_triangles.clear();
}

//Iso_Parm methods

void  Iso_Parm::set_dist_du_du_inter(float len){
    this->dist_du_du_inter = len;
}

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
void  Iso_Parm::set_dist_du_du_lnk(float dist){
    this->dist_du_du_lnk=dist;
}
void  Iso_Parm::set_dist_du_du_scf(float dist){
    this->dist_du_du_scf=dist;
}

void  Iso_Parm::set_diff_num_atoms(int diff){
    this->diff_num_atoms = diff;
}
void  Iso_Parm::set_iso_num_top(int num){
    this->iso_num_top = num;
}
void Iso_Parm::set_iso_rank_score_sel(std::string score){
    this->iso_rank_score_sel = score;
}

void  Iso_Parm::set_rank(bool if_rank){
    this->rank = if_rank;
}

void  Iso_Parm::set_write_libraries(bool if_write){
    this->write_libraries = if_write;
}

void Iso_Parm::set_iso_fraglib(bool if_frag,std::string PATH){


    this->iso_fraglib      = if_frag;
    this->iso_fraglib_path = PATH;
    //if (if_frag){
    //    
    //}


}
void Iso_Parm::set_iso_score_sel(std::string score){
    this->iso_score_sel = score;
}

void Iso_Parm::set_iso_rank_reverse(bool reverse){
    this->iso_rank_reverse = reverse;
}
void Iso_Parm::set_iso_write_freq_cutoff(int cutoff){
    this->iso_write_freq_cutoff = cutoff;
};

void Iso_Parm::set_iso_cos_score_cutoff(float cutoff){
    this->iso_cos_score_cutoff = cutoff;
};

float  Iso_Parm::get_dist_du_du_inter(){
    return this->dist_du_du_inter;
}
float Iso_Parm::get_dist_tol_sid(){
    return this->dist_tol_sid;
};
float Iso_Parm::get_dist_tol_lnk(){
    return this->dist_tol_lnk;
}
float Iso_Parm::get_dist_du_du_lnk(){
    return this->dist_du_du_lnk;
}
float Iso_Parm::get_dist_du_du_scf(){
    return this->dist_du_du_scf;
}

bool Iso_Parm::get_rank(){
    return this->rank;
}

bool Iso_Parm::get_write_libraries(){
    return this->write_libraries;
}

int Iso_Parm::get_diff_num_atoms(){
    return this->diff_num_atoms;
}
int Iso_Parm::get_iso_num_top(){
    return this->iso_num_top;
}
float Iso_Parm::get_iso_cos_score_cutoff(){
    return this->iso_cos_score_cutoff;
}


std::pair<bool,
      std::vector<std::string>> Iso_Parm::get_iso_fraglib(){


    std::vector<std::string> tmp {};

    DIR *OPEN_DIR_PATH; 
    struct dirent *dp;
    OPEN_DIR_PATH = opendir((*this).iso_fraglib_path.c_str());
    if (OPEN_DIR_PATH == NULL) {
        printf ("Cannot open directory '%s'\n", (*this).iso_fraglib_path.c_str());
        exit(0);
    }
    while ((dp = readdir(OPEN_DIR_PATH)) != NULL){
            std::string filename(dp->d_name);
            tmp.push_back(filename);
    }
    closedir(OPEN_DIR_PATH);
    
    return std::make_pair(this->iso_fraglib,tmp);
}
std::string Iso_Parm::get_iso_score_sel(){
    return this->iso_score_sel;
}

std::string Iso_Parm::get_iso_rank_score_sel(){
    return this->iso_rank_score_sel;
}

bool Iso_Parm::get_iso_rank_reverse(){
    return this->iso_rank_reverse;
}
int Iso_Parm::get_iso_write_freq_cutoff(){
    return this->iso_write_freq_cutoff;
};

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

Iso_Score::Domain::~Domain(){

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


void Iso_Score::Score::Analytical_method(Fragment & reffrag, Fragment & testfrag)
{
    int         i,
                j;
    float       pi = 3.14,
                d_sqr,
                overlap_temp;


    int reffrag_heavy_atoms =0;

    //create obj for finger print calc
    Fingerprint fing;

    //create obj for hungar rmsd
    Hungarian_RMSD h;

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


    for (i = 0; i < reffrag.mol.num_atoms; i++){
        if (reffrag.mol.atom_active_flags[i]) {
            for (j = 0; j < reffrag.mol.num_atoms; j++){
                if (reffrag.mol.atom_active_flags[j]) {
                    d_sqr = (reffrag.mol.x[i] - reffrag.mol.x[j]) * 
                            (reffrag.mol.x[i] - reffrag.mol.x[j]) + 
                            (reffrag.mol.y[i] - reffrag.mol.y[j]) * 
                            (reffrag.mol.y[i] - reffrag.mol.y[j]) + 
                            (reffrag.mol.z[i] - reffrag.mol.z[j]) * 
                            (reffrag.mol.z[i] - reffrag.mol.z[j]);
                    if (d_sqr < (reffrag.mol.amber_at_radius[i] + reffrag.mol.amber_at_radius[j]) *
                        (reffrag.mol.amber_at_radius[i] + reffrag.mol.amber_at_radius[j])) {
                        if (d_sqr > 0.00001)
                            overlap_temp = pi / 12 * 
                                (reffrag.mol.amber_at_radius[i] + reffrag.mol.amber_at_radius[j] - sqrt(d_sqr)) * 
                                (reffrag.mol.amber_at_radius[i] + reffrag.mol.amber_at_radius[j] - sqrt(d_sqr)) * 
                                (sqrt(d_sqr) + 2 * (reffrag.mol.amber_at_radius[i] + reffrag.mol.amber_at_radius[j]) - 3 / 
                                 sqrt(d_sqr) * (reffrag.mol.amber_at_radius[i] - reffrag.mol.amber_at_radius[j])       * 
                                (reffrag.mol.amber_at_radius[i] - reffrag.mol.amber_at_radius[j]));
                        else
                            overlap_temp = pi / 12 * 2 * 
                                (reffrag.mol.amber_at_radius[i] + reffrag.mol.amber_at_radius[j]) * 
                                (reffrag.mol.amber_at_radius[i] + reffrag.mol.amber_at_radius[j]) *  
                                (reffrag.mol.amber_at_radius[i] + reffrag.mol.amber_at_radius[j]);

                        ref_o.set_overlap(ref_o.get_overlap() + overlap_temp);
                        if (reffrag.mol.atom_types[i] != "H" && reffrag.mol.atom_types[j] != "H")

                            ref_o.set_overlap_hvy(ref_o.get_overlap_hvy() + overlap_temp);

                        if (reffrag.mol.chem_types[i] == "hydrophobic" && reffrag.mol.chem_types[j] == "hydrophobic")
                             ref_o.set_overlap_pho(ref_o.get_overlap_pho() + overlap_temp); 
                        if ((reffrag.mol.chem_types[i] == "donar"    || 
                             reffrag.mol.chem_types[i] == "acceptor" || 
                             reffrag.mol.chem_types[i] == "polar")   && 
                            (reffrag.mol.chem_types[j] == "donar"    ||  
                             reffrag.mol.chem_types[j] == "acceptor" || 
                             reffrag.mol.chem_types[j] == "polar"))

                            ref_o.set_overlap_phi(ref_o.get_overlap_phi() + overlap_temp);
                        if (reffrag.mol.charges[i] < 0 && reffrag.mol.charges[j] < 0)

                            ref_o.set_overlap_neg(ref_o.get_overlap_neg() + overlap_temp);  

                        if (reffrag.mol.charges[i] > 0 && reffrag.mol.charges[j] > 0)

                            ref_o.set_overlap_pos(ref_o.get_overlap_pos() + overlap_temp);
                    }
                }
            }
        }
    }

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
                    if (d_sqr > 0.00001)
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
                        if (d_sqr > 0.00001)
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



    //Hungarian scoring
    std::pair <float, int> h_result = h.calc_Hungarian_RMSD_dissimilar(reffrag.mol,testfrag.mol);


    (*this).hrmsd_score = h_result.first;

    for (int i=0; i<reffrag.mol.num_atoms; i++){
        if (reffrag.mol.atom_types[i].compare("H") != 0 && reffrag.mol.atom_types[i].compare("Du") != 0){ 
            reffrag_heavy_atoms++;
        }   
    }
    (*this).hms_score   = -5 * (reffrag_heavy_atoms -h_result.second) / reffrag_heavy_atoms +
                           1 * (h_result.first);

    (*this).tani_score  = fing.compute_tanimoto(reffrag.mol, testfrag.mol);

    (*this).vos_score = (total_component+heavy_atom_component+
                     negative_component+positive_component+
                     hydrophobic_component+hydrophilic_component)/6;

    (*this).total_component       = total_component;
    (*this).heavy_atom_component  = heavy_atom_component;
    (*this).negative_component    = negative_component;
    (*this).positive_component    = positive_component;
    (*this).hydrophobic_component = hydrophobic_component;
    (*this).hydrophilic_component = hydrophilic_component;
    
}





//Iso_Score class Methods


void Iso_Score::Score::operator=(const Iso_Score::Score& original){


    this->hrmsd_score           = original.hrmsd_score;
    this->hms_score             = original.hms_score;
    this->tani_score            = original.tani_score;
    this->vos_score             = original.vos_score;

    this->score                 = original.score;
    this->total_component       = original.total_component;
    this->heavy_atom_component  = original.heavy_atom_component;
    this->negative_component    = original.negative_component;
    this->positive_component    = original.positive_component;
    this->hydrophobic_component = original.hydrophobic_component;
    this->hydrophilic_component = original.hydrophilic_component;
}
Iso_Score::Score::Score(){
  
    this->hrmsd_score           =9999.0;
    this->hms_score             =9999.0;
    this->tani_score            =0.0; 
    this->vos_score             =0.0; 
    this->score                 =0.0; 
    this->total_component       =0.0; 
    this->heavy_atom_component  =0.0; 
    this->negative_component    =0.0; 
    this->positive_component    =0.0; 
    this->hydrophobic_component =0.0; 
    this->hydrophilic_component =0.0; 

}
Iso_Score::Score::~Score(){
  
    this->hrmsd_score           =9999.0;
    this->hms_score             =9999.0;
    this->tani_score            =0.0; 
    this->vos_score             =0.0; 
    this->score                 =0.0; 
    this->total_component       =0.0; 
    this->heavy_atom_component  =0.0; 
    this->negative_component    =0.0; 
    this->positive_component    =0.0; 
    this->hydrophobic_component =0.0; 
    this->hydrophilic_component =0.0; 

}
Iso_Score::Score::Score( const Iso_Score::Score& original){
  
    this->hrmsd_score           = original.hrmsd_score;
    this->hms_score             = original.hms_score;
    this->tani_score            = original.tani_score;
    this->vos_score             = original.vos_score;
    this->score                 = original.score;
    this->total_component       = original.total_component;
    this->heavy_atom_component  = original.heavy_atom_component;
    this->negative_component    = original.negative_component;
    this->positive_component    = original.positive_component;
    this->hydrophobic_component = original.hydrophobic_component;
    this->hydrophilic_component = original.hydrophilic_component;
}

float Iso_Score::Score::get_score(std::string score_type){

    std::vector<std::string> list_of_scores = {
        "hrmsd_score",
        "hms_score",
        "tanimoto",
        "vos_score",
        "score",
        "total_component",
        "heavy_atom_component",
        "negative_component",
        "positive_component",
        "hydrophobic_component",
        "hydrophilic_component" 
    };
    float hyb_score = 0.0;

    if(score_type == "hrmsd_score"){
        return this->hrmsd_score;
    }
    else if(score_type == "hms_score"){
        return this->hms_score;
    }
    else if(score_type == "tanimoto"){
        return this->tani_score;
    }
    else if(score_type == "vos_score"){
        return this->vos_score;
    }
    else if(score_type == "score"){
        return this->score;
    }
    else if(score_type == "total_component"){
        return this->total_component;
    }
    else if(score_type == "heavy_atom_component"){
        return this->heavy_atom_component;
    }
    else if(score_type == "negative_component"){
        return this->negative_component;
    }
    else if(score_type == "positive_component"){
        return this->positive_component;
    }
    else if(score_type == "hydrophobic_component"){
        return this->hydrophobic_component;
    }
    else if(score_type =="hydrophilic_component"){
        return this->hydrophilic_component;
    } 
    // If not any of these, it is probably a hybrid
    else{

       std::vector <std::string> tokens;
       Tokenizer(score_type,tokens,*";");
       for (int i=0; i<tokens.size(); i++){
           std::vector <std::string> subtokens;
           Tokenizer(tokens[i],subtokens,*",");
 
           if (subtokens.size() != 2) {
               std::cout << "Must either specify weights or score between ';'" << std::endl;
               std::cout << "Or you have more than 2 strings in a score component" << std::endl;
               std::cout << "usage: [weight0],[score_type0];[weight1],[score_type1]" << std::endl;
               std::cout << "example: -1,hms_score;-5,heavy_atom_component" << std::endl;
               exit(1);
           }
 
           if ( !(std::find(list_of_scores.begin(), list_of_scores.end(), subtokens[1]) != list_of_scores.end()) ){
               std::cout << "One of your rank score component '" << subtokens[1] <<"' is not available"  << std::endl;
               exit(1);
           }
           std::string subscore_type = subtokens[1];
           float       subscore = 0;
           float       weight   = std::stof(subtokens[0]);
           if(subscore_type == "hrmsd_score"){
               subscore = this->hrmsd_score;
           }
           else if(subscore_type == "hms_score"){
               subscore = this->hms_score;
           }
           else if(subscore_type == "tanimoto"){
               subscore = this->tani_score;
           }
           else if(subscore_type == "vos_score"){
               subscore = this->vos_score;
           }
           else if(subscore_type == "score"){
               subscore = this->score;
           }
           else if(subscore_type == "total_component"){
               subscore = this->total_component;
           }
           else if(subscore_type == "heavy_atom_component"){
               subscore = this->heavy_atom_component;
           }
           else if(subscore_type == "negative_component"){
               subscore = this->negative_component;
           }
           else if(subscore_type == "positive_component"){
               subscore = this->positive_component;
           }
           else if(subscore_type == "hydrophobic_component"){
               subscore = this->hydrophobic_component;
           }
           else if(subscore_type =="hydrophilic_component"){
               subscore = this->hydrophilic_component;
           } 
           hyb_score = hyb_score + (weight * subscore);
       }
       return hyb_score;

    }

   
    std::cout<<"You didn't select the right scoring function for ranking. exiting..." << std::endl;
    exit(0);
}

void Iso_Score::Score::calc_score(Fragment& reffrag,Fragment& testfrag){
 
    (*this).Analytical_method(reffrag,testfrag);
}


//Iso_Score class Methods END


//A_Pair class methods

Iso_Acessory::A_Pair::A_Pair(){
  
    this->title = "";
    this->cos_sim = 0.0;
    (*this).first = 0;
    (*this).second = 0;

}

Iso_Acessory::A_Pair::A_Pair( const A_Pair & original){

    (*this).title = original.title; 
    (*this).cos_sim = original.cos_sim; 
    (*this).first = original.first;
    (*this).second = original.second;


}

void Iso_Acessory::A_Pair::operator=( const A_Pair & original){

    (*this).title = original.title; 
    (*this).cos_sim = original.cos_sim; 
    (*this).first = original.first;
    (*this).second = original.second;

}

Iso_Acessory::A_Pair::~A_Pair(){
  
    this->title = "";
    this->cos_sim = 0.0;
    (*this).first = 0;
    (*this).second = 0;

}
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
    this->hrmsd = 9999;
    this->volume_overlap = 9998;
    (*this).three_pairs = {};
}

Iso_Acessory::Scored_Triangle::Scored_Triangle( const Scored_Triangle & original){

    (*this).clear_triangle();

    (*this).alignable = original.alignable; 
    (*this).hrmsd = original.hrmsd; 
    (*this).volume_overlap = original.volume_overlap; 
    
    for (const A_Pair & one_pair : original.three_pairs){
        (*this).three_pairs.push_back(one_pair);
    }
    if ((*this).three_pairs.size() > 3) {
        std::cout << "During assigning Scored_Triangle, the threepairs have more than 3 pairs. \n exiting..." << std::endl;
        exit(0);
    }

}

void Iso_Acessory::Scored_Triangle::operator=( const Scored_Triangle & original){

    (*this).clear_triangle();

    (*this).alignable = original.alignable; 
    (*this).hrmsd = original.hrmsd; 
    (*this).volume_overlap = original.volume_overlap; 

    for (const A_Pair & one_pair : original.three_pairs){
        (*this).three_pairs.push_back(one_pair);
    }
    if ((*this).three_pairs.size() > 3) {
        std::cout << "During assigning Scored_Triangle, the threepairs have more than 3 pairs. \n exiting..." << std::endl;
        exit(0);
    }

}

Iso_Acessory::Scored_Triangle::~Scored_Triangle(){
  
    this->alignable = false;
    this->hrmsd = 9999;
    this->volume_overlap = 9998;
    (*this).three_pairs.clear();
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
        float cos = (*this).three_pairs[i].get_cos_sim();
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



void Iso_Acessory::Scored_Triangle::clear_triangle(){
    (*this).three_pairs.clear();
    this->hrmsd = 9999;
    this->volume_overlap = 9998;
    this->alignable = false;

}

void Iso_Acessory::Scored_Triangle::is_not_alignable(){

    this->alignable = false;

}

void Iso_Acessory::Scored_Triangle::is_alignable(){

    this->alignable = true;

}

void Iso_Acessory::Scored_Triangle::set_hrmsd(float score){
    this->hrmsd = score;
}

float Iso_Acessory::Scored_Triangle::get_hrmsd(){

    if (this->hrmsd == 9999){
        std::cout << "Iso_Acessory::Scored_Triangle::get_hrmsd: you are trying to access an unscored iso frag" << std::endl;
        exit(0);
    }
    return this->hrmsd;
}

void Iso_Acessory::Scored_Triangle::set_vos(float score){
    this->volume_overlap = score;
}

float Iso_Acessory::Scored_Triangle::get_vos(){

    if (this->volume_overlap== 9998){
        std::cout << "Iso_Acessory::Scored_Triangle::get_vos: you are trying to access an unscored iso frag" << std::endl;
        exit(0);
    }
    return this->volume_overlap;
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
//checks num of H relative to heavy atoms
bool Iso_Align::check_H_vs_atom (Fragment & reffrag, Fragment & testfrag){
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

//calculates dot_product
float Iso_Align::dot_product(std::vector<float> vec1,std::vector<float>vec2) {
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
std::vector<float> Iso_Align::cross_prod(const std::vector<float> & v1, const std::vector<float> & v2) 
{
    std::vector <float> vec(3,0.0);

    vec[0] = v1[1] * v2[2] - v1[2] * v2[1];
    vec[1] = -v1[0] * v2[2] + v1[2] * v2[0];
    vec[2] = v1[0] * v2[1] - v1[1] * v2[0];

    return vec;
}

//calcualtes the length of a vector
float Iso_Align::length_vec(std::vector<float> len_vec)
{
    float           len;

    len = sqrt( ( len_vec[0] * len_vec[0] )   
             +  ( len_vec[1] * len_vec[1] )
             +  ( len_vec[2] * len_vec[2] ));  
    return len;
}

//substracts a vector
std::vector<float> Iso_Align::subtract_vec(const std::vector<float> &v1, const std::vector<float> &v2)
{
    std::vector<float> diff {};

    diff.push_back(v1[0] - v2[0]);
    diff.push_back(v1[1] - v2[1]);
    diff.push_back(v1[2] - v2[2]);

    return diff;
}

//you have to pass through by reference or the molecules will be distorted
void Iso_Align::normalize_vec(std::vector <float> &vec_to_norm)
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
float Iso_Align::magnitude(std::vector<float> mag_vec){
    float mag_result = 0.0;
    for (unsigned int i =0; i<mag_vec.size(); i++){
        mag_result += (mag_vec[i] * mag_vec[i]);
    }

    return mag_result;
}

float Iso_Align::calc_2atoms_length(Fragment  frag, int origin, int target)
{
    float           len;

    len = sqrt( ( (frag.mol.x[origin] - frag.mol.x[target]) * (frag.mol.x[origin] - frag.mol.x[target]) )   
             +  ( (frag.mol.y[origin] - frag.mol.y[target]) * (frag.mol.y[origin] - frag.mol.y[target]) )
             +  ( (frag.mol.z[origin] - frag.mol.z[target]) * (frag.mol.z[origin] - frag.mol.z[target]) ) ); 
    return len;
};
//calculates the vector angle
float Iso_Align::get_vector_angle(const std::vector<float> & v1, const std::vector<float> & v2)
{
    float           mag,
                    prod;
    float           result;
    float           pi = 3.1415926536;

    mag = length_vec(v1) * length_vec(v2);

    if (mag == 0) {
        std::cout << "A magnitude value is 0. It will return calc a -nan value. exiting..." <<
            std::endl;
        exit(0);
    }

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
//This is a different torsion angle calculator
//than the one we have in the code base
float Iso_Align::get_torsion_angle(const std::vector<float> v1
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


//calculates cosine similarity
float Iso_Align::calc_cos_similarity( Fragment &reffrag, Fragment & testfrag,int atomref, int atomtest){

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

float 
Iso_Align::threeD_length(std::vector<float> first, std::vector<float> second)
{
    float           len;

    len = sqrt( ( (first[0] - second[0]) * (first[0] - second[0] ) )   
             +  ( (first[1] - second[1]) * (first[1] - second[1] ) ) 
             +  ( (first[2] - second[2]) * (first[2] - second[2] ) ));  
    return len;

}

float
Iso_Align::get_dist_diff_atat(DOCKMol & refmol, DOCKMol & testmol, int ref_atomnum, int test_atomnum)
{

    float return_dist = 0;
    std::vector<float> ref_vec{};
    std::vector<float> test_vec{};

    ref_vec.push_back(refmol.x[ref_atomnum]);
    ref_vec.push_back(refmol.y[ref_atomnum]);
    ref_vec.push_back(refmol.z[ref_atomnum]);
    
    test_vec.push_back(testmol.x[test_atomnum]);
    test_vec.push_back(testmol.y[test_atomnum]);
    test_vec.push_back(testmol.z[test_atomnum]);

    return_dist = threeD_length(ref_vec,test_vec);

    return return_dist;
}
//gets the difference between two atom positions in length
bool 
Iso_Align::get_diff_dist_two_mol (Fragment & reffrag, 
		       Fragment & testfrag,
		       Iso_Acessory::Scored_Triangle triangle,
		       float cutoff,
		        bool sidechains_or_scaffold)
{

    bool diff         = false;
    float ref_0_1_dist  = 0.0; 
    float ref_0_2_dist  = 0.0;
    float ref_1_2_dist  = 0.0;

    float test_0_1_dist = 0.0;
    float test_0_2_dist = 0.0;
    float test_1_2_dist = 0.0;



    //If linker
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
std::vector<std::vector<float>> Iso_Align::get_centroid(DOCKMol & refmol
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
void Iso_Align::target_translation (DOCKMol & trans_mol,std::vector<float> trans_vec) {

    for (unsigned int i = 0; i<trans_mol.num_atoms; i++){
        trans_mol.x[i]+=trans_vec[0];
        trans_mol.y[i]+=trans_vec[1];
        trans_mol.z[i]+=trans_vec[2];
    }
}

//calculates and returns rotation matrix
//that swivels the moving vec to the ref vec
std::vector<std::vector<double>> Iso_Align::get_rotation_mat( std::vector<float>& vec1_test, std::vector<float>& vec2) {
    std::vector<std::vector<double>> final_mat_vec;

    double dot_test        = 0.0; 
    double vec1_magsq_test = 0.0; 
    double vec2_magsq      = 0.0; 

    double cos_theta_test  = 0.0; 
    double sin_theta_test  = 0.0; 
   
    double finalmat_test[3][3];

    dot_test        =dot_product(vec1_test,vec2);

    vec1_magsq_test =magnitude(vec1_test);
    vec2_magsq      =magnitude(vec2); 

    cos_theta_test = dot_test/ (sqrt(vec1_magsq_test* vec2_magsq));

    if (cos_theta_test > 1.0){
        cos_theta_test = 1.0;
        sin_theta_test = 0.0;
    }
    else if (cos_theta_test < -1.0){
        cos_theta_test = -1.0;
        sin_theta_test = 0.0;
    }
    else{
        sin_theta_test = sqrt (1 - (cos_theta_test * cos_theta_test));
    }

    if (cos_theta_test == -1){
        finalmat_test[0][0] =  -1; finalmat_test[0][1] = 0; finalmat_test[0][2] = 0;
        finalmat_test[1][0] =  0; finalmat_test[1][1] = -1; finalmat_test[1][2] = 0;
        finalmat_test[2][0] =  0; finalmat_test[2][1] = 0; finalmat_test[2][2] = -1;
    }else if (cos_theta_test == 1){
        finalmat_test[0][0] =  1; finalmat_test[0][1] = 0; finalmat_test[0][2] = 0;
        finalmat_test[1][0] =  0; finalmat_test[1][1] = 1; finalmat_test[1][2] = 0;
        finalmat_test[2][0] =  0; finalmat_test[2][1] = 0; finalmat_test[2][2] = 1;


    }else if(cos_theta_test != 1){
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

          // make inverse matrix of coordRot matrix - since coorRot is an orthonal matrix 
          // the inverse if its transpose, coorRot^T
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
    }

    //convert double 3 by 3 into vector mat 
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

void Iso_Align::align_molcentroid (DOCKMol& refmol, DOCKMol& testmol,std::vector<std::vector<float>> centroid_vec,Iso_Acessory::Scored_Triangle three_atom_pairs){

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

//fn that asks...are the Du vectors pointing at 
//neighboring atom overlapping?
bool 
Iso_Align::are_Du_axis_overlapped(DOCKMol &refmol, DOCKMol &testmol
                       ,Iso_Acessory::Scored_Triangle atom_pairs
                       ,float angle_cutoff,float dist_cutoff)
{
    DOCKMol tmprefmol;
    DOCKMol tmptestmol;

    std::vector<float> test_trans_vec {};
    std::vector<float> ref_trans_vec  {};

    std::vector<int> refnbrs {};
    std::vector<int> testnbrs {};
    std::vector <float> all_angles {};
    std::vector <float> all_distances {};
    

    //Get distanace difference betwen two fragments bases on Du-Du pair
    for (unsigned int i = 0; i<atom_pairs.three_pairs.size(); i++){
        if (refmol.atom_types[atom_pairs.three_pairs[i].get_ref()] == "Du" &&
            testmol.atom_types[atom_pairs.three_pairs[i].get_test()] == "Du"){
            
            int refmol_atomnum = atom_pairs.three_pairs[i].get_ref();
            int testmol_atomnum = atom_pairs.three_pairs[i].get_test();

            if (i==0){all_distances.push_back(get_dist_diff_atat(refmol,testmol,refmol_atomnum,testmol_atomnum));}
            if (i==1){all_distances.push_back(get_dist_diff_atat(refmol,testmol,refmol_atomnum,testmol_atomnum));} 
            if (i==2){all_distances.push_back(get_dist_diff_atat(refmol,testmol,refmol_atomnum,testmol_atomnum));}
        }
    }


    //Then translate the du-du atoms to origin, then calcualte the angle for the neighboring atoms.
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
    for (unsigned int i=0;i<all_distances.size();i++){
        if (all_distances[i] < dist_cutoff){
            true_counter++;
        }
    }
 
    for (unsigned int i=0;i<all_angles.size();i++){
        if (all_angles[i] < angle_cutoff){
            true_counter++;
        }
    }
    if (true_counter == (all_angles.size() + all_distances.size()) ){return true;}

    return false;
}


//gets hrmsd scores in a pair (first: if Du-du-angles are tolerable, second: another PAIR where 
//FIRST: hrmsd, SECOND: num of dissimilar atoms)
//the score is calcualted between refmol and testmol based on the
//the testing three_atom_pairs
std::pair<bool,std::pair <Iso_Score::Score, int>>
Iso_Align::get_aliscore_and_du_angles(Fragment reffrag, Fragment testfrag
                       ,Iso_Acessory::Scored_Triangle atom_pairs,float ang_tol, float dist_tol){
    float angle_tolerance = ang_tol;
    float du_du_dist_tolerance = dist_tol;


    ////instantiate template DOCKMol objects
    //DOCKMol tmprefmol;
    //DOCKMol tmptestmol;
  
    Iso_Score::Score tmp_score;

    //initialize variables to get hrmsd 
    std::vector<std::vector<float>> tmp_centroid_coor {};
    Hungarian_RMSD h;
    std::pair<bool,std::pair <Iso_Score::Score, int>> result;
    result.second.second = 0; 

    //get coordinates of the centroid for both tmprefmol and tmptestmol
    tmp_centroid_coor = get_centroid(reffrag.mol,testfrag.mol,atom_pairs,true);

    //alignment the molecules based on centroid coordinates
    align_molcentroid(reffrag.mol,testfrag.mol,tmp_centroid_coor,atom_pairs);
    
    ///get the results for hungarian 
    //result.second  = h.calc_Hunrarian_RMSD_dissimilar(reffrag.mol,testfrag.mol);
    
    tmp_score.calc_score(reffrag,testfrag);

    result.second.first = tmp_score;

    //get result for du_axi_overlap
    result.first = are_Du_axis_overlapped(reffrag.mol,testfrag.mol,atom_pairs
                                          ,angle_tolerance,du_du_dist_tolerance); 

    return result;
}


Iso_Acessory::Scored_Triangle Iso_Align::get_three_atoms_pairs(Fragment& reffrag
                                                   ,Fragment& testfrag, Iso_Parm parameters){


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
    std::vector<std::pair<float
                          ,Iso_Acessory::Scored_Triangle>>  sorting_stack;


    Iso_Acessory::Scored_Triangle notalignable;

    //create a pair that can accept a bool and  hrmsd scores
    //from the hungarian function.
    std::pair<bool,std::pair <Iso_Score::Score, int>> should_align {};
   
    //pairwise calculations to calcualte cos_sim_score 
    //and get the du_du filled vectors.
    for (unsigned int i = 0; i < reffrag.mol.num_atoms; i++){
        for (unsigned int m = 0; m < testfrag.mol.num_atoms; m++){
            cos_score = calc_cos_similarity(reffrag,testfrag, i, m);
            Iso_Acessory::A_Pair one_pair;
            one_pair.set_cos_sim(cos_score);
            one_pair.set_ref(i);
            one_pair.set_test(m);

            //Where you set your cutoff of cos_score cutoff
            if (cos_score < parameters.get_iso_cos_score_cutoff()){

                //std::cout << "REJ:"<< cos_score << " " << parameters.get_iso_cos_score_cutoff() << std::endl;
                one_pair.clear();
                continue;
            } 

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

    //std::cout << "heavy_and_du_pairs.size() " << heavy_and_du_pairs.size() << std::endl; 
    //If there are no three pairs to align with, return not alignable
    if (heavy_and_du_pairs.size() == 0) { return notalignable; }
 
    //If there are no three pairs for du_dus, return not alignable
    if (du_du_pairs.size() == 0) { return notalignable; }

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

    //if the differences of num_atoms are too much...
    if (abs(reffrag.mol.num_atoms - testfrag.mol.num_atoms) > parameters.get_diff_num_atoms()){
	return notalignable;
    }

    //if the freq_num cutoff is not met...
    if (testfrag.freq_num < parameters.get_iso_write_freq_cutoff()){
        return notalignable;
    }

    //std::cout <<reffrag.get_num_du() << " num_du  " <<testfrag.get_num_du() << std::endl;
    //std::cout << " du_du_size " <<du_du_pairs.size() << std::endl;
    //for (int z = 0;z<du_du_pairs.size(); z++){
    //    du_du_pairs[z].print_all(reffrag,testfrag);
    //}
    //if ref and test fragments are SIDECHAINS...
    if (reffrag.get_num_du() == 1 && testfrag.get_num_du() == 1 && du_du_pairs.size()==1){

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

	//get info on neighboring atoms
	std::vector <int> ref_neighbhors_du = 
	    reffrag.mol.get_atom_neighbors(temp_triangle.three_pairs[0].get_ref());
	std::vector <int> test_neighbhors_du = 
	    testfrag.mol.get_atom_neighbors(temp_triangle.three_pairs[0].get_test());
  
        // Filter out if the bond types for the sidechains between du and heavy atom  
        // are the same. if not return nothign 
        int bond_id_ref =  reffrag.mol.get_bond(temp_triangle.three_pairs[0].get_ref(),
                                                ref_neighbhors_du[0]);
        int bond_id_test=  testfrag.mol.get_bond(temp_triangle.three_pairs[0].get_test(),
                                                 test_neighbhors_du[0]);

        std::string bond_type_ref = reffrag.mol.bond_types[bond_id_ref];
        std::string bond_type_test= testfrag.mol.bond_types[bond_id_test];
        // check if the bond_types are same
        if (bond_type_ref != bond_type_test) {return notalignable;}
 
	bool if_nbr_exit = false;

        for (unsigned int j = 0; j < heavy_and_du_pairs.size(); j++){
	    //check if atoms are Du, if it is, skip
            if (reffrag.mol.atom_types[heavy_and_du_pairs[j].get_ref()] == "Du" ||
                testfrag.mol.atom_types[heavy_and_du_pairs[j].get_test()] == "Du"){continue;}
	    
	    //get neighboring atoms of the Du atom, and set as second atom pair
	    if (reffrag.mol.atom_names[ref_neighbhors_du[0]] == reffrag.mol.atom_names[heavy_and_du_pairs[j].get_ref()]   &&
		testfrag.mol.atom_names[test_neighbhors_du[0]] == testfrag.mol.atom_names[heavy_and_du_pairs[j].get_test()]
		){
		if_nbr_exit = true;
		temp_triangle.push_back(heavy_and_du_pairs[j]);
	    }
	}

	if (if_nbr_exit == false){
	    //std::cout << "ERROR: during calculations Du atoms show now neighboring atoms. Exiting..." 
	    //    << std::endl;
	    return notalignable;
	}

	//Where you do the real alignment
	for (unsigned int j = 0; j < heavy_and_du_pairs.size(); j++){
            if (reffrag.mol.atom_types[heavy_and_du_pairs[j].get_ref()] == "Du" ||
                testfrag.mol.atom_types[heavy_and_du_pairs[j].get_test()] == "Du"){continue;}

	    temp_triangle.push_back(heavy_and_du_pairs[j]);

            if (temp_triangle.three_pairs.size() == 3 
                && temp_triangle.check_if_redundant()){
              temp_triangle.pop_back(); 
                continue;
            }

            if (temp_triangle.three_pairs.size() == 2 && temp_triangle.check_if_redundant()){
                //std::cout<< "ERROR: du-du first pair and second non du pair"
                //         << "have atom redundancies, when they should not. check your mol2 file. exiting.."<< endl;
                return notalignable;
            }

            should_align.second.second = 0; 
            more_than    = false;  
     
            //get if the distance differences of corresponding triangle sides
            //between refmol and testmol are too much, return true. 
            more_than = get_diff_dist_two_mol(reffrag,testfrag,temp_triangle
                                              ,parameters.get_dist_tol_sid(),false);
     
            //if true, must skip because you will waste time aligning
            if (more_than == true){temp_triangle.pop_back(); continue;}
     
            //If false, do the alignment
            if (more_than == false && temp_triangle.three_pairs.size() == 3 && 
                !temp_triangle.check_if_redundant()){
     
                //should these molecules be aligned?
                should_align  = get_aliscore_and_du_angles(reffrag,testfrag
                          ,temp_triangle,parameters.get_bond_angle_tol_sid()
                          ,parameters.get_dist_du_du_inter());
                
                //if du_du angles are too far away, SKIP;      
                if (!should_align.first){
                    temp_triangle.pop_back();
                    continue;
                }                        

     
                //once you get the hrmsd score, push back the three atom pairs
                float iso_score = should_align.second.first.get_score(parameters.get_iso_score_sel());
                sorting_stack.push_back(std::make_pair(iso_score,temp_triangle));
     
                //get rid of the top pair on top of LIFO stack
                temp_triangle.pop_back();
                continue;
            }
	}


        //sort tmp stack to get the best hrmsd score
        try {

            //deactivate each heavy atom 
            //for (unsigned int i =0; i<reffrag.mol.num_atoms;i++){
            //    reffrag.mol.atom_active_flags[i] = false;
            //}
            //for (unsigned int i =0; i<testfrag.mol.num_atoms;i++){
            //    testfrag.mol.atom_active_flags[i] = false;
            //}

            //sometimes the tmp_stack will have nothing
            //if so, throw a 0 value, so it can be caught
            if (sorting_stack.empty()) { throw 0;}

            //sort and get the best hrmsd scoring three_atom_pair
            //and return

            if (parameters.get_iso_score_sel() == "hrmsd_score" || parameters.get_iso_score_sel() == "hms_score"){
                std::sort(sorting_stack.begin(),sorting_stack.end()
                          ,[](std::pair<float,Iso_Acessory::Scored_Triangle>  a
                          ,std::pair<float,Iso_Acessory::Scored_Triangle> b)
                          {return a.first < b.first;});
            } else {

                std::sort(sorting_stack.begin(),sorting_stack.end()
                          ,[](std::pair<float,Iso_Acessory::Scored_Triangle>  a
                          ,std::pair<float,Iso_Acessory::Scored_Triangle> b)
                          {return a.first > b.first;});
            }

	    if (parameters.get_rank()){
		//set the hrmsd score for later ranking
                float tmpscore =  sorting_stack[0].first;
		sorting_stack[0].second.set_vos(tmpscore);
	    }
	
            return sorting_stack[0].second;

        } catch(int size_of_stack) {
            //std::cout << testfrag.mol.title  
            //          << " is not alignable to reference. Size of stack is: " 
            //          << size_of_stack 
            //          << std::endl;
            return notalignable;
        }

    //IF LINKER
    }else if ( reffrag.get_num_du() == 2 && 
               testfrag.get_num_du() == 2 && 
               du_du_pairs.size() == 4 ){

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
  
	//get info on neighboring atoms
	std::vector <int> ref_neighbhors_du_1 = 
	    reffrag.mol.get_atom_neighbors(pair_of_Du_Dus.first.three_pairs[0].get_ref());
	std::vector <int> test_neighbhors_du_1 = 
	    testfrag.mol.get_atom_neighbors(pair_of_Du_Dus.first.three_pairs[0].get_test());


	std::vector <int> ref_neighbhors_du_2 = 
	    reffrag.mol.get_atom_neighbors(pair_of_Du_Dus.first.three_pairs[1].get_ref());
	std::vector <int> test_neighbhors_du_2 = 
	    testfrag.mol.get_atom_neighbors(pair_of_Du_Dus.first.three_pairs[1].get_test());


        // Filter out if the bond types for the sidechains between du and heavy atom  
        // are the same. if not return nothign 

        int bond_id_ref_first_1 = reffrag.mol.get_bond(pair_of_Du_Dus.first.three_pairs[0].get_ref(),
                                                       ref_neighbhors_du_1[0]);
        int bond_id_test_first_1 = testfrag.mol.get_bond(pair_of_Du_Dus.first.three_pairs[0].get_test(),
                                                         test_neighbhors_du_1[0]);

        int bond_id_ref_first_2 = reffrag.mol.get_bond(pair_of_Du_Dus.first.three_pairs[1].get_ref(),
                                                       ref_neighbhors_du_2[0]);
        int bond_id_test_first_2 = testfrag.mol.get_bond(pair_of_Du_Dus.first.three_pairs[1].get_test(),
                                                         test_neighbhors_du_2[0]);


        std::vector<std::string> ref_bond_types {reffrag.mol.bond_types[bond_id_ref_first_1],
                                                 reffrag.mol.bond_types[bond_id_ref_first_2]};


        std::vector<std::string> test_bond_types {testfrag.mol.bond_types[bond_id_test_first_1],
                                                  testfrag.mol.bond_types[bond_id_test_first_2]};


        // if there is absolutely no match in bond types for
        // at least one bond, you must return unalignable.
        for (unsigned int i =0; i<ref_bond_types.size(); i++){
            int check = 0;
            for (unsigned int j =0; j<test_bond_types.size(); j++){
                if (ref_bond_types[i] == test_bond_types[j]){
                    check++;
                }
            }
            if (check == 0){
               return notalignable;
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
            du_more_than = get_diff_dist_two_mol(reffrag,testfrag,temp_triangle,parameters.get_dist_du_du_lnk(),true);

            //if there is a big difference, skip for loop.
            if(du_more_than == true){ continue;}

            //this is where we loop thru to get the THIRD
            //atom-atom pair that excludes Du atoms
            for (unsigned int j = 0; j < heavy_and_du_pairs.size(); j++){

                //assign values to set up a fresh alignment
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
                more_than = get_diff_dist_two_mol(reffrag,testfrag,temp_triangle,parameters.get_dist_tol_lnk(),false);
    
                //if tolerance is broken, popback and skip 
                if (more_than == true){temp_triangle.pop_back();  continue;}

                //if tolerance is within boundaries, align and calc hrmsd.
                if (more_than == false && temp_triangle.three_pairs.size() == 3){

                    should_align  = get_aliscore_and_du_angles(reffrag,testfrag
                                    ,temp_triangle,parameters.get_bond_angle_tol_lnk()
                                    ,parameters.get_dist_du_du_inter());
    
                    //if du_du angles are too far away, SKIP; 
                    if (!should_align.first){ temp_triangle.pop_back();  continue;}
                    float iso_score = should_align.second.first.get_score(parameters.get_iso_score_sel());
                    sorting_stack.push_back(std::make_pair(iso_score,temp_triangle));
                    temp_triangle.pop_back();
                }
            }

        } 

        //sort tmp stack to get the best hrmsd score
        try {
            //deactivate each heavy atom 
            //for (unsigned int i =0; i<reffrag.mol.num_atoms;i++){
            //    reffrag.mol.atom_active_flags[i] = false;
            //}
            //for (unsigned int i =0; i<testfrag.mol.num_atoms;i++){
            //    testfrag.mol.atom_active_flags[i] = false;
            //}

            //sometimes the tmp_stack will have nothing
            //if so, throw a 0 value, so it can be caught 
            if (sorting_stack.empty()) { 
                throw 0;
            }

            //sort and get the best hrmsd scoring three_atom_pair
            //and return
            if (parameters.get_iso_score_sel() == "hrmsd_score" || 
                    parameters.get_iso_score_sel() == "hms_score"){
                std::sort(sorting_stack.begin(),sorting_stack.end()
                          ,[](std::pair<float,Iso_Acessory::Scored_Triangle>  a
                          ,std::pair<float,Iso_Acessory::Scored_Triangle> b)
                          {return a.first < b.first;});
            } else {

                std::sort(sorting_stack.begin(),sorting_stack.end()
                          ,[](std::pair<float,Iso_Acessory::Scored_Triangle>  a
                          ,std::pair<float,Iso_Acessory::Scored_Triangle> b)
                          {return a.first > b.first;});
            }

            ////PRINT
            //std::cout << sorting_stack[0].first << " ";
            //sorting_stack[0].second.print(reffrag,testfrag);

            for (int si = 0; si<sorting_stack.size(); si++){

	        if (parameters.get_rank()){
                    float tmpscore =  sorting_stack[si].first;
	            sorting_stack[si].second.set_vos(tmpscore);
	        }

                // Check if the bonds that are overlapped have the same bond types
                Iso_Acessory::Scored_Triangle see_if_bond_overlap = sorting_stack[si].second;

	        ref_neighbhors_du_1 = 
	            reffrag.mol.get_atom_neighbors(see_if_bond_overlap.three_pairs[0].get_ref());
	        test_neighbhors_du_1 = 
	            testfrag.mol.get_atom_neighbors(see_if_bond_overlap.three_pairs[0].get_test());


	        ref_neighbhors_du_2 = 
	            reffrag.mol.get_atom_neighbors(see_if_bond_overlap.three_pairs[1].get_ref());
	        test_neighbhors_du_2 = 
	            testfrag.mol.get_atom_neighbors(see_if_bond_overlap.three_pairs[1].get_test());
       
                bond_id_ref_first_1 = reffrag.mol.get_bond(see_if_bond_overlap.three_pairs[0].get_ref(),
                                                               ref_neighbhors_du_1[0]);
                bond_id_test_first_1 = testfrag.mol.get_bond(see_if_bond_overlap.three_pairs[0].get_test(),
                                                                 test_neighbhors_du_1[0]);

                bond_id_ref_first_2 = reffrag.mol.get_bond(see_if_bond_overlap.three_pairs[1].get_ref(),
                                                               ref_neighbhors_du_2[0]);
                bond_id_test_first_2 = testfrag.mol.get_bond(see_if_bond_overlap.three_pairs[1].get_test(),
                                                             test_neighbhors_du_2[0]);

                 std::vector<std::string> ref_bond_types {reffrag.mol.bond_types[bond_id_ref_first_1],
                                                          reffrag.mol.bond_types[bond_id_ref_first_2]};


                 std::vector<std::string> test_bond_types {testfrag.mol.bond_types[bond_id_test_first_1],
                                                           testfrag.mol.bond_types[bond_id_test_first_2]};

                // If they do not have the same bond types, return unalignable
                if (ref_bond_types[0] != test_bond_types[0]) {
                    continue;
                }

                if (ref_bond_types[1] != test_bond_types[1]) {
                    continue;
                }
                return sorting_stack[si].second;
            } 
            int size_of_vec =  sorting_stack.size();
            throw size_of_vec;

        } catch(int error_code) {
            //if (error_code == 0) {
            //    std::cout <<  "ERROR_"<< error_code <<": " 
            //              << testfrag.mol.title
            //              <<" is not alignable to referece."
            //              <<" Size of the linker return stack is: "<<error_code<<"." 
            //              << std::endl;
	    //}
            return notalignable; 
        }

    //if SCAFFOLDS for 3 du atoms and above
    }else if ( reffrag.get_num_du() == testfrag.get_num_du() && 
               du_du_pairs.size() == (reffrag.get_num_du() * testfrag.get_num_du()) ){

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
        single_du_du.clear_triangle();

	//get info on neighboring atoms
	std::vector <int> ref_neighbhors_du_1 = 
	     reffrag.mol.get_atom_neighbors(vec_three_dus_dus[0].three_pairs[0].get_ref());
	std::vector <int> test_neighbhors_du_1 = 
	    testfrag.mol.get_atom_neighbors(vec_three_dus_dus[0].three_pairs[0].get_test());

	std::vector <int> ref_neighbhors_du_2 = 
	     reffrag.mol.get_atom_neighbors(vec_three_dus_dus[0].three_pairs[1].get_ref());
	std::vector <int> test_neighbhors_du_2 = 
	    testfrag.mol.get_atom_neighbors(vec_three_dus_dus[0].three_pairs[1].get_test());

        std::vector <int> ref_neighbhors_du_3 = 
             reffrag.mol.get_atom_neighbors(vec_three_dus_dus[0].three_pairs[2].get_ref());
        std::vector <int> test_neighbhors_du_3 = 
            testfrag.mol.get_atom_neighbors(vec_three_dus_dus[0].three_pairs[2].get_test());
        
        // Filter out if the bond types for the sidechains between du and heavy atom  
        // are the same. if not return nothign 

        int bond_id_ref_first_1  =  reffrag.mol.get_bond(vec_three_dus_dus[0].three_pairs[0].get_ref(),
                                                         ref_neighbhors_du_1[0]);
        int bond_id_test_first_1 = testfrag.mol.get_bond(vec_three_dus_dus[0].three_pairs[0].get_test(),
                                                         test_neighbhors_du_1[0]);

        int bond_id_ref_first_2  =  reffrag.mol.get_bond(vec_three_dus_dus[0].three_pairs[1].get_ref(), 
                                                         ref_neighbhors_du_2[0]);
        int bond_id_test_first_2 = testfrag.mol.get_bond(vec_three_dus_dus[0].three_pairs[1].get_test(),
                                                         test_neighbhors_du_2[0]);

        int bond_id_ref_first_3  =  reffrag.mol.get_bond(vec_three_dus_dus[0].three_pairs[2].get_ref(),
                                                         ref_neighbhors_du_3[0]);
        int bond_id_test_first_3 = testfrag.mol.get_bond(vec_three_dus_dus[0].three_pairs[2].get_test(),
                                                         test_neighbhors_du_3[0]);


        std::vector<std::string> ref_bond_types {reffrag.mol.bond_types[bond_id_ref_first_1],
                                                 reffrag.mol.bond_types[bond_id_ref_first_2],
                                                 reffrag.mol.bond_types[bond_id_ref_first_3]};


        std::vector<std::string> test_bond_types {testfrag.mol.bond_types[bond_id_test_first_1],
                                                  testfrag.mol.bond_types[bond_id_test_first_2],
                                                  testfrag.mol.bond_types[bond_id_test_first_3]};

        // if there is absolutely no match in bond types for
        // at least one bond, you must return unalignable.
        for (unsigned int i =0; i<ref_bond_types.size(); i++){
            int check = 0;
            for (unsigned int j =0; j<test_bond_types.size(); j++){
                if (ref_bond_types[i] == test_bond_types[j]){
                    check++;
                }
            }
            if (check == 0){return notalignable;}
        }
    
        //loop through vec_three_du_du where we can assess if the 
        //distances of corresponding du-du distances between two mol
        //are tolerable 
        std::vector <bool> vec_of_bool {};
        for (unsigned int i = 0;i < vec_three_dus_dus.size(); i++){
            //too check if the du-du distances are at least matching
            du_more_than = get_diff_dist_two_mol(reffrag,testfrag,vec_three_dus_dus[i]
						,parameters.get_dist_du_du_scf(),false);
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
                //std::cout << "ERROR "<< error_code << ": "
                //          <<" Aligning " <<testfrag.mol.title << " and" << reffrag.mol.title
                //          <<" are giving an error...REASON:\n" 
                //          <<" std::vector<std::pair<float,std::pair<int,int>>> Iso_Align::select_three_atoms_pairs():"  
                //          <<" error in matching vector size between 'vec_of_bool' and 'vec_three_dus_dus'."
                //          <<" returning empty scaffold vector..." 
                //          << std::endl;
                return notalignable;
            }else if (error_code == 1){
                //std::cout << "ALIGNABLE " << error_code <<": " 
                //          << testfrag.mol.title 
                //          <<" All Du-Du distances between ref and testmol"
                //          <<" are tolerable and are alignable."
                //          << std::endl;
            }else if (error_code == 2){
                //std::cout << "WARNING " << error_code <<": " 
                //          << testfrag.mol.title << " is could be differet from " 
                //          << reffrag.mol.title  << "."
                //          <<" Some Du-Du distances could be unalignable."
                //          <<" Double check on the alignments."
                //          << std::endl; 
            }else if (error_code == 3){
                //std::cout << "ERROR " << error_code <<": " 
                //          << testfrag.mol.title << " is too differet from " 
                //          << reffrag.mol.title  << "."
                //          <<" Du-Du distances are too different. Can't be aligned."
                //          <<" Returning empty vector."
                //          << std::endl;  
                return notalignable;
            }
        }



        //this is where we align and check.
        for (unsigned int i = 0; i < vec_three_dus_dus.size(); i++){

            //if the element is false, align!
            if (!vec_of_bool[i]){

                should_align.first         = false;
                should_align.second.second = 0;

                should_align   = get_aliscore_and_du_angles(reffrag,testfrag,vec_three_dus_dus[i]
                                                    ,parameters.get_bond_angle_tol_scf()
                                                    ,parameters.get_dist_du_du_inter());
            
                //if du_du angles are too faraway SKIP!
                if (!should_align.first){ 
                    vec_three_dus_dus[i].is_not_alignable();
                    continue;
                }
                vec_three_dus_dus[i].is_alignable();
                float iso_score = should_align.second.first.get_score(parameters.get_iso_score_sel());
                sorting_stack.push_back(std::make_pair(iso_score,vec_three_dus_dus[i]));
            }
            vec_three_dus_dus[i].is_not_alignable();
        }


        //try check, if the sorting_stack is 0, throw 0
        try {

            //deactivate each heavy atom 
            //for (unsigned int i =0; i<reffrag.mol.num_atoms;i++){
            //    reffrag.mol.atom_active_flags[i] = false;
            //}
            //for (unsigned int i =0; i<testfrag.mol.num_atoms;i++){
            //    testfrag.mol.atom_active_flags[i] = false;
            //}

            if (sorting_stack.empty()) { throw 0;}

            //sort tmp stack to get the best hrmsd score
            if (parameters.get_iso_score_sel() == "hrmsd_score" || parameters.get_iso_score_sel() == "hms_score"){
                std::sort(sorting_stack.begin(),sorting_stack.end()
                          ,[](std::pair<float,Iso_Acessory::Scored_Triangle>  a
                          ,std::pair<float,Iso_Acessory::Scored_Triangle> b)
                          {return a.first < b.first;});
            } else {
                std::sort(sorting_stack.begin(),sorting_stack.end()
                          ,[](std::pair<float,Iso_Acessory::Scored_Triangle>  a
                          ,std::pair<float,Iso_Acessory::Scored_Triangle> b)
                          {return a.first > b.first;});
            }

            for (int si = 0; si < sorting_stack.size(); si++){

	        if (parameters.get_rank()){
	            //set the hrmsd score for later ranking
                    float tmpscore = sorting_stack[si].first;
	            sorting_stack[si].second.set_vos(tmpscore);
	        }

                // Check if the bonds that are overlapped have the same bond types
                Iso_Acessory::Scored_Triangle see_if_bond_overlap = sorting_stack[si].second;

	        ref_neighbhors_du_1 = 
	            reffrag.mol.get_atom_neighbors(see_if_bond_overlap.three_pairs[0].get_ref());
	        test_neighbhors_du_1 = 
	            testfrag.mol.get_atom_neighbors(see_if_bond_overlap.three_pairs[0].get_test());

	        ref_neighbhors_du_2 = 
	            reffrag.mol.get_atom_neighbors(see_if_bond_overlap.three_pairs[1].get_ref());
	        test_neighbhors_du_2 = 
	            testfrag.mol.get_atom_neighbors(see_if_bond_overlap.three_pairs[1].get_test());

	        ref_neighbhors_du_3 = 
	            reffrag.mol.get_atom_neighbors(see_if_bond_overlap.three_pairs[2].get_ref());
	        test_neighbhors_du_3 = 
	            testfrag.mol.get_atom_neighbors(see_if_bond_overlap.three_pairs[2].get_test());
       
                bond_id_ref_first_1 = reffrag.mol.get_bond(see_if_bond_overlap.three_pairs[0].get_ref(),
                                                               ref_neighbhors_du_1[0]);
                bond_id_test_first_1 = testfrag.mol.get_bond(see_if_bond_overlap.three_pairs[0].get_test(),
                                                                 test_neighbhors_du_1[0]);

                bond_id_ref_first_2 = reffrag.mol.get_bond(see_if_bond_overlap.three_pairs[1].get_ref(),
                                                               ref_neighbhors_du_2[0]);
                bond_id_test_first_2 = testfrag.mol.get_bond(see_if_bond_overlap.three_pairs[1].get_test(),
                                                             test_neighbhors_du_2[0]);

                bond_id_ref_first_3 = reffrag.mol.get_bond(see_if_bond_overlap.three_pairs[2].get_ref(),
                                                               ref_neighbhors_du_3[0]);
                bond_id_test_first_3 = testfrag.mol.get_bond(see_if_bond_overlap.three_pairs[2].get_test(),
                                                             test_neighbhors_du_3[0]);

                std::vector<std::string> ref_bond_types {reffrag.mol.bond_types[bond_id_ref_first_1],
                                                         reffrag.mol.bond_types[bond_id_ref_first_2],
                                                         reffrag.mol.bond_types[bond_id_ref_first_3]};


                std::vector<std::string> test_bond_types {testfrag.mol.bond_types[bond_id_test_first_1],
                                                          testfrag.mol.bond_types[bond_id_test_first_2],
                                                          testfrag.mol.bond_types[bond_id_test_first_3]};

                // If they do not have the same bond types on the bond overlap, return unalignable
                if (ref_bond_types[0] != test_bond_types[0] ) {continue;}

                if (ref_bond_types[1] != test_bond_types[1] ) {continue;}

                if (ref_bond_types[2] != test_bond_types[2] ) {continue;}

                return sorting_stack[si].second;
            }

            int size_of_vec =  sorting_stack.size();
            throw size_of_vec;

        }catch (int size_of_stack){
            //std::cout << testfrag.mol.title 
            //          << " is not alignable to reference. Size of stack is: " 
            //          << size_of_stack 
            //          << std::endl;
            return notalignable;  
        }

    }else if (reffrag.get_num_du() != testfrag.get_num_du()){
        //std::cout << testfrag.mol.title 
        //          << " is not alignable to reference due to different num of attachment points. " 
        //          << "refmol num of Du: " << reffrag.get_num_du() << " testmol num of Du: " << testfrag.get_num_du()  
        //          << " Returning an empty vector."
        //          << std::endl; 
        return notalignable;
    } 


    return notalignable;
}

//get the legnths for the twoD atom plots
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
    b1 = twoD_length({half_sq[0].first,y2}, {half_sq[0].first, 0});
    
    //b2 is the lenght of the right of the trapzoid
    b2 = twoD_length({half_sq[3].first,y1},{half_sq[3].first, 0}); 
    area1 = (b1 + b2)/2 * h1;

    //calculate the area of the trianlge that sits on top of the trapezoid
    b3 = twoD_length(half_sq[1],{bin,y1});
    h2 = twoD_length(half_sq[1],{bin-1,y2});  
    area2 = (b3 * h2)/2;

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
    } 
    return half_sq_coordinates;
}

//Where the intial step of iso_align protocol is called.
//1st step -> get radial distances
//2nd step -> get the 2D atom plots for each atom within the internal structure
//3rd step -> calculate the similarity between each atom freq_dist 
//4th step -> align
std::vector <float> get_twoD_atom_plots(std::vector<float> & vec_of_centers){

    //these are set min and max bins
    float max_bin = 6.0;
    float min_bin = 0.0;

    //structures to hold triangle areas
    std::pair <float,float> tmp_pair {};

    //The vector that contians the freq distributions 
    std::vector<float> atom_plot(int(max_bin),0.0);

    std::vector<std::pair <float,float>> single_tri_coords {};

    
    std::vector<float> vec_of_center_temp = vec_of_centers;

    for (unsigned int i = 0; i<vec_of_center_temp.size(); i++){
        single_tri_coords = get_plot_triangle_coord(vec_of_center_temp[i],min_bin,max_bin);  
        for (int z=0; z<int(max_bin); z++){
            if (single_tri_coords[1].first == z) {
                atom_plot[z-1] = atom_plot[z-1] + 0.5;         
                atom_plot[z] = atom_plot[z] + 0.5;
                continue;
            } 

            if (single_tri_coords[1].first > z-1 && single_tri_coords[1].first < z ) {
                float  left_area=0;
                float   mid_area=0;
                float right_area=0;

		//calculate the mid_area that is shaped as a pentagon.
                mid_area =  area_pentagon(single_tri_coords,z);  

                //check if there is a triangle in the bin before or after the current bin number
                tmp_pair   = area_triangle(single_tri_coords,z);
                left_area  = tmp_pair.first;
                right_area = tmp_pair.second;
 
                atom_plot[z-1] = atom_plot[z-1] + mid_area;
		if (z != 1){
		    atom_plot[z-2] = atom_plot[z-2] + left_area;
		}
                atom_plot[z]   =   atom_plot[z] + right_area;

            }
        }
	single_tri_coords.clear();
    }

    single_tri_coords.clear();
    vec_of_center_temp.clear();

    return atom_plot; 
};


std::vector<float> 
Iso_Align::get_a_radial_dist(Fragment & frag, unsigned int atom_num){

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
    
    list_of_radial_dist.clear();
    return return_radial_areas_vec;
}

void
Iso_Align::get_all_radial_dist(Fragment& frag)
{

    std::vector<float> list_of_radial_dist {};
    float tmplength = 0;
    std::vector<float> tmp_radial{};


    for (int origin=0; origin<frag.mol.num_atoms;origin++){
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
Iso_Align::align(Fragment& reffrag, std::vector<Fragment>& vec_testfrags, 
                 Iso_Parm  parameters){

    int num_frags =                 vec_testfrags.size();
    std::vector<Iso_Acessory::Scored_Triangle>    three_atoms_pairs {};
    std::vector<std::vector<float>> tmp_centroid_coor {}; 

    //calculate the radial dist of ref frag
    if (!reffrag.is_it_radial_set() && !reffrag.alloc_set){
        reffrag.allocate_radial_dist_distri();
        reffrag.calc_num_du();
        reffrag.calc_radial_dist_distri();
    }

    //calculate the radial dist in test frags
    for (unsigned int i = 0; i<vec_testfrags.size(); i++){
        if (!vec_testfrags[i].is_it_radial_set() && !vec_testfrags[i].alloc_set){
            vec_testfrags[i].allocate_radial_dist_distri();
            vec_testfrags[i].calc_num_du();
            vec_testfrags[i].calc_radial_dist_distri();
        }
    }

    //returns best triagnle from the three atom pairs
    for (unsigned int i =0; i<num_frags; i++){
        Iso_Acessory::Scored_Triangle best_triangle;
        best_triangle = get_three_atoms_pairs(reffrag,vec_testfrags[i],parameters);


	//If empty, there are no viable alignments, skip
        if(best_triangle.three_pairs.empty()) {
            vec_testfrags[i].is_not_iso_aligned();
	    vec_testfrags[i].set_iso_score(9998);
            continue;
        }    
	//marks the frags as iso alignable
        vec_testfrags[i].is_iso_aligned();

	//calculates the centroids
        tmp_centroid_coor = get_centroid(reffrag.mol
                                        ,vec_testfrags[i].mol
                                        ,best_triangle
                                        ,true);
	//3D aligns the molecules
        align_molcentroid(reffrag.mol,vec_testfrags[i].mol,tmp_centroid_coor,best_triangle);

	//if rank is turned on, set teh hrmsd score
	if (parameters.get_rank()){
            vec_testfrags[i].set_iso_score(best_triangle.get_vos());	
	}

        best_triangle.clear_triangle();
    }

    if (parameters.get_rank()){
	//sort the pairs based on cos_similarity
	//sort the fragments based on their iso_frag hrmsd
    	Iso_Align::frag_sort(vec_testfrags
    	         ,Iso_Align::fragment_sort); 
    }

};

void
Iso_Align::align(Fragment& reffrag, std::vector<Fragment>& vec_testfrags, 
                 Iso_Parm  parameters, std::vector<Iso_Acessory::Scored_Triangle> & best_triangles){

    int num_frags =                 vec_testfrags.size();
    std::vector<Iso_Acessory::Scored_Triangle>    three_atoms_pairs {};
    std::vector<std::vector<float>> tmp_centroid_coor {}; 

    //calculate the radial dist of ref frag
    if (!reffrag.is_it_radial_set() && !reffrag.alloc_set){
        reffrag.allocate_radial_dist_distri();
        reffrag.calc_num_du();
        reffrag.calc_radial_dist_distri();
    }

    //calculate the radial dist in test frags
    for (unsigned int i = 0; i<vec_testfrags.size(); i++){
        if (!vec_testfrags[i].is_it_radial_set() && !vec_testfrags[i].alloc_set){
            vec_testfrags[i].allocate_radial_dist_distri();
            vec_testfrags[i].calc_num_du();
            vec_testfrags[i].calc_radial_dist_distri();
        }
    }
    //returns best triagnle from the three atom pairs
    for (unsigned int i =0; i<num_frags; i++){
        Iso_Acessory::Scored_Triangle best_triangle;
        best_triangle = get_three_atoms_pairs(reffrag,vec_testfrags[i],parameters);

        vec_testfrags[i].best_tri = best_triangle;
        best_triangles.push_back(best_triangle);


	//If empty, there are no viable alignments, skip
        if(best_triangle.three_pairs.empty()) {
            vec_testfrags[i].is_not_iso_aligned();
	    vec_testfrags[i].set_iso_score(9998);
            continue;
        }    
	//marks the frags as iso alignable
        vec_testfrags[i].is_iso_aligned();

	//calculates the centroids
        tmp_centroid_coor = get_centroid(reffrag.mol
                                        ,vec_testfrags[i].mol
                                        ,best_triangle
                                        ,true);
	//3D aligns the molecules
        align_molcentroid(reffrag.mol,vec_testfrags[i].mol,tmp_centroid_coor,best_triangle);

	//set teh score
        vec_testfrags[i].set_iso_score(best_triangle.get_vos());	

        best_triangle.clear_triangle();
    }

    if (parameters.get_rank()){

        Iso_Score::Score tmp_score;

        for (unsigned int i=0; i<vec_testfrags.size(); i++){
            tmp_score.calc_score(reffrag,vec_testfrags[i]);
            vec_testfrags[i].set_iso_score(
                tmp_score.get_score(parameters.get_iso_rank_score_sel()));
        }       

        if (parameters.get_iso_rank_reverse()==false){
            if (parameters.get_iso_rank_score_sel()== "hrmsd_score" || parameters.get_iso_rank_score_sel()== "hms_score"){
	        //sort the pairs based on cos_similarity
	        //sort the fragments based on their iso_frag hrmsd
    	        Iso_Align::frag_sort(vec_testfrags
                         ,Iso_Align::fragment_sort_reverse); 
            }else{
    	        Iso_Align::frag_sort(vec_testfrags
    	                 ,Iso_Align::fragment_sort); 
            }
        }else{
            if (parameters.get_iso_rank_score_sel()== "hrmsd_score" || parameters.get_iso_rank_score_sel()== "hms_score"){
	        //sort the pairs based on cos_similarity
	        //sort the fragments based on their iso_frag hrmsd
    	        Iso_Align::frag_sort(vec_testfrags
                         ,Iso_Align::fragment_sort); 
            }else{
    	        Iso_Align::frag_sort(vec_testfrags
    	                 ,Iso_Align::fragment_sort_reverse); 
            }
        }
    }
  
};
// +++++++++++++++++++++++++++++++++++++++++
// Activate all of the atoms AND bonds in a given DOCKMol
void
activate_mol( DOCKMol & mol )
{
     // Iterate through all atoms and set atom_active_flag to true
     for (int i=0; i<mol.num_atoms; i++){
               mol.atom_active_flags[i] = true;
     }    

     // Iterate through all bonds and set bond_active_flag to true
     for (int i=0; i<mol.num_bonds; i++){
               mol.bond_active_flags[i] = true;
     }    

    return;

} // end DN_Build::activate_mol()

// Function to find the possible
// permutations
void 
Iso_Align::permutations( std::vector<std::vector<Iso_Acessory::A_Pair> >& res,
                   std::vector< Iso_Acessory::A_Pair > nums, int l, int h)
{
    // Base case
    // Add the vector to result and return
    if (l == h) {
        res.push_back(nums);
        return;
    }
 
    // Permutations made
    for (int i = l; i <= h; i++) {
 
        // Swapping
        swap(nums[l], nums[i]);
 
        // Calling permutations for
        // next greater value of l
        permutations(res, nums, l + 1, h);
 
        // Backtracking
        swap(nums[l], nums[i]);
    }
}

// Function to get the permutations
std::vector<std::vector<Iso_Acessory::A_Pair> > 
Iso_Align::permute( std::vector< Iso_Acessory::A_Pair >& nums)
{
    // Declaring result variable
    std::vector<std::vector<Iso_Acessory::A_Pair> > res;
    int x = nums.size() - 1;
 
    // Calling permutations for the first
    // time by passing l
    // as 0 and h = nums.size()-1
    permutations(res, nums, 0, x);
    return res;
}


Fragment
Iso_Align::align_two_frags(Fragment reffrag, Fragment testfrag){

    Hungarian_RMSD h;    
    std::vector<std::vector<Iso_Acessory::A_Pair> > all_poss_pairs
        = permute(testfrag.best_tri.three_pairs);

    std::vector<Fragment> poss_frags {};
    std::vector<float> poss_hrmsd {};

    for ( int i =0; i < all_poss_pairs.size(); i++){
 
        Fragment copy_ref = reffrag;
        Fragment copy_test = testfrag;
 
        Iso_Acessory::Scored_Triangle BEST_TRI;
        BEST_TRI.three_pairs = all_poss_pairs[i];

        //Grab both centroids
        std::vector<std::vector<float>> 
            centroid_coor = get_centroid(reffrag.mol
                                    ,testfrag.mol
                                    ,BEST_TRI
                                    ,true);

        align_molcentroid(copy_ref.mol,copy_test.mol,centroid_coor, BEST_TRI);

        std::pair <float, int> h_result = 
            h.calc_Hungarian_RMSD_dissimilar(reffrag.mol,copy_test.mol); 

        poss_frags.push_back( copy_test );
        poss_hrmsd.push_back( h_result.first );
  
    }

    std::vector<int> indices(poss_hrmsd.size());
    std::iota(indices.begin(),indices.end(),0);
    std::sort(indices.begin(), indices.end(),
              [&](int A, int B) -> bool {
                   return poss_hrmsd[A]
                        < poss_hrmsd[B];
               });  

    return poss_frags[indices[0]];

}
