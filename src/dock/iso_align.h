#ifndef ISO_ALIGN_H
#define ISO_ALIGN_H

#include "fragment.h"

#include <vector>
#include <iostream>
#include <string>

class Fragment;


namespace Iso_Score{

    namespace{

        class Domain{
            private:
                float          overlap;
                float      overlap_hvy;
                float      overlap_pho;
                float      overlap_phi;
                float      overlap_neg;
                float      overlap_pos;
            public:
                float      get_overlap();
                float      get_overlap_hvy();
                float      get_overlap_pho();
                float      get_overlap_phi();
                float      get_overlap_neg();
                float      get_overlap_pos();

                void           set_overlap(float);
                void       set_overlap_hvy(float);
                void       set_overlap_pho(float);
                void       set_overlap_phi(float);
                void       set_overlap_neg(float);
                void       set_overlap_pos(float); 
                Domain();
  
        };

        class Score{
            private:
                float                   score;
                float                vo_score;

            public:
                float         get_score();
                float         get_vo_score();
                void          set_score(float);
                void          set_vo_score(float);
                void          calc_score(Fragment&,Fragment&);
        };
    }
    float    Analytical_method(Fragment&, Fragment&);

}


namespace Iso_Acessory{
    namespace{
        class A_Pair: private std::pair<int,int>{
        
            public:
        
            private:
                std::string      title;
                     float      cos_sim;
        
            public:
                int              get_ref();
                int              get_test();
                void             set_ref(int);
                void             set_test(int);
                std::string      get_title();
                void             set_title(std::string);
        	void             print_all(Fragment&, Fragment&);
                void             print_all();
                void            print_ref();
                void            print_test();
                float              get_cos_sim();
                void                set_cos_sim(float);
                void             clear();
        };
        
        class Scored_Triangle{
        
            private:
                bool                alignable;
        
            public:
                void                push_back(A_Pair&);
                void                pop_back();
                void                print();
                void                print(Fragment&, Fragment &);
		void		    clear();
                bool                check_if_redundant();
                void                is_not_alignable();
                void                is_alignable();
                std::vector<A_Pair>  three_pairs;
                Scored_Triangle();
        };


    }


}

class Iso_Parm {

    private:
        float               bond_angle_tol_sid;
        float               bond_angle_tol_lnk;
        float               bond_angle_tol_scf;
        float               dist_tol_sid;
        float               dist_tol_lnk;
        float               dist_tol_scf;
    public:
        void                set_bond_angle_tol_sid(float);
        void                set_bond_angle_tol_lnk(float);
        void                set_bond_angle_tol_scf(float);
        float               get_bond_angle_tol_sid();
        float               get_bond_angle_tol_lnk();
        float               get_bond_angle_tol_scf();

        void                set_dist_tol_sid(float);
        void                set_dist_tol_lnk(float);
        void                set_dist_tol_scf(float);
        float               get_dist_tol_sid();
        float               get_dist_tol_lnk();
        float               get_dist_tol_scf();

};

class Iso_Align {

    public: 
        void                align(Fragment&, std::vector<Fragment>&, Iso_Parm);
        void                get_all_radial_dist(Fragment&);
        std::vector<float>  get_a_radial_dist(Fragment &, unsigned int);
        
};

float                                calc_2atoms_length(Fragment, int , int );
std::vector <float>                  get_twoD_atom_plot(std::vector<float>);
std::vector<std::pair <float,float>> get_plot_triangle_coord(float,float,float);
std::pair<float,float>               area_triangle(std::vector<std::pair<float,float>>,float);
float                                area_pentagon(std::vector<std::pair <float,float>> ,float);
float                                twoD_length(std::pair<float,float>,std::pair<float, float>);
float                                calc_cos_similarity(Fragment&,Fragment&,int,int);

void                                 align_molcentroid (Fragment, Fragment testmol
                                                        ,std::vector<std::vector<float>>
                                                        ,Iso_Acessory::Scored_Triangle);

std::vector<std::vector<float>>      get_centroid(Fragment&,Fragment&
                                                  ,Iso_Acessory::Scored_Triangle,bool);
Iso_Acessory::Scored_Triangle        get_three_atoms_pairs(Fragment&,Fragment&); 

bool                                 get_diff_dist_two_mol(Fragment&, Fragment&
                                                          ,Iso_Acessory::Scored_Triangle,float,bool);


std::pair<bool,std::pair<double,int>> get_hrmsd_and_du_angles(Fragment&, Fragment &
                                                                ,Iso_Acessory::Scored_Triangle);

float                                  magnitude(std::vector<float>);
float                                  dot_product(std::vector<float>,std::vector<float>);
bool                                   check_H_vs_atom(Fragment&,Fragment&);
float                                twoD_length(std::pair<float,float>, std::pair<float, float>);
float                                area_pentagon(std::vector<std::pair <float,float>>, float);
float                                calc_2atoms_length(Fragment , int , int);
bool          are_Du_axis_overlapped(DOCKMol &, DOCKMol &
                                     ,Iso_Acessory::Scored_Triangle
                                     ,float );
float get_vector_angle(const std::vector<float> & , const std::vector<float> & );






#endif
