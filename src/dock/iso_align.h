#ifndef ISO_ALIGN_H
#define ISO_ALIGN_H

#include "fragment.h"
#include <vector>
#include <iostream>
#include <string>
#include <functional>


class Fragment;

// Atomic weights from General Chemistry 3rd Edition Darrell D. Ebbing
#define ATOMIC_WEIGHT_H 1.00794
#define ATOMIC_WEIGHT_C 12.011
#define ATOMIC_WEIGHT_N 14.00647
#define ATOMIC_WEIGHT_O 15.9994
#define ATOMIC_WEIGHT_S 32.066
#define ATOMIC_WEIGHT_P 30.973762
#define ATOMIC_WEIGHT_F 18.9984032
#define ATOMIC_WEIGHT_Cl 35.4527
#define ATOMIC_WEIGHT_Br 79.904
#define ATOMIC_WEIGHT_I 126.90447


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
                ~Domain();
  
        };

        class Score{
            private:
                float                 hrmsd_score;
                float                   hms_score;
                float                  tani_score;
                float                   vos_score;
                float                       score;
                float             total_component;
                float        heavy_atom_component; 
                float          negative_component;    
                float          positive_component;    
                float       hydrophobic_component; 
                float       hydrophilic_component; 

            public:
                float         get_score(std::string);

                void          set_score(float);

                void          calc_score(Fragment&,Fragment&);
                Score(const Score&);
                Score();
                ~Score();
                void          operator=(const Iso_Score::Score&);
                void          Analytical_method(Fragment&, Fragment&);
  
        };

        std::string   SCORECLUSTTYPE = "heavy_atom_component";
    }

}


namespace Iso_Acessory{

        class A_Pair: private std::pair<int,int>{
        
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
                void             print_ref();
                void             print_test();
                float            get_cos_sim();
                void             set_cos_sim(float);
                void             clear();

                A_Pair();
                ~A_Pair();
                A_Pair(const A_Pair&);
                void operator=(const A_Pair&);
        };

        class Scored_Triangle{
        
            private:
                bool                alignable;
		float               hrmsd;
                float               volume_overlap;
        
            public:
                void                push_back(A_Pair&);
                void                pop_back();
                void                print();
                void                print(Fragment&, Fragment &);
		void		    clear_triangle();
                bool                check_if_redundant();
                void                is_not_alignable();
                void                is_alignable();
                void                set_hrmsd(float);
                float               get_hrmsd();
                void                set_vos(float);
                float               get_vos();
                std::vector<A_Pair> three_pairs;

                Scored_Triangle();
                ~Scored_Triangle();
                Scored_Triangle(const Scored_Triangle&);
                void operator=(const Scored_Triangle&);
        };

}


namespace Iso_Table{

    class Iso_Tab: private std::vector<std::pair<std::string,
				       std::vector<Fragment>>> {

	private:
            std::vector<std::vector<Iso_Acessory::Scored_Triangle>>
                                                       best_triangles;
            std::vector<std::string>  head_names;
            std::vector<Fragment>     head_frags;
      
	public:
            void                                       write_table(std::string,int);
	    void                                       set(Fragment,std::vector<Fragment>&); 
	    void 	                               print_table();
	    void                                       set_tri(Fragment,std::vector<Iso_Acessory::Scored_Triangle>);
	    int                                        get_size();
	    bool                                       check_if(Fragment);
            std::string                                get_string(int);
            std::vector<Fragment>                      get(Fragment,int);
            Fragment                                   get_head(Fragment);
            std::vector<Iso_Acessory::Scored_Triangle> get_tri(Fragment);
            void                                       clear_table();
              int                                      get_size_total();
            Iso_Tab();
            ~Iso_Tab();
            Iso_Tab(const Iso_Tab&);
            void operator=(const Iso_Tab&);
            
    };
}


class Iso_Parm {

    private:
        float               bond_angle_tol_sid;
        float               bond_angle_tol_lnk;
        float               bond_angle_tol_scf;
        float               dist_du_du_inter;
        float               dist_tol_sid;
        float               dist_tol_lnk;
	float               dist_du_du_lnk;
        float               dist_du_du_scf;
	int                 diff_num_atoms;
	bool                rank;
        bool                write_libraries;
        bool                iso_fraglib;
        std::string         iso_fraglib_path;
        int                 iso_num_top;
        std::string         iso_rank_score_sel;
        std::string         iso_score_sel;
        bool                iso_rank_reverse;
        int                 iso_write_freq_cutoff;
        float               iso_cos_score_cutoff;

    public:
        void                set_bond_angle_tol_sid(float);
        void                set_bond_angle_tol_lnk(float);
        void                set_bond_angle_tol_scf(float);
        void                set_dist_du_du_inter(float);
        void                set_dist_tol_sid(float);
        void                set_dist_tol_lnk(float);
	void 		    set_dist_du_du_lnk(float);
        void                set_dist_du_du_scf(float);
	void                set_diff_num_atoms(int);
	void                set_rank(bool);
        void                set_write_libraries(bool);
        void                set_iso_fraglib(bool,std::string);
        void                set_iso_num_top(int);
        void                set_iso_score_sel(std::string);
        void                set_iso_rank_score_sel(std::string);
        void                set_iso_rank_reverse(bool);
        void                set_iso_write_freq_cutoff(int);
        void                set_iso_cos_score_cutoff(float); 


        float               get_bond_angle_tol_sid();
        float               get_bond_angle_tol_lnk();
        float               get_bond_angle_tol_scf();
        float               get_dist_du_du_inter();
        float               get_dist_tol_sid();
        float               get_dist_tol_lnk();
	float               get_dist_du_du_lnk();
        float               get_dist_du_du_scf();
	bool                get_rank();
        bool                get_write_libraries();
	int                 get_diff_num_atoms();
        int                 get_iso_num_top();
        std::pair<bool,     
              std::vector<std::string>> get_iso_fraglib();
        std::string         get_iso_score_sel();
        std::string         get_iso_rank_score_sel();
        bool                get_iso_rank_reverse();
        int                 get_iso_write_freq_cutoff();
        float               get_iso_cos_score_cutoff(); 
};

class Iso_Align {
    private:
         static bool                     fragment_sort(  Fragment & a,  Fragment & b);
         static bool                     fragment_sort_reverse(  Fragment & a,  Fragment & b);
         void                            frag_sort(std::vector<Fragment> &, std::function<bool(Fragment&,Fragment&)>);
    public: 
        void                             align(Fragment&, std::vector<Fragment>&, Iso_Parm);
        void                             align(Fragment&, std::vector<Fragment>&, Iso_Parm,
                                               std::vector<Iso_Acessory::Scored_Triangle>&);
        Fragment                         align_two_frags(Fragment,Fragment);
        void                             get_all_radial_dist(Fragment&);
        std::vector<float>               get_a_radial_dist(Fragment &, unsigned int);
        
       float                             calc_2atoms_length(Fragment, int , int );
       std::vector <float>               get_twoD_atom_plot(std::vector<float> &);
       std::vector<std::pair <float,
                             float>>     get_plot_triangle_coord(float,float,float);
       std::pair<float,float>            area_triangle(std::vector<std::pair<float,float>>,float);
       float                             area_pentagon(std::vector<std::pair <float,float>>, float);
       float                             twoD_length(std::pair<float,float>, std::pair<float, float>);
       float                             threeD_length(std::vector<float>, std::vector<float>);

       float                             calc_cos_similarity(Fragment&,Fragment&,int,int);
 
       void                              permutations(std::vector<std::vector<Iso_Acessory::A_Pair>>& ,std::vector<Iso_Acessory::A_Pair>,
                                                      int , int );

       std::vector<std::vector<Iso_Acessory::A_Pair>> 
                                         permute( std::vector<Iso_Acessory::A_Pair> &);

       void                              align_molcentroid (DOCKMol&, DOCKMol& 
                                                            ,std::vector<std::vector<float>>
                                                            ,Iso_Acessory::Scored_Triangle);

       std::vector<std::vector<float>>   get_centroid(DOCKMol&,DOCKMol&
                                                         ,Iso_Acessory::Scored_Triangle,bool);
       Iso_Acessory::Scored_Triangle     get_three_atoms_pairs(Fragment&,Fragment&,Iso_Parm); 
       
       bool                              get_diff_dist_two_mol(Fragment&, Fragment&
                                                              ,Iso_Acessory::Scored_Triangle,float,bool);
       float                             get_dist_diff_atat(DOCKMol&, DOCKMol&
                                                            ,int,int);
       
       std::pair<bool,
                 std::pair<
                  Iso_Score::Score,
                           int>>         get_aliscore_and_du_angles(Fragment, Fragment
       						     ,Iso_Acessory::Scored_Triangle, float, float);
       
       float                             magnitude(std::vector<float>);
       float                             dot_product(std::vector<float>,std::vector<float>);
       bool                              check_H_vs_atom(Fragment&,Fragment&);
       bool                              are_Du_axis_overlapped(DOCKMol &, DOCKMol &
                                                                ,Iso_Acessory::Scored_Triangle
                                                                ,float ,float);
       float                            get_vector_angle(const std::vector<float> &,
                                                         const std::vector<float> & );
       float                            get_torsion_angle(const std::vector<float> 
                                                          ,const std::vector<float> 
                                                          ,const std::vector<float> 
                                                          ,const std::vector<float>);
       std::vector<float>               subtract_vec(const std::vector<float> &, const std::vector<float> &);
       std::vector<float>               cross_prod(const std::vector<float> &, const std::vector<float> &);
       void                             normalize_vec(std::vector <float> &);
       float                            length_vec(std::vector<float>);
       std::vector<std::vector<double>> get_rotation_mat(std::vector<float>&, std::vector<float>&);
       void                             target_translation (DOCKMol&,std::vector<float>);
};

#endif
