#ifndef FRAGMENT_H
#define FRAGMENT_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "amber_typer.h"
#include "dockmol.h"
#include "master_score.h"
#include "utils.h"
#include "iso_align.h"

// +++++++++++++++++++++++++++++++++++++++++
// Attachment points are the dummy atom + heavy atom it is connected to
class           AttPoint {
   public:
       int    dummy_atom;
       int    heavy_atom;
       std::string frag_name;
      AttPoint();
      ~AttPoint();
      AttPoint(const AttPoint&);
      void operator=(const AttPoint&);
};


// +++++++++++++++++++++++++++++++++++++++++
// Scaffolds, Linkers, Sidechains, and Rigids are condensed into one class called Fragment
class           Fragment {

    private:
        std::vector <float>             *radial_dist_distri;
        int                             num_du;
        bool                            iso_aligned;
	float                           iso_score;
        bool                            radial_set;

    public:
        void                            is_iso_aligned();
        void                            is_not_iso_aligned();
        bool                            is_it_iso_aligned();
        void                            calc_radial_dist_distri();
        void                            set_radial_dist_distri(int,std::vector<float>);
        std::vector<float>              get_radial_dist_distri(int);
        void                            print_radial_dist_distri();
        void                            allocate_radial_dist_distri();
        void                            clear_radial_dist_distri();
        void                            calc_num_du();
        int                             get_num_du();
	void                            set_iso_score(float);
	float                           get_iso_score();
        bool                            is_it_radial_set();
        bool                            alloc_set;
        Iso_Acessory::Scored_Triangle   best_tri;

        Fragment(const Fragment&);
        void operator=(const Fragment&);

   public:
  
       DOCKMol                          mol;
       std::vector <AttPoint>           aps;
       bool                             used;
       float                            tanimoto;
       int                              last_ap_heavy;
       int                              scaffolds_this_layer;
       std::vector < std::pair < int, int > >  torenv_recheck_indices;
       std::vector < DOCKMol >                 frag_growth_tree;

       void 				calc_mol_wt();
       // Need to update the DN code to utilize this information when comparing bonds - CS: 09/19/16
       std::vector < std::pair < int, std::string > > aps_bonds_type; // Bond number and type for each aps

       // Parameters for scaffold hopping
       int                              size;
       int                              ring_size;
       std::vector <float>              aps_cos;
       std::vector <std::vector <float> > aps_dist; // Pairwise distances between all aps
        
       // To keep track of mutations
       int                              mut_type;
       //
       int                              freq_num;

       // Variables and data for isoswap_GA
       double                           iso_ori_mat[3][3];
       float                            iso_tors_turned;
       int                              iso_targeted_AP;
       int                              iso_frag_att_num_APs;
       int                              iso_head_attached_ind;
       DOCKVector                       iso_centers_com;
       DOCKVector                       iso_spheres_com;

       // To debug issues w/ fragments
       void print(int index, std:: string label);

       // Initialize fragment object
       Fragment              read_mol( Fragment &);
       
       Fragment();
       ~Fragment();
};

#endif //FRAGMENT_H
