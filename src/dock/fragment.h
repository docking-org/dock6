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

// +++++++++++++++++++++++++++++++++++++++++
// Attachment points are the dummy atom + heavy atom it is connected to
class           AttPoint {
   public:
       int    dummy_atom;
       int    heavy_atom;

      AttPoint();
      ~AttPoint();
};


// +++++++++++++++++++++++++++++++++++++++++
// Scaffolds, Linkers, Sidechains, and Rigids are condensed into one class called Fragment
class           Fragment {

    private:
        std::vector <float>     *radial_dist_distri;
        int                     num_du;
        bool                    iso_aligned;

    public:
        void                    is_iso_aligned();
        void                    is_not_iso_aligned();
        bool                    is_it_iso_aligned();
        void                    calc_radial_dist_distri();
        void                    set_radial_dist_distri(int,std::vector<float>);
        std::vector<float>      get_radial_dist_distri(int);
        void                    print_radial_dist_distri();
        void                    allocate_radial_dist_distri();
        void                    clear_radial_dist_distri();
        void                    calc_num_du();
        int                     get_num_du();



   public:
  
       DOCKMol                          mol;
       std::vector <AttPoint>           aps;
       bool                             used;
       float                            tanimoto;
       int                              last_ap_heavy;
       int                              scaffolds_this_layer;
       std::vector < std::pair < int, int > >  torenv_recheck_indices;
       std::vector < DOCKMol >                 frag_growth_tree;
       // Need to update the DN code to utilize this information when comparing bonds - CS: 09/19/16
       std::vector < std::pair < int, std::string > > aps_bonds_type; // Bond number and type for each aps

       // Parameters for scaffold hopping
       int                              size;
       int                              ring_size;
       std::vector <float>              aps_cos;
       std::vector <std::vector <float> > aps_dist; // Pairwise distances between all aps
        
       // To keep track of mutations
       int                              mut_type;

       /*#ifdef BUILD_DOCK_WITH_RDKIT
       // Descriptor-driven de novo variables
       bool                             fail_clogp;
       bool                             fail_esol;
       bool                             fail_qed;
       bool                             fail_sa;
       bool                             fail_stereo;
       int                              dropped_at_layer;
       #endif
       */

       // Initialize fragment object
       Fragment              read_mol( Fragment &);
       
       Fragment();
       ~Fragment();
};

#endif //FRAGMENT_H
