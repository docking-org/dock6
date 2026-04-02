#ifndef RDDRIVE_H
#define RDDRIVE_H

#include "dockmol.h"
#include "rdtyper.h"
 
class RDTYPER;

namespace RD_Parm{
    class RD_Drive_Parm{

    private:
        std::map<unsigned int, double>    fragMap;
        bool                              create_smiles;
        std::map<std::string,std::string> PAINSmap;
        int                               upper_stereocenter;

        float                             lower_tpsa;
        float                             upper_tpsa;
        float                             tpsa_std_dev;

        float                             lower_clogp;
        float                             upper_clogp;
        float                             clogp_std_dev;

        float                             lower_esol;
        float                             upper_esol;
        float                             esol_std_dev;

        float                             lower_qed;
        float                             qed_std_dev;

        float                             upper_sa;
        float                             sa_std_dev;

        int                               upper_pains;

    public:
        void                              set_fragMap( std::map<unsigned int, double> ); 
        void                              set_createsmiles( bool );  
        void                              set_PAINSmap( std::map<std::string,std::string> );
        void                              set_upper_stereocenter( int );

        void                              set_lower_tpsa( float );
        void                              set_upper_tpsa( float );
        void                              set_tpsa_std_dev( float );

        void                              set_lower_clogp( float );
        void                              set_upper_clogp( float );
        void                              set_clogp_std_dev( float );

        void                              set_lower_esol( float );
        void                              set_upper_esol( float );
        void                              set_esol_std_dev( float );

        void                              set_lower_qed( float );
        void                              set_qed_std_dev( float );

        void                              set_upper_sa( float );
        void                              set_sa_std_dev( float );

        void                              set_upper_pains( int );

        std::map<unsigned int, double>    get_fragMap();
        bool                              get_createsmiles();
        std::map<std::string,std::string> get_PAINSmap();
        int                               get_upper_stereocenter();

        float                             get_lower_tpsa();
        float                             get_upper_tpsa();
        float                             get_tpsa_std_dev();

        float                             get_lower_clogp();
        float                             get_upper_clogp();
        float                             get_clogp_std_dev();

        float                             get_lower_esol();
        float                             get_upper_esol();
        float                             get_esol_std_dev();

        float                             get_lower_qed();
        float                             get_qed_std_dev();

        float                             get_upper_sa();
        float                             get_sa_std_dev();

        int                               get_upper_pains( );
    
                                          RD_Drive_Parm();
                                          ~RD_Drive_Parm();
    };
};
namespace RD_Calc{
    namespace{
        class Calc{
            private:
                RDTYPER rd_typer;

            public:
                std::map<unsigned int, double>    defFragMap;
                bool                              defCreate_smiles;
                std::map<std::string,std::string> defPAINSmap;


                void calc_des( DOCKMol &, RD_Parm::RD_Drive_Parm );
                Calc();
                ~Calc();

        };

        class Drive{

            private:
                bool         stereocenter_cutoff( DOCKMol, int = 2);

                bool                 tpsa_cutoff( DOCKMol, float = 28.53, 
                                                           float = 113.20, 
                                                           float = 42.33);

                bool                clogp_cutoff( DOCKMol, float = -0.30,
                                                           float = 3.75,
                                                           float = 2.02);

                bool                 esol_cutoff( DOCKMol, float = -5.23,
                                                           float = -1.35, 
                                                           float = 1.94);

                bool                  qed_cutoff( DOCKMol, float = 0.61,
                                                           float = 0.19);

                bool                   sa_cutoff( DOCKMol, float = 3.34,
                                                           float = 0.9);

                bool                pains_cutoff( DOCKMol, float = 1 );
     
            public:
                bool                       drive_growth( DOCKMol, RD_Parm::RD_Drive_Parm );
 

        };

    }
};


class RD_Drive {

    private:
        RD_Calc::Calc     calc_des;
        RD_Calc::Drive    drive_des;

    public:
        bool              drive( DOCKMol);
        bool              drive( DOCKMol, RD_Parm::RD_Drive_Parm&);



};

#endif
