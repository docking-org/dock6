/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
9/96
*/

typedef struct score_grid_struct
{
  int flag;			/* Flag for precomputed grid use */
  int flag37;			/* Flag for using grid generated with dock3.7 */
  int flagqnifft;		/* Flag for using QNIFFT grid (True) or DelPhi grid (False)  */
  int init_flag;		/* Flag for grid initialization */
  float version;		/* Chemgrid version format */

  int size;			/* Number of grid points */
  int psize;            /* Number of Delphi grid points */
  int nsize;            /* Number of Delphi grid points in one dimension */
  int span[3];			/* Number of points along each grid edge */
  int pspan[3];	            /* Number of points along each Delphi grid edge */
  float origin[3];		/* Coordinates of origin of grids */
  float porigin[3];           /* Coordinates of center of Delphi grids */
  float spacing;		/* Spacing in angstroms between grid points */
  float pspacing;       /* Spacing in angstroms between Delphi grid points */
  float distance;		/* Longest interaction distance */

  FILE_NAME file_prefix;	/* Filename prefix for grid files */
  FILE_NAME receptor_file_name;	/* Filename of receptor atoms */

  MOLECULE receptor;		/* Receptor atoms for when grids not used */
  SLINT **atom;			/* 3D grid of receptor atoms */

  int grid_cutoff;              /* Cutoff converted to integer grid units */

} SCORE_GRID;


typedef struct score_bump_struct
{
  int flag;			/* Flag for bump checking */
  int init_flag;		/* Flag for bump check initialization */
  unsigned char *grid;		/* Bump grid */

  float clash_energy;		/* Maximum VDW energy for grid point */
  float clash_overlap;		/* VDW overlap allowed of probe with receptor */
  float distance;		/* Longest distance between bumping atoms */

  int maximum;			/* Maximum bumps allowed for orientation */

} SCORE_BUMP;


typedef struct score_contact_struct
{
  int flag;			/* Flag for contact score use */
  int init_flag;		/* Flag for contact score initialization */
  short int *grid;		/* Contact score grid */
  float distance;         	/* Contact distance between heavy atoms */
  float clash_overlap;		/* VDW overlap allowed of probe with receptor */
  float clash_penalty;		/* Contact score penalty for each atom clash */

} SCORE_CONTACT;


typedef struct score_chemical_struct
{
  int flag;			/* Flag for chemical score use */
  int init_flag;		/* Flag for chemical score initialization */

  float **grid;			/* Chemical score grids */
  float distance;         	/* Longest chemical interaction distance */

} SCORE_CHEMICAL;


typedef struct score_energy_struct
{
  int flag;			/* Flag for energy score use */
  int init_flag;		/* Flag for energy score initialization */
  int decomp_flag;		/* Flag to decompose energy by atom */

  float *avdw, *bvdw, *es;	/* VDW and electrostatic potential grids */
#ifdef BUILD_FOR_RDSOL
  float *dslb, *dslx;         /* Bulk and explicit desolvation grids */
#endif
  float *phi;                 /* Delphi electrostatic potential grid */

  float distance;          	/* Cutoff distance for NB interaction */

  int distance_dielectric;      /* Flag for distance dependent dielectric */
  float dielectric_factor;      /* Factor to multiply dielectric function */

  int repulsive_exponent;	/* Exponent of repulsive LJ-pot'l term */
  int attractive_exponent;	/* Exponent of attractive LJ-pot'l term */
  int atom_model;		/* Flag for all atom or united atom model */

  float scale_electro;		/* Scaling factor for electrostatic term */
  float scale_vdw;		/* Scaling factor for Van der Waals term */

  int vdw_init_flag;		/* Flag for vdw initialization */
  float *vdwA;			/* Array of VDW repulsive terms */
  float *vdwB;			/* Array of VDW repulsive terms */

} SCORE_ENERGY;

typedef struct score_rmsd_struct
{
  int flag;			/* Flag for rmsd score use */
  int init_flag;		/* Flag for rmsd score initialization */

} SCORE_RMSD;


/* //////////////////// Score Routines ////////////////////////// */

void make_receptor_grid
(
  SCORE_GRID            *grid,
  SCORE_ENERGY          *energy,
  LABEL                 *label
);

void free_receptor_grid (SCORE_GRID *grid);

int get_grid_index
(
  SCORE_GRID    *grid,
  XYZ           coord
);

int get_grid_coordinate
(
  SCORE_GRID    *grid,
  XYZ           coord,
  int           grid_coord[3]
);

void calc_inter_score_cont
(
  SCORE_GRID    *grid,
  void          *score,
  float         distance_cutoff,
  void          calc_inter_score
                  (SCORE_GRID *, void *, LABEL *, MOLECULE *,
                  int, int, SCORE_PART *),
  LABEL         *label,
  MOLECULE      *molecule,
  int           atom,
  SCORE_PART    *inter
);

int check_bump
(
  SCORE_GRID *,
  SCORE_BUMP *,
  LABEL *,
  MOLECULE *
);

void calc_inter_contact
(
  SCORE_GRID *,
  SCORE_BUMP *,
  SCORE_CONTACT *,
  LABEL *,
  MOLECULE *,
  int,
  SCORE_PART *
);

void calc_inter_contact_grid
(
  SCORE_GRID *,
  SCORE_BUMP *,
  SCORE_CONTACT *,
  LABEL *,
  MOLECULE *,
  int,
  SCORE_PART *
);

void calc_inter_contact_cont
(
  SCORE_GRID *,
  SCORE_CONTACT *,
  LABEL *,
  MOLECULE *,
  int,
  int,
  SCORE_PART *
);

void calc_pairwise_contact
(
  SCORE_CONTACT *contact,
  LABEL         *label,
  MOLECULE      *origin,
  MOLECULE      *target,
  int           origin_atom,
  int           target_atom,
  float		*score
);

void calc_inter_energy
(
  SCORE_GRID *,
  SCORE_ENERGY *,
  LABEL *,
  MOLECULE *,
  int,
  SCORE_PART *
);

void initialize_vdw_energy
(
  SCORE_ENERGY *energy,
  LABEL_VDW *vdw
);

void free_vdw_energy (SCORE_ENERGY *energy);

void calc_inter_energy_grid
(
  SCORE_GRID *,
  SCORE_ENERGY *,
  LABEL *,
  MOLECULE *,
  int,
  SCORE_PART *
);

void calc_inter_energy_cont
(
  SCORE_GRID *,
  SCORE_ENERGY *,
  LABEL *,
  MOLECULE *,
  int,
  int,
  SCORE_PART *
);

void calc_pairwise_energy
(
  SCORE_ENERGY  *energy,
  LABEL         *label,
  MOLECULE      *origin,
  MOLECULE      *target,
  int           origin_atom,
  int           target_atom,
  float         *vdwA,
  float         *vdwB,
  float         *electro
);

void calc_inter_chemical
(
  SCORE_GRID *,
  SCORE_ENERGY *,
  SCORE_CHEMICAL *,
  LABEL *,
  MOLECULE *,
  int,
  SCORE_PART *
);

void calc_inter_chemical_grid
(
  SCORE_GRID *,
  SCORE_ENERGY *,
  SCORE_CHEMICAL *,
  LABEL *,
  MOLECULE *,
  int,
  SCORE_PART *
);

void calc_inter_chemical_cont
(
  SCORE_GRID *,
  SCORE_ENERGY *,
  LABEL *,
  MOLECULE *,
  int,
  int,
  SCORE_PART *
);

void calc_inter_rmsd
(
  SCORE_GRID *,
  MOLECULE *,
  int,
  SCORE_PART *
);

float calc_rmsd
(
  MOLECULE      *mol_ori,
  MOLECULE      *mol_ref
);

float calc_layer_rmsd
(
  MOLECULE      *mol_ori,
  MOLECULE      *mol_ref
);

float calc_segment_rmsd
(
  MOLECULE      *mol_ori,
  MOLECULE      *mol_ref,
  int		segment
);


void calc_intra_contact
(
  SCORE_CONTACT *contact,
  LABEL		*label,
  MOLECULE	*molecule,
  int		atomi,
  int		atomj,
  SCORE_PART	*intra
);

void calc_intra_energy
(
  SCORE_ENERGY	*energy,
  LABEL		*label,
  MOLECULE	*molecule,
  int		atomi,
  int		atomj,
  SCORE_PART	*intra
);

void calc_intra_chemical
(
  SCORE_ENERGY	*energy,
  LABEL		*label,
  MOLECULE	*molecule,
  int		atomi,
  int		atomj,
  SCORE_PART	*intra
);

void sum_contact
(
  SCORE_PART	*sum,
  SCORE_PART	*increment
);

void sum_energy
(
  SCORE_PART	*sum,
  SCORE_PART	*increment
);

void sum_chemical
(
  SCORE_PART	*sum,
  SCORE_PART	*increment
);

void sum_rmsd
(
  SCORE_PART	*sum,
  SCORE_PART	*increment
);


