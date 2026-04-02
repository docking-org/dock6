/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
9/96
*/
#define NONE            0
#define CONTACT         1
#define CHEMICAL        2
#define ENERGY          3
#define RMSD            4
#define SCORE_TOTAL     5


typedef struct score_type_struct
{
  int		flag;			/* Flag for this type of scoring */
  STRING20	name;			/* Name of scoring */
  STRING5	abbrev;			/* Abbreviated name of scoring */
  float		maximum;		/* Maximum score filter */
  float		size_penalty;		/* Size penalty to compare molecules */
  int		minimize;		/* Flag for if minimization performed */
  float		convergence;		/* Convergence criteria for minimizer */
  float		termination;		/* Termination criteria for minimizer */
  FILE_NAME	file_name;		/* Output file name */
  FILE		*file;			/* Output file stream */
  int		number_written;		/* Number of orientations written */

} SCORE_TYPE;

typedef struct minimize_struct
{
  int		flag;			/* Flag for minimization */
  int		torsion_flag;		/* Flag for torsion minimization */

  int		iteration;		/* Maximum simplex iterations */
  int		cycle;			/* Maximum simplex cycles */
  float		cycle_converge;		/* Simplex cycle convergence */
  float		rotation;		/* Initial rotation step size */
  float		translation;		/* Initial translation step size */
  float		torsion;		/* Initial torsion step size */

  int		call_total;		/* Total calls */
  int		call_sub_total;		/* Total calls per molecule */
  int		call_min;		/* Minimum calls per molecule */
  int		call_max;		/* Maximum calls per molecule */
  int		vertex_total;		/* Total vertices */
  int		vertex_min;		/* Minimum vertices per call */
  int		vertex_max;		/* Maximum vertices per call */
  int		cycle_total;		/* Total cycles */
  int		cycle_min;		/* Minimum cycles per call */
  int		cycle_max;		/* Maximum cycles per call */
  int		iteration_total;	/* Total iterations */
  int		iteration_min;		/* Minimum iterations per call */
  int		iteration_max;		/* Maximum iterations per call */
  float		delta_total;		/* Total score change */
  float		delta_min;		/* Minimum score change per call */
  float		delta_max;		/* Maximum score change per call */

} MINIMIZE;

typedef struct near_struct
{
  int		**flag;			/* Array of nearby atoms */
  int		total;			/* Total flags */
  char		*name;			/* Name of molecule described */

} NEAR;

typedef struct score_struct
{                              
  int		flag;			/* Flag to use ANY type of scoring */
  int		intra_flag;		/* Flag to get intra-ligand energy */
  int		inter_flag;		/* Flag to get ligand-rec score */

  SCORE_GRID	grid;			/* Grid scoring related info */
  SCORE_BUMP	bump;			/* Bump checking info */
  SCORE_CONTACT	contact;		/* Contact scoring info */
  SCORE_CHEMICAL chemical;		/* Chemical scoring info */
  SCORE_ENERGY	energy;			/* Energy scoring info */
  SCORE_RMSD	rmsd;			/* RMSD scoring info */

  SCORE_TYPE	type[SCORE_TOTAL];	/* Info for each type of scoring */
  int		key[SCORE_TOTAL];	/* Key to scoring types selected */

  MINIMIZE	minimize;		/* Minimization info */

  NEAR		near;			/* Atom proximity info */
  float		rmsd_override;		/* RMS cutoff to force write */
  float		time;			/* Time spent getting energy of pair */

} SCORE;

typedef struct list_struct
{
  int lite_flag;			/* Flag to not store coordinates */
  int total[SCORE_TOTAL];		/* Length of each list */
  int max[SCORE_TOTAL];			/* Maximum length of each list */
  MOLECULE **member[SCORE_TOTAL];	/* Sorted lists of molecules */

} LIST;


void score_setup	(SCORE *, LABEL *, MOLECULE *);

int get_anchor_score
(
  DOCK		*dock,
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*mol_ref,
  MOLECULE	*mol_init,
  MOLECULE	*mol_ori,
  MOLECULE	*mol_min,
  LIST		*best_anchors,
  int		orient_id,
  int		write_flag
);

void get_peripheral_score
(
  DOCK		*dock,
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*mol_ref,
  MOLECULE	*mol_init,
  MOLECULE	*mol_conf,
  MOLECULE	*mol_ori,
  MOLECULE	*mol_min,
  LIST		*best_anchors,
  LIST		*best_orient,
  int		write_flag,
  int		anchor
);

void write_periph_structures
(
  DOCK		*dock,
  SCORE		*score,
  LIST		*list,
  MOLECULE	*mol_ref,
  MOLECULE	*mol_ori,
  int		type,
  int		anchor,
  int		layer,
  int		segment,
  int		shrink_flag
);

void minimize_ligand
(
  DOCK          *dock,
  LABEL         *label,
  SCORE         *score,
  MOLECULE      *mol_ref,
  MOLECULE      *mol_ori,
  MOLECULE      *mol_score,
  int		rigid_flag,
  int		layer_initial,
  int		layer_final
);

float simplex_score (void *simplex, float *vertex);

float calc_score
(
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*molecule,
  int		layer_final
);

void calc_inter_score
(
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*molecule,
  int		atom,
  SCORE_PART	*inter
);

void calc_intra_score
(
  LABEL         *label,
  SCORE         *score,
  MOLECULE      *molecule,
  int           atomi,
  int           atomj,
  SCORE_PART    *intra
);

void sum_score
(
  int           score_type,
  SCORE_PART    *sum,
  SCORE_PART    *increment
);

void output_score_info
(
  DOCK		*dock,
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*molecule,
  LIST		*list,
  int		number_docked,
  float		time
);

void  allocate_lists	(SCORE *, LIST *, int, int);
void  allocate_list	(LIST *, int);
void  reset_lists	(SCORE *, LIST *);
void  reset_list	(LIST *, int);
void  free_lists	(SCORE *, LIST *);
void  free_list		(LIST *, int);
void  reallocate_lists	(SCORE *, LIST *);
void  reallocate_list	(LIST *, int);
void  copy_lists	(SCORE *, LIST *, LIST *);
void  copy_list		(LIST *, LIST *, int);
void  inter_lists	(SCORE *, LIST *, MOLECULE *);
void  merge_lists	(SCORE *, LIST *, LIST *);
void  update_list	(LIST *, MOLECULE *);
int   compare_score	(MOLECULE *, MOLECULE *);
void  sort_list		(LIST *, int type);
void  copy_member	(int, MOLECULE *, MOLECULE *);
int   compare_member	(MOLECULE **, MOLECULE **);
int   print_lists	(SCORE *, LIST *, FILE *);
int   print_list	(SCORE *, LIST *, int, FILE *);
int   save_lists	(SCORE *, LIST *, FILE *);
int   load_lists	(SCORE *, LIST *, FILE *);

void shrink_list
(
  LIST          *full,
  LIST          *shrunk,
  int           type,
  float		rmsd_cutoff
);

void update_shrunk_list
(
  LIST          *shrunk,
  MOLECULE	*molecule
);

void initialize_near	(NEAR *near, MOLECULE *molecule);
void allocate_near	(NEAR *near);
void free_near		(NEAR *near);

void free_scores	(LABEL *label, SCORE *score);
