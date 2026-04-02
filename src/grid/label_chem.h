/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/

/* Structures to store atom labelling definitions */

typedef struct chemical_member_struct
{
  STRING20	name;			/* Member name */
  NODE		*definition;		/* Member definitions */
  int		definition_total;	/* Number of definitions */

  float		radius;			/* Interaction radius */
  float		tolerance;		/* Interaction tolerance */

} CHEMICAL_MEMBER;

typedef struct chemical_screen_struct
{
  int		flag;			/* Flag for chemical screen */
  int		process_flag;		/* Flag to screen ligands */
  int		construct_flag;		/* Flag for keying database */
  int		pharmaco_flag;		/* Flag for pharmacophore search */
  int		similar_flag;		/* Flag for similarity search */
  float		dissimilar_maximum;	/* Similarity cutoff */
  int		fold_flag;		/* Flag for folding label keys */
  float		distance_minimum;	/* Minimum distance of interest */
  float		distance_maximum;	/* Maximum distance of interest */
  float		distance_interval;	/* Distance interval between keys */
  int		interval_total;		/* Number of distance keys */

} CHEMICAL_SCREEN;

typedef struct label_chemical_struct
{
  int		flag;			/* Flag for chemical labeling */
  int		init_flag;		/* Flag for chem label initialization */
  CHEMICAL_MEMBER *member;		/* Label members */
  int		total;			/* Number of members */
  NODE		*definition;		/* Label definitions */

  CHEMICAL_SCREEN screen;		/* Chemical screen parameters */

  float		**match_table;		/* Match compatibility table */
  float		**score_table;		/* Score weight table */
  float		**screen_table;		/* Screen weight table */

  FILE_NAME	file_name;		/* File containing label data */
  FILE_NAME	match_file_name;	/* File with match table */
  FILE_NAME	score_file_name;	/* File with score table */
  FILE_NAME	screen_file_name;	/* File with screen table */

} LABEL_CHEMICAL;

/* Routines to manipulate atom label structures */

void get_chemical_labels	(LABEL_CHEMICAL *);
void free_chemical_labels	(LABEL_CHEMICAL *);
int  assign_chemical_labels	(LABEL_CHEMICAL *, MOLECULE *);
void read_chemical_labels	(LABEL_CHEMICAL *, MOLECULE *,
				FILE_NAME, FILE *);
void get_table			(LABEL_CHEMICAL *, FILE_NAME, float ***);
void free_table			(LABEL_CHEMICAL *, float ***);

