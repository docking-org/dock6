/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/96
*/

/* Structures to store atom labelling definitions */

typedef struct flex_member_struct
{
  STRING20	name;			/* Member name */
  NODE		*definition;		/* Member definitions */
  int		definition_total;	/* Number of definitions */

  int		drive_id;		/* Torsion search label */
  int		torsion_total;		/* Number of torsion positions */
  float		*torsion;		/* Torsion values */
  int		minimize_flag;		/* Flag for rotation during minim'n */
  float		penalty;		/* Penalty value (unused) */

} FLEX_MEMBER;

typedef struct label_flex_struct
{
  int		flag;			/* Flag for flex labeling */
  int		init_flag;		/* Flag for flex label initialization */

  int		drive_flag;		/* Flag for torsion driver search */
  float		clash_overlap;		/* Clash vdw overlap threshold */
  int		max_conforms;		/* Maximum conformations per layer */
  float		diversity_factor;	/* Scaling of diversity when pruning */

  int		minimize_flag;		/* Flag for torsional minimization */

  int		anchor_flag;		/* Flag for anchor search */
  int		multiple_anchors;	/* Flag for multiple anchors */
  int		anchor_size;		/* Mininum size of anchor segment */
  int		write_flag;		/* Flag to write partial structures */

  int		periph_flag;		/* Flag for peripheral search */
  float		min_rmsd;		/* Minimum rmsd allowed for search */
  int		minimize_anchor_flag;	/* Flag to minimize anchor */
  int		reminimize_layers;	/* Previous layers to reminimize */
  int		reminimize_anchor_flag;	/* Flag to reminimize anchor */
  int		reminimize_ligand_flag;	/* Flag to reminimize ligand */

  int		max_torsions;		/* Maximum torsions per molecule */

  FLEX_MEMBER	*member;		/* Label members */
  int		total;			/* Number of label members */
  NODE		*definition;		/* Member definitions */
  FILE_NAME	file_name;		/* File containing flex data */
  FILE_NAME	search_file_name;	/* File containing search data */

} LABEL_FLEX;

/* Routines to manipulate atom label structures */

int get_flex_labels (LABEL_FLEX *flex);

void assign_flex_labels
(
  LABEL_FLEX	*label_flex,
  MOLECULE	*molecule
);

void free_flex_labels (LABEL_FLEX *label_flex);

int get_flex_search (LABEL_FLEX *);

int check_peripheral_torsion
(
  MOLECULE      *molecule,
  int           torsion_id
);

int get_anchor
(
  LABEL_FLEX	*label_flex,
  MOLECULE	*mol_init,
  MOLECULE	*mol_anch,
  int		anchor
);

void get_segments (LABEL_FLEX *label_flex, MOLECULE *mol_init);

int get_anchor_segment
(
  LABEL_FLEX	*label_flex,
  MOLECULE	*molecule,
  int		anchor_count
);

void get_layers
(
  LABEL_FLEX	*label_flex,
  MOLECULE	*mol_anch
);

void get_single_anchor
(
  LABEL_FLEX	*label_flex,
  MOLECULE	*mol_anch
);

void initialize_segments
(
  LABEL_FLEX	*label_flex,
  MOLECULE	*mol_anch
);

void initialize_anchor_layer
(
  LABEL_FLEX    *label_flex,
  MOLECULE      *molecule
);

void initialize_layer
(
  LABEL_FLEX    *label_flex,
  MOLECULE      *molecule,
  int		layer
);

