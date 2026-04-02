/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/

/* Structures to store atom labelling definitions */

typedef struct vdw_member_struct
{
  STRING20	name;			/* Member name */
  NODE		*definition;		/* Atom definition */
  int		definition_total;	/* Number of definitions */

  int		atom_model;		/* All-atom or united-atom */
  int		bump_id;		/* Identifier used in bump grid */
  int		heavy_flag;		/* Flag for heavy (non-hydrogen) atom */
  int		valence;		/* Number of allowed bonds */
  float		radius;			/* VDW radius */
  float		well_depth;		/* VDW well-depth */

} VDW_MEMBER;

typedef struct label_vdw_struct
{
  int		flag;			/* Flag for vdw labeling */
  int		init_flag;		/* Flag for vdw initialization */
  VDW_MEMBER	*member;		/* Label members */
  int		total;			/* Number of label members */
  NODE		*definition;		/* Member definitions */
  FILE_NAME	file_name;		/* File containing label parameters */

} LABEL_VDW;


/* Routines to manipulate atom label structures */

void get_vdw_labels (LABEL_VDW *vdw);
void free_vdw_labels (LABEL_VDW *vdw);

int  assign_vdw_labels
(
  LABEL_VDW	*vdw,
  int		atom_model,
  MOLECULE	*molecule
);

int count_heavies
(
  MOLECULE	*molecule,
  int		segment_id
);

