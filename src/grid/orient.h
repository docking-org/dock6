/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/

typedef struct orient_struct
{
  int flag;			/* Flag for orienting molecule */
  int init_flag;		/* Flag for initializing orient variables */
  int random_flag;		/* Flag for random orientation sampling */
  int max;			/* Maximum number of orientations */

  XYZ span;			/* Half width of box */
  XYZ center;			/* Center of box */

  int *orient_list;
  int *target_list;

  MOLECULE orient;
  MOLECULE target;

  MATCH match;

} ORIENT;


/*
Routines defined in transform.c, that are called by outside functions
*/

int get_orientation
(
  DOCK          *dock,
  ORIENT        *orient,
  LABEL         *label,
  MOLECULE      *mol_ref,
  MOLECULE      *mol_conf,
  MOLECULE      *mol_ori,
  int           molecule_id,
  int           conformation_id,
  int           orientation_id
);

void free_orients
(
  LABEL         *label,
  ORIENT        *orient
);

int get_orient
(
  ORIENT        *orient,
  LABEL         *label,
  MOLECULE      *mol_ref,
  MOLECULE      *mol_ori,
  int           molecule,
  int           conformation,
  int           orientation
);

void free_orient
(
  ORIENT        *orient
);

