/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
12/96
*/


int get_anchor_conformation
(
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*mol_init,
  MOLECULE	*mol_anch,
  int		conformer
);

int get_layer_conformation
(
  LABEL		*label,
  SCORE		*score,
  MOLECULE      *molecule,
  int           layer,
  int		conformer
);

int get_segment_conformation
(
  LABEL		*label,
  SCORE		*score,
  MOLECULE      *molecule,
  int           layer,
  int           segment
);

int update_segment_torsions
(
  LABEL_FLEX    *label_flex,
  MOLECULE      *molecule,
  int           segment,
  int           conformer
);

int check_segment_conformation
(
  LABEL		*label,
  SCORE		*score,
  MOLECULE      *molecule,
  int           layer,
  int           segment
);

