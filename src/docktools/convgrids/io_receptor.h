/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/

/* Routines used by dock to read receptor data */

int read_receptor
(
  SCORE_ENERGY	*energy,
  LABEL		*label,
  MOLECULE	*molecule,
  FILE_NAME	molecule_file_name,
  FILE		*molecule_file,
  int		need_bonds,
  int		label_chemical,
  int		label_vdw
);


