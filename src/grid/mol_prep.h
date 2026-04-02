/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
2/96
*/

/* Routines used to read and write molecule data */

int prepare_molecule
(
  MOLECULE	*molecule,
  FILE_NAME	molecule_file_name,
  FILE		*molecule_file,
  LABEL		*label,
  int		atom_model,
  int		need_bonds,
  int		assign_chemical,
  int		assign_vdw
);

void deduce_bonds	(LABEL_VDW *, int, MOLECULE *);

