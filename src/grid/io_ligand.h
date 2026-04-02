/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/

/* Routines used by dock to read and write ligand data */

int get_ligand
(
  DOCK		*dock,
  SCORE		*score,
  LABEL		*label,
  MOLECULE	*lig_orig,
  MOLECULE	*lig_init,
  int		need_bonds,
  int		chem_label,
  int		vdw_label
);

int write_ligand
(
  DOCK		*dock,
  SCORE		*score,
  MOLECULE	*molecule,
  FILE_NAME	molecule_file_name,
  FILE		*molecule_file
);

void write_ligand_info
(
  DOCK		*dock,
  SCORE	        *score,
  MOLECULE	*molecule,
  FILE_NAME	molecule_file_name,
  FILE		*molecule_file
);

