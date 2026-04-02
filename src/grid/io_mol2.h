/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/

/* Routines used to read and write molecule data */

int read_mol2
(
  MOLECULE  *molecule,
  FILE_NAME molecule_file_name,
  FILE      *molecule_file
);

int write_mol2
(
  MOLECULE  *molecule,
  FILE      *molecule_file
);

