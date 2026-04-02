/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/

/* Routines used to read and write molecule data */

int read_sph
(
  MOLECULE  *molecule,
  FILE_NAME molecule_file_name,
  FILE      *molecule_file
);

int write_sph
(
  MOLECULE  *molecule,
  FILE      *molecule_file
);

