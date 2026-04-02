/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/

/* Routines used to read and write molecule data */

int read_ptr
(
  MOLECULE  *mol_ref,
  MOLECULE  *mol_init,
  FILE_NAME in_file_name,
  FILE      *in_file,
  int       read_source
);

int write_ptr
(
  MOLECULE  *molecule,
  FILE_NAME in_file_name,
  FILE      *out_file
);

