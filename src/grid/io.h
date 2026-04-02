/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/

/* Routines used to read and write molecule data */

enum FILE_FORMAT {Unknown, Mol2, Pdb, Xpdb, Ptr, Sph};

enum FILE_FORMAT check_file_extension (FILE_NAME file_name, int report_error);

int read_molecule
(
  MOLECULE  *mol_orig,
  MOLECULE  *mol_init,
  FILE_NAME in_file_name,
  FILE      *in_file,
  int       read_source
);

int write_molecule
(
  MOLECULE  *molecule,
  FILE_NAME in_file_name,
  FILE_NAME out_file_name,
  FILE      *out_file
);

