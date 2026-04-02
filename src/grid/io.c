/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/
#include "define.h"
#include "utility.h"
#include "mol.h"
#include "global.h"
#include "io.h"
#include "io_mol2.h"
#include "io_pdb.h"
#include "io_ptr.h"
#include "io_sph.h"
#include "transform.h"
#include "rotrans.h"

/* ////////////////////////////////////////////////////////////// */

enum FILE_FORMAT check_file_extension
(
  FILE_NAME file_name,
  int report_error
)
{
  enum FILE_FORMAT format;
  char *extension;

  format = Unknown;

  if (strchr (file_name, '.'))
  {
    extension = strrchr (file_name, '.') + 1;

    if (!strcmp (extension, "mol2"))
      format = Mol2;

    else if (!strcmp (extension, "pdb"))
      format = Pdb;

    else if (!strcmp (extension, "xpdb"))
      format = Xpdb;

    else if (!strcmp (extension, "ptr"))
      format = Ptr;

    else if (!strcmp (extension, "sph"))
      format = Sph;
  }

  if ((format == Unknown) && report_error)
  {
    fprintf (global.outfile,
      "WARNING check_file_extension: \n"
      "        Coordinate files must have a recognized file extension.\n"
      "        Recognized file extensions are: mol2, pdb, xpdb, ptr, sph\n");
  }

  return format;
}


/* ////////////////////////////////////////////////////////////////// */

int read_molecule
(
  MOLECULE *mol_ref,
  MOLECULE *mol_init,
  FILE_NAME in_file_name,
  FILE *in_file,
  int read_source
)
{
  enum FILE_FORMAT format;
  int return_value;

/*
* Determine format of input file
* 2/96 te
*/
  format = check_file_extension (in_file_name, TRUE);

  switch (format)
  {
    case Mol2:
      return_value = 
        read_mol2
        (
          mol_ref,
          in_file_name,
          in_file
        );

      break;

    case Pdb:
    case Xpdb:
      return_value = 
        read_pdb
        (
          mol_ref,
          in_file_name,
          in_file
        );

      break;

    case Ptr:
      return_value = 
        read_ptr
        (
          mol_ref,
          mol_init,
          in_file_name,
          in_file,
          read_source
        );

      break;

    case Sph:
      return_value = 
        read_sph
        (
          mol_ref,
          in_file_name,
          in_file
        );

      break;


    default:
      return FALSE;
  }

/*
* Initialize accessory molecule information
* 10/96 te
*/
  if ((return_value == TRUE) && (format != Ptr) && (mol_ref != NULL))
  {
    atom_neighbors (mol_ref);

    center_of_mass
    (
      mol_ref->coord,
      mol_ref->total.atoms,
      mol_ref->transform.com
    );

    if (mol_init != NULL)
      copy_molecule (mol_init, mol_ref);
  }

  return return_value;
}

/* ////////////////////////////////////////////////////////////////// */

int write_molecule
(
  MOLECULE *molecule,
  FILE_NAME in_file_name,
  FILE_NAME out_file_name,
  FILE *out_file
)
{
  enum FILE_FORMAT format;

/*
* Determine format of output file
* 2/96 te
*/
  format = check_file_extension (out_file_name, TRUE);

  switch (format)
  {
    case Mol2:
      return
        write_mol2
        (
          molecule,
          out_file
        );

    case Pdb:
    case Xpdb:
      return
        write_pdb
        (
          molecule,
          out_file_name,
          out_file
        );

    case Ptr:
      return
        write_ptr
        (
          molecule,
          in_file_name,
          out_file
        );

    case Sph:
      return
        write_sph
        (
          molecule,
          out_file
        );

    default:
      return FALSE;
  }
}


