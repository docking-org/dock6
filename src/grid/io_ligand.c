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
#include "dock.h"
#include "search.h"
#include "label.h"
#include "score.h"
#include "score_dock.h"
#include "io.h"
#include "mol_prep.h"
#include "io_ligand.h"
#include "transform.h"
#include "rotrans.h"

/* ////////////////////////////////////////////////////////////// */

int get_ligand
(
  DOCK *dock,
  SCORE *score,
  LABEL *label,
  MOLECULE *lig_ref,
  MOLECULE *lig_init,
  int need_bonds,
  int chemical_label,
  int vdw_label
)
{
  int return_value;
  FILE_NAME ready_file_name;
  FILE *ready_file, *quit_file;

/*
* Perform parallel client operations, if necessary
* 10/95 te
*/
  if (dock->parallel.flag && !dock->parallel.server)
  {
/*
*   If the ligand file has already been open, then close it now
*   11/96 te
*/
    if (dock->ligand_file != NULL)
    {
      fclose (dock->ligand_file);
      remove (dock->ligand_file_name);
    }

/*
*   Announce that we are ready to accept a job
*   10/95 te
*/
    sprintf (dock->ligand_file_name, "%s_%s.ptr",
      dock->parallel.server_name,
      dock->parallel.client_name[0]);
    sprintf (ready_file_name, "%s.ready", dock->parallel.client_name[0]);
    ready_file = rfopen (ready_file_name, "w", global.outfile);
    fclose (ready_file);
  
/*
*   Wait until the ready file has been removed by the server
*   10/95 te
*/
    while (ready_file = fopen (ready_file_name, "r"))
    {
      fclose (ready_file);

      if (quit_file = fopen (dock->quit_file_name, "r"))
      {
        fclose (quit_file);
        remove (ready_file_name);
        remove (dock->ligand_file_name);
        return EOF;
      }
    }

    dock->ligand_file = rfopen (dock->ligand_file_name, "r", global.outfile);
  }

/*
* Reset molecule data structure
* 2/96 te
*/
  reset_molecule (lig_ref);
  lig_ref->info.file_position = ftell (dock->ligand_file);

/*
* Read in coordinates of molecule
* 6/95 te
*/
  return_value =
    read_molecule
    (
      lig_ref,
      lig_init,
      dock->ligand_file_name,
      dock->ligand_file,
      !label->chemical.screen.process_flag
    );

  if (return_value != TRUE)
    return return_value;

/*
* Prepare molecule for docking
* 6/95 te
*/
  if (!label->chemical.screen.process_flag &&
    ((return_value =
      prepare_molecule
      (
        lig_init,
        dock->ligand_file_name,
        dock->ligand_file,
        label,
        score->energy.atom_model,
        need_bonds,
        chemical_label,
        vdw_label
      )) != TRUE))
    return return_value;

/*
* Check number of heavy atoms
* 6/96 te
*/
  if (dock->multiple_ligands && vdw_label &&
    ((lig_init->transform.heavy_total > dock->max_heavies) ||
    (lig_init->transform.heavy_total < dock->min_heavies)))
  {
    if (global.output_volume != 't')
      fprintf (global.outfile, "Skipped %s (%d heavy atoms).\n",
        lig_init->info.name, lig_init->transform.heavy_total);

    return NULL;
  }

  return TRUE;
}


/* ////////////////////////////////////////////////////////////////// */

int write_ligand
(
  DOCK *dock,
  SCORE *score,
  MOLECULE *molecule,
  FILE_NAME input_molecule_file_name,
  FILE *molecule_file
)
{
  int i;
  char *molecule_file_name;
  FILE_NAME ready_file_name, delegate_file_name;
  FILE *ready_file, *quit_file;

/*
* Perform parallel server operations
* 10/95 te
*/
  if (dock->parallel.server)
  {
/*
*   Look for a client that is ready to process a molecule
*   10/95 te
*/
    for (molecule_file_name = NULL; !molecule_file_name;)
    {
      for (i = 0; i < dock->parallel.client_total; i++)
      {
        sprintf (ready_file_name, "%s.ready", dock->parallel.client_name[i]);
  
        if (ready_file = fopen (ready_file_name, "r"))
        {
          fclose (ready_file);
          sprintf
            (delegate_file_name, "%s_%s.ptr",
            dock->parallel.server_name, dock->parallel.client_name[i]);
          molecule_file_name = delegate_file_name;

          molecule_file = rfopen (molecule_file_name, "w", global.outfile);
          break;
        }
      }

      if (quit_file = fopen (dock->quit_file_name, "r"))
      {
        fclose (quit_file);
        return EOF;
      }
    }
  }

  else
    molecule_file_name = input_molecule_file_name;

/*
* Write out coordinates of molecule
* 7/95 te
*/
  write_ligand_info
  (
    dock,
    score,
    molecule,
    molecule_file_name,
    molecule_file
  );

/*
* Perform server duties
* 10/95 te
*/
  if (dock->parallel.server)
  {
    efclose (&molecule_file);
    rremove (ready_file_name, global.outfile);
  }

  return TRUE;
}


/* ////////////////////////////////////////////////////////////////// */

void write_ligand_info
(
  DOCK *dock,
  SCORE *score,
  MOLECULE *molecule,
  FILE_NAME molecule_file_name,
  FILE *molecule_file
)
{
  enum FILE_FORMAT format;
  char *flag = NULL;
  void write_header
    (DOCK *, SCORE *, MOLECULE *, FILE *, char *);
  void write_trailer
    (DOCK *, SCORE *, MOLECULE *, FILE *);

  molecule->info.output_id =
    ++(score->type[molecule->score.type].number_written);

/*
* Determine format of output file
* 2/96 te
*/
  format = check_file_extension (molecule_file_name, TRUE);

/*
* Assign the proper flag to identify comment info in coordinate file
* 2/96 te
*/
  switch (format)
  {
    case Mol2:
      flag = "##########";
      break;

    case Pdb:
    case Xpdb:
      flag = "REMARK    ";
      break;

    case Ptr:
    case Sph:
      flag = "";
      break;
  }

/*
* Write out comment information in proper format
* 2/96 te
*/
  switch (format)
  {
    case Mol2:
    case Pdb:
    case Xpdb:
    {
/*
*     Score information is put in header before coordinates
*     2/96 te
*/
      write_header
      (
        dock,
        score,
        molecule,
        molecule_file,
        flag
      );

      write_molecule
      (
        molecule,
        dock->ligand_file_name,
        molecule_file_name,
        molecule_file
      );

      break;
    }

    case Ptr:
    {
/*
*     Score information is put in trailer on same line as pointer information
*     2/96 te
*/
      write_molecule
      (
        molecule,
        dock->ligand_file_name,
        molecule_file_name,
        molecule_file
      );

      write_trailer
      (
        dock,
        score,
        molecule,
        molecule_file
      );

      break;
    }

    case Sph:
    {
/*
*     No score header information is written
*     2/96 te
*/
      write_molecule
      (
        molecule,
        dock->ligand_file_name,
        molecule_file_name,
        molecule_file
      );

      break;
    }
  }
}



/* ////////////////////////////////////////////////////////////////// */

void write_header
(
  DOCK *dock,
  SCORE *score,
  MOLECULE *molecule,
  FILE *molecule_file,
  char *flag
)
{
  fprintf (molecule_file, "%s %-11s : %d\n", flag, "Number",
    molecule->info.output_id);

  if (molecule->info.input_id != NEITHER)
    fprintf (molecule_file, "%s %-11s : %d\n", flag, "Source num",
      molecule->info.input_id);

  fprintf
    (molecule_file, "%s %-11s : %s\n", flag, "Name", molecule->info.name);
  fprintf
    (molecule_file, "%s %-11s : %s\n", flag, "Description",
    molecule->info.comment);

  if (molecule->transform.refl_flag != NEITHER)
    fprintf (molecule_file, "%s %-11s : %d\n", flag, "Reflect",
      molecule->transform.refl_flag);

  if (score->bump.flag)
    fprintf (molecule_file, "%s %-40s : %10d\n", flag, "Number of bumps",
      molecule->score.bumpcount);

  if (molecule->score.total != INITIAL_SCORE)
  {
    if (molecule->score.type == NONE)
      fprintf (molecule_file, "%s %-40s : %10.2f\n", flag, "Score",
        molecule->score.total);

    else if (molecule->score.type == CONTACT)
      fprintf (molecule_file, "%s %-40s : %10.2f\n", flag, "Contact score",
        molecule->score.total);

    else if (molecule->score.type == CHEMICAL)
    {
      fprintf (molecule_file, "%s %-40s : %10.2f\n", flag, "Chemical score",
        molecule->score.total);

      if (score->inter_flag)
      {
        fprintf (molecule_file, "%s %-40s   %10.2f\n", flag,
          "  intermolecular van der Waals", molecule->score.inter.vdw);

        fprintf (molecule_file, "%s %-40s   %10.2f\n", flag,
          "  intermolecular electrostatic", molecule->score.inter.electro);
      }

      if (score->intra_flag)
      {
        fprintf (molecule_file, "%s %-40s   %10.2f\n", flag,
          "  intramolecular van der Waals", molecule->score.intra.vdw);

        fprintf (molecule_file, "%s %-40s   %10.2f\n", flag,
          "  intramolecular electrostatic", molecule->score.intra.electro);
      }
    }
    
    else if (molecule->score.type == ENERGY)
    {
      fprintf (molecule_file, "%s %-40s : %10.2f\n", flag, "Energy score",
        molecule->score.total);

      if (score->inter_flag)
      {
        fprintf (molecule_file, "%s %-40s   %10.2f\n", flag,
          "  intermolecular van der Waals", molecule->score.inter.vdw);

        fprintf (molecule_file, "%s %-40s   %10.2f\n", flag,
          "  intermolecular electrostatic", molecule->score.inter.electro);
      }

      if (score->intra_flag)
      {
        fprintf (molecule_file, "%s %-40s   %10.2f\n", flag,
          "  intramolecular van der Waals", molecule->score.intra.vdw);

        fprintf (molecule_file, "%s %-40s   %10.2f\n", flag,
          "  intramolecular electrostatic", molecule->score.intra.electro);
      }
    }

    else if (molecule->score.type == RMSD)
      fprintf (molecule_file, "%s %-40s : %10.2f\n", flag, "RMSD score",
        molecule->score.total);
  }

  if (dock->multiple_conforms ||
    dock->multiple_orients ||
    score->type[molecule->score.type].minimize)
    fprintf (molecule_file, "%s %-40s : %10.2f\n",
      flag, "RMSD from input orientation (A)", molecule->transform.rmsd);
}



/* ////////////////////////////////////////////////////////////////// */

void write_trailer
(
  DOCK *dock,
  SCORE *score,
  MOLECULE *molecule,
  FILE *molecule_file
)
{
/*
* Move back one byte to overwrite the previous new line character
* 2/96 te
*/
  fseek (molecule_file, -7, SEEK_END);

/*
* Output score information
* 6/95 te
*/
  if (score->bump.flag)
    fprintf (molecule_file, " <BMP> %d", molecule->score.bumpcount);

  if (molecule->score.total != INITIAL_SCORE)
  {
    if (molecule->score.type == NONE)
      fprintf (molecule_file, " <SCORE> %.2f", molecule->score.total);

    else if (molecule->score.type == CONTACT)
      fprintf (molecule_file, " <CNT> %.2f", molecule->score.total);

    else if (molecule->score.type == CHEMICAL)
    {
      fprintf (molecule_file, " <CHM> %.2f", molecule->score.total);

      if (score->inter_flag)
      {
        fprintf (molecule_file, " <INTER_VDW> %.2f",
          molecule->score.inter.vdw);
        fprintf (molecule_file, " <INTER_ELE> %.2f",
          molecule->score.inter.electro);
      }

      if (score->intra_flag)
      {
        fprintf (molecule_file, " <INTRA_VDW> %.2f",
          molecule->score.intra.vdw);
        fprintf (molecule_file, " <INTRA_ELE> %.2f",
          molecule->score.intra.electro);
      }
    }

    else if (molecule->score.type == ENERGY)
    {
      fprintf (molecule_file, " <NRG> %.2f", molecule->score.total);

      if (score->inter_flag)
      {
        fprintf (molecule_file, " <INTER_VDW> %.2f",
          molecule->score.inter.vdw);
        fprintf (molecule_file, " <INTER_ELE> %.2f",
          molecule->score.inter.electro);
      }

      if (score->intra_flag)
      {
        fprintf (molecule_file, " <INTRA_VDW> %.2f",
          molecule->score.intra.vdw);
        fprintf (molecule_file, " <INTRA_ELE> %.2f",
          molecule->score.intra.electro);
      }
    }

    else if (molecule->score.type == RMSD)
      fprintf (molecule_file, " <RMS> %.2f", molecule->score.total);
  }

  if (dock->multiple_orients || score->type[molecule->score.type].minimize)
    fprintf (molecule_file, " <RMSD> %.2f", molecule->transform.rmsd);

  fprintf (molecule_file, " <END>\n");
}


