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
#include "io_ligand.h"
#include "transform.h"
#include "rank.h"


/* ////////////////////////////////////////////////////////////// */

int write_info
(
  DOCK *dock,
  SCORE *score,
  LIST *list,
  int ligand_read_num,
  int ligand_dock_num,
  int ligand_skip_num,
  float time
)
{
  int i, j;
  FILE *infofile = NULL;

  infofile = rfopen (dock->info_file_name, "w", global.outfile);

  fprintf (infofile, "%-35s : %10d\n", "Compounds read", ligand_read_num);
  fprintf (infofile, "%-35s : %10d\n", "Compounds docked", ligand_dock_num);
  fprintf (infofile, "%-35s : %10d\n", "Compounds skipped", ligand_skip_num);
  fprintf (infofile, "%-35s : %10.2f\n", "Elapsed CPU time (sec)", time);
  fprintf (infofile, "%-35s : %10.2f\n", "Time per docked compound (sec)",
    time / (float) ligand_dock_num);

  if (dock->rank_ligands)
    for (i = 0; i < SCORE_TOTAL; i++)
    {
      if (score->type[i].flag)
      {
        fprintf (infofile,
          "\n\nCurrent best %s scorers: \n", score->type[i].name);

        for (j = 0; j < list->total[i]; j++)
          fprintf (infofile, "%2d: %7.2f %s %s\n",
            j + 1,
            list->member[i][j]->score.total,
            list->member[i][j]->info.name,
            (list->member[i][j]->transform.refl_flag == 1
              ? "(REFLECTED)" : ""));
      }
    }

  fclose (infofile);
  return TRUE;
}

/* ////////////////////////////////////////////////////////////// */

int write_topscorers
(
  DOCK *dock,
  SCORE *score,
  LIST *list,
  MOLECULE *mol_ref,
  MOLECULE *mol_ori
)
{
  int i, j;
  long file_position;

/*
* Save the current position in the ligand input file
* 7/95 te
*/
  if (!dock->parallel.flag || dock->parallel.server)
    file_position = ftell (dock->ligand_file);

/*
* Loop through all requested scoring types
* 7/95 te
*/
  for (i = 0; i < SCORE_TOTAL; i++)
  {
    if (!score->type[i].flag)
      continue;

/*
*   Open up the output file and reset output numbering
*   6/97 te
*/
    if (dock->rank_ligands)
      score->type[i].file =
        rfopen (score->type[i].file_name, "w", global.outfile);

    for (j = 0; j < list->total[i]; j++)
    {
/*
*     If coordinates are needed for output (output file is NOT ptr format),
*     AND the list doesn't contain coordinates,
*     then reread them and transform them
*     9/97 te
*/
      if ((check_file_extension (score->type[i].file_name, '.') != Ptr) &&
        (list->lite_flag == TRUE))
      {
        reset_molecule (mol_ref);

        fseek
          (dock->ligand_file,
          list->member[i][j]->info.file_position,
          SEEK_SET);

        read_molecule
          (mol_ref, NULL, dock->ligand_file_name, dock->ligand_file, TRUE);

        if (strcmp
          (list->member[i][j]->info.name, mol_ref->info.name) != 0)
        {
          fprintf (global.outfile,
            "ERROR write_topscorers: "
            "Read incorrect ligand from input file:\n");
          fprintf (global.outfile, "      Intended : %s\n",
            list->member[i][j]->info.name);
          fprintf (global.outfile, "      Actual   : %s\n",
            mol_ref->info.name);
          exit (EXIT_FAILURE);
        }

        copy_molecule (mol_ori, mol_ref);
        copy_member (list->lite_flag, mol_ori, list->member[i][j]);

        if ((mol_ori->total.torsions > 0) && (mol_ori->transform.tors_flag))
        {
          revise_atom_neighbors (mol_ori);
          get_torsion_neighbors (mol_ori);
        }

        if (mol_ori->transform.flag)
          transform_molecule (mol_ori, mol_ref);

        write_ligand
        (
          dock,
          score,
          mol_ori,
          score->type[i].file_name,
          score->type[i].file
        );
      }

      else
        write_ligand
        (
          dock,
          score,
          list->member[i][j],
          score->type[i].file_name,
          score->type[i].file
        );
    }

    if (dock->rank_ligands)
      fclose (score->type[i].file);
  }

  if (!dock->parallel.flag || dock->parallel.server)
    fseek (dock->ligand_file, file_position, SEEK_SET);

  return TRUE;
}

/* ////////////////////////////////////////////////////////////////// */

int write_restartinfo
(
  DOCK *dock,
  SCORE *score,
  LIST *list,
  int ligand_read_num,
  int ligand_dock_num,
  int ligand_skip_num,
  float clock_elapsed
)
{
  long file_position;
  FILE *restart_file;

  void efwrite (void *, size_t, size_t, FILE *);

  file_position = ftell (dock->ligand_file);

  restart_file = efopen (dock->restart_file_name, "w", global.outfile);

  efwrite (&file_position, sizeof (long), 1, restart_file);
  efwrite (&ligand_read_num, sizeof (int), 1, restart_file);
  efwrite (&ligand_dock_num, sizeof (int), 1, restart_file);
  efwrite (&ligand_skip_num, sizeof (int), 1, restart_file);
  efwrite (&clock_elapsed, sizeof (float), 1, restart_file);

  save_lists (score, list, restart_file);

  fclose (restart_file);
  return TRUE;
}

/* ////////////////////////////////////////////////////////////////// */

int read_restartinfo
(
  DOCK *dock,
  SCORE *score,
  LIST *list,
  int *ligand_read_num,
  int *ligand_dock_num,
  int *ligand_skip_num,
  float *clock_elapsed
)
{
  long file_position;
  FILE *restart_file;

  void efread (void *, size_t, size_t, FILE *);

  restart_file = efopen (dock->restart_file_name, "r", global.outfile);

  efread (&file_position, sizeof (long), 1, restart_file);
  efread (ligand_read_num, sizeof (int), 1, restart_file);
  efread (ligand_dock_num, sizeof (int), 1, restart_file);
  efread (ligand_skip_num, sizeof (int), 1, restart_file);
  efread (clock_elapsed, sizeof (float), 1, restart_file);

  if (!load_lists (score, list, restart_file))
    exit (fprintf (global.outfile,
      "ERROR read_restartinfo: Error reading restart information.\n"));

  fclose (restart_file);
  fseek (dock->ligand_file, file_position, SEEK_SET);
  return TRUE;
}

