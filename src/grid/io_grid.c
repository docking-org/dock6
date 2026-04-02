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
#include "search.h"
#include "label.h"
#include "score.h"
#include "io_grid.h"
#include "io_pdb.h"
#include "extsymbols.h"

void read_grids
(
  SCORE_GRID		*grid,
  SCORE_BUMP		*bump,
  SCORE_CONTACT		*contact,
  SCORE_CHEMICAL	*chemical,
  SCORE_ENERGY		*energy,
  LABEL_CHEMICAL	*label_chemical
)
{
  STRING100 grid_file_name;
  FILE *grid_file = NULL;
  int i, j, jj, grid_size_check;
  int input_var;
  STRING20 *name = NULL;
  int gridnum;
  int *grid2label = NULL;

  void EXTERNAL_READGRID (char *, int *, float *, int [], float [],
    float [], float [], float [], char []);

  grid->init_flag = TRUE;

/*
* Check if old grids requested
* 3/96 te
*/
  if (grid->version < 3.99)
  {
    emalloc
    (
      (void **) &bump->grid,
      grid->size * sizeof (char),
      "bump array",
      global.outfile
    );

    emalloc
    (
      (void **) &energy->avdw,
      grid->size * sizeof (float),
      "energy grid array",
      global.outfile
    );

    emalloc
    (
      (void **) &energy->bvdw,
      grid->size * sizeof (float),
      "energy grid array",
      global.outfile
    );

    emalloc
    (
      (void **) &energy->es,
      grid->size * sizeof (float),
      "energy grid array",
      global.outfile
    );

    if (strlen (grid->file_prefix) < sizeof (FILE_NAME))
      grid->file_prefix[strlen (grid->file_prefix)] = ' ';

    else
      exit (fprintf (global.outfile,
        "ERROR read_grids: Grid file prefix too long.\n"));

    EXTERNAL_READGRID
    (
      grid->file_prefix,
      &grid->size,
      &grid->spacing,
      grid->span,
      grid->origin,
      energy->avdw,
      energy->bvdw,
      energy->es,
      bump->grid
    );

/*
    energy->atom_model = 'a';
*/
    energy->attractive_exponent = 6;
    energy->repulsive_exponent = 12;

    return;
  }

/*
* Read in bump grid
* 6/95 te
*/
  strcat (strcpy (grid_file_name, grid->file_prefix), ".bmp");
  grid_file = efopen (grid_file_name, "r", global.outfile);
  fprintf
    (global.outfile, "Reading general grid info from %s\n", grid_file_name);

  efread (&grid->size, sizeof (int), 1, grid_file);
  efread (&grid->spacing, sizeof (float), 1, grid_file);
  efread (grid->origin, sizeof (float), 3, grid_file);
  efread (grid->span, sizeof (int), 3, grid_file);

  if (bump->flag || contact->flag)
  {
    fprintf (global.outfile, "Reading bump grid from %s\n", grid_file_name);
    fflush (global.outfile);

    emalloc
    (
      (void **) &bump->grid,
      grid->size * sizeof (char),
      "bump array",
      global.outfile
    );

    efread
    (
      bump->grid,
      sizeof (char),
      grid->size,
      grid_file
    );
  }

  fclose (grid_file);

/*
* Read in contact grid(s)
* 6/95 te
*/
  if (contact->flag)
  {
    strcat (strcpy (grid_file_name, grid->file_prefix), ".cnt");
    fprintf (global.outfile, "Reading contact grid from %s\n", grid_file_name);
    grid_file = efopen (grid_file_name, "r", global.outfile);
    fflush (global.outfile);

    efread (&grid_size_check, sizeof (int), 1, grid_file);
    if ((grid->size) && (grid->size != grid_size_check))
      exit (fprintf (global.outfile,
        "ERROR read_grids: Contact grid incompatible with bump grid.\n"));

    emalloc
    (
      (void **) &contact->grid,
      grid->size * sizeof (short int),
      "contact array",
      global.outfile
    );

    efread
    (
      contact->grid,
      sizeof (short),
      grid->size,
      grid_file
    );

    fclose (grid_file);
  }

/*
* Read in chemical grid(s)
* 12/96 te
*/
  if (chemical->flag)
  {
    strcat (strcpy (grid_file_name, grid->file_prefix), ".chm");
    fprintf
      (global.outfile, "Reading chemical grids from %s\n", grid_file_name);
    grid_file = efopen (grid_file_name, "r", global.outfile);
    fflush (global.outfile);

    efread (&grid_size_check, sizeof (int), 1, grid_file);
    if (grid->size && (grid->size != grid_size_check))
      exit (fprintf (global.outfile,
        "ERROR read_grids: Chemical grids incompatible with bump grid.\n"));

/*
*   Read in number of chemical grids
*   6/95 te
*/
    efread (&gridnum, sizeof (int), 1, grid_file);

    emalloc
    (
      (void **) &name,
      gridnum * sizeof (STRING20),
      "Label names read in from chemical grid",
      global.outfile
    );

    for (i = 0; i < gridnum; i++)
      efread (name[i], sizeof (STRING20), 1, grid_file);

/*
*   Construct conversion of grid number to chemical_label number
*   6/95 te
*/
    ecalloc
    (
      (void **) &grid2label,
      gridnum,
      sizeof (int),
      "grid number to chemical_label number conversion list",
      global.outfile
    );

    for (i = 0; i < gridnum; i++)
    {
      for (j = 0, jj = -1; j < label_chemical->total; j++)
        if (!strcmp (name[i], label_chemical->member[j].name))
          jj = j;

      if (jj == -1)
        exit (fprintf (global.outfile,
          "ERROR read_grids: Label in chemical grid undefined (%s).\n",
          name[i]));

      grid2label[i] = jj;
    }

    ecalloc
    (
      (void **) &chemical->grid,
      label_chemical->total,
      sizeof (float *),
      "chemical grid arrays",
      global.outfile
    );

    for (i = 0; i < gridnum; i++)
    {
      fprintf (global.outfile, "  Reading %s chemical grid.\n",
        label_chemical->member[grid2label[i]].name);
      fflush (global.outfile);

      emalloc
      (
        (void **) &chemical->grid[grid2label[i]],
        grid->size * sizeof (float),
        "chemical array",
        global.outfile
      );

      efread
      (
        chemical->grid[grid2label[i]],
        sizeof (float),
        grid->size,
        grid_file
      );
    }

    fclose (grid_file);

    efree ((void **) &grid2label);
    efree ((void **) &name);
  }

/*
* Read in energy grids
* 6/95 te
*/
  if (energy->flag || chemical->flag)
  {
    strcat (strcpy (grid_file_name, grid->file_prefix), ".nrg");
    fprintf (global.outfile, "Reading energy grids from %s\n", grid_file_name);
    grid_file = efopen (grid_file_name, "r", global.outfile);
    fflush (global.outfile);

    efread (&grid_size_check, sizeof (int), 1, grid_file);
    if ((grid->size) && (grid->size != grid_size_check))
      exit (fprintf (global.outfile,
        "ERROR read_grid: Energy grids incompatible with bump grid.\n"));

    efread (&input_var, sizeof (int), 1, grid_file);

    if ((energy->atom_model != '0') && (energy->atom_model != input_var))
      fprintf (global.outfile,
        "WARNING read_grid: atom_model from grid (%c) \n"
        "                   incompatible with user specified value (%c)\n",
        input_var, energy->atom_model);

    else if ((input_var != 'a') && (input_var != 'u'))
      exit (fprintf (global.outfile,
        "ERROR read_grids: Energy grids are incompatible.\n"));

    else
      energy->atom_model = input_var;

    efread (&input_var, sizeof (int), 1, grid_file);

    if
    (
      (energy->attractive_exponent != 0) &&
      (energy->attractive_exponent != input_var)
    )
      fprintf (global.outfile,
        "WARNING read_grid: attractive_exponent from grid (%d) \n"
        "                   incompatible with user specified value (%d)\n",
        input_var, energy->attractive_exponent);

    else
      energy->attractive_exponent = input_var;

    efread (&input_var, sizeof (int), 1, grid_file);

    if
    (
      (energy->repulsive_exponent != 0) &&
      (energy->repulsive_exponent != input_var)
    )
      fprintf (global.outfile,
        "WARNING read_grid: repulsive_exponent from grid (%d) \n"
        "                   incompatible with user specified value (%d)\n",
        input_var, energy->repulsive_exponent);

    else
      energy->repulsive_exponent = input_var;

    fprintf (global.outfile,
      "  VDW grids use a %1d-%1d Lennard-Jones potential with %s atom model.\n",
      energy->attractive_exponent, energy->repulsive_exponent,
      energy->atom_model == 'a' ? "an all" : "a united");

    emalloc
    (
      (void **) &energy->avdw,
      grid->size * sizeof (float),
      "energy grid array",
      global.outfile
    );

    emalloc
    (
      (void **) &energy->bvdw,
      grid->size * sizeof (float),
      "energy grid array",
      global.outfile
    );

    emalloc
    (
      (void **) &energy->es,
      grid->size * sizeof (float),
      "energy grid array",
      global.outfile
    );

    fprintf (global.outfile, "  Reading attractive VDW energy grid.\n");
    fflush (global.outfile);
    efread
    (
      energy->bvdw,
      sizeof (float),
      grid->size,
      grid_file
    );

    fprintf (global.outfile, "  Reading repulsive VDW energy grid.\n");
    fflush (global.outfile);
    efread
    (
      energy->avdw,
      sizeof (float),
      grid->size,
      grid_file
    );

    fprintf (global.outfile, "  Reading electrostatic energy grid.\n");
    fflush (global.outfile);
    efread
    (
      energy->es,
      sizeof (float),
      grid->size,
      grid_file
    );

    fclose (grid_file);
  }

  else
  {
    energy->atom_model = 'a';
    energy->attractive_exponent = 6;
    energy->repulsive_exponent = 12;
  }
}


/* ////////////////////////////////////////////////////////////// */

void write_grids
(
  SCORE_GRID		*grid,
  SCORE_BUMP		*bump,
  SCORE_CONTACT		*contact,
  SCORE_CHEMICAL	*chemical,
  SCORE_ENERGY		*energy,
  LABEL_CHEMICAL	*label_chemical
)
{
  int i;
  STRING100 grid_file_name;
  FILE *grid_file;

/*
* Write out bump grid
* 6/95 te
*/
  strcat (strcpy (grid_file_name, grid->file_prefix), ".bmp");
  fprintf (global.outfile, "Writing general grid info to %s\n", grid_file_name);
  grid_file = efopen (grid_file_name, "w", global.outfile);

  efwrite (&grid->size, sizeof (int), 1, grid_file);
  efwrite (&grid->spacing, sizeof (float), 1, grid_file);
  efwrite (grid->origin, sizeof (float), 3, grid_file);
  efwrite (grid->span, sizeof (int), 3, grid_file);

  if (bump->flag || contact->flag || chemical->flag)
  {
    fprintf (global.outfile, "Writing bump grid to %s\n", grid_file_name);
    fflush (global.outfile);

    efwrite
    (
      bump->grid,
      sizeof (char),
      grid->size,
      grid_file
    );
  }

  fclose (grid_file);

/*
* Write out contact grid
* 6/95 te
*/
  if (contact->flag)
  {
    strcat (strcpy (grid_file_name, grid->file_prefix), ".cnt");
    fprintf (global.outfile, "Writing contact grid to %s\n", grid_file_name);
    grid_file = efopen (grid_file_name, "w", global.outfile);
    fflush (global.outfile);

    efwrite (&grid->size, sizeof (int), 1, grid_file);
    efwrite
    (
      contact->grid,
      sizeof (short),
      grid->size,
      grid_file
    );

    fclose (grid_file);
  }
    
/*
* Write out chemical grids
* 12/96 te
*/
  if (chemical->flag)
  {
    strcat (strcpy (grid_file_name, grid->file_prefix), ".chm");
    fprintf (global.outfile, "Writing chemical grids to %s\n", grid_file_name);
    grid_file = efopen (grid_file_name, "w", global.outfile);
    fflush (global.outfile);

    efwrite (&grid->size, sizeof (int), 1, grid_file);
    efwrite (&label_chemical->total, sizeof (int), 1, grid_file);

    for (i = 0; i < label_chemical->total; i++)
      efwrite (label_chemical->member[i].name, sizeof (STRING20), 1, grid_file);

    for (i = 0; i < label_chemical->total; i++)
    {
      fprintf (global.outfile, "  Writing %s chemical grid\n",
        label_chemical->member[i].name);
      fflush (global.outfile);

      efwrite
      (
        chemical->grid[i],
        sizeof (float),
        grid->size,
        grid_file
      );
    }

    fclose (grid_file);
  }
    
/*
* Write out energy grids
* 6/95 te
*/
  if (energy->flag)
  {
    strcat (strcpy (grid_file_name, grid->file_prefix), ".nrg");
    fprintf (global.outfile, "Writing energy grids to %s\n", grid_file_name);
    grid_file = efopen (grid_file_name, "w", global.outfile);
    fflush (global.outfile);

    efwrite (&grid->size, sizeof (int), 1, grid_file);
    efwrite (&energy->atom_model, sizeof (int), 1, grid_file);
    efwrite (&energy->attractive_exponent, sizeof (int), 1, grid_file);
    efwrite (&energy->repulsive_exponent, sizeof (int), 1, grid_file);

    fprintf (global.outfile, "  Writing attractive VDW energy grid\n");
    fflush (global.outfile);
    efwrite
    (
      energy->bvdw,
      sizeof (float),
      grid->size,
      grid_file
    );

    fprintf (global.outfile, "  Writing repulsive VDW energy grid\n");
    fflush (global.outfile);
    efwrite
    (
      energy->avdw,
      sizeof (float),
      grid->size,
      grid_file
    );

    fprintf (global.outfile, "  Writing electrostatic energy grid\n");
    fflush (global.outfile);
    efwrite
    (
      energy->es,
      sizeof (float),
      grid->size,
      grid_file
    );

    fclose (grid_file);
  }
}

/* ///////////////////////////////////////////////////////////////// */

int read_box (char *filename,
              XYZ com,
              XYZ dimension)
{
  STRING80 line;
  FILE *file;

/*
* Open file
* 6/95 te
*/
  file = efopen (filename, "rb", global.outfile);

/*
* Read in data
* 6/95 te
*/
  fgets (line, 80, file);
  fgets (line, 80, file);

  if (sscanf (&line[25], "%f %f %f", &com[0], &com[1], &com[2]) != 3)
  {
    fprintf (global.outfile, "Error reading center of mass from box file.\n");
    return FALSE;
  }

  fprintf (global.outfile, "%-40s: %8.3f %8.3f %8.3f\n",
    "Box center of mass", com[0], com[1], com[2]);

  fgets (line, 80, file);

  if (sscanf (&line[29], "%f %f %f",
    &dimension[0], &dimension[1], &dimension[2]) != 3)
  {
    fprintf (global.outfile, "Error reading dimensions from box file.\n");
    return FALSE;
  }

  fprintf (global.outfile, "%-40s: %8.3f %8.3f %8.3f\n",
    "Box dimensions", dimension[0], dimension[1], dimension[2]);

  fclose (file);
  return TRUE;
}



/* ///////////////////////////////////////////////////////////////// */

void free_grids
(
  SCORE_GRID		*grid,
  SCORE_BUMP		*bump,
  SCORE_CONTACT		*contact,
  SCORE_CHEMICAL	*chemical,
  SCORE_ENERGY		*energy,
  LABEL_CHEMICAL	*label_chemical
)
{
  int i;

  if (grid->init_flag)
    free_receptor_grid (grid);

  if (bump->grid)
    efree ((void **) &bump->grid);

  if (contact->flag)
    efree ((void **) &contact->grid);

  if (chemical->flag)
  {
    for (i = 0; i < label_chemical->total; i++)
      efree ((void **) &chemical->grid[i]);

    efree ((void **) &chemical->grid);
  }

  if (energy->flag || chemical->flag)
  {
    efree ((void **) &energy->avdw);
    efree ((void **) &energy->bvdw);
    efree ((void **) &energy->es);
  }
}

