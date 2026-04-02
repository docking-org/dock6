/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
9/96
*/
#include <time.h>
#include "define.h"
#include "utility.h"
#include "mol.h"
#include "global.h"
#include "label.h"
#include "score.h"
#include "io.h"
#include "io_receptor.h"
#include "io_grid.h"
#include "grid.h"
#include "score_grid.h"
#include "parm_grid.h"

GLOBAL global = {0};

int main (int argc, char *argv[])
{
/*
* General variables (see header files for globally defined variables
* 10/95 te
*/
  int i;
  GRID grid = {0};
  SCORE_GRID score_grid = {0};
  SCORE_BUMP score_bump = {0};
  SCORE_CONTACT score_contact = {0};
  SCORE_CHEMICAL score_chemical = {0};
  SCORE_ENERGY score_energy = {0};
  LABEL label = {0};
  MOLECULE receptor;

/*
* Functions used in the main routine
* 10/95 te
*/
  void write_program_header (void);
  void set_memory_limit     (void);

/*
* Process command line arguments
* 10/95 te
*/
  process_commands (argc, argv);

/* 
* Zero the memory
* 2/99 zou
*/
  bzero(&receptor, sizeof(receptor));

/*
* Output program header and user information
* 10/95 te
*/
  write_program_header ();

/*
* Set ceiling on how memory can be dynamically allocated for this run
* 11/95 te
*/
  set_memory_limit ();

/*
* Read in parameters from input file
* 10/95 te
*/
  if (!get_parameters
    (
      &grid,
      &score_grid,
      &score_bump,
      &score_contact,
      &score_chemical,
      &score_energy,
      &label
    ))
  {
    fprintf (stdout, "Error reading control parameters.\n");
    exit (EXIT_FAILURE);
  }

/*
* Read in vdw definitions
* 10/95 te
*/
  if (label.vdw.flag)
    get_vdw_labels (&label.vdw);

/*
* Read in label definitions
* 10/95 te
*/

  if (label.chemical.flag)
    get_chemical_labels (&label.chemical);

/*
* Read in receptor coordinates
* 10/95 te
*/
  fprintf (global.outfile, "\nReading in coordinates of receptor.\n");
  fflush (global.outfile);
  grid.in_file = efopen (grid.in_file_name, "r", global.outfile);

  if (read_receptor
    (&score_energy, &label, &receptor, grid.in_file_name, grid.in_file,
    label.vdw.flag && label.chemical.flag,
    label.chemical.flag, label.vdw.flag)
    != TRUE)
  {
    fprintf (global.outfile, "Error reading in receptor.\n");
    exit (EXIT_FAILURE);
  }

  fclose (grid.in_file);

/*
* Write out receptor coordinates
* 10/95 te
*/
  if (grid.output_molecule)
  {
    fprintf (global.outfile, "\nWriting out processed receptor.\n");
    fflush (global.outfile);
    grid.out_file = efopen (grid.out_file_name, "w", global.outfile);

    write_molecule
    (
      &receptor,
      grid.in_file_name,
      grid.out_file_name,
      grid.out_file
    );

    fclose (grid.out_file);
  }

  if (score_grid.flag)
  {
/*
*   Read in box coordinates which define the boundaries of the grid
*   10/95 te
*/
    fprintf (global.outfile, "\nReading in grid box information.\n");
    if (!read_box
      (grid.box_file_name, grid.box_com, grid.box_dimension))
    {
      fprintf (global.outfile, "Error reading box file %s.\n",
        grid.box_file_name);
      exit (EXIT_FAILURE);
    }

/*
*   Calculate grid variables
*   10/95 te
*/
    for (i = 0, score_grid.size = 1; i < 3; i++)
    {
      score_grid.origin[i] =
        grid.box_com[i] - grid.box_dimension[i] / 2.0;
      score_grid.span[i] =
        (int) (grid.box_dimension[i] / score_grid.spacing + 1.0) + 1;
      score_grid.size *= score_grid.span[i];
    }

    fprintf (global.outfile, "%-40s: %8d %8d %8d\n",
      "Number of grid points per side [x y z]",
      score_grid.span[0], score_grid.span[1], score_grid.span[2]);

    fprintf (global.outfile, "%-40s: %8d\n",
      "Total number of grid points", score_grid.size);
    fflush (global.outfile);

/*
*   Allocate memory for grids
*   10/95 te
*/
    if (score_bump.flag)
    {
      emalloc
      (
        (void **) &score_bump.grid,
        score_grid.size * sizeof (char),
        "chemgrid arrays",
        global.outfile
      );
      memset (score_bump.grid, ~((char) 0), score_grid.size);
    }

    if (score_contact.flag)
      ecalloc
      (
        (void **) &score_contact.grid,
        score_grid.size,
        sizeof (short),
        "chemgrid arrays",
        global.outfile
      );

    if (score_chemical.flag)
    {
      ecalloc
      (
        (void **) &score_chemical.grid,
        label.chemical.total,
        sizeof (float *),
        "chemgrid arrays",
        global.outfile
      );

      for (i = 0; i < label.chemical.total; i++)
        ecalloc
        (
          (void **) &score_chemical.grid[i],
          score_grid.size,
          sizeof (float),
          "chemgrid arrays",
          global.outfile
        );
    }

    if (score_energy.flag)
    {
      ecalloc
      (
        (void **) &score_energy.avdw,
        score_grid.size,
        sizeof (float),
        "chemgrid arrays",
        global.outfile
      );

      ecalloc
      (
        (void **) &score_energy.bvdw,
        score_grid.size,
        sizeof (float),
        "chemgrid arrays",
        global.outfile
      );

      ecalloc
      (
        (void **) &score_energy.es,
        score_grid.size,
        sizeof (float),
        "chemgrid arrays",
        global.outfile
      );
    }

/*
*   Calculate score grids
*   10/95 te
*/
    fprintf (global.outfile, "\nGenerating scoring grids.\n");
    fflush (global.outfile);

    make_grids
    (
      &score_grid,
      &score_bump,
      &score_contact,
      &score_chemical,
      &score_energy,
      &label,
      &receptor
    );

/*
*   Write out grids
*   10/95 te
*/
    write_grids
    (
      &score_grid,
      &score_bump,
      &score_contact,
      &score_chemical,
      &score_energy,
      &label.chemical
    );
  }

  fprintf (global.outfile, "\nFinished calculation.\n");

  return EXIT_SUCCESS;
}


/* /////////////////////////////////////////////////////////////////////// */

void write_program_header ()
{
  extern GLOBAL global;

  if (global.outfile != stdout)
    fprintf (global.outfile, "\n\n"
      "        UUUUUUU    CCCCCCCC    SSSSSSS/  FFFFFFFFF  \n"
      "      UU/    UU/  CC/    CC/     SS/     FF/    FF/ \n"
      "     UU/         CCCCCCCC/      SS/     FF/    FF/  \n"
      "    UU/  UUUU/  CC/ CC\\        SS/     FF/    FF/   \n"
      "   UU/    UU/  CC/   CC\\      SS/     FF/    FF/    \n"
      "   UUUUUUU/   CC/     CC\\  SSSSSS/  FFFFFFFFF/      \n\n\n");

  else
    fprintf (global.outfile, "\n\n"
      "        U[4mUUUU[0mUU    CC[4mCCCC[0mCC    [4mSS[0mSSS[4mSS[0m/  FFF[4mFFFF[0mFF  \n"
      "      UU/    UU/  CC/    CC/     SS/     FF/    FF/ \n"
      "     UU/         CCCCCCCC/      SS/     FF/    FF/  \n"
      "    UU/  [4mU[0mUUU/  CC/ CC\\        SS/     FF/    FF/   \n"
      "   UU/    UU/  CC/   CC\\      SS/     FF/    FF/    \n"
      "   [4mUUUUUUU[0m/   [4mCC[0m/     [4mCC\\[0m  [4mSSSSSS[0m/  [4mFFFFFFFFF[0m/      \n\n\n");

  fprintf (global.outfile,
    "University of California at San Francisco, DOCK %s\n", DOCK_VERSION);

  fflush (global.outfile);
}
