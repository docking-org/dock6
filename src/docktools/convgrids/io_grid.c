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

#ifdef BUILD_FOR_RDSOL
  void readgrid_ (char *, int *, float *, int [], float [],
    float [], float [], float [], float [], float [], unsigned char [],
    int *, float *, int [], float [], float []);
#else 
  void readgrid_ (char *, int *, float *, int [], float [],
    float [], float [], float [], unsigned char [], int *, float *, 
    int [], float [], float []);
#endif

  grid->init_flag = TRUE;

/*
* Check if old grids requested
* 3/96 te
/
  if (grid->version < 3.99)
  {
*/
    emalloc
    (
      (void **) &bump->grid,
      grid->size * sizeof (unsigned char),
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

#ifdef BUILD_FOR_RDSOL
    emalloc
    (
      (void **) &energy->dslb,
      grid->size * sizeof (float),
      "energy grid array",
      global.outfile
    );

    emalloc
    (
      (void **) &energy->dslx,
      grid->size * sizeof (float),
      "energy grid array",
      global.outfile
    );
#endif

    emalloc
    (
      (void **) &energy->phi,
	grid->size * sizeof (float),
	"Delphi grid array",
	global.outfile
    );

    if (strlen (grid->file_prefix) < sizeof (FILE_NAME))
      grid->file_prefix[strlen (grid->file_prefix)] = ' ';

    else
      exit (fprintf (global.outfile,
        "ERROR read_grids: Grid file prefix too long.\n"));

    readgrid_
    (
      grid->file_prefix,
      &grid->size,
      &grid->spacing,
      grid->span,
      grid->origin,
      energy->avdw,
      energy->bvdw,
      energy->es,
#ifdef BUILD_FOR_RDSOL
      energy->dslb,
      energy->dslx,
#endif
      bump->grid,
      &grid->psize,
      &grid->pspacing,
      grid->pspan,
      grid->porigin,
      energy->phi
    );


    energy->atom_model = 'a';

    energy->attractive_exponent = 6;
    energy->repulsive_exponent = 12;

    return;

}




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
  int i,j;
  STRING100 grid_file_name;
  FILE *grid_file;
  char new_fname[]="chem";

/*
* Write out bump grid
* 6/95 te


  int file_pref_len = strlen(grid->file_prefix);
  char new_fname[file_pref_len-2],ch;
  j=0;
  while (j < file_pref_len-2)
	    j++;
	    ch = grid->file_prefix[j];
	    new_fname[j]=ch;
	    printf("%d",j);
  
  
//printf("%c",grid->file_prefix[file_pref_len-2]);
*/
//  strcat (strcpy (grid_file_name, grid->file_prefix), ".bmp5.2");
  strcat (strcpy (grid_file_name, new_fname), "52.bmp");
  fprintf (global.outfile, "Writing general grid info to %s\n", grid_file_name);
  grid_file = efopen (grid_file_name, "w", global.outfile);

  fwrite (&grid->size, sizeof (int), 1, grid_file);
  fwrite (&grid->spacing, sizeof (float), 1, grid_file);
  fwrite (grid->origin, sizeof (float), 3, grid_file);
  fwrite (grid->span, sizeof (int), 3, grid_file);
  
  //printf(" origin %f \n", grid->origin[0]);

  if (bump->flag || contact->flag || chemical->flag)
  {
    fprintf (global.outfile, "Writing bump grid to %s\n", grid_file_name);
    fflush (global.outfile);

    fwrite
    (
      bump->grid,
      sizeof (unsigned char),
      grid->size,
      grid_file
    );
  }

  fclose (grid_file);

/* 
 * Write out Delphi map in Dock5 readable format -- kxr 03/05  
*/
    strcat (strcpy (grid_file_name, new_fname), "52.phi");
    fprintf (global.outfile, "Writing Delphi grid info to %s\n", grid_file_name);
    grid_file = efopen (grid_file_name, "w", global.outfile);
	                                                                                                                             
    fwrite (&grid->psize, sizeof (int), 1, grid_file);
    fwrite (&grid->pspacing, sizeof (float), 1, grid_file);
    fwrite (grid->porigin, sizeof (float), 3, grid_file);
    fwrite (grid->pspan, sizeof (int), 3, grid_file);

    fprintf (global.outfile, "Writing Delphi grid to %s\n", grid_file_name);
    fflush (global.outfile);

/* valgrind detects memory problems here; see the cvs log; srb apr 30, 2009 */

    fwrite
    (
      energy->phi,
	sizeof (float),
	grid->size,
	grid_file
    );

    fclose(grid_file);



/*
* Write out contact grid
* 6/95 te
*/

  if (contact->flag)
  {
    fprintf (global.outfile, "Not converting the contact grid\n" );
    fprintf (global.outfile, "It should already be in DOCK readable format\n" );
  }
    
/*
* Write out chemical grids
* 12/96 te
*/
  if (chemical->flag)
  {
//  strcat (strcpy (grid_file_name, grid->file_prefix), ".chm");
    strcat (strcpy (grid_file_name, new_fname), "52.chm");
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
//  strcat (strcpy (grid_file_name, grid->file_prefix), ".nrg");
    strcat (strcpy (grid_file_name, new_fname), "52.cmg");
    fprintf (global.outfile, "Writing energy grids to %s\n", grid_file_name);
    grid_file = efopen (grid_file_name, "w", global.outfile);
    fflush (global.outfile);

    fwrite (&grid->size, sizeof (int), 1, grid_file);
    fwrite (&grid->spacing, sizeof (float), 1, grid_file);
    fwrite (grid->origin, sizeof (float), 3, grid_file);
    fwrite (grid->span, sizeof (int), 3, grid_file);

//    fwrite (&grid->size, sizeof (int), 1, grid_file);
    fwrite (&energy->atom_model, sizeof (int), 1, grid_file);
    fwrite (&energy->attractive_exponent, sizeof (int), 1, grid_file);
    fwrite (&energy->repulsive_exponent, sizeof (int), 1, grid_file);

    fprintf (global.outfile, "  Writing attractive VDW energy grid\n");
    fflush (global.outfile);
    fwrite
    (
      energy->bvdw,
      sizeof (float),
      grid->size,
      grid_file
    );

    fprintf (global.outfile, "  Writing repulsive VDW energy grid\n");
    fflush (global.outfile);
    fwrite
    (
      energy->avdw,
      sizeof (float),
      grid->size,
      grid_file
    );

    fprintf (global.outfile, "  Writing electrostatic energy grid\n");
    fflush (global.outfile);
    fwrite
    (
      energy->es,
      sizeof (float),
      grid->size,
      grid_file
    );

    fclose (grid_file);
  }

#ifdef BUILD_FOR_RDSOL

// Write occupancy based desolvation grids
  if(energy->flag)
  {
//  strcat (strcpy (grid_file_name, grid->file_prefix), ".nrg");
    strcat (strcpy (grid_file_name, new_fname), "52.dsl");
    fprintf (global.outfile, "Writing desolvation grids to %s\n", grid_file_name);
    grid_file = efopen (grid_file_name, "w", global.outfile);
    fflush (global.outfile);

    fwrite (&grid->size, sizeof (int), 1, grid_file);
    fwrite (&grid->spacing, sizeof (float), 1, grid_file);
    fwrite (grid->origin, sizeof (float), 3, grid_file);
    fwrite (grid->span, sizeof (int), 3, grid_file);

    fprintf (global.outfile, "  Writing bulk desolvation energy grid\n");
    fflush (global.outfile);
    fwrite
    (
      energy->dslb,
      sizeof (float),
      grid->size,
      grid_file
    );

    fprintf (global.outfile, "  Writing explicit desolvation energy grid\n");
    fflush (global.outfile);
    fwrite
    (
      energy->dslx,
      sizeof (float),
      grid->size,
      grid_file
    );

    fclose (grid_file);
  }
#endif

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

/*  if (grid->init_flag)
    free_receptor_grid (grid);
*/
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
#ifdef BUILD_FOR_RDSOL
    efree ((void **) &energy->dslb);
    efree ((void **) &energy->dslx);
#endif
  }
}

