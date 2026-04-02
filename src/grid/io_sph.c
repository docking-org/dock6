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
#include "io_sph.h"

int read_sph
(
  MOLECULE *molecule,
  FILE_NAME molecule_file_name,
  FILE *molecule_file
)
{
  int i;
  int flag;
  char *line = NULL;
  int return_value = TRUE;
  STRING20 temp;

  if (find_record (&line, "cluster", molecule_file) == NULL)
    return EOF;

  vstrcpy (&molecule->info.name, "site_points");
  vstrcpy (&molecule->info.comment, molecule_file_name);

  memset (temp, 0, sizeof (STRING20));
  molecule->total.atoms = atoi (strncpy (temp, &line[45], 5));

/*
* Allocate space for molecule components
* 2/96 te
*/
  reallocate_atoms (molecule);

/*
* Store site point information
* 6/95 te
*/
  for (i = 0; i < molecule->total.atoms; i++)
  {
/*
*   Read in line of file
*   6/95 te
*/
    if (vfgets (&line, molecule_file) == NULL)
      exit (fprintf (global.outfile,
        "ERROR read_sph: Insufficient sphere records in %s\n",
        molecule_file_name));

    if (strlen (line) < 52)
      exit (fprintf (global.outfile,
        "ERROR read_sph: Sphere file data line too short in %s\n%s",
        molecule_file_name, line));

    molecule->atom[i].number = i + 1;

    memset (temp, 0, sizeof (STRING20));
    sscanf (&line[0], "%s", temp);
    vstrcpy (&molecule->atom[i].name, temp);

    vstrcpy (&molecule->atom[i].type, "Du");

    memset (temp, 0, sizeof (STRING20));
    molecule->coord[i][0] = atof (strncpy (temp, &line[5], 10));

    memset (temp, 0, sizeof (STRING20));
    molecule->coord[i][1] = atof (strncpy (temp, &line[15], 10));

    memset (temp, 0, sizeof (STRING20));
    molecule->coord[i][2] = atof (strncpy (temp, &line[25], 10));

    memset (temp, 0, sizeof (STRING20));
    molecule->atom[i].charge = atof (strncpy (temp, &line[35], 8));

    memset (temp, 0, sizeof (STRING20));
    molecule->atom[i].subst_id = atoi (strncpy (temp, &line[48], 2)) - 1;

    memset (temp, 0, sizeof (STRING20));
    molecule->atom[i].chem_id = atoi (strncpy (temp, &line[50], 3));
  }

/*
* Convert zero-value critical id numbers to one above the maximum
* 10/95 te
*/
  for (i = 0, flag = FALSE; i < molecule->total.atoms; i++)
    if (molecule->atom[i].subst_id + 1 > molecule->total.substs)
      molecule->total.substs = molecule->atom[i].subst_id + 1;

  for (i = 0, flag = FALSE; i < molecule->total.atoms; i++)
    if (molecule->atom[i].subst_id < 0)
    {
      flag = TRUE;
      molecule->atom[i].subst_id = molecule->total.substs;
    }

  if (flag)
    molecule->total.substs++;

/*
* Allocate substructure records
* 10/95 te
*/
  reallocate_substs (molecule);

  for (i = 0; i < molecule->total.substs; i++)
  {
    molecule->subst[i].number = i + 1;
    vstrcpy (&molecule->subst[i].name, "SPH");
  }

  efree ((void **) &line);
  return return_value;
}


/* ////////////////////////////////////////////////////// */

int write_sph
(
  MOLECULE *molecule,
  FILE *molecule_file
)
{
  int i, j;
  STRING80 line;

/*
* Write header
* 6/95 te
*/
  fprintf (molecule_file, "cluster %5d   ", molecule->info.output_id);
  fprintf (molecule_file, "number of spheres in cluster %5d\n",
    molecule->total.atoms);

/*
* Write out atoms
* 6/95 te
*/
  for (i = 0; i < molecule->total.atoms; i++)
  {
/*
*   Print out site point info
*   6/95 te
*/
    memset (line, 0, sizeof (STRING80));

    if (molecule->atom[i].name)
      sprintf (&line[0], "%5d", i + 1);

    else
      sprintf (&line[0], "%5d", i + 1);

    sprintf (&line[5], "%10.5f", molecule->coord[i][0]);
    sprintf (&line[15], "%10.5f", molecule->coord[i][1]);
    sprintf (&line[25], "%10.5f", molecule->coord[i][2]);
    sprintf (&line[35], "%8.3f", molecule->atom[i].charge);
    sprintf (&line[48], "%2d", molecule->atom[i].subst_id + 1);
    sprintf (&line[50], "%3d", molecule->atom[i].chem_id);

/*
*   Replace null characters with spaces
*   6/95 te
*/
    for (j = 0; j < sizeof (STRING80); j++)
      if (line[j] == '\0') line[j] = ' ';
    line[sizeof (STRING80) - 2] = '\0';

    fprintf (molecule_file, "%s\n", line);
  }

  return TRUE;
}


