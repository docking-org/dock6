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
#include "label_node.h"
#include "label_vdw.h"

/* ==================================================================== */
void get_vdw_labels (LABEL_VDW *label_vdw)
{
  int i, definition_total, definition_count;
  STRING100 line, model;
  FILE *vdw_file;

  label_vdw->init_flag = TRUE;

  vdw_file = efopen (label_vdw->file_name, "r", global.outfile);

/*
* Count up the number of label_vdw->member declarations and definitions,
* then allocate memory
* 6/95 te
*/
  label_vdw->total = 1;
  definition_count = 0;

  while (fgets (line, 100, vdw_file) != NULL)
  {
    if (!strncmp (line, "name", 4))
      label_vdw->total++;
    if (!strncmp (line, "definition", 10))
      definition_count++;
  }

  rewind (vdw_file);
  definition_total = definition_count;

  ecalloc
  (
    (void **) &label_vdw->member,
    label_vdw->total,
    sizeof (VDW_MEMBER),
    "label_vdw->member structures",
    global.outfile
  );

  ecalloc
  (
    (void **) &label_vdw->definition,
    definition_count,
    sizeof (NODE),
    "atom label_vdw->member definitions",
    global.outfile
  );

/*
* Read in the atom label_vdw->member definitions
* 6/95 te
*/
  strcpy (label_vdw->member[0].name, "null");

  for (i = 0; i < label_vdw->total; i++)
    label_vdw->member[i].valence = INT_MAX;

  label_vdw->total = 1;
  definition_count = definition_total = 0;
  while (fgets (line, 100, vdw_file) != NULL)
  {
/*
*   Process label_vdw->member declaration
*   6/95 te
*/
    if (!strncmp (line, "name", 4))
    {
      label_vdw->member[label_vdw->total - 1].definition_total =
        definition_total;
      definition_total = 0;
      label_vdw->member[label_vdw->total].definition =
        &label_vdw->definition[definition_count];

      if (sscanf (line, "%*s %s", label_vdw->member[label_vdw->total].name) < 1)
      {
        fprintf (global.outfile, "Incomplete vdw member declaration.\n");
        exit (EXIT_FAILURE);
      }

      label_vdw->total++;
    }

    else if (!strncmp (line, "atom_model", 10))
    {
      if (sscanf (line, "%*s %s", model) != 1)
      {
        fprintf (global.outfile, "Incomplete atom_model specification.\n");
        exit (EXIT_FAILURE);
      }

      label_vdw->member[label_vdw->total - 1].atom_model = tolower (model[0]);

      if ((label_vdw->member[label_vdw->total - 1].atom_model != 'a') &&
        (label_vdw->member[label_vdw->total - 1].atom_model != 'u') &&
        (label_vdw->member[label_vdw->total - 1].atom_model != 'e'))
      {
        fprintf (global.outfile,
          "Atom_model specification restricted to ALL, UNITED, or EITHER.\n");
        exit (EXIT_FAILURE);
      }
    }

    else if (!strncmp (line, "heavy_flag", 8))
    {
      if (sscanf (line, "%*s %d",
        &label_vdw->member[label_vdw->total - 1].heavy_flag) != 1)
      {
        fprintf (global.outfile, "Incomplete heavy_flag specification.\n");
        exit (EXIT_FAILURE);
      }
    }

    else if (!strncmp (line, "radius", 6))
    {
      if (sscanf (line, "%*s %f",
        &label_vdw->member[label_vdw->total - 1].radius) != 1)
      {
        fprintf (global.outfile, "Incomplete radius specification.\n");
        exit (EXIT_FAILURE);
      }
    }

    else if (!strncmp (line, "well_depth", 10))
    {
      if (sscanf (line, "%*s %f",
        &label_vdw->member[label_vdw->total - 1].well_depth) != 1)
      {
        fprintf (global.outfile, "Incomplete well_depth specification.\n");
        exit (EXIT_FAILURE);
      }
    }

    else if (!strncmp (line, "valence", 7))
    {
      if (sscanf (line, "%*s %d",
        &label_vdw->member[label_vdw->total - 1].valence) != 1)
      {
        fprintf (global.outfile, "Incomplete valence specification.\n");
        exit (EXIT_FAILURE);
      }
    }

/*
*   Process label_vdw->member definition
*   6/95 te
*/
    else if (!strncmp (line, "definition", 10))
    {
      strtok (white_line (line), " ");

      if (!assign_node
        (&label_vdw->definition[definition_count], TRUE))
      {
        fprintf (global.outfile,
          "Error assigning vdw member definitions.\n");
        exit (EXIT_FAILURE);
      }

      definition_count++;
      definition_total++;
    }
  }

/*
* Update last label_vdw->member info also
* 6/95 te
*/
  label_vdw->member[label_vdw->total - 1].definition_total = definition_total;

  efclose (&vdw_file);

/*
* Calculate parameters
* 6/95 te
*/
  for (i = 0; i < label_vdw->total; i++)
  {
    if (label_vdw->member[i].heavy_flag)
      label_vdw->member[i].bump_id =
        NINT (10.0 * label_vdw->member[i].radius);

    else
      label_vdw->member[i].bump_id = 0;
  }

/*
* Print out the label_vdw->member and their definitions
* 6/95 te
*/

/*
  if (global.output_volume == 'v')
  {
    fprintf (global.outfile, "\n____VDW_Label_Definitions____\n\n");
    for (i = 1; i < label_vdw->total; i++)
    {
      MFPUTC ('_', 40);
      fprintf (global.outfile, "\n");
      fprintf (global.outfile, "%-20s%s\n", "name",
        label_vdw->member[i].name);
      fprintf (global.outfile, "%-20s%8c\n", "atom_model",
        label_vdw->member[i].atom_model);
      fprintf (global.outfile, "%-20s%8.3f\n", "radius",
        label_vdw->member[i].radius);
      fprintf (global.outfile, "%-20s%8.3f\n", "well_depth",
        label_vdw->member[i].well_depth);
      fprintf (global.outfile, "%-20s%8d\n", "heavy_flag",
        label_vdw->member[i].heavy_flag);
      fprintf (global.outfile, "%-20s%8d\n", "bump_id",
        label_vdw->member[i].bump_id);

      fprintf (global.outfile, "\n");

      for (j = 0; j < label_vdw->member[i].definition_total; j++)
      {
        fprintf (global.outfile, "%-20s", "definition");
        print_node (&label_vdw->member[i].definition[j], 0);
        fprintf (global.outfile, "\n");
      }

      fprintf (global.outfile, "\n");
    }

    fprintf (global.outfile, "\n\n");
  }
*/
}


/* ////////////////////////////////////////////////////////////////////// */

void free_vdw_labels (LABEL_VDW *label_vdw)
{
  int i, j;

  for (i = 0; i < label_vdw->total; i++)
    for (j = 0; j < label_vdw->member[i].definition_total; j++)
      free_node (&label_vdw->member[i].definition[j]);

  efree ((void **) &label_vdw->member);
  efree ((void **) &label_vdw->definition);
}

/* ////////////////////////////////////////////////////////////////////// */

int assign_vdw_labels
(
  LABEL_VDW *label_vdw,
  int atom_model,
  MOLECULE *molecule
)
{
  int i, j, k;
  int vdw_assigned = FALSE;

  if (!label_vdw->init_flag)
    get_vdw_labels (label_vdw);

  for (i = 0; i < molecule->total.atoms; i++)
  {
    for (j = 1; j < label_vdw->total; j++)
    {
      if ((atom_model == 'a') &&
        (label_vdw->member[j].atom_model == 'u'))
          continue;

      if ((atom_model == 'u') &&
        (label_vdw->member[j].atom_model == 'a'))
          continue;

      for (k = 0; k < label_vdw->member[j].definition_total; k++)
        if (check_atom
          (molecule, i, &label_vdw->member[j].definition[k]))
        {
          molecule->atom[i].vdw_id = j;
          vdw_assigned = TRUE;
        }
    }

    if (!vdw_assigned)
    {
      fprintf (global.outfile,
        "WARNING assign_vdw_labels: "
        "No vdw parameters for %1d %s %s %1d %s %s.\n",
        i + 1, molecule->atom[i].name, molecule->atom[i].type,
        molecule->atom[i].subst_id,
        molecule->subst[molecule->atom[i].subst_id].name,
        molecule->info.name);
      return FALSE;
    }

    if (molecule->atom[i].neighbor_total > 
      label_vdw->member[molecule->atom[i].vdw_id].valence)
    {
      for (j = k = 0; j < molecule->atom[i].neighbor_total; j++)
        if (strcmp (molecule->atom[molecule->atom[i].neighbor[j].id].type, "LP")
          &&
          strcmp (molecule->atom[molecule->atom[i].neighbor[j].id].type, "Du"))
          k++;

      if (k > label_vdw->member[molecule->atom[i].vdw_id].valence)
      {
        fprintf (global.outfile,
          "WARNING assign_vdw_labels: Atom valence violated for "
          "%s %s %s %d %s %d.\n",
          molecule->info.name,
          molecule->info.comment,
          molecule->subst[molecule->atom[i].subst_id].name,
          molecule->atom[i].subst_id + 1,
          molecule->atom[i].name, i + 1);

        return FALSE;
      }
    }

/*
*   Check to see if partial charge should be transferred from an atom
*   that is excluded from the current model (ie. aliphatic hydrogen in
*   a united atom model) to one that is
*   6/95 te
*/
    if ((atom_model == 'u') &&
      (ABS (label_vdw->member[molecule->atom[i].vdw_id].well_depth) < 0.00001))
    {
      if (molecule->atom[i].neighbor_total == 1)
      {
        molecule->atom[molecule->atom[i].neighbor[0].id].charge +=
          molecule->atom[i].charge;
        molecule->atom[i].charge = 0.0;
      }

      else
      {
        fprintf (global.outfile,
          "WARNING assign_vdw_labels: "
          "Unable to transfer partial charge away from "
          "%s %s %s %d %s %d.\n",
          molecule->info.name,
          molecule->info.comment,
          molecule->subst[molecule->atom[i].subst_id].name,
          molecule->atom[i].subst_id + 1,
          molecule->atom[i].name, i + 1);

        molecule->atom[i].charge = 0.0;
      }
    }
  }

  for (i = molecule->transform.heavy_total = 0; i < molecule->total.atoms; i++)
  {
    molecule->atom[i].heavy_flag =
      label_vdw->member[molecule->atom[i].vdw_id].heavy_flag;

    if (molecule->atom[i].heavy_flag == TRUE)
      molecule->transform.heavy_total++;
  }

  if (global.output_volume == 'v')
  {
    fprintf (global.outfile, "____VDW_Assignments____\n");
    fprintf (global.outfile, "%s:\n\n", molecule->info.name);
    for (i = 0; i < molecule->total.atoms; i++)
      fprintf (global.outfile, "%6d %4s %5s %-20s %d %d\n", i + 1,
        molecule->atom[i].name,
        molecule->atom[i].type,
        label_vdw->member[molecule->atom[i].vdw_id].name,
        molecule->atom[i].heavy_flag,
        label_vdw->member[molecule->atom[i].vdw_id].bump_id);
    fprintf (global.outfile, "\n\n");
  }

  return TRUE;
}


/* ////////////////////////////////////////////////////////////////////// */

int count_heavies
(
  MOLECULE	*molecule,
  int		segment_id
)
{
  int i;
  int heavy_total;

  if (segment_id == NEITHER)
  {
    for (i = heavy_total = 0; i < molecule->total.atoms; i++)
      if (molecule->atom[i].heavy_flag)
        heavy_total++;

    return heavy_total;
  }

  else if ((segment_id >= 0) && (segment_id < molecule->total.segments))
  {
    for
    (
      i = heavy_total = 0;
      i < molecule->segment[segment_id].atom_total;
      i++
    )
      if (molecule->atom[molecule->segment[segment_id].atom[i]].heavy_flag)
        heavy_total++;

    return heavy_total;
  }

  else
    return NEITHER;
}

