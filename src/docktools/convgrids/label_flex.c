/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
12/96
*/
#include "define.h"
#include "utility.h"
#include "mol.h"
#include "global.h"
#include "search.h"
#include "label_node.h"
#include "label_vdw.h"
#include "label_flex.h"
#include "transform.h"


/* ==================================================================== */
int get_flex_labels (LABEL_FLEX *label_flex)
{
  int i, definition_total, definition_count;
  STRING100 line;
  FILE *flex_file;

  label_flex->init_flag = TRUE;

  flex_file = efopen (label_flex->file_name, "r", global.outfile);

/*
* Count up the number of flex label declarations and definitions,
* then allocate memory
* 6/95 te
*/
  for (label_flex->total = definition_count = 0; fgets (line, 100, flex_file);)
  {
    if (!strncmp (line, "name", 4))
      label_flex->total++;
    if (!strncmp (line, "definition", 10))
      definition_count++;
  }

  rewind (flex_file);

  definition_total = definition_count;

  ecalloc
  (
    (void **) &label_flex->member,
    label_flex->total + 1,
    sizeof (FLEX_MEMBER),
    "flexible bond labels",
    global.outfile
  );

  ecalloc
  (
    (void **) &label_flex->definition,
    definition_total,
    sizeof (NODE),
    "flexible bond label definitions",
    global.outfile
  );

  strcpy (label_flex->member[0].name, "null");

/*
* Read in the atom label definitions
* 6/95 te
*/
  for (label_flex->total = 1, definition_count = definition_total = 0;
    fgets (line, 100, flex_file);)
  {
/*
*   Process flexible label declaration
*   6/95 te
*/
    if (!strncmp (line, "name", 4))
    {
      if (label_flex->total > 1)
      {
        if (label_flex->member[label_flex->total - 1].drive_id == NEITHER)
          exit (fprintf (global.outfile,
            "ERROR get_flex_labels: Missing drive_id parameter in %s\n",
            label_flex->file_name));

        else if (label_flex->member[label_flex->total - 1].minimize_flag ==
          NEITHER)
          exit (fprintf (global.outfile,
            "ERROR get_flex_labels: Missing minimize parameter in %s\n",
            label_flex->file_name));

        else if (definition_total != 2)
          exit (fprintf (global.outfile,
            "ERROR get_flex_labels: Incorrect number of definitions in %s\n",
            label_flex->file_name));
      }

      label_flex->member[label_flex->total - 1].definition_total =
        definition_total;
      label_flex->member[label_flex->total].definition =
        &label_flex->definition[definition_count];
      definition_total = 0;

      if (sscanf
        (line, "%*s %s", label_flex->member[label_flex->total].name) < 1)
      {
        fprintf (global.outfile,
          "Incomplete label_flex->member declaration.\n");
        exit (EXIT_FAILURE);
      }

/*
*     Convert flex label name to lowercase
*     6/95 te
*/
      for (i = 0; i < strlen (label_flex->member[label_flex->total].name); i++)
        label_flex->member[label_flex->total].name[i] =
          (char) tolower (label_flex->member[label_flex->total].name[i]);

      label_flex->total++;

      label_flex->member[label_flex->total - 1].drive_id = NEITHER;
      label_flex->member[label_flex->total - 1].minimize_flag = NEITHER;
    }

/*
*   Process label_flex->member search identifier
*   10/95 te
*/
    else if (!strncmp (line, "drive_id", 6))
      sscanf
        (line, "%*s %d", &label_flex->member[label_flex->total - 1].drive_id);

/*
*   Process flex label minimize flag
*   10/95 te
*/
    else if (!strncmp (line, "minimize", 8))
      sscanf
        (line, "%*s %d",
        &label_flex->member[label_flex->total - 1].minimize_flag);

/*
*   Process flex label definition
*   6/95 te
*/
    else if (!strncmp (line, "definition", 10))
    {
      strtok (white_line (line), " ");

      if (!assign_node
        (&label_flex->definition[definition_count], TRUE))
      {
        fprintf (global.outfile, "Error assigning flex label definitions.\n");
        exit (EXIT_FAILURE);
      }

      definition_count++;
      definition_total++;
    }
  }

/*
* Update last flex label info also
* 6/95 te
*/
  label_flex->member[label_flex->total - 1].definition_total =
    definition_total;

  efclose (&flex_file);

/*
* Assign search torsion parameters if requested
* 10/96 te
*/
  if (label_flex->drive_flag)
    get_flex_search (label_flex);

/*
* Print out the flex label and their definitions
* 6/95 te

  if (global.output_volume == 'v')
  {
    fprintf (global.outfile, "\n____Flexible_Bond_Label_Definitions____\n\n");
    for (i = 0; i < label_flex->total; i++)
    {
      fprintf (global.outfile, "\n");
      fprintf (global.outfile, "%-20s%s\n", "name",
        label_flex->member[i].name);
      fprintf (global.outfile, "%-20s%d\n", "drive_id",
        label_flex->member[i].drive_id);
      fprintf (global.outfile, "%-20s%d\n", "minimize",
        label_flex->member[i].minimize_flag);

      fprintf (global.outfile, "%-20s%d\n", "positions",
        label_flex->member[i].torsion_total);

      fprintf (global.outfile, "%-20s", "torsions");
      for (j = 0; j < label_flex->member[i].torsion_total; j++)
        fprintf (global.outfile, "%g ", label_flex->member[i].torsion[j]);

      fprintf (global.outfile, "\n");
      for (j = 0; j < label_flex->member[i].definition_total; j++)
      {
        fprintf (global.outfile, "%-20s", "definition");
        print_node (&label_flex->member[i].definition[j], 0);
        fprintf (global.outfile, "\n");
      }

      fprintf (global.outfile, "\n");
    }

    fprintf (global.outfile, "\n\n");
  }
*/

  return TRUE;
}


/* ////////////////////////////////////////////////////////////////////// */

void free_flex_labels (LABEL_FLEX *label_flex)
{
  int i, j;

  for (i = 0; i < label_flex->total; i++)
  {
    efree ((void **) &label_flex->member[i].torsion);

    for (j = 0; j < label_flex->member[i].definition_total; j++)
      free_node (&label_flex->member[i].definition[j]);
  }

  efree ((void **) &label_flex->member);
  efree ((void **) &label_flex->definition);
}


/* ///////////////////////////////////////////////////////////////////

Program to read in torsion.defn
11/95 Yax Sun
10/96 edited by Todd Ewing

  //////////////////////////////////////////////////////////////////// */

int get_flex_search (LABEL_FLEX *label_flex)
{
  int i, j;			/* Counter variables */
  STRING200 line;		/* String used to compose output */
  char *token;
  FILE  *search_file;

  int   drive_id;
  int   torsion_total;
  float *torsion = NULL;	/* temporarily hold the torsion values */

  search_file = efopen (label_flex->search_file_name, "r", global.outfile);

  while (fgets (line, 200, search_file))
  {
    token = strtok (white_line (line), " ");

    if (!strcmp (token, "drive_id"))
    {
      if (token = strtok (NULL, " "))
        drive_id = atoi (token);

      else
        exit (fprintf (global.outfile,
          "ERROR get_flex_search: Search_id value not specified in %s\n",
          label_flex->search_file_name));

      if (!fgets (line, 200, search_file) ||
        !(token = strtok (white_line (line), " ")) ||
        strcmp (token, "positions"))
        exit (fprintf (global.outfile,
          "ERROR get_flex_search: Positions field doesn't follow Id in %s\n",
          label_flex->search_file_name));

      if (token = strtok (NULL, " "))
        torsion_total = atoi (token);

      else
        exit (fprintf (global.outfile,
          "ERROR get_flex_search: Postions value not specified in %s\n",
          label_flex->search_file_name));

      ecalloc
      (
        (void **) &torsion,
        torsion_total,
        sizeof (float),
        "torsion sample values",
        global.outfile
      );

      if (!fgets (line, 200, search_file) ||
        !(token = strtok (white_line (line), " ")) ||
        strcmp (token, "torsions"))
        exit (fprintf (global.outfile,
          "ERROR get_flex_search: Torsions doesn't follow Positions in %s\n",
          label_flex->search_file_name));

      for (i = 0; i < torsion_total; i++)
      {
        if (token = strtok (NULL, " "))
          torsion[i] = atof (token);

        else
          exit (fprintf (global.outfile,
            "ERROR get_flex_search: Insufficient number of torsions in %s\n",
            label_flex->search_file_name));
      }

      for (i = 0; i < label_flex->total; i++)
      {
        if (label_flex->member[i].drive_id == drive_id)
        {
          label_flex->member[i].torsion_total = torsion_total;

          ecalloc
          (
            (void **) &label_flex->member[i].torsion,
            label_flex->member[i].torsion_total,
            sizeof (float),
            "stored torsion sample values",
            global.outfile
          );

          for (j = 0; j < label_flex->member[i].torsion_total; j++)
            label_flex->member[i].torsion[j] = torsion[j];

/*
          fprintf (global.outfile,
            "search %d, %s\n", drive_id, label_flex->member[i].name);
*/
        }
      }

      torsion_total = 0;
      efree ((void **) &torsion);
    }
  }

  for (i = 1; i < label_flex->total; i++)
    if (label_flex->member[i].torsion_total < 1)
      exit (fprintf (global.outfile,
        "ERROR get_flex_search: Missing torsion parameters in %s\n",
        label_flex->search_file_name));

  fclose (search_file);
  return TRUE;
}


/* //////////////////////////////////////////////////////////////////// */

void assign_flex_labels
(
  LABEL_FLEX	*label_flex,
  MOLECULE	*molecule
)
{
  int i, j;
  int bond_id, bond_found;
  int flex_id;
  float conf_total;
  static MOLECULE temporary = {0};

  void detect_rings (MOLECULE *);

  if (!label_flex->init_flag)
    get_flex_labels (label_flex);

  reset_molecule (&temporary);

/*
* Make a copy of the original set of torsion bonds
* 10/96 te
*/
  copy_torsions (&temporary, molecule);

/*
* Flag all ring bonds
* 10/96 te
*/
  detect_rings (molecule);

/*
* Identify all bonds with flex labels
* 10/96 te
*/
  for (i = molecule->total.torsions = 0; i < molecule->total.bonds; i++)
  {
    molecule->bond[i].flex_id = 0;

    if (molecule->bond[i].ring_flag == TRUE)
      continue;

    if ((molecule->atom[molecule->bond[i].origin].neighbor_total <= 1) ||
      (molecule->atom[molecule->bond[i].target].neighbor_total <= 1))
      continue;

    for (j = 1; j < label_flex->total; j++)
    {
      if
      (
        check_atom
          (molecule, molecule->bond[i].origin,
          &label_flex->member[j].definition[0]) &&
        check_atom
          (molecule, molecule->bond[i].target,
          &label_flex->member[j].definition[1])
      )
        molecule->bond[i].flex_id = j;

      else if
      (
        check_atom
          (molecule, molecule->bond[i].origin,
          &label_flex->member[j].definition[1]) &&
        check_atom
          (molecule, molecule->bond[i].target,
          &label_flex->member[j].definition[0])
      )
        molecule->bond[i].flex_id = j;
    }

    if (molecule->bond[i].flex_id)
      molecule->total.torsions++;
  }

/*
* Check for sufficient space to store new and old torsions
* 10/96 te
*/
  if ((molecule->total.torsions + temporary.total.torsions) >
    molecule->max.torsions)
  {
    molecule->total.torsions += temporary.total.torsions;
    reallocate_torsions (molecule);
  }

  else
    reset_torsions (molecule);

/*
* Store flexible torsions
* 10/96 te
*/
  for (i = molecule->total.torsions = 0; i < molecule->total.bonds; i++)
    if (molecule->bond[i].flex_id)
    {
      molecule->torsion[molecule->total.torsions].bond_id = i;
      molecule->torsion[molecule->total.torsions].flex_id =
        molecule->bond[i].flex_id;
      molecule->total.torsions++;
    }

/*
* Append original torsions
* 10/96 te
*/
  for (i = 0; i < temporary.total.torsions; i++)
  {
    for (j = 0, bond_found = FALSE; j < molecule->total.torsions; j++)
      if (temporary.torsion[i].bond_id == molecule->torsion[j].bond_id)
      {
        molecule->torsion[j].target_angle = temporary.torsion[i].target_angle;
        bond_found = TRUE;
      }

    if (!bond_found)
      copy_torsion
        (&molecule->torsion[molecule->total.torsions++],
        &temporary.torsion[i]);
  }

/*
* Compute current torsion angles and conformation total
* 10/96 te
*/
  get_torsion_neighbors (molecule);

  for (i = 0, conf_total = 1.0; i < molecule->total.torsions; i++)
  {
    molecule->torsion[i].current_angle =
      molecule->torsion[i].target_angle =
      compute_torsion (molecule, i);

    molecule->torsion[i].periph_flag =
      check_peripheral_torsion (molecule, i);

    bond_id = molecule->torsion[i].bond_id;

    if (label_flex->drive_flag)
    {
      flex_id = molecule->bond[bond_id].flex_id;

      if (flex_id > 0)
      {
        if (label_flex->member[flex_id].torsion_total < 1)
          exit (fprintf (global.outfile,
            "ERROR assign_flex_labels: bond with no positions detected\n"));

        conf_total *= (float) label_flex->member[flex_id].torsion_total;
      }
    }
  }

/*
* Print out flex label assignments
* 10/95 te
*/
  if (global.output_volume == 'v')
  {
    fprintf (global.outfile, "____Flexible_Bond_Label_Assignments____\n");
    fprintf (global.outfile, "%s:\n\n", molecule->info.name);
    fprintf (global.outfile,
      "%-4s %-4s | %-4s %-4s %-4s %-4s | %7s %3s %3s %s\n",
      "tors", "bond", "orig", "orig", "targ", "targ", "input", "pos", "drv", "flex");
    fprintf (global.outfile,
      "%-4s %-4s | %-4s %-4s %-4s %-4s | %7s %3s %3s %s\n",
      " id", " id", "nhbr", "", "", "nhbr", "angle", "num", "id", "label");
    fprintf (global.outfile,
      "-----------------------------------------"
      "---------------------------\n");

    for (i = 0; i < molecule->total.torsions; i++)
    {
      bond_id = molecule->torsion[i].bond_id;

      fprintf
      (
        global.outfile,
        "%-4d %-4d | %-4s %-4s %-4s %-4s | %7.2f %3d %3d %s\n",
        i + 1, bond_id + 1,
        molecule->atom[molecule->torsion[i].origin_neighbor].name,
        molecule->atom[molecule->torsion[i].origin].name,
        molecule->atom[molecule->torsion[i].target].name,
        molecule->atom[molecule->torsion[i].target_neighbor].name,
        molecule->torsion[i].current_angle / PI * 180,
        label_flex->member[molecule->bond[bond_id].flex_id].torsion_total,
        molecule->bond[bond_id].flex_id ?
          label_flex->member[molecule->bond[bond_id].flex_id].drive_id :
          0,
        molecule->bond[bond_id].flex_id ?
          label_flex->member[molecule->bond[bond_id].flex_id].name :
          "INFLEXIBLE"
      );
    }

    if (label_flex->drive_flag)
      fprintf (global.outfile,
        "\nTotal conformations (ignoring symmetry and clashes) for %s %.1g\n",
        molecule->info.name, conf_total);
  }
}
 

/* //////////////////////////////////////////////////////////////////////

Subroutine to identify all bonds in rings.

11/96 te

////////////////////////////////////////////////////////////////////// */

void detect_rings (MOLECULE *molecule)
{
  int i, j;

  int next_atom_level (int, int, int, MOLECULE *);

/*
* Reset atom flags and bond ring flags
* 4/97 te
*/
  for (i = 0; i < molecule->total.atoms; i++)
    molecule->atom[i].flag = 0;

  for (i = 0; i < molecule->total.bonds; i++)
    molecule->bond[i].ring_flag = 0;

/*
* Identify ring bonds
* 4/97 te
*/
  next_atom_level (1, 0, 0, molecule);

/*
* Also identify bonds specified in RIGID set, if present
* 4/97 te
*/
  for (i = 0; i < molecule->total.sets; i++)
  {
    if ((molecule->set[i].name != NULL) &&
      !strcmp (molecule->set[i].name, "RIGID"))
    {
      if ((molecule->set[i].type != NULL) &&
        !strcmp (molecule->set[i].type, "STATIC"))
      {
        if ((molecule->set[i].obj_type != NULL) &&
          !strcmp (molecule->set[i].obj_type, "BONDS"))
        {
          if (molecule->set[i].member_total > 0)
          {
            for (j = 0; j < molecule->set[i].member_total; j++)
            {
              if ((molecule->set[i].member[j] >= 0) &&
                (molecule->set[i].member[j] < molecule->total.bonds))
                molecule->bond[molecule->set[i].member[j]].ring_flag = TRUE;

              else
                exit (fprintf (global.outfile,
                  "ERROR detect_rings: RIGID set contains invalid bonds\n"));
            }
          }

          else
            exit (fprintf (global.outfile,
              "ERROR detect_rings: RIGID set empty\n"));
        }

        else
          exit (fprintf (global.outfile,
            "ERROR detect_rings: RIGID set must be BONDS obj_type\n"));
      }

      else
        exit (fprintf (global.outfile,
          "ERROR detect_rings: RIGID set must be STATIC type\n"));
    }
  }

/*
  for (i = 0; i < molecule->total.bonds; i++)
    if (molecule->bond[i].ring_flag)
      fprintf (global.outfile, "Ring bond %d %s %s\n", i + 1,
        molecule->atom[molecule->bond[i].origin].name,
        molecule->atom[molecule->bond[i].target].name);
*/
}

 
/* //////////////////////////////////////////////////////////////////////

Recursive subroutine that traverses all atom linkage paths in a molecule.
The level of recursion is recorded for each atom as it is traversed.
The recursion level is reported back, unless an atom is re-encountered
during a higher level of recursion; then, the level of the previously seen
atom is reported.

A ring bond is identified when an atom is linked to an atom previously seen.

4/97 te

////////////////////////////////////////////////////////////////////// */

int next_atom_level
(
  int level,
  int current_atom,
  int previous_atom,
  MOLECULE *molecule
)
{
  int i;
  int bond_id;			/* Bond id linking current neighbor */
  int neighbor_flag;		/* Flag of neighbor atom */

/*
* If the current atom hasn't been flagged, then check the flags
* of its neighbors
* 4/97 te
*/
  if (molecule->atom[current_atom].flag == 0)
  {
    molecule->atom[current_atom].flag = level;        

    for (i = 0; i < molecule->atom[current_atom].neighbor_total; i++)
    {
      if (molecule->atom[current_atom].neighbor[i].id != previous_atom) 
      {
        neighbor_flag =
          next_atom_level
          (
            level + 1,
            molecule->atom[current_atom].neighbor[i].id,
            current_atom,
            molecule
          );

        if (neighbor_flag <= level)
        {
          bond_id = molecule->atom[current_atom].neighbor[i].bond_id;
          molecule->bond[bond_id].ring_flag = TRUE;
        }

        if (neighbor_flag < molecule->atom[current_atom].flag)
          molecule->atom[current_atom].flag = neighbor_flag;
      }
    }
  }

  return molecule->atom[current_atom].flag;
}


/* /////////////////////////////////////////////////////////////////// */

int check_peripheral_torsion
(
  MOLECULE	*molecule,
  int		torsion_id
)
{
  int atom_id;

  atom_id = molecule->torsion[torsion_id].origin_neighbor;
  if (molecule->atom[atom_id].heavy_flag == FALSE)
    return TRUE;

  atom_id = molecule->torsion[torsion_id].target_neighbor;
  if (molecule->atom[atom_id].heavy_flag == FALSE)
    return TRUE;

  return FALSE;
}


/* /////////////////////////////////////////////////////////////////// */
/*
Routine to assign torsions to segments and segments to layers.
It also loops through all anchor fragments.

Return values:
  TRUE	successful identification of anchor segment
  EOF	unable to identify any more anchor segments

12/96 te

*/
/* /////////////////////////////////////////////////////////////////// */

int get_anchor
(
  LABEL_FLEX	*label_flex,
  MOLECULE	*mol_init,
  MOLECULE	*mol_anch,
  int		anchor
) 
{
/*
* Initialize variables
* 12/96 te
*/
  if (anchor == 0)
  {
/*
*   Assign flexible labels
*   12/96 te
*/
    if (label_flex->flag)
    {
      mol_init->transform.tors_flag = FALSE;
      assign_flex_labels (label_flex, mol_init);

      if (mol_init->total.torsions > label_flex->max_torsions)
      {
        if (global.output_volume != 't')
          fprintf(global.outfile, "Skipped %s (%d rotatable bonds).\n",
            mol_init->info.name, mol_init->total.torsions);

        return EOF;
      }
    }

    get_segments (label_flex, mol_init);
    copy_molecule (mol_anch, mol_init);
  }

/*
* Exit if return visit, but multiple anchors not requested
* 12/96 te
*/
  else if (!label_flex->multiple_anchors)
    return EOF;

  else
    copy_segments (mol_anch, mol_init);

/*
* Find next anchor
* 12/96 te
*/
  if (get_anchor_segment (label_flex, mol_anch, anchor) != TRUE)
    return EOF;

  initialize_segments (label_flex, mol_anch);
  initialize_anchor_layer (label_flex, mol_anch);

  return TRUE;
}


/* /////////////////////////////////////////////////////////////////////

Routine to divide a molecule into rigid segments.

11/96 te

///////////////////////////////////////////////////////////////////// */

void get_segments (LABEL_FLEX *label_flex, MOLECULE *mol_init)
{
  int i, j;
  int atom_id;			/* Atom id */
  int bond_id;			/* Bond id */
  int torsion_id;		/* Torsion id */
  int origin;			/* Origin atom id */
  int target;			/* Target atom id */
  int segment;			/* Current segment id */
  int neighbor;			/* Neighboring segment id */
  int neighbor_id;		/* Current position in neighbor list */

  static SEARCH search = {0};

/*
* Allocate more than enough space for segments
* 11/96 te
*/
  mol_init->total.segments = mol_init->total.torsions + 1;
  reallocate_segments (mol_init);

  for (i = 0; i < mol_init->total.segments; i++)
  {
    mol_init->segment[i].atom_total = mol_init->total.atoms;
    mol_init->segment[i].neighbor_total = mol_init->total.torsions;

    reallocate_segment_atoms (&mol_init->segment[i]);
    reallocate_segment_neighbors (&mol_init->segment[i]);

    reset_segment (&mol_init->segment[i]);
  }

/*
* Identify rigid segments separated by rotatable bonds
* 11/96 te
*/
  if (mol_init->total.torsions == 0)
  {
    for (i = 0; i < mol_init->total.atoms; i++)
    {
      mol_init->segment[0].atom[i] = i;
      mol_init->atom[i].segment_id = 0;
    }

    mol_init->segment[0].atom_total = mol_init->total.atoms;
  }

  else
  {
    for
    (
      i = mol_init->total.segments = 0;
      breadth_search
      (
        &search,
        mol_init->atom,
        mol_init->total.atoms,
        get_atom_neighbor,
        NULL,
        &i, 1, NEITHER, i
      ) != EOF;
      i++
    )
    {
      origin = get_search_origin (&search, NEITHER, NEITHER, NEITHER);
      target = get_search_target (&search, NEITHER, NEITHER, NEITHER);

/*
*     If this is the first atom, then initialize segment list only
*     11/96 te
*/
      if (i == 0)
      {
        mol_init->atom[target].segment_id = 0;
        mol_init->segment[0].atom[0] = target;
        mol_init->segment[0].atom_total = 1;
        mol_init->total.segments = 1;
        continue;
      }

/*
*     Check if this linkage is a rotatable bond
*     11/96 te
*/
      for
      (
        j = 0, torsion_id = NEITHER;
        (j < mol_init->total.torsions) && (torsion_id == NEITHER);
        j++
      )
        if (mol_init->torsion[j].flex_id > 0)
        {
          bond_id = mol_init->torsion[j].bond_id;

          if (((mol_init->bond[bond_id].origin == origin) &&
            (mol_init->bond[bond_id].target == target)) ||
            ((mol_init->bond[bond_id].origin == target) &&
            (mol_init->bond[bond_id].target == origin)))
            torsion_id = j;
        }

/*
*     If the link is flexible, then increment the segment total,
*     update the segment neighbor lists and atom list, and atom segment id
*     11/96 te
*/
      if (torsion_id != NEITHER)
      {
        segment = mol_init->total.segments++;
        neighbor = mol_init->atom[origin].segment_id;

        neighbor_id = mol_init->segment[segment].neighbor_total;
        mol_init->segment[segment].neighbor[neighbor_id].id = neighbor;
        mol_init->segment[segment].neighbor_total++;

        neighbor_id = mol_init->segment[neighbor].neighbor_total;
        mol_init->segment[neighbor].neighbor[neighbor_id].id = segment;
        mol_init->segment[neighbor].neighbor_total++;

        atom_id = mol_init->segment[segment].atom_total;
        mol_init->segment[segment].atom[atom_id] = target;
        mol_init->segment[segment].atom_total++;

        mol_init->atom[target].segment_id = segment;

        mol_init->segment[segment].periph_flag =
          mol_init->torsion[torsion_id].periph_flag;
      }

/*
*     Otherwise, update the segment atom list and atom segment id
*     11/96 te
*/
      else
      {
        segment = mol_init->atom[origin].segment_id;

        atom_id = mol_init->segment[segment].atom_total;
        mol_init->segment[segment].atom[atom_id] = target;
        mol_init->segment[segment].atom_total++;

        mol_init->atom[target].segment_id = segment;
      }
    }
  }

/*
* Determine the size of each segment
* 1/97 te
*/
  for (i = 0; i < mol_init->total.segments; i++)
  {
    mol_init->segment[i].heavy_total = 0;

    for (j = 0; j < mol_init->segment[i].atom_total; j++)
    {
      atom_id = mol_init->segment[i].atom[j];

      if (mol_init->atom[atom_id].heavy_flag == TRUE)
        mol_init->segment[i].heavy_total++;
    }
  }

  if (label_flex->flag && (global.output_volume == 'v'))
  {
    fprintf (global.outfile, "\nInitial segments\n");

    for (i = 0; i < mol_init->total.segments; i++)
    {
      fprintf (global.outfile, "  segment %d\n", i + 1);
      fprintf (global.outfile, "    neighbors:");

      for (j = 0; j < mol_init->segment[i].neighbor_total; j++)
        fprintf (global.outfile, " %d",
          mol_init->segment[i].neighbor[j].id + 1);

      fprintf (global.outfile, "\n    atoms    :");

      for (j = 0; j < mol_init->segment[i].atom_total; j++)
        fprintf (global.outfile, " %s",
          mol_init->atom[mol_init->segment[i].atom[j]].name);

      fprintf (global.outfile, "\n");
    }
  }
}


/* /////////////////////////////////////////////////////////////////////

Routine to identify segments to use as anchors.

Return values:
	TRUE:	next anchor found
	EOF:	unable to find any more anchors

12/96 te

///////////////////////////////////////////////////////////////////// */

int get_anchor_segment
(
  LABEL_FLEX	*label_flex,
  MOLECULE	*molecule,
  int		anchor_count
)
{
  int i;
  int size;
  int anchor_found = FALSE;	/* Flag for whether anchor found */

  if (anchor_count > 0)
  {
    if (label_flex->multiple_anchors == FALSE)
      return EOF;
  }

  else
  {
    if (label_flex->anchor_flag == FALSE)
      return TRUE;

    else if (molecule->transform.anchor_atom != NEITHER)
      exit (fprintf (global.outfile,
        "ERROR get_anchor_segment: Molecule %s %s\n"
        "  This molecule has an anchor atom and transformation from input.\n"
        "  Convert input file to non-PTR format before processing.\n",
        molecule->info.name, molecule->info.comment));
  }

/*
* Check if an anchor was specified in input
* 1/97 te
*/
  for (i = 0; i < molecule->total.sets; i++)
    if (!strcmp (molecule->set[i].name, "ANCHOR"))
    {
      if ((strcmp (molecule->set[i].type, "STATIC")) ||
        (strcmp (molecule->set[i].obj_type, "ATOMS")))
        exit (fprintf (global.outfile,
          "ERROR get_anchor_segment: Molecule %s %s\n"
          "  This molecule has an ANCHOR set in input file.\n"
          "  ANCHOR set must be of type STATIC ATOMS.\n",
          molecule->info.name, molecule->info.comment));

      if (label_flex->multiple_anchors == TRUE)
        exit (fprintf (global.outfile,
          "ERROR get_anchor_segment: Molecule %s %s\n"
          "  This molecule has an ANCHOR atom set in input file.\n"
          "  Delete ANCHOR set before processing with multiple_anchors.\n",
          molecule->info.name, molecule->info.comment));

      if
      (
        (molecule->set[i].member_total > 0) &&
        (molecule->set[i].member[0] >= 0) &&
        (molecule->set[i].member[0] < molecule->total.atoms)
      )
      {
        molecule->transform.anchor_atom = molecule->set[i].member[0];
        return TRUE;
      }

      else
        exit (fprintf (global.outfile,
          "ERROR get_anchor_segment: Molecule %s %s\n"
          "  This molecule has an ANCHOR atom set in input file.\n"
          "  The ANCHOR set contains no atoms or incorrect atoms.\n",
          molecule->info.name, molecule->info.comment));
    }

/*
* If multiple anchors allowed, then find the next segment
* that satisfies the anchor size cutoff
* 1/97 te
*/
  if ((anchor_found == FALSE) &&
    (label_flex->multiple_anchors == TRUE))
  {
    for
    (
      anchor_found = FALSE,
        i = (anchor_count > 0) ?
          molecule->atom[molecule->transform.anchor_atom].segment_id + 1 : 0;
      (anchor_found == FALSE) &&
        (i < molecule->total.segments);
      i++
    )
      if ((molecule->segment[i].heavy_total +
        molecule->segment[i].neighbor_total) >=
        label_flex->anchor_size)
      {
        molecule->transform.anchor_atom = molecule->segment[i].atom[0];
        anchor_found = TRUE;
      }
  }

/*
* Identify the largest segment (counting all heavy atoms and one atom
* from each neighbor)
* 1/97 te
*/
  if
  (
    (anchor_found == FALSE) &&
    (((label_flex->multiple_anchors == TRUE) && (anchor_count == 0)) ||
    (label_flex->multiple_anchors == FALSE))
  )
  {
    for
    (
      i = 0, size = INT_MIN;
      i < molecule->total.segments;
      i++
    )
      if ((molecule->segment[i].heavy_total +
        molecule->segment[i].neighbor_total) > size)
      {
        size = molecule->segment[i].heavy_total +
          molecule->segment[i].neighbor_total;
        molecule->transform.anchor_atom = molecule->segment[i].atom[0];
      }

    anchor_found = TRUE;
  }

  if (anchor_found == TRUE)
    return TRUE;

  else
    return EOF;
}


/* /////////////////////////////////////////////////////////////////////

Routine to initialize segments by
1. updating atom segment ids
2. assigning torsions to each segment and
3. transfering bond target atoms from the target segment back to the
   origin segment (since they do not move during bond rotation).

11/96 te

///////////////////////////////////////////////////////////////////// */

void initialize_segments
(
  LABEL_FLEX	*label_flex,
  MOLECULE	*molecule
)
{
  int i, j, k;
  int atom_id;			/* Atom id */
  int origin;			/* Origin atom id */
  int target;			/* Target atom id */
  int segment;			/* Current segment id */
  int anchor_atom;		/* Anchor atom of molecule */

  static SEARCH search = {0};
  static MOLECULE temporary = {0};

/*
* Quit this routine if only one segment (no rotatable bonds
* 3/98 te
*/
  if (molecule->total.segments == 1)
    return;

/*
* Reset atom segment records and neighbor flags
* 12/96 te
*/
  for (i = 0; i < molecule->total.segments; i++)
    for (j = 0; j < molecule->segment[i].atom_total; j++)
    {
      atom_id = molecule->segment[i].atom[j];
      molecule->atom[atom_id].segment_id = i;
    }

/*
* Allocate extra space in temporary structure and copy from molecule
* 12/96 te
*/
  temporary.total.segments = molecule->total.segments;
  reallocate_segments (&temporary);

  for (i = 0; i < molecule->total.segments; i++)
  {
    temporary.segment[i].atom_total = molecule->total.atoms;
    reallocate_segment_atoms (&temporary.segment[i]);
  }

  copy_atoms (&temporary, molecule);
  copy_bonds (&temporary, molecule);
  copy_segments (&temporary, molecule);

/*
* Sort atoms radially from anchor atom
* 1/97 te
*/
  if (molecule->transform.anchor_atom != NEITHER)
    anchor_atom = molecule->transform.anchor_atom;

  else
    anchor_atom = 0;

  for
  (
    i = 0;
    breadth_search
    (
      &search,
      molecule->atom,
      molecule->total.atoms,
      get_atom_neighbor,
      flag_atom_neighbor,
      &anchor_atom, 1,
      NEITHER, i
    ) != EOF;
    i++
  );

/*
* Assign each flexible bond to segments and transfer target atom
* over to origin segment, since it does not move upon bond rotation.
* 11/96 te
*/
  for (i = 0; i < molecule->total.torsions; i++)
  {
    if (molecule->torsion[i].flex_id)
    {
/*

      fprintf (global.outfile, "tors %d orig %s (%d) targ %s (%d)\n",
        i, molecule->atom[molecule->torsion[i].origin].name,
        molecule->atom[molecule->torsion[i].origin].segment_id,
        molecule->atom[molecule->torsion[i].target].name,
        molecule->atom[molecule->torsion[i].target].segment_id);
*/

      origin = molecule->torsion[i].origin;
      target = molecule->torsion[i].target;

      if (get_search_radius (&search, origin, NEITHER) >
        get_search_radius (&search, target, NEITHER))
        reverse_torsion (&molecule->torsion[i]);

      origin = molecule->torsion[i].origin;
      target = molecule->torsion[i].target;

/*
*     Update segment torsion records
*     3/97 te
*/
      segment = temporary.atom[target].segment_id;
      molecule->torsion[i].segment_id = segment;

      if (temporary.segment[segment].torsion_id == NEITHER)
        temporary.segment[segment].torsion_id = i;

      else
        exit (fprintf (global.outfile,
          "ERROR initialize_segments: torsion assigned twice to segment\n"));

      temporary.segment[segment].conform_total =
        label_flex->member[molecule->torsion[i].flex_id].torsion_total;

/*
*     Update segment atom records
*     3/97 te
*/
      if (temporary.atom[origin].segment_id !=
        temporary.atom[target].segment_id)
      {
/*
*       Delete target atom from target segment
*       11/96 te
*/
        for (j = 0; j < temporary.segment[segment].atom_total; j++)
          if (temporary.segment[segment].atom[j] == target)
          {
            for (k = j + 1; k < temporary.segment[segment].atom_total; k++)
              temporary.segment[segment].atom[k - 1] = 
                temporary.segment[segment].atom[k];

            temporary.segment[segment].atom_total--;

            if (temporary.atom[target].heavy_flag)
              temporary.segment[segment].heavy_total--;

            break;
          }

/*
*       Add target atom to origin segment
*       11/96 te
*/
        segment = temporary.atom[origin].segment_id;
        atom_id = temporary.segment[segment].atom_total++;
        temporary.segment[segment].atom[atom_id] = target;
        molecule->atom[target].segment_id = segment;

        if (temporary.atom[target].heavy_flag)
          temporary.segment[segment].heavy_total++;
      }

      else
        exit (fprintf (global.outfile,
          "ERROR initialize_segments: torsion atoms in same segment\n"));
    }
  }

  copy_segments (molecule, &temporary);

  if (label_flex->flag && (global.output_volume == 'v'))
  {
    fprintf (global.outfile, "\nFinal segments\n");

    for (i = 0; i < molecule->total.segments; i++)
    {
      fprintf (global.outfile, "  segment    : %d\n", i + 1);
      fprintf (global.outfile, "    conforms : %d\n",
        molecule->segment[i].conform_total);
      fprintf (global.outfile, "    neighbors:");

      for (j = 0; j < molecule->segment[i].neighbor_total; j++)
        fprintf (global.outfile, " %d",
          molecule->segment[i].neighbor[j].id + 1);

      fprintf (global.outfile, "\n    atoms    :");

      for (j = 0; j < molecule->segment[i].atom_total; j++)
        fprintf (global.outfile, " %s",
          molecule->atom[molecule->segment[i].atom[j]].name);

      fprintf (global.outfile, "\n");
    }
  }
}


/* /////////////////////////////////////////////////////////////////// */

void initialize_anchor_layer
(
  LABEL_FLEX	*label_flex,
  MOLECULE	*molecule
)
{
  int segment;
  int segment_id;
  int segment_anchor;
  int layer;
  static SEARCH search = {0};
  SORT_INT *list = NULL;

/*
* If no anchor search, initialize first layer to contain all segments
* 3/97 te
*/
  if (label_flex->anchor_flag == FALSE)
  {
    molecule->total.layers = 1;
    reallocate_layers (molecule);

    molecule->layer[0].segment_total = molecule->total.segments;
    reallocate_layer_segments (&molecule->layer[0]);

    if (molecule->transform.anchor_atom != NEITHER)
      segment_anchor =
        molecule->atom[molecule->transform.anchor_atom].segment_id;

    else
      segment_anchor = 0;

    molecule->layer[0].segment[0] = segment_anchor;
    molecule->layer[0].conform_total = 1;

    for
    (
      segment = 0, segment_id = 1;
      segment < molecule->total.segments;
      segment++
    )
    {
      if (segment != segment_anchor)
        molecule->layer[0].segment[segment_id++] = segment;

      if (label_flex->drive_flag == TRUE)
      {
        if
        (
          molecule->layer[0].conform_total <
          INT_MAX / molecule->segment[segment].conform_total
        )
          molecule->layer[0].conform_total *=
            molecule->segment[segment].conform_total;

        else
          molecule->layer[0].conform_total = INT_MAX;
      }

      molecule->segment[segment].layer_id = 0;
      molecule->segment[segment].active_flag = TRUE;
      molecule->segment[segment].min_flag = TRUE;
    }

    return;
  }

/*
* Otherwise, allocate enough layers to hold one segment each AND
* each of these layers to hold all segments (if necessary).
* 3/97 te
*/
  molecule->total.layers = molecule->total.segments + 1;
  reallocate_layers (molecule);

  for (layer = 0; layer < molecule->total.layers; layer++)
  {
    molecule->layer[layer].segment_total = molecule->total.segments;
    reallocate_layer_segments (&molecule->layer[layer]);
  }

  molecule->total.layers = 1;

/*
* Assign anchor segment to first layer and initialize layer_id of all segments
* 3/97 te
*/
  segment_anchor = molecule->atom[molecule->transform.anchor_atom].segment_id;
  molecule->layer[0].segment_total = 1;
  molecule->layer[0].segment[0] = segment_anchor;
  molecule->layer[0].conform_total = 1;

  molecule->segment[segment_anchor].layer_id = 0;
  molecule->segment[segment_anchor].active_flag = TRUE;
  molecule->segment[segment_anchor].min_flag = TRUE;

  for (segment = 0; segment < molecule->total.segments; segment++)
    if (segment != segment_anchor)
      molecule->segment[segment].layer_id = NEITHER;

/*
* Flag the segment neighbors radially from the anchor segment
* 3/97 te

  NOTE: This section of the code has been disabled.
        It allows additional segments to be added to the anchor layer
        to satisfy the anchor_size parameter.  Since this adds flexibility
        to the anchor, I have decided to temporarily disable it until
        it can be proven to be necessary.
        3/97 te

  for
  (
    segment = 0;
    breadth_search
    (
      &search,
      molecule->segment,
      molecule->total.segments,
      get_segment_neighbor,
      flag_segment_neighbor,
      &segment_anchor, 1,
      NEITHER, segment
    ) != EOF;
    segment++
  );


* Build the anchor layer from anchor segment until it is sufficiently large.
* Merge the largest adjacent segment at each iteration.
* 3/97 te

  for
  (
    anchor_size =
      molecule->segment[segment_anchor].heavy_total;
    anchor_size < label_flex->anchor_size;
  )
  {
    for
    (
      segment_id = 0, segment_anchor = NEITHER, heavy_max = INT_MIN;
      segment_id < molecule->layer[0].segment_total;
      segment_id++
    )
    {
      segment = molecule->layer[0].segment[segment_id];

      for
      (
        neighbor_id = 0;
        neighbor_id < molecule->segment[segment].neighbor_total;
        neighbor_id++
      )
      {
        if (molecule->segment[segment].neighbor[neighbor_id].out_flag == TRUE)
        {
          neighbor = molecule->segment[segment].neighbor[neighbor_id].id;

          if
          (
            (molecule->segment[neighbor].layer_id == NEITHER) &&
            (molecule->segment[neighbor].heavy_total > heavy_max)
          )
          {
            heavy_max = molecule->segment[neighbor].heavy_total;
            segment_anchor = neighbor;
          }
        }
      }
    }

    if (segment_anchor != NEITHER)
    {
      molecule->layer[0].segment[molecule->layer[0].segment_total++] =
        segment_anchor;
      molecule->layer[0].conform_total *=
        molecule->segment[segment_anchor].conform_total;
      molecule->segment[segment_anchor].layer_id = 0;
      molecule->segment[segment_anchor].active_flag = TRUE;
      molecule->segment[segment_anchor].min_flag = TRUE;
      anchor_size += molecule->segment[segment_anchor].heavy_total;
    }

    else
      break;
  }
*/

/*
* Regroup the segments radially from the anchor segment
* 3/97 te
*/
  for
  (
    segment = 0;
    breadth_search
    (
      &search,
      molecule->segment,
      molecule->total.segments,
      get_segment_neighbor,
      flag_segment_neighbor,
      molecule->layer[0].segment, molecule->layer[0].segment_total,
      NEITHER, segment
    ) != EOF;
    segment++
  );

/*
* Assign segments to outer layers
* 3/97 te
*/
  ecalloc
  (
    (void **) &list,
    molecule->total.segments,
    sizeof (SORT_INT),
    "initial layer sort list",
    global.outfile
  );

  molecule->total.layers = search.radius + 1;

  for
  (
    layer = 1;
    layer <= search.radius;
    layer++
  )
  {
    molecule->layer[layer].segment_total = search.total[layer];
    molecule->layer[layer].conform_total = 1;

    for
    (
      segment_id = 0;
      segment_id < search.total[layer];
      segment_id++
    )
    {
      segment = search.target[layer][segment_id];

      list[segment_id].id = segment;
      list[segment_id].value = molecule->segment[segment].heavy_total;
    }

    qsort
    (
      list,
      search.total[layer],
      sizeof (SORT_INT),
      (int (*)()) compare_int_descend
    );

    for
    (
      segment_id = 0;
      segment_id < search.total[layer];
      segment_id++
    )
    {
      segment = list[segment_id].id;
      molecule->segment[segment].layer_id = layer;
      molecule->layer[layer].segment[segment_id] = segment;
      molecule->layer[layer].conform_total *=
        molecule->segment[segment].conform_total;
    }
  }

  efree ((void **) &list);

  if (global.output_volume == 'v')
  {
    fprintf (global.outfile, "Molecule Layers\n");
    for (layer = 0; layer < molecule->total.layers; layer++)
    {
      fprintf (global.outfile, "layer %2d: conforms %4d, segments",
        layer + 1, molecule->layer[layer].conform_total);

      for
      (
        segment = 0;
        segment < molecule->layer[layer].segment_total;
        segment++
      )
        fprintf (global.outfile, " %2d ",
          molecule->layer[layer].segment[segment] + 1);

      fprintf (global.outfile, "\n");
    }
  }
}

