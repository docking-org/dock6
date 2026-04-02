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
#include "label_chem.h"

/* ==================================================================== */
void get_chemical_labels (LABEL_CHEMICAL *label_chemical)
{
  int i, definition_total, definition_count;
  int values_read;
  STRING100 line;
  FILE *chemical_file;

  label_chemical->init_flag = TRUE;

  chemical_file = efopen (label_chemical->file_name, "r", global.outfile);

/*
* Count up the number of chemical label declarations and definitions,
* then allocate memory
* 6/95 te
*/
  label_chemical->total = 1;
  definition_count = 0;
  while (fgets (line, 100, chemical_file) != NULL)
  {
    if (!strncmp (line, "name", 4))
      label_chemical->total++;
    if (!strncmp (line, "definition", 10))
      definition_count++;
  }
  rewind (chemical_file);
  definition_total = definition_count;

  ecalloc
  (
    (void **) &label_chemical->member,
    label_chemical->total,
    sizeof (CHEMICAL_MEMBER),
    "chemical labels",
    global.outfile
  );

  ecalloc
  (
    (void **) &label_chemical->definition,
    definition_count,
    sizeof (NODE),
    "chemical label definitions",
    global.outfile
  );

/*
* Read in the atom label definitions
* 6/95 te
*/
  strcpy (label_chemical->member[0].name, "null");
  label_chemical->member[0].definition_total = 0;
  label_chemical->member[0].definition = &label_chemical->definition[0];
  strcpy (label_chemical->member[0].definition[0].type, "*");
  label_chemical->member[0].radius = label_chemical->member[0].tolerance = 0;

  label_chemical->total = 1;
  definition_count = definition_total = 0;
  while (fgets (line, 100, chemical_file) != NULL)
  {
/*
*   Process chemical label declaration
*   6/95 te
*/
    if (!strncmp (line, "name", 4))
    {
      label_chemical->member[label_chemical->total - 1].definition_total =
        definition_total;
      label_chemical->member[label_chemical->total].definition =
        &label_chemical->definition[definition_count];
      definition_total = 0;

      if (sscanf
        (line, "%*s %s", label_chemical->member[label_chemical->total].name) <
        1)
        exit (fprintf (global.outfile,
          "ERROR get_chemical_labels: Incomplete label declaration in %s\n",
          label_chemical->file_name));
/*
*     Convert label_chemical->member name to lowercase
*     6/95 te
*/
      for (i = 0; i < strlen
        (label_chemical->member[label_chemical->total].name); i++)
        label_chemical->member[label_chemical->total].name[i] =
          (char) tolower
            (label_chemical->member[label_chemical->total].name[i]);

/*
*     Make sure that "null" label was not in input
*     6/95 te
*/
      if (!strncmp
        (label_chemical->member[label_chemical->total].name, "null", 7))
        exit (fprintf (global.outfile,
          "ERROR get_chemical_labels: <null> specified if %s\n",
          label_chemical->file_name));
 
      label_chemical->total++;
    }

/*
*   Process label_chemical->member radius
*   10/95 te
*/
    else if (!strncmp (line, "radius", 5))
    {
      values_read = sscanf (line, "%*s %f %f",
        &label_chemical->member[label_chemical->total - 1].radius,
        &label_chemical->member[label_chemical->total - 1].tolerance);

      if (values_read < 1)
        exit (fprintf (global.outfile,
          "ERROR get_chemical_labels: Incomplete radius specification in %s\n",
          label_chemical->file_name));

      if (values_read < 2)
        label_chemical->member[label_chemical->total - 1].tolerance = 0.0;
    }

/*
*   Process label_chemical->member definition
*   6/95 te
*/
    else if (!strncmp (line, "definition", 10))
    {
      strtok (white_line (line), " ");

      if (!assign_node (&label_chemical->definition[definition_count],
        TRUE))
        exit (fprintf (global.outfile,
          "ERROR get_chemical_labels: Improper label definition in %s\n",
          label_chemical->file_name));

      label_chemical->definition[definition_count].weight = 1.0;
      definition_count++;
      definition_total++;
    }

    else if (!strncmp (line, "weight", 6))
      if (sscanf (line, "%*s %f",
        &label_chemical->definition[definition_count - 1].weight) < 1)
        label_chemical->definition[definition_count - 1].weight = 1.0;
  }

/*
* Update last label_chemical->member info also
* 6/95 te
*/
  label_chemical->member[label_chemical->total - 1].definition_total =
    definition_total;

  efclose (&chemical_file);


/*
* Print out the label_chemical->members and their definitions
* 6/95 te

  if (global.output_volume == 'v')
  {
    fprintf (global.outfile, "\n____Chemical_Label_Definitions____\n\n");
    for (i = 0; i < label_chemical->total; i++)
    {
      MFPUTC ('_', 40);
      fprintf (global.outfile, "\n");
      fprintf (global.outfile, "%-20s%s\n", "name",
        label_chemical->member[i].name);
      fprintf (global.outfile, "%-20s%-7.2f", "radius",
        label_chemical->member[i].radius);
      if (label_chemical->member[i].tolerance)
        fprintf (global.outfile, "%-7.2f", label_chemical->member[i].tolerance);
      fprintf (global.outfile, "\n");

      for (j = 0; j < label_chemical->member[i].definition_total; j++)
      {
        fprintf (global.outfile, "%-20s", "definition");
        print_node (&label_chemical->member[i].definition[j], 0);
        fprintf (global.outfile, "\n");
        fprintf (global.outfile, "%-20s%-7.2f\n", "weight",
            label_chemical->member[i].definition[j].weight);
        fprintf (global.outfile, "\n");
      }
      fprintf (global.outfile, "\n");
    }
    fprintf (global.outfile, "\n\n");
  }
*/
}


/* ////////////////////////////////////////////////////////////////////// */

void free_chemical_labels (LABEL_CHEMICAL *label_chemical)
{
  int i, j;

  for (i = 0; i < label_chemical->total; i++)
    for (j = 0; j < label_chemical->member[i].definition_total; j++)
      free_node (&label_chemical->member[i].definition[j]);

  efree ((void **) &label_chemical->member);
  efree ((void **) &label_chemical->definition);
}


/* ============================================================= */

void get_table
(
  LABEL_CHEMICAL	*label_chemical,
  FILE_NAME		table_file_name,
  float			***table
)
{
  int i, j;
  int continue_loop;
  STRING20 receptor_axis = "RECEPTOR", value;
  FILE *table_file;
  STRING100 line;
  STRING20 *table_label = NULL;
  int label_count, label_match, *label_conversion = NULL;
  char *token, *token_arg;

  if (!label_chemical->init_flag)
    get_chemical_labels (label_chemical);

/*
* Allocate memory (and set elements to zero) for chemical matching table
* 6/95 te
*/
  ecalloc
  (
    (void **) table,
    label_chemical->total,
    sizeof (float *),
    "chemical table",
    global.outfile
  );

  for (i = 0; i < label_chemical->total; i++)
    ecalloc
    (
      (void **) &(*table)[i],
      label_chemical->total,
      sizeof (float),
      "chemical table",
      global.outfile
    );
 
/*
* Read in interaction table
* 10/95 te
*/
  table_file = efopen (table_file_name, "r", global.outfile);

  for (label_count = 0; fgets (line, 100, table_file);)
    if (!strncmp (line, "label", 5))
      label_count++;

  rewind (table_file);

  emalloc
  (
    (void **) &label_conversion,
    label_count * sizeof (int),
    "label conversion array",
    global.outfile
  );

  ecalloc
  (
    (void **) &table_label,
    label_count,
    sizeof (STRING20),
    "label name array",
    global.outfile
  );

  for (label_count = 0; fgets (line, 100, table_file);)
    if (!strncmp (line, "label", 5))
    {
      if (sscanf (line, "%*s %s", table_label[label_count]))
      {
        for (i = 0, label_match = FALSE; i < label_chemical->total; i++)
        {
          if (!strcmp
            (table_label[label_count], label_chemical->member[i].name))
          {
            label_conversion[label_count] = i;
            label_match = TRUE;
          }
        }

        if (label_match)
          label_count++;

        else
          exit (fprintf (global.outfile,
            "ERROR get_table: %s found in %s, but not defined in %s\n",
            table_label[label_count], table_file_name,
            label_chemical->file_name));
      }

      else
        exit (fprintf (global.outfile,
          "ERROR get_table: Empty line encountered in %s\n",
          table_file_name));
    }

  rewind (table_file);

  while (strncmp (line, "table", 5))
    if (!fgets (line, 100, table_file))
      exit (fprintf (global.outfile,
        "ERROR get_table: No table found in %s\n", table_file_name));

  for (i = 0; i < label_count; i++)
  {
    if (!fgets (line, 100, table_file))
      exit (fprintf (global.outfile,
        "ERROR get_table: Table empty in %s\n", table_file_name));

/*
*   Replace tab characters in input line with spaces
*   6/95 te
*/
    white_line (line);

    for (j = 0, token_arg = line; j <= i; j++, token_arg = NULL)
    {
      if (token = strtok (token_arg, " "))
        (*table)[label_conversion[i]][label_conversion[j]] =
          (*table)[label_conversion[j]][label_conversion[i]] =
            atof (token);

      else
        exit (fprintf (global.outfile,
          "ERROR get_table: Insufficient table entries in %s\n",
          table_file_name));
    }
  }

  fclose (table_file);

  efree ((void *) &label_conversion);
  efree ((void *) &table_label);

/*
* Print out match table 
* 6/95 te
*/

  fprintf (global.outfile, "\n____________Chemical_Table____________\n\n");


  fprintf (global.outfile, "   _");
  for (i = 0; i < label_chemical->total; i++)
    fprintf (global.outfile, "%6s", "______");
  fprintf (global.outfile, "_\n");

  fprintf (global.outfile, "  | ");
  for (i = 0; i < label_chemical->total; i++)
    fprintf (global.outfile, "%6s", "");
  fprintf (global.outfile, " |  LIGAND\n");

  for (i = 0; i < label_chemical->total; i++)
  {
    fprintf (global.outfile, "  | ");
    for (j = 0; j < label_chemical->total; j++)
    {
      sprintf (value, "%6.2f", (*table)[i][j]);

      if (value[5] == '0')
      {
        value[5] = ' ';
        if (value[4] == '0')
        {
          value[4] = ' ';
          value[3] = ' ';

          if (atoi (value) == 0)
            strcpy (value, "");
        }
      }

      fprintf (global.outfile, "%6s", value);
    }

    fprintf (global.outfile, " |  %s\n", label_chemical->member[i].name);
  }

  fprintf (global.outfile, "  |_");
  for (i = 0; i < label_chemical->total; i++)
    fprintf (global.outfile, "%6s", "______");
  fprintf (global.outfile, "_|\n\n");

  for (i = 0, continue_loop = TRUE;
    (i < sizeof (STRING20)) && continue_loop; i++)
  {
    if (receptor_axis[i])
    {
      fprintf (global.outfile, "  %c", receptor_axis[i]);
      continue_loop = TRUE;
    }
    else
      fprintf (global.outfile, "%3s", "");

    for (j = continue_loop = 0; j < label_chemical->total; j++)
    {
      if (label_chemical->member[j].name[i])
      {
        fprintf (global.outfile, "   %c  ", label_chemical->member[j].name[i]);
        continue_loop = TRUE;
      }
      else
        fprintf (global.outfile, "%6s", "");
    }
    fprintf (global.outfile, "\n");
  }

  fprintf (global.outfile, "\n\n");
}


/* ============================================================= */

void free_table
(
  LABEL_CHEMICAL	*label_chemical,
  float			***table
)
{
  int i;

  for (i = 0; i < label_chemical->total; i++)
    efree ((void **) &(*table)[i]);
 
  efree ((void **) table);
}


/* //////////////////////////////////////////////////////////////////// */

int assign_chemical_labels
(
  LABEL_CHEMICAL *label_chemical,
  MOLECULE *molecule
)
{
  int i, j, k;

  if (!label_chemical->init_flag)
    get_chemical_labels (label_chemical);

  molecule->info.assign_chem = TRUE;

  for (i = 0; i < molecule->total.atoms; i++)
    molecule->atom[i].chem_id = 0;

  for (i = 0; i < molecule->total.atoms; i++)
    for (j = 1; j < label_chemical->total; j++)
      for (k = 0; k < label_chemical->member[j].definition_total; k++)
        if (check_atom
          (molecule, i, &label_chemical->member[j].definition[k]))
          molecule->atom[i].chem_id = j;

/*
* Print out chemical label assignments
* 10/95 te
*/
  if (global.output_volume == 'v')
  {
    fprintf (global.outfile, "____Chemical_Assignments____\n");
    fprintf (global.outfile, "%s:\n\n", molecule->info.name);
    for (i = 0; i < molecule->total.atoms; i++)
      fprintf (global.outfile, "%6d %4s %5s   %-20s\n", i + 1,
        molecule->atom[i].name,
        molecule->atom[i].type,
        label_chemical->member[molecule->atom[i].chem_id].name);
    fprintf (global.outfile, "\n\n");
  }

  return TRUE;
}
 

/* //////////////////////////////////////////////////////////////// */

void read_chemical_labels
(
  LABEL_CHEMICAL *label_chemical,
  MOLECULE *molecule,
  FILE_NAME molecule_file_name,
  FILE *molecule_file
)
{
  int i, j;
  long file_position;
  int label_count, label_match, *label_conversion = NULL;
  char buff[200];

/*
* Record current file position
* 2/96 te
*/
  file_position = ftell (molecule_file);
  rewind (molecule_file);

/*
* Read in header line of file
* 6/95 te
*/
  if (!fgets (buff, 199, molecule_file))
    exit (fprintf (global.outfile,
      "ERROR read_chemical_labels: No data in SPH file %s\n",
      molecule_file_name));

  if (!strncmp (buff, "cluster", 7))
    exit (fprintf (global.outfile,
      "ERROR read_chemical_labels: No chemical labels are in SPH file %s\n",
      molecule_file_name));

/*
* Construct list to convert chemical labels read here to proper
* chem_id values
* 6/95 te
*/
  ecalloc
  (
    (void **) &label_conversion,
    label_chemical->total,
    sizeof (int),
    "chemical label conversion list",
    global.outfile
  );

  for (i = 1; i < label_chemical->total; i++)
    label_conversion[i] = -1;

  if (!fgets (buff, 199, molecule_file))
    exit (fprintf (global.outfile,
      "ERROR read_chemical_labels: No chemical labels are in SPH file %s\n",
      molecule_file_name));

  label_count = 0;
  while (strncmp (buff, "cluster", 7))
  {
    label_match = FALSE;

    for (i = 0; i < label_chemical->total; i++)
      if (!strncmp
        (buff, label_chemical->member[i].name,
        strlen (label_chemical->member[i].name)))
      {
        label_conversion[i] = ++label_count;
        label_match = TRUE;
        break;
      }

    if (!label_match)
      exit (fprintf (global.outfile,
        "ERROR read_chemical_labels: No definition for %s in SPH file %s\n",
        buff, molecule_file_name));

    if (!fgets (buff, 199, molecule_file))
      exit (fprintf (global.outfile,
        "ERROR read_chemical_labels: Insufficient labels in SPH file %s\n",
        molecule_file_name));
  }

  for (i = 0; i < molecule->total.atoms; i++)
  {
/*
*   Convert chemical label integer
*   6/95 te
*/
    if (label_chemical->total && label_conversion)
    {
      if (molecule->atom[i].chem_id > label_count + 1)
        exit (fprintf (global.outfile,
          "ERROR read_chemical_labels: "
          "Improper chemical label for entry %s in %s\n",
          molecule->atom[i].name, molecule_file_name));

      for (j = 0; j < label_chemical->total; j++)
        if (label_conversion[j] == molecule->atom[i].chem_id)
        {
          molecule->atom[i].chem_id = j;
          break;
        }
    }
  }

  efree ((void **) &label_conversion);

/*
* Return to original file position
* 2/96 te
*/
  fseek (molecule_file, file_position, SEEK_SET);
}


