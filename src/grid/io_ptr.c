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
#include "io_ptr.h"
#include "vector.h"
#include "transform.h"
#include "rotrans.h"


int read_ptr
(
  MOLECULE *mol_ref,
  MOLECULE *mol_init,
  FILE_NAME in_file_name,
  FILE *in_file,
  int read_source
)
{
  int i, j;
  int read_flag = FALSE;
  char *line = NULL;
  char *token;
  STRING20 field;
  FILE *source_file = NULL;
  static FILE_NAME previous_file_name;
  static FILE *previous_file = NULL;
  MOLECULE temporary = {0};

  reset_molecule (&temporary);

/*
* Check if termination requested
* 5/97 te
*/
  if ((mol_ref == NULL) && (mol_init == NULL))
  {
    if (previous_file != NULL)
      efclose (&previous_file);

    return EOF;
  }

/*
* Read in line of data from file
* 6/96 te
*/
  if (!vfgets (&line, in_file))
    return EOF;

/*
* Process each data item in line, one at a time
* 10/96 te
*/
  for (token = strtok (white_line (line), " "); token;
    token = strtok (NULL, " "))
  {
/*
*   Read in file name containing source molecule
*   10/96 te
*/
    if (!strcmp (token, "<FILE>"))
    {
      if (token = strtok (NULL, " "))
        vstrcpy (&mol_ref->info.source_file, token);

      else
        exit (fprintf (global.outfile,
          "ERROR read_ptr: Missing <FILE> data in %s\n", in_file_name));
    }

/*
*   Get file position
*   10/96 te
*/
    else if (!strcmp (token, "<FPOS>"))
    {
      if (token = strtok (NULL, " "))
        mol_ref->info.source_position = atol (token);

      else
        exit (fprintf (global.outfile,
          "ERROR read_ptr: Missing <FPOS> data in %s\n", in_file_name));
    }

/*
*   If the source file is not to be read (chemical screen run),
*   then be sure to read in the molecule name, description and screen keys
*   11/96 te
*/
    if (!read_source)
    {
/*
*     Read in file name containing source molecule
*     11/96 te
*/
      if (!strcmp (token, "<NAME>"))
      {
        if (token = strtok (NULL, " "))
          vstrcpy (&mol_ref->info.name, token);

        else
          exit (fprintf (global.outfile,
            "ERROR read_ptr: Missing <NAME> data in %s\n", in_file_name));
      }

/*
*     Read in file name containing source molecule
*     11/96 te
*/
      else if (!strcmp (token, "<DESCR>"))
      {
        if (token = strtok (NULL, " "))
          vstrcpy (&mol_ref->info.comment, token);

        else
          exit (fprintf (global.outfile,
            "ERROR read_ptr: Missing <DESCR> data in %s\n",
            in_file_name));
      }

/*
*     Check if chemical key information is present
*     10/96 te
*/
      else if (!strcmp (token, "<KEY>"))
      {
        if (token = strtok (NULL, " "))
          mol_ref->total.keys = atoi (token);

        else
          exit (fprintf (global.outfile,
            "* * * Missing <KEY> data in pointer file.\n"));

        reallocate_keys (mol_ref);

        if (!(token = strtok (NULL, " ")) || strcmp (token, "<KFOLD>"))
          exit (fprintf (global.outfile,
            "ERROR read_ptr: Missing <KFOLD> field in %s\n",
            in_file_name));

        if (token = strtok (NULL, " "))
          mol_ref->transform.fold_flag = atoi (token);

        else
          exit (fprintf (global.outfile,
            "ERROR read_ptr: Missing <KFOLD> data in %s\n",
            in_file_name));

        if (mol_ref->transform.fold_flag != TRUE)
        {
          mol_ref->total.atoms =
            mol_ref->transform.heavy_total =
            mol_ref->total.keys;
          reallocate_atoms (mol_ref);
        }

        for (i = 0; i < mol_ref->total.keys; i++)
        {
          sprintf (field, "<KI%d>", i + 1);

          if (!(token = strtok (NULL, " ")) || strcmp (token, field))
            exit (fprintf (global.outfile,
              "ERROR read_ptr: Missing %s field in %s\n",
              field, in_file_name));

          if (token = strtok (NULL, " "))
          {
            if (mol_ref->transform.fold_flag != TRUE)
              mol_ref->atom[i].chem_id = atoi (token);
          }

          else
            exit (fprintf (global.outfile,
              "ERROR read_ptr: Missing %s data in %s\n",
              field, in_file_name));

          for (j = i; j < mol_ref->total.keys; j++)
          {
            sprintf (field, "<KJ%d>", j + 1);

            if (!(token = strtok (NULL, " ")) || strcmp (token, field))
              exit (fprintf (global.outfile,
                "ERROR read_ptr: Missing %s field in %s\n",
                field, in_file_name));

            if (token = strtok (NULL, " "))
              mol_ref->key[i][j].count = atoi (token);

            else
              exit (fprintf (global.outfile,
                "ERROR read_ptr: Missing %s data in %s\n",
                field, in_file_name));

            if (token = strtok (NULL, " "))
              mol_ref->key[i][j].distance = strtoul (token, NULL, 16);

            else
              exit (fprintf (global.outfile,
                "ERROR read_ptr: Missing %s data in %s\n",
                field, in_file_name));
          }
        }
      }

      continue;
    }

/*
*   Skip other initialization info if no data structure was passed to hold it
*   10/96 te
*/
    if (mol_init == NULL)
      continue;

/*
*   Check if rot/trans information is present
*   10/96 te
*/
    else if (!strcmp (token, "<TRANS>"))
    {
      temporary.transform.flag = TRUE;
      temporary.transform.trans_flag = TRUE;

      for (i = 0; i < 3; i++)
      {
        if (token = strtok (NULL, " "))
          temporary.transform.translate[i] = atof (token);

        else
          exit (fprintf (global.outfile,
            "ERROR read_ptr: Missing <TRANS> data in %s\n",
            in_file_name));
      }
    }

    else if (!strcmp (token, "<ROT>"))
    {
      temporary.transform.flag = TRUE;
      temporary.transform.rot_flag = TRUE;

      for (i = 0; i < 3; i++)
      {
        if (token = strtok (NULL, " "))
          temporary.transform.rotate[i] = atof (token);

        else
          exit (fprintf (global.outfile,
            "ERROR read_ptr: Missing <ROT> data in %s\n",
            in_file_name));
      }
    }

/*
*   Check if torsion information is present
*   10/96 te
*/
    else if (!strcmp (token, "<TORS>"))
    {
      temporary.transform.flag = TRUE;

      if (token = strtok (NULL, " "))
        temporary.total.torsions = atoi (token);

      else
        exit (fprintf (global.outfile,
          "ERROR read_ptr: Missing <TORS> data in %s\n",
          in_file_name));

      if (!(token = strtok (NULL, " ")) || strcmp (token, "<TANCHOR>"))
        exit (fprintf (global.outfile,
          "ERROR read_ptr: Missing <TANCHOR> field in %s\n",
          in_file_name));

      if (token = strtok (NULL, " "))
        temporary.transform.anchor_atom = atoi (token) - 1;

      else
        exit (fprintf (global.outfile,
          "ERROR read_ptr: Missing <TANCHOR> data in %s\n",
          in_file_name));

      if (temporary.total.torsions > 0)
      {
        temporary.transform.tors_flag = TRUE;
        reallocate_torsions (&temporary);

        for (i = 0; i < temporary.total.torsions; i++)
        {
          sprintf (field, "<T%d>", i + 1);

          if (!(token = strtok (NULL, " ")) || strcmp (token, field))
            exit (fprintf (global.outfile,
              "ERROR read_ptr: Missing %s field in %s\n",
              field, in_file_name));

          if (token = strtok (NULL, " "))
            temporary.torsion[i].bond_id = atoi (token) - 1;

          else
            exit (fprintf (global.outfile,
              "ERROR read_ptr: Missing %s data in %s\n",
              field, in_file_name));

          if (token = strtok (NULL, " "))
            temporary.torsion[i].target_angle = atof (token) * PI / 180.0;

          else
            exit (fprintf (global.outfile,
              "ERROR read_ptr: Missing %s data in %s\n",
              field, in_file_name));
        }
      }
    }

/*
*   Check if reflection information is present
*   10/96 te
*/
    else if (!strcmp (token, "<REFL>"))
    {
      if (token = strtok (NULL, " "))
        temporary.transform.refl_flag = atoi (token);

      else
        exit (fprintf (global.outfile,
          "ERROR read_ptr: Missing <REFL> data in %s\n",
          in_file_name));

      if (temporary.transform.refl_flag == TRUE)
        temporary.transform.flag = TRUE;
    }
  }

  efree ((void **) &line);

/*
* Check if molecule file and position were read
* 10/96 te
*/
  if (!strcmp (mol_ref->info.source_file, "") ||
    (mol_ref->info.source_position < 0))
    exit (fprintf (global.outfile,
      "ERROR read_ptr: Missing <FILE> and/or <FPOS> field in %s\n",
      in_file_name));

  if (read_source)
  {
/*
*   Check to see if the current file is still open
*   6/96 te
*/
    if (!strcmp (mol_ref->info.source_file, previous_file_name))
      source_file = previous_file;

    else
    {
      source_file = efopen (mol_ref->info.source_file, "r", global.outfile);
      strcpy (previous_file_name, mol_ref->info.source_file);

      if (previous_file != NULL)
        efclose (&previous_file);

      previous_file = source_file;
    }

/*
*   Read in molecule
*   10/96 te
*/
    fseek (source_file, mol_ref->info.source_position, SEEK_SET);

    read_flag = read_molecule
      (mol_ref, mol_init, mol_ref->info.source_file, source_file, FALSE);

    if ((read_flag != TRUE) || (mol_init == NULL))
      return read_flag;

    if (temporary.transform.flag == TRUE)
    {
      copy_transform (mol_init, &temporary);

      if (temporary.total.torsions > 0)
      {
        copy_torsions (mol_init, &temporary);

        if (mol_init->transform.anchor_atom != NEITHER)
          revise_atom_neighbors (mol_init);

        get_torsion_neighbors (mol_init);
      }

      center_of_mass
      (
        mol_init->coord,
        mol_init->total.atoms,
        mol_init->transform.com
      );

/*
*     Generate initial coordinates based on transformation specified in input
*     7/96 te
*/
      transform_molecule (mol_init, mol_ref);
    }
  }

  else if (mol_init != NULL)
    copy_molecule (mol_init, mol_ref);

  free_molecule (&temporary);
  return TRUE;
}


/* ////////////////////////////////////////////////////////////////// */

int write_ptr
(
  MOLECULE *molecule,
  FILE_NAME in_file_name,
  FILE *out_file
)
{
  int i, j;
  char *word = NULL;

/*
* Output molecule identifiers
* 10/96 te
*/
  fprintf (out_file, "<ID> %d", molecule->info.output_id);

  if (molecule->info.input_id != NEITHER)
    fprintf (out_file, " <SRCID> %d", molecule->info.input_id);

  fprintf (out_file, " <NAME> %s",
    subst_char (vstrcpy (&word, molecule->info.name), ' ', '_'));
  fprintf (out_file, " <DESCR> %s",
    subst_char (vstrcpy (&word, molecule->info.comment), ' ', '_'));

  efree ((void **) &word);

/*
* Output source file name and position
* 6/95 te
*/
  if (molecule->info.source_file)
  {
    fprintf (out_file, " <FILE> %s", molecule->info.source_file);
    fprintf (out_file, " <FPOS> %ld", molecule->info.source_position);
  }

  else
  {
    fprintf (out_file, " <FILE> %s", in_file_name);
    fprintf (out_file, " <FPOS> %ld", molecule->info.file_position);
  }

/*
* Output chemical keys
* 1/97 te
*/
  if (molecule->total.keys > 0)
  {
    fprintf (out_file, " <KEY> %d", molecule->total.keys);
    fprintf (out_file, " <KFOLD> %d", molecule->transform.fold_flag);

/*
*   Write out each chemical key
*   11/96 te
*/
    for (i = 0; i < molecule->total.keys; i++)
    {
      fprintf (out_file, " <KI%d>", i + 1);

      if (molecule->transform.fold_flag == TRUE)
        fprintf (out_file, " %d", i);

      else
        fprintf (out_file, " %d", molecule->atom[i].chem_id);

      for (j = i; j < molecule->total.keys; j++)
      {
        fprintf (out_file, " <KJ%d>", j + 1);
        fprintf (out_file, " %d", molecule->key[i][j].count);
        fprintf (out_file, " %lx", molecule->key[i][j].distance);
      }
    }

    fprintf (out_file, " <END>\n");
    return TRUE;
  }

/*
* If a transformation has occurred, then output rotation/translation info
* 10/96 te
*/
  if (molecule->transform.trans_flag == TRUE)
  {
    fprintf (out_file, " <TRANS>");

    for (i = 0; i < 3; i++)
      fprintf (out_file, " %.6g", molecule->transform.translate[i]);
  }

  if (molecule->transform.rot_flag == TRUE)
  {
    fprintf (out_file, " <ROT>");

    for (i = 0; i < 3; i++)
      fprintf (out_file, " %.6g", molecule->transform.rotate[i]);
  }

/*
* Signal flexible bond rotation
* 10/96 te
*/
  if ((molecule->transform.tors_flag != NEITHER) &&
    (molecule->transform.anchor_atom != NEITHER))
  {
    fprintf (out_file, " <TORS> %d", molecule->total.torsions);
    fprintf (out_file, " <TANCHOR> %d", molecule->transform.anchor_atom + 1);

/*
*   Write out each rotatable bond and angle.
*   10/96 te
*/
    for (i = 0; i < molecule->total.torsions; i++)
    {
      fprintf (out_file, " <T%d>", i + 1);
      fprintf (out_file, " %d", molecule->torsion[i].bond_id + 1);
      fprintf (out_file, " %.6g",
        molecule->torsion[i].current_angle / PI * 180.0);
    }
  }

/*
* Write out reflection state
* 10/96 te
  fprintf (out_file, " <HEAVY> %d", molecule->transform.heavy_total);
*/
  if (molecule->transform.refl_flag != NEITHER)
    fprintf (out_file, " <REFL> %d", molecule->transform.refl_flag);

  fprintf (out_file, " <END>\n");
  return TRUE;
}

