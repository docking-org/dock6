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
#include "io_mol2.h"

int read_mol2
(
  MOLECULE *molecule,
  FILE_NAME molecule_file_name,
  FILE *molecule_file
)
{
  int i, j;			/* Iteration variables */
  int token_error = FALSE;	/* Flag for error reading tokens */
  char *line = NULL;		/* Line of data */
  char *token;			/* Pointer to identify each token */
  int return_value = TRUE;	/* Value to return upon completion */

/*
* Find the next molecule data record
* 6/96 te
*/
  if (find_record (&line, "@<TRIPOS>MOLECULE", molecule_file) == NULL)
  {
    return_value = EOF;
    goto terminate;
  }

/*
* Record molecule name.  If no name, then record default name.
* 10/96 te
*/
  if (vfgets (&line, molecule_file) == NULL)
  {
    return_value = EOF;
    goto terminate;
  }

  for (token = strtok (white_line (line), " "); token;
    token = strtok (NULL, " "))
  {
    if (molecule->info.name)
      vstrcat (&molecule->info.name, " ");
    vstrcat (&molecule->info.name, token);
  }
    
  if (!molecule->info.name || !strcmp (molecule->info.name, ""))
    vstrcpy (&molecule->info.name, "****");

/*
* Record the total molecule components and types
* 10/96 te
*/
  if (vfgets (&line, molecule_file) == NULL)
  {
    return_value = EOF;
    goto terminate;
  }

  if (token = strtok (white_line (line), " "))
    molecule->total.atoms = atoi (token);

  else
  {
    fprintf (global.outfile,
      "WARNING read_mol2: incorrect MOLECULE record for %s in %s\n",
      molecule->info.name, molecule_file_name);
    return_value = FALSE;
    goto terminate;
  }

  if (token = strtok (NULL, " "))
  {
    molecule->total.bonds = atoi (token);

    if (token = strtok (NULL, " "))
    {
      molecule->total.substs = atoi (token);

      token = strtok (NULL, " ");
        /* Skip number features data record */

      if (token = strtok (NULL, " "))
        molecule->total.sets = atoi (token);

      else
        molecule->total.sets = 0;
    }

    else
    {
      molecule->total.substs = 0;
      molecule->total.sets = 0;
    }
  }

  else
  {
    molecule->total.bonds = 0;
    molecule->total.substs = 0;
    molecule->total.sets = 0;
  }

  if (vfgets (&line, molecule_file) == NULL)
  {
    return_value = EOF;
    goto terminate;
  }

  vstrcpy (&molecule->info.molecule_type, strtok (white_line (line), " "));

  if (vfgets (&line, molecule_file) == NULL)
  {
    return_value = EOF;
    goto terminate;
  }

  vstrcpy (&molecule->info.charge_type, strtok (white_line (line), " "));

  if (vfgets (&line, molecule_file) == NULL)
  {
    return_value = EOF;
    goto terminate;
  }

  if (strncmp (line, "@<TRIPOS>", 9))
  {
    vstrcpy (&molecule->info.status_bits, strtok (white_line (line), " "));

/*
*   Record molecule comment.  If no comment was read,
*   then record default comment.
*   10/96 te
*/
    if (vfgets (&line, molecule_file) == NULL)
  {
    return_value = EOF;
    goto terminate;
  }

    for (token = strtok (white_line (line), " "); token;
      token = strtok (NULL, " "))
    {
      if (molecule->info.comment)
        vstrcat (&molecule->info.comment, " ");
      vstrcat (&molecule->info.comment, token);
    }
    
    if (!molecule->info.comment || !strcmp (molecule->info.comment, ""))
      vstrcpy (&molecule->info.comment, "****");
  }

  else
  {
    vstrcpy (&molecule->info.status_bits, "****");
    vstrcpy (&molecule->info.comment, "****");
    fseek (molecule_file, (long) -strlen (line), SEEK_CUR);
  }

/*
* Allocate space for molecule components
* 2/96 te
*/
  reallocate_atoms (molecule);
  reallocate_bonds (molecule);
  reallocate_substs (molecule);
  reallocate_sets (molecule);

/*
* Read in atoms
* 6/95 te
*/
  if (find_record (&line, "@<TRIPOS>ATOM", molecule_file) == NULL)
  {
    fprintf (global.outfile,
      "WARNING read_mol2: Unable to find ATOM record for %s in %s\n",
      molecule->info.name, molecule_file_name);
    return_value = FALSE;
    goto terminate;
  }

  for (i = 0; i < molecule->total.atoms; i++)
  {
    if (vfgets (&line, molecule_file) == NULL)
    {
      fprintf (global.outfile,
        "WARNING read_mol2: Incomplete ATOM record for %s in %s\n",
        molecule->info.name, molecule_file_name);
      return_value = FALSE;
      goto terminate;
    }

/*
*   Parse each data line into space-separated tokens and check that each
*   token is found.
*   9/95 te
*/
    molecule->atom[i].number = i + 1;

    if (token = strtok (white_line (line), " "));
      /* Skip the atom id field */
    else
      token_error = TRUE;

    if (!token_error && (token = strtok (NULL, " ")))
      vstrcpy (&molecule->atom[i].name, token);
    else
      token_error = TRUE;

    if (!token_error && (token = strtok (NULL, " ")))
      molecule->coord[i][0] = atof (token);
    else
      token_error = TRUE;

    if (!token_error && (token = strtok (NULL, " ")))
      molecule->coord[i][1] = atof (token);
    else
      token_error = TRUE;

    if (!token_error && (token = strtok (NULL, " ")))
      molecule->coord[i][2] = atof (token);
    else
      token_error = TRUE;

    if (!token_error && (token = strtok (NULL, " ")))
      vstrcpy (&molecule->atom[i].type, token);
    else
      token_error = TRUE;

    if (!token_error && (token = strtok (NULL, " ")))
      molecule->atom[i].subst_id = atoi (token) - 1;
    else
      token_error = TRUE;

    if ((molecule->atom[i].subst_id < 0) ||
      (molecule->atom[i].subst_id >= molecule->total.substs))
    {
      if (molecule->total.substs == 1)
        molecule->atom[i].subst_id = 0;

      else
      {
        fprintf (global.outfile,
          "WARNING read_mol2: "
          "Improper ATOM record substructure id for %s in %s.\n",
          molecule->info.name, molecule_file_name);
        return_value = FALSE;
        goto terminate;
      }
    }

    if (!token_error && (token = strtok (NULL, " ")));
      /* Skip the substructure name field */
    else
      token_error = TRUE;

    if (!token_error && (token = strtok (NULL, " ")))
      molecule->atom[i].charge = atof (token);
    else
      token_error = TRUE;

    if (token_error)
    {
      fprintf (global.outfile,
        "WARNING read_mol2: Error reading ATOM record for %s in %s.\n",
        molecule->info.name, molecule_file_name);
      return_value = FALSE;
      goto terminate;
    }
  }

/*
* Read in bonds
* 6/95 te
*/
  if (molecule->total.bonds > 0)
  {
    if (find_record (&line, "@<TRIPOS>BOND", molecule_file) == NULL)
    {
      fprintf (global.outfile,
        "WARNING read_mol2: Unable to fine BOND record for %s in %s\n",
        molecule->info.name, molecule_file_name);
      return_value = FALSE;
      goto terminate;
    }

    for (i = 0; i < molecule->total.bonds; i++)
    {
      if (vfgets (&line, molecule_file) == NULL)
      {
        fprintf (global.outfile,
          "WARNING read_mol2: Incomplete BOND record for %s in %s\n",
          molecule->info.name, molecule_file_name);
        return_value = FALSE;
        goto terminate;
      }

/*
*     Parse each data line into space-separated tokens and check that each
*     token is found.
*     9/95 te
*/
      if (token = strtok (white_line (line), " "));
        /* Skip the bond id field */
      else
        token_error = TRUE;

      if (!token_error && (token = strtok (NULL, " ")))
        molecule->bond[i].origin = atoi (token) - 1;
      else
        token_error = TRUE;

      if (!token_error && (token = strtok (NULL, " ")))
        molecule->bond[i].target = atoi (token) - 1;
      else
        token_error = TRUE;

      if (!token_error && (token = strtok (NULL, " ")))
        vstrcpy (&molecule->bond[i].type, strip_char (token, '\n'));
      else
        token_error = TRUE;

      if (token_error)
      {
        fprintf (global.outfile,
          "WARNING read_mol2: Error in BOND record of %s in %s\n",
          molecule->info.name, molecule_file_name);
        return_value = FALSE;
        goto terminate;
      }

      if ((molecule->bond[i].origin >= molecule->total.atoms) ||
        (molecule->bond[i].origin < 0) ||
        (molecule->bond[i].target >= molecule->total.atoms) ||
        (molecule->bond[i].target < 0))
      {
        fprintf (global.outfile,
          "WARNING read_mol2: Improper BOND origin/target for %s in %s.\n",
          molecule->info.name, molecule_file_name);
        return_value = FALSE;
        goto terminate;
      }
    }
  }

/*
* Read in substructure information
* 9/95 te
*/
  if (molecule->total.substs > 0)
  {
    if (find_record (&line, "@<TRIPOS>SUBSTRUCTURE", molecule_file) == NULL)
    {
      fprintf (global.outfile,
        "WARNING read_mol2: No SUBSTRUCTURE record for %s in %s\n",
        molecule->info.name, molecule_file_name);
      return_value = FALSE;
      goto terminate;
    }

    for (i = 0; i < molecule->total.substs; i++)
    {
      if (vfgets (&line, molecule_file) == NULL)
      {
        fprintf (global.outfile,
          "WARNING read_mol2: Incomplete SUBSTRUCTURE record for %s in %s\n",
          molecule->info.name, molecule_file_name);
        return_value = FALSE;
        goto terminate;
      }
/*
*     Parse each data line into space-separated tokens and check that each
*     token is found.
*     9/95 te
*/
      molecule->subst[i].number = i + 1;
#ifdef DEBUG
      printf( "data line:%s:end data line\n", line );
#endif

      if (token = strtok (white_line (line), " "));
        /* Skip the substructure id field */
      else
        token_error = TRUE;

      if (!token_error && (token = strtok (NULL, " ")))
        vstrcpy (&molecule->subst[i].name, token);
      else
        token_error = TRUE;

      if (!token_error && (token = strtok (NULL, " ")))
        molecule->subst[i].root_atom = atoi (token) - 1;
      else
        token_error = TRUE;

/*
* The rest of the fields in the substructure data line are optional.
* http://www.tripos.com/mol2/mol2_format33.html#1370
* Initialization of the optional fields, notably sub_type, avoids
* problems in read_receptor (io_receptor.c:89) for MOLDEN generated
* mol2 files.
* 01/2006 Scott Brozell
*/
      vstrcpy (&molecule->subst[i].type, "");
      molecule->subst[i].dict_type = -1;
      vstrcpy (&molecule->subst[i].chain, "");
      vstrcpy (&molecule->subst[i].sub_type, "");
      molecule->subst[i].inter_bonds = -1;
      vstrcpy (&molecule->subst[i].status, "");

      if (!token_error && (token = strtok (NULL, " ")))
      {
        vstrcpy (&molecule->subst[i].type, strip_char (token, '\n'));

        if (token = strtok (NULL, " "))
        {
          molecule->subst[i].dict_type = atoi (token) - 1;

          if (token = strtok (NULL, " "))
          {
            vstrcpy (&molecule->subst[i].chain, strip_char (token, '\n'));

            if (token = strtok (NULL, " "))
            {
              vstrcpy (&molecule->subst[i].sub_type, strip_char (token, '\n'));

              if (token = strtok (NULL, " "))
              {
                molecule->subst[i].inter_bonds = atoi (token);

                if (token = strtok (NULL, " "))
                  vstrcpy (&molecule->subst[i].status, strip_char (token, '\n'));
              }
            }
          }
        }
      }

      if (token_error)
      {
        fprintf (global.outfile,
          "WARNING read_mol2: Error in SUBSTRUCTURE record for %s in %s\n",
          molecule->info.name, molecule_file_name);
        return_value = FALSE;
        goto terminate;
      }
    }

/*
*   Update substructure information if more that one substructure is present.
*   (Assume molecule is a macromolecule.)  Extract out the residue number
*   from the substructure name only if the beginning of the name is the
*   same as the sub_type.
*   10/95 te
*/
    if (molecule->total.substs > 1)
    {
      for (i = 0, token_error = TRUE; i < molecule->total.substs; i++)
        if ((strlen (molecule->subst[i].name) <= 3) ||
          !isdigit (molecule->subst[i].name[3]) ||
          (atoi (&molecule->subst[i].name[3]) <= 0))
          token_error = FALSE;
     
      if (token_error)
        for (i = 0; i < molecule->total.substs; i++)
          molecule->subst[i].number = atoi
            (&molecule->subst[i].name[3]);
    }
  }

  else
  {
    molecule->total.substs = 1;
    reallocate_substs (molecule);

    molecule->subst[0].number = 1;
    vstrcpy (&molecule->subst[0].name, "****");
    molecule->subst[0].root_atom = 0;

    for (i = 0; i < molecule->total.atoms; i++)
      molecule->atom[i].subst_id = 0;
  }

/*
* Read in set information
* 10/96 te
*/
  if (molecule->total.sets > 0)
  {
    token_error = FALSE;

    if (find_record (&line, "@<TRIPOS>SET", molecule_file) == NULL)
    {
      fprintf (global.outfile,
        "WARNING read_mol2: No SET record for %s in %s\n",
        molecule->info.name, molecule_file_name);
      return_value = FALSE;
      goto terminate;
    }

    for (i = 0; i < molecule->total.sets; i++)
    {
      if (vfgets (&line, molecule_file) == NULL)
      {
        fprintf (global.outfile,
          "WARNING read_mol2: Incomplete SET record for %s in %s\n",
          molecule->info.name, molecule_file_name);
        return_value = FALSE;
        goto terminate;
      }

      if (token = strtok (white_line (line), " "))
        vstrcpy (&molecule->set[i].name, token);
      else
        token_error = TRUE;

      if (!token_error && (token = strtok (NULL, " ")))
        vstrcpy (&molecule->set[i].type, token);
      else
        token_error = TRUE;

      if (!token_error && (token = strtok (NULL, " ")))
        vstrcpy (&molecule->set[i].obj_type, token);
      else
        token_error = TRUE;

      if (!token_error && (token = strtok (NULL, " ")))
      {
        vstrcpy (&molecule->set[i].sub_type, token);

        if (token = strtok (NULL, " "))
        {
          vstrcpy (&molecule->set[i].status, token);

          while (token = strtok (NULL, " "))
          {
            if (molecule->set[i].comment)
              vstrcat (&molecule->set[i].comment, " ");
            vstrcat (&molecule->set[i].comment, token);
          }
        }
      }

      if (token_error)
      {
        fprintf (global.outfile,
          "WARNING read_mol2: Error in SET record for %s in %s\n",
          molecule->info.name, molecule_file_name);
        return_value = FALSE;
        goto terminate;
      }

      if (vfgets (&line, molecule_file) == NULL)
      {
        fprintf (global.outfile,
          "WARNING read_mol2: Incomplete SET record for %s in %s\n",
          molecule->info.name, molecule_file_name);
        return_value = FALSE;
        goto terminate;
      }

      if (token = strtok (white_line (line), " "))
        molecule->set[i].member_total = atoi (token);
      else
        token_error = TRUE;

      emalloc
      (
        (void **) &molecule->set[i].member,
        molecule->set[i].member_total * sizeof (int),
        "molecule set members",
        global.outfile
      );

      for (j = 0; !token_error && (j < molecule->set[i].member_total); j++)
      {
        if (token = strtok (NULL, " "))
        {
          if (strcmp (token, "\\"))
            molecule->set[i].member[j] = atoi (token) - 1;

          else if ((vfgets (&line, molecule_file) == NULL) ||
            ((token = strtok (white_line (line), " ")) == NULL))
            token_error = TRUE;
        }

        else
          token_error = TRUE;
      }

      if (token_error)
      {
        fprintf (global.outfile,
          "WARNING read_mol2: Error in SET record for %s in %s\n",
          molecule->info.name, molecule_file_name);
        return_value = FALSE;
        goto terminate;
      }
    }
  }

/*
* Upon termination, free up space
* 11/96 te
*/
  terminate:

  efree ((void **) &line);

  return return_value;
}


/* ========================================================== */

int write_mol2
(
  MOLECULE *molecule,
  FILE *molecule_file
)
{
  int i, j, k;			/* Iteration variables */
  STRING20 subst_name;		/* String to construct substructure name */
  int segment;			/* Segment id */

  int *atom_key = NULL;		/* Key for old to new atom ids */
  int *bond_key = NULL;		/* Key for old to new bond ids */
  int *subst_key = NULL;	/* Key for old to new subst ids */

  int atom_total = 0;		/* Total atoms in new set */
  int bond_total = 0; 		/* Total bonds in new set */
  int subst_total = 0; 		/* Total substs in new set */

  long header_position;		/* File position of header information */

/*
* Check if partial structure
* 2/98 te
  partial_flag = FALSE;

  if (molecule_in->total.layers > 0)
    for (i = 0; i < molecule_in->total.atoms; i++)
      if (((segment = molecule->atom[atom].segment_id) != NEITHER) &&
        (molecule->segment[segment].active_flag == FALSE))
      {
        partial_flag == TRUE;
        break;
      }

  if (partial_flag == TRUE)
  {
    ecalloc
    (
      (void **) &molecule,
      1,
      sizeof (MOLECULE),
      "write_mol2 temporary molecule",
      global.outfile
    );

    molecule->total = molecule_in->total;
    reallocate_molecule (molecule);
    reset_molecule (molecule);
  }

  else
    molecule = molecule_in;
*/

  fprintf (molecule_file, "\n");
  fprintf (molecule_file, "@<TRIPOS>MOLECULE\n");
  fprintf (molecule_file, "%s\n",
    molecule->info.name ? molecule->info.name : "****");

  header_position = ftell (molecule_file);

  fprintf (molecule_file, "%-5d %-5d %-5d %-5d %-5d\n",
    molecule->total.atoms,
    molecule->total.bonds,
    molecule->total.substs, 0,
    molecule->total.sets);

  fprintf (molecule_file, "%s\n",
    molecule->info.molecule_type ? molecule->info.molecule_type : "****");
  fprintf (molecule_file, "%s\n",
    molecule->info.charge_type ? molecule->info.charge_type : "****");

  if (molecule->info.status_bits)
  {
    fprintf (molecule_file, "%s\n", molecule->info.status_bits);
    fprintf (molecule_file, "%s\n",
      molecule->info.comment ? molecule->info.comment : "****");
  }

/*
* Print out atoms
* 6/95 te
*/
  fprintf (molecule_file, "@<TRIPOS>ATOM\n");
  for (i = 0; i < molecule->total.atoms; i++)
  {
/*
*   Check if only a sub-structure will be written out
*   2/98 te
*/
    if
    (
      (molecule->total.layers > 0) &&
      (
        ((segment = molecule->atom[i].segment_id) != NEITHER) &&
        (molecule->segment[segment].active_flag == FALSE)
      )
    )
    {
      if (atom_key == NULL)
      {
        ecalloc
        (
          (void **) &atom_key,
          molecule->total.atoms,
          sizeof (int),
          "write_mol2 atom_key array",
          global.outfile
        );

        ecalloc
        (
          (void **) &bond_key,
          molecule->total.bonds,
          sizeof (int),
          "write_mol2 bond_key array",
          global.outfile
        );

        ecalloc
        (
          (void **) &subst_key,
          molecule->total.substs,
          sizeof (int),
          "write_mol2 subst_key array",
          global.outfile
        );

        for (j = 0; j < i; j++)
          atom_key[j] = j;

        atom_total = i;
      }

      atom_key[i] = NEITHER;
      continue;
    }

    if (atom_key != NULL)
      atom_key[i] = atom_total++;

    if ((molecule->total.substs > 0) &&
      (molecule->atom[i].subst_id < molecule->total.substs) &&
      molecule->subst[molecule->atom[i].subst_id].name &&
      strcmp (molecule->subst[molecule->atom[i].subst_id].name, "****"))
      sprintf
        (subst_name, "%s", molecule->subst[molecule->atom[i].subst_id].name);
    else
      sprintf (subst_name, "<%d>", molecule->atom[i].subst_id + 1);

    fprintf (molecule_file, "%-5d %-5s %9.4f %9.4f %9.4f %-5s %3d %4s %8.4f\n",
      atom_key == NULL ? i + 1 : atom_key[i] + 1,
      molecule->atom[i].name ? molecule->atom[i].name : "****",
      molecule->coord[i][0],
      molecule->coord[i][1],
      molecule->coord[i][2],
      molecule->atom[i].type ? molecule->atom[i].type : "****",
      molecule->atom[i].subst_id + 1, subst_name,
      molecule->atom[i].charge);

    fflush (molecule_file);
  }

/*
* Print out bonds
* 6/95 te
*/
  fprintf (molecule_file, "@<TRIPOS>BOND\n");

  if (atom_key == NULL)
    for (i = 0; i < molecule->total.bonds; i++)
    {
      fprintf (molecule_file, "%-5d %-5d %-5d %s\n",
        i + 1, molecule->bond[i].origin + 1, molecule->bond[i].target + 1,
        molecule->bond[i].type ? molecule->bond[i].type : "****");
    }

  else
  {
    for (i = 0; i < molecule->total.bonds; i++)
    {
      if ((atom_key[molecule->bond[i].origin] != NEITHER) &&
        (atom_key[molecule->bond[i].target] != NEITHER))
      {
        bond_key[i] = bond_total++;

        fprintf (molecule_file, "%-5d %-5d %-5d %s\n",
          bond_total,
          atom_key[molecule->bond[i].origin] + 1,
          atom_key[molecule->bond[i].target] + 1,
          molecule->bond[i].type ? molecule->bond[i].type : "****");
      }

      else
        bond_key[i] = NEITHER;
    }
  }

/*
* Print out substructures
* 6/95 te
*/
  if (molecule->total.substs > 0)
  {
    fprintf (molecule_file, "@<TRIPOS>SUBSTRUCTURE\n");
    for (i = 0; i < molecule->total.substs; i++)
    {
      if (atom_key == NULL)
        fprintf (molecule_file, "%-3d %-6s %-5d",
          i + 1, molecule->subst[i].name, molecule->subst[i].root_atom + 1);

      else
      {
        if (atom_key[molecule->subst[i].root_atom] != NEITHER)
          fprintf (molecule_file, "%-3d %-6s %-5d",
            i + 1, molecule->subst[i].name,
            atom_key[molecule->subst[i].root_atom] + 1);

        else
        {
          for (j = 0; j < molecule->total.atoms; j++)
            if ((atom_key[j] != NEITHER) && (molecule->atom[j].subst_id == i))
            {
              fprintf (molecule_file, "%-3d %-6s %-5d",
                i + 1, molecule->subst[i].name,
                atom_key[j] + 1);
              break;
            }

/*
*         Check if no root atoms written for substructure
*         2/98 te
*/
          if (j >= molecule->total.atoms)
          {
            subst_key[i] = NEITHER;
            continue;
          }
        }

        subst_key[i] = subst_total++;
      }

      if (molecule->subst[i].type)
      {
        fprintf (molecule_file, " %-10s", molecule->subst[i].type);
        fprintf (molecule_file, " %-3d", molecule->subst[i].dict_type + 1);
        if (molecule->subst[i].chain)
        {
          fprintf (molecule_file, " %-5s", molecule->subst[i].chain);
          if (molecule->subst[i].sub_type)
          {
            fprintf (molecule_file, " %-5s", molecule->subst[i].sub_type);
            if (molecule->subst[i].inter_bonds)
            {
              fprintf (molecule_file, " %1d", molecule->subst[i].inter_bonds);
              if (molecule->subst[i].status)
                fprintf (molecule_file, " %s", molecule->subst[i].status);
            }
          }
        }
      }

      fprintf (molecule_file, "\n");
    }
  }

/*
* Print out sets
* 10/96 te
*/
  if (molecule->total.sets > 0)
  {
    fprintf (molecule_file, "@<TRIPOS>SET\n");

    for (i = 0; i < molecule->total.sets; i++)
    {
      fprintf (molecule_file, "%s %s %s",
        molecule->set[i].name,
        molecule->set[i].type,
        molecule->set[i].obj_type);

      if (molecule->set[i].sub_type)
      {
        fprintf (molecule_file, " %s", molecule->set[i].sub_type);

        if (molecule->set[i].status)
        {
          fprintf (molecule_file, " %s", molecule->set[i].status);

          if (molecule->set[i].comment)
            fprintf (molecule_file, " %s", molecule->set[i].comment);
        }
      }

      if (atom_key == NULL)
      {
        fprintf (molecule_file, "\n%d", molecule->set[i].member_total);

        for (j = 0; j < molecule->set[i].member_total; j++)
          fprintf (molecule_file, " %d", molecule->set[i].member[j] + 1);
      }

/*
*     Edit set for non-printed members
*     2/98 te
*/
      else
      {
        if (!strcmp (molecule->set[i].type, "ATOMS"))
        {
          for (j = k = 0; j < molecule->set[i].member_total; j++)
            if ((molecule->set[i].member[j] < molecule->total.atoms) &&
              (atom_key[molecule->set[i].member[j]] != NEITHER))
              k++;

          fprintf (molecule_file, "\n%d", k);

          for (j = 0; j < molecule->set[i].member_total; j++)
            if ((molecule->set[i].member[j] < molecule->total.atoms) &&
              (atom_key[molecule->set[i].member[j]] != NEITHER))
              fprintf (molecule_file, " %d",
                atom_key[molecule->set[i].member[j]] + 1);
        }

        if (!strcmp (molecule->set[i].type, "BONDS"))
        {
          for (j = k = 0; j < molecule->set[i].member_total; j++)
            if ((molecule->set[i].member[j] < molecule->total.bonds) &&
              (bond_key[molecule->set[i].member[j]] != NEITHER))
              k++;

          fprintf (molecule_file, "\n%d", k);

          for (j = 0; j < molecule->set[i].member_total; j++)
            if ((molecule->set[i].member[j] < molecule->total.bonds) &&
              (bond_key[molecule->set[i].member[j]] != NEITHER))
              fprintf (molecule_file, " %d",
                bond_key[molecule->set[i].member[j]] + 1);
        }

        if (!strcmp (molecule->set[i].type, "SUBSTS"))
        {
          for (j = k = 0; j < molecule->set[i].member_total; j++)
            if ((molecule->set[i].member[j] < molecule->total.substs) &&
              (subst_key[molecule->set[i].member[j]] != NEITHER))
              k++;

          fprintf (molecule_file, "\n%d", k);

          for (j = 0; j < molecule->set[i].member_total; j++)
            if ((molecule->set[i].member[j] < molecule->total.substs) &&
              (subst_key[molecule->set[i].member[j]] != NEITHER))
              fprintf (molecule_file, " %d",
                subst_key[molecule->set[i].member[j]] + 1);
        }
      }


      fprintf (molecule_file, "\n");
    }
  }

  fprintf (molecule_file, "\n\n");

  if (atom_key != NULL)
  {
    fseek (molecule_file, header_position, SEEK_SET);

    fprintf (molecule_file, "%-5d %-5d %-5d %-5d %-5d\n",
      atom_total,
      bond_total,
      subst_total, 0,
      molecule->total.sets);

    fseek (molecule_file, 0, SEEK_END);

    efree ((void **) &atom_key);
    efree ((void **) &bond_key);
    efree ((void **) &subst_key);
  }

  return TRUE;
}

