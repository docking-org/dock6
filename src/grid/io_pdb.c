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
#include "io_pdb.h"

typedef struct link_atom
{
  ATOM atom;
  XYZ coord;
  struct link_atom *previous, *next;
} LINK_ATOM;

typedef struct link_bond
{
  BOND bond;
  struct link_bond *previous, *next;
} LINK_BOND;

typedef struct link_subst
{
  SUBST subst;
  struct link_subst *previous, *next;
} LINK_SUBST;


int read_pdb
(
  MOLECULE *molecule,
  FILE_NAME molecule_file_name,
  FILE *molecule_file
)
{
  int i, j, origin, target;
  char buff[200], buff2[200];
  STRING20 temp, subst_name;
  int subst_number;
  int found_atom;
  int new_substructure;
  LINK_ATOM *first_atom, *previous_atom, *current_atom;
  LINK_BOND *first_bond, *previous_bond, *current_bond;
  LINK_SUBST *first_subst, *previous_subst, *current_subst;

  void assign_atom_type (ATOM *);

/*
* Initialize molecule info
* 6/95 te
*/
  vstrcpy (&molecule->info.name, "****");
  vstrcpy (&molecule->info.status_bits, "****");
  vstrcpy (&molecule->info.comment, "****");

/*
* Read in molecule info if found before ATOM records
* 6/95 te
*/
  if (!fgets (buff, 199, molecule_file))
    return EOF;

  while (strncmp (buff, "ATOM", 4) && strncmp (buff, "HETATM", 6))
  {
    if (!strncmp (buff, "HEADER", 6))
    {
      sscanf (buff, "%*s %s", buff2);
      vstrcpy (&molecule->info.name, buff2);
    }

    if (!strncmp (buff, "COMPND", 6))
    {
      for (i = 6; (i < strlen (buff)) && isspace (buff[i]); i++);
      vstrcpy (&molecule->info.comment, strip_newline (&buff[i]));
      if (strchr (molecule->info.comment, ' '))
        strcpy (strchr (molecule->info.comment, ' '), "");
    }

    if (!fgets (buff, 199, molecule_file))
      return EOF;
  }

/*
* Initialize linked lists and counters
* 6/95 te
*/
  first_atom = current_atom = previous_atom = NULL;
  first_subst = current_subst = previous_subst = NULL;
  first_bond = current_bond = previous_bond = NULL;

/*
* Read in atom information and store in a linked list
* 6/95 te
*/
  while (!strncmp (buff, "ATOM", 4) || !strncmp (buff, "HETATM", 6))
  {
/*
*   Allocate memory for this ATOM link
*   6/95 te
*/
    previous_atom = current_atom;
    current_atom = NULL;

    ecalloc
    (
      (void **) &current_atom,
      1,
      sizeof (LINK_ATOM),
      "linked atom list",
      global.outfile
    );

/*
*   Update pointers between ATOM links, initialize first link
*   6/95 te
*/
    current_atom->previous = previous_atom;

    if (previous_atom)
      previous_atom->next = current_atom;

    else
      first_atom = current_atom;

/*
*   Deposit atom info in this ATOM link
*   6/95 te
*/
    memset (temp, 0, sizeof (STRING20));
    current_atom->atom.number = atoi (strncpy (temp, &buff[6], 5));

    memset (buff2, 0, 199);
    strncpy (buff2, &buff[12], 4);
    sscanf (buff2, "%s", temp);
    vstrcpy (&current_atom->atom.name, temp);

    memset (temp, 0, sizeof (STRING20));
    strncpy (subst_name, &buff[17], 3);
    subst_name[3] = 0;

    memset (temp, 0, sizeof (STRING20));
    subst_number = atoi (strncpy (temp, &buff[22], 4));

    memset (temp, 0, sizeof (STRING20));
    current_atom->coord[0] = atof (strncpy (temp, &buff[30], 8));

    memset (temp, 0, sizeof (STRING20));
    current_atom->coord[1] = atof (strncpy (temp, &buff[38], 8));

    memset (temp, 0, sizeof (STRING20));
    current_atom->coord[2] = atof (strncpy (temp, &buff[46], 8));

    if (!strcmp (strrchr (molecule_file_name, '.'), ".xpdb"))
    {
      memset (temp, 0, sizeof (STRING20));
      current_atom->atom.charge = atof (strncpy (temp, &buff[55], 8));

      memset (temp, 0, sizeof (STRING20));
      strncpy (current_atom->atom.type, &buff[71], 5);
    }

    else
    {
      current_atom->atom.charge = 0.0;
      assign_atom_type (&current_atom->atom);
    }

/*
*   Store substructure information if the current residue is new
*   6/95 te
*/
    if (!previous_atom || !current_subst ||
      (subst_number != current_subst->subst.number))
    {
/*
*     Check to see if this substructure has been seen before
*     6/95 te
*/
      for (i = 0, current_subst = first_subst, new_substructure = TRUE;
        current_subst != NULL;
        i++, previous_subst = current_subst,
        current_subst = current_subst->next)
      {
        if (current_subst->subst.number == subst_number)
        {
          new_substructure = FALSE;
          current_atom->atom.subst_id = i;
        }
      }

/*
*     If the current atom is a member of a new residue, then add a new link
*     10/95 te
*/
      if (new_substructure)
      {
        current_subst = NULL;

        ecalloc
        (
          (void **) &current_subst,
          1,
          sizeof (LINK_SUBST),
          "linked subst list",
          global.outfile
        );

/*
*       Update pointers between SUBST links, initialize first link
*       6/95 te
*/
        current_subst->previous = previous_subst;

        if (previous_subst)
          previous_subst->next = current_subst;
        else
          first_subst = current_subst;

/*
*       Deposit substructure info in this SUBST link
*       6/95 te
*/
        current_subst->subst.number = subst_number;
        vstrcpy (&current_subst->subst.name, subst_name);
        current_subst->subst.root_atom = molecule->total.atoms;

        current_atom->atom.subst_id = molecule->total.substs++;
      }

/*
*     Otherwise, set current subst pointer to the last link added
*     10/95 te
*/
      else
        current_subst = previous_subst;
    }

    else
      current_atom->atom.subst_id = molecule->total.substs - 1;

/*
*   Update atom info
*   10/95 te
*/
    molecule->total.atoms++;

/*
*   Read in next line of file
*   6/95 te
*/
    if (!fgets (buff, 199, molecule_file))
      break;
  }

/*
* Advance to TER record, but stop if CONECT records are found
* 6/95 te
*/
  while (strncmp (buff, "TER", 3))
  {
    if (!strncmp (buff, "CONECT", 6))
      break;

    if (!fgets (buff, 199, molecule_file))
      break;
  }

/*
* Read in connectivity information
* 6/95 te
*/
  while (!strncmp (buff, "CONECT", 6))
  {
    memset (temp, 0, sizeof (STRING20));
    origin = atoi (strncpy (temp, &buff[6], 5));

/*
*   Scan the data in the CONECT record
*   6/95 te
*/
    for (i = 1; target = atoi (strncpy (temp, &buff[i * 5 + 6], 5)); i++)
    {
      if (target > origin)
      {
/*
*       Allocate memory for this BOND link
*       6/95 te
*/
        previous_bond = current_bond;
        current_bond = NULL;

        ecalloc
        (
          (void **) &current_bond,
          1,
          sizeof (LINK_BOND),
          "linked bond list",
          global.outfile
        );

/*
*       Update pointers between BOND links, initialize first link
*       6/95 te
*/
        current_bond->previous = previous_bond;

        if (previous_bond)
          previous_bond->next = current_bond;

        else
          first_bond = current_bond;

/*
*       Deposit bond info in this BOND link
*       6/95 te
*/
        current_bond->bond.origin = origin;
        current_bond->bond.target = target;
        vstrcpy (&current_bond->bond.type, "1");
        molecule->total.bonds++;
      }
    }

/*
*   Read in next line of file
*   6/95 te
*/
    if (!fgets (buff, 199, molecule_file))
      break;
  }

/*
* Advance to TER record
* 6/95 te
*/
  while (strncmp (buff, "TER", 3))
    if (!fgets (buff, 199, molecule_file))
      break;

/*
* Allocate space for molecule components
* 2/96 te
*/
  reallocate_atoms (molecule);
  reallocate_bonds (molecule);
  reallocate_substs (molecule);

/*
* Copy atom info into molecule structure
* 6/95 te
*/
  for (i = 0, current_atom = first_atom;
    (i < molecule->total.atoms) && (current_atom != NULL);
    i++, current_atom = current_atom->next)
  {
    copy_atom (&molecule->atom[i], &current_atom->atom);
    copy_coord (molecule->coord[i], current_atom->coord);
  }

  if (i != molecule->total.atoms)
    exit (fprintf (global.outfile,
      "ERROR read_pdb: Error in atom linked list for %s in %s.\n",
      molecule->info.name, molecule_file_name));

/*
* Copy bond info into molecule structure
* 6/95 te
*/
  for (i = 0, current_bond = first_bond;
    (i < molecule->total.bonds) && (current_bond != NULL);
    i++, current_bond = current_bond->next)
  {
    copy_bond (&molecule->bond[i], &current_bond->bond);

/*
*   Convert bond reference from atom number to position in atom array
*   10/95 te
*/
    for (j = 0, found_atom = FALSE; j < molecule->total.atoms; j++)
    {
      if (molecule->bond[i].origin == molecule->atom[j].number)
      {
        molecule->bond[i].origin = j;
        found_atom = TRUE;
      }
    }

    if (!found_atom)
      exit (fprintf (global.outfile,
        "ERROR read_pdb: Unknown atom in CONECT record for %s in %s\n",
        molecule->info.name, molecule_file_name));

    for (j = 0, found_atom = FALSE; j < molecule->total.atoms; j++)
    {
      if (molecule->bond[i].target == molecule->atom[j].number)
      {
        molecule->bond[i].target = j;
        found_atom = TRUE;
      }
    }

    if (!found_atom)
      exit (fprintf (global.outfile,
        "ERROR read_pdb: Unknown atom in CONECT record for %s in %s\n",
        molecule->info.name, molecule_file_name));
  }

  if (i != molecule->total.bonds)
    exit (fprintf (global.outfile,
      "ERROR read_pdb: Error in bond linked list for %s in %s\n",
      molecule->info.name, molecule_file_name));

/*
* Copy substructure info into molecule structure
* 6/95 te
*/
  for (i = 0, current_subst = first_subst;
    (i < molecule->total.substs) && (current_subst != NULL);
    i++, current_subst = current_subst->next)
    copy_subst (&molecule->subst[i], &current_subst->subst);

  if (i != molecule->total.substs)
    exit (fprintf (global.outfile,
      "ERROR read_pdb: Error in subst linked list for %s in %s\n",
      molecule->info.name, molecule_file_name));

/*
* Update substructure information if more than one substructure read
* (Assume molecule is a macromolecule.)
* 10/95 te
*/
  if (molecule->total.substs > 1)
  {
    for (i = 0; i < molecule->total.substs; i++)
    {
      vstrcpy (&molecule->subst[i].type, "RESIDUE");
      molecule->subst[i].dict_type = 0;
      vstrcpy (&molecule->subst[i].chain, "A");
      vstrcpy (&molecule->subst[i].sub_type, molecule->subst[i].name);

      sprintf (temp, "%d", molecule->subst[i].number);
      vstrcat (&molecule->subst[i].name, temp);
    }
  }

/*
* Free up memory allocated to linked lists
* 6/95 te
*/
  for (previous_atom = first_atom;
    previous_atom != NULL;
    previous_atom = current_atom)
  {
    current_atom = previous_atom->next;
    free_atom (&previous_atom->atom);
    efree ((void **) &previous_atom);
  }

  for (previous_bond = first_bond;
    previous_bond != NULL;
    previous_bond = current_bond)
  {
    current_bond = previous_bond->next;
    free_bond (&previous_bond->bond);
    efree ((void **) &previous_bond);
  }

  for (previous_subst = first_subst;
    previous_subst != NULL;
    previous_subst = current_subst)
  {
    current_subst = previous_subst->next;
    free_subst (&previous_subst->subst);
    efree ((void **) &previous_subst);
  }

  return TRUE;
}


/* ////////////////////////////////////////////////////// */

int write_pdb
(
  MOLECULE	*molecule,
  FILE_NAME	molecule_file_name,
  FILE		*molecule_file
)
{
  int i, j;
  int atom;
  int segment;
  int neighbor;
  STRING80 line;

  fprintf (molecule_file, "HEADER     %s\n", molecule->info.name);
  fprintf (molecule_file, "COMPND     %s\n", molecule->info.comment);
  fprintf (molecule_file, "AUTHOR     Generated by %s version %s\n",
    global.executable, DOCK_VERSION);

/*
* Write out atoms
* 12/96 te
*/
  for (atom = 0; atom < molecule->total.atoms; atom++)
  {
    if
    (
      (molecule->total.layers > 0) &&
      (
        ((segment = molecule->atom[atom].segment_id) != NEITHER) &&
        (molecule->segment[segment].active_flag == FALSE)
      )
    )
      continue;

    memset (line, 0, sizeof (STRING80));
    sprintf (&line[0], "%4s", "ATOM");
    sprintf (&line[6], "%5d", molecule->atom[atom].number);
    sprintf (&line[12], "%-4s",
      molecule->atom[atom].name ? molecule->atom[atom].name : "UNK");
    sprintf (&line[17], "%3.3s",
      (molecule->total.substs > molecule->atom[atom].subst_id) &&
      molecule->subst[molecule->atom[atom].subst_id].name &&
        strcmp (molecule->subst[molecule->atom[atom].subst_id].name, "****") ?
        molecule->subst[molecule->atom[atom].subst_id].name : "UNK");
    sprintf (&line[22], "%4d", molecule->total.substs > 1 ?
      molecule->subst[molecule->atom[atom].subst_id].number :
      molecule->info.output_id);
    sprintf (&line[30], "%8.3f", molecule->coord[atom][0]);
    sprintf (&line[38], "%8.3f", molecule->coord[atom][1]);
    sprintf (&line[46], "%8.3f", molecule->coord[atom][2]);

/*
*   If extended pdb requested, then print out charge and atom type
*   6/95 te
*/
    if (!strcmp (strrchr (molecule_file_name, '.') + 1, "xpdb"))
    {
      sprintf (&line[55], "%10.4f", molecule->atom[atom].charge);
      sprintf (&line[71], "%-5s", molecule->atom[atom].type);
    }

/*
*   Replace null characters with spaces
*   6/95 te
*/
    for (i = 0; i < sizeof (STRING80); i++)
      if (line[i] == '\0') line[i] = ' ';
    line[sizeof (STRING80) - 2] = '\0';

    fprintf (molecule_file, "%s\n", line);
  }

/*
* Write out connectivity
* 12/96 te
*/
  if (molecule->total.bonds)
  {
    for (atom = 0; atom < molecule->total.atoms; atom++)
    {
      if
      (
        (molecule->total.layers > 0) &&
        (
          ((segment = molecule->atom[atom].segment_id) == NEITHER) ||
          (molecule->segment[segment].active_flag == FALSE)
        )
      )
        continue;

      memset (line, 0, sizeof (STRING80));
      sprintf (&line[0], "%-5s", "CONECT");
      sprintf (&line[6], "%5d", molecule->atom[atom].number);

      for (i = j = 0; i < molecule->atom[atom].neighbor_total; i++)
      {
        neighbor = molecule->atom[atom].neighbor[i].id;

        if
        (
          (molecule->total.layers > 0) &&
          (
            ((segment = molecule->atom[neighbor].segment_id) == NEITHER) ||
            (molecule->segment[segment].active_flag == FALSE)
          )
        )
          continue;

        sprintf (&line[j++ * 5 + 11], "%5d",
          molecule->atom[neighbor].number);
      }

/*
*     Replace null characters with spaces
*     6/95 te
*/
      for (i = 0; i < sizeof (STRING80); i++)
        if (line[i] == '\0') line[i] = ' ';
      line[sizeof (STRING80) - 2] = '\0';

      fprintf (molecule_file, "%s\n", line);
    }
  }

  fprintf (molecule_file, "TER\n");

  return TRUE;
}


/* ////////////////////////////////////////////////////////////////////// */

void assign_atom_type (ATOM *atom)
{
  int i;
  STRING5 atom_name[13] =
    {"C", "N", "O", "S", "P", "H", "F", "CL", "BR", "SI", "I", "DU", "LP"};
  STRING5 atom_type[13] =
    {"C.3", "N.3", "O.3", "S.3", "P.3", "H", "F", "Cl", "Br", "Si", "I", "Du", "LP"};

  for (i = 0; i < 13; i++)
    if (strstr (atom->name, atom_name[i]))
      vstrcpy (&atom->type, atom_type[i]);

  if (atom->type == NULL)
    vstrcpy (&atom->type, "Du");
}

