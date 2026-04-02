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


/* ==================================================================== */
int assign_node
(
  NODE *node,
  int include
)
{
  char *branch;
  STRING5 temp;

  (*node).next_total = 0;
  strcpy (temp, strtok (NULL, " "));

  if (isdigit (temp[0]))
  {
    (*node).multiplicity = atoi (temp);
    if ((*node).multiplicity < 0)
    {
      fprintf (global.outfile,
        "Cannot specify negative multiplicity in definition.\n");
      return FALSE;
    }
    strcpy ((*node).type, strtok (NULL, " "));
  }
  else
  {
    (*node).multiplicity = 0;
    strcpy ((*node).type, temp);
  }

  (*node).include = include;

  while (branch = strtok (NULL, " "))
  {
    if ((!strncmp (branch, "(", 1)) || (!strncmp (branch, "[", 1)))
    {
      if ((*node).next_total >= 6)
      {
        fprintf (global.outfile,
          "Cannot exceed 6 substituents for each atom in definition.\n");
        return FALSE;
      }

      ecalloc
      (
        (void **) &(*node).next[(*node).next_total],
        1,
        sizeof (NODE),
        "next label node",
        global.outfile
      );

      if ((!strncmp (branch, "(", 1)) &&
        (!assign_node ((*node).next[(*node).next_total], TRUE)))
        return FALSE;
      if ((!strncmp (branch, "[", 1)) &&
        (!assign_node ((*node).next[(*node).next_total], FALSE)))
        return FALSE;
      (*node).next_total++;
    }

    else if ((!strncmp (branch, ")", 1)) || (!strncmp (branch, "]", 1)))
      return TRUE;

    else
      return FALSE;
  }
  return TRUE;
}

/* ==================================================================== */
void free_node (NODE *node)
{
  int i;

  for (i = 0; i < node->next_total; i++)
  {
    free_node (node->next[i]);
    efree ((void **) &node->next[i]);
  }
}


/* ==================================================================== */
void print_node (NODE *node, int level)
{
  int i;

  if (level)
  {
    if ((*node).include)
      fprintf (global.outfile, " (");
    else
      fprintf (global.outfile, " [");
  }

  if ((*node).multiplicity)
    fprintf (global.outfile, "%d ", (*node).multiplicity);

  fprintf (global.outfile, "%s", (*node).type);

  for (i = 0; i < (*node).next_total; i++)
    print_node ((*node).next[i], level + 1);

  if (level)
  {
    if ((*node).include)
      fprintf (global.outfile, ")");
    else
      fprintf (global.outfile, "]");
  }
}


/* =================================================================== */
int check_type (char *candidate, char *reference)
{
  if ((strstr (candidate, reference)) || (reference[0] == '*'))
    return TRUE;
  else
    return FALSE;
}


/* =================================================================== */
int check_atom
(
  MOLECULE *mol,
  int currentatom,
  NODE *node
)
{
  int i;
  int match;

  if (match = check_type (mol->atom[currentatom].type, (*node).type))
    for (i = 0; i < (*node).next_total; i++)
      if (!check_bonded_atoms
        (mol, currentatom, currentatom, (*node).next[i]))
        match = FALSE;

  return match; 
}

/* =================================================================== */
int check_bonded_atoms
(
  MOLECULE *mol,
  int currentatom,
  int previousatom,
  NODE *node
)
{
  int i, j, nextatom;
  int match, match_count = 0;

  for (i = 0; i < mol->atom[currentatom].neighbor_total; i++)
  {
    nextatom = mol->atom[currentatom].neighbor[i].id;
    match = FALSE;

    if ((nextatom != previousatom) &&
      (match = check_type (mol->atom[nextatom].type, (*node).type)))
    {
      for (j = 0; j < (*node).next_total; j++)
        if (!check_bonded_atoms
          (mol, nextatom, currentatom, (*node).next[j]))
          match = FALSE;
    }

    if (match)
      match_count++;
  }

  if ((*node).multiplicity)
  {
    if ((*node).multiplicity == match_count)
      match = TRUE;
    else
      match = FALSE;
  }
  else
  {
    if (match_count)
      match = TRUE;

    else
      match = FALSE;
  }

  if (match == (*node).include)
    return TRUE; 
  else
    return FALSE; 
}


