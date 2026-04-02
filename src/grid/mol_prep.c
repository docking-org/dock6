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
#include "mol_prep.h"

typedef struct link_bond
{
  BOND bond;
  struct link_bond *previous, *next;
} LINK_BOND;


int prepare_molecule
(
  MOLECULE *molecule,
  FILE_NAME molecule_file_name,
  FILE *molecule_file,
  LABEL *label,
  int atom_model,
  int need_bonds,
  int assign_chemical,
  int assign_vdw
)
{
  int return_value;

/*
* If bonds are needed, but were not read in, then deduce them
* 10/95 te
*/
  if (need_bonds && (molecule->total.bonds == 0))
    deduce_bonds (&label->vdw, atom_model, molecule);

/*
* Update atom neighbor information
* 9/95 te
*/
  if (molecule->total.bonds > 0)
    atom_neighbors (molecule);

/*
  fprintf (global.outfile, "Breadth search\n");
  for
  (
    i = 0;
    (atom_id = breadth_search
      (
        &search,
        molecule->atom,
        molecule->total.atoms,
        get_atom_neighbor,
        NULL, &i, 1,
        NEITHER,
        i
      )) != EOF;
    i++
  )
  {
    origin = get_search_origin (&search, NEITHER, NEITHER, NEITHER);
    target = get_search_target (&search, NEITHER, NEITHER, NEITHER);
    thread = get_search_thread (&search, NEITHER, NEITHER);

    fprintf (global.outfile, "%d %d %s %d %s %d %s %d\n",
      get_search_radius (&search, NEITHER, NEITHER),
      atom_id, molecule->atom[atom_id].name,
      origin, (origin >= 0) ? molecule->atom[origin].name : "UNK",
      target, (target >= 0) ? molecule->atom[target].name : "UNK",
      thread);
  }
*/

/*
  fprintf (global.outfile, "Depth search\n");
  for
  (
    i = 0;
    (atom_id = depth_search
      (
        &search,
        molecule->atom,
        molecule->total.atoms,
        get_atom_neighbor,
        0,
        NEITHER,
        i
      )) != EOF;
    i++
  )
  {
    origin = get_search_origin (&search, NEITHER, NEITHER, NEITHER);
    target = get_search_target (&search, NEITHER, NEITHER, NEITHER);

    fprintf (global.outfile, "%d %d %s %d %s %d %s\n",
      get_search_radius (&search, NEITHER, NEITHER),
      atom_id, molecule->atom[atom_id].name,
      origin, (origin >= 0) ? molecule->atom[origin].name : "UNK",
      target, (target >= 0) ? molecule->atom[target].name : "UNK");
  }
*/

/*
* Assign new chemical labels.  If the input file was not SPH format
* 10/95 te
*/
  if (assign_chemical)
  {
    if (strcmp (strrchr (molecule_file_name, '.'), ".sph"))
    {
      if ((return_value =
        assign_chemical_labels (&label->chemical, molecule)) != TRUE)
        return return_value;
    }

    else
      read_chemical_labels
       (&label->chemical, molecule, molecule_file_name, molecule_file);
  }

/*
* Assign vdw label types to atoms
* 10/95 te
*/
  if (assign_vdw)
    if ((return_value =
      assign_vdw_labels (&label->vdw, atom_model, molecule)) != TRUE)
      return return_value;

  return TRUE;
}

/* ////////////////////////////////////////////////////////////////////// */

void deduce_bonds
(
  LABEL_VDW *label_vdw,
  int atom_model,
  MOLECULE *molecule
)
{
  int i, j;
  LINK_BOND *first_bond = NULL;
  LINK_BOND *previous_bond = NULL;
  LINK_BOND *current_bond = NULL;
  float cutoff;

  float square_distance (XYZ, XYZ);

/*
* Make preliminary vdw type assignments
* 10/95 te
*/
  if (!assign_vdw_labels (label_vdw, atom_model, molecule))
    exit (EXIT_FAILURE);

/*
* Loop through all possible pairs of atoms, checking interatomic distance.
* If the distance is below some cutoff, then assume the two atoms are bonded.
* 10/95 te
*/
  molecule->total.bonds = 0;

  for (i = 0; i < molecule->total.atoms; i++)
    for (j = i + 1; j < molecule->total.atoms; j++)
    {
      cutoff = 0.53 *
        (label_vdw->member[molecule->atom[i].vdw_id].radius +
        label_vdw->member[molecule->atom[j].vdw_id].radius);

      if (square_distance (molecule->coord[i], molecule->coord[j]) <
        SQR (cutoff))
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
        current_bond->bond.origin = i;
        current_bond->bond.target = j;
        vstrcpy (&current_bond->bond.type, "1");
        molecule->total.bonds++;
      }
    }

  molecule->max.bonds = molecule->total.bonds;
  allocate_bonds (molecule);

/*
* Copy bond info into molecule structure
* 6/95 te
*/
  for (i = 0, current_bond = first_bond;
    (i < molecule->total.bonds) && (current_bond != NULL);
    i++, current_bond = current_bond->next)
    copy_bond (&molecule->bond[i], &current_bond->bond);

  if (i != molecule->total.bonds)
    exit (fprintf (global.outfile,
      "ERROR deduce_bonds: Error reading bond linked list.\n"));

/*
* Free up memory allocated to linked list
* 6/95 te
*/
  for (previous_bond = first_bond;
    previous_bond != NULL;
    previous_bond = current_bond)
  {
    current_bond = previous_bond->next;
    efree ((void **) &previous_bond);
  }
}

