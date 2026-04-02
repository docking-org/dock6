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
#include "io.h"
#include "mol_prep.h"
#include "io_receptor.h"

MOLECULE temporary = {0};

/* ////////////////////////////////////////////////////////////// */

int read_receptor
(
  SCORE_ENERGY *energy,
  LABEL *label,
  MOLECULE *molecule,
  FILE_NAME molecule_file_name,
  FILE *molecule_file,
  int need_bonds,
  int label_chemical,
  int label_vdw
)
{
  int i, j, return_value;
  int root_atom, root_subst_id, neighbor_subst_id;
  int *old2new = NULL;
  STRING80 line;
  float residue_charge, molecule_charge;
  int *order = NULL;

  int compare_substs ();
  int compare_bonds ();

/*
* Read in coordinates of molecule
* 2/96 te
*/
  return_value =
    read_molecule
    (
      &temporary,
      molecule,
      molecule_file_name,
      molecule_file,
      TRUE
    );

  free_molecule (&temporary);

  if (return_value != TRUE)
    return return_value;

/*
* Prepare molecule for grid calculation
* 2/96 te
*/
  if ((return_value =
    prepare_molecule
    (
      molecule,
      molecule_file_name,
      molecule_file,
      label,
      energy->atom_model,
      need_bonds,
      label_chemical,
      label_vdw
    )) != TRUE)
    return return_value;

/*
* If capping groups were added by SYBYL, then merge them with the substructure
* to which they are attached.
* 10/95 te
*/
  for (i = 0; i < molecule->total.substs; i++)
    if ((strstr (molecule->subst[i].sub_type, "AMN")) ||
      (strstr (molecule->subst[i].sub_type, "CXL")))
    {
      root_atom = molecule->subst[i].root_atom;

      root_subst_id = neighbor_subst_id = molecule->atom[root_atom].subst_id;

      for (j = 0;
        (neighbor_subst_id == root_subst_id) &&
        (j < molecule->atom[root_atom].neighbor_total); j++)
      {
        neighbor_subst_id =
          molecule->atom[molecule->atom[root_atom].neighbor[j].id].subst_id;
      }

      if (neighbor_subst_id != root_subst_id)
      {
        fprintf (global.outfile, "Merging %s %d cap residue with %s %d residue.\n",
          molecule->subst[root_subst_id].name, root_subst_id + 1,
          molecule->subst[neighbor_subst_id].name, neighbor_subst_id + 1);

/*
*       Copy over merged residue info to root atom in cap
*/
        molecule->atom[root_atom].subst_id = neighbor_subst_id;

/*
*       Copy over merged residue info to other atoms in cap
*/
        for (j = 0; j < molecule->atom[root_atom].neighbor_total; j++)
          if (molecule->atom[molecule->atom[root_atom].neighbor[j].id].subst_id
            == root_subst_id)
            molecule->atom[molecule->atom[root_atom].neighbor[j].id].subst_id =
              neighbor_subst_id;

/*
*       Edit neighbor substructure info
*/
        molecule->subst[neighbor_subst_id].inter_bonds = 1;
      }

      else
        exit (fprintf (global.outfile, "Unable to merge substructure %s %d.\n",
          molecule->subst[i].name, i + 1));
    }

/*
* Remove caps from atom and subst records
*/
  for (i = molecule->total.substs - 1; i >= 0; i--)
    if ((strstr (molecule->subst[i].sub_type, "AMN")) ||
      (strstr (molecule->subst[i].sub_type, "CXL")))
    {
/*
*     Reorder substructure ids stored in atom records
*/
      for (j = 0; j < molecule->total.atoms; j++)
        if (molecule->atom[j].subst_id > i)
          molecule->atom[j].subst_id--;

/*
*     Delete cap residues from substructure list
*/
      for (j = i; j < molecule->total.substs - 1; j++)
        molecule->subst[j] = molecule->subst[j + 1];
      molecule->total.substs--;
    }

/*
* Allocate space for use while shuffling molecular components
* 10/95 te
*/
  temporary.max = molecule->max;
  allocate_molecule (&temporary);
  copy_substs (&temporary, molecule);

/*
* Reshuffle substructures so that they are in ascending residue number order
* 10/95 te
*/
  ecalloc
  (
    (void **) &old2new,
    molecule->total.substs,
    sizeof (int),
    "old2new",
    global.outfile
  );

  emalloc
  (
    (void **) &order,
    molecule->total.substs * sizeof (int),
    "order",
    global.outfile
  );

  for (i = 0; i < molecule->total.substs; i++)
    order[i] = i;

  qsort (order, molecule->total.substs, sizeof (int), compare_substs);

  for (i = 0; i < molecule->total.substs; i++)
  {
    old2new[order[i]] = i;

    copy_subst
      (&molecule->subst[i], &temporary.subst[order[i]]);
  }

  for (i = 0; i < molecule->total.atoms; i++)
    molecule->atom[i].subst_id = old2new[molecule->atom[i].subst_id];
  
  efree ((void **) &order);
  efree ((void **) &old2new);

/*
* Reshuffle atoms so that they are also in ascending residue number order.
* Also, report any charged residues.
* 10/95 te
*/
  ecalloc
  (
    (void **) &old2new,
    molecule->max.atoms,
    sizeof (int),
    "old2new",
    global.outfile
  );

  for (i = temporary.total.atoms = 0, molecule_charge = 0.0;
    i < molecule->total.substs; i++)
  {
    for (j = 0, residue_charge = 0.0; j < molecule->total.atoms; j++)
      if (molecule->atom[j].subst_id == i)
      {
        old2new[j] = temporary.total.atoms;
        residue_charge += molecule->atom[j].charge;

        copy_atom
          (&temporary.atom[temporary.total.atoms],
          &molecule->atom[j]);

        temporary.atom[temporary.total.atoms].number =
          temporary.total.atoms + 1;

        copy_coord
          (temporary.coord[temporary.total.atoms],
          molecule->coord[j]);

        temporary.total.atoms++;
      }

    if (ABS (residue_charge) > 0.00001)
    {
      //added this line 3/14 to prevent grid generation for receptors with residues of non integer charge (with a tolerance of 0.01)
      if ( !energy->non_integral_charges && 
      (ABS (residue_charge) - rint(ABS (residue_charge)) >= 0.001) ) 
      {
        fprintf(global.outfile, "%s has a non-integer charge, indicating problems with your receptor file. \n  Quitting without making grids! \n", molecule->subst[i].name);
        return FALSE;
      }
      sprintf (line, "CHARGED RESIDUE %s", molecule->subst[i].name);
      fprintf (global.outfile, "%-40s: %8.3f\n", line, residue_charge);

      molecule_charge += residue_charge;
    }
  }

  sprintf (line, "Total charge on %s", molecule->info.name);
  fprintf (global.outfile, "\n%-40s: %8.3f\n", line, molecule_charge);

  copy_atoms (molecule, &temporary);

/*
* Update bond records to new atom numbering
* 10/95 te
*/
  for (i = 0; i < molecule->total.bonds; i++)
  {
    molecule->bond[i].origin = old2new[molecule->bond[i].origin];
    molecule->bond[i].target = old2new[molecule->bond[i].target];
  }

/*
* Update substructure records to new atom numbering
* 10/95 te
*/
  for (i = 0; i < molecule->total.substs; i++)
    molecule->subst[i].root_atom = old2new[molecule->subst[i].root_atom];

  free_molecule (&temporary);
  efree ((void **) &old2new);

/*
* Reorder bond records also
* 10/95 te
*/
  for (i = 0; i < molecule->total.bonds; i++)
    if (molecule->bond[i].origin > molecule->bond[i].target)
    {
      j = molecule->bond[i].origin;
      molecule->bond[i].origin = molecule->bond[i].target;
      molecule->bond[i].target = j;
    }

  qsort (molecule->bond, molecule->total.bonds, sizeof (BOND), compare_bonds);

  atom_neighbors (molecule);

  return TRUE;
}

  
/* ////////////////////////////////////////////////////////////////// */

int compare_substs (int *subst1, int *subst2)
{
  extern MOLECULE temporary;

  if (temporary.subst[*subst1].chain && temporary.subst[*subst2].chain)
  {
    if (!strcmp
      (temporary.subst[*subst1].chain, temporary.subst[*subst2].chain))
    {
      if (temporary.subst[*subst1].number >
        temporary.subst[*subst2].number)
        return 1;

      else
        return -1;
    }

    else
      return strcmp
        (temporary.subst[*subst1].chain, temporary.subst[*subst2].chain);
  }

  else
  {
    if (temporary.subst[*subst1].number >
      temporary.subst[*subst2].number)
      return 1;

    else
      return -1;
  }
}
  
/* ////////////////////////////////////////////////////////////////// */

int compare_bonds (BOND *bond1, BOND *bond2)
{
  if (bond1->origin == bond2->origin)
  {
    if (bond1->target > bond2->target)
      return 1;

    else
      return -1;
  }

  else if (bond1->origin > bond2->origin)
    return 1;

  else
    return -1;
}
  
