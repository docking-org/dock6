/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
1/97
*/
#include "define.h"
#include "utility.h"
#include "mol.h"
#include "global.h"
#include "label.h"
#include "screen.h"
#include "vector.h"

/* ////////////////////////////////////////////////////////////////////

Routine to update the distance keys for the current conformation.

//////////////////////////////////////////////////////////////////// */

int update_keys
(
  LABEL_CHEMICAL	*label_chemical,
  MOLECULE		*molecule,
  int			conform_id
)
{
  if (label_chemical->screen.fold_flag)
    return
      update_folded_keys (label_chemical, molecule, conform_id);

  else
    return
      update_unfolded_keys (label_chemical, molecule, conform_id);
}


/* ////////////////////////////////////////////////////////////////////

Routine to update the distance keys between labeled atoms
for the current conformation.  The distance keys for all atoms with
the same label are folded on top of each other.

Return values:
	TRUE:	no problems
	FALSE:	big problem (unable to assign chemical labels)

//////////////////////////////////////////////////////////////////// */

int update_folded_keys
(
  LABEL_CHEMICAL	*label_chemical,
  MOLECULE		*molecule,
  int			conform_id
)
{
  int i, j;
  int li, lj;

  if (molecule->info.assign_chem == FALSE)
    if (assign_chemical_labels (label_chemical, molecule) != TRUE)
      return FALSE;

  if (conform_id == 0)
  {
    molecule->transform.fold_flag = TRUE;

    reset_keys (molecule);
    molecule->total.keys = label_chemical->total;
    reallocate_keys (molecule);

    for (i = 0; i < molecule->total.keys; i++)
      for (j = 0; j < molecule->total.keys; j++)
        molecule->key[i][j].count = 0;

    for (i = 0; i < molecule->total.atoms; i++)
    {
      if (molecule->atom[i].heavy_flag != TRUE)
        continue;

      li = molecule->atom[i].chem_id;

      for (j = i + 1; j < molecule->total.atoms; j++)
      {
        if (molecule->atom[j].heavy_flag != TRUE)
          continue;

        lj = molecule->atom[j].chem_id;

        if (li < lj)
          molecule->key[li][lj].count++;

        else
          molecule->key[lj][li].count++;
      }
    }
  }

  for (i = 0; i < molecule->total.atoms; i++)
  {
    if (molecule->atom[i].heavy_flag != TRUE)
      continue;

    li = molecule->atom[i].chem_id;

    for (j = i + 1; j < molecule->total.atoms; j++)
    {
      if (molecule->atom[j].heavy_flag != TRUE)
        continue;

      lj = molecule->atom[j].chem_id;

/*
      if ((li == 2) && (lj == 2))
        fprintf (global.outfile, "atoms %s %s: dist %g\n",
          molecule->atom[i].name,
          molecule->atom[j].name,
          molecule->distance[i][j]);
*/

      if (li < lj)
        molecule->key[li][lj].distance |= 
          get_mask
          (
            &label_chemical->screen,
            dist3 (molecule->coord[i], molecule->coord[j])
          );

      else
        molecule->key[lj][li].distance |=
          get_mask
          (
            &label_chemical->screen,
            dist3 (molecule->coord[i], molecule->coord[j])
          );
    }
  }

  return TRUE;
}


/* ////////////////////////////////////////////////////////////////////

Routine to update the distance keys between ATOMS for the current conformation.

Return values:
	TRUE:	no problems
	FALSE:	big problem (unable to assign chemical labels)

//////////////////////////////////////////////////////////////////// */

int update_unfolded_keys
(
  LABEL_CHEMICAL	*label_chemical,
  MOLECULE		*molecule,
  int			conform_id
)
{
  int i, j;		/* Atom iterators */

  if (molecule->info.assign_chem == FALSE)
    if (assign_chemical_labels (label_chemical, molecule) != TRUE)
      return FALSE;

  if (conform_id == 0)
  {
    molecule->transform.fold_flag = FALSE;

    reset_keys (molecule);
    molecule->total.keys = molecule->total.atoms;
    reallocate_keys (molecule);

    for (i = 0; i < molecule->total.atoms; i++)
      if (molecule->atom[i].heavy_flag == TRUE)
        for (j = i + 1; j < molecule->total.atoms; j++)
          if (molecule->atom[j].heavy_flag == TRUE)
            molecule->key[i][j].count = 1;
  }

  for (i = 0; i < molecule->total.atoms; i++)
    if (molecule->atom[i].heavy_flag == TRUE)
      for (j = i + 1; j < molecule->total.atoms; j++)
        if (molecule->atom[j].heavy_flag == TRUE)
          molecule->key[i][j].distance |=
            get_mask
            (
              &label_chemical->screen,
              dist3 (molecule->coord[i], molecule->coord[j])
            );

  return TRUE;
}


/* ////////////////////////////////////////////////////////////////////

Routine to construct a bit mask encoding the current distance.
11/96 te

//////////////////////////////////////////////////////////////////// */

MASK get_mask
(
  CHEMICAL_SCREEN	*screen,
  float			distance
)
{
  int key;

  key = NINT ((distance - screen->distance_minimum) / 
    screen->distance_interval) + 1;

  key = MAX (key, 0);
  key = MIN (key, screen->interval_total - 1);

  return (MASK) 1 << key;
}


/* ////////////////////////////////////////////////////////////////////

Routine to read unfolded keys from one molecule and record
them as folded keys in another molecule.

//////////////////////////////////////////////////////////////////// */

void fold_keys
(
  LABEL_CHEMICAL	*label_chemical,
  MOLECULE		*copy,
  MOLECULE		*original
)
{
  int i, j;
  int li, lj;

  copy->transform.fold_flag = TRUE;

  if (original->transform.fold_flag == TRUE)
  {
    copy_keys (copy, original);
    return;
  }

  reset_keys (copy);
  copy->total.keys = label_chemical->total;
  reallocate_keys (copy);

  for (i = 0; i < original->total.atoms; i++)
  {
    if (original->atom[i].heavy_flag != TRUE)
      continue;

    li = original->atom[i].chem_id;

    for (j = i + 1; j < original->total.atoms; j++)
    {
      if (original->atom[j].heavy_flag != TRUE)
        continue;

      lj = original->atom[j].chem_id;

      if (li <= lj)
      {
        copy->key[li][lj].count++;

        copy->key[li][lj].distance |=
          original->key[i][j].distance;
      }

      else
      {
        copy->key[lj][li].count++;

        copy->key[lj][li].distance |=
          original->key[i][j].distance;
      }
    }
  }
}


/* ////////////////////////////////////////////////////////////////////

Routine to check if a candidate molecule might include a target.

Return values:
	TRUE:	candidate MAY include target
	FALSE:	candidate CANNOT include target
11/96 te
//////////////////////////////////////////////////////////////////// */

int check_pharmacophore
(
  LABEL_CHEMICAL	*label_chemical,
  float			uncertainty,
  MOLECULE		*target,
  MOLECULE		*candidate
)
{
  static MOLECULE tmp_targ = {0};
  static MOLECULE tmp_cand = {0};

  reset_keys (&tmp_cand);
  reset_keys (&tmp_targ);

  fold_keys (label_chemical, &tmp_targ, target);
  fold_keys (label_chemical, &tmp_cand, candidate);

/*
  copy_keys (&tmp_targ, target);
  copy_keys (&tmp_cand, candidate);
*/

  if (label_chemical->match_table == NULL)
    get_table
    (
      label_chemical,
      label_chemical->match_file_name,
      &label_chemical->match_table
    );

  update_equivalency (label_chemical, &tmp_cand);
  update_uncertainty (label_chemical, uncertainty, &tmp_cand);

  return
    mask_keys (&tmp_targ, &tmp_cand);
}


/* //////////////////////////////////////////////////////////////

Routine to supplement keys to include off-diagonal equivalencies
in the chemical match table.
1/97 te

////////////////////////////////////////////////////////////// */

void update_equivalency
(
  LABEL_CHEMICAL	*label_chemical,
  MOLECULE		*molecule
)
{
  int i, j, k;
  static MOLECULE temporary = {0};

  if (molecule->transform.fold_flag != TRUE)
    exit (fprintf (global.outfile, "update_equivalency: keys not folded\n"));

  copy_keys (&temporary, molecule);

  for (i = 0; i < molecule->total.keys; i++)
    for (j = i + 1; j < molecule->total.keys; j++)
      if (label_chemical->match_table[i][j])
      {
/*
*       Update keys involving only equivalent labels
*       1/97 te
*/
        molecule->key[i][i].count +=
          temporary.key[i][j].count + temporary.key[j][j].count;

        molecule->key[i][j].count +=
          temporary.key[i][i].count + temporary.key[j][j].count;

        molecule->key[j][j].count +=
          temporary.key[i][i].count + temporary.key[i][j].count;

        molecule->key[i][i].distance |=
          temporary.key[i][j].distance | temporary.key[j][j].distance;

        molecule->key[i][j].distance |=
          temporary.key[i][i].distance | temporary.key[j][j].distance;

        molecule->key[j][j].distance |=
          temporary.key[i][i].distance | temporary.key[i][j].distance;

/*
*       Update keys involving equivalent labels and another label
*       1/97 te
*/
        for (k = 0; k < molecule->total.keys; k++)
          if ((k != i) && (k != j))
          {
            if (j < k)
            {
              molecule->key[i][k].count += temporary.key[j][k].count;
              molecule->key[i][k].distance |= temporary.key[j][k].distance;

              molecule->key[j][k].count += temporary.key[i][k].count;
              molecule->key[j][k].distance |= temporary.key[i][k].distance;
            }

            else if (i < k)
            {
              molecule->key[i][k].count += temporary.key[k][j].count;
              molecule->key[i][k].distance |= temporary.key[k][j].distance;

              molecule->key[k][j].count += temporary.key[i][k].count;
              molecule->key[k][j].distance |= temporary.key[i][k].distance;
            }

            else
            {
              molecule->key[k][i].count += temporary.key[k][j].count;
              molecule->key[k][i].distance |= temporary.key[k][j].distance;

              molecule->key[k][j].count += temporary.key[k][i].count;
              molecule->key[k][j].distance |= temporary.key[k][i].distance;
            }
          }
      }
}


/* //////////////////////////////////////////////////////////////

Routine to smear keys according to the uncertainty in distance.
11/96 te

////////////////////////////////////////////////////////////// */

void update_uncertainty
(
  LABEL_CHEMICAL	*label_chemical,
  float			uncertainty,
  MOLECULE		*molecule
)
{
  int i, j, k;
  int smear;

  smear = (int) (uncertainty /
    label_chemical->screen.distance_interval + .9999);

  for (i = 0; i < molecule->total.keys; i++)
    for (j = i; j < molecule->total.keys; j++)
      for (k = 0; k < smear; k++)
        molecule->key[i][j].distance |=
          molecule->key[i][j].distance << 1 |
          molecule->key[i][j].distance >> 1;
}


/* //////////////////////////////////////////////////////////////

Routine to check if candidate molecule contains all keys of target molecule.

Return values:
	TRUE:	candidate molecule contains ALL keys
	FALSE:	candidate molecule doesn't contain all keys
11/96 te

////////////////////////////////////////////////////////////// */

int mask_keys
(
  MOLECULE      *target,
  MOLECULE      *candidate
)
{
  int i, j;

/*
* Verify that a comparison can be made
* 11/96 te
*/
  if ((target->transform.fold_flag != TRUE) ||
    (candidate->transform.fold_flag != TRUE))
    exit (fprintf (global.outfile, "ERROR mask_keys: keys not folded.\n"));

/*
* Check chemical label compositions
* 11/96 te
*/
  for (i = 0; i < target->total.keys; i++)
    for (j = i; j < target->total.keys; j++)
    {
      if (candidate->key[i][j].count < target->key[i][j].count)
        return FALSE;

      if (target->key[i][j].distance !=
        (target->key[i][j].distance & candidate->key[i][j].distance))
        return FALSE;
    }

  return TRUE;
}


/* ////////////////////////////////////////////////////////////////////

Routine to evaluate the similarity of a candidate molecule with
a target molecule.

Return value: FLOAT dissimilarity
11/96 te
//////////////////////////////////////////////////////////////////// */

float check_dissimilarity
(
  LABEL_CHEMICAL	*label_chemical,
  MOLECULE		*target,
  MOLECULE		*candidate
)
{
  static MOLECULE tmp_targ = {0};
  static MOLECULE tmp_cand = {0};

  if (label_chemical->screen_table == NULL)
    get_table
    (
      label_chemical,
      label_chemical->screen_file_name,
      &label_chemical->screen_table
    );

  reset_keys (&tmp_cand);
  reset_keys (&tmp_targ);

  fold_keys (label_chemical, &tmp_targ, target);
  fold_keys (label_chemical, &tmp_cand, candidate);

  return
    1.0 - compare_keys (label_chemical, &tmp_targ, &tmp_cand);
}


/* //////////////////////////////////////////////////////////////

Routine to compare the similarity of two molecules.

Return values:
	TRUE:	candidate molecule contains ALL keys
	FALSE:	candidate molecule doesn't contain all keys
11/96 te

////////////////////////////////////////////////////////////// */

float compare_keys
(
  LABEL_CHEMICAL	*label_chemical,
  MOLECULE      	*target,
  MOLECULE      	*candidate
)
{
  int i, j;
  float numerator = 0;
  float denominator = 0;

/*
* Verify that a comparison can be made
* 11/96 te
*/
  if ((target->transform.fold_flag != TRUE) ||
    (candidate->transform.fold_flag != TRUE))
    exit (fprintf (global.outfile, "ERROR compare_keys: keys not folded.\n"));

/*
* Check chemical label compositions
* 1/97 te
*/
  for (i = 0; i < target->total.keys; i++)
    for (j = i; j < target->total.keys; j++)
    {
      numerator +=
        label_chemical->screen_table[i][j] * (float)
        (MIN (candidate->key[i][j].count, target->key[i][j].count) *
        bit_count (target->key[i][j].distance & candidate->key[i][j].distance));

      denominator +=
        label_chemical->screen_table[i][j] * (float)
        (MAX (candidate->key[i][j].count, target->key[i][j].count) *
        bit_count (target->key[i][j].distance | candidate->key[i][j].distance));
    }

  if (denominator > 0.0)
    return numerator / denominator;

  else
    return 1;
}


/* //////////////////////////////////////////////////////////////

Routine to count the number of bits turned ON in a mask.

Return values:	INT	the number of ON bits
11/96 te

////////////////////////////////////////////////////////////// */

int bit_count (MASK mask)
{
  int i;

  for (i = 0; mask; i++)
    mask >>= 1;

  return i;
}
/* //////////////////////////////////////////////////////////////

Routine to print key.

5/98 te

////////////////////////////////////////////////////////////// */

int print_keys (MOLECULE *mol)
{
  int i, j;

/*
* Check chemical label compositions
* 11/96 te
*/
  for (i = 0; i < mol->total.keys; i++)
    for (j = i; j < mol->total.keys; j++)
      fprintf (global.outfile, "key %d %d %d %x\n", i+1, j+1,
        mol->key[i][j].count, mol->key[i][j].distance);

  return TRUE;
}


