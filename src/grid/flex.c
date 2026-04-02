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
#include "label.h"
#include "transform.h"
#include "score.h"
#include "dock.h"
#include "score_dock.h"
#include "vector.h"
#include "flex.h"


/* /////////////////////////////////////////////////////////////////// */
/*
Routine to loop through all anchor fragment conformations.

Return values:
  TRUE	successful generation of anchor conformation
  EOF	unable to generate any more anchor conformations

12/96 te

*/
/* /////////////////////////////////////////////////////////////////// */

int get_anchor_conformation 
(
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*mol_init,
  MOLECULE	*mol_anch,
  int		conformer
) 
{
  int i;
  int conform_found;

/*
* Check whether to continue loop
* 12/96 te
*/
  if
  (
    !label->flex.drive_flag ||
    (mol_anch->layer[0].conform_total == 1)
  )
  {
    if (conformer == 0)
      return TRUE;

    else
      return EOF;
  }

  conform_found =
    get_layer_conformation
    (
      label,
      score,
      mol_anch,
      0,
      conformer
    );

/*
* Check whether conformation was formed.
* 12/96 te
*/
  if (conform_found == TRUE)
    return TRUE;

/*
* If not -- and this is the first requested -- then use input conformation.
* 12/96 te
*/
  else if (conformer == 0)
  {
    for (i = 0; i < mol_anch->total.torsions; i++)
      mol_anch->torsion[i].target_angle = mol_init->torsion[i].target_angle;

    copy_coords (mol_anch, mol_init);

    return TRUE;
  }

  else
    return EOF;
}


/* /////////////////////////////////////////////////////////////////// */

int get_layer_conformation
(
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*molecule,
  int		layer,
  int		conformer
) 
{
  int segment_id;
  int segment;
  int max_conforms;

/*
* If the total conformations is outside the limit, then:
* 1. Assign a random search seed to each segment;
* 2. Initialize the starting segment as always the first
* 3/97 te
*/
  if (layer == 0)
  {
    if (molecule->layer[layer].segment_total == 1)
      max_conforms = 1;

    else
      max_conforms = label->flex.max_conforms *
        (molecule->layer[layer].segment_total - 1);
  }

  else
    max_conforms = label->flex.max_conforms *
      molecule->layer[layer].segment_total;

  if (molecule->layer[layer].conform_total > max_conforms)
  {
    if (conformer >= max_conforms)
      return EOF;

    for
    (
      segment_id = 0;
      segment_id < molecule->layer[layer].segment_total;
      segment_id++
    )
    {
      segment = molecule->layer[layer].segment[segment_id];
      molecule->segment[segment].conform_seed = rand ();
      molecule->segment[segment].conform_count = 0;
    }

    segment_id = 0;
  }

/*
* Otherwise, initialize the starting segment as the first or the last
* 3/97 te
*/
  else
  {
    if (conformer == 0)
    {
      segment_id = 0;
      segment = molecule->layer[layer].segment[segment_id];
      molecule->segment[segment].conform_count = 0;
    }

    else
      segment_id = molecule->layer[layer].segment_total - 1;
  }


/*
* Loop until a complete layer conformation is found, or no more possible
* 3/97 te
*/
  for (;;)
  {
    if
    (
      get_segment_conformation
      (
        label,
        score,
        molecule,
        layer,
        segment_id
      ) == TRUE
    )
    {
      if (++segment_id >= molecule->layer[layer].segment_total)
        return TRUE;

      segment = molecule->layer[layer].segment[segment_id];
      molecule->segment[segment].conform_count = 0;
    }

    else if (--segment_id < 0)
      return EOF;
  }
}


/* /////////////////////////////////////////////////////////////////// */

int get_segment_conformation
(
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*molecule,
  int		layer,
  int		segment_id
) 
{
  int segment;

  segment = molecule->layer[layer].segment[segment_id];

  if (label->flex.drive_flag == TRUE)
  {
    for (;;)
    {
      if
      (
        update_segment_torsions
        (
          &label->flex,
          molecule,
          segment,
          molecule->segment[segment].conform_count++
        ) != TRUE
      )
        return EOF;

      if (check_segment_conformation
        (label, score, molecule, layer, segment_id)
        == TRUE)
        return TRUE;
    }
  }

  else
  {
    if (molecule->segment[segment].conform_count++ == 0)
      return TRUE;

    else
      return EOF;
  }
}


/* /////////////////////////////////////////////////////////////////// */

int update_segment_torsions
(
  LABEL_FLEX	*label_flex,
  MOLECULE	*molecule,
  int		segment,
  int		conformer
) 
{
  int torsion;
  int flex_id;
  int index;

/*
* Assign new torsion angle and apply it
* 3/97 te
*/
  torsion = molecule->segment[segment].torsion_id;

  if (torsion == NEITHER)
    return (conformer == 0) ? TRUE : EOF;

  flex_id = molecule->torsion[torsion].flex_id;

  if (conformer >= label_flex->member[flex_id].torsion_total)
    return EOF;

  else
  {
    index = (molecule->segment[segment].conform_seed + conformer) %
      label_flex->member[flex_id].torsion_total;

    molecule->torsion[torsion].target_angle =
      label_flex->member[flex_id].torsion[index] / 180.0 * PI;

    molecule->transform.tors_flag = TRUE;

    torsion_transform (molecule, torsion);
    return TRUE;
  }
}


/* ///////////////////////////////////////////////////////////////////

Routine to check all inter-segment distances to update distance matrix
and check for clashes

Return values:
	TRUE	segment conformation ok
	FALSE	clash detected
1/97 te

/////////////////////////////////////////////////////////////////// */

int check_segment_conformation
(
  LABEL		*label,
  SCORE		*score,
  MOLECULE	*molecule,
  int		li,			/* Input layer */
  int		sli			/* Input layer segment_id */
) 
{
  int lj;			/* Comparison layer */
  int slj;			/* Segment id in layer list */
  int si, sj;			/* Segment id */
  int asi, asj;			/* Atom id in segment list */
  int ai, aj;			/* Atom id */
  float distance;		/* Square of distance between atoms */
  float reference;		/* Cutoff distance allowed */


/*
* Update atom connection info
* 3/97 te
*/
  initialize_near (&score->near, molecule);

  si = molecule->layer[li].segment[sli];

/*
* Loop through atoms in si segment
* 1/97 te
*/
  for (asi = molecule->segment[si].atom_total - 1; asi >= 0; asi--)
  {
    ai = molecule->segment[si].atom[asi];

/*
*   Loop from current layer to all inner layers
*   1/97 te
*/
    for (lj = li; lj >= 0; lj--)
    {
/*
*     Loop through segments in layer
*       For current layer, only consider previous segments
*       For inner layers, consider all segments
*     1/97 te
*/
      for
      (
        slj = (li == lj ? sli - 1 : molecule->layer[lj].segment_total - 1);
        slj >= 0;
        slj--
      )
      {
        sj = molecule->layer[lj].segment[slj];

/*
*       Loop through atoms in current segment
*       1. Check if unflagged atoms clash with atom ai
*       2. Update portions of distance matrix that have changed
*       1/97 te
*/
        for (asj = molecule->segment[sj].atom_total - 1; asj >= 0; asj--)
        {
          aj = molecule->segment[sj].atom[asj];

          distance = square_distance
            (molecule->coord[ai], molecule->coord[aj]);

          if (!score->near.flag[ai][aj])
          {
            reference = label->flex.clash_overlap *
              (label->vdw.member[molecule->atom[ai].vdw_id].radius +
              label->vdw.member[molecule->atom[aj].vdw_id].radius);

            if (distance < SQR (reference))
              return FALSE;
          }
        } /* End of asj loop */
      } /* End of slj loop */
    } /* End of lj loop */
  } /* End of asi loop */

  return TRUE;
}


