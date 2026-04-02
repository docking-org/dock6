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
#include "vector.h"
#include "search.h"
#include "transform.h"
#include "dock.h"
#include "label.h"
#include "score.h"
#include "score_dock.h"
#include "match.h"
#include "orient.h"


/* ///////////////////////////////////////////////////////////// */

int get_orientation
(
  DOCK		*dock,
  ORIENT	*orient,
  LABEL		*label,
  MOLECULE	*mol_ref,
  MOLECULE	*mol_conf,
  MOLECULE	*mol_ori,
  int		molecule_id,
  int		conformation_id,
  int		orientation_id
)
{
/*
* Match molecule to site points
* 10/96 te
*/
  if (orient->match.flag)
  {
    if (orient->random_flag)
      get_random_match
      (
        &orient->match,
        label,
        mol_conf,
        molecule_id,
        conformation_id,
        orientation_id
      );
    
    else if (get_match
      (
        dock,
        &orient->match,
        label,
        mol_conf,
        conformation_id,
        orientation_id
      ) == EOF)
      return EOF;

    extract_clique (&orient->match, mol_ori);

    orient_molecule
    (
      &orient->match.receptor_clique,
      &orient->match.ligand_clique,
      mol_ref, mol_conf, mol_ori
    );

    return TRUE;
  }

  else if (orient->random_flag)
  {
    return
      get_orient
      (
        orient,
        label,
        mol_ref,
        mol_ori,
        molecule_id,
        conformation_id,
        orientation_id
      );
  }

  else
  {
    if (orientation_id == 0)
      return TRUE;

    else
      return EOF;
  }
}


/* ///////////////////////////////////////////////////////////// */

void free_orients
(
  LABEL		*label,
  ORIENT	*orient
)
{
  if (orient->match.flag)
  {
    if (orient->random_flag)
      free_random_matches (&orient->match);

    else
      free_matches (label, &orient->match);
  }

  else if (orient->random_flag)
    free_orient (orient);
}


/* /////////////////////////////////////////////////////////////////

Routine to randomly orient a ligand in a rectangular box
2/97 te

///////////////////////////////////////////////////////////////// */

int get_orient
(
  ORIENT	*orient,
  LABEL		*label,
  MOLECULE	*mol_ref,
  MOLECULE	*mol_ori,
  int		molecule,
  int		conformation,
  int		orientation
)
{
  int i, j;
  XYZ max = {FLT_MIN, FLT_MIN, FLT_MIN};
  XYZ min = {FLT_MAX, FLT_MAX, FLT_MAX};

  if (!molecule && !conformation && !orientation)
  {
    get_site (&orient->match, label);
    get_centers (&orient->match, label);

    for (i = 0; i < orient->match.receptor_site.total.atoms; i++)
    {
      for (j = 0; j < 3; j++)
      {
        min[j] = MIN (min[j], orient->match.receptor_site.coord[i][j]);
        max[j] = MAX (max[j], orient->match.receptor_site.coord[i][j]);
      }
    }

    for (j = 0; j < 3; j++)
    {
      orient->center[j] = (max[j] + min[j]) / 2.0;
      orient->span[j] = (max[j] - min[j]) / 2.0;
    }
  }

  for (i = 0; i < 3; i++)
  {
    mol_ori->transform.translate[i] =
      orient->center[i] +
      orient->span[i] * (1.0 - 2.0 * (float) rand() / (float) RAND_MAX) -
      (orient->match.centers_flag
        ? orient->match.ligand_center.transform.com[i]
        : mol_ref->transform.com[i]);

    mol_ori->transform.rotate[i] =
      1.0 - 2.0 * (float) rand() / (float) RAND_MAX;
  }

  mol_ori->transform.trans_flag =
    mol_ori->transform.rot_flag = TRUE;

  transform_molecule (mol_ori, mol_ref);

  return TRUE;
}

void free_orient (ORIENT *orient)
{
  free_molecule (&orient->match.receptor_site);
  free_molecule (&orient->match.ligand_center);
}
