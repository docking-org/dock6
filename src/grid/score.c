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
#include "label.h"
#include "score.h"
#include "io_grid.h"
#include "io_receptor.h"
#include "vector.h"
#include "extsymbols.h"

/* ///////////////////////////////////////////////////////////////

Routine to read in receptor atoms and partition them into a
3-D grid for rapid evaluation of continuous scores.

2/97 te

/////////////////////////////////////////////////////////////// */

void make_receptor_grid
(
  SCORE_GRID		*grid,
  SCORE_ENERGY		*energy,
  LABEL			*label
)
{
  FILE *receptor_file;
  XYZ span;
  int i, j, index, grid_coord[3];
  SLINT **grid_current = NULL;
  SLINT **grid_previous = NULL;

  grid->init_flag = TRUE;

/*
* Open and read receptor atom file
* 9/96 te
*/
  receptor_file = efopen (grid->receptor_file_name, "r", global.outfile);

  read_receptor
  (
    energy,
    label,
    &grid->receptor,
    grid->receptor_file_name,
    receptor_file,
    label->vdw.flag || label->chemical.flag,
    label->chemical.flag,
    label->vdw.flag
  );

/*
* Determine size of box enclosing receptor
* 9/96 te
*/
  for (j = 0; j < 3; j++)
  {
    grid->origin[j] = grid->receptor.coord[0][j];
    span[j] = grid->receptor.coord[0][j];
  }

  for (i = 1; i < grid->receptor.total.atoms; i++)
    for (j = 0; j < 3; j++)
    {
      if (grid->receptor.coord[i][j] < grid->origin[j])
        grid->origin[j] = grid->receptor.coord[i][j];

      if (grid->receptor.coord[i][j] > span[j])
        span[j] = grid->receptor.coord[i][j];
    }

  grid->size = 1;

  for (j = 0; j < 3; j++)
  {
    grid->span[j] = NINT ((span[j] - grid->origin[j]) / grid->spacing) + 1;
    grid->size *= grid->span[j];
  }

  ecalloc
  (
    (void **) &grid->atom,
    grid->size,
    sizeof (SLINT *),
    "receptor atom grid",
    global.outfile
  );

  ecalloc
  (
    (void **) &grid_current,
    grid->size,
    sizeof (SLINT *),
    "receptor atom ptr grid",
    global.outfile
  );

  ecalloc
  (
    (void **) &grid_previous,
    grid->size,
    sizeof (SLINT *),
    "receptor atom ptr grid",
    global.outfile
  );

  for (i = 0; i < grid->receptor.total.atoms; i++)
  {
    for (j = 0; j < 3; j++)
      grid_coord[j] = NINT ((grid->receptor.coord[i][j] - grid->origin[j])
        / grid->spacing);

    index =
      grid->span[0] * grid->span[1] * grid_coord[2] +
      grid->span[0] * grid_coord[1] +
      grid_coord[0];

    ecalloc
    (
      (void **) &grid_current[index],
      1,
      sizeof (SLINT),
      "next atom in receptor atom grid",
      global.outfile
    );

    grid_current[index]->value = i;

    if (grid->atom[index])
      grid_previous[index]->next = grid_current[index];

    else
      grid->atom[index] = grid_current[index];

    grid_previous[index] = grid_current[index];
    grid_current[index] = NULL;
  }

  efree ((void **) &grid_current);
  efree ((void **) &grid_previous);
}


/* ///////////////////////////////////////////////////////////////

Routine to free receptor grid

5/97 te

/////////////////////////////////////////////////////////////// */

void free_receptor_grid (SCORE_GRID *grid)
{
  int i;
  SLINT *previous = NULL;
  SLINT *current = NULL;

  free_molecule (&grid->receptor);

  if (grid->atom)
  {
    for (i = 0; i < grid->size; i++)
    {
      for
      (
        previous = grid->atom[i];
        previous != NULL;
        previous = current
      )
      {
        current = previous->next;
        efree ((void **) &previous);
      }
    }

    efree ((void **) &grid->atom);
  }
}


/* ///////////////////////////////////////////////////////////////

Routine to compute the 1D array index given 3D coordinates.

Return values:
	POSITIVE INTEGER	grid index
	NEITHER (-1)		coordinates outside grid

2/97 te

/////////////////////////////////////////////////////////////// */

int get_grid_index
(
  SCORE_GRID	*grid,
  XYZ		coord
)
{
  int grid_coord[3];

  if (get_grid_coordinate (grid, coord, grid_coord))
    return
      grid->span[0] * grid->span[1] * grid_coord[2] +
      grid->span[0] * grid_coord[1] +
      grid_coord[0];

  else
    return NEITHER;
}


/* ///////////////////////////////////////////////////////////////

Routine to compute the 3D integer grid coordinates given
3D real coordinates.

Return values:
	TRUE		coordinates inside grid
	FALSE		coordinates outside grid

2/97 te

/////////////////////////////////////////////////////////////// */

int get_grid_coordinate
(
  SCORE_GRID	*grid,
  XYZ		coord,
  int		grid_coord[3]
)
{
  int i;
  int inside_flag = TRUE;

  for (i = 0; i < 3; i++)
  {
    grid_coord[i] = NINT ((coord[i] - grid->origin[i]) / grid->spacing);

    if
    (
      (grid_coord[i] < 0) ||
      (grid_coord[i] > grid->span[i] - 1)
    )
      inside_flag = FALSE;
  }

  return inside_flag;
}


/* ///////////////////////////////////////////////////////////////

Routine to identify all receptor atoms near the ligand atom
and compute a continuous intermolecular score.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_inter_score_cont
(
  SCORE_GRID	*grid,
  void		*score,
  float		distance_cutoff,
  void		calc_inter_score
                  (SCORE_GRID *, void *, LABEL *, MOLECULE *,
                  int, int, SCORE_PART *),
  LABEL		*label,
  MOLECULE	*molecule,
  int		atom,
  SCORE_PART	*inter
)
{
  int i, j, k;
  int ilo, jlo, klo;
  int ihi, jhi, khi;
  int grid_cutoff;
  int grid_coord[3];
  int index;
  int rec_atom;
  SLINT *ptr;

  get_grid_coordinate (grid, molecule->coord[atom], grid_coord);
  grid_cutoff = (int) (distance_cutoff / grid->spacing) + 1,

  ilo = MAX (0, (grid_coord[0] - grid_cutoff));
  jlo = MAX (0, (grid_coord[1] - grid_cutoff));
  klo = MAX (0, (grid_coord[2] - grid_cutoff));
  ihi = MIN (grid->span[0], (grid_coord[0] + grid_cutoff + 1));
  jhi = MIN (grid->span[1], (grid_coord[1] + grid_cutoff + 1));
  khi = MIN (grid->span[2], (grid_coord[2] + grid_cutoff + 1));

  for (i = ilo; i < ihi; i++)
    for (j = jlo; j < jhi; j++)
      for (k = klo; k < khi; k++)
      {
        index =
          grid->span[0] * grid->span[1] * k +
          grid->span[0] * j + i;

        for (ptr = grid->atom[index]; ptr; ptr = ptr->next)
        {
          rec_atom = ptr->value;

          calc_inter_score
          (
            grid,
            score,
            label,
            molecule,
            atom,
            rec_atom,
            inter
          );
        }
      }
}


/* ///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

Routine to evaluate the number of intermolecular bumps between
a ligand and receptor

2/97 te

/////////////////////////////////////////////////////////////// */

int check_bump
(
  SCORE_GRID *grid,
  SCORE_BUMP *bump,
  LABEL *label,
  MOLECULE *molecule
)
{
  int atom;
  int segment;
  int index;

  if (bump->grid == NULL)
    exit (fprintf (global.outfile, "ERROR check_bump: No bump grid loaded\n"));

  for (atom = molecule->score.bumpcount = 0;
    (molecule->score.bumpcount <= bump->maximum) && 
    (atom < molecule->total.atoms); atom++)
  {
    if
    (
      ((segment = molecule->atom[atom].segment_id) != NEITHER) &&
      (segment < molecule->total.segments) &&
      (molecule->segment[segment].active_flag != TRUE)
    )
      continue;

    if (molecule->atom[atom].heavy_flag == TRUE)
    {
      index = get_grid_index (grid, molecule->coord[atom]);

      if (index != NEITHER)
      {
        if (grid->version < 3.99)
        {
          if (bump->grid[index] == 'T')
            molecule->score.bumpcount++;
        }

        else
        {
          if (label->vdw.member[molecule->atom[atom].vdw_id].bump_id >=
            bump->grid[index])
            molecule->score.bumpcount++;
        }
      }
    }
  }

  return molecule->score.bumpcount;
}


/* ///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

Routine to compute the intermolecular contact score between
a ligand atom and the receptor.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_inter_contact
(
  SCORE_GRID	*grid,
  SCORE_BUMP	*bump,
  SCORE_CONTACT	*contact,
  LABEL		*label,
  MOLECULE	*molecule,
  int		atom,
  SCORE_PART	*inter
)
{
  if (molecule->atom[atom].heavy_flag == TRUE)
  {
    if (grid->flag)
      calc_inter_contact_grid
        (grid, bump, contact, label, molecule, atom, inter);

    else
      calc_inter_score_cont
      (
        grid,
        (void *) contact,
        contact->distance,
        (void (*)()) calc_inter_contact_cont,
        label,
        molecule,
        atom,
        inter
      );
  }
}


/* ///////////////////////////////////////////////////////////////

Routine to compute the intermolecular contact score using a 
precomputed grid.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_inter_contact_grid
(
  SCORE_GRID	*grid,
  SCORE_BUMP	*bump,
  SCORE_CONTACT	*contact,
  LABEL		*label,
  MOLECULE	*molecule,
  int		atom,
  SCORE_PART	*inter
)
{
  int index;

  index = get_grid_index (grid, molecule->coord[atom]);

  if (index != NEITHER)
  {
    if (label->vdw.member[molecule->atom[atom].vdw_id].bump_id >=
      bump->grid[index])
      inter->total += contact->clash_penalty;

    else
      inter->total += (float) contact->grid[index];
  }
}


/* ///////////////////////////////////////////////////////////////

Routine to compute the intermolecular contact score in a continuous
fashion given a ligand atom and a receptor atom.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_inter_contact_cont
(
  SCORE_GRID	*grid,
  SCORE_CONTACT	*contact,
  LABEL		*label,
  MOLECULE	*molecule,
  int		atom,
  int		rec_atom,
  SCORE_PART	*inter
)
{
  float score;

  if (grid->receptor.atom[rec_atom].heavy_flag == TRUE)
  {
    calc_pairwise_contact
    (
      contact,
      label,
      molecule,
      &grid->receptor,
      atom,
      rec_atom,
      &score
    );

    inter->total += score;
  }
}


/* ///////////////////////////////////////////////////////////////

Routine to compute the intramolecular contact score
between two ligand atoms.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_intra_contact
(
  SCORE_CONTACT *contact,
  LABEL		*label,
  MOLECULE	*molecule,
  int		atomi,
  int		atomj,
  SCORE_PART	*intra
)
{
  float score;

  calc_pairwise_contact
  (
    contact,
    label,
    molecule,
    molecule,
    atomi,
    atomj,
    &score
  );

  intra->total += score;
}


/* ///////////////////////////////////////////////////////////////

Routine to compute the contact score between any two atoms.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_pairwise_contact
(
  SCORE_CONTACT	*contact,
  LABEL		*label,
  MOLECULE	*origin,
  MOLECULE	*target,
  int		origin_atom,
  int		target_atom,
  float		*score
)
{
  float reference;
  float sq_distance;

  if
  (
    origin->atom[origin_atom].heavy_flag &&
    target->atom[target_atom].heavy_flag &&
    ((sq_distance = square_distance
      (origin->coord[origin_atom], target->coord[target_atom]))
        < SQR (contact->distance))
  )
  {
    reference = contact->clash_overlap *
      (label->vdw.member[origin->atom[origin_atom].vdw_id].radius +
      label->vdw.member[target->atom[target_atom].vdw_id].radius);

    if (sq_distance < SQR (reference))
      *score = contact->clash_penalty;

    else
      *score = -1;
  }

  else
    *score = 0;
}


/* ///////////////////////////////////////////////////////////////

Routine to add contact score components.

2/97 te

/////////////////////////////////////////////////////////////// */

void sum_contact
(
  SCORE_PART	*sum,
  SCORE_PART	*increment
)
{
  sum->total += increment->total;
}


/* ///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

Routine to compute the intermolecular energy score between
a ligand atom and the receptor.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_inter_energy
(
  SCORE_GRID	*grid,
  SCORE_ENERGY	*energy,
  LABEL		*label,
  MOLECULE	*molecule,
  int		atom,
  SCORE_PART	*inter
)
{
  if (!energy->vdw_init_flag)
    initialize_vdw_energy (energy, &label->vdw);

  if (label->vdw.member[molecule->atom[atom].vdw_id].well_depth != 0.0)
  {
    if (grid->flag)
      calc_inter_energy_grid (grid, energy, label, molecule, atom, inter);

    else
      calc_inter_score_cont
      (
        grid,
        (void *) energy,
        energy->distance,
        (void (*)()) calc_inter_energy_cont,
        label,
        molecule,
        atom,
        inter
      );
  }
}


/* ///////////////////////////////////////////////////////////////

Routine to initialize vdw parameters.

2/97 te

/////////////////////////////////////////////////////////////// */

void initialize_vdw_energy
(
  SCORE_ENERGY	*energy,
  LABEL_VDW	*label_vdw
)
{
  int i;

  if (!label_vdw->init_flag)
    get_vdw_labels (label_vdw);

  ecalloc
  (
    (void **) &energy->vdwA,
    label_vdw->total,
    sizeof (float),
    "energy vdwA terms",
    global.outfile
  );

  ecalloc
  (
    (void **) &energy->vdwB,
    label_vdw->total,
    sizeof (float),
    "energy vdwB terms",
    global.outfile
  );

  for (i = 0; i < label_vdw->total; i++)
  {
    energy->vdwA[i] = sqrt (label_vdw->member[i].well_depth *
      (float) energy->attractive_exponent /
      (float) (energy->repulsive_exponent - energy->attractive_exponent) *
      pow (2 * label_vdw->member[i].radius, energy->repulsive_exponent));

    energy->vdwB[i] = sqrt (label_vdw->member[i].well_depth *
      (float) energy->repulsive_exponent /
      (float) (energy->repulsive_exponent - energy->attractive_exponent) *
      pow (2 * label_vdw->member[i].radius, energy->attractive_exponent));
  }

  energy->dielectric_factor = 332.0 / energy->dielectric_factor;
  energy->vdw_init_flag = TRUE;
}


/* /////////////////////////////////////////////////////////////// */

void free_vdw_energy (SCORE_ENERGY *energy)
{
  efree ((void **) &energy->vdwA);
  efree ((void **) &energy->vdwB);
}


/* ///////////////////////////////////////////////////////////////

Routine to compute the intermolecular energy score using a 
precomputed grid.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_inter_energy_grid
(
  SCORE_GRID	*grid,
  SCORE_ENERGY	*energy,
  LABEL		*label,
  MOLECULE	*molecule,
  int		atom,
  SCORE_PART	*inter
)
{
  int index[8];
  int out_flag;
  XYZ cube_coord;

  void EXTERNAL_GET_INDEX (XYZ, XYZ, float *, int *, int *, int *, float *);
  float EXTERNAL_GET_VALUE (float *, XYZ, int *);

  if (label->vdw.member[molecule->atom[atom].vdw_id].well_depth != 0.0)
  {
    EXTERNAL_GET_INDEX
    (
      molecule->coord[atom],
      grid->origin,
      &grid->spacing,
      grid->span,
      &out_flag,
      index,
      cube_coord
    );

    if (out_flag == 0)
    {
      inter->vdw +=
        (
          energy->vdwA[molecule->atom[atom].vdw_id] *
            EXTERNAL_GET_VALUE(energy->avdw, cube_coord, index) -
          energy->vdwB[molecule->atom[atom].vdw_id] *
            EXTERNAL_GET_VALUE(energy->bvdw, cube_coord, index)
        ) * energy->scale_vdw;

      inter->electro +=
        molecule->atom[atom].charge *
        EXTERNAL_GET_VALUE(energy->es, cube_coord, index) *
        energy->scale_electro;

      inter->total = inter->vdw + inter->electro;
    }
  }
}


/* ///////////////////////////////////////////////////////////////

Routine to compute the intermolecular energy score in a continuous
fashion given a ligand atom and a receptor atom.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_inter_energy_cont
(
  SCORE_GRID	*grid,
  SCORE_ENERGY	*energy,
  LABEL		*label,
  MOLECULE	*molecule,
  int		atom,
  int		rec_atom,
  SCORE_PART	*inter
)
{
  float vdwA, vdwB, electro;

  if (label->vdw.member[grid->receptor.atom[rec_atom].vdw_id].well_depth != 0.0)
  {
    calc_pairwise_energy
    (
      energy,
      label,
      molecule,
      &grid->receptor,
      atom,
      rec_atom,
      &vdwA,
      &vdwB,
      &electro
    );

    inter->vdw += (vdwA - vdwB) * energy->scale_vdw;
    inter->electro += electro * energy->scale_electro;
    inter->total = inter->vdw + inter->electro;
  }
}


/* ///////////////////////////////////////////////////////////////

Routine to compute the intramolecular energy score
between two ligand atoms.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_intra_energy
(
  SCORE_ENERGY	*energy,
  LABEL		*label,
  MOLECULE	*molecule,
  int		atomi,
  int		atomj,
  SCORE_PART	*intra
)
{
  float vdwA, vdwB, electro;

  calc_pairwise_energy
  (
    energy,
    label,
    molecule,
    molecule,
    atomi,
    atomj,
    &vdwA,
    &vdwB,
    &electro
  );

  intra->vdw += (vdwA - vdwB) * energy->scale_vdw;
  intra->electro += electro * energy->scale_electro;
  intra->total = intra->vdw + intra->electro;

/*
  fprintf (global.outfile, "%s %s %f\n",
    molecule->atom[atomi].name,
    molecule->atom[atomj].name,
    vdwA - vdwB + electro);
*/
}


/* ///////////////////////////////////////////////////////////////

Routine to compute the energy score between any two atoms.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_pairwise_energy
(
  SCORE_ENERGY	*energy,
  LABEL		*label,
  MOLECULE	*origin,
  MOLECULE	*target,
  int		origin_atom,
  int		target_atom,
  float		*vdwA,
  float		*vdwB,
  float		*electro
)
{
  int power_exponent;
  float power_distance;
  float distance;
  int square_flag = TRUE;

  if (!energy->vdw_init_flag)
    initialize_vdw_energy (energy, &label->vdw);

  if
  (
    (label->vdw.member[origin->atom[origin_atom].vdw_id].well_depth != 0.0) &&
    (label->vdw.member[target->atom[target_atom].vdw_id].well_depth != 0.0) &&
    ((distance = square_distance
      (origin->coord[origin_atom], target->coord[target_atom]))
      <= SQR (energy->distance))
  )
  {
    distance = 1.0 / MAX (distance, DISTANCE_MIN);

    if (energy->repulsive_exponent % 2)
    {
      power_exponent = energy->repulsive_exponent;
      distance = sqrt (distance);
      square_flag = FALSE;
    }

    else
      power_exponent = energy->repulsive_exponent / 2;

    POWER (distance, power_exponent, power_distance);

    *vdwA =
      energy->vdwA[origin->atom[origin_atom].vdw_id] *
      energy->vdwA[target->atom[target_atom].vdw_id] *
      power_distance;

    if (square_flag == TRUE)
    {
      if (energy->repulsive_exponent % 2)
      {
        power_exponent = energy->attractive_exponent;
        distance = sqrt (distance);
        square_flag = FALSE;
      }

      else
        power_exponent = energy->attractive_exponent / 2;
    }

    else
      power_exponent = energy->attractive_exponent;

    POWER (distance, power_exponent, power_distance);

    *vdwB =
      energy->vdwB[origin->atom[origin_atom].vdw_id] *
      energy->vdwB[target->atom[target_atom].vdw_id] *
      power_distance;

    *electro =
      energy->dielectric_factor *
      origin->atom[origin_atom].charge *
      target->atom[target_atom].charge *
      (energy->distance_dielectric == TRUE
          ? (square_flag == TRUE ? distance : SQR (distance))
          : (square_flag == TRUE ? sqrt (distance) : distance));
  }

  else
    *vdwA = *vdwB = *electro = 0;
}


/* ///////////////////////////////////////////////////////////////

Routine to add energy score components.

2/97 te

/////////////////////////////////////////////////////////////// */

void sum_energy
(
  SCORE_PART	*sum,
  SCORE_PART	*increment
)
{
  sum->vdw += increment->vdw;
  sum->electro += increment->electro;
  sum->total = sum->vdw + sum->electro;
}


/* ///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

Routine to compute the intermolecular chemical score between
a ligand atom and the receptor.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_inter_chemical
(
  SCORE_GRID	*grid,
  SCORE_ENERGY	*energy,
  SCORE_CHEMICAL *chemical,
  LABEL		*label,
  MOLECULE	*molecule,
  int		atom,
  SCORE_PART	*inter
)
{
  if (!energy->vdw_init_flag)
    initialize_vdw_energy (energy, &label->vdw);

  if (label->chemical.score_table == NULL)
    get_table
    (
      &label->chemical,
      label->chemical.score_file_name,
      &label->chemical.score_table
    );

  if (label->vdw.member[molecule->atom[atom].vdw_id].well_depth != 0.0)
  {
    if (grid->flag)
      calc_inter_chemical_grid
        (grid, energy, chemical, label, molecule, atom, inter);

    else
      calc_inter_score_cont
      (
        grid,
        (void *) energy,
        energy->distance,
        (void (*)()) calc_inter_chemical_cont,
        label,
        molecule,
        atom,
        inter
      );
  }
}


/* ///////////////////////////////////////////////////////////////

Routine to compute the intermolecular chemical score using a 
precomputed grid.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_inter_chemical_grid
(
  SCORE_GRID	*grid,
  SCORE_ENERGY	*energy,
  SCORE_CHEMICAL *chemical,
  LABEL		*label,
  MOLECULE	*molecule,
  int		atom,
  SCORE_PART	*inter
)
{
  int i;
  int index[8];
  int out_flag;
  XYZ cube_coord;
  float vdw, electro;

  void EXTERNAL_GET_INDEX (XYZ, XYZ, float *, int *, int *, int *, float *);
  float EXTERNAL_GET_VALUE (float *, XYZ, int *);

  if (label->vdw.member[molecule->atom[atom].vdw_id].well_depth != 0.0)
  {
    EXTERNAL_GET_INDEX
    (
      molecule->coord[atom],
      grid->origin,
      &grid->spacing,
      grid->span,
      &out_flag,
      index,
      cube_coord
    );

    if (out_flag == 0)
    {
      vdw =
        energy->vdwA[molecule->atom[atom].vdw_id] *
        EXTERNAL_GET_VALUE(energy->avdw, cube_coord, index);

      for (i = 0; i < label->chemical.total; i++)
      {
        if (chemical->grid[i] != NULL)
          vdw -=
            label->chemical.score_table[molecule->atom[atom].chem_id][i] *
            energy->vdwB[molecule->atom[atom].vdw_id] *
            EXTERNAL_GET_VALUE(chemical->grid[i], cube_coord, index);
      }

      electro =
        molecule->atom[atom].charge *
        EXTERNAL_GET_VALUE(energy->es, cube_coord, index);
    }

    else
    {
      vdw = 10.0;
      electro = 10.0;
    }

    inter->vdw += vdw * energy->scale_vdw;
    inter->electro += electro * energy->scale_electro;
    inter->total = inter->vdw + inter->electro;
  }
}


/* ///////////////////////////////////////////////////////////////

Routine to compute the intermolecular chemical score in a continuous
fashion given a ligand atom and a receptor atom.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_inter_chemical_cont
(
  SCORE_GRID	*grid,
  SCORE_ENERGY	*energy,
  LABEL		*label,
  MOLECULE	*molecule,
  int		atom,
  int		rec_atom,
  SCORE_PART	*inter
)
{
  float vdwA, vdwB, electro;

  if (label->vdw.member[grid->receptor.atom[rec_atom].vdw_id].well_depth != 0.0)
  {
    calc_pairwise_energy
    (
      energy,
      label,
      molecule,
      &grid->receptor,
      atom,
      rec_atom,
      &vdwA,
      &vdwB,
      &electro
    );

    inter->vdw += (vdwA - vdwB *
      label->chemical.score_table
        [molecule->atom[atom].chem_id]
        [grid->receptor.atom[rec_atom].chem_id]) * energy->scale_vdw;

    inter->electro += electro * energy->scale_electro;
    inter->total = inter->vdw + inter->electro;
  }
}


/* ///////////////////////////////////////////////////////////////

Routine to compute the intramolecular chemical score
between two ligand atoms.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_intra_chemical
(
  SCORE_ENERGY	*energy,
  LABEL		*label,
  MOLECULE	*molecule,
  int		atomi,
  int		atomj,
  SCORE_PART	*intra
)
{
  float vdwA, vdwB, electro;

  calc_pairwise_energy
  (
    energy,
    label,
    molecule,
    molecule,
    atomi, atomj,
    &vdwA,
    &vdwB,
    &electro
  );

  intra->vdw += (vdwA - vdwB *
    label->chemical.score_table
      [molecule->atom[atomi].chem_id]
      [molecule->atom[atomj].chem_id]) * energy->scale_vdw;

  intra->electro += electro * energy->scale_electro;
  intra->total = intra->vdw + intra->electro;
}


/* ///////////////////////////////////////////////////////////////

Routine to add chemical score components.

2/97 te

/////////////////////////////////////////////////////////////// */

void sum_chemical
(
  SCORE_PART	*sum,
  SCORE_PART	*increment
)
{
  sum->vdw += increment->vdw;
  sum->electro += increment->electro;
  sum->total = sum->vdw + sum->electro;
}


/* ///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

Routine to compute the intermolecular rmsd score between
a ligand atom and the receptor.

2/97 te

/////////////////////////////////////////////////////////////// */

void calc_inter_rmsd
(
  SCORE_GRID	*grid,
  MOLECULE	*molecule,
  int		atom,
  SCORE_PART	*inter
)
{
/*
* Check if input is alright
* 1/97 te
*/
  if (molecule->total.atoms != grid->receptor.total.atoms)
    exit (fprintf (global.outfile,
      "ERROR calc_rmsd_score: molecules have different number of atoms\n"));

  if (molecule->atom[atom].heavy_flag == TRUE)
  {
    inter->vdw += 1;
    inter->electro +=
      square_distance (molecule->coord[atom], grid->receptor.coord[atom]);
  }
}


/* ///////////////////////////////////////////////////////////////

Routine to add rmsd score components.

2/97 te

/////////////////////////////////////////////////////////////// */

void sum_rmsd
(
  SCORE_PART	*sum,
  SCORE_PART	*increment
)
{
  sum->vdw += increment->vdw;
  sum->electro += increment->electro;

  if (sum->vdw > 0)
    sum->total = sqrt (sum->electro / sum->vdw);

  else
    sum->total = 0;
}


/* ///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

Routine to compute the rms deviation between two configurations 
of the same molecule.

2/97 te

/////////////////////////////////////////////////////////////// */

float calc_rmsd
(
  MOLECULE	*mol_ori,
  MOLECULE	*mol_ref
)
{
  int atom;
  int layer_ref;
  int layer_ori;
  int heavy_total = 0;
  float rmsd = 0;

  if (mol_ori->total.atoms != mol_ref->total.atoms)
    exit (fprintf (global.outfile,
      "ERROR calc_rmsd: molecules have different number of atoms\n"));

/*
* If both molecules have an anchor layer, but the anchors are
* different, then set the rmsd arbitrarily high
* 2/97 te
  if ((mol_ori->total.layers > 1) &&
    (mol_ref->total.layers > 1))
  {
    if (mol_ori->layer[0].segment_total != mol_ref->layer[0].segment_total)
      return FLT_MAX;

    for (segment = 0; segment < mol_ref->layer[0].segment_total; segment++)
      if (mol_ori->layer[0].segment[segment] !=
        mol_ref->layer[0].segment[segment])
        return FLT_MAX;
  }
*/

/*
* Loop through atoms
* 2/97 te
*/
  for (atom = 0; atom < mol_ori->total.atoms; atom++)
  {
    if (mol_ori->atom[atom].heavy_flag != TRUE)
      continue;

    if ((mol_ori->total.segments > 0) &&
      (mol_ori->segment[mol_ori->atom[atom].segment_id].active_flag == FALSE))
      continue;
    
    if ((mol_ref->total.layers > 1) && (mol_ori->total.layers > 1))
    {
      layer_ref = mol_ref->segment[mol_ref->atom[atom].segment_id].layer_id;
      layer_ori = mol_ori->segment[mol_ori->atom[atom].segment_id].layer_id;

      if
      (
        ((layer_ref == NEITHER) || (layer_ori == NEITHER)) &&
        (layer_ref != layer_ori)
      )
        return FLT_MAX;
    }

    rmsd += square_distance (mol_ori->coord[atom], mol_ref->coord[atom]);
    heavy_total++;
  }

  if (heavy_total > 0)
    rmsd = sqrt (rmsd / (float) heavy_total);

  else
    rmsd = 0;

  return rmsd;
}


/* ///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

Routine to compute the layer-weighted rms deviation between two configurations 
of the same molecule.

3/97 te

/////////////////////////////////////////////////////////////// */

float calc_layer_rmsd
(
  MOLECULE	*mol_ori,
  MOLECULE	*mol_ref
)
{
  int atom;
  int layer_ref;
  int layer_ori;
  int heavy_total = 0;
  float rmsd = 0;

  if (mol_ori->total.atoms != mol_ref->total.atoms)
    exit (fprintf (global.outfile,
      "ERROR calc_rmsd: molecules have different number of atoms\n"));

/*
* Loop through atoms
* 3/97 te
*/
  for (atom = 0; atom < mol_ori->total.atoms; atom++)
  {
/*
    if (mol_ori->atom[atom].heavy_flag != TRUE)
      continue;
*/

    if ((mol_ori->total.segments > 0) &&
      (mol_ori->segment[mol_ori->atom[atom].segment_id].active_flag == FALSE))
      continue;
    
    layer_ref = mol_ref->segment[mol_ref->atom[atom].segment_id].layer_id;
    layer_ori = mol_ori->segment[mol_ori->atom[atom].segment_id].layer_id;

    if
    (
      ((layer_ref == NEITHER) || (layer_ori == NEITHER)) &&
      (layer_ref != layer_ori)
    )
      return FLT_MAX;

    rmsd +=
      square_distance (mol_ori->coord[atom], mol_ref->coord[atom]) *
      (layer_ori + 1);
    heavy_total += layer_ori + 1;
  }

  if (heavy_total > 0)
    rmsd = sqrt (rmsd / (float) heavy_total);

  else
    rmsd = 0;

  return rmsd;
}


/* ///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

Routine to compute the rms deviation of a segment from two configurations 
of the same molecule.

3/97 te

/////////////////////////////////////////////////////////////// */

float calc_segment_rmsd
(
  MOLECULE	*mol_ori,
  MOLECULE	*mol_ref,
  int		segment
)
{
  int atom;
  int atom_id;
  int atom_total = 0;
  float rmsd = 0;

  if (mol_ori->total.atoms != mol_ref->total.atoms)
    exit (fprintf (global.outfile,
      "ERROR calc_segment_rmsd: molecules have different number of atoms\n"));

/*
* Loop through atoms
* 2/97 te
*/
  for (atom_id = 0; atom_id < mol_ori->segment[segment].atom_total; atom_id++)
  {
    atom = mol_ori->segment[segment].atom[atom_id];

    if ((atom < 0) || (atom >= mol_ori->total.atoms))
      exit (fprintf (global.outfile,
        "ERROR calc_segment_rmsd: segment contains bad atom\n"));

    rmsd += square_distance (mol_ori->coord[atom], mol_ref->coord[atom]);
    atom_total++;                                                       
  }                                                                      

  if (atom_total > 0)
    rmsd = sqrt (rmsd / (float) atom_total);

  else                                       
    rmsd = 0;

  return rmsd;
}             

