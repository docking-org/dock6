/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/
#include <stdio.h>
#include "define.h"
#include "utility.h"
#include "global.h"
#include "search.h"
#include "mol.h"
#include "vector.h"

/* ///////////////////////////////////////////////////////////// */

void allocate_molecule (MOLECULE *molecule)
{
  allocate_info		(molecule);
  allocate_transform	(molecule);
  allocate_score	(&molecule->score);
  allocate_atoms	(molecule);
  allocate_bonds	(molecule);
  allocate_torsions	(molecule);
  allocate_substs	(molecule);
  allocate_segments	(molecule);
  allocate_layers	(molecule);
  allocate_sets		(molecule);
  allocate_keys		(molecule);
}

/* //////////////////////////////////////////////////////////// */

void reset_molecule (MOLECULE *molecule)
{
  reset_info		(molecule);
  reset_transform 	(molecule);
  reset_score		(&molecule->score);
  reset_atoms		(molecule);
  reset_bonds		(molecule);
  reset_torsions	(molecule);
  reset_substs		(molecule);
  reset_segments	(molecule);
  reset_layers		(molecule);
  reset_sets		(molecule);
  reset_keys		(molecule);
}

/* //////////////////////////////////////////////////////////// */

void free_molecule (MOLECULE *molecule)
{
  free_info		(molecule);
  free_atoms		(molecule);
  free_bonds		(molecule);
  free_torsions		(molecule);
  free_substs		(molecule);
  free_segments		(molecule);
  free_layers		(molecule);
  free_sets		(molecule);
  free_keys		(molecule);
}

/* ///////////////////////////////////////////////////////////// */

void reallocate_molecule (MOLECULE *molecule)
{
  reallocate_atoms	(molecule);
  reallocate_bonds	(molecule);
  reallocate_torsions	(molecule);
  reallocate_substs	(molecule);
  reallocate_segments	(molecule);
  reallocate_layers	(molecule);
  reallocate_sets	(molecule);
  reallocate_keys	(molecule);
}

/* ///////////////////////////////////////////////////////////// */

void copy_molecule (MOLECULE *copy, MOLECULE *original)
{
  copy_info		(copy, original);
  copy_transform	(copy, original);
  copy_score		(&copy->score, &original->score);
  copy_atoms		(copy, original);
  copy_bonds		(copy, original);
  copy_torsions		(copy, original);
  copy_substs		(copy, original);
  copy_segments		(copy, original);
  copy_layers		(copy, original);
  copy_sets		(copy, original);
  copy_keys		(copy, original);
}

/* ///////////////////////////////////////////////////////////// */

void save_molecule (MOLECULE *molecule, FILE *file)
{
  efwrite (&molecule->max, sizeof (ARRAYSIZE), 1, file);
  efwrite (&molecule->total, sizeof (ARRAYSIZE), 1, file);

  save_info		(molecule, file);
  save_transform	(molecule, file);
  save_score		(&molecule->score, file);
  save_atoms		(molecule, file);
  save_bonds		(molecule, file);
  save_torsions		(molecule, file);
  save_substs		(molecule, file);
  save_segments		(molecule, file);
  save_layers		(molecule, file);
  save_sets		(molecule, file);
  save_keys		(molecule, file);
}

/* ///////////////////////////////////////////////////////////// */

void load_molecule (MOLECULE *molecule, FILE *file)
{
  free_molecule (molecule);

  efread (&molecule->max, sizeof (ARRAYSIZE), 1, file);
  efread (&molecule->total, sizeof (ARRAYSIZE), 1, file);

  allocate_molecule	(molecule);

  load_info		(molecule, file);
  load_transform	(molecule, file);
  load_score		(&molecule->score, file);
  load_atoms		(molecule, file);
  load_bonds		(molecule, file);
  load_torsions		(molecule, file);
  load_substs		(molecule, file);
  load_segments		(molecule, file);
  load_layers		(molecule, file);
  load_sets		(molecule, file);
  load_keys		(molecule, file);
}

/* ///////////////////////////////////////////////////////////// */

int get_identifier (void)
{
  static int identifier = 0;
  return ++identifier;
}

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

void allocate_info (MOLECULE *molecule)
{
  reset_info (molecule);
}

/* ///////////////////////////////////////////////////////////// */

void reset_info (MOLECULE *molecule)
{
  efree ((void **) &molecule->info.name);
  efree ((void **) &molecule->info.comment);
  efree ((void **) &molecule->info.molecule_type);
  efree ((void **) &molecule->info.charge_type);
  efree ((void **) &molecule->info.status_bits);
  efree ((void **) &molecule->info.source_file);

  memset (&molecule->info, 0, sizeof (INFO));
  molecule->info.input_id = NEITHER;
}

/* ///////////////////////////////////////////////////////////// */

void free_info (MOLECULE *molecule)
{
  efree ((void **) &molecule->info.name);
  efree ((void **) &molecule->info.comment);
  efree ((void **) &molecule->info.molecule_type);
  efree ((void **) &molecule->info.charge_type);
  efree ((void **) &molecule->info.status_bits);
  efree ((void **) &molecule->info.source_file);

  memset (&molecule->info, 0, sizeof (INFO));
  molecule->info.input_id = NEITHER;
}

/* ///////////////////////////////////////////////////////////// */

void copy_info (MOLECULE *copy, MOLECULE *original)
{
  copy->info.input_id = original->info.input_id;
  copy->info.output_id = original->info.output_id;

  vstrcpy (&copy->info.name, original->info.name);
  vstrcpy (&copy->info.comment, original->info.comment);
  vstrcpy (&copy->info.molecule_type, original->info.molecule_type);
  vstrcpy (&copy->info.charge_type, original->info.charge_type);
  vstrcpy (&copy->info.status_bits, original->info.status_bits);
  vstrcpy (&copy->info.source_file, original->info.source_file);

  copy->info.source_position = original->info.source_position;
  copy->info.file_position = original->info.file_position;

  copy->info.allocated = original->info.allocated;
  copy->info.assign_chem = original->info.assign_chem;
  copy->info.assign_flex = original->info.assign_flex;
  copy->info.assign_vdw = original->info.assign_vdw;
}

/* ///////////////////////////////////////////////////////////// */

void save_info (MOLECULE *molecule, FILE *file)
{
  efwrite (&molecule->info.input_id, sizeof (int), 1, file);
  efwrite (&molecule->info.output_id, sizeof (int), 1, file);

  save_string (&molecule->info.name, file);
  save_string (&molecule->info.comment, file);
  save_string (&molecule->info.molecule_type, file);
  save_string (&molecule->info.charge_type, file);
  save_string (&molecule->info.status_bits, file);
  save_string (&molecule->info.source_file, file);

  efwrite (&molecule->info.source_position, sizeof (long), 1, file);
  efwrite (&molecule->info.file_position, sizeof (long), 1, file);

  efwrite (&molecule->info.allocated, sizeof (int), 1, file);
  efwrite (&molecule->info.assign_chem, sizeof (int), 1, file);
  efwrite (&molecule->info.assign_vdw, sizeof (int), 1, file);
  efwrite (&molecule->info.assign_flex, sizeof (int), 1, file);
}

/* ///////////////////////////////////////////////////////////// */

void load_info (MOLECULE *molecule, FILE *file)
{
  efread (&molecule->info.input_id, sizeof (int), 1, file);
  efread (&molecule->info.output_id, sizeof (int), 1, file);

  load_string (&molecule->info.name, file);
  load_string (&molecule->info.comment, file);
  load_string (&molecule->info.molecule_type, file);
  load_string (&molecule->info.charge_type, file);
  load_string (&molecule->info.status_bits, file);
  load_string (&molecule->info.source_file, file);

  efread (&molecule->info.source_position, sizeof (long), 1, file);
  efread (&molecule->info.file_position, sizeof (long), 1, file);

  efread (&molecule->info.allocated, sizeof (int), 1, file);
  efread (&molecule->info.assign_chem, sizeof (int), 1, file);
  efread (&molecule->info.assign_vdw, sizeof (int), 1, file);
  efread (&molecule->info.assign_flex, sizeof (int), 1, file);
}

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

void allocate_transform (MOLECULE *molecule)
{
  reset_transform (molecule);
}

/* ///////////////////////////////////////////////////////////// */

void reset_transform (MOLECULE *molecule)
{
  memset (&molecule->transform, 0, sizeof (TRANSFORM));
  molecule->transform.refl_flag = NEITHER;
  molecule->transform.tors_flag = NEITHER;
  molecule->transform.anchor_segment = NEITHER;
  molecule->transform.anchor_atom = NEITHER;
}

/* ///////////////////////////////////////////////////////////// */

void copy_transform (MOLECULE *copy, MOLECULE *original)
{
  copy->transform.flag = original->transform.flag;
  copy->transform.trans_flag = original->transform.trans_flag;
  copy->transform.rot_flag = original->transform.rot_flag;
  copy->transform.refl_flag = original->transform.refl_flag;
  copy->transform.tors_flag = original->transform.tors_flag;

  copy->transform.translate[0] = original->transform.translate[0];
  copy->transform.translate[1] = original->transform.translate[1];
  copy->transform.translate[2] = original->transform.translate[2];

  copy->transform.rotate[0] = original->transform.rotate[0];
  copy->transform.rotate[1] = original->transform.rotate[1];
  copy->transform.rotate[2] = original->transform.rotate[2];

  copy->transform.conf_total = original->transform.conf_total;
  copy->transform.anchor_segment = original->transform.anchor_segment;
  copy->transform.anchor_atom = original->transform.anchor_atom;
  copy->transform.active_layer = original->transform.active_layer;
  copy->transform.fold_flag = original->transform.fold_flag;

  copy->transform.com[0] = original->transform.com[0];
  copy->transform.com[1] = original->transform.com[1];
  copy->transform.com[2] = original->transform.com[2];

  copy->transform.rmsd = original->transform.rmsd;
  copy->transform.heavy_total = original->transform.heavy_total;
}

/* ///////////////////////////////////////////////////////////// */

void save_transform (MOLECULE *molecule, FILE *file)
{
  efwrite (&molecule->transform, sizeof (TRANSFORM), 1, file);
}

/* ///////////////////////////////////////////////////////////// */

void load_transform (MOLECULE *molecule, FILE *file)
{
  efread (&molecule->transform, sizeof (TRANSFORM), 1, file);
}

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

void allocate_score (MOL_SCORE *score)
{
  reset_score (score);
}

/* ///////////////////////////////////////////////////////////// */

void reset_score (MOL_SCORE *score)
{
  score->type = 0;
  score->bumpcount = 0;

  reset_score_parts (&score->intra);
  reset_score_parts (&score->inter);

  score->inter.total = INITIAL_SCORE;
  score->total = INITIAL_SCORE;
}

/* ///////////////////////////////////////////////////////////// */

void reset_score_parts (SCORE_PART *part)
{
  part->current_flag = FALSE;
  part->total = 0;
  part->electro = 0;
  part->vdw = 0;
}

/* ///////////////////////////////////////////////////////////// */

void copy_score (MOL_SCORE *copy, MOL_SCORE *original)
{
  copy->type = original->type;
  copy->bumpcount = original->bumpcount;
  copy->intra = original->intra;
  copy->inter = original->inter;
  copy->total = original->total;
}

/* ///////////////////////////////////////////////////////////// */

void add_score_parts (SCORE_PART *sum, SCORE_PART *increment)
{
  sum->vdw += increment->vdw;
  sum->electro += increment->electro;
  sum->total += increment->total;
}

/* ///////////////////////////////////////////////////////////// */

void save_score (MOL_SCORE *score, FILE *file)
{
  efwrite (score, sizeof (MOL_SCORE), 1, file);
}

/* ///////////////////////////////////////////////////////////// */

void load_score (MOL_SCORE *score, FILE *file)
{
  efread (score, sizeof (MOL_SCORE), 1, file);
}

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

void allocate_atoms (MOLECULE *molecule)
{
  int i;

  if (molecule->max.atoms > 0)
  {
    ecalloc
    (
      (void **) &molecule->atom,
      molecule->max.atoms,
      sizeof (ATOM),
      molecule->info.name,
      global.outfile
    );

    ecalloc
    (
      (void **) &molecule->coord,
      molecule->max.atoms,
      sizeof (XYZ),
      molecule->info.name,
      global.outfile
    );

    for (i = 0; i < molecule->total.atoms; i++)
      reset_atom (&molecule->atom[i]);
  }
}

/* ///////////////////////////////////////////////////////////// */

void allocate_atom_neighbors (ATOM *atom)
{
  if (atom->neighbor_max > 0)
  {
    ecalloc
    (
      (void **) &atom->neighbor,
      atom->neighbor_max,
      sizeof (LINK),
      "atom neighbors",
      global.outfile
    );

    reset_atom_neighbors (atom);
  }
}

/* ///////////////////////////////////////////////////////////// */

void reset_atoms (MOLECULE *molecule)
{
  int i;

  if (molecule->max.atoms > 0)
  {
    for (i = 0; i < molecule->max.atoms; i++)
      reset_atom (&molecule->atom[i]);

    memset (molecule->coord, 0, molecule->max.atoms * sizeof (XYZ));
  }

  molecule->total.atoms = 0;
}

/* ///////////////////////////////////////////////////////////// */

void reset_atom (ATOM *atom)
{
  efree ((void **) &atom->name);
  efree ((void **) &atom->type);

  atom->number = 0;
  atom->subst_id = 0;
  atom->chem_id = 0;
  atom->vdw_id = 0;
  atom->heavy_flag = FALSE;
  atom->centrality = NEITHER;
  atom->segment_id = NEITHER;
  atom->flag = 0;
  atom->charge = 0;

  reset_atom_neighbors (atom);
}

/* ///////////////////////////////////////////////////////////// */

void reset_atom_neighbors (ATOM *atom)
{
  int i;

  for (i = 0; i < atom->neighbor_max; i++)
  {
    atom->neighbor[i].id = NEITHER;
    atom->neighbor[i].bond_id = NEITHER;
    atom->neighbor[i].out_flag = FALSE;
  }

  atom->neighbor_total = 0;
}

/* ///////////////////////////////////////////////////////////// */

void free_atoms (MOLECULE *molecule)
{
  int i;

  for (i = 0; i < molecule->max.atoms; i++)
    free_atom (&molecule->atom[i]);

  efree ((void **) &molecule->atom);
  efree ((void **) &molecule->coord);
  molecule->max.atoms = 0;
}

/* ///////////////////////////////////////////////////////////// */

void free_atom (ATOM *atom)
{
  efree ((void **) &atom->name);
  efree ((void **) &atom->type);

  free_atom_neighbors (atom);
}

/* ///////////////////////////////////////////////////////////// */

void free_atom_neighbors (ATOM *atom)
{
  efree ((void **) &atom->neighbor);
  atom->neighbor_max = 0;
}

/* ///////////////////////////////////////////////////////////// */

void reallocate_atoms (MOLECULE *molecule)
{
  if (molecule->total.atoms > molecule->max.atoms)
  {
    free_atoms (molecule);
    molecule->max.atoms = molecule->total.atoms;
    allocate_atoms (molecule);
  }
}

/* ///////////////////////////////////////////////////////////// */

void reallocate_atom_neighbors (ATOM *atom)
{
  if (atom->neighbor_total > atom->neighbor_max)
  {
    free_atom_neighbors (atom);
    atom->neighbor_max = atom->neighbor_total;
    allocate_atom_neighbors (atom);
  }
}

/* ///////////////////////////////////////////////////////////// */

void copy_atoms (MOLECULE *copy, MOLECULE *original)
{
  int i;

  copy->total.atoms = original->total.atoms;
  reallocate_atoms (copy);

  for (i = 0; i < original->total.atoms; i++)
  {
    copy_atom (&copy->atom[i], &original->atom[i]);
    copy_coord (copy->coord[i], original->coord[i]);
  }
}

/* ///////////////////////////////////////////////////////////// */

void copy_atom (ATOM *copy, ATOM *original)
{
  copy->number = original->number;
  copy->subst_id = original->subst_id;
  copy->chem_id = original->chem_id;
  copy->vdw_id = original->vdw_id;
  copy->heavy_flag = original->heavy_flag;
  copy->centrality = original->centrality;
  copy->segment_id = original->segment_id;
  copy->flag = original->flag;
  copy->charge = original->charge;

  vstrcpy (&copy->name, original->name);
  vstrcpy (&copy->type, original->type);

  copy_atom_neighbors (copy, original);
}

/* ///////////////////////////////////////////////////////////// */

void copy_atom_neighbors (ATOM *copy, ATOM *original)
{
  copy->neighbor_total = original->neighbor_total;
  reallocate_atom_neighbors (copy);

  memcpy
    (copy->neighbor, original->neighbor,
    original->neighbor_total * sizeof (LINK));

  copy->neighbor_total = original->neighbor_total;
}

/* ///////////////////////////////////////////////////////////// */

void copy_coords (MOLECULE *copy, MOLECULE *original)
{
  int i;

  if (copy->info.input_id != original->info.input_id)
    exit (fprintf (global.outfile,
      "ERROR copy_coords: "
      "Error copying coordinates to out-of-date molecule structure.\n"));

  for (i = 0; i < original->total.atoms; i++)
    copy_coord (copy->coord[i], original->coord[i]);
}

/* ///////////////////////////////////////////////////////////// */

void copy_coord (XYZ copy, XYZ original)
{
  copy[0] = original[0];
  copy[1] = original[1];
  copy[2] = original[2];
}

/* ///////////////////////////////////////////////////////////// */

void save_atoms (MOLECULE *molecule, FILE *file)
{
  int i;

  if (molecule->max.atoms > 0)
  {
    for (i = 0; i < molecule->max.atoms; i++)
      save_atom (&molecule->atom[i], file);

    efwrite
    (
      molecule->coord,
      sizeof (XYZ),
      molecule->max.atoms,
      file
    );
  }
}

/* ///////////////////////////////////////////////////////////// */

void save_atom (ATOM *atom, FILE *file)
{
  efwrite (atom, sizeof (ATOM), 1, file);

  save_string (&atom->name, file);
  save_string (&atom->type, file);

  save_atom_neighbors (atom, file);
}

/* ///////////////////////////////////////////////////////////// */

void save_atom_neighbors (ATOM *atom, FILE *file)
{
  efwrite (&atom->neighbor_total, sizeof (int), 1, file);
  efwrite (&atom->neighbor_max, sizeof (int), 1, file);
  efwrite (atom->neighbor, sizeof (LINK), atom->neighbor_total, file);
}

/* ///////////////////////////////////////////////////////////// */

void load_atoms (MOLECULE *molecule, FILE *file)
{
  int i;

  if (molecule->max.atoms > 0)
  {
    for (i = 0; i < molecule->max.atoms; i++)
      load_atom (&molecule->atom[i], file);

    efread
    (
      molecule->coord,
      sizeof (XYZ),
      molecule->max.atoms,
      file
    );
  }
}

/* ///////////////////////////////////////////////////////////// */

void load_atom (ATOM *atom, FILE *file)
{
  efread (atom, sizeof (ATOM), 1, file);

  atom->name = NULL;
  atom->type = NULL;
  atom->neighbor = NULL;

  load_string (&atom->name, file);
  load_string (&atom->type, file);

  load_atom_neighbors (atom, file);
}

/* ///////////////////////////////////////////////////////////// */

void load_atom_neighbors (ATOM *atom, FILE *file)
{
  efread (&atom->neighbor_total, sizeof (int), 1, file);
  efread (&atom->neighbor_max, sizeof (int), 1, file);

  allocate_atom_neighbors (atom);

  efread (atom->neighbor, sizeof (LINK), atom->neighbor_total, file);
}

/* ///////////////////////////////////////////////////////////// */

int get_atom_neighbor
(
  void		*atom,
  int		atom_id,
  int		neighbor_id
)
{
  if (neighbor_id < ((ATOM *) atom)[atom_id].neighbor_total)
    return ((ATOM *) atom)[atom_id].neighbor[neighbor_id].id;

  else
    return EOF;
}

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

void allocate_bonds (MOLECULE *molecule)
{
  int i;

  if (molecule->max.bonds > 0)
    ecalloc
    (
      (void **) &molecule->bond,
      molecule->max.bonds,
      sizeof (BOND),
      molecule->info.name,
      global.outfile
    );

  for (i = 0; i < molecule->max.bonds; i++)
    reset_bond (&molecule->bond[i]);
}

/* ///////////////////////////////////////////////////////////// */

void reset_bonds (MOLECULE *molecule)
{
  int i;

  for (i = 0; i < molecule->max.bonds; i++)
    reset_bond (&molecule->bond[i]);

  molecule->total.bonds = 0;
}

/* ///////////////////////////////////////////////////////////// */

void reset_bond (BOND *bond)
{
  bond->id = 0;
  bond->origin = NEITHER;
  bond->target = NEITHER;
  bond->ring_flag = 0;
  bond->flex_id = 0;
  efree ((void **) &bond->type);
}

/* ///////////////////////////////////////////////////////////// */

void free_bonds (MOLECULE *molecule)
{
  int i;

  for (i = 0; i < molecule->max.bonds; i++)
    free_bond (&molecule->bond[i]);

  efree ((void **) &molecule->bond);
  molecule->max.bonds = 0;
}

/* ///////////////////////////////////////////////////////////// */

void free_bond (BOND *bond)
{
  efree ((void **) &bond->type);
}

/* ///////////////////////////////////////////////////////////// */

void reallocate_bonds (MOLECULE *molecule)
{
  if (molecule->total.bonds > molecule->max.bonds)
  {
    free_bonds (molecule);
    molecule->max.bonds = molecule->total.bonds;
    allocate_bonds (molecule);
  }
}

/* ///////////////////////////////////////////////////////////// */

void copy_bonds (MOLECULE *copy, MOLECULE *original)
{
  int i;

  copy->total.bonds = original->total.bonds;
  reallocate_bonds (copy);

  for (i = 0; i < original->total.bonds; i++)
    copy_bond (&copy->bond[i], &original->bond[i]);

  copy->total.bonds = original->total.bonds;
}

/* ///////////////////////////////////////////////////////////// */

void copy_bond (BOND *copy, BOND *original)
{
  copy->id = original->id;
  copy->origin = original->origin;
  copy->target = original->target;
  copy->ring_flag = original->ring_flag;
  copy->flex_id = original->flex_id;
  vstrcpy (&copy->type, original->type);
}

/* ///////////////////////////////////////////////////////////// */

void save_bonds (MOLECULE *molecule, FILE *file)
{
  int i;

  for (i = 0; i < molecule->max.bonds; i++)
    save_bond (&molecule->bond[i], file);
}

/* ///////////////////////////////////////////////////////////// */

void save_bond (BOND *bond, FILE *file)
{
  efwrite (bond, sizeof (BOND), 1, file);
  save_string (&bond->type, file);
}

/* ///////////////////////////////////////////////////////////// */

void load_bonds (MOLECULE *molecule, FILE *file)
{
  int i;

  for (i = 0; i < molecule->max.bonds; i++)
    load_bond (&molecule->bond[i], file);
}

/* ///////////////////////////////////////////////////////////// */

void load_bond (BOND *bond, FILE *file)
{
  efread (bond, sizeof (BOND), 1, file);
  bond->type = NULL;
  load_string (&bond->type, file);
}

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

void allocate_torsions (MOLECULE *molecule)
{
  int i;

  if (molecule->max.torsions > 0)
    ecalloc
    (
      (void **) &molecule->torsion,
      molecule->max.torsions,
      sizeof (TORSION),
      molecule->info.name,
      global.outfile
    );

  for (i = 0; i < molecule->max.torsions; i++)
    reset_torsion (&molecule->torsion[i]);
}

/* ///////////////////////////////////////////////////////////// */

void reset_torsions (MOLECULE *molecule)
{
  int i;

  for (i = 0; i < molecule->max.torsions; i++)
    reset_torsion (&molecule->torsion[i]);

  molecule->total.torsions = 0;
}

/* ///////////////////////////////////////////////////////////// */

void reset_torsion (TORSION *torsion)
{
  torsion->flex_id = 0;
  torsion->bond_id = NEITHER;
  torsion->segment_id = NEITHER;
  torsion->periph_flag = 0;
  torsion->reverse_flag = 0;
  torsion->origin = NEITHER;
  torsion->target = NEITHER;
  torsion->origin_neighbor = NEITHER;
  torsion->target_neighbor = NEITHER;
  torsion->current_angle = 0;
  torsion->target_angle = 0;
}

/* ///////////////////////////////////////////////////////////// */

void free_torsions (MOLECULE *molecule)
{
  efree ((void **) &molecule->torsion);
  molecule->max.torsions = 0;
}

/* ///////////////////////////////////////////////////////////// */

void reallocate_torsions (MOLECULE *molecule)
{
  if (molecule->total.torsions > molecule->max.torsions)
  {
    free_torsions (molecule);
    molecule->max.torsions = molecule->total.torsions;
    allocate_torsions (molecule);
  }
}

/* ///////////////////////////////////////////////////////////// */

void copy_torsions (MOLECULE *copy, MOLECULE *original)
{
  int i;

  copy->total.torsions = original->total.torsions;
  reallocate_torsions (copy);

  for (i = 0; i < original->total.torsions; i++)
    copy_torsion (&copy->torsion[i], &original->torsion[i]);

  copy->total.torsions = original->total.torsions;
}

/* ///////////////////////////////////////////////////////////// */

void copy_torsion (TORSION *copy, TORSION *original)
{
  copy->flex_id = original->flex_id;
  copy->bond_id = original->bond_id;
  copy->segment_id = original->segment_id;
  copy->periph_flag = original->periph_flag;
  copy->reverse_flag = original->reverse_flag;
  copy->origin = original->origin;
  copy->target = original->target;
  copy->origin_neighbor = original->origin_neighbor;
  copy->target_neighbor = original->target_neighbor;
  copy->current_angle = original->current_angle;
  copy->target_angle = original->target_angle;
}

/* ///////////////////////////////////////////////////////////// */

void save_torsions (MOLECULE *molecule, FILE *file)
{
  if (molecule->max.torsions > 0)
    efwrite (molecule->torsion, sizeof (TORSION), molecule->max.torsions, file);
}

/* ///////////////////////////////////////////////////////////// */

void load_torsions (MOLECULE *molecule, FILE *file)
{
  if (molecule->max.torsions > 0)
    efread (molecule->torsion, sizeof (TORSION), molecule->max.torsions, file);
}

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

void allocate_substs (MOLECULE *molecule)
{
  int i;

  if (molecule->max.substs > 0)
    ecalloc
    (
      (void **) &molecule->subst,
      molecule->max.substs,
      sizeof (SUBST),
      molecule->info.name,
      global.outfile
    );

  for (i = 0; i < molecule->max.substs; i++)
    reset_subst (&molecule->subst[i]);
}

/* ///////////////////////////////////////////////////////////// */

void reset_substs (MOLECULE *molecule)
{
  int i;

  for (i = 0; i < molecule->max.substs; i++)
    reset_subst (&molecule->subst[i]);

  molecule->total.substs = 0;
}

/* ///////////////////////////////////////////////////////////// */

void reset_subst (SUBST *subst)
{
  efree ((void **) &subst->name);
  efree ((void **) &subst->type);
  efree ((void **) &subst->chain);
  efree ((void **) &subst->sub_type);
  efree ((void **) &subst->status);

  memset (subst, 0, sizeof (SUBST));
}

/* ///////////////////////////////////////////////////////////// */

void free_substs (MOLECULE *molecule)
{
  int i;

  for (i = 0; i < molecule->max.substs; i++)
    free_subst (&molecule->subst[i]);

  efree ((void **) &molecule->subst);
  molecule->max.substs = 0;
}

/* ///////////////////////////////////////////////////////////// */

void free_subst (SUBST *subst)
{
  efree ((void **) &subst->name);
  efree ((void **) &subst->type);
  efree ((void **) &subst->chain);
  efree ((void **) &subst->sub_type);
  efree ((void **) &subst->status);
}

/* ///////////////////////////////////////////////////////////// */

void reallocate_substs (MOLECULE *molecule)
{
  if (molecule->total.substs > molecule->max.substs)
  {
    free_substs (molecule);
    molecule->max.substs = molecule->total.substs;
    allocate_substs (molecule);
  }
}

/* ///////////////////////////////////////////////////////////// */

void copy_substs (MOLECULE *copy, MOLECULE *original)
{
  int i;

  copy->total.substs = original->total.substs;
  reallocate_substs (copy);

  for (i = 0; i < original->total.substs; i++)
    copy_subst (&copy->subst[i], &original->subst[i]);

  copy->total.substs = original->total.substs;
}

/* ///////////////////////////////////////////////////////////// */

void copy_subst (SUBST *copy, SUBST *original)
{
  copy->number = original->number;
  copy->root_atom = original->root_atom;
  copy->dict_type = original->dict_type;
  copy->inter_bonds = original->inter_bonds;

  vstrcpy (&copy->name, original->name);
  vstrcpy (&copy->type, original->type);
  vstrcpy (&copy->chain, original->chain);
  vstrcpy (&copy->sub_type, original->sub_type);
  vstrcpy (&copy->status, original->status);
}

/* ///////////////////////////////////////////////////////////// */

void save_substs (MOLECULE *molecule, FILE *file)
{
  int i;

  for (i = 0; i < molecule->max.substs; i++)
    save_subst (&molecule->subst[i], file);
}

/* ///////////////////////////////////////////////////////////// */

void save_subst (SUBST *subst, FILE *file)
{
  efwrite (subst, sizeof (SUBST), 1, file);

  save_string (&subst->name, file);
  save_string (&subst->type, file);
  save_string (&subst->chain, file);
  save_string (&subst->sub_type, file);
  save_string (&subst->status, file);
}

/* ///////////////////////////////////////////////////////////// */

void load_substs (MOLECULE *molecule, FILE *file)
{
  int i;

  for (i = 0; i < molecule->max.substs; i++)
    load_subst (&molecule->subst[i], file);
}

/* ///////////////////////////////////////////////////////////// */

void load_subst (SUBST *subst, FILE *file)
{
  efread (subst, sizeof (SUBST), 1, file);

  subst->name = NULL;
  subst->type = NULL;
  subst->chain = NULL;
  subst->sub_type = NULL;
  subst->status = NULL;

  load_string (&subst->name, file);
  load_string (&subst->type, file);
  load_string (&subst->chain, file);
  load_string (&subst->sub_type, file);
  load_string (&subst->status, file);
}

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

void allocate_segments (MOLECULE *molecule)
{
  int i, j;

  ecalloc
  (
    (void **) &molecule->segment,
    molecule->max.segments,
    sizeof (SEGMENT),
    molecule->info.name,
    global.outfile
  );

  for (i = 0; i < molecule->max.segments; i++)
    reset_segment (&molecule->segment[i]);
}

/* ///////////////////////////////////////////////////////////// */

void allocate_segment_atoms (SEGMENT *segment)
{
  if (segment->atom_max > 0)
  {
    ecalloc
    (
      (void **) &segment->atom,
      segment->atom_max,
      sizeof (int),
      "segment atoms",
      global.outfile
    );
  }
}

/* ///////////////////////////////////////////////////////////// */

void allocate_segment_neighbors (SEGMENT *segment)
{
  int i;

  if (segment->neighbor_max > 0)
  {
    ecalloc
    (
      (void **) &segment->neighbor,
      segment->neighbor_max,
      sizeof (LINK),
      "segment neighbors",
      global.outfile
    );

    for (i = 0; i < segment->neighbor_max; i++)
    {
      segment->neighbor[i].id = NEITHER;
      segment->neighbor[i].bond_id = NEITHER;
      segment->neighbor[i].out_flag = FALSE;
    }
  }
}

/* ///////////////////////////////////////////////////////////// */

void reset_segments (MOLECULE *molecule)
{
  int i;

  for (i = 0; i < molecule->max.segments; i++)
    reset_segment (&molecule->segment[i]);

  molecule->total.segments = 0;
}

/* ///////////////////////////////////////////////////////////// */

void reset_segment (SEGMENT *segment)
{
  segment->id = NEITHER;
  segment->torsion_id = NEITHER;
  segment->layer_id = NEITHER;
  segment->periph_flag = NEITHER;
  segment->heavy_total = 0;
  segment->active_flag = FALSE;
  segment->min_flag = FALSE;
  segment->conform_count = 0;
  segment->conform_total = 1;
  reset_score (&segment->score);

  memset (segment->atom, 0, segment->atom_max * sizeof (int));
  segment->atom_total = 0;

  reset_segment_neighbors (segment);
}

/* ///////////////////////////////////////////////////////////// */

void reset_segment_neighbors (SEGMENT *segment)
{
  int i;

  for (i = 0; i < segment->neighbor_max; i++)
  {
    segment->neighbor[i].id = NEITHER;
    segment->neighbor[i].bond_id = NEITHER;
    segment->neighbor[i].out_flag = FALSE;
  }

  segment->neighbor_total = 0;
}

/* ///////////////////////////////////////////////////////////// */

void free_segments (MOLECULE *molecule)
{
  int i;

  for (i = 0; i < molecule->max.segments; i++)
  {
    free_segment_atoms (&molecule->segment[i]);
    free_segment_neighbors (&molecule->segment[i]);
  }

  efree ((void **) &molecule->segment);
  molecule->max.segments = 0;
}

/* ///////////////////////////////////////////////////////////// */

void free_segment_atoms (SEGMENT *segment)
{
  efree ((void **) &segment->atom);
  segment->atom_max = 0;
}

/* ///////////////////////////////////////////////////////////// */

void free_segment_neighbors (SEGMENT *segment)
{
  efree ((void **) &segment->neighbor);
  segment->neighbor_max = 0;
}

/* ///////////////////////////////////////////////////////////// */

void reallocate_segments (MOLECULE *molecule)
{
  if (molecule->total.segments > molecule->max.segments)
  {
    free_segments (molecule);
    molecule->max.segments = molecule->total.segments;
    allocate_segments (molecule);
  }
}

/* ///////////////////////////////////////////////////////////// */

void reallocate_segment_atoms (SEGMENT *segment)
{
  if (segment->atom_total > segment->atom_max)
  {
    free_segment_atoms (segment);
    segment->atom_max = segment->atom_total;
    allocate_segment_atoms (segment);
  }
}

/* ///////////////////////////////////////////////////////////// */

void reallocate_segment_neighbors (SEGMENT *segment)
{
  if (segment->neighbor_total > segment->neighbor_max)
  {
    free_segment_neighbors (segment);
    segment->neighbor_max = segment->neighbor_total;
    allocate_segment_neighbors (segment);
  }
}

/* ///////////////////////////////////////////////////////////// */

void copy_segments (MOLECULE *copy, MOLECULE *original)
{
  int i, j;
  int atom;

  copy->total.segments = original->total.segments;
  reallocate_segments (copy);

  for (i = 0; i < original->total.segments; i++)
  {
    copy_segment (&copy->segment[i], &original->segment[i]);

    for (j = 0; j < copy->segment[i].atom_total; j++)
    {
      atom = copy->segment[i].atom[j];

      if ((atom >= 0) && (atom < copy->total.atoms))
        copy->atom[atom].segment_id = i;

      else
        exit (fprintf (global.outfile,
          "ERROR copy_segments: unable to copy atoms\n"));
    }
  }

  if ((copy->total.atoms > 0) &&
    (copy->total.segments > 1) &&
    (copy->transform.anchor_atom != NEITHER))
    revise_atom_neighbors (copy);
/*
*/
}

/* ///////////////////////////////////////////////////////////// */

void copy_segment (SEGMENT *copy, SEGMENT *original)
{
  copy->id = original->id;
  copy->torsion_id = original->torsion_id;
  copy->layer_id = original->layer_id;
  copy->periph_flag = original->periph_flag;
  copy->heavy_total = original->heavy_total;
  copy->active_flag = original->active_flag;
  copy->min_flag = original->min_flag;
  copy->conform_count = original->conform_count;
  copy->conform_total = original->conform_total;
  copy_score (&copy->score, &original->score);

  copy->atom_total = original->atom_total;
  reallocate_segment_atoms (copy);
  memcpy
    (copy->atom, original->atom,
    original->atom_total * sizeof (int));

  copy->neighbor_total = original->neighbor_total;
  reallocate_segment_neighbors (copy);
  memcpy
    (copy->neighbor, original->neighbor,
    original->neighbor_total * sizeof (LINK));
}

/* ///////////////////////////////////////////////////////////// */

void save_segments (MOLECULE *molecule, FILE *file)
{
  int i;

  for (i = 0; i < molecule->max.segments; i++)
    save_segment (&molecule->segment[i], file);
}

/* ///////////////////////////////////////////////////////////// */

void save_segment (SEGMENT *segment, FILE *file)
{
  efwrite (&segment->id, sizeof (int), 1, file);
  efwrite (&segment->torsion_id, sizeof (int), 1, file);
  efwrite (&segment->layer_id, sizeof (int), 1, file);
  efwrite (&segment->heavy_total, sizeof (int), 1, file);
  efwrite (&segment->periph_flag, sizeof (int), 1, file);
  efwrite (&segment->active_flag, sizeof (int), 1, file);
  efwrite (&segment->min_flag, sizeof (int), 1, file);
  efwrite (&segment->conform_count, sizeof (int), 1, file);
  efwrite (&segment->conform_total, sizeof (int), 1, file);
  save_score (&segment->score, file);

  efwrite (&segment->atom_total, sizeof (int), 1, file);
  efwrite (&segment->atom_max, sizeof (int), 1, file);
  efwrite (segment->atom, sizeof (int), segment->atom_max, file);

  efwrite (&segment->neighbor_total, sizeof (int), 1, file);
  efwrite (&segment->neighbor_max, sizeof (int), 1, file);
  efwrite (segment->neighbor, sizeof (LINK), segment->neighbor_max, file);
}

/* ///////////////////////////////////////////////////////////// */

void load_segments (MOLECULE *molecule, FILE *file)
{
  int i;

  for (i = 0; i < molecule->max.segments; i++)
    load_segment (&molecule->segment[i], file);
}

/* ///////////////////////////////////////////////////////////// */

void load_segment (SEGMENT *segment, FILE *file)
{
  efread (&segment->id, sizeof (int), 1, file);
  efread (&segment->torsion_id, sizeof (int), 1, file);
  efread (&segment->layer_id, sizeof (int), 1, file);
  efread (&segment->heavy_total, sizeof (int), 1, file);
  efread (&segment->periph_flag, sizeof (int), 1, file);
  efread (&segment->active_flag, sizeof (int), 1, file);
  efread (&segment->min_flag, sizeof (int), 1, file);
  efread (&segment->conform_count, sizeof (int), 1, file);
  efread (&segment->conform_total, sizeof (int), 1, file);
  load_score (&segment->score, file);

  efread (&segment->atom_total, sizeof (int), 1, file);
  efread (&segment->atom_max, sizeof (int), 1, file);
  allocate_segment_atoms (segment);
  efread (segment->atom, sizeof (int), segment->atom_max, file);

  efread (&segment->neighbor_total, sizeof (int), 1, file);
  efread (&segment->neighbor_max, sizeof (int), 1, file);
  allocate_segment_neighbors (segment);
  efread (segment->neighbor, sizeof (LINK), segment->neighbor_max, file);
}

/* ///////////////////////////////////////////////////////////// */

int get_segment_neighbor
(
  void		*segment,
  int		segment_id,
  int		neighbor_id
)
{
  if (neighbor_id < ((SEGMENT *) segment)[segment_id].neighbor_total)
    return ((SEGMENT *) segment)[segment_id].neighbor[neighbor_id].id;

  else
    return EOF;
}

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

void allocate_layers (MOLECULE *molecule)
{
  int i;

  if (molecule->max.layers > 0)
  {
    ecalloc
    (
      (void **) &molecule->layer,
      molecule->max.layers,
      sizeof (LAYER),
      molecule->info.name,
      global.outfile
    );
  }

  for (i = 0; i < molecule->max.layers; i++)
    reset_layer (molecule, i);
}

/* ///////////////////////////////////////////////////////////// */

void allocate_layer_segments (LAYER *layer)
{
  if (layer->segment_max > 0)
  {
    ecalloc
    (
      (void **) &layer->segment,
      layer->segment_max,
      sizeof (int),
      "layer atoms",
      global.outfile
    );
  }
}

/* ///////////////////////////////////////////////////////////// */

void reset_layers (MOLECULE *molecule)
{
  int i;

  for (i = 0; i < molecule->max.layers; i++)
    reset_layer (molecule, i);

  molecule->total.layers = 0;
}

/* ///////////////////////////////////////////////////////////// */

void reset_layer (MOLECULE *molecule, int layer)
{
  int segment;
  int segment_id;

  for
  (
    segment_id = 0;
    segment_id < molecule->layer[layer].segment_total;
    segment_id++
  )
  {
    segment = molecule->layer[layer].segment[segment_id];
    molecule->segment[segment].layer_id = NEITHER;
  }

  molecule->layer[layer].id = NEITHER;
  molecule->layer[layer].active_flag = FALSE;
  molecule->layer[layer].conform_total = 1;
  reset_score (&molecule->layer[layer].score);

  molecule->layer[layer].segment_total = 0;
  memset
  (
    molecule->layer[layer].segment,
    0,
    molecule->layer[layer].segment_max * sizeof (int)
  );
}

/* ///////////////////////////////////////////////////////////// */

void free_layers (MOLECULE *molecule)
{
  int i;

  for (i = 0; i < molecule->max.layers; i++)
    free_layer_segments (&molecule->layer[i]);

  efree ((void **) &molecule->layer);
  molecule->max.layers = 0;
}

/* ///////////////////////////////////////////////////////////// */

void free_layer_segments (LAYER *layer)
{
  efree ((void **) &layer->segment);
  layer->segment_max = 0;
}

/* ///////////////////////////////////////////////////////////// */

void reallocate_layers (MOLECULE *molecule)
{
  if (molecule->total.layers > molecule->max.layers)
  {
    free_layers (molecule);
    molecule->max.layers = molecule->total.layers;
    allocate_layers (molecule);
  }
}

/* ///////////////////////////////////////////////////////////// */

void reallocate_layer_segments (LAYER *layer)
{
  if (layer->segment_total > layer->segment_max)
  {
    free_layer_segments (layer);
    layer->segment_max = layer->segment_total;
    allocate_layer_segments (layer);
  }
}

/* ///////////////////////////////////////////////////////////// */

void copy_layers (MOLECULE *copy, MOLECULE *original)
{
  int i, j;
  int segment;

  copy->total.layers = original->total.layers;
  reallocate_layers (copy);

  for (i = 0; i < original->total.layers; i++)
  {
    copy_layer (&copy->layer[i], &original->layer[i]);

    for (j = 0; j < copy->layer[i].segment_total; j++)
    {
      segment = copy->layer[i].segment[j];

      if ((segment >= 0) && (segment < copy->total.segments))
        copy->segment[segment].layer_id = i;

      else
        exit (fprintf (global.outfile,
          "ERROR copy_layers: unable to update segments\n"));
    }
  }
}

/* ///////////////////////////////////////////////////////////// */

void copy_layer (LAYER *copy, LAYER *original)
{
  copy->id = original->id;
  copy->active_flag = original->active_flag;
  copy->conform_total = original->conform_total;
  copy_score (&copy->score, &original->score);

  copy->segment_total = original->segment_total;
  reallocate_layer_segments (copy);
  memcpy
    (copy->segment, original->segment,
    original->segment_total * sizeof (int));
}

/* ///////////////////////////////////////////////////////////// */

void save_layers (MOLECULE *molecule, FILE *file)
{
  int i;

  for (i = 0; i < molecule->max.layers; i++)
    save_layer (&molecule->layer[i], file);
}

/* ///////////////////////////////////////////////////////////// */

void save_layer (LAYER *layer, FILE *file)
{
  efwrite (&layer->id, sizeof (int), 1, file);
  efwrite (&layer->active_flag, sizeof (int), 1, file);
  efwrite (&layer->conform_total, sizeof (int), 1, file);
  save_score (&layer->score, file);

  efwrite (&layer->segment_total, sizeof (int), 1, file);
  efwrite (&layer->segment_max, sizeof (int), 1, file);
  efwrite (layer->segment, sizeof (int), layer->segment_max, file);
}

/* ///////////////////////////////////////////////////////////// */

void load_layers (MOLECULE *molecule, FILE *file)
{
  int i;

  for (i = 0; i < molecule->max.layers; i++)
    load_layer (&molecule->layer[i], file);
}

/* ///////////////////////////////////////////////////////////// */

void load_layer (LAYER *layer, FILE *file)
{
  efread (&layer->id, sizeof (int), 1, file);
  efread (&layer->active_flag, sizeof (int), 1, file);
  efread (&layer->conform_total, sizeof (int), 1, file);
  load_score (&layer->score, file);

  efread (&layer->segment_total, sizeof (int), 1, file);
  efread (&layer->segment_max, sizeof (int), 1, file);
  allocate_layer_segments (layer);
  efread (layer->segment, sizeof (int), layer->segment_max, file);
}


/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

void allocate_sets (MOLECULE *molecule)
{
  int i;

  if (molecule->max.sets > 0)
    ecalloc
    (
      (void **) &molecule->set,
      molecule->max.sets,
      sizeof (SET),
      molecule->info.name,
      global.outfile
    );

  for (i = 0; i < molecule->max.sets; i++)
    reset_set (&molecule->set[i]);
}

/* ///////////////////////////////////////////////////////////// */

void reset_sets (MOLECULE *molecule)
{
  int i;

  for (i = 0; i < molecule->max.sets; i++)
    reset_set (&molecule->set[i]);

  molecule->total.sets = 0;
}

/* ///////////////////////////////////////////////////////////// */

void reset_set (SET *set)
{
  efree ((void **) &set->name);
  efree ((void **) &set->type);
  efree ((void **) &set->obj_type);
  efree ((void **) &set->sub_type);
  efree ((void **) &set->status);
  efree ((void **) &set->comment);
  efree ((void **) &set->member);

  memset (set, 0, sizeof (SET));
}

/* ///////////////////////////////////////////////////////////// */

void free_sets (MOLECULE *molecule)
{
  int i;

  for (i = 0; i < molecule->max.sets; i++)
    free_set (&molecule->set[i]);

  efree ((void **) &molecule->set);
  molecule->max.sets = 0;
}

/* ///////////////////////////////////////////////////////////// */

void free_set (SET *set)
{
  efree ((void **) &set->name);
  efree ((void **) &set->type);
  efree ((void **) &set->obj_type);
  efree ((void **) &set->sub_type);
  efree ((void **) &set->status);
  efree ((void **) &set->comment);
  efree ((void **) &set->member);

  set->member_total = 0;
}

/* ///////////////////////////////////////////////////////////// */

void reallocate_sets (MOLECULE *molecule)
{
  if (molecule->total.sets > molecule->max.sets)
  {
    free_sets (molecule);
    molecule->max.sets = molecule->total.sets;
    allocate_sets (molecule);
  }
}

/* ///////////////////////////////////////////////////////////// */

void copy_sets (MOLECULE *copy, MOLECULE *original)
{
  int i;

  copy->total.sets = original->total.sets;
  reallocate_sets (copy);

  for (i = 0; i < original->total.sets; i++)
    copy_set (&copy->set[i], &original->set[i]);

  copy->total.sets = original->total.sets;
}

/* ///////////////////////////////////////////////////////////// */

void copy_set (SET *copy, SET *original)
{
  int i;

  vstrcpy (&copy->name, original->name);
  vstrcpy (&copy->type, original->type);
  vstrcpy (&copy->obj_type, original->obj_type);
  vstrcpy (&copy->sub_type, original->sub_type);
  vstrcpy (&copy->status, original->status);
  vstrcpy (&copy->comment, original->comment);

  if (copy->member_total != original->member_total)
  {
    efree ((void **) &copy->member);
    ecalloc
    (
      (void **) &copy->member,
      original->member_total,
      sizeof (int),
      "set members",
      global.outfile
    );
    copy->member_total = original->member_total;
  }

  for (i = 0; i < original->member_total; i++)
    copy->member[i] = original->member[i];
}

/* ///////////////////////////////////////////////////////////// */

void save_sets (MOLECULE *molecule, FILE *file)
{
  int i;

  for (i = 0; i < molecule->max.sets; i++)
    save_set (&molecule->set[i], file);
}

/* ///////////////////////////////////////////////////////////// */

void save_set (SET *set, FILE *file)
{
  efwrite (set, sizeof (SET), 1, file);

  save_string (&set->name, file);
  save_string (&set->type, file);
  save_string (&set->obj_type, file);
  save_string (&set->sub_type, file);
  save_string (&set->status, file);
  save_string (&set->comment, file);

  efwrite (set->member, sizeof (int), set->member_total, file);
}

/* ///////////////////////////////////////////////////////////// */

void load_sets (MOLECULE *molecule, FILE *file)
{
  int i;

  for (i = 0; i < molecule->max.sets; i++)
    load_set (&molecule->set[i], file);
}

/* ///////////////////////////////////////////////////////////// */

void load_set (SET *set, FILE *file)
{
  efread (set, sizeof (SET), 1, file);

  set->name = NULL;
  set->type = NULL;
  set->obj_type = NULL;
  set->sub_type = NULL;
  set->status = NULL;
  set->comment = NULL;
  set->member = NULL;

  load_string (&set->name, file);
  load_string (&set->type, file);
  load_string (&set->obj_type, file);
  load_string (&set->sub_type, file);
  load_string (&set->status, file);
  load_string (&set->comment, file);

  emalloc
  (
    (void **) &set->member,
    set->member_total * sizeof (int),
    "loaded set members",
    global.outfile
  );

  efread (set->member, sizeof (int), set->member_total, file);
}

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

void allocate_keys (MOLECULE *molecule)
{
  int i;

  if (molecule->max.keys > 0)
  {
    ecalloc
    (
      (void **) &molecule->key,
      molecule->max.keys,
      sizeof (KEY *),
      molecule->info.name,
      global.outfile
    );

    for (i = 0; i < molecule->max.keys; i++)
      ecalloc
      (
        (void **) &molecule->key[i],
        molecule->max.keys,
        sizeof (KEY),
        molecule->info.name,
        global.outfile
      );
  }
}

/* ///////////////////////////////////////////////////////////// */

void reset_keys (MOLECULE *molecule)
{
  int i, j;

  for (i = 0; i < molecule->max.keys; i++)
    for (j = 0; j < molecule->max.keys; j++)
    {
      molecule->key[i][j].count = 0;
      molecule->key[i][j].distance = (MASK) 0;
    }

  molecule->total.keys = 0;
}

/* ///////////////////////////////////////////////////////////// */

void free_keys (MOLECULE *molecule)
{
  int i;

  for (i = 0; i < molecule->max.keys; i++)
    efree ((void **) &molecule->key[i]);

  efree ((void **) &molecule->key);
  molecule->max.keys = 0;
}

/* ///////////////////////////////////////////////////////////// */

void reallocate_keys (MOLECULE *molecule)
{
  if (molecule->total.keys > molecule->max.keys)
  {
    free_keys (molecule);
    molecule->max.keys = molecule->total.keys;
    allocate_keys (molecule);
  }
}

/* ///////////////////////////////////////////////////////////// */

void copy_keys (MOLECULE *copy, MOLECULE *original)
{
  int i;

  copy->total.keys = original->total.keys;
  reallocate_keys (copy);

  for (i = 0; i < original->total.keys; i++)
    memcpy
    (
      copy->key[i],
      original->key[i],
      original->total.keys * sizeof (KEY)
    );
}

/* ///////////////////////////////////////////////////////////// */

void save_keys (MOLECULE *molecule, FILE *file)
{
  int i;

  for (i = 0; i < molecule->max.keys; i++)
    efwrite (molecule->key[i], sizeof (KEY), molecule->max.keys, file);
}

/* ///////////////////////////////////////////////////////////// */

void load_keys (MOLECULE *molecule, FILE *file)
{
  int i;

  for (i = 0; i < molecule->max.keys; i++)
    efread (molecule->key[i], sizeof (KEY), molecule->max.keys, file);
}

/* //////////////////////////////////////////////////////////// */
/*
Routine to save the contents of a dynamically allocated character
string.
2/96 te
*/

void save_string (char **string, FILE *file)
{
  int length;

/*
* If the string has been allocated, then find the length + 1 (for
* the null character)
* 2/96 te
*/
  if (*string)
    length = strlen (*string) + 1;

  else
    length = 0;

  efwrite (&length, sizeof (int), 1, file);

  if (length > 0)
    efwrite (*string, sizeof (char), length, file);
}

/* //////////////////////////////////////////////////////////// */
/*
Routine to load the contents of a dynamically allocated character
string.
2/96 te
*/

void load_string (char **string, FILE *file)
{
  int length;

  efread (&length, sizeof (int), 1, file);

  if (length > 0)
  {
    ecalloc
    (
      (void **) string,
      length,
      sizeof (char),
      "molecule string",
      global.outfile
    );

    efread (*string, sizeof (char), length, file);
  }
}


/* ///////////////////////////////////////////////////////////// */

void atom_neighbors (MOLECULE *molecule)
{
  int i;
  int origin;
  int target;
  int neighbor;

  for (i = 0; i < molecule->total.atoms; i++)
  {
    molecule->atom[i].neighbor_total = 0;
    molecule->atom[i].flag = TRUE;
  }

  if (molecule->total.bonds > 0)
  {
/*
*   Compute the total neighbors for each atom
*   12/96 te
*/
    for (i = 0; i < molecule->total.bonds; i++)
    {
      molecule->atom[molecule->bond[i].origin].neighbor_total++;
      molecule->atom[molecule->bond[i].target].neighbor_total++;
    }

/*
*   Allocate space for atom neighbor lists and reset totals
*   12/96 te
*/
    for (i = 0; i < molecule->total.atoms; i++)
    {
      reallocate_atom_neighbors (&molecule->atom[i]);
      molecule->atom[i].neighbor_total = 0;
    }

/*
*   Record bond information into atom neighbor lists
*   12/96 te
*/
    for (i = 0; i < molecule->total.bonds; i++)
    {
      origin = molecule->bond[i].origin;
      target = molecule->bond[i].target;

      if ((origin < 0) || (origin >= molecule->total.atoms))
        exit (fprintf (global.outfile,
          "ERROR atom_neighbors: origin atom incorrect for bond %d\n", i + 1));

      if ((target < 0) || (target >= molecule->total.atoms))
        exit (fprintf (global.outfile,
          "ERROR atom_neighbors: target atom incorrect for bond %d\n", i + 1));

      neighbor = molecule->atom[origin].neighbor_total;
      molecule->atom[origin].neighbor[neighbor].id = target;
      molecule->atom[origin].neighbor[neighbor].bond_id = i;

      if (molecule->atom[target].flag == TRUE)
        molecule->atom[origin].neighbor[neighbor].out_flag = TRUE;

      molecule->atom[target].flag = FALSE;
      molecule->atom[origin].neighbor_total++;

      neighbor = molecule->atom[target].neighbor_total;
      molecule->atom[target].neighbor[neighbor].id = origin;
      molecule->atom[target].neighbor[neighbor].bond_id = i;
      molecule->atom[target].neighbor[neighbor].out_flag = FALSE;
      molecule->atom[target].neighbor_total++;
    }
  }
}


/* ///////////////////////////////////////////////////////////// */

void flag_atom_neighbor
(
  void	*atom,
  int	atom_id,
  int	neighbor_id,
  int	flag
)
{
  ((ATOM *) atom)[atom_id].neighbor[neighbor_id].out_flag = flag;
}


/* ///////////////////////////////////////////////////////////// */

void flag_segment_neighbor
(
  void	*segment,
  int	segment_id,
  int	neighbor_id,
  int	flag
)
{
  ((SEGMENT *) segment)[segment_id].neighbor[neighbor_id].out_flag = flag;
}


/* ///////////////////////////////////////////////////////////// */

void revise_atom_neighbors (MOLECULE *molecule)
{
  SEARCH search = {0};
  int i, j;
  int bond;

/*
* Perform breadth-first atom search, starting from anchor segment.
* Also, flag atom neighbors during this search run.
* 11/96 te
*/
  for
  (
    i = 0;
    breadth_search
    (
      &search,
      molecule->atom,
      molecule->total.atoms,
      get_atom_neighbor,
      flag_atom_neighbor,
      &molecule->transform.anchor_atom, 1,
      NEITHER, i
    ) != EOF;
    i++
  );

/*
* Update the origin and target entries in bond records
* 1/97 te
*/
  for (i = 0; i < molecule->total.atoms; i++)
    for (j = 0; j < molecule->atom[i].neighbor_total; j++)
      if (molecule->atom[i].neighbor[j].out_flag == TRUE)
      {
        bond = molecule->atom[i].neighbor[j].bond_id;

        if (molecule->bond[bond].target == i)
        {
          molecule->bond[bond].origin = i;
          molecule->bond[bond].target = molecule->atom[i].neighbor[j].id;
        }
      }

  free_search (&search);
}

/* //////////////////////////////////////////////////////////////////////

Subroutine to compute the centrality of an atom.  This is done by performing
a breadth-first search to identify the bond distance of the most distant atom
from each atom.

11/96 te

////////////////////////////////////////////////////////////////////// */

int get_atom_centrality (MOLECULE *molecule, int atom_id)
{
  int i;
  static SEARCH search = {0};

  if (molecule->atom[atom_id].centrality != NEITHER)
    return molecule->atom[atom_id].centrality;

  else
  {
    for
    (
      i = molecule->atom[atom_id].centrality = 0;
      breadth_search
      (
        &search,
        molecule->atom,
        molecule->total.atoms,
        get_atom_neighbor,
        NULL,
        &atom_id, 1,
        NEITHER,
        i
      ) != EOF;
      i++
    )
      molecule->atom[atom_id].centrality += search.radius;

    return molecule->atom[atom_id].centrality;
/*
    return
      molecule->atom[atom_id].centrality =
        get_search_radius (&search, NEITHER, NEITHER);
*/
  }
}


/* //////////////////////////////////////////////////////////////////////

Subroutine that determines the neighboring atoms to the origin and target
atoms for each torsion.

11/96 te

////////////////////////////////////////////////////////////////////// */

void get_torsion_neighbors (MOLECULE *molecule)
{
  int i, j, k;
  int bond_id;
  int origin;
  int target;
  int neighbor;
  int central_neighbor;
  int central;
  static SEARCH search = {0};

  for (i = 0; i < molecule->total.torsions; i++)
  {
    bond_id = molecule->torsion[i].bond_id;
    origin = molecule->torsion[i].origin = molecule->bond[bond_id].origin;
    target = molecule->torsion[i].target = molecule->bond[bond_id].target;

/*
*   Determine which neighbor of the origin atom is most central
*   11/96 te
*/
    if (molecule->atom[origin].neighbor_total > 1)
    {
      for
      (
        j = central_neighbor = 0, central = INT_MIN;
        (neighbor = get_atom_neighbor (molecule->atom, origin, j)) != EOF;
        j++
      )
      {
        if (neighbor != target)
        {
          for
          (
            k = 0;
            breadth_search
            (
              &search,
              molecule->atom,
              molecule->total.atoms,
              get_atom_neighbor,
              NULL,
              &neighbor, 1,
              origin,
              k
            ) != EOF;
            k++
          );

          if (search.radius > central)
          {
            central = search.radius;
            central_neighbor = neighbor;
          }
        }
      }

      molecule->torsion[i].origin_neighbor = central_neighbor;
    }

    else
      molecule->torsion[i].origin_neighbor = NEITHER;

/*
*   Determine which neighbor of the target atom is most central
*   11/96 te
*/
    if (molecule->atom[target].neighbor_total > 1)
    {
      for
      (
        j = central_neighbor = 0, central = INT_MIN;
        (neighbor = get_atom_neighbor (molecule->atom, target, j)) != EOF;
        j++
      )
      {
        if (neighbor != origin)
        {
          for
          (
            k = 0;
            breadth_search
            (
              &search,
              molecule->atom,
              molecule->total.atoms,
              get_atom_neighbor,
              NULL,
              &neighbor, 1,
              target,
              k
            ) != EOF;
            k++
          );

          if (search.radius > central)
          {
            central = search.radius;
            central_neighbor = neighbor;
          }
        }
      }

      molecule->torsion[i].target_neighbor = central_neighbor;
    }

    else
      molecule->torsion[i].target_neighbor = NEITHER;
  }
}


/* //////////////////////////////////////////////////////////////////////

Subroutine that switches the origin and target atoms in a torsion.

11/96 te

////////////////////////////////////////////////////////////////////// */

void reverse_torsion (TORSION *torsion)
{
  int atom;

  atom = torsion->origin;
  torsion->origin = torsion->target;
  torsion->target = atom;

  atom = torsion->origin_neighbor;
  torsion->origin_neighbor = torsion->target_neighbor;
  torsion->target_neighbor = atom;
}
 

/* //////////////////////////////////////////////////////////////////////

Recursive subroutine that counts up the number of atoms directly
and indirectly bonded to a given atom.

10/95 te

////////////////////////////////////////////////////////////////////// */

int bonded_atoms
(
  MOLECULE *molecule,
  int current_atom
)
{
  int i;
  int total;		/* Number of atoms seen on linkage path */

  if (molecule->atom[current_atom].flag)
    return 0;

  else
  {
    molecule->atom[current_atom].flag = TRUE;
    total = 1;

    for (i = 0; i < molecule->atom[current_atom].neighbor_total; i++)
      total +=
        bonded_atoms
        (
          molecule,
          molecule->atom[current_atom].neighbor[i].id
        );

    return total;
  }
}


/* //////////////////////////////////////////////////////////////

Routine to construct a distance matrix.
3/97 te

////////////////////////////////////////////////////////////// */

float calculate_distances (MOLECULE *molecule, float ***distances, int *size)
{
  int i, j;
  float distance_max;

/*
* Either allocate space for distance matrix, or reset previous space
* 3/97 te
*/
  if (*size < molecule->total.atoms)
  {
    free_distances (distances, size);

    *size = molecule->total.atoms;

    ecalloc
    (
      (void **) distances,
      *size,
      sizeof (float *),
      "distance matrix",
      global.outfile
    );

    for (i = 0; i < *size; i++)
      ecalloc
      (
        (void **) &(*distances)[i],
        *size,
        sizeof (float),
        "distance matrix",
        global.outfile
      );
  }

  else
  {
    for (i = 0; i < *size; i++)
      memset ((*distances)[i], 0, *size * sizeof (float));
  }

/*
* Compute the distance between all pairs of atoms
* 3/97 te
*/
  for (i = distance_max = 0; i < molecule->total.atoms; i++)
  {
    for (j = i + 1; j < molecule->total.atoms; j++)
    {
      (*distances)[i][j] =
        (*distances)[j][i] =
          dist3 (molecule->coord[i], molecule->coord[j]);

      if ((*distances)[i][j] > distance_max)
        distance_max = (*distances)[i][j];
    }
  }

  return distance_max;
}

/* //////////////////////////////////////////////////////////////

Routine to free distance matrix.
5/97 te

////////////////////////////////////////////////////////////// */

void free_distances (float ***distances, int *size)
{
  int i;

  for (i = 0; i < *size; i++)
    efree ((void **) &(*distances)[i]);

  efree ((void **) distances);
  *size = 0;
}

