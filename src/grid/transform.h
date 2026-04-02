/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/

#include "extsymbols.h"

/*
Routines defined in transform.c, that are called by outside functions
*/

int orient_molecule
(
  MOLECULE *target,
  MOLECULE *current,
  MOLECULE *mol_ref,
  MOLECULE *mol_conf,
  MOLECULE *mol_ori
);

void transform_molecule
(
  MOLECULE *mol_ori,
  MOLECULE *mol_ref
);

void rigid_transform
(
  MOLECULE *mol_ori,
  MOLECULE *mol_ref
);

void torsion_transform
(
  MOLECULE *molecule,
  int torsion_id
);

void rotate_atoms
(
  MOLECULE *molecule,
  XYZ matrix[3],
  int origin,
  int atom
);

int get_matrix_from_quaternion (XYZ [3], XYZ, int);
int get_quaternion_from_matrix (XYZ, int *, XYZ [3]);
float compute_torsion (MOLECULE *, int);
void overall_rotation (MOLECULE *, MOLECULE *);
void overall_translation (MOLECULE *, MOLECULE *);

/*
Fortran routines in transformf.f
*/
int EXTERNAL_ORIENT_GK(int *, XYZ *, XYZ *, XYZ, XYZ *, XYZ, int *);
void EXTERNAL_TRANSFORM(int *, XYZ *, XYZ, XYZ *, XYZ, XYZ *);
void EXTERNAL_TRANSFORM_ATOM(XYZ, XYZ *, XYZ);

