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
#include "transform.h"


/* ///////////////////////////////////////////////////////////// */

int orient_molecule
(
  MOLECULE *target,
  MOLECULE *current,
  MOLECULE *mol_ref,
  MOLECULE *mol_conf,
  MOLECULE *mol_ori
)
{
  static XYZ rotation_matrix[3];

/*
* Calculate rotation matrix and translation vector.
* 6/95 te
*/
  if (!EXTERNAL_ORIENT_GK
  (
    &target->total.atoms,
    target->coord,
    current->coord,
    mol_ori->transform.com,
    rotation_matrix,
    mol_ori->transform.translate,
    &mol_ori->transform.refl_flag
  ))
    return FALSE;

/*
* Update transform records
* 10/96 te
*/
  get_quaternion_from_matrix
  (
    mol_ori->transform.rotate,
    &mol_ori->transform.refl_flag,
    rotation_matrix
  );

  mol_ori->transform.trans_flag = TRUE;
  mol_ori->transform.rot_flag = TRUE;

/*
* If a previous transformation is present, then update this info
* 10/96 te
*/
  if (mol_conf->transform.rot_flag)
    overall_rotation (mol_ori, mol_conf);

  if (mol_conf->transform.trans_flag)
    overall_translation (mol_ori, mol_conf);

/*
* Transform the molecule coordinates
* 10/96 te
*/
  transform_molecule (mol_ori, mol_ref);

  return TRUE;
}


/* ///////////////////////////////////////////////////////////// */

void transform_molecule
(
  MOLECULE *mol_trans,
  MOLECULE *mol_ref
)
{
  int i;

  mol_trans->transform.flag = TRUE;

/*
* Perform a rigid-body rotation/translation
* 12/96 te
*/
  if ((mol_trans->transform.trans_flag == TRUE) ||
    (mol_trans->transform.rot_flag == TRUE))
    rigid_transform (mol_trans, mol_ref);

/*
* Perform a rotation about each rotatable bond vector
* 12/96 te
*/
  if (mol_trans->transform.tors_flag == TRUE)
    for (i = 0; i < mol_trans->total.torsions; i++)
      torsion_transform (mol_trans, i);
}


/* ///////////////////////////////////////////////////////////// */

void rigid_transform
(
  MOLECULE *mol_trans,
  MOLECULE *mol_ref
)
{
  static XYZ rotation[3];

  mol_trans->transform.flag = TRUE;

/*
* Calculate a rotation matrix from the quaternion
* 12/96 te
*/
  get_matrix_from_quaternion
  (
    rotation,
    mol_trans->transform.rotate,
    mol_trans->transform.refl_flag
  );

  EXTERNAL_TRANSFORM
  (
    &mol_trans->total.atoms,
    mol_ref->coord,
    mol_trans->transform.com,
    rotation,
    mol_trans->transform.translate,
    mol_trans->coord
  );
}


/* ///////////////////////////////////////////////////////////// */

void torsion_transform
(
  MOLECULE	*molecule,
  int		torsion_id
)
{
  int i;
  static XYZ bond_matrix[3];
  static XYZ bond_vector;
  static float difference;

  void rotateaxis (float, XYZ, XYZ *);

  molecule->transform.flag = TRUE;

/*
* Compute current torsion angle and difference
* 12/96 te
*/
  molecule->torsion[torsion_id].current_angle =
    compute_torsion (molecule, torsion_id);

  difference =
    molecule->torsion[torsion_id].target_angle -
    molecule->torsion[torsion_id].current_angle;

/*
* If difference is significant, then rotate bond
* 12/96 te
*/
  if (NON_ZERO (difference))
  {
    for (i = 0; i < 3; i++)
      bond_vector[i] =
        molecule->coord[molecule->torsion[torsion_id].target][i] -
        molecule->coord[molecule->torsion[torsion_id].origin][i];

    rotateaxis (difference, bond_vector, bond_matrix);

    rotate_atoms
    (
      molecule,
      bond_matrix,
      molecule->torsion[torsion_id].target,
      molecule->torsion[torsion_id].target
    );

    molecule->torsion[torsion_id].current_angle =
      molecule->torsion[torsion_id].target_angle;
  }
}


/* //////////////////////////////////////////////////////////////////////

Subroutine that calculate a rotation matrix, given a vector axis
and an angle of rotation.

10/95 ys

////////////////////////////////////////////////////////////////////// */

void rotateaxis (float phi, XYZ vect, XYZ rot[3])
{
  float del1,del2,del3;
  float vect_length;
  float vcosx,vcosy,vcosz;

  del1=1-cos(phi);
  del2=sin(phi);
  del3=cos(phi);

  vect_length=sqrt(vect[0]*vect[0]+
                   vect[1]*vect[1]+
                   vect[2]*vect[2]);
  vcosx= vect[0]/vect_length;
  vcosy= vect[1]/vect_length;
  vcosz= vect[2]/vect_length;

  rot[0][0]=vcosx*vcosx*del1 + del3;
  rot[0][1]=vcosx*vcosy*del1 + vcosz*del2;
  rot[0][2]=vcosx*vcosz*del1 - vcosy*del2;
  rot[1][0]=vcosy*vcosx*del1 - vcosz*del2;
  rot[1][1]=vcosy*vcosy*del1 + del3;
  rot[1][2]=vcosy*vcosz*del1 + vcosx*del2;
  rot[2][0]=vcosz*vcosx*del1 + vcosy*del2;
  rot[2][1]=vcosz*vcosy*del1 - vcosx*del2;
  rot[2][2]=vcosz*vcosz*del1 + del3;

  return;
}

/* /////////////////////////////////////////////////////////////////// */

void rotate_atoms
(
  MOLECULE	*molecule,
  XYZ		bond_matrix[3],
  int		origin,
  int		atom
)
{
  int i;

/*
  fprintf (global.outfile, "orig %s atom %s\n",
    molecule->atom[origin].name,
    molecule->atom[atom].name);
*/

  EXTERNAL_TRANSFORM_ATOM
  (
    molecule->coord[atom],
    bond_matrix,
    molecule->coord[origin]
  );

  for (i = 0; i < molecule->atom[atom].neighbor_total; i++)
    if (molecule->atom[atom].neighbor[i].out_flag == TRUE)
      rotate_atoms
        (molecule, bond_matrix, origin, molecule->atom[atom].neighbor[i].id);
}


/* /////////////////////////////////////////////////////////////////// */
/*
Function to compute a rotation matrix from quaternion values.
For discussion of this technique, see "Computer Simulation of Liquids"
by M.P. Allen and D.J. Tildesley from Oxford Science Publishers, 1987
10/96 te
*/

int get_matrix_from_quaternion
(
  XYZ	m[3],		/* rotation matrix */
  XYZ	qin,		/* input independent quaternion elements */
  int	refl		/* reflection flag */
)
{
  int	i;		/* iterater */
  float	qn;		/* dependent quaternion element, q-naught */
  float	qn2;		/* square of q-naught */
  XYZ	q;		/* independent quaternion elements */
  XYZ	q2;		/* square of independent quaternion elements */
  float	sum2;		/* sum of squares of independent q values */
  float	sum;		/* sum of independent q values */

/*
* Check that each q-value is between -1.0 and 1.0.
* If not, remap into this range using wrap-around.
* Compute q-squared values and their sum
* 10/96 te
*/
  for (i = 0, sum2 = 0.0; i < 3; i++)
  {
    q[i] = qin[i];

    if (q[i] > 1.0)
      q[i] = fmod (q[i] + 1.0, 2.0) - 1.0;

    else if (q[i] < -1.0)
      q[i] = fmod (q[i] - 1.0, 2.0) + 1.0;

    sum2 += q2[i] = SQR (q[i]);
  }

/*
* If the sum-of-squares is less than 1.0, compute q-naught
* 10/96 te
*/
  if (sum2 < 1.0)
  {
    qn2 = 1.0 - sum2;
    qn = sqrt (qn2);
  }

/*
* If the sum is 1.0, set q-naught to zero
* 10/96 te
*/
  else if (sum2 == 1.0)
    qn = qn2 = 0.0;

/*
* Otherwise, renormalize q-values and set q-naught to zero
* 10/96 te
*/
  else
  {
    for (i = 0; i < 3; i++)
      q2[i] /= sum2;

    for (i = 0, sum = sqrt (sum2); i < 3; i++)
      q[i] /= sum;

    qn = qn2 = 0.0;
  }

/*
* Compute rotation matrix elements
* 10/96 te
*/
  m[0][0] = qn2 + q2[0] - q2[1] - q2[2];
  m[0][1] = 2.0 * (q[0] * q[1] + qn * q[2]);
  m[0][2] = 2.0 * (q[0] * q[2] - qn * q[1]);

  m[1][0] = 2.0 * (q[0] * q[1] - qn * q[2]);
  m[1][1] = qn2 - q2[0] + q2[1] - q2[2];
  m[1][2] = 2.0 * (q[1] * q[2] + qn * q[0]);

  m[2][0] = 2.0 * (q[0] * q[2] + qn * q[1]);
  m[2][1] = 2.0 * (q[1] * q[2] - qn * q[0]);
  m[2][2] = qn2 - q2[0] - q2[1] + q2[2];

/*
* If a reflection is required, then update the matrix
* 10/96 te
*/
  if (refl == TRUE)
    for (i = 0; i < 3; i++)
      m[i][2] *= -1.0;

  return TRUE;
}

/* /////////////////////////////////////////////////////////////////// */
/*
Function to compute the quaternion values from a rotation matrix
For discussion of this technique, see "Computer Simulation of Liquids"
by M.P. Allen and D.J. Tildesley from Oxford Science Publishers, 1987
10/96 te
*/

int get_quaternion_from_matrix
(
  XYZ	q,		/* independent quaternion elements */
  int	*refl,		/* reflection flag */
  XYZ	m[3]		/* rotation matrix */
)
{
  int	i;		/* iterater */
  int	max;		/* quaternion element with greatest magnitude */
  float	max_val;	/* value of greatest quaternion element */
  float	qn;		/* dependent quaternion element, q-naught */
  float	qn2;		/* square of q-naught */
  XYZ	q2;		/* square of independent quaternion elements */
  float	det;		/* determinant of input matrix */

/*
* Compute determinant of matrix
* 10/96 te
*/
  det =
    m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
    m[0][1] * (m[1][0] * m[2][2] - m[2][0] * m[1][2]) +
    m[0][2] * (m[1][0] * m[2][1] - m[2][0] * m[1][1]);

/*
* If the determinant is negative, then set the reflection flag and
* invert matrix
* 10/96 te
*/
  if (det < 0.0)
  {
    *refl = TRUE;

    for (i = 0; i < 3; i++)
      m[i][2] *= -1.0;
  }

  else if (*refl != NEITHER)
    *refl = FALSE;

/*
* Calculate q-squared values (times 4)
* 10/96 te
*/
  qn2   = 1.0 + m[0][0] + m[1][1] + m[2][2];
  q2[0] = 1.0 + m[0][0] - m[1][1] - m[2][2];
  q2[1] = 1.0 - m[0][0] + m[1][1] - m[2][2];
  q2[2] = 1.0 - m[0][0] - m[1][1] + m[2][2];

/*
* Identify the largest q-squared value
* 10/96 te
*/
  for (i = max = 0, max_val = qn2; i < 3; i++)
    if (q2[i] > max_val)
    {
      max_val = q2[i];
      max = i + 1;
    }

/*
* Compute the other q-values based on the positive root
* of the largest q-squared value
* 10/96 te
*/
  switch (max)
  {
    case 0:
    {
      qn   = 0.5 * sqrt (qn2);

      q[0] = 0.25 * (m[1][2] - m[2][1]) / qn;
      q[1] = 0.25 * (m[2][0] - m[0][2]) / qn;
      q[2] = 0.25 * (m[0][1] - m[1][0]) / qn;

      break;
    }

    case 1:
    {
      q[0] = 0.5 * sqrt (q2[0]);

      qn   = 0.25 * (m[1][2] - m[2][1]) / q[0];
      q[1] = 0.25 * (m[0][1] + m[1][0]) / q[0];
      q[2] = 0.25 * (m[0][2] + m[2][0]) / q[0];

      break;
    }

    case 2:
    {
      q[1] = 0.5 * sqrt (q2[1]);

      qn   = 0.25 * (m[2][0] - m[0][2]) / q[1];
      q[0] = 0.25 * (m[1][0] + m[0][1]) / q[1];
      q[2] = 0.25 * (m[2][1] + m[1][2]) / q[1];

      break;
    }

    case 3:
    {
      q[2] = 0.5 * sqrt (q2[2]);

      qn   = 0.25 * (m[0][1] - m[1][0]) / q[2];
      q[0] = 0.25 * (m[2][0] + m[0][2]) / q[2];
      q[1] = 0.25 * (m[2][1] + m[1][2]) / q[2];

      break;
    }

    default:
    {
      printf ("Unable to find q2 value > 0\n");
      return FALSE;
    }
  }

/*
* Reset q-values so that q-naught is always positive
* 10/96 te
*/
  if (qn < 0.0)
    for (i = 0; i < 3; i++)
      q[i] *= -1.0;

  return TRUE;
}


/* /////////////////////////////////////////////////////////////////
Compute current angle of specified bond
10/96 te
///////////////////////////////////////////////////////////////// */

float compute_torsion (MOLECULE *molecule, int torsion_id)
{
  if
  (
    (molecule->torsion[torsion_id].origin_neighbor == NEITHER) ||
    (molecule->torsion[torsion_id].origin == NEITHER) ||
    (molecule->torsion[torsion_id].target == NEITHER) ||
    (molecule->torsion[torsion_id].target_neighbor == NEITHER)
  )
    exit (fprintf (global.outfile,
      "ERROR compute_torsion: Torsion atoms undefined in %s\n",
      molecule->info.name));

  return
    dihed
    (
      molecule->coord[molecule->torsion[torsion_id].origin_neighbor],
      molecule->coord[molecule->torsion[torsion_id].origin],
      molecule->coord[molecule->torsion[torsion_id].target],
      molecule->coord[molecule->torsion[torsion_id].target_neighbor]
    );
}


/* /////////////////////////////////////////////////////////////////
Compute a new overall rotation
10/96 te
////////////////////////////////////////////////////////////////// */

void overall_rotation
(
  MOLECULE *mol_new,
  MOLECULE *mol_orig
)
{
  static XYZ rot_orig[3], rot_new[3], rot_temp[3];

  get_matrix_from_quaternion
    (rot_orig, mol_orig->transform.rotate, mol_orig->transform.refl_flag);
  get_matrix_from_quaternion
    (rot_new, mol_new->transform.rotate, mol_new->transform.refl_flag);
  mmult3 (rot_orig, rot_new, rot_temp);
  get_quaternion_from_matrix
    (mol_new->transform.rotate, &mol_new->transform.refl_flag, rot_temp);
}


/* /////////////////////////////////////////////////////////////////
Compute a new overall translation
10/96 te
////////////////////////////////////////////////////////////////// */

void overall_translation
(
  MOLECULE *mol_new,
  MOLECULE *mol_orig
)
{
  int i;

  for (i = 0; i < 3; i++)
    mol_new->transform.translate[i] += mol_orig->transform.translate[i];
}

