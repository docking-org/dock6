/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

#include "define.h"
#include "vector.h"
#include "rotrans.h"
float x_axis[] = {1, 0, 0};
float y_axis[] = {0, 1, 0};
float z_axis[] = {0, 0, 1};

/*
 * =====================================================================
 * This file contains a series of useful rotation/translation routines
 *
 * Diana Roe 11/91
 * =====================================================================
 */

/*
 * =====================================================================
 * Function: center_of_mass
 * Output: com = array containing x-,y-, and z-coords of center of mass.
 * ---------------------------------------------------------------------
 */
void center_of_mass (float coord[][3], int atomnum, float com[3])
{
  int i, j;

  for(j = 0; j < 3; j++)
    com[j] = 0;
  for (i = 0; i < atomnum; i++)
    for(j = 0; j < 3; j++)
      com[j] += coord[i][j] / (float) atomnum;
  return;
}
	
/*
 * =====================================================================
 * Function:  minimum_angle
 * Purpose:  Given any rotation matrix, it calculates the minimum rotation
 *		angle between the 2 objects described by the rotation matrix.
 *		It uses the formula 2cos(theta)=trace-1 (where trace is the
 *		sum of the diangonals.
 * Output: returns the angle in radians.
 * ---------------------------------------------------------------------
 */

float minimum_angle (float matrix[3][3])
{
	float angle,t;
	t=(trace(matrix)-1.0)/2.0;		
	angle=(float) acos((float)(t));		
	return angle;
}
/*
 * =====================================================================
 * Function: make_rot_matrix
 * Purpose: this function takes an angle and an axis and writes out the
 *		corresponding rotation matrix.
 * Input: angle==in radians. axis-- 0:  x-axis, 1: y-axis, 2: z-axis. 
 *		(Use constants X_AXIS,Y_AXIS,Z_AXIS)
 * Output: rotation matrix in [row][column]
 * ---------------------------------------------------------------------
 */
int make_rot_matrix( float angle, int axis, float matrix[3][3] )
{
	int i,j;
	float cost,sint;
	cost=(float)cos((float)angle);
	sint=(float)sin((float)angle);
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			matrix[i][j]=0;
	switch (axis){
		case X_AXIS:    
			matrix[1][1]=cost;
			matrix[2][2]=cost;
			matrix[1][2]=-sint;
			matrix[2][1]=sint;
			break;
		case Y_AXIS:    
			matrix[0][0]=cost;
			matrix[2][2]=cost;
			matrix[0][2]=sint;
			matrix[2][0]=-sint;
			break;
		case Z_AXIS:    
			matrix[0][0]=cost;
			matrix[1][1]=cost;
			matrix[0][1]=-sint;
			matrix[1][0]=sint;
			break;
		default:
			return FALSE;
	}
	matrix[axis][axis]=1;
	return TRUE;
}
/*
 * =====================================================================
 * Function: Allign_to_z_axis
 * Purpose: takes coords of 1 atoms and pt of origin and finds the rotation
 *		matrix around origin to place first atom on the z-axis 
 *Input: coord1== coordinates of 2 atoms.  Origin== origin through
 *	    which to rotate.
 *Output: matrix =rotation matrix for rotation about origin. 
 *		trans= translation --to go from origin as origin to 0 as origin.
 *		NOTE: Can either subtract origin from coordinates and then
 * 		multiply by matrix to get coordinates with origin as center of
 *		origin.  Or multiply unsubtracted coordinates by matrix and add
 *		trans to get coordinates with 0 as center of origin.
 * ---------------------------------------------------------------------
 */
void allign_to_z_axis(float coord1[3],float origin[3], float matrix[3][3],
		float trans[3])
{
/*	static float x[3]={1,0,0},y[3]={0,1,0},z[3]={0,0,1};*/
	float v1[3],vxy[3],v1rot[3];
	float t1matrix[3][3],t2matrix[3][3];
	int i;
	for (i=0;i<3;i++){
		v1[i]=coord1[i]-origin[i];
	}
/* 
 * 	First place atom1 on z-axis.  This is done by taking the projection
 *   of atom1 on the xy plane and rotating it about z onto the y- axis, 
 *	which will rotate atom1 onto the y-z plane,followed by a rotation about 
 *	the x-axis onto the z-axis.  Result in matrix.
 */ 
	vxy[0]=v1[0];vxy[1]=v1[1];vxy[2]=0;
	vector_rotate(vxy, y_axis, Z_AXIS, t1matrix);
	mvmult3(t1matrix,v1,v1rot);
	vector_rotate(v1rot,z_axis,X_AXIS,t2matrix);
	mmult3(t2matrix,t1matrix,matrix);
	mvmult3(matrix,v1,v1rot);
/*
 *  Finally calculate translation 
 */	
	mvmult3(matrix,origin,trans);
	for (i=0;i<3;i++)
		trans[i]=origin[i]-trans[i];
	return;
}
/*
 * =====================================================================
 * Function: Allign_to_xz_axes
 * Purpose: takes coords of 2 atoms and pt of origin and finds the rotation
 *		matrix around origin to place first atom on the z-axis and align 
 *		second atom to the x-z plane.
 *Input: coord1,coord2 == coordinates of 2 atoms.  Origin== origin through
 *	    which to rotate.
 *Output: matrix =rotation matrix for rotation about origin. 
 *		trans= translation --to go from origin as origin to 0 as origin.
 *		NOTE: Can either subtract origin from coordinates and then
 * 		multiply by matrix to get coordinates with origin as center of
 *		origin.  Or multiply unsubtracted coordinates by matrix and add
 *		trans to get coordinates with 0 as center of origin.
 * ---------------------------------------------------------------------
 */
void allign_to_xz_axes(float coord1[3], float coord2[3], float origin[3],
		float matrix[3][3], float trans[3] )
{
/*	static float x[3]={1,0,0},y[3]={0,1,0},z[3]={0,0,1};*/
	float v1[3],v2[3],v2rot[3],vxy[3],v1rot[3];
	float t1matrix[3][3],t2matrix[3][3],t3matrix[3][3];
	float v1rot2[3],v2rot2[3];
	int i;
	for (i=0;i<3;i++){
		v1[i]=coord1[i]-origin[i];
		v2[i]=coord2[i]-origin[i];
	}
/* 
 * 	First place atom1 on z-axis.  This is done by taking the projection
 *   of atom1 on the xy plane and rotating it about z onto the y- axis, 
 *	which will rotate atom1 onto the y-z plane,followed by a rotation about 
 *	the x-axis onto the z-axis.  Result in t3matrix.
 */ 
	vxy[0]=v1[0];vxy[1]=v1[1];vxy[2]=0;
	vector_rotate(vxy,y_axis,Z_AXIS,t1matrix);
	mvmult3(t1matrix,v1,v1rot);
	vector_rotate(v1rot,z_axis,X_AXIS,t2matrix);
	mmult3(t2matrix,t1matrix,t3matrix);
	mvmult3(t3matrix,v1,v1rot);
	mvmult3(t3matrix,v2,v2rot);
/*
 *   Now  find matrix to rotate rotated v2 onto xz plane while maintaining
 *	v1 on z.  This is done by a rotation of projection of v2 on xy plane
 *	onto the x-axis about the z-axis .
 */
	vxy[0]=v2rot[0];vxy[1]=v2rot[1];
	vector_rotate(vxy,x_axis,Z_AXIS,t1matrix);
/*
 *  Finally multiply matricies to get total rotation/translation matrix
 */	
	mvmult3(t1matrix,v1rot,v1rot2);
	mvmult3(t1matrix,v2rot,v2rot2);
	mmult3(t1matrix,t3matrix,matrix);
	mvmult3(matrix,v1,v1rot2);
	mvmult3(matrix,v2,v2rot2);
	mvmult3(matrix,origin,trans);
	for (i=0;i<3;i++)
		trans[i]=origin[i]-trans[i];
	
	return;
}
/*
 * =====================================================================
 * Function: vector_rotate
 * Purpose: Calculate a rotation matrix to rotate vector 1 to vector 2,
 *		about a specified axis.
 * Input: v1= first vector (rotatable vector). 
 *		v2=fixed vector.  
 *		axis=axis to rotate. 0=x-axis, 1=y-axis, 2=z-axis.
 *		(Use constants X_AXIS,Y_AXIS,Z_AXIS)
 * Output: matrix = 3X3 matrix to rotate v1 onto v2
 * ---------------------------------------------------------------------
 */
void vector_rotate( float v1[3], float v2[3], int axis, float matrix[3][3])
{
	float phi,ref1[3],ref2[3];
	int i;
	if ((v1[0]==0) && (v1[1] == 0) && (v1[2] ==0))
		phi=0;
	else
		phi=vangle(v1,v2);
/*
 *	check sign of phi by comparing to cross of vectors to axis vector.  
 */	
	for (i=0;i<3;i++)
		ref1[i]=0;
	ref1[axis]=1;
	vcrossv3(v1,v2,ref2);
	if (vdotv3(ref1,ref2)>0)
		phi=(-phi);
	make_rot_matrix(phi,axis,matrix);
	return;
}
		

/*
 * =====================================================================
 * Function: RotatePhiAboutAxis
 * Purpose: Calculate a rotation matrix to rotate by phi 
 *		about specified axis.  
 * Input: phi = amount to rotate (in radians) 
 *		origin = origin of rotation for  vector.
 *		axis = endpt of rotation axis (with origin as other end pt) 
 * ---------------------------------------------------------------------
 */
void RotatePhiAboutAxis( float phi, float origin[3], 
	float axis[3], float matrix[3][3],float trans[3])
{

	float t1matrix[3][3], rotmatrix[3][3], t1invmatrix[3][3];
	float t3matrix[3][3];
	float uaxis[3], taxis[3], saxis[3];
	int i;

	for (i = 0; i < 3; i++)
		uaxis[i] = axis[i] - origin[i];
/*
 *	find orthonormal axis system including uaxis 
 */
	normalize3(uaxis);
	get_norm3(uaxis, saxis);
	vcrossv3(saxis, uaxis, taxis);
	for (i = 0; i < 3; i++){
		t1matrix[0][i] = saxis[i];
		t1matrix[1][i] = taxis[i];
		t1matrix[2][i] = uaxis[i];
	}
	transpose3(t1matrix, t1invmatrix);
/*
 *	rotate by angle, then tranform back
 */

	make_rot_matrix(phi, Z_AXIS, rotmatrix);
	mmult3(rotmatrix, t1matrix, t3matrix);
	mmult3(t1invmatrix, t3matrix, matrix);
	mvmult3(matrix,origin,trans);
     for ( i = 0; i<3; i++ )
          trans[i] = origin[i] - trans[i];

	return;
}
/*
 * ==========================================================================
 * Function: RotateV2toAngleFromV1 
 * Purpose:  Rotates V2 to specified angle phi from v1 (w/random angle psi)
 * Input: coord1 = fixed point
 *		coord2 = rotatable point
 *		origin = origin  for 2 points.
 *		target_angle = target angle about axis  for vectors.
 * Output:
 *		matrix, trans
 * --------------------------------------------------------------------------
 */ 
void RotateV2toAngleFromV1( float coord1[3], float coord2[3], float origin[3], 
		float target_angle, float matrix[3][3], 
		float trans[3] )
{
	float axis[3];
	float v1[3], v2[3], v2rot[3], v1rot[3];
	float matrix1[3][3], trans1[3];
	float matrix2[3][3];
	float invmatrix[3][3], tmatrix[3][3];
	float act_angle, phi;
	int i;
/*
 *	use cross-product of points as rotation axis
 */
	for (i = 0; i < 3; i ++){
		v1[i] = coord1[i] - origin[i];
		v2[i] = coord2[i] - origin[i];
	}
	vcrossv3(v1, v2, axis);
	translate3(axis, origin, 1);
/*
 *	rotate axis onto z-axis
 */
	allign_to_z_axis( axis, origin, matrix1, trans1);
	mvmult3( matrix1, v1, v1rot);
	mvmult3( matrix1, v2, v2rot);
/*
 *	calculate rotation angle phi 
 *	(note: no sign worries since axis is cross product of vectors)
 */
	act_angle = vangle(v1rot, v2rot);
	phi = act_angle - target_angle;
	make_rot_matrix( phi, Z_AXIS, matrix2);

/*
 *	rotate by angle, then rotate back.	
 */
	transpose3(matrix1, invmatrix);
	mmult3(matrix2, matrix1, tmatrix);
	mmult3(invmatrix, tmatrix, matrix);
	mvmult3(matrix, origin, trans);
	for (i = 0; i < 3; i++)
		trans[i] = origin[i] - trans[i];
	
	return;
}

