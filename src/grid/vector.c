/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

#include "define.h"
#include "vector.h"
/*
 * ====================================================================
 *  All routines in this file assume fixed size arrays.  See vectors2
 *  for routines with variable indexed arrays.  The notation for
 *  matricies is [row][column]
 * ====================================================================
 */
/*
 * copies v2 from v1
 */
void copy_vector3(float v2[3], float v1[3] )
{
	register int i;
	for (i =0; i < 3; i++)
		v2[i] = v1[i];
	return;
}
/*
 * copies m2 from m1
 */
void copy_matrix3(float m2[3][3], float m1[3][3] ) 
{
	register int i,j;
	for (i =0; i < 3; i++)
		for (j =0; j < 3; j++)
		m2[i][j] = m1[i][j];
	return;
}
float dist3(float a[3], float b[3])
{
	return sqrt(SQR(a[0] - b[0]) + SQR(a[1] - b[1]) + SQR(a[2] - b[2]));
}

float square_distance (float a[3], float b[3])
{
	return
		(a[0] - b[0]) * (a[0] - b[0]) +
		(a[1] - b[1]) * (a[1] - b[1]) +
		(a[2] - b[2]) * (a[2] - b[2]);
}
void vector_add(float a[3], float b[3], float c[3], int sign)
{
	int i;
	for (i = 0; i < 3; i++)
		c[i] = a[i] + ((float) sign) * b[i];
}

void translate3(float a[3], float b[3], int sign)
{
	int i;
	for (i = 0; i < 3; i++)
		a[i] += ((float) sign) * b[i];
}
/*
*/
/*
 * =====================================================================
 * Function: rotate3
 * Purpose: rotate given coordinates by a given 3X3 matrix and trans
 *        vector
 * ---------------------------------------------------------------------
 */
void rotate3(float coords[3], float matrix[3][3], float trans[3])
{
	float v1[3];
	copy_vector3(v1, coords);
	mvmult3(matrix, v1, coords);
	translate3(coords, trans, 1);
	return;
}


/* 
 * Function vangle(a,b)
 * Purpose: calculate the angle between two vectors in radians.
 */
float vangle (float a[3], float b[3])
{
  float theta, costheta;

  if (!NON_ZERO (a[0] - b[0]) &&
    !NON_ZERO (a[1] - b[1]) &&
    !NON_ZERO (a[2] - b[2]))
    return 0;

  costheta = norm3 (a) * norm3 (b);

  if (!NON_ZERO (costheta))
    return 0;

  costheta = vdotv3 (a, b) / costheta;

  if (costheta > 1.0)
    costheta = 1.0;

  if (costheta < -1.0)
    costheta = -1.0;

  theta = acos (costheta);
  return theta;
}

/*
 * Function: axis_angle
 * Purpose: calculate the angle and axis of rotation between two vectors
 */
float axis_angle (float axis[3], float a[3], float b[3])
{
  float theta;

  theta = vangle (a, b);

  if (!NON_ZERO (theta))
    return 0.0;

  vcrossv3 (a, b, axis);

  return theta;
}

/* 
 * Function uangle(a,b)
 * Purpose: calculate the angle between 2 unit vectors a,b in radians
 * 11/95 te
 */
float uangle (float a[3], float b[3])
{
  float theta, costheta;

  if (((a[0]==0) && (a[1] == 0) && (a[2] == 0)) ||
    ((b[0]==0) && (b[1] == 0) && (b[2] == 0)))
  return 0;

  costheta = vdotv3 (a, b);
  theta = acos (costheta);
  return theta;
}

float angle(float x[3],float y[3],float z[3]) 
/* compute the angle between points x, y and z. */
{
        register int    i;
        register float acc;
        register float d1, d2;

        acc = 0;
        for (i = 0; i < 3; i++)
                acc += (y[i] - x[i]) * (y[i] - z[i]);
        d1 = dist3(x, y);
        d2 = dist3(z, y);
        if (d1 <= 0 || d2 <= 0)
                return 0.0;
        acc /= (d1 * d2);
        if (acc > 1)
                acc = 1;
        else if (acc < -1)
                acc = -1;
        return acos(acc);
}

float dihed (float a[3],float b[3],float c[3],float d[3])
{
        float q[3],r[3],s[3],t[3],u[3];
        float v[3];
        static float z[3] = {0,0,0};
        register int i;
        register float acc;
        register float w;

        for (i=0;i<3;i++) {
                q[i] = b[i] - a[i];
                r[i] = b[i] - c[i];
                s[i] = c[i] - d[i];
        }
        vcrossv3 (q,r,t);
        vcrossv3 (s,r,u);
        vcrossv3 (u,t,v);
        w = vdotv3(v,r);
        acc = angle (t,z,u);
        if (w<0)  acc = -acc;
        return (acc);
}

float vdotv3 (float x[3], float y[3])
{
        register float acc;
        register int i;

        acc = 0;
        for (i=0; i<3; i++)
                acc += x[i] * y[i];
        return(acc);
}

/*
 * cross-product
 *  (c <- a x b)
 */
void vcrossv3 (float a[3], float b[3], float c[3])
{
        c[0] = a[1] * b[2] - a[2] * b[1];
        c[1] = a[2] * b[0] - a[0] * b[2];
        c[2] = a[0] * b[1] - a[1] * b[0];
}

/*
 * normalize:
 *   Normalize a vector
 */
float normalize3(float v[3])
{
     float     d;

     d = norm3(v);
     v[0] /= d;
     v[1] /= d;
     v[2] /= d;
     return d;
}
/*
 * norm3:
 *   Find the norm of a  vector
 */
float norm3(float v[3])
{
     float     d;

     d = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	return d;
}
/*
 * get_norm3
 *	Find a normal (y) to a given vector (x)
 */
void get_norm3(float x[3], float y[3])
{
	if ( x[0] == 0.0 && x[1] == 0.0 ){
		y[0] = 1.0;
		y[1] = 0.0;
		y[2] = 0.0;
	}
	else if ( x[0] == 0.0 && x[2] == 0.0 ){
		y[0] = 0.0;
		y[1] = 0.0;
		y[2] = 1.0;
	}
	else if ( x[1] == 0.0 && x[2] == 0.0 ){
		y[0] = 0.0;
		y[1] = 1.0;
		y[2] = 0.0;
	}
	else{
		y[0] = 0.0;
		y[1] = -x[2];
		y[2] = x[1];
		normalize3(y);
	}
}
/*
 * transpose3:
 * take the transpose of a 3X3 matrix
 */
void transpose3(float matrix[3][3],float tmatrix[3][3])
{
	int i,j;
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
			tmatrix[i][j]=matrix[j][i];
	return;
}
/*
 * Function: mvmult3
 * Purpose:
 *   This function multiplies vector v1 by rotation 
 *   matrix (dimension 3X3) and stores result in vector v2.
 *	v2=matrix*v1
 */
void mvmult3(float matrix[3][3],float v1[3],float v2[3])
{
        int i,j;
        v2[0]=0;v2[1]=0;v2[2]=0;
        for (i=0;i<3;i++) {
                for (j=0;j<3;j++)
                        v2[i]=v2[i]+v1[j]*matrix[i][j];
        }
        return;
}
/*
 * Function: mmult3
 * Purpose:
 *   This function multiplies (matrix1)(matrix2)  
 *	and stores the result in outmatrix.
 */
void mmult3(float matrix1[3][3],float matrix2[3][3],float outmatrix[3][3])
{
  int i,j,k;
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      outmatrix[i][j]=0;

  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      for (k=0;k<3;k++)
        outmatrix[i][j]+=matrix1[i][k]*matrix2[k][j];

  return;
}
/*
 * Function: vmmult3
 * Purpose:
 *   This function multiplies vector v1 by rotation 
 *   matrix (dimension 3X3) and stores result in vector v2 .
 *	m2=v1*matrix
 */
void vmmult3(float v1[3],float matrix[3][3],float v2[3])
{
        int i;
   /*     for (i=0;i<3;i++) 
			v2[i]=0;
        for (i=0;i<3;i++) 
                for (j=0;j<3;j++)
                        v2[i]+=v1[j]*matrix[j][i];*/
	for (i = 0; i<3; i++)
		v2[i] = v1[0]*matrix[0][i] + v1[1]*matrix[1][i] +
			v1[2]*matrix[2][i];
        return;
}
/*
 * Function: mvmult
 * Purpose:
 *   This function effectively multiplies vector v1 by rotation translation
 *   matrix and stores result in vector v2.
 */
void mvmult(float matrix[4][4],float v1[3],float v2[3])
{
        int i,j;
        v2[0]=0;v2[1]=0;v2[2]=0;
        for (i=0;i<3;i++) {
                for (j=0;j<3;j++)
                        v2[i]=v2[i]+v1[j]*matrix[i][j];
                v2[i]=v2[i]+matrix[i][3];
        }
        return;
}


/*
 * Function: mvinmult
 * Purpose:
 *   This function effectively multiplies vector v1 by the transpose of the
 *   rotation translation matrix and stores result in vector v2.
 */
void mvinmult(float matrix[4][4],float v1[3],float v2[3])
{
     int i,j;
     v2[0]=0;v2[1]=0;v2[2]=0;
        for (i=0;i<3;i++) {
                v1[i]=v1[i]-matrix[i][3];
                for (j=0;j<3;j++)
                        v2[j]=v2[j]+v1[i]*matrix[i][j];
        }
     return;
}

/*
 * Function:trace 
 * Purpose: Calculates the trace (the sum of the diagonals) of a 
 *		3 dimensional matrix 
 *
 */
float trace(float matrix[3][3])
{
	int i;
	float t=0;
    	for (i=0;i<3;i++)
		t+=matrix[i][i];
	return t;	
}
/*
 * Function: read_matrix
 * Purpose: read  a 4X4 dimensional matrix from the specified infile
 *		into matrix.
 * Returns: False if problem reading in matrix
 */
int read_matrix(FILE *infile,float matrix[4][4])
{
	int i;
	for (i=0;i<4;i++)
		if (fscanf(infile,"%f %f %f %f",&matrix[i][0],&matrix[i][1],
			&matrix[i][2], &matrix[i][3]) != 4)  
			return FALSE;
	return TRUE;
}

