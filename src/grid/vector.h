/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/
/*
 * header file for vector.c
 * written by Diana Roe 
 */
#ifndef VECTORS_INCLUDE
#define VECTORS_INCLUDE
typedef float VECTOR3[3];
#include <stdio.h>
#define RADTODEG	57.2957795
#define DEGTORAD    0.0174533

void copy_vector3 (float v2[3], float v1[3]); 		/* v2 <- v1*/

void copy_matrix3 (float m2[3][3], float m1[3][3]);		/* m2 <- m1*/

float dist3 (float a[3], float b[3]);

float square_distance (float a[3], float b[3]);

void translate3 (float a[3], float b[3], int sign); 

void rotate3 (float coords[3], float matrix[3][3], float trans[3]);

float vangle(float a[3], float b[3]);		   

float angle(float x[3], float y[3], float z[3]) ;

float dihed (float a[3], float b[3], float c[3], float d[3]);

float vdotv3(float x[3], float y[3]);

void   vcrossv3(float a[3], float b[3], float c[3]); 

float   normalize3(float v[3]);		 /*normalizes a vector*/

float norm3(float v[3]);			 /*    calculates norm*/

void    get_norm3(float x[3], 		 /* finds y normal to x */
				    float y[3]); 	

void transpose3(float matrix[3][3], float tmatrix[3][3]);

void mvmult3(float matrix[3][3], float v1[3], float v2[3]); /* v2=matrix*v1 */

void vmmult3(float v1[3], float matrix[3][3], float v2[3]); /* m2=v1*matrix */

void mvmult(float matrix[4][4], float v1[3], float v2[3]);
	/* effective multplies v1 by rot/trans matrix => v2*/

void mmult3(float matrix1[3][3], float matrix2[3][3], float outmatrix[3][3]);
	/*(matrix1)(matrix2) = outmatrix*/

void mvinmult(float matrix[4][4], float v1[3], float v2[3]);
	/*v2 = v1*(matrix)^Transpose*/

float trace(float matrix[3][3]);

int read_matrix(FILE *infile, float matrix[4][4]);

#endif
