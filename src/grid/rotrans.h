/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
 * header file  for rotrans.c
 */
#ifndef ROTRANS_INCLUDE 
#define ROTRANS_INCLUDE 
#define X_AXIS 0 
#define Y_AXIS 1
#define Z_AXIS 2
#define NUMSYMOLS 20

void center_of_mass(float coord[][3], int number, float com[3]);
float minimum_angle( float matrix[3][3] );
int make_rot_matrix( float angle, int axis, float matrix[3][3] );
void allign_to_z_axis( float coord1[3],float origin[3], 
	float matrix[3][3], float trans[3] );
void allign_to_yz_axes( float coord1[3], float coord2[3], 
	float origin[3], float matrix[3][3], float trans[3] );
void vector_rotate( float v1[3], float v2[3], int axis, 
		float matrix[3][3] );
void RotatePhiAboutAxis( float phi, float origin[3],
     float vaxis[3], float matrix[3][3],float trans[3]);
void RotateV2toAngleFromV1( float v1[3], float v2[3], float
		origin[3],  float target_angle, 
		float matrix[3][3], float trans[3] );
#endif
