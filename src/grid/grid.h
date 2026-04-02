/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
9/96
*/

typedef struct grid_struct
{
  int output_molecule;		/* Flag for whether to write receptor */

  int visualize;		/* Write out grids in MIDAS format */

  XYZ box_com;			/* Center of mass of box enclosing grids */
  XYZ box_dimension;		/* Dimensions of box enclosing grids */

  FILE_NAME in_file_name;	/* Input receptor file name */
  FILE_NAME out_file_name;	/* Output receptor file name */

  FILE *in_file;		/* Input receptor file pointer */
  FILE *out_file;		/* Output receptor file pointer */

  FILE_NAME box_file_name;	/* Grid-enclosing box file */

} GRID;

