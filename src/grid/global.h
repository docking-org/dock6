/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Header file for structure containing global variables

Include the following statements before the main() routine:

#include "global.h"
GLOBAL global = {0};

Written by Todd Ewing
10/95
*/

typedef struct global_struct
{
  FILE_NAME executable;		/* Name of executable used for run */
  FILE_NAME job_name;		/* Name identifying particular execution */
  int output_volume;		/* Amount of output information
                                   (NULL = normal, v = verbose, t = terse) */
  FILE *infile;			/* Input file stream */
  FILE *outfile;		/* Output file stream */

} GLOBAL;

extern GLOBAL global;
