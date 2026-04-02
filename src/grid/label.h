/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
10/95
*/
#include "label_node.h"
#include "label_vdw.h"
#include "label_chem.h"
#include "label_flex.h"

/* Structures to store atom labelling definitions */

typedef struct label_struct
{
  LABEL_VDW vdw;
  LABEL_CHEMICAL chemical;
  LABEL_FLEX flex;

  FILE_NAME definition_path;

} LABEL;


void free_labels (LABEL *);
