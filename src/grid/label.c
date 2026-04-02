/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

#include "define.h"
#include "mol.h"
#include "label.h"

void free_labels (LABEL *label)
{
  free_vdw_labels (&label->vdw);
  free_chemical_labels (&label->chemical);
  free_flex_labels (&label->flex);
}

