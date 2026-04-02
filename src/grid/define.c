/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*
Written by Todd Ewing
2/97
*/
#include "define.h"

int compare_int_descend (SORT_INT *a, SORT_INT *b)
{
  return a->value == b->value ? 0 : (a->value < b->value ? 1 : -1);
}

int compare_int_ascend (SORT_INT *a, SORT_INT *b)
{
  return a->value == b->value ? 0 : (a->value > b->value ? 1 : -1);
}
