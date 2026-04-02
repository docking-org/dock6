/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*

Written by Todd Ewing
5/97

*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>

/* Preprocessor constants, strings, and macros */

#define DOCK_VERSION "4.0.1"
#define RELEASING_DATE "May 17, 1998"
#define PARAMETER_PATH "PARAMETER_PATH_NOT_SET/"
#define TRUE 1
#define FALSE 0
#define NEITHER -1
#define INITIAL_SCORE 1000
#define DISTANCE_MIN 1.0
#define PI 3.141592654
#define PRINTVAR(var, fmt) printf (#var " = " #fmt "\n", var)
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define NINT(x) (int) ((x) > 0 ? ((x) + 0.5) : ((x) - 0.5))
#define SIGN(x) ((x) < 0 ? -1 : 1)
#define NON_ZERO(x) (ABS (x) < 0.00001 ? FALSE : TRUE)
#define POWER(b, e, p) \
{ \
  int word; \
  float run; \
  POWER_CORE(b, e, p) \
}

#define POWER_CORE(base, exponent, product) \
  word = exponent; product = 1.0; run = base; \
  while (word) \
  { \
    if (word & 1) \
      product *= run; \
    run *= run; \
    word >>= 1; \
  }

/* String variable types */

typedef char STRING5[6];
typedef char STRING10[11];
typedef char STRING20[21];
typedef char STRING40[41];
typedef char STRING50[51];
typedef char STRING60[61];
typedef char STRING80[81];
typedef char STRING100[101];
typedef char STRING200[201];
typedef STRING80 FILE_NAME;
typedef float X;
typedef float XYZ[3];
typedef unsigned long MASK;
#define MASK_LENGTH 32

typedef struct sort_int
{
  int id;
  int value;
} SORT_INT;

int compare_int_descend (SORT_INT *a, SORT_INT *b);
int compare_int_ascend (SORT_INT *a, SORT_INT *b);

typedef struct singly_linked_integer
{
  int value;
  struct singly_linked_integer *next;
} SLINT;

typedef struct singly_linked_integer_pair
{
  int i, j;
  struct singly_linked_integer_pair *next;
} SLINT2;

typedef struct singly_linked_float
{
  int value;
  struct singly_linked_float *next;
} SLFLT;

