/*
Originally transcribed by Dan Gschwend.
Converted to C, and generalized for any optimization problem by Todd Ewing.
3/96

Multidimensional minimization of the function funk(x) where x is an
n-dimensional vector, by the downhill simplex method of Nelder 
and Mead.  In this case, funk is the score function.
Input is a matrix p whose n+1 rows are n-dimensional
vectors which are the vertices of the starting simplex.  Also input 
is the vector y of length n+1, whose  components must be 
pre-initialized to the values of funk evaluated at the n+1 
vertices (rows) of p.
On output, p and y will have been reset to n+1 new points all 
within score_converge of a minimum function value; iteration gives the number
of iterations taken.

Taken from Numerical Recipes:  The Art of Scientific
Computing, by Press, Flannery, Teukolsky, and Vetterling.
1986 by Cambridge University Press, C edition
pp. 292-293.

References and further reading:
  Nelder, J.A., and Mead, R. 1965, Computer Journal, vol. 7,
      p. 308.
  Yarbro, L.A., and Deming, S.N. 1974, Analytica Chim. Acta,
      vol. 73, p. 391.
  Jacoby, S.L.S., Kowalik, J.S., and Pizzo, J.T. 1972, 
      Iterative Methods for Nonlinear Optimization Problems
      (Englewood Cliffs, N.J.:  Prentice-Hall).
*/

#include "define.h"
#include "utility.h"
#include "global.h"

float simplex_optimize
(
  void *simplex,		/* Structure of info passed to scoring fn */
  float *solution,		/* Array of variables to optimize */
  int size,			/* Number of variables */
  float score_converge,		/* Convergence criteria */
  int *iteration,		/* Number of simplex iterations */
  int iteration_max,		/* Maximum number of simplex iterations */
  float (*simplex_score)	/* Scoring function */
    (void *, float *),
  float *delta			/* Minimizer score improvement (init - fin) */
)
{
  int i, j;
  int ilo, ihi, inhi;
  float **p = NULL;
  float *pr = NULL;
  float *prr = NULL;
  float *pbar = NULL;
  float *y = NULL;
  float ypr = 0;
  float yprr = 0;
  float alpha = 1.0;	/* range: 0=no extrap, 1=unit step extrap, higher OK */
  float beta = 0.5;	/* range: 0=no contraction, 1=full contraction */
  float optimum;
  int replace_flag;	/* flag for whether bad vertex replaced */

/*
* Allocate space for arrays
*/
  ecalloc
  (
    (void **) &p,
    (size + 1),
    sizeof (float *),
    "simplex array",
    global.outfile
  );

  for (i = 0; i < size + 1; i++)
    ecalloc
    (
      (void **) &p[i],
      size,
      sizeof (float),
      "simplex array",
      global.outfile
    );

  ecalloc
  (
    (void **) &y,
    size + 1,
    sizeof (float),
    "simplex array",
    global.outfile
  );

  ecalloc
  (
    (void **) &pr,
    size,
    sizeof (float),
    "simplex array",
    global.outfile
  );

  ecalloc
  (
    (void **) &prr,
    size,
    sizeof (float),
    "simplex array",
    global.outfile
  );

  ecalloc
  (
    (void **) &pbar,
    size,
    sizeof (float),
    "simplex array",
    global.outfile
  );

/*
* Begin the simplex optimization
*/
  *iteration = 0;

  do
  {
    if (*iteration == 0)
    {
/*
*     Initialize simplex array.
*       The first row contains values initialized externally.
*       All other rows contain perturbations from the first row,
*         by amounts between -1.0 and 1.0.
*     1/97 te
*/
      for (j = 0; j < size; j++)
        p[0][j] = solution[j];

      for (i = 1; i < size + 1; i++)
        for (j = 0; j < size; j++)
          p[i][j] =
            solution[j] +
            2.0 * ((float) rand () / (float) RAND_MAX - 0.5);

      for (i = 0; i < size + 1; i++)
      {
        for (j = 0; j < size; j++)
          solution[j] = p[i][j];

        y[i] = simplex_score (simplex, solution);
      }

      *delta = y[0];
    }

    else
    {
/*
*     Begin a new iteration.  Compute the vector average of all points
*     except the highest, i.e. the center of the "face" of the simplex
*     across from the high point.  We will subsequently explore the
*     ray from the high point through that center.
*/
      for (j = 0; j < size; j++)
        pbar[j] = 0.0;

      for (i = 0; i < size + 1; i++)
        if (i != ihi)
          for (j = 0; j < size; j++)
            pbar[j] += p[i][j];

/*
*     Extrapolate by a factor alpha through the face, i.e. reflect the
*     simplex from the high point.
*/
      for (j = 0; j < size; j++)
      {
        pbar[j] /= (float) size;
        solution[j] = pr[j] = (1.0 + alpha) * pbar[j] - alpha * p[ihi][j];
      }

/*
*     Evaluate the function at the reflected point.  
*/
      ypr = simplex_score (simplex, solution);

      if (ypr <= y[ilo])
      {
/*
*       Gives a result better than the best point, so try an 
*       additional extrapolation by a factor alpha.
*/
        for (j = 0; j < size; j++)
          solution[j] = prr[j] = (1.0 + alpha) * pr[j] - alpha * pbar[j];

/*
*       Check out the function there.
*/
        yprr = simplex_score (simplex, solution);

        if (yprr < y[ilo])
        {

/*
*         The additional extrapolation succeeded, and replaces 
*         the high point.
*/
          for (j = 0; j < size; j++)
            p[ihi][j] = prr[j];

          y[ihi] = yprr;
        }

        else
        {

/*
*         The additional extrapolation failed, but we can still use
*         the reflected point
*/
          for (j = 0; j < size; j++)
            p[ihi][j] = pr[j];

          y[ihi] = ypr;
        }
      }

      else if (ypr >= y[inhi])
      {
/*
*       The reflected point is worse than the second-highest.  
*       If it's better than the highest, then replace the highest.
*/
        replace_flag = FALSE;

        if (ypr < y[ihi])
        {
          for (j = 0; j < size; j++)
            p[ihi][j] = pr[j];

          y[ihi] = ypr;
          replace_flag = TRUE;
        }

/*
*       But look for an intermediate lower point, in other words,
*       perform a contraction of the simplex along one dimension, then
*       evaluate the function.  
*/
        for (j = 0; j < size; j++)
          solution[j] = prr[j] = beta * p[ihi][j] + (1.0 - beta) * pbar[j];

        yprr = simplex_score (simplex, solution);

        if (yprr < y[ihi])
        {

/*
*         Contraction gives an improvement, so accept it.
*/
          for (j = 0; j < size; j++)
            p[ihi][j] = prr[j];

          y[ihi] = yprr;
          replace_flag = TRUE;
        }

        if (replace_flag == FALSE)
        {
/*
*         Can't seem to get rid of that high point.  Better contract
*         around the lowest (best) point.
*/
          for (i = 0; i < size + 1; i++)
            if (i != ilo)
          {
            for (j = 0; j < size; j++)
              solution[j] = p[i][j] = pr[j] = 0.5 * (p[i][j] + p[ilo][j]);

            y[i] = simplex_score (simplex, solution);
          }
        }
      }

      else
      {
/*
*       We arrive here if the original reflection gives a middling point.
*       Replace the old high point and continue.
*/
        for (j = 0; j < size; j++)
          p[ihi][j] = pr[j];

        y[ihi] = ypr;
      }
    }
/*
*   Identify the best and worst vertices in the current simplex.
*   We must determine which point is the highest (worst), next-
*   highest, and lowest (best).
*/
    ilo = 0;

    if (y[0] > y[1])
    {
      ihi = 0;
      inhi = 1;
    }

    else
    {
      ihi = 1;
      inhi = 0;
    }

/*
*   Loop over points in the simplex
*/
/*
    if (size == 6)
    fprintf (global.outfile, "  simplex scores:");
*/

    for (i = 0; i < size + 1; i++)
    {
/*
    if (size == 6)
      fprintf (global.outfile, " %g", y[i]);
*/

      if (y[i] < y[ilo])
        ilo = i;

      if (y[i] > y[ihi])
      {
        inhi = ihi;
        ihi = i;
      }

      else if (y[i] > y[inhi])
      {
        if (i != ihi)
          inhi = i;
      }
    }
/*
    if (size == 6)
    fprintf (global.outfile, "\n");
*/
  } while
      ((++(*iteration) < iteration_max) &&
      (ABS (y[ihi] - y[ilo]) > score_converge));

/*
* Store the best scoring vertex of the simplex
*/ 
  for (i = 0; i < size; i++)
    solution[i] = p[ilo][i];

  optimum = y[ilo];

  *delta -= optimum;

/*
* Free up space allocated for simplex arrays
*/
  for (i = 0; i < size + 1; i++)
    efree ((void **) &p[i]);

  efree ((void **) &p);
  efree ((void **) &y);
  efree ((void **) &pr);
  efree ((void **) &prr);
  efree ((void **) &pbar);

/*
* Return the best scoring vertex
  fprintf (global.outfile, "Finished simplex, best score %g\n", optimum);
*/
  return optimum;
}

