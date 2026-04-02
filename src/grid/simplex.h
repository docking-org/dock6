/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

/*

General-purpose simplex optimizer

Define simplex_struct and simplex_score() in file containing calling routine.
11/96 te

*/

float simplex_optimize
(
  void	*simplex,		/* Structure of info passed to scoring fn */
  float	*solution,		/* Array of variables to optimize */
  int	size,			/* Number of variables */
  float	score_converge,		/* Convergence criteria */
  int	*iteration,		/* Number of simplex iterations */
  int	iteration_max,		/* Maximum number of simplex iterations */
  float	(*simplex_score)	/* Scoring function */
    (void *, float *),
  float *delta			/* Minimizer score improvement (init - fin) */
);
