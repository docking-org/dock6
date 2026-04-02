#include <stdio.h>
#include <math.h>
#include	"defreal.h"

#define	IM1	2147483563
#define	IM2	2147483399
#define	AM	( 1.0 / IM1 )
#define	IMM1	( IM1 - 1 )
#define	IA1	40014
#define	IA2	40692
#define	IQ1	53668
#define	IQ2	52774
#define	IR1	12211
#define	IR2	3791
#define	NTAB	32
#define	NDIV	( 1 + IMM1 / NTAB )
#define	EPS	1.2e-13
#define	RNMX	( 1.0 - EPS )

static	int	idum2 = 0;
static  int idum3 = -1;   /* surrogate for idum  */
static	int	iy = 0;
static	int	iv[ NTAB ];

static REAL_T x;

/*
     Get a pseudo-random number in the range 0 -> 1.  If "idum" is
     negative, reset the seeds so that the same sequence can be
     re-created.

     Note that the "driver" routine, rand2(), which is provided at the 
     bottom of this file, ignores its argument if the argument is positive;
     if the argument is negative, the sequence is reset to a new starting
     position, and the first pseudo-random number of that sequence is 
     returned.

*/

static
REAL_T	rand2a( int *idum )
{
	int		j, k;
	REAL_T	temp;

	if( *idum <= 0 ){
		if( -*idum < 1 )
			*idum = 1;
		else
			*idum = -*idum;
		idum2 = *idum;
		for( j = NTAB + 7; j >= 0; j-- ){
			k = *idum / IQ1;
			*idum = IA1 * ( *idum - k * IQ1 ) - k * IR1;
			if( *idum < 0 )
				*idum += IM1;
			if( j < NTAB )
				iv[ j ] = *idum;
		}
		iy = iv[ 0 ];
	}
	k = *idum / IQ1;
	*idum = IA1 * ( *idum - k * IQ1 ) - k * IR1;
	if( *idum < 0 )
		*idum += IM1;
	k = idum2 / IQ2;
	idum2 = IA2 * ( idum2 - k * IQ2 ) - k * IR2;
	if( idum2 < 0 )
		idum2 += IM2; 
	j = iy / NDIV;
	iy = iv[ j ] - idum2;
	iv[ j ] = *idum;
	if( iy < 1 )
		iy += IMM1;
	if( ( temp = AM * iy ) > RNMX )
		return( RNMX );
	else
		return( temp );
}

/*
   Generate a pseudo-random sequence of numbers with a given mean
   and standard deviation.  Use the Box & Mueller method, but only
   use the first value, since the two values are correlated.
*/

static
REAL_T gaussa( REAL_T *mean, REAL_T *sd, int *idum )
{
	REAL_T fac,gdev1,rsq,s1,s2;

		do {
			s1 = 2.*rand2a(idum) - 1.;
			s2 = 2.*rand2a(idum) - 1.;
			rsq = s1*s1 + s2*s2;
		} while ( rsq >= 1. || rsq == 0.0 );
		fac = sqrt(-2.*log(rsq)/rsq);
		gdev1 = s1*fac;

		return( *sd*gdev1 + *mean );

}

/*   
    Driver routine for randa(), so that state information is kept within
    this routine.  Same for gaussa().
*/

REAL_T  rand2( int *idum )
{
    if( *idum <= 0 ) {
        x = rand2a( idum ); 
        idum3 = *idum;
        return x;
    } else {
        return rand2a( &idum3 );
    }
}

REAL_T gauss( REAL_T *mean, REAL_T *sd, int *idum )
{
    return gaussa( mean, sd, &idum3 );
}
