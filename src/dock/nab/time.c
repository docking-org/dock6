#include <time.h>
#include	"defreal.h"

/* NAB interface to date and time. Note, repeated calls will leak memory! */

char	*date( void )
{

	static char string[11];
	size_t smax=11;
	time_t now;

	now = time( NULL );
	strftime( string, smax, "%m/%d/%Y", localtime( &now ) );

	return( string );

}

char	*timeofday( void )
{

	static char string[9];
	size_t smax=9;
	time_t now;

	now = time( NULL );
	strftime( string, smax, "%H:%M:%S", localtime( &now ) );

	return( string );

}

char	*ftime( char *fmt )
{

/*   NAB interface to system routine strftime   */

	static char string[50];
	size_t smax=50;
	time_t now;

	now = time( NULL );
	strftime( string, smax, fmt, localtime( &now ) );
	string[49] = '\0';  /* in case of overflow, no explicit checks here */

	return( string );

}

REAL_T	second( void )
{

#ifdef DIFFTIME
	/* 
	 *  Here use the standard C difftime() function to get calendar times.
	 *  The disadvantage here is that you just get seconds, not fractions
	 *  of a second, but is easier to interpret in parallel runs.
	 *  Based on code in Harbisson & Steele, 5th ed., p. 448
	 */

	struct tm ref_struct = {0};
	static time_t tm_ref;
	static int first=1;

	if( first ){                  /* get time for April 15, 2000  */
		ref_struct.tm_year = 100;
		ref_struct.tm_mon  = 3;
		ref_struct.tm_mday  = 15;
		tm_ref = mktime( &ref_struct);
		first = 0;
	}
	
	return difftime( time(NULL), tm_ref );

#else
	/*
	 *  This is an interface to clock(), which returns
	 *  processor time, not calendar (wallclock) time
	 */
	REAL_T rv;
#  if defined CLOCKS_PER_SEC
	rv = clock();
	rv /= CLOCKS_PER_SEC;
#  else
	rv = 0.0;
#  endif
	return( rv );

#endif

}
