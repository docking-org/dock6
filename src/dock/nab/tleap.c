#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "nab.h"

/*
 * Call tleap and execute the commands in the string leap_cmds.
 * The path of teLeap is $DOCK_HOME/bin.
 * The leap files used are in $DOCK_HOME/parameters/.
 */

int tleap( char* leap_cmds )
{

/* the leap command string, lcmd, can get really big if DOCK_HOME is big */
#define MAXSTRING 3067

    char lcmd[ MAXSTRING ];
    char *nabhome;
    FILE *fp;
    char lcfile[] = "tleap.in";
    const char leap_invocation[] = 
"%s/bin/teLeap -s -f %s -I%s/parameters/leap/cmd -I%s/parameters/leap/parm -I%s/parameters/leap/prep -I%s/parameters/leap/lib > tleap.out";

    nabhome = getenv( "DOCK_HOME" ); 
    if( nabhome == NULL ){
        fprintf( stderr, "Warning: DOCK_HOME is not defined.\n" );
        fprintf( stderr, "   Using the relative path: ../../..\n" );
    }

    fp = fopen( lcfile, "w" );

    fprintf( fp, "%s\n",  leap_cmds );
    fprintf( fp, "quit\n" );
    fclose( fp );

    sprintf( lcmd, leap_invocation,
             nabhome, lcfile, nabhome, nabhome, nabhome, nabhome );

    fprintf( nabout, "\nRunning: %s\n", lcmd );
    system( lcmd );

/* do not delete the leap input file if debugging */
#ifdef NDEBUG
    unlink( lcfile );
#else
#endif

    return( 0 );
}

