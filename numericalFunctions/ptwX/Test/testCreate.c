/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>

#include <nf_utilities.h>
#include <ptwX.h>

#define nXs 23

static int verbose = 0;

void printMsg( const char *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0;
    int64_t i, errors = 0;
    double xs[nXs], d;
    nfu_status status;
    ptwXPoints *ptwX;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfu_printErrorMsg( "Error %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    for( i = 0; i < nXs; i++ ) xs[i] = -10. + i * i;
    
    ptwX = ptwX_create( 10, nXs, xs, &status );

    for( i = 0; i < nXs; i++ ) {
        d = *ptwX_getPointAtIndex( ptwX, i );
        if( d != xs[i] ) {
            errors++;
            nfu_printMsg( "Error %s: bad ptwX_getPointAtIndex return value = %g, xs[%d] = %g", __FILE__, d, i, xs[i] );
        }
    }

    exit( errors );
}
