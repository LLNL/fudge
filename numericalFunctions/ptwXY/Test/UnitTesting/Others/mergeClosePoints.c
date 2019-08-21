/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include <ptwXY.h>
#include <nf_utilities.h>

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int checkMerge( ptwXYPoints *XYs, double epsilon, int n );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;
    nfu_status status;
    ptwXYPoints *XYs;
    double x, epsilon = 1e-10;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    if( ( XYs = ptwXY_new( ptwXY_interpolationLinLin, 4, 1.e-3, 10, 10, &status, 0 ) ) == NULL ) 
            nfu_printErrorMsg( "ERROR %s: XYs creation, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    for( i = 0, x = 1; i < 5; i++, x += epsilon / 8 ) {
        if( ( status = ptwXY_setValueAtX( XYs, x, ( x - 1 ) / epsilon ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    for( ; i < 50; i++, x += 2 * epsilon ) {
        if( ( status = ptwXY_setValueAtX( XYs, x, ( x - 1 ) / epsilon ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkMerge( XYs, epsilon, 45 );

    for( i = 0, x = 1; i < 45; i++, x += 2 * epsilon ) {
        if( ( status = ptwXY_setValueAtX( XYs, x, ( x - 1 ) / epsilon ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    for( ; i < 50; i++, x += epsilon / 8 ) {
        if( ( status = ptwXY_setValueAtX( XYs, x, ( x - 1 ) / epsilon ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkMerge( XYs, epsilon, 46 );

    for( i = 0, x = 1; i < 15; i++, x += 2 * epsilon ) {
        if( ( status = ptwXY_setValueAtX( XYs, x, ( x - 1 ) / epsilon ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    for( ; i < 30; i++, x += epsilon / 8 ) {
        if( ( status = ptwXY_setValueAtX( XYs, x, ( x - 1 ) / epsilon ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    for( ; i < 50; i++, x += 2 * epsilon ) {
        if( ( status = ptwXY_setValueAtX( XYs, x, ( x - 1 ) / epsilon ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkMerge( XYs, epsilon, 36 );

    for( i = 0, x = 1; i < 100; i++, x += epsilon / 8 ) {
        if( ( status = ptwXY_setValueAtX( XYs, x, i ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkMerge( XYs, epsilon, 13 );

    ptwXY_free( XYs );
    exit( errCount );
}
/*
************************************************************
*/
static int checkMerge( ptwXYPoints *XYs, double epsilon, int n ) {

    int errCount = 0;
    nfu_status status;

    if( verbose ) printf( "# epsilon = %e\n", epsilon );
    printIfVerbose( XYs );
    if( ( status = ptwXY_mergeClosePoints( XYs, epsilon ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: ptwXY_mergeClosePoints 1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    printIfVerbose( XYs );
    if( ptwXY_length( XYs ) != n ) {
        nfu_printMsg( "ERROR %s: ptwXY_length( XYs ) = %d != n = %d", __FILE__, ptwXY_length( XYs ), n );
        errCount++;
    }
    ptwXY_clear( XYs );
    return( errCount );
}
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    printf( "# length = %d\n", (int) data->length );
    ptwXY_simpleWrite( data, stdout, fmtXY );
    printf( "\n\n" );
}
