/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ptwXY.h>
#include <nf_utilities.h>

#define nSame 6

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int checkDullEdges( ptwXYPoints *data, double lowerEps, double upperEps );
static int checkDullEdges2( ptwXYPoints *data, double lowerEps, double upperEps, int positiveXOnly );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;
    nfu_status status;
    ptwXYPoints *XY;

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

    if( ( XY = ptwXY_new( ptwXY_interpolationLinLin, 4, 1.e-3, 10, 10, &status, 0 ) ) == NULL ) 
            nfu_printErrorMsg( "ERROR %s: XY new, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    for( i = 0; i < nSame; i++ ) {
        if( ( status = ptwXY_setValueAtX( XY, 0.2 * i - .5, 0.7 + i + .1 ) ) != nfu_Okay )
                nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX 1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkDullEdges( XY, 1e-10, 1e-10 );
    ptwXY_neg( XY );
    errCount += checkDullEdges( XY, 1e-10, 1e-10 );

    ptwXY_clear( XY );
    for( i = 0; i < nSame; i++ ) {
        if( ( status = ptwXY_setValueAtX( XY, 0.2 * i, 0.7 + i + .1 ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX 1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkDullEdges( XY, 1e-10, 1e-10 );

    ptwXY_clear( XY );
    for( i = 0; i < nSame; i++ ) {
        if( ( status = ptwXY_setValueAtX( XY, -0.2 * i, 0.7 + i + .1 ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX 1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkDullEdges( XY, 1e-10, 1e-10 );

    ptwXY_free( XY );

    exit( errCount );
}
/*
************************************************************
*/
static int checkDullEdges( ptwXYPoints *data, double lowerEps, double upperEps ) {

    int errCount;

    errCount = checkDullEdges2( data, lowerEps, upperEps, 0 );
    errCount += checkDullEdges2( data,  lowerEps, upperEps, 1 );
    errCount += checkDullEdges2( data, -lowerEps, upperEps, 0 );
    errCount += checkDullEdges2( data, -lowerEps, upperEps, 1 );
    errCount += checkDullEdges2( data,  lowerEps, -upperEps, 0 );
    errCount += checkDullEdges2( data,  lowerEps, -upperEps, 1 );
    errCount += checkDullEdges2( data, -lowerEps, -upperEps, 0 );
    errCount += checkDullEdges2( data, -lowerEps, -upperEps, 1 );
    return( errCount );
}
/*
************************************************************
*/
static int checkDullEdges2( ptwXYPoints *data, double lowerEps, double upperEps, int positiveXOnly ) {

    int errCount = 0;
    ptwXYPoints *dullEdges;
    nfu_status status;

    printIfVerbose( data );
    if( verbose ) {
        printf( "# lowerEps = %.14e\n", lowerEps );
        printf( "# upperEps = %.14e\n", upperEps );
        printf( "# positiveXOnly = %d\n", positiveXOnly );
    }
    if( ( dullEdges = ptwXY_clone( data, &status ) ) == NULL )
            nfu_printErrorMsg( "ERROR %s: data, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_dullEdges( dullEdges, lowerEps, upperEps, positiveXOnly ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_dullEdges, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    printIfVerbose( dullEdges );

    ptwXY_free( dullEdges );

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
