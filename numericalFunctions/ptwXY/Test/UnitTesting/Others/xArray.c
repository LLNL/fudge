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

#define nXs 99

static int verbose = 0;
static char *fmtX = "%19.12e\n", *fmtXY = "%19.12e %19.12e\n";
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0;
    nfu_status status;
    ptwXYPoints *XYs;
    ptwXPoints *Xs, *xArray;
    double x;

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

    if( ( Xs = ptwX_new( nXs, &status ) ) == NULL ) 
            nfu_printErrorMsg( "ERROR %s: Xs creation, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    if( ( XYs = ptwXY_new( ptwXY_interpolationLinLin, 4, 1.e-3, nXs / 2, 10, &status, 0 ) ) == NULL ) 
            nfu_printErrorMsg( "ERROR %s: XYs creation, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    for( i = 0; i < nXs; i++ ) {
        x = 2.1 * i + 0.1;
        Xs->points[i] = x;
        Xs->length++;
        if( ( status = ptwXY_setValueAtX( XYs, x, i * i ) ) != nfu_Okay ) 
            nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }

    if( verbose ) {
        printf( "# length = %d\n", (int) XYs->length );
        ptwXY_simpleWrite( XYs, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) Xs->length );
        for( i = 0; i < ptwX_length( Xs ); i++ ) printf( fmtX, ptwX_getPointAtIndex_Unsafely( Xs, i ) );
        printf( "\n\n" );
    }

    if( ( xArray = ptwXY_getXArray( XYs, &status ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: ptwXY_getXArray call, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    if( verbose ) {
        printf( "# length = %d\n", (int) nXs );
        for( i = 0; i < nXs; i++ ) printf( fmtX, ptwX_getPointAtIndex_Unsafely( xArray, i ) );
    }

    if( ptwX_length( Xs ) != ptwX_length( xArray ) ) nfu_printErrorMsg( "ERROR %s: length not the same", __FILE__ );
    for( i = 0; i < nXs; i++ ) {
        if( Xs->points[i] != xArray->points[i] ) nfu_printErrorMsg( "ERROR %s: at index %d, %e != %e", __FILE__, i, Xs->points[i], xArray->points[i] );
    }

    ptwXY_free( XYs );
    ptwX_free( Xs );
    ptwX_free( xArray );

    exit( EXIT_SUCCESS );
}
