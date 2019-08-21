/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>


static int verbose = 0;
static char *fmtXY = "%25.15e %25.15e\n";
static FILE *infoF;

static void flatInterpolationToLinear( ptwXYPoints *XYs );
static void flatInterpolationToLinear2( ptwXYPoints *XYs, double epsm, double epsp );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0, i;
    ptwXYPoints *fineXYs, *coarseXYs;
    nfu_status status;
    double *fineYs, *coarseYs, accuracy = 1e-3;
    double fineXs[] = { -2.0, -1.5, -1.4, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 2.2, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
    double coarseXs[] = { -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0 };
    int nFineXs= sizeof( fineXs ) / ( sizeof( double ) ), nCoarseXs = sizeof( coarseXs ) / ( sizeof( double ) );

    infoF = stdout;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) fprintf( stderr, "%s\n", __FILE__ );

    if( verbose ) fprintf( infoF, "# accuracy = %e\n", accuracy );

    if( ( fineYs = nfu_malloc( nFineXs * sizeof( double ) ) ) == NULL ) nfu_printErrorMsg( "ERROR %s: nfu_malloc-ing fineXYs", __FILE__ );
    for( i = 0; i < nFineXs; i++ ) fineYs[i] = 2 * i + 3;
    fineYs[nFineXs-1] = fineYs[nFineXs-2];
    if( ( fineXYs = ptwXY_new( ptwXY_interpolationFlat, 0., accuracy, 10, 10, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: dataXY creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_setXYDataFromXsAndYs( fineXYs, nFineXs, fineXs, fineYs ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: ptwXY_setXYDataFromXsAndYs status = %d for fineXYs: %s", __FILE__, status, nfu_statusMessage( status ) );
    flatInterpolationToLinear( fineXYs );

    if( ( coarseYs = nfu_malloc( nCoarseXs * sizeof( double ) ) ) == NULL ) nfu_printErrorMsg( "ERROR %s: nfu_malloc-ing nCoarseXs", __FILE__ );
    for( i = 0; i < nCoarseXs ; i++ ) coarseYs[i] = -i + 10;
    coarseYs[nCoarseXs-1] = coarseYs[nCoarseXs-2];
    if( ( coarseXYs = ptwXY_new( ptwXY_interpolationFlat, 0., accuracy, 10, 10, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: coarseXYs creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_setXYDataFromXsAndYs( coarseXYs, nCoarseXs, coarseXs, coarseYs ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: ptwXY_setXYDataFromXsAndYs status = %d for coarseXYs: %s", __FILE__, status, nfu_statusMessage( status ) );

    free( fineYs );
    free( coarseYs );
    ptwXY_free( fineXYs );
    ptwXY_free( coarseXYs );
    exit( errCount );
}
/*
************************************************************
*/
static void flatInterpolationToLinear( ptwXYPoints *XYs ) {

    double dx = 1e-2;

    flatInterpolationToLinear2( XYs, 0., dx );
    flatInterpolationToLinear2( XYs, dx, 0. );
    flatInterpolationToLinear2( XYs, dx, dx );
}
/*
************************************************************
*/
static void flatInterpolationToLinear2( ptwXYPoints *XYs, double epsm, double epsp ) {

    nfu_status status;
    ptwXYPoints *linearXYs;

    if( verbose ) {
        fprintf( infoF, "# epsm = %e\n", epsm );
        fprintf( infoF, "# epsp = %e\n", epsp );
    }
    printIfVerbose( XYs );
    if( ( linearXYs = ptwXY_flatInterpolationToLinear( XYs, epsm, epsp, &status ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: ptwXY_flatInterpolationToLinear: status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    printIfVerbose( linearXYs );
    ptwXY_free( linearXYs );
}
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    fprintf( infoF, "# length = %d\n", (int) ptwXY_length( data ) );
    ptwXY_simpleWrite( data, infoF, fmtXY );
    fprintf( infoF, "\n\n" );
}
