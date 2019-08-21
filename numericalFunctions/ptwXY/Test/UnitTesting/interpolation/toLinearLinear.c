/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include <ptwXY.h>
#include <ptwXY_utilities.h>

static int verbose = 0;
char fmt[] = "%22.14e %22.14e\n";

void printMsg( const char *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, xlog = 0, ylog = 0, errCount;
    ptwXYPoints *pSparse, *pDense, *pLinear;
    double accuracy = 1e-3, xys[2*2] = { 1.0, 2.0, 10.0, 100.0 };
    nfu_status status;
    FILE *ff;
    ptwXY_interpolation interpolation = ptwXY_interpolationLinLin;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-xlog", argv[iarg] ) == 0 ) {
            xlog = 1; }
        else if( strcmp( "-ylog", argv[iarg] ) == 0 ) {
            ylog = 1; }
        else {
            printMsg( "Error %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s%s%s\n", __FILE__, xlog ? " -xlog" : "", ylog ? " -ylog" : "" );
    
    nfu_setMemoryDebugMode( 0 );

    if( xlog ) {
        interpolation = ptwXY_interpolationLogLin;
        if( ylog ) interpolation = ptwXY_interpolationLogLog; }
    else if( ylog ) {
        interpolation = ptwXY_interpolationLinLog;
    }

    if( ( pSparse = ptwXY_create( interpolation, 5, accuracy, 10, 10,    2, xys, &status, 0 ) ) == NULL ) 
        printMsg( "pSparse creation: status = %d: %s", status, nfu_statusMessage( status ) );
    if( ( pDense = ptwXY_clone( pSparse, &status ) ) == NULL ) printMsg( "pDense creation: status = %d: %s", status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_thicken( pDense, 4000, 0., 1 + 2.e-4 ) ) != nfu_Okay ) printMsg( "thicken: status = %d: %s", status, nfu_statusMessage( status ) );

    ff = fopen( "curve_sparse.dat", "w" );
    fprintf( ff, "# xlog = %d\n", xlog );
    fprintf( ff, "# ylog = %d\n", ylog );
    fprintf( ff, "# accuracy = %e\n", accuracy );
    fprintf( ff, "# length = %d\n", (int) pSparse->length );
    ptwXY_simpleWrite( pSparse, ff, fmt );
    fclose( ff );

    ff = fopen( "curve_dense.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) pDense->length );
    ptwXY_simpleWrite( pDense, ff, fmt );
    fclose( ff );

    if( ( pLinear = ptwXY_toOtherInterpolation( pSparse, ptwXY_interpolationLinLin, accuracy, &status ) ) == NULL ) 
        printMsg( "pLinear creation: status = %d: %s", status, nfu_statusMessage( status ) );

    ff = fopen( "curve_linear.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) pLinear->length );
    ptwXY_simpleWrite( pLinear, ff, fmt );
    fclose( ff );

    if( ( errCount = nfu_ptwXY_cmp( pDense, pLinear, verbose, accuracy ) ) )
        nfu_printMsg( "Error %s: nfu_ptwXY_cmp found %d differences", __FILE__, errCount );

    ptwXY_free( pSparse );
    ptwXY_free( pDense );
    ptwXY_free( pLinear );

    exit( errCount ? EXIT_FAILURE : EXIT_SUCCESS );
}
/*
****************************************************************
*/
void printMsg( const char *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    fprintf( stderr, "\n" );
    va_end( args );
    exit( EXIT_FAILURE );
}
