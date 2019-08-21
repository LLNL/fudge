/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include <ptwXY.h>

#define size 1001

static int verbose = 0;
char fmt[] = "%22.14e %22.14e\n";

void printMsg( const char *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0;
    ptwXYPoints *l2, *lLog, *lLin;
    double xy[2*size], xy2[2*2], a, x, y, f;
    nfu_status status;
    FILE *ff;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else {
            printMsg( "Error %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );
    
    nfu_setMemoryDebugMode( 0 );

    xy2[0] = 1.0;
    xy2[1] = 1.0;
    xy2[2] = 100.0;
    xy2[3] = 10.0;

    a = log( xy2[3] / xy2[1] ) / log( xy2[2] / xy2[0] );
    f = pow( xy2[2] / xy2[0], 1. / ( size - 1 ) );
    x = xy2[0];

    xy[0] = xy2[0];
    xy[1] = xy2[1];
    for( i = 1; i < size - 1; i++ ) {
        x *= f;
        y = pow( x / xy2[0], a );
        xy[2*i] = x;
        xy[2*i+1] = y;
    }
    xy[2 * size - 2] = xy2[2];
    xy[2 * size - 1] = xy2[3];

    if( ( l2    = ptwXY_create( ptwXY_interpolationLogLog, 5, 1e-3, 10, 10,    2, xy2, &status, 0 ) ) == NULL ) printMsg( "l2 creation: status = %d: %s", 
        status, nfu_statusMessage( status ) );
    if( ( lLog = ptwXY_create( ptwXY_interpolationLogLog, 5, 1e-3, 25, 10, size,  xy, &status, 0 ) ) == NULL ) printMsg( "lLog creation: status = %d: %s", 
        status, nfu_statusMessage( status ) );

    ff = fopen( "curve_u_loglog.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) l2->length );
    ptwXY_simpleWrite( l2, ff, fmt );
    fclose( ff );

    ff = fopen( "curve_u_loglog_dense.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) lLog->length );
    ptwXY_simpleWrite( lLog, ff, fmt );
    fclose( ff );

    if( ( lLin = ptwXY_toOtherInterpolation( l2, ptwXY_interpolationLinLin, 1e-3, &status ) ) == NULL ) 
        printMsg( "lLin creation: status = %d: %s", status, nfu_statusMessage( status ) );

    ff = fopen( "curve_u_interpolatedToLinear.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) lLin->length );
    ptwXY_simpleWrite( lLin, ff, fmt );
    fclose( ff );

    ptwXY_free( l2 );
    ptwXY_free( lLog );
    ptwXY_free( lLin );

    exit( EXIT_SUCCESS );
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
