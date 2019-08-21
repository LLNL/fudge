/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

/*
    This routine tests the insertion of points where a point already exists.
    The number of initial points, currently 7, must be a little smaller than the space allocated for 
    points (currently 10) and overflowPoints (currently 10).
*/
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>

#include <ptwXY.h>

static int verbose = 0;

void printMsg( const char *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    ptwXYPoint *point;
    ptwXYPoints *u;
    double uXY[] = { 3.732, 0, 5.598, 0, 5.679, 0, 6.459, 0, 6.557, 0, 6.838, 0, 7.1138, 0 };
    int i, j, n, nuXY = sizeof( uXY ) / ( 2 * sizeof( double ) ), iarg, echo = 0;
    nfu_status status;
    FILE *ff;
    char *fmt = "%7.4f, %5.1f\n";

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            printMsg( "Error %s: unsupported option = '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    nfu_setMemoryDebugMode( 0 );

    if( ( u = ptwXY_create( ptwXY_interpolationLinLin, 5, 1e-3, 10, 10, nuXY, uXY, &status, 0 ) ) == NULL ) printMsg( "u creation: status = %d: %s", 
        status, nfu_statusMessage( status ) );

    if( ( ff = fopen( "curve_u.dat", "w" ) ) == NULL ) printMsg( "Error %s: could not open output file", __FILE__ );
    ptwXY_simpleWrite( u, ff, fmt );
    fprintf( ff, "\n" );
    if( verbose ) ptwXY_simpleWrite( u, stdout, fmt );
    if( verbose ) printf( "\n" );

/*
 Overwrite existing y values for each point. This should not change the number of points.
*/
    for( i = 0; i < nuXY; i++ ) ptwXY_setValueAtX( u, uXY[2 * i], -i - 1. );
    if( ptwXY_length( u ) != nuXY ) printMsg( "Error %s: ptwXY_length( u ) = %d != nuXY = %d", __FILE__, (int) ptwXY_length( u ), nuXY );
    ptwXY_simpleWrite( u, ff, fmt );
    fprintf( ff, "\n" );
    if( verbose ) ptwXY_simpleWrite( u, stdout, fmt );
    if( verbose ) printf( "\n" );

/*
 Add a point before xMin, between each point and after xMax. This should change the number of points to 2 * nuXY + 1.
*/
    ptwXY_setValueAtX( u, uXY[0] - .1, 4 );
    for( i = 0; i < nuXY - 1; i++ ) ptwXY_setValueAtX( u, 0.5 * ( uXY[2 * i] + uXY[2 * i + 2] ), i + 5. );
    ptwXY_setValueAtX( u, uXY[2 * nuXY - 2] + .1, 11 );
    if( ptwXY_length( u ) != 2 * nuXY + 1 ) printMsg( "Error %s: ptwXY_length( u ) = %d != nuXY = %d", __FILE__, (int) ptwXY_length( u ), nuXY );
    ptwXY_simpleWrite( u, ff, fmt );
    fprintf( ff, "\n" );
    if( verbose ) ptwXY_simpleWrite( u, stdout, fmt );
    if( verbose ) printf( "\n" );
    if( verbose ) ptwXY_showInteralStructure( u, stdout, 1 );

/*
 Adds some points and after each add overwrite all existing y values for each point. The overwrite should not change the number of points.
*/
    for( i = 0; i < nuXY - 1; i++ ) {
        ptwXY_setValueAtX( u, 0.75 * uXY[2 * i] + 0.25 * uXY[2 * i + 2], i + 55. );
        n = ptwXY_length( u );
        for( j = 0; j < n; j++ ) {
            point = ptwXY_getPointAtIndex( u, j );
            ptwXY_setValueAtX( u, point->x, i + j );
            if( ptwXY_length( u ) != n ) printMsg( "Error %s: ptwXY_length( u ) = %d != n = %d", __FILE__, (int) ptwXY_length( u ), n );
        }
    }
    ptwXY_simpleWrite( u, ff, fmt );
    if( verbose ) ptwXY_simpleWrite( u, stdout, fmt );
    if( verbose ) printf( "\n" );
    if( verbose ) ptwXY_showInteralStructure( u, stdout, 1 );

    ptwXY_free( u );
    fclose( ff );

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
