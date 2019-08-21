/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>

#include <ptwXY.h>

#define size 1001

static int verbose = 0;

double getDouble( const char const *s );
void printMsg( const char *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int i, n, iarg, echo = 0;
    ptwXYPoints *u, *v, *y, *e;
    double xy2[2*2], x, dx, u1, v1;
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

    xy2[0] = 0.0;
    xy2[1] = 1.0;
    xy2[2] = 1.0;
    xy2[3] = -0.2;
    if( argc == 5 ) {
        xy2[0] = getDouble( argv[1] );
        xy2[1] = getDouble( argv[2] );
        xy2[2] = getDouble( argv[3] );
        xy2[3] = getDouble( argv[4] );
    }

    if( ( u = ptwXY_create( ptwXY_interpolationLinLin, 5, 1e-3, 10, 10, 2, xy2, &status, 0 ) ) == NULL ) printMsg( "u creation: status = %d: %s", 
        status, nfu_statusMessage( status ) );
    xy2[3] = -1.0;
    if( ( v = ptwXY_create( ptwXY_interpolationLinLin, 5, 1e-3, 10, 10, 2, xy2, &status, 0 ) ) == NULL ) printMsg( "v creation: status = %d: %s", 
        status, nfu_statusMessage( status ) );

    ff = fopen( "curve_u.dat", "w" );
    ptwXY_simpleWrite( u, ff, "%g, %g\n" );
    fclose( ff );

    ff = fopen( "curve_v.dat", "w" );
    ptwXY_simpleWrite( v, ff, "%g, %g\n" );
    fclose( ff );

    if( ( y = ptwXY_mul2_ptwXY( u, v, &status ) ) == NULL ) printMsg( "y creation: status = %d: %s", status, nfu_statusMessage( status ) );

    ff = fopen( "u_times_v.dat", "w" );
    ptwXY_simpleWrite( y, ff, "%g, %g\n" );
    fclose( ff );

    n = 1000;
    x = ptwXY_getXMin( y );
    dx = ( ptwXY_getXMax( y ) - x ) / n;
    if( ( e = ptwXY_new( ptwXY_interpolationLinLin, 5, 1e-3, 10, 10, &status, 0 ) ) == NULL ) printMsg( "e creation: status = %d: %s", 
        status, nfu_statusMessage( status ) );
    for( i = 0; i < n; i++ ) {
        ptwXY_getValueAtX( u, x, &u1 );
        ptwXY_getValueAtX( v, x, &v1 );
        ptwXY_setValueAtX( e, x, u1 * v1 );
        x += dx;
    }
    if( ptwXY_getXMax( e ) < ptwXY_getXMax( y ) ) {
        x = ptwXY_getXMax( y );
        ptwXY_getValueAtX( u, x, &u1 );
        ptwXY_getValueAtX( v, x, &v1 );
        ptwXY_setValueAtX( e, x, u1 * v1 );
    }
    ff = fopen( "exact.dat", "w" );
    ptwXY_simpleWrite( e, ff, "%g, %g\n" );
    fclose( ff );


    ptwXY_free( u );
    ptwXY_free( v );
    ptwXY_free( y );
    ptwXY_free( e );

    exit( EXIT_SUCCESS );
}
/*
****************************************************************
*/
double getDouble( const char const *s ) {

    double d;
    char *e;

    errno = 0;
    d = strtod( s, &e );
    if( ( *e != 0 ) || ( errno != 0 ) ) printMsg( "could not convert '%s' to double, err = %d, e = %s", s, errno, e );
    return( d );
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
