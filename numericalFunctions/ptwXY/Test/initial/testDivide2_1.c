/*
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>

#include <nfut_utilities.h>
#include <ptwXY.h>

#define size 1001

static int verbose = 0;

double getDouble( char const *s );
ptwXYPoints *randomUV( statusMessageReporting *smr );
void printMsg( char const *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int doRandom = 0, iarg, echo = 0;
    int64_t i, j, n;
    ptwXYPoints *u, *v, *y, *e;
    double xyPoints[2*2], x, dx, u1, v1;
    ptwXYPoint *xy1, *xy2;
    nfu_status status;
    FILE *ff;
    char fmt[] = "%.14e %.14e\n";
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-r", argv[iarg] ) == 0 ) {
            doRandom = 1; }
        else {
            printMsg( "Error %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    nfu_setMemoryDebugMode( 0 );

    xyPoints[0] = 0.0;
    xyPoints[1] = 1.0;
    xyPoints[2] = 1.0;
    xyPoints[3] = -0.2;
    if( argc == 5 ) {
        xyPoints[0] = getDouble( argv[1] );
        xyPoints[1] = getDouble( argv[2] );
        xyPoints[2] = getDouble( argv[3] );
        xyPoints[3] = getDouble( argv[4] );
    }

    if( doRandom ) {
        u = randomUV( &smr );
        v = randomUV( &smr ); }
    else {
        if( ( u = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10, 1e-3, 10, 10, 2, xyPoints, 0 ) ) == NULL )
            nfut_printSMRErrorExit2p( &smr, "Via." );
        xyPoints[3] = -1.0;
        if( ( v = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10, 1e-3, 10, 10, 2, xyPoints, 0 ) ) == NULL )
            nfut_printSMRErrorExit2p( &smr, "Via." );
    }

    ff = fopen( "curve_u.dat", "w" );
    ptwXY_simpleWrite( u, ff, fmt );
    fclose( ff );

    ff = fopen( "curve_v.dat", "w" );
    ptwXY_simpleWrite( v, ff, fmt );
    fclose( ff );

    if( ( y = ptwXY_div_ptwXY( &smr, u, v, 1 ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( status = ptwXY_simpleCoalescePoints( &smr, y ) ) != nfu_Okay )
        nfut_printSMRErrorExit2p( &smr, "Via." );

    ff = fopen( "u_divide_v.dat", "w" );
    ptwXY_simpleWrite( y, ff, fmt );
    fclose( ff );

    if( ( e = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 5, 1e-3, 10, 10, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    n = ptwXY_length( &smr, y );
    xy1 = ptwXY_getPointAtIndex_Unsafely( y, 0 );
    for( i = 1; i < n; i++ ) {
        xy2 = ptwXY_getPointAtIndex_Unsafely( y, i );
        x = xy1->x;
        if( xy1->y * xy2->y < 0. ) {
            if( ptwXY_getValueAtX( &smr, u, x, &u1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
            if( ptwXY_getValueAtX( &smr, v, x, &v1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
            if( v1 != 0. ) {
                if( ptwXY_setValueAtX( &smr, e, x, u1 / v1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
            } }
        else {
            dx = ( xy2->x - xy1->x ) / 5;
            for( j = 0; j < 5; j++ ) {
                if( ptwXY_getValueAtX( &smr, u, x, &u1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
                if( ptwXY_getValueAtX( &smr, v, x, &v1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
                if( v1 != 0. ) {
                    if( ptwXY_setValueAtX( &smr, e, x, u1 / v1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
                }
                x += dx;
            }
        }
        xy1 = xy2;
    }
    if( ptwXY_getValueAtX( &smr, u, xy2->x, &u1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_getValueAtX( &smr, v, xy2->x, &v1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( v1 != 0. ) {
        if( ptwXY_setValueAtX( &smr, e, xy2->x, u1 / v1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    ff = fopen( "exactDivide.dat", "w" );
    ptwXY_simpleWrite( e, ff, fmt );
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
double getDouble( char const *s ) {

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
ptwXYPoints *randomUV( statusMessageReporting *smr ) {

    int64_t i;
    double x, y;
    nfu_status status;
    ptwXYPoints *f;

    if( ( f = ptwXY_new( smr, ptwXY_interpolationLinLin, NULL, 5, 1e-3, 10, 10, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( smr, "Via." );
    for( i = 0, x = 0, y = 0; i < size; i++ ) {
        x += drand48( );
        y += drand48( ) - 0.5;
        if( ( status = ptwXY_setValueAtX( smr, f, x, y ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    }
    return( f );
}
/*
****************************************************************
*/
void printMsg( char const *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    fprintf( stderr, "\n" );
    va_end( args );
    exit( EXIT_FAILURE );
}
