/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
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
void printMsg( char const *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int i, n, iarg, echo = 0;
    ptwXYPoints *u, *v, *y, *e;
    double xy2[2*2], x, dx, u1, v1, domainMaxY, domainMaxE;
    FILE *ff;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

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

    if( ( u = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 5, 1e-3, 10, 10, 2, xy2, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    xy2[3] = -1.0;
    if( ( v = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 5, 1e-3, 10, 10, 2, xy2, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );

    ff = fopen( "curve_u.dat", "w" );
    ptwXY_simpleWrite( u, ff, "%g, %g\n" );
    fclose( ff );

    ff = fopen( "curve_v.dat", "w" );
    ptwXY_simpleWrite( v, ff, "%g, %g\n" );
    fclose( ff );

    if( ( y = ptwXY_mul2_ptwXY( &smr, u, v ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );

    ff = fopen( "u_times_v.dat", "w" );
    ptwXY_simpleWrite( y, ff, "%g, %g\n" );
    fclose( ff );

    n = 1000;
    if( ptwXY_domainMin( &smr, y, &x  ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_domainMax( &smr, y, &domainMaxY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    dx = ( domainMaxY - x ) / n;
    if( ( e = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 5, 1e-3, 10, 10, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i = 0; i < n; i++ ) {
        if( ptwXY_getValueAtX( &smr, u, x, &u1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
        if( ptwXY_getValueAtX( &smr, v, x, &v1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
        if( ptwXY_setValueAtX( &smr, e, x, u1 * v1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
        x += dx;
    }
    if( ptwXY_domainMax( &smr, e, &domainMaxE ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( domainMaxE < domainMaxY ) {
        x = domainMaxY;
        if( ptwXY_getValueAtX( &smr, u, x, &u1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
        if( ptwXY_getValueAtX( &smr, v, x, &v1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
        if( ptwXY_setValueAtX( &smr, e, x, u1 * v1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
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
void printMsg( char const *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    fprintf( stderr, "\n" );
    va_end( args );
    exit( EXIT_FAILURE );
}
