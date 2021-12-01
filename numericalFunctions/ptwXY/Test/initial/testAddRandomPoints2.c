/*
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

/*
    This routine test putting random x, y values into a ptwXYPoints instance. Because primarySize = 7 and secondarySize = 13, 
    this routine can take some time.
*/
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdarg.h>

#include <nfut_utilities.h>
#include <ptwXY.h>

static int verbose = 0;

int xCompare( void const *, void const * );
void printMsg( char const *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int i3, n = 10 * 1000, i1, i2, iarg, echo = 0;
    unsigned short seed16v[3] = { 1242, 14213, 543 };
    double accuracy = 1e-3, biSectionMax = 3., xMin = -100, xMax = 100, yMin = 0, yMax = 10, r, x, y, *points, *p;
    nfu_status status;
    ptwXYPoints *f;
    ptwXYPoint *point;
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

    seed48( seed16v );

    if( ( points = malloc( 2 * sizeof( double ) * n ) ) == NULL ) printMsg( "Error allocating points\n" );
    for( i1 = ptwXY_minimumSize; i1 < ptwXY_minimumSize + 10; i1++ ) {
        for( i2 = ptwXY_minimumOverflowSize; i2 < ptwXY_minimumOverflowSize + 10; i2++ ) {
            if( ( f = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, biSectionMax, accuracy, i1, i2, 0 ) ) == NULL ) 
                nfut_printSMRErrorExit2p( &smr, "Via." );

            for( i3 = 0, p = points; i3 < n; i3++ ) {
                r = drand48( );
                x = r * xMin + ( 1. - r ) * xMax;
                *(p++) = x;
                r = drand48( );
                y = r * yMin + ( 1. - r ) * yMax;
                *(p++) = y;
                if( ( status = ptwXY_setValueAtX( &smr, f, x, y ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
            }
            qsort( points, n, 2 * sizeof( double ), xCompare );
            for( i3 = 0, p = points; i3 < n; i3++ ) {
                point = ptwXY_getPointAtIndex_Unsafely( f, i3 );
                x = *(p++);
                y = *(p++);
                if( ( x != point->x ) || ( y != point->y ) ) 
                    printMsg( "Error x,y = %16e, %16e != point = %16e, %16e: i1 = %d, i2 = %d, i3 = %d\n", 
                            x, y, point->x, point->y, i1, i2, i3 );
            }
            ptwXY_free( f );
        }
    }

    exit( EXIT_SUCCESS );
}
/*
****************************************************************
*/
int xCompare( void const *p1, void const *p2 ) {

    double *x1 = (double *) p1, *x2 = (double *) p2;

    if( *x1 > *x2 ) return( 1 );
    if( *x1 == *x2 ) return( 0 );
    return( -1 );
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
