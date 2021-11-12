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

void printMsg( char const *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int i, n = 100 * 1000, primarySize = 7, secondarySize = 13, iarg, echo = 0;
    unsigned short seed16v[3] = { 1242, 14213, 543 };
    double accuracy = 1e-3, biSectionMax = 3., xMin = -100, xMax = 100, yMin = 0, yMax = 10, r, x, y;
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

    if( ( f = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, biSectionMax, accuracy, primarySize, secondarySize, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );

    for( i = 0; i < n; i++ ) {
        r = drand48( );
        x = r * xMin + ( 1. - r ) * xMax;
        r = drand48( );
        y = r * yMin + ( 1. - r ) * yMax;
        if( ( status = ptwXY_setValueAtX( &smr, f, x, y ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }

    for( i = 0; i < n; i++ ) {
        point = ptwXY_getPointAtIndex_Unsafely( f, i );
        if( i > 0 ) {
            if( x >= point->x ) printMsg( "Error x[%d] = %16e > x[%d] = %16e\n", i - 1, x, i, point->x );
        }
        x = point->x;
    }
    ptwXY_free( f );
    exit( EXIT_SUCCESS );
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
