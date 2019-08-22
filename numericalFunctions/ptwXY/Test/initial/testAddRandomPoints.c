/*
# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
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

#include <ptwXY.h>

static int verbose = 0;

void printMsg( const char *fmt, ... );
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

    if( ( f = ptwXY_new( ptwXY_interpolationLinLin, biSectionMax, accuracy, primarySize, secondarySize, &status, 0 ) ) == NULL ) 
        printMsg( "u creation: status = %d: %s", status, nfu_statusMessage( status ) );

    for( i = 0; i < n; i++ ) {
        r = drand48( );
        x = r * xMin + ( 1. - r ) * xMax;
        r = drand48( );
        y = r * yMin + ( 1. - r ) * yMax;
        if( ( status = ptwXY_setValueAtX( f, x, y ) ) != nfu_Okay ) printMsg( "Error setting x, y = %16e, %16e at index %d with status = %d\n", 
            x, y, i, status );
    }

    for( i = 0; i < n; i++ ) {
        point = ptwXY_getPointAtIndex( f, i );
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
void printMsg( const char *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    fprintf( stderr, "\n" );
    va_end( args );
    exit( EXIT_FAILURE );
}
