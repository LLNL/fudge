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

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include <ptwXY.h>
#include <nf_utilities.h>

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int checkTrim( ptwXYPoints *data );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, n, iarg, echo = 0, errCount = 0;
    nfu_status status;
    ptwXYPoints *XYs;
    double x;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    if( ( XYs = ptwXY_new( ptwXY_interpolationLinLin, 4, 1.e-3, 10, 10, &status, 0 ) ) == NULL ) 
            nfu_printErrorMsg( "ERROR %s: XYs creation, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    n = 5;
    for( i = 0; i < n; i++ ) {
        x = .2 * i + .5;
        if( ( status = ptwXY_setValueAtX( XYs, x, 0. ) ) != nfu_Okay )
                    nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkTrim( XYs );

    n += 7;
    for( ; i < n; i++ ) {
        x = .2 * i + .5;
        if( ( status = ptwXY_setValueAtX( XYs, x, 2 * i - 20. ) ) != nfu_Okay )
                    nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    n += 3;
    for( ; i < n; i++ ) {
        x = .2 * i + .5;
        if( ( status = ptwXY_setValueAtX( XYs, x, 0. ) ) != nfu_Okay )
                    nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkTrim( XYs );

    ptwXY_clear( XYs );
    for( i = 0; i < 12; i++ ) {
        x = .2 * i + .5;
        if( ( status = ptwXY_setValueAtX( XYs, x, 2 * i - 20. ) ) != nfu_Okay )
                    nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkTrim( XYs );

    ptwXY_free( XYs );

    exit( errCount );
}
/*
************************************************************
*/
static int checkTrim( ptwXYPoints *data ) {

    int allPointsZero = 0;
    int64_t i, i1, errCount;
    nfu_status status = nfu_Okay;
    ptwXYPoints *trimmed;
    ptwXYPoint *point1, *point2;

    printIfVerbose( data );
    if( ( trimmed = ptwXY_clone( data, &status ) ) == NULL )
            nfu_printErrorMsg( "ERROR %s: data clone, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_trim( trimmed ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_trim, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    printIfVerbose( trimmed );

    for( i1 = 0; i1 < data->length; i1++ ) {
        point1 = ptwXY_getPointAtIndex_Unsafely( data, i1 );
        if( point1->y != 0. ) break;
    }
    if( i1 == data->length ) {      /* All points are zero. */
        allPointsZero = 1;
        i1 = 0; }
    else if( i1 > 0 ) {
        i1--;
    }
    for( i = 0; ( i < trimmed->length ) && ( i1 < data->length ); i++, i1++ ) {
        point1 = ptwXY_getPointAtIndex_Unsafely( data, i1 );
        point2 = ptwXY_getPointAtIndex_Unsafely( trimmed, i );
        if( ( point1->x != point2->x ) || ( point1->y != point2->y ) ) {
                nfu_printErrorMsg( "ERROR %s: ptwXY_trim at index %d, (x,y) = (%e %e) != (%e %e)", __FILE__, (int) i, 
                    point1->x, point1->y, point2->x, point2->y );
                errCount++;
        }
        if( allPointsZero ) i1 = data->length - 2;
    }

    ptwXY_free( trimmed );

    return( status );
}
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    printf( "# length = %d\n", (int) data->length );
    ptwXY_simpleWrite( data, stdout, fmtXY );
    printf( "\n\n" );
}
