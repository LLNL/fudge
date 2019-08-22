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

#define allocatedSize 100

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int compareXYs( ptwXYPoints *XY1, ptwXYPoints *XY2 );
static int compareXYsToCList( ptwXYPoints *XY1, int64_t nPoints, double *xy );
static int compareValues( int64_t i, double x1, double y1, double x2, double y2 );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;
    int64_t returnedPoints;
    double x, xy[2 * allocatedSize];
    nfu_status status;
    ptwXYPoints *XYSrc, *XYDesc;

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

    if( ( XYSrc = ptwXY_new( ptwXY_interpolationLinLin, 4, 1.e-3, 10, 10, &status, 0 ) ) == NULL ) 
            nfu_printErrorMsg( "ERROR %s: XYSrc creation, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    for( i = 0, x = 1; i < 45; i++, x += 1.1 ) {
        if( ( status = ptwXY_setValueAtX( XYSrc, x, x * x ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }

    if( ( XYDesc = ptwXY_new( ptwXY_interpolationLinLin, 4, 1.e-3, 10, 10, &status, 0 ) ) == NULL ) 
            nfu_printErrorMsg( "ERROR %s: XYDesc creation, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );


    if( ( status = ptwXY_copy( XYDesc, XYSrc ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    errCount += compareXYs( XYDesc, XYSrc );
    ptwXY_free( XYDesc );

    if( ( status = ptwXY_copyToC_XY( XYSrc, 0, ptwXY_length( XYSrc ), allocatedSize, &returnedPoints, xy ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_copyToC_XY, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    errCount += compareXYsToCList( XYSrc, returnedPoints, xy );

    ptwXY_free( XYSrc );
    exit( errCount );
}
/*
************************************************************
*/
static int compareXYs( ptwXYPoints *XY1, ptwXYPoints *XY2 ) {

    int errCount = 0;
    int64_t i, n = ptwXY_length( XY1 );
    nfu_status status;
    double x1, y1, x2, y2;

    printIfVerbose( XY1 );
    printIfVerbose( XY2 );
    if( n != ptwXY_length( XY2 ) ) {
        errCount++;
        nfu_printMsg( "ERROR %s: compareXYs, len( XY1 ) = %d != len( XY2 ) = %d", __FILE__, (int) n, (int) ptwXY_length( XY2 ) ); }
    else {
        for( i = 0; i < n; i++ ) {
            if( ( status = ptwXY_getXYPairAtIndex( XY1, i, &x1, &y1 ) ) != nfu_Okay )
                nfu_printMsg( "ERROR %s: ptwXY_getXYPairAtIndex( XY1 ), status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
            if( ( status = ptwXY_getXYPairAtIndex( XY2, i, &x2, &y2 ) ) != nfu_Okay )
                nfu_printMsg( "ERROR %s: ptwXY_getXYPairAtIndex( XY2 ), status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
            errCount += compareValues( i, x1, y1, x2, y2 );
        }
    }

    return( errCount );
}
/*
************************************************************
*/
static int compareXYsToCList( ptwXYPoints *XY1, int64_t nPoints, double *xy ) {

    int errCount = 0;
    int64_t i, n = ptwXY_length( XY1 );
    nfu_status status;
    double x1, y1, *p;

    if( n != nPoints ) {
        errCount++;
        nfu_printMsg( "ERROR %s: compareXYsToCList, len( XY1 ) = %d != nPoints = %d", __FILE__, (int) n, (int) nPoints ); }
    else {
        for( i = 0, p = xy; i < n; i++, p += 2 ) {
            if( ( status = ptwXY_getXYPairAtIndex( XY1, i, &x1, &y1 ) ) != nfu_Okay )
                nfu_printMsg( "ERROR %s: ptwXY_getXYPairAtIndex( XY1 ), status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
            errCount += compareValues( i, x1, y1, *p, p[1] );
        }
    }

    return( errCount );
}
/*
************************************************************
*/
static int compareValues( int64_t i, double x1, double y1, double x2, double y2 ) {

    if( ( x1 != x2 ) || ( y1 != y2 ) ) {
        nfu_printMsg( "ERROR %s: at index %3d ( x1 = %.17e != x2 = %.17e ) or ( y1 = %.17e != y2 = %.17e )", __FILE__, (int) i, x1, x2, y1, y2 );
        return( 1 );
    }
    return( 0 );
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
