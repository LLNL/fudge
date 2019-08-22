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
#include <math.h>

#include <ptwXY.h>
#include <nf_utilities.h>

#define nSame 6

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int checkDullEdges( ptwXYPoints *data, double lowerEps, double upperEps );
static int checkDullEdges2( ptwXYPoints *data, double lowerEps, double upperEps, int positiveXOnly );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;
    nfu_status status;
    ptwXYPoints *XY;

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

    if( ( XY = ptwXY_new( ptwXY_interpolationLinLin, 4, 1.e-3, 10, 10, &status, 0 ) ) == NULL ) 
            nfu_printErrorMsg( "ERROR %s: XY new, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    for( i = 0; i < nSame; i++ ) {
        if( ( status = ptwXY_setValueAtX( XY, 0.2 * i - .5, 0.7 + i + .1 ) ) != nfu_Okay )
                nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX 1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkDullEdges( XY, 1e-10, 1e-10 );
    ptwXY_neg( XY );
    errCount += checkDullEdges( XY, 1e-10, 1e-10 );

    ptwXY_clear( XY );
    for( i = 0; i < nSame; i++ ) {
        if( ( status = ptwXY_setValueAtX( XY, 0.2 * i, 0.7 + i + .1 ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX 1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkDullEdges( XY, 1e-10, 1e-10 );

    ptwXY_clear( XY );
    for( i = 0; i < nSame; i++ ) {
        if( ( status = ptwXY_setValueAtX( XY, -0.2 * i, 0.7 + i + .1 ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX 1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkDullEdges( XY, 1e-10, 1e-10 );

    ptwXY_free( XY );

    exit( errCount );
}
/*
************************************************************
*/
static int checkDullEdges( ptwXYPoints *data, double lowerEps, double upperEps ) {

    int errCount;

    errCount = checkDullEdges2( data, lowerEps, upperEps, 0 );
    errCount += checkDullEdges2( data,  lowerEps, upperEps, 1 );
    errCount += checkDullEdges2( data, -lowerEps, upperEps, 0 );
    errCount += checkDullEdges2( data, -lowerEps, upperEps, 1 );
    errCount += checkDullEdges2( data,  lowerEps, -upperEps, 0 );
    errCount += checkDullEdges2( data,  lowerEps, -upperEps, 1 );
    errCount += checkDullEdges2( data, -lowerEps, -upperEps, 0 );
    errCount += checkDullEdges2( data, -lowerEps, -upperEps, 1 );
    return( errCount );
}
/*
************************************************************
*/
static int checkDullEdges2( ptwXYPoints *data, double lowerEps, double upperEps, int positiveXOnly ) {

    int errCount = 0;
    ptwXYPoints *dullEdges;
    nfu_status status;

    printIfVerbose( data );
    if( verbose ) {
        printf( "# lowerEps = %.14e\n", lowerEps );
        printf( "# upperEps = %.14e\n", upperEps );
        printf( "# positiveXOnly = %d\n", positiveXOnly );
    }
    if( ( dullEdges = ptwXY_clone( data, &status ) ) == NULL )
            nfu_printErrorMsg( "ERROR %s: data, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_dullEdges( dullEdges, lowerEps, upperEps, positiveXOnly ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_dullEdges, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    printIfVerbose( dullEdges );

    ptwXY_free( dullEdges );

    return( errCount );
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
