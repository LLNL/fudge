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

static int verbose = 0;
static char *fmtXY = "%17.8e%17.8e\n";

static int checkMutualify( ptwXYPoints *data );
static int checkMutualify2( ptwXYPoints *data, int64_t i1, int64_t i2 );
static int checkMutualify3( ptwXYPoints *d1, ptwXYPoints *d2, int64_t i1, int64_t i2 );
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
    for( i = 0; i < 10; i++ ) {
        if( ( status = ptwXY_setValueAtX( XY, 0.2 * i - .5, 0.7 + i + .1 ) ) != nfu_Okay )
                nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX 1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkMutualify( XY );
    ptwXY_neg( XY );
    errCount += checkMutualify( XY );

    ptwXY_free( XY );

    exit( errCount );
}
/*
************************************************************
*/
static int checkMutualify( ptwXYPoints *data ) {

    int errCount = 0;
    int64_t n = data->length - 1;

    errCount += checkMutualify2( data, 2, n );
    errCount += checkMutualify2( data, 0, n - 3 );
    return( errCount );
}
/*
************************************************************
*/
static int checkMutualify2( ptwXYPoints *data, int64_t i1, int64_t i2 ) {

    int errCount = 0;
    ptwXYPoints *clone, *sliced;
    nfu_status status;

    if( ( clone = ptwXY_clone( data, &status ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: cloning data, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( sliced = ptwXY_slice( data, i1, i2, 0, &status ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: slicing data, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    errCount += checkMutualify3( clone, sliced, i1, i2 );

    if( ( clone = ptwXY_clone( data, &status ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: cloning data, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( sliced = ptwXY_slice( data, i1, i2, 0, &status ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: slicing data, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    errCount += checkMutualify3( sliced, clone, i1, i2 );

    return( errCount );
}
/*
************************************************************
*/
static int checkMutualify3( ptwXYPoints *d1, ptwXYPoints *d2, int64_t i1, int64_t i2 ) {

    int errCount = 0, positiveXOnly = 1;
    double lowerEps = 1e-6, upperEps = 1e-6;
    nfu_status status;

    if( verbose ) {
        printf( "# i1 = %d\n", (int) i1 );
        printf( "# i2 = %d\n", (int) i2 );
        printf( "# lowerEps = %.14e\n", lowerEps );
        printf( "# upperEps = %.14e\n", upperEps );
        printf( "# positiveXOnly = %d\n", positiveXOnly );
    }
    printIfVerbose( d1 );
    printIfVerbose( d2 );

    if( ( status = ptwXY_mutualifyDomains( d1, lowerEps, upperEps, positiveXOnly, d2, lowerEps, upperEps, positiveXOnly ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: ptwXY_mutualifyDomains, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    printIfVerbose( d1 );
    printIfVerbose( d2 );

    if( ( status = ptwXY_areDomainsMutual( d1, d2 ) ) != nfu_Okay ) {
        errCount++;
        nfu_printMsg( "ERROR %s: ptwXY_MutualifyDomains, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }

    ptwXY_free( d1 );
    ptwXY_free( d2 );

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
