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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

    if( ( XY = ptwXY_new( ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 10, 10, &status, 0 ) ) == NULL ) 
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
