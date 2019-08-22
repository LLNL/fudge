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
#include <time.h>

#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>


static int verbose = 0;
static char *fmtXY = "%25.15e %25.15e\n";
static FILE *infoF;

static void flatInterpolationToLinear( ptwXYPoints *XYs );
static void flatInterpolationToLinear2( ptwXYPoints *XYs, double epsm, double epsp );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0, i;
    ptwXYPoints *fineXYs, *coarseXYs;
    nfu_status status;
    double *fineYs, *coarseYs, accuracy = 1e-3;
    double fineXs[] = { -2.0, -1.5, -1.4, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 2.2, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
    double coarseXs[] = { -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0 };
    int nFineXs= sizeof( fineXs ) / ( sizeof( double ) ), nCoarseXs = sizeof( coarseXs ) / ( sizeof( double ) );

    infoF = stdout;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) fprintf( stderr, "%s\n", __FILE__ );

    if( verbose ) fprintf( infoF, "# accuracy = %e\n", accuracy );

    if( ( fineYs = nfu_malloc( nFineXs * sizeof( double ) ) ) == NULL ) nfu_printErrorMsg( "ERROR %s: nfu_malloc-ing fineXYs", __FILE__ );
    for( i = 0; i < nFineXs; i++ ) fineYs[i] = 2 * i + 3;
    fineYs[nFineXs-1] = fineYs[nFineXs-2];
    if( ( fineXYs = ptwXY_new( ptwXY_interpolationFlat, NULL, 0., accuracy, 10, 10, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: dataXY creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_setXYDataFromXsAndYs( fineXYs, nFineXs, fineXs, fineYs ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: ptwXY_setXYDataFromXsAndYs status = %d for fineXYs: %s", __FILE__, status, nfu_statusMessage( status ) );
    flatInterpolationToLinear( fineXYs );

    if( ( coarseYs = nfu_malloc( nCoarseXs * sizeof( double ) ) ) == NULL ) nfu_printErrorMsg( "ERROR %s: nfu_malloc-ing nCoarseXs", __FILE__ );
    for( i = 0; i < nCoarseXs ; i++ ) coarseYs[i] = -i + 10;
    coarseYs[nCoarseXs-1] = coarseYs[nCoarseXs-2];
    if( ( coarseXYs = ptwXY_new( ptwXY_interpolationFlat, NULL, 0., accuracy, 10, 10, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: coarseXYs creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_setXYDataFromXsAndYs( coarseXYs, nCoarseXs, coarseXs, coarseYs ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: ptwXY_setXYDataFromXsAndYs status = %d for coarseXYs: %s", __FILE__, status, nfu_statusMessage( status ) );

    free( fineYs );
    free( coarseYs );
    ptwXY_free( fineXYs );
    ptwXY_free( coarseXYs );
    exit( errCount );
}
/*
************************************************************
*/
static void flatInterpolationToLinear( ptwXYPoints *XYs ) {

    double dx = 1e-2;

    flatInterpolationToLinear2( XYs, 0., dx );
    flatInterpolationToLinear2( XYs, dx, 0. );
    flatInterpolationToLinear2( XYs, dx, dx );
}
/*
************************************************************
*/
static void flatInterpolationToLinear2( ptwXYPoints *XYs, double epsm, double epsp ) {

    nfu_status status;
    ptwXYPoints *linearXYs;

    if( verbose ) {
        fprintf( infoF, "# epsm = %e\n", epsm );
        fprintf( infoF, "# epsp = %e\n", epsp );
    }
    printIfVerbose( XYs );
    if( ( linearXYs = ptwXY_flatInterpolationToLinear( XYs, epsm, epsp, &status ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: ptwXY_flatInterpolationToLinear: status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    printIfVerbose( linearXYs );
    ptwXY_free( linearXYs );
}
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    fprintf( infoF, "# length = %d\n", (int) ptwXY_length( data ) );
    ptwXY_simpleWrite( data, infoF, fmtXY );
    fprintf( infoF, "\n\n" );
}
