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
static char *fmtXY = "%19.12e %19.12e\n";

static int checkIntegration( ptwXYPoints *data, double xMin, double xMax, double expectedSum );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;
    nfu_status status;
    ptwXYPoints *XY, *expXY, *mulXY;
    double x, XYs[8] = { 2., 2., 4., 4., 6., 2., 8., 6. };

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

    if( ( XY = ptwXY_create( ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 10, 10, 4, XYs, &status, 0 ) ) == NULL ) 
            nfu_printErrorMsg( "ERROR %s: XY new, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    errCount += checkIntegration( XY, ptwXY_domainMin( XY ), ptwXY_domainMax( XY ), 20. );
    errCount += checkIntegration( XY, 3., ptwXY_domainMax( XY ), 17.5 );
    errCount += checkIntegration( XY, 5., ptwXY_domainMax( XY ), 10.5 );
    errCount += checkIntegration( XY, 7., ptwXY_domainMax( XY ), 5. );
    errCount += checkIntegration( XY, ptwXY_domainMin( XY ), 7., 15. );
    errCount += checkIntegration( XY, ptwXY_domainMin( XY ), 5., 9.5 );
    errCount += checkIntegration( XY, ptwXY_domainMin( XY ), 3., 2.5 );
    errCount += checkIntegration( XY, 3, 7., 12.5 );

    ptwXY_clear( XY );
    for( i = 0; i < 501; i++ ) {
        x = i * M_PI / 50;
        if( ( status = ptwXY_setValueAtX( XY, x, sin( x ) ) ) != nfu_Okay )
                nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX 3, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }

    XYs[0] = ptwXY_domainMin( XY );
    XYs[1] = 0.;
    XYs[2] = ptwXY_domainMax( XY );
    XYs[3] = 1.;
    if( ( expXY = ptwXY_create( ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 100, 10, 2, XYs, &status, 0 ) ) == NULL )
            nfu_printErrorMsg( "ERROR %s: XYs create, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    printIfVerbose( expXY );
    if( ( status = ptwXY_exp( expXY, 1. ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: ptwXY_exp, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    printIfVerbose( expXY );
    if( ( mulXY = ptwXY_mul_ptwXY( XY, expXY, &status ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: ptwXY_mul_ptwXY, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    errCount += checkIntegration( mulXY, ptwXY_domainMin( mulXY ), ptwXY_domainMax( mulXY ), -1.71786 );
    errCount += checkIntegration( mulXY, ptwXY_domainMin( mulXY ) - 100, ptwXY_domainMax( mulXY ) + 100, -1.71786 );

    ptwXY_free( XY );
    ptwXY_free( expXY );
    ptwXY_free( mulXY );

    exit( errCount );
}
/*
************************************************************
*/
static int checkIntegration( ptwXYPoints *data, double xMin, double xMax, double expectedSum ) {

    int errCount = 0;
    double sum, invSum;
    nfu_status status;
    ptwXYPoints *normed;

    if( verbose ) {
        printf( "# xMin = %.12e\n", xMin );
        printf( "# xMax = %.12e\n", xMax );
    }
    printIfVerbose( data );
    invSum = ptwXY_integrate( data, xMax, xMin, &status );
    if( status != nfu_Okay ) nfu_printErrorMsg( "ERROR %s: data inv. integration, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    sum = ptwXY_integrate( data, xMin, xMax, &status );
    if( status != nfu_Okay ) nfu_printErrorMsg( "ERROR %s: data integration, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( fabs( sum - expectedSum ) > 1e-6 * ( fabs( sum ) + fabs( expectedSum ) ) ) {
        nfu_printMsg( "ERROR %s: sum = %.8e != expectedSum = %.8e, sum - expectedSum = %e", __FILE__, sum, expectedSum, sum - expectedSum );
        errCount += 1;
    }
    if( fabs( sum + invSum ) > 1e-12 * ( fabs( sum ) + fabs( invSum ) ) ) {
        nfu_printMsg( "ERROR %s: sum + invSum != 0, sum = %e  invSum = %e   sum + invSum = %e", __FILE__, sum, invSum, sum + invSum );
        errCount += 1;
    }
    if( verbose ) {
        printf( "# sum = %.12e  invSum = %.12e, dSum = %.12e\n", sum, invSum, sum + invSum );
    }

    if( ( normed = ptwXY_clone( data, &status ) ) == NULL )
        nfu_printMsg( "ERROR %s: cloning, %d, %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_normalize( normed ) ) != nfu_Okay )
        nfu_printMsg( "ERROR %s: norm, %d, %s", __FILE__, status, nfu_statusMessage( status ) );
    sum = ptwXY_integrateDomain( normed, &status );
    if( status != nfu_Okay ) nfu_printErrorMsg( "ERROR %s: domain integration, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    printIfVerbose( normed );
    if( fabs( 1. - sum ) > 1e-14 ) {
        nfu_printMsg( "ERROR %s: norm sum = %.14e != 1", __FILE__, sum );
        errCount += 1;
    }
    ptwXY_free( normed );

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
