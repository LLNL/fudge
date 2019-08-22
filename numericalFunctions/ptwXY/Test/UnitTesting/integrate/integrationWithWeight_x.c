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

typedef double (* integrator)( ptwXYPoints *, double, double, nfu_status * );

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int checkIntegration( ptwXYPoints *data, double xMin, double xMax, double expectedIntegral, integrator integrator_ );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0, nXYs;
    nfu_status status;
    ptwXYPoints *XY, *y_is_x;
    double sum, XYs[] = { 1., 1., 2., 2., 3., 4., 5., 1., 6., -1. };
    double sumRegion1x = 7. / 3., sumRegion2x = 23. / 3., sumRegion3x = 19., sumRegion4x = -1. / 6.;
    double sumLowerRegion1x = 19. / 24, sumLowerRegion2x = 17. / 6., sumLowerRegion3x = 45. / 4., sumLowerRegion4x = 31. / 24.;
    ptwXPoints *groupBoundaries, *groupIntegrals;

    nXYs = sizeof( XYs ) / ( 2 * sizeof( XYs[0] ) );

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

    if( ( XY = ptwXY_create( ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 10, 10, nXYs, XYs, &status, 0 ) ) == NULL ) 
            nfu_printErrorMsg( "ERROR %s: XY create, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    errCount += checkIntegration( XY, XYs[0], XYs[2], sumRegion1x, ptwXY_integrateWithWeight_x );
    errCount += checkIntegration( XY, XYs[2], XYs[4], sumRegion2x, ptwXY_integrateWithWeight_x );
    errCount += checkIntegration( XY, XYs[4], XYs[6], sumRegion3x, ptwXY_integrateWithWeight_x );
    errCount += checkIntegration( XY, XYs[6], XYs[8], sumRegion4x, ptwXY_integrateWithWeight_x );

    errCount += checkIntegration( XY, ptwXY_domainMin( XY ), ptwXY_domainMax( XY ), sumRegion1x + sumRegion2x + sumRegion3x + sumRegion4x, 
        ptwXY_integrateWithWeight_x );
    errCount += checkIntegration( XY, ptwXY_domainMin( XY ) - 10, ptwXY_domainMax( XY ) + 10, sumRegion1x + sumRegion2x + sumRegion3x + sumRegion4x, 
        ptwXY_integrateWithWeight_x );

    errCount += checkIntegration( XY, XYs[0], 0.5 * ( XYs[0] + XYs[2] ), sumLowerRegion1x, ptwXY_integrateWithWeight_x );
    errCount += checkIntegration( XY, XYs[2], 0.5 * ( XYs[2] + XYs[4] ), sumLowerRegion2x, ptwXY_integrateWithWeight_x );
    errCount += checkIntegration( XY, XYs[4], 0.5 * ( XYs[4] + XYs[6] ), sumLowerRegion3x, ptwXY_integrateWithWeight_x );
    errCount += checkIntegration( XY, XYs[6], 0.5 * ( XYs[6] + XYs[8] ), sumLowerRegion4x, ptwXY_integrateWithWeight_x );

    errCount += checkIntegration( XY, XYs[0], 0.5 * ( XYs[6] + XYs[8] ), 
        sumRegion1x + sumRegion2x + sumRegion3x + sumLowerRegion4x, ptwXY_integrateWithWeight_x );
    errCount += checkIntegration( XY, 0.5 * ( XYs[0] + XYs[2] ), 0.5 * ( XYs[6] + XYs[8] ), 
        ( sumRegion1x - sumLowerRegion1x ) + sumRegion2x + sumRegion3x + sumLowerRegion4x, ptwXY_integrateWithWeight_x );
    errCount += checkIntegration( XY, 0.5 * ( XYs[2] + XYs[4] ), 0.5 * ( XYs[6] + XYs[8] ), 
        ( sumRegion2x - sumLowerRegion2x ) + sumRegion3x + sumLowerRegion4x, ptwXY_integrateWithWeight_x );
    errCount += checkIntegration( XY, 0.5 * ( XYs[2] + XYs[4] ), 0.5 * ( XYs[4] + XYs[6] ), 
        ( sumRegion2x - sumLowerRegion2x ) + sumLowerRegion3x, ptwXY_integrateWithWeight_x );
    errCount += checkIntegration( XY, 0.5 * ( XYs[0] + XYs[2] ), 0.5 * ( XYs[4] + XYs[6] ), 
        ( sumRegion1x - sumLowerRegion1x ) + sumRegion2x + sumLowerRegion3x, ptwXY_integrateWithWeight_x );

    sum = ptwXY_integrateDomainWithWeight_x( XY, &status );
    if( status != nfu_Okay ) nfu_printErrorMsg( "ERROR %s: ptwXY_integrateWithWeight_x, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( fabs( sum - ( sumRegion1x + sumRegion2x + sumRegion3x + sumRegion4x ) ) > 1e-12 * sum )
        nfu_printErrorMsg( "ERROR %s: ptwXY_integrateDomainWithWeight_x = %.12e != sum = %.12e", __FILE__, sum,
            sumRegion1x + sumRegion2x + sumRegion3x + sumRegion4x );

    if( verbose ) printf( "# groupTwoFunctions testing\n" );
    if( ( groupBoundaries = ptwX_new( nXYs, &status ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: x ptwX_new, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    for( i = 0; i < nXYs; i++ ) {
        if( ( status = ptwX_setPointAtIndex( groupBoundaries, i, XYs[2 * i] ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: x ptwX_setPointAtIndex, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }

    if( ( y_is_x = ptwXY_new( ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 10, 10, &status, 0 ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: x ptwXY_new, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_setValueAtX( y_is_x, XYs[0], XYs[0] ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_setValueAtX( y_is_x, XYs[8], XYs[8] ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( groupIntegrals = ptwXY_groupTwoFunctions( XY, y_is_x, groupBoundaries, ptwXY_group_normType_none, NULL, &status ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: ptwXY_groupTwoFunctions, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    printIfVerbose( y_is_x );

    ptwXY_free( y_is_x );
    ptwX_free( groupBoundaries );
    ptwX_free( groupIntegrals );

    ptwXY_free( XY );

/*
*   Testing for flat interpolations.
*/
    if( ( XY = ptwXY_create( ptwXY_interpolationFlat, NULL, 4, 1.e-3, 10, 10, nXYs, XYs, &status, 0 ) ) == NULL ) 
            nfu_printErrorMsg( "ERROR %s: XY create, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    errCount += checkIntegration( XY, XYs[0], XYs[8], 44., ptwXY_integrateWithWeight_x );
    errCount += checkIntegration( XY, 0.5 * ( XYs[0] + XYs[2] ), 0.5 * ( XYs[4] + XYs[6] ), 39.75 / 2., ptwXY_integrateWithWeight_x );

    ptwXY_free( XY );

    exit( errCount );
}
/*
************************************************************
*/
static int checkIntegration( ptwXYPoints *data, double xMin, double xMax, double expectedIntegral, integrator integrator_ ) {

    int errCount = 0;
    double sum, invSum;
    nfu_status status;

    if( verbose ) {
        printf( "# xMin = %.12e\n", xMin );
        printf( "# xMax = %.12e\n", xMax );
    }
    printIfVerbose( data );
    invSum = integrator_( data, xMax, xMin, &status );
    if( status != nfu_Okay ) nfu_printErrorMsg( "ERROR %s: data inv. integration, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    sum = integrator_( data, xMin, xMax, &status );
    if( status != nfu_Okay ) nfu_printErrorMsg( "ERROR %s: data integration, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( fabs( sum - expectedIntegral ) > 1e-12 * ( fabs( sum ) + fabs( expectedIntegral ) ) ) {
        nfu_printMsg( "ERROR %s: sum = %.8e != expectedIntegral = %.8e, sum - expectedIntegral = %e", __FILE__, sum, expectedIntegral, sum - expectedIntegral );
        errCount += 1;
    }
    if( fabs( sum + invSum ) > 1e-12 * ( fabs( sum ) + fabs( invSum ) ) ) {
        nfu_printMsg( "ERROR %s: sum + invSum != 0, sum = %e  invSum = %e   sum + invSum = %e", __FILE__, sum, invSum, sum + invSum );
        errCount += 1;
    }
    if( verbose ) {
        printf( "# sum = %.12e  invSum = %.12e, dSum = %.12e\n", sum, invSum, sum + invSum );
    }

    return( errCount );
}
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    printf( "# length = %d\n", (int) data->length );
    printf( "# interpolation = %d\n", (int) data->interpolation );
    ptwXY_simpleWrite( data, stdout, fmtXY );
    printf( "\n\n" );
}
