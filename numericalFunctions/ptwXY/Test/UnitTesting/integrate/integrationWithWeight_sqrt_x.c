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

typedef double (* integrator)( ptwXYPoints *, double, double, nfu_status * );

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int checkIntegration( ptwXYPoints *data, double xMin, double xMax, double expectedIntegral, integrator integrator_ );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errCount = 0, nXYs;
    nfu_status status;
    ptwXYPoints *XY;
    double sum, XYs[] = { 1., 1., 2., 2., 3., 4., 5., 1., 6., -1. };
    double sumRegion1x = 2. * ( 4. * sqrt( 2. ) - 1. ) / 5., sumRegion2x = 8. * ( 6. * sqrt( 3. ) - sqrt( 2. ) ) / 15., 
        sumRegion3x = 40. * sqrt( 5. ) / 3. - 58. * sqrt( 3. ) / 5., sumRegion4x = 76. * sqrt( 6. ) / 5. - 50. * sqrt( 5. ) / 3.;
    double sumLowerRegion1x = ( 9. * sqrt( 6. ) - 8. ) / 20., sumLowerRegion2x = ( 25. * sqrt( 5. ) - 16. ) / 15. / sqrt( 2. ), 
        sumLowerRegion3x = 392. / 15. - 58. * sqrt( 3. ) / 5., sumLowerRegion4x = 121. * sqrt( 22. ) / 15. - 50. * sqrt( 5. ) / 3.;

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

    if( ( XY = ptwXY_create( ptwXY_interpolationLinLin, 4, 1.e-3, 10, 10, nXYs, XYs, &status, 0 ) ) == NULL ) 
            nfu_printErrorMsg( "ERROR %s: XY create, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    errCount += checkIntegration( XY, XYs[0], XYs[2], sumRegion1x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( XY, XYs[2], XYs[4], sumRegion2x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( XY, XYs[4], XYs[6], sumRegion3x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( XY, XYs[6], XYs[8], sumRegion4x, ptwXY_integrateWithWeight_sqrt_x );

    errCount += checkIntegration( XY, ptwXY_getXMin( XY ), ptwXY_getXMax( XY ), sumRegion1x + sumRegion2x + sumRegion3x + sumRegion4x, 
        ptwXY_integrateWithWeight_sqrt_x );

    errCount += checkIntegration( XY, ptwXY_getXMin( XY ) - 10, ptwXY_getXMax( XY ) + 10, sumRegion1x + sumRegion2x + sumRegion3x + sumRegion4x, 
        ptwXY_integrateWithWeight_sqrt_x );

    errCount += checkIntegration( XY, XYs[0], 0.5 * ( XYs[0] + XYs[2] ), sumLowerRegion1x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( XY, XYs[2], 0.5 * ( XYs[2] + XYs[4] ), sumLowerRegion2x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( XY, XYs[4], 0.5 * ( XYs[4] + XYs[6] ), sumLowerRegion3x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( XY, XYs[6], 0.5 * ( XYs[6] + XYs[8] ), sumLowerRegion4x, ptwXY_integrateWithWeight_sqrt_x );

    errCount += checkIntegration( XY, XYs[0], 0.5 * ( XYs[6] + XYs[8] ), 
        sumRegion1x + sumRegion2x + sumRegion3x + sumLowerRegion4x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( XY, 0.5 * ( XYs[0] + XYs[2] ), 0.5 * ( XYs[6] + XYs[8] ), 
        ( sumRegion1x - sumLowerRegion1x ) + sumRegion2x + sumRegion3x + sumLowerRegion4x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( XY, 0.5 * ( XYs[2] + XYs[4] ), 0.5 * ( XYs[6] + XYs[8] ), 
        ( sumRegion2x - sumLowerRegion2x ) + sumRegion3x + sumLowerRegion4x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( XY, 0.5 * ( XYs[2] + XYs[4] ), 0.5 * ( XYs[4] + XYs[6] ), 
        ( sumRegion2x - sumLowerRegion2x ) + sumLowerRegion3x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( XY, 0.5 * ( XYs[0] + XYs[2] ), 0.5 * ( XYs[4] + XYs[6] ), 
        ( sumRegion1x - sumLowerRegion1x ) + sumRegion2x + sumLowerRegion3x, ptwXY_integrateWithWeight_sqrt_x );

    sum = ptwXY_integrateDomainWithWeight_sqrt_x( XY, &status );
    if( status != nfu_Okay ) nfu_printErrorMsg( "ERROR %s: ptwXY_integrateWithWeight_sqrt_x, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( fabs( sum - ( sumRegion1x + sumRegion2x + sumRegion3x + sumRegion4x ) ) > 1e-12 * sum )
        nfu_printErrorMsg( "ERROR %s: ptwXY_integrateDomainWithWeight_sqrt_x = %.12e != sum = %.12e", __FILE__, sum,
            sumRegion1x + sumRegion2x + sumRegion3x + sumRegion4x );

    ptwXY_free( XY );

/*
*   Testing for flat interpolations.
*/
    if( ( XY = ptwXY_create( ptwXY_interpolationFlat, 4, 1.e-3, 10, 10, nXYs, XYs, &status, 0 ) ) == NULL )
            nfu_printErrorMsg( "ERROR %s: XY create, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    errCount += checkIntegration( XY, XYs[0], XYs[8], 22.678150766024306, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( XY, 0.5 * ( XYs[0] + XYs[2] ), 0.5 * ( XYs[4] + XYs[6] ), 11.294767148502110, ptwXY_integrateWithWeight_sqrt_x );

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
