/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_utilities.h>

typedef nfu_status (* integrator)( statusMessageReporting *smr, ptwXYPoints *, double, double, double * );

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int checkIntegration( statusMessageReporting *smr, ptwXYPoints *data, double xMin, double xMax, double expectedIntegral, integrator integrator_ );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errCount = 0, nXYs;
    ptwXYPoints *XY;
    double domainMin, domainMax;
    double sum, XYs[] = { 1., 1., 2., 2., 3., 4., 5., 1., 6., -1. };
    double sumRegion1x = 2. * ( 4. * sqrt( 2. ) - 1. ) / 5., sumRegion2x = 8. * ( 6. * sqrt( 3. ) - sqrt( 2. ) ) / 15., 
        sumRegion3x = 40. * sqrt( 5. ) / 3. - 58. * sqrt( 3. ) / 5., sumRegion4x = 76. * sqrt( 6. ) / 5. - 50. * sqrt( 5. ) / 3.;
    double sumLowerRegion1x = ( 9. * sqrt( 6. ) - 8. ) / 20., sumLowerRegion2x = ( 25. * sqrt( 5. ) - 16. ) / 15. / sqrt( 2. ), 
        sumLowerRegion3x = 392. / 15. - 58. * sqrt( 3. ) / 5., sumLowerRegion4x = 121. * sqrt( 22. ) / 15. - 50. * sqrt( 5. ) / 3.;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

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

    if( ( XY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 10, 10, nXYs, XYs, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );

    errCount += checkIntegration( &smr, XY, XYs[0], XYs[2], sumRegion1x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( &smr, XY, XYs[2], XYs[4], sumRegion2x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( &smr, XY, XYs[4], XYs[6], sumRegion3x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( &smr, XY, XYs[6], XYs[8], sumRegion4x, ptwXY_integrateWithWeight_sqrt_x );

    if( ptwXY_domainMin( &smr, XY, &domainMin ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_domainMax( &smr, XY, &domainMax ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkIntegration( &smr, XY, domainMin, domainMax, sumRegion1x + sumRegion2x + sumRegion3x + sumRegion4x, 
        ptwXY_integrateWithWeight_sqrt_x );

    errCount += checkIntegration( &smr, XY, domainMin - 10, domainMax + 10, sumRegion1x + sumRegion2x + sumRegion3x + sumRegion4x, 
        ptwXY_integrateWithWeight_sqrt_x );

    errCount += checkIntegration( &smr, XY, XYs[0], 0.5 * ( XYs[0] + XYs[2] ), sumLowerRegion1x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( &smr, XY, XYs[2], 0.5 * ( XYs[2] + XYs[4] ), sumLowerRegion2x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( &smr, XY, XYs[4], 0.5 * ( XYs[4] + XYs[6] ), sumLowerRegion3x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( &smr, XY, XYs[6], 0.5 * ( XYs[6] + XYs[8] ), sumLowerRegion4x, ptwXY_integrateWithWeight_sqrt_x );

    errCount += checkIntegration( &smr, XY, XYs[0], 0.5 * ( XYs[6] + XYs[8] ), 
        sumRegion1x + sumRegion2x + sumRegion3x + sumLowerRegion4x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( &smr, XY, 0.5 * ( XYs[0] + XYs[2] ), 0.5 * ( XYs[6] + XYs[8] ), 
        ( sumRegion1x - sumLowerRegion1x ) + sumRegion2x + sumRegion3x + sumLowerRegion4x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( &smr, XY, 0.5 * ( XYs[2] + XYs[4] ), 0.5 * ( XYs[6] + XYs[8] ), 
        ( sumRegion2x - sumLowerRegion2x ) + sumRegion3x + sumLowerRegion4x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( &smr, XY, 0.5 * ( XYs[2] + XYs[4] ), 0.5 * ( XYs[4] + XYs[6] ), 
        ( sumRegion2x - sumLowerRegion2x ) + sumLowerRegion3x, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( &smr, XY, 0.5 * ( XYs[0] + XYs[2] ), 0.5 * ( XYs[4] + XYs[6] ), 
        ( sumRegion1x - sumLowerRegion1x ) + sumRegion2x + sumLowerRegion3x, ptwXY_integrateWithWeight_sqrt_x );

    if( ptwXY_integrateDomainWithWeight_sqrt_x( &smr, XY, &sum ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( fabs( sum - ( sumRegion1x + sumRegion2x + sumRegion3x + sumRegion4x ) ) > 1e-12 * sum )
        nfu_printErrorMsg( "ERROR %s: ptwXY_integrateDomainWithWeight_sqrt_x = %.12e != sum = %.12e", __FILE__, sum,
            sumRegion1x + sumRegion2x + sumRegion3x + sumRegion4x );

    ptwXY_free( XY );

/*
*   Testing for flat interpolations.
*/
    if( ( XY = ptwXY_create( &smr, ptwXY_interpolationFlat, NULL, 4, 1.e-3, 10, 10, nXYs, XYs, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );

    errCount += checkIntegration( &smr, XY, XYs[0], XYs[8], 22.678150766024306, ptwXY_integrateWithWeight_sqrt_x );
    errCount += checkIntegration( &smr, XY, 0.5 * ( XYs[0] + XYs[2] ), 0.5 * ( XYs[4] + XYs[6] ), 11.294767148502110, ptwXY_integrateWithWeight_sqrt_x );

    ptwXY_free( XY );

    exit( errCount );
}
/*
************************************************************
*/
static int checkIntegration( statusMessageReporting *smr, ptwXYPoints *data, double xMin, double xMax, double expectedIntegral, integrator integrator_ ) {

    int errCount = 0;
    double sum, invSum;

    if( verbose ) {
        printf( "# xMin = %.12e\n", xMin );
        printf( "# xMax = %.12e\n", xMax );
    }
    printIfVerbose( data );
    if( integrator_( smr, data, xMax, xMin, &invSum ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( integrator_( smr, data, xMin, xMax, &sum    ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
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
