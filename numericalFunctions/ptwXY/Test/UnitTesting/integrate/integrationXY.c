/*
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
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

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int checkIntegration( statusMessageReporting *smr, ptwXYPoints *data, double xMin, double xMax, double expectedSum );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;
    ptwXYPoints *XY, *expXY, *mulXY;
    double domainMin, domainMax;
    double x, XYs[8] = { 2., 2., 4., 4., 6., 2., 8., 6. };
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

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

    if( ( XY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 10, 10, 4, XYs, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_domainMin( &smr, XY, &domainMin ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_domainMax( &smr, XY, &domainMax ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkIntegration( &smr, XY, domainMin, domainMax, 20. );
    errCount += checkIntegration( &smr, XY, 3., domainMax, 17.5 );
    errCount += checkIntegration( &smr, XY, 5., domainMax, 10.5 );
    errCount += checkIntegration( &smr, XY, 7., domainMax, 5. );
    errCount += checkIntegration( &smr, XY, domainMin, 7., 15. );
    errCount += checkIntegration( &smr, XY, domainMin, 5., 9.5 );
    errCount += checkIntegration( &smr, XY, domainMin, 3., 2.5 );
    errCount += checkIntegration( &smr, XY, 3, 7., 12.5 );

    if( ptwXY_clear( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i = 0; i < 501; i++ ) {
        x = i * M_PI / 50;
        if( ptwXY_setValueAtX( &smr, XY, x, sin( x ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }

    if( ptwXY_domainMin( &smr, XY, &domainMin ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_domainMax( &smr, XY, &domainMax ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    XYs[0] = domainMin;
    XYs[1] = 0.;
    XYs[2] = domainMax;
    XYs[3] = 1.;
    if( ( expXY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 100, 10, 2, XYs, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    printIfVerbose( expXY );
    if( ptwXY_exp( &smr, expXY, 1. ) != nfu_Okay ) nfut_printSMRError2p( &smr, "Via." );
    printIfVerbose( expXY );
    if( ( mulXY = ptwXY_mul_ptwXY( &smr, XY, expXY ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_domainMin( &smr, mulXY, &domainMin ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_domainMax( &smr, mulXY, &domainMax ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkIntegration( &smr, mulXY, domainMin, domainMax, -1.71786 );
    errCount += checkIntegration( &smr, mulXY, domainMin - 100, domainMax + 100, -1.71786 );

    ptwXY_free( XY );
    ptwXY_free( expXY );
    ptwXY_free( mulXY );

    exit( errCount );
}
/*
************************************************************
*/
static int checkIntegration( statusMessageReporting *smr, ptwXYPoints *data, double xMin, double xMax, double expectedSum ) {

    int errCount = 0;
    double sum, invSum;
    nfu_status status;
    ptwXYPoints *normed;

    if( verbose ) {
        printf( "# xMin = %.12e\n", xMin );
        printf( "# xMax = %.12e\n", xMax );
    }
    printIfVerbose( data );
    if( ( status = ptwXY_integrate( smr, data, xMax, xMin, &invSum ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ( status = ptwXY_integrate( smr, data, xMin, xMax, &sum ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
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

    if( ( normed = ptwXY_clone( smr, data ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_normalize( smr, normed ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_integrateDomain( smr, normed, &sum ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
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
