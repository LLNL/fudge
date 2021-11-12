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

#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>
#include <nfut_utilities.h>

static int verbose = 0;
static int Bins[] = { 2, 3, 4, 5, 7, 10, 20, 100 };
static int nBins = sizeof( Bins ) / sizeof( Bins[0] );

static int epb( statusMessageReporting *smr, ptwXY_interpolation interpolation, int nXYs, double *XYs );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0;
    double XYs1[] = { 0, 2.643700e-2, 0.4, 9.808100e-1, 0.8, 1.435800e+0, 1.2, 5.697500e-2, 1.6, 0 };
    double XYs2[] = { 0, 9.164300e-4, 0.4, 3.272800e-2, 0.8, 9.592600e-1, 1.2, 1.250000e+0, 1.6, 2.571200e-1, 2, 0 };
    statusMessageReporting smr;
    int nXYs1 = sizeof( XYs1 ) / sizeof( XYs1[0] ) / 2;
    int nXYs2 = sizeof( XYs2 ) / sizeof( XYs2[0] ) / 2;

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

    errCount += epb( &smr, ptwXY_interpolationLinLin, nXYs1, XYs1 );
    errCount += epb( &smr, ptwXY_interpolationFlat, nXYs1, XYs1 );
    errCount += epb( &smr, ptwXY_interpolationLinLin, nXYs2, XYs2 );
    errCount += epb( &smr, ptwXY_interpolationFlat, nXYs2, XYs2 );

    exit( errCount );
}
/*
************************************************************
*/
static int epb( statusMessageReporting *smr, ptwXY_interpolation interpolation, int nXYs, double *XYs ) {

    int i1, i2, errCount = 0;
    ptwXYPoints *ptwXY;
    ptwXPoints *ptwX;
    double integral;

    if( ( ptwXY = ptwXY_create( smr, interpolation, NULL, 4, 1.e-3, 10, 10, nXYs, XYs, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( smr, "Via." );

    if( ptwXY_integrateDomain( smr, ptwXY, &integral ) != nfu_Okay )
        nfut_printSMRErrorExit2p( smr, "Via." );

    nfu_printXYDataOnVerbosity( verbose, ptwXY );
    for( i1 = 0; i1 < nBins; ++i1 ) {
        if( verbose ) printf( "# numberOfBins = %d\n", Bins[i1] );
        if( ( ptwX = ptwXY_equalProbableBins( smr, ptwXY, Bins[i1] ) ) == NULL ) {
            errCount++; }
        else {
            nfu_printXDataOnVerbosity( verbose, ptwX );
            for( i2 = 0; i2 < Bins[i1]; ++i2 ) {
                double x1 = ptwX_getPointAtIndex_Unsafely( ptwX, i2 );
                double x2 = ptwX_getPointAtIndex_Unsafely( ptwX, i2 + 1 );
                double subIntegral;
                nfu_status status = ptwXY_integrate( smr, ptwXY, x1, x2, &subIntegral );

                if( status != nfu_Okay ) {
                    nfu_printSMRError2p( smr, "Via" );
                    break;
                }
                if( nfu_cmpDoubles( subIntegral, integral / Bins[i1], 1e-13 ) ) {
                    fprintf( stderr, "Boundary error: nXYs = %d  i1 = %d  number of bins = %d\n", nXYs, i1, Bins[i1] );
                    ++errCount;
                    break;
                }
            }
            ptwX_free( ptwX );
        }
    }

    ptwXY_free( ptwXY );
    return( errCount );
}
