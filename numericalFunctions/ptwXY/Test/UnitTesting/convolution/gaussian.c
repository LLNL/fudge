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

#include <ptwXY.h>
#include <nf_utilities.h>
#include <nfut_utilities.h>
#include <ptwXY_utilities.h>

static int verbose = 0;
static char *fmtXY = "%.12e %.12e\n";

static void printIfVerbose( const char * const label, ptwXYPoints *data, double accuracy, double offset, double sigma, double amplitude, 
    double xMin, double xMax, double dullEps );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0;
    double accuracy = 1e-5, area, offset, amplitude = 5., sigma = 10., xMin, xMax;
    ptwXYPoints *gaussian;
    nfu_status status;
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

    if( ( gaussian = ptwXY_createGaussianCenteredSigma1( &smr, accuracy ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    printIfVerbose( "ptwXY_createGaussianCenteredSigma1", gaussian, accuracy, 0., 1., 1., 0., 0., 0. );

    if( ( status = ptwXY_integrateDomain( &smr, gaussian, &area ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "# area = %.12e  area err = %e\n", area, area / 2.50662827463 - 1 );
    if( fabs( area / 2.50662827463 - 1 ) > 0.2 * accuracy )
        nfu_printErrorMsg( "ERROR %s: area err = %e, area = %.12e", __FILE__, area / 2.50662827463 - 1, area );
    ptwXY_free( gaussian );

    offset = 110.;
    xMin = offset - 66.;
    xMax = offset + 66.;
    if( ( gaussian = ptwXY_createGaussian( &smr, accuracy, offset, sigma, amplitude, xMin, xMax, 0. ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    printIfVerbose( "ptwXY_createGaussian", gaussian, accuracy, offset, sigma, amplitude, xMin, xMax, 0. );
    if( ( status = ptwXY_integrateDomain( &smr, gaussian, &area ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "# area = %.12e  area err = %e\n", area, area / ( sigma * amplitude * 2.50662827463 ) - 1 );
    if( fabs( area / ( sigma * amplitude * 2.50662827463 ) - 1 ) > 0.2 * accuracy )
        nfu_printErrorMsg( "ERROR %s: area err = %e, area = %.12e", __FILE__, area / ( sigma * amplitude * 2.50662827463 ) - 1, area );
    ptwXY_free( gaussian );

    return( 0 );
}
/*
************************************************************
*/
static void printIfVerbose( const char * const label, ptwXYPoints *data, double accuracy, double offset, double sigma, double amplitude, 
    double xMin, double xMax, double dullEps ) {

    if( !verbose ) return;
    printf( "# %s\n", label );
    printf( "# accuracy = %g\n", accuracy );
    printf( "# offset = %g\n", offset );
    printf( "# sigma = %g\n", sigma );
    printf( "# amplitude = %g\n", amplitude );
    printf( "# xMin = %.17e\n", xMin );
    printf( "# xMax = %.17e\n", xMax );
    printf( "# dullEps = %g\n", dullEps );
    printf( "# length = %d\n", (int) data->length );
    ptwXY_simpleWrite( data, stdout, fmtXY );
    printf( "\n\n" );
}
