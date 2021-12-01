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
#include <time.h>

#include <ptwXY.h>
#include <nfut_utilities.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>


static int verbose = 0;
static char *fmtXY = "%25.15e %25.15e\n";
static FILE *infoF;

static void flatInterpolationToLinear( statusMessageReporting *smr, ptwXYPoints *XYs );
static void flatInterpolationToLinear2( statusMessageReporting *smr, ptwXYPoints *XYs, double epsm, double epsp );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0;
    double accuracy = 1e-3, XYData[] = { 0, 1.841839e-8, 5.7e6, 6.22331e-10, 5.9e6, 8.33984e-11, 6099999.939, 4.42672e-7, 6.1e6, 0 };
    int nXYs= sizeof( XYData ) / ( sizeof( double ) ) / 2;
    ptwXYPoints *XYs;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

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

    if( ( XYs = ptwXY_create( &smr, ptwXY_interpolationFlat, NULL, 0., accuracy, 10, 10, nXYs, XYData, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    flatInterpolationToLinear( &smr, XYs );

    ptwXY_free( XYs );
    exit( errCount );
}
/*
************************************************************
*/
static void flatInterpolationToLinear( statusMessageReporting *smr, ptwXYPoints *XYs ) {

    double dx = 1e-9;

    flatInterpolationToLinear2( smr, XYs, 0., dx );
    flatInterpolationToLinear2( smr, XYs, dx, 0. );
    flatInterpolationToLinear2( smr, XYs, dx, dx );
}
/*
************************************************************
*/
static void flatInterpolationToLinear2( statusMessageReporting *smr, ptwXYPoints *XYs, double epsm, double epsp ) {

    ptwXYPoints *linearXYs;

    if( verbose ) {
        fprintf( infoF, "# epsm = %e\n", epsm );
        fprintf( infoF, "# epsp = %e\n", epsp );
    }
    printIfVerbose( XYs );
    if( ( linearXYs = ptwXY_flatInterpolationToLinear( smr, XYs, epsm, epsp ) ) == NULL )
        nfut_printSMRError2p( smr, "Via." );
    printIfVerbose( linearXYs );
    ptwXY_free( linearXYs );
}
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    fprintf( infoF, "# length = %d\n", (int) ptwXY_length( NULL, data ) );
    ptwXY_simpleWrite( data, infoF, fmtXY );
    fprintf( infoF, "\n\n" );
}
