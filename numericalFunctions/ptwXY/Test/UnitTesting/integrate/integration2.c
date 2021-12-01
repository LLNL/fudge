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
#include <ptwXY_utilities.h>

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int integrate( statusMessageReporting *smr, int n1, double *xys1, double xMin, double xMax );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0;
    double XYs1[8] = { 0.00000000e+00, 0.00000000e+00, 1.95312500e+04, 9.59428000e-08, 2.14843750e+05, 2.76211000e-07, 2.00000000e+07, 1.58192490e-12 };
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

    errCount += integrate( &smr, 4, XYs1, -1e-5, 1e+3 );
    errCount += integrate( &smr, 4, XYs1, 1e-5, 1e+3 );
    errCount += integrate( &smr, 4, XYs1, 1e-5, XYs1[2] );
    errCount += integrate( &smr, 4, XYs1, 1e-5, 0.5 * ( XYs1[2] + XYs1[4] ) );
    errCount += integrate( &smr, 4, XYs1, 0.5 * ( XYs1[2] + XYs1[4] ), 0.6 * ( XYs1[2] + XYs1[4] ) );
    errCount += integrate( &smr, 4, XYs1, 0.5 * ( XYs1[2] + XYs1[4] ), 0.6 * ( XYs1[4] + XYs1[6] ) );

    exit( errCount );
}
/*
************************************************************
*/
static int integrate( statusMessageReporting *smr, int n1, double *xys1, double xMin, double xMax ) {

    int errCount = 0;
    nfu_status status;
    double integral;
    ptwXYPoints *ptwXY1;

    if( ( ptwXY1 = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, 8, 1e-5, 10, 10, n1, xys1, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( smr, "Via." );

    if( verbose ) {
        printf( "# xMin = %.8e\n", xMin );
        printf( "# xMax = %.8e\n", xMax );
    }
    printIfVerbose( ptwXY1 );

    if( ( status = ptwXY_integrate( smr, ptwXY1, xMin, xMax, &integral ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( verbose ) printf( "# integral = %.8e\n\n\n", integral );
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
