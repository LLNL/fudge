/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

/*
*   2011 July 11: A bug was discovered in the ptwXY_union in the fill section. This is a test of that the bug is fixed.
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

struct XYData {
    int points;
    double *XYs;
};

static int verbose = 0;
static double biSectionMax = 5., accuracy = 1e-3;
static char *fmtXY = "%18.11e %18.11e\n";

/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0;
    int64_t i;
    double x, y;
    double XY0[] = { 1.423095e+07, 0.0, 1.440000e+07, 1.745688e-05, 1.460000e+07, 2.285425e-04, 1.480000e+07, 7.085735e-04 };
    double XY1[] = { 1.423095e+07, 0.0, 1.423100e+07, 0.000000e+00, 1.440000e+07, 0.000000e+00, 1.460000e+07, 0.000000e+00, 1.480000e+07, 6.155894e-05 };
    int nXY0 = sizeof( XY0 ) / ( 2 * sizeof( double ) ), nXY1 = sizeof( XY1 ) / ( 2 * sizeof( double ) );
    ptwXYPoints *ptwXY0, *ptwXY1, *add1, *add2, *diff;
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
    if(echo ) printf( "%s\n", __FILE__ );

    if( ( ptwXY0 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, biSectionMax, accuracy, 10, 10, nXY0, XY0, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( ptwXY1 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, biSectionMax, accuracy, 10, 10, nXY1, XY1, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );

    if( ( add1 = ptwXY_add_ptwXY( &smr, ptwXY0, ptwXY1 ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( add2 = ptwXY_add_ptwXY( &smr, ptwXY1, ptwXY0 ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( diff = ptwXY_sub_ptwXY( &smr, add1, add2 ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );

    if( verbose ) {
        printf( "# biSectionMax = %.12e\n", biSectionMax );
        printf( "# accuracy = %.12e\n", accuracy );
        printf( "# length = %d\n", (int) ptwXY_length( &smr, ptwXY0 ) );
        ptwXY_simpleWrite( ptwXY0, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( &smr, ptwXY1 ) );
        ptwXY_simpleWrite( ptwXY1, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( &smr, add1 ) );
        ptwXY_simpleWrite( add1, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( &smr, add2 ) );
        ptwXY_simpleWrite( add2, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( &smr, diff ) );
        ptwXY_simpleWrite( diff, stdout, fmtXY );
        printf( "\n\n" );
    }

    for( i = 0; i < diff->length; i++ ) {
        if( ( status = ptwXY_getXYPairAtIndex( &smr, diff, i, &x, &y ) ) != nfu_Okay )
            nfut_printSMRErrorExit2p( &smr, "Via." );
        if( y != 0 ) {
            nfu_printMsg( "ERROR %s: value at index %d = %e is not 0.", __FILE__, (int) i, y );
            errCount++;
        }
    }

    ptwXY_free( ptwXY0 );
    ptwXY_free( ptwXY1 );
    ptwXY_free( add1 );
    ptwXY_free( add2 );
    ptwXY_free( diff );

    exit( errCount );
}
