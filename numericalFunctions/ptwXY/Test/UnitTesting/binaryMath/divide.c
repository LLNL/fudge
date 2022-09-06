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
#include <ptwXY_utilities.h>

struct XYData {
    int points;
    double *XYs;
};

static int verbose = 0;
static double biSectionMax = 5., accuracy = 1e-3;
static char *fmtXY = "%18.11e %18.11e\n";

static int divide( statusMessageReporting *smr, int n1, double *XY1, int n2, double *XY2 );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0, i, j, nDatas;
    double XY1[6] = { 2., 0., 2.2, 3., 2.6, 0. }, XY2[4] = { 2., 0., 2.6, 3. }, XY3[4] = { 2., 2., 2.6, 0. }, XY4[4] = { 2., 2., 2.6, 1. };
    double XY5[10] = { 2., 0., 2.2, 3., 2.3, 0., 2.5, 0., 2.6, 0. }, XY6[10] = { 2., 0., 2.2, 0., 2.3, 0., 2.5, 3., 2.6, 0. };
    double XY7[4] = { 2., 1.0, 2.6, -1.0 }, XY8[6] = { 2., 1., 2.2, 0., 2.6, 2. };
    struct XYData XYDatas[] = { { sizeof( XY1 ) / ( 2 * sizeof( double ) ), XY1 }, { sizeof( XY2 ) / ( 2 * sizeof( double ) ), XY2 },
                                { sizeof( XY3 ) / ( 2 * sizeof( double ) ), XY3 }, { sizeof( XY4 ) / ( 2 * sizeof( double ) ), XY4 },
                                { sizeof( XY5 ) / ( 2 * sizeof( double ) ), XY5 }, { sizeof( XY6 ) / ( 2 * sizeof( double ) ), XY6 },
                                { sizeof( XY7 ) / ( 2 * sizeof( double ) ), XY7 }, { sizeof( XY8 ) / ( 2 * sizeof( double ) ), XY8 } };

    nDatas = sizeof( XYDatas ) / sizeof( struct XYData );
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

    for( i = 0; i < nDatas; i++ ) {
        for( j = 0; j < nDatas; j++ ) {
            errCount += divide( &smr, XYDatas[i].points, XYDatas[i].XYs, XYDatas[j].points, XYDatas[j].XYs );
        }
    }

    exit( errCount );
}
/*
************************************************************
*/
static int divide( statusMessageReporting *smr, int n1, double *XY1, int n2, double *XY2 ) {

    int errCount = 0;
    ptwXYPoints *ptwXY1, *ptwXY2, *division;

    if( ( ptwXY1 = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, biSectionMax, accuracy, 10, 10, n1, XY1, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( smr, "Via." );
    if( ( ptwXY2 = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, biSectionMax, accuracy, 10, 10, n2, XY2, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( smr, "Via." );

    if( verbose ) {
        printf( "# biSectionMax = %.12e\n", biSectionMax );
        printf( "# accuracy = %.12e\n", accuracy );
        printf( "# length = %d\n", (int) ptwXY_length( smr, ptwXY1 ) );
        ptwXY_simpleWrite( ptwXY1, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( smr, ptwXY2 ) );
        ptwXY_simpleWrite( ptwXY2, stdout, fmtXY );
        printf( "\n\n" );
    }

    if( ( division = ptwXY_div_ptwXY( smr, ptwXY1, ptwXY2, 1 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );

    if( verbose ) {
        printf( "# length = %d\n", (int) ptwXY_length( smr, division ) );
        ptwXY_simpleWrite( division, stdout, fmtXY );
        printf( "\n\n" );
    }

    ptwXY_free( ptwXY1 );
    ptwXY_free( ptwXY2 );
    ptwXY_free( division );

    return( errCount );
}
