/*
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <ptwXY_utilities.h>

static int verbose = 0;
char fmt[] = "%22.14e %22.14e\n";

#if 0
nfu_status chargedParticleGetValue( void *argList, double x, double *y, double x1, double y1, double x2, double y2 );
#endif
void printMsg( const char *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errCount = 0, nPoints;
    ptwXYPoints *pSparse, *pDense;
    double accuracy = 1e-2;
    double xys1[] = { /* Following is data from ENDF/B-VII.1 reaction H2(H2,n)He3. */
        1e2, 1.1231e-58,          1e3, 5.5988e-18,          1e4, 8.8216e-6,           2e4, 2.7734e-4,           3e4, 1.1747e-3, 
        4e4, 2.6785e-3,           5e4, 4.6105e-3,           6e4, 6.805e-3,            7e4, 9.1418e-3,           8e4, 0.011541, 
        9e4, 0.01395,             1e5, 0.016339,            1.5e5, 0.027495,          2e5, 0.037113,            2.5e5, 0.045358 };
    FILE *ff;
#include "chargedParticle.h"
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else {
            printMsg( "Error %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );
    
    nfu_setMemoryDebugMode( 0 );

    nPoints = sizeof( xys1 ) / sizeof( xys1[0] ) / 2;
    if( ( pSparse = ptwXY_create( &smr, ptwXY_interpolationOther, "charged-particle", 5, accuracy, 10, 10, nPoints, xys1, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    ff = fopen( "curve_sparse.dat", "w" );
    fprintf( ff, "# accuracy = %e\n", accuracy );
    fprintf( ff, "# length = %d\n", (int) pSparse->length );
    ptwXY_simpleWrite( pSparse, ff, fmt );
    fclose( ff );

    nPoints = sizeof( xys1_answer ) / sizeof( xys1_answer[0] ) / 2;
    if( ( pDense = ptwXY_create( &smr, ptwXY_interpolationLinLin, "charged-particle", 5, accuracy, 10, 10, nPoints, xys1_answer, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    ff = fopen( "curve_dense.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) pDense->length );
    ptwXY_simpleWrite( pDense, ff, fmt );
    fclose( ff );

    ptwXY_free( pSparse );
    ptwXY_free( pDense );

    exit( errCount ? EXIT_FAILURE : EXIT_SUCCESS );
}
#if 0
/*
****************************************************************
*/
nfu_status chargedParticleGetValue( void *argList, double x, double *y, double x1, double y1, double x2, double y2 ) {

    double A, B, T = *((double *) argList);

    B = log( x2 * y2 / ( x1 * y1 ) ) / ( 1. / sqrt( x1 - T ) - 1. / sqrt( x2 - T ) );
    A = x1 * y1 * exp( B / sqrt( x1 - T ) );
    *y = A * exp( - B / sqrt( x - T ) ) / x;
    return( nfu_Okay );
}
#endif
/*
****************************************************************
*/
void printMsg( const char *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    fprintf( stderr, "\n" );
    va_end( args );
    exit( EXIT_FAILURE );
}
