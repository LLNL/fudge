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
#include <time.h>
#include <float.h>

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>


static int verbose = 0, timing = 0;
static char *fmtXY = "%25.17e %25.17e\n";
static FILE *infoF;
static double epsilon = 1e-10;

static int checkThinning( statusMessageReporting *smr, const char * const label, double epsilon, int n1, double *input, int n2, double *output );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0;
    statusMessageReporting smr;
#include "thinDomain.h"

    smr_initialize( &smr, smr_status_Ok );

    infoF = stdout;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else if( strcmp( "-t", argv[iarg] ) == 0 ) {
            timing = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    checkThinning( &smr, "test 1.1", epsilon,      n_test1_input, test1_input, n_test1_1output, test1_1output );
    checkThinning( &smr, "test 1.2", epsilon / 3., n_test1_input, test1_input, n_test1_2output, test1_2output );
    checkThinning( &smr, "test 1.3", 3 * epsilon,  n_test1_input, test1_input, n_test1_3output, test1_3output );
/*
*   Test for point close to x[0] but far from x[2].
*/
    checkThinning( &smr, "test 2.1", epsilon,      n_test2_input, test2_input, n_test2_1output, test2_1output );
    checkThinning( &smr, "test 2.2", epsilon / 3., n_test2_input, test2_input, n_test2_2output, test2_2output );
    checkThinning( &smr, "test 2.3", 3 * epsilon,  n_test2_input, test2_input, n_test2_3output, test2_3output );
/*
*   Test for point close to x[2] but far from x[0].
*/
    checkThinning( &smr, "test 3.1", epsilon,      n_test3_input, test3_input, n_test3_1output, test3_1output );
    checkThinning( &smr, "test 3.2", epsilon / 3., n_test3_input, test3_input, n_test3_2output, test3_2output );
    checkThinning( &smr, "test 3.3", 3 * epsilon,  n_test3_input, test3_input, n_test3_3output, test3_3output );

    exit( errCount );
}
/*
************************************************************
*/
static int checkThinning( statusMessageReporting *smr, const char * const label, double epsilon, int n1, double *input, int n2, double *output ) {

    int errCount = 0;
    ptwXYPoint *points1, *points2;
    ptwXYPoints *inputXY, *outputXY, *thinned;
    clock_t time0;
    int64_t i1;
    double x1, x2;

    if( ( inputXY = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, 10., 1e-3, n1, 10, n1, input, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( smr, "Via." );
    if( ( outputXY = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, 10., 1e-3, n2, 10, n2, output, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( smr, "Via." );

    if( verbose ) fprintf( infoF, "# label = %s\n# epsilon = %e\n", label, epsilon );
    printIfVerbose( inputXY );
    if( strcmp( "test 3.3", label ) ) {
        
    }
    printIfVerbose( outputXY );
    time0 = clock( );
    if( ( thinned = ptwXY_thinDomain( smr, inputXY, epsilon ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( timing ) printf( "# %s: time = %.4f sec", label, ( clock( ) - time0 ) / ( (double) CLOCKS_PER_SEC ) );
    if( timing ) printf( "  length = %d\n", (int) thinned->length );
    printIfVerbose( thinned );

    points1 = thinned->points;
    x1 = points1->x;
    for( i1 = 1, ++points1; i1 < thinned->length; ++i1, ++points1 ) {
        x2 = points1->x;
        if( ( x2 - x1 ) < ( 0.5 * epsilon * ( fabs( x1 ) + fabs( x2 ) ) ) ) break;
        x1 = x2;
    }
    if( thinned->length != outputXY->length ) {
        fprintf( infoF, "ERROR %s: %s compare, length not equal: %d vs. %d\n", __FILE__, label, (int) thinned->length, (int) outputXY->length ); }
    else {
        points1 = thinned->points;
        points2 = outputXY->points;
        for( i1 = 0; i1 < thinned->length; ++i1, ++points1, ++points2 ) {
            x1 = points1->x;
            x2 = points2->x;
            if( fabs( x2 - x1 ) > ( 5 * DBL_EPSILON *( fabs( x1 ) + fabs( x2 ) ) ) ) {
                fprintf( infoF, "ERROR %s: %s compare, answer and results differ at index %d: %.17e vs. %.17e\n", 
                        __FILE__, label, (int) i1, x1, x2 );
                break;
            }
        }
    }

    ptwXY_free( thinned );
    ptwXY_free( inputXY );
    ptwXY_free( outputXY );
    return( errCount );
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
