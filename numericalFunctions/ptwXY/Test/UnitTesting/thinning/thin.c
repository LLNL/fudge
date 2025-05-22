/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>


static int verbose = 0, timing = 0;
static char *fmtXY = "%25.17e %25.17e\n";
static FILE *infoF;

static int checkThinning( statusMessageReporting *smr, const char * const label, ptwXYPoints *ptwXY1, double accuracy, int deleteXYs );
static void compareDoubles( double d1, double d2, double eps, const char * const label );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0, i, n;
    ptwXYPoints *dataXY;
    double x, xMax;
    double zeros[] = { -1.5, 0., -1.49999999, 0., 0.5, 0., 2., 0., 2.2, 1, 2.20000001, 0 };
    double triangle[] = { -1.0, 0., 0., 1., 1.0, 0. };
    int nZeros = sizeof( zeros ) / ( 2 * sizeof( double ) ), nTriangles = sizeof( triangle ) / ( 2 * sizeof( double ) );
    statusMessageReporting smr;

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

    if( ( dataXY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., 1e-3, 10, 10, nZeros, zeros, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    checkThinning( &smr, "zeros", dataXY, 1e-3, 1 );

    if( ( dataXY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., 1e-3, 10, 10, nTriangles, triangle, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    checkThinning( &smr, "triangle", dataXY, 1e-3, 0 );
    if( ptwXY_thicken( &smr, dataXY, 100, 1e-5, 1. + 1e-5 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    checkThinning( &smr, "triangle thickened", dataXY, 1e-3, 1 );

    n = 1001;
    if( ( dataXY = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 5, 1e-5, n, 40, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i = 0; i < n; i++ ) {
        x = i * M_PI / ( n - 1 );
        dataXY->points[i].x = x;
        dataXY->points[i].y = sin( x );
    }
    dataXY->length = n;
    checkThinning( &smr, "sin", dataXY, 1e-3, 1 );

    n = 100001;
    if( ( dataXY = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 5, 1e-5, n, 40, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i = 0; i < n; i++ ) {
        x = i * M_PI / ( n - 1 );
        dataXY->points[i].x = x;
        dataXY->points[i].y = sin( x );
    }
    dataXY->length = n;
    checkThinning( &smr, "sin with many points", dataXY, 1e-3, 1 );

    n = 1000001;
    if( ( dataXY = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 5, 1e-5, n, 40, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    xMax = 25 * M_PI;
    for( i = 0; i < n; i++ ) {
        x = i  * xMax / ( n - 1 );
        dataXY->points[i].x = x;
        dataXY->points[i].y = exp( 2. * ( x - 0.2 * xMax ) / xMax ) * sin( x );
    }
    dataXY->length = n;
    checkThinning( &smr, "exp * sin with many points", dataXY, 1e-2, 1 );

    exit( errCount );
}
/*
************************************************************
*/
static int checkThinning( statusMessageReporting *smr, const char * const label, ptwXYPoints *ptwXY1, double accuracy, int deleteXYs ) {

    int errCount = 0;
    ptwXYPoints *thinned;
    clock_t time0;
    int64_t i;
    double x, y;

    if( verbose ) fprintf( infoF, "# label = %s\n# accuracy = %e\n", label, accuracy );
    printIfVerbose( ptwXY1 );
    time0 = clock( );
    if( ( thinned = ptwXY_thin( smr, ptwXY1, accuracy ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( timing ) printf( "# %s: time = %.4f sec", label, ( clock( ) - time0 ) / ( (double) CLOCKS_PER_SEC ) );
    if( timing ) printf( "  length = %d\n", (int) thinned->length );
    printIfVerbose( thinned );
    for( i = 0; i < ptwXY1->length; i++ ) {     /* This logic requires that ptwXY_thin coalesed ptwXY1. */
        x = ptwXY1->points[i].x;
        if( ptwXY_getValueAtX( smr, thinned, x, &y ) != nfu_Okay )
            nfut_printSMRErrorExit2p( smr, "Via." );
        compareDoubles( y, ptwXY1->points[i].y, accuracy, label );
    }
    ptwXY_free( thinned );

    if( deleteXYs ) ptwXY_free( ptwXY1 );
    return( errCount );
}
/*
************************************************************
*/
static void compareDoubles( double d1, double d2, double eps, const char * const label ) {

    double s, d, r;

    s = 0.5 * ( fabs( d1 ) + fabs( d2 ) );
    d = d2 - d1;
    r = d;
    if( s != 0 ) r /= s;
    if( fabs( r ) > eps ) fprintf( infoF, "ERROR %s: %s compare, %e %e %e %e %e\n", __FILE__, label, d1, d2, s, d, r );
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
