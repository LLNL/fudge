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

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>

#define M_SQRT2PI ( 2. * M_SQRT2 / M_2_SQRTPI )

static int verbose = 0, mode = 0, doBruteForce = 0, timing = 0;
static char *fmtXY = "%25.17e %25.17e\n";
static FILE *infoF;

static int checkConvolution( statusMessageReporting *smr, const char * const label, ptwXYPoints *ptwXY1, 
        ptwXYPoints* ptwXY2, int deleteXYs, double accuracy, double area );
static void bruteForce( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, double accuracy );
static void compareDoubles( double d1, double d2, double eps, const char * const funcName );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0;
    ptwXYPoints *wedgeXY, *box1XY, *box2XY, *triangle1XY, *triangle2XY, *gaussian1XY, *gaussian2XY;
    double sigma, accuracy = 1e-3;
/*                      x1               x2       x3      x4       x5             x6 */
    double box1[] = { -1.5, 0., -1.49999999, 1., 0.5, 1., 2., 1., 2.2, 1, 2.20000001, 0 };
    int nBox1 = sizeof( box1 ) / ( 2 * sizeof( double ) );
/*                      x1               x2       x3              x4  x6 */
    double box2[] = { -0.5, 0., -0.499999999999, 1., 0.5, 1., 0.500000000001, 0. };
    int nBox2 = sizeof( box2 ) / ( 2 * sizeof( double ) );
    double wedge[] = { 0.0, 0., 1.0, 1., 1.0 + 1e-15, 0. };
    double triangle[] = { -1.0, 0., 0., 1., 1.0, 0. };
    double triangle1[] = { 1.0, 0., 2.0, 2., 3.0, 0. }, triangle2[] = { -2, 0., -0.0, 3., 2., 0. };
    double triangle3[] = { 1.0, 0., 1.5, 2., 3.0, 0. }, triangle4[] = { -5, 0., -0.5, 3., 2., 0. };
    int nTriangles = 3;
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
        else if( strcmp( "+1", argv[iarg] ) == 0 ) {
            mode = 1; }
        else if( strcmp( "-1", argv[iarg] ) == 0 ) {
            mode = -1; }
        else if( strcmp( "-b", argv[iarg] ) == 0 ) {
            doBruteForce = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    if( ( box1XY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, 10, 10, nBox1, box1, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( box2XY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, 10, 10, nBox2, box2, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    checkConvolution( &smr, "box 1 & 2", box1XY, box2XY, 1, accuracy, 3.7 );
exit( 0 );

    if( ( box1XY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, 10, 10, nBox1, box1, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( wedgeXY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, 10, 10, 3, wedge, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    checkConvolution( &smr, "box 1 + wedge", box1XY, wedgeXY, 1, accuracy, 1.85 );

    if( ( wedgeXY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, 10, 10, 3, wedge, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    checkConvolution( &smr, "wedge + wedge", wedgeXY, wedgeXY, 1, accuracy, 0.25 );

    if( ( triangle1XY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, 10, 10, nTriangles, triangle, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( box1XY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, 10, 10, nBox1, box1, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    checkConvolution( &smr, "box 1 + triangle", triangle1XY, box1XY, 1, accuracy, 3.70 );

    if( ( triangle1XY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, 10, 10, nTriangles, triangle, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( wedgeXY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, 10, 10, 3, wedge, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    checkConvolution( &smr, "triangle + wedge", triangle1XY, wedgeXY, 1, accuracy, 0.5 );

    if( ( triangle1XY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, 10, 10, nTriangles, triangle1, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( triangle2XY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, 10, 10, nTriangles, triangle2, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    checkConvolution( &smr, "triangles 1 & 2", triangle1XY, triangle2XY, 0, accuracy, 12. );
    if( ptwXY_thicken( &smr, triangle1XY, 100, 1e-5, 1. + 1e-5 ) != nfu_Okay )
        nfut_printSMRErrorExit2p( &smr, "Via." );
   if(  ptwXY_thicken( &smr, triangle2XY, 100, 1e-5, 1. + 1e-5 ) != nfu_Okay ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    checkConvolution( &smr, "triangles 1 & 2 thickened", triangle1XY, triangle2XY, 1, accuracy, 12. );

    if( ( triangle1XY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, 10, 10, nTriangles, triangle3, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( triangle2XY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, 10, 10, nTriangles, triangle4, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    checkConvolution( &smr, "triangles 3 & 4", triangle2XY, triangle1XY, 0, accuracy, 21.0 );
    if( ptwXY_thicken( &smr, triangle1XY, 100, 1e-5, 1. + 1e-5 ) != nfu_Okay ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_thicken( &smr, triangle2XY, 100, 1e-5, 1. + 1e-5 ) != nfu_Okay ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    checkConvolution( &smr, "triangles 3 & 4 thickened", triangle1XY, triangle2XY, 1, accuracy, 21.0 );

    sigma = 3.;
    accuracy = 2e-3;
    if( ( gaussian1XY = ptwXY_createGaussian( &smr, accuracy, 0., sigma, 1.0 / ( M_SQRT2PI * sigma ), -4. * sigma, 4. * sigma, 0. ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    sigma = 4.;
    if( ( gaussian2XY = ptwXY_createGaussian( &smr, accuracy, 0., sigma, 1.0 / ( M_SQRT2PI * sigma ), -4. * sigma, 4. * sigma, 0. ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    checkConvolution( &smr, "Gaussian", gaussian1XY, gaussian2XY, 1, accuracy, 1. );

    exit( errCount );
}
/*
************************************************************
*/
static int checkConvolution( statusMessageReporting *smr, const char * const label, ptwXYPoints *ptwXY1, 
        ptwXYPoints *ptwXY2, int deleteXYs, double accuracy, double area ) {

    int errCount = 0;
    ptwXYPoints *convolute;
    nfu_status status;
    double calculatedArea;
    clock_t time0;

    if( verbose ) {
        fprintf( infoF, "# label = %s\n", label );
        fprintf( infoF, "# mode = %d\n", mode );
        fprintf( infoF, "# accuracy = %e\n", accuracy );
    }
    printIfVerbose( ptwXY1 );
    printIfVerbose( ptwXY2 );
    time0 = clock( );
    convolute = ptwXY_convolution( smr, ptwXY1, ptwXY2, mode );
    if( timing ) printf( "# %s: time = %.4f sec", label, ( clock( ) - time0 ) / ( (double) CLOCKS_PER_SEC ) );
    if( convolute == NULL ) {
        if( timing ) printf( "\n" );
        errCount++;
        nfut_printSMRError2p( smr, "Via." ); }
    else {
        if( timing ) printf( ", length = %d\n", (int) convolute->length );
        printIfVerbose( convolute );
        if( ( status = ptwXY_integrateDomain( smr, convolute, &calculatedArea ) ) != nfu_Okay ) {
            errCount++;
            nfut_printSMRError2p( smr, "Via." );
        }
        if( verbose ) fprintf( infoF, "# Area = %.6g\n", calculatedArea );
        compareDoubles( area, calculatedArea, 1e-3, label );
        ptwXY_free( convolute );
    }

    if( doBruteForce ) bruteForce( smr, ptwXY1, ptwXY2, accuracy );

    if( deleteXYs ) {
        if( ptwXY2 != ptwXY1 ) ptwXY_free( ptwXY1 );
        ptwXY_free( ptwXY2 );
    }
    return( errCount );
}
/*
************************************************************
*/
static void bruteForce( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, double accuracy ) {

    int64_t i, j, k, n = 101;
    double s, ds, sMin, sMax, x11, x12, t, v;
    double x1Min, x1Max, x2Min, x2Max;
    nfu_status status;
    ptwXYPoints *f2, *clone, *mulXY, *bfXY, *d1, *d2;

    if( ptwXY_domainMin( smr, ptwXY1, &x1Min ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_domainMax( smr, ptwXY1, &x1Max ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_domainMin( smr, ptwXY2, &x2Min ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_domainMax( smr, ptwXY2, &x2Max ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );

    sMin = x1Min + x2Min;
    sMax = x1Max + x2Max;
    ds = ( sMax - sMin ) / ( n - 1 );
    if( verbose ) {
        fprintf( infoF, "\n\n" );
        fprintf( infoF, "# x1Min = %12.5e\n", x1Min );
        fprintf( infoF, "# x1Max = %12.5e\n", x1Max );
        fprintf( infoF, "# x2Min = %12.5e\n", x2Min );
        fprintf( infoF, "# x2Max = %12.5e\n", x2Max );
        fprintf( infoF, "# sMin = %12.5e\n", sMin );
        fprintf( infoF, "# sMax = %12.5e\n", sMax );
    }


    if( ( d1 = ptwXY_clone( smr, ptwXY1 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ( status = ptwXY_dullEdges( smr, d1, 1e-14, 1e-14, 0 ) ) != nfu_Okay )
        nfut_printSMRErrorExit2p( smr, "Via." );
    if( ( d2 = ptwXY_clone( smr, ptwXY2 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ( status = ptwXY_dullEdges( smr, d2, 1e-14, 1e-14, 0 ) ) != nfu_Okay )
        nfut_printSMRErrorExit2p( smr, "Via." );

    if( ( bfXY = ptwXY_new( smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, n, 10, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( smr, "Via." );

    if( ptwXY_setValueAtX( smr, bfXY, sMin, 0. ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    for( i = 1; i < n - 1; i++ ) {
        s = sMin + i * ds;
        if( s > sMax ) break;

        x11 = s - x2Max;
        if( x11 < x1Min ) x11 = x1Min;
        x12 = x1Min + ( s - sMin );
        if( x12 > x1Max ) x12 = x1Max;

        if( ( clone = ptwXY_clone( smr, d2 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
        for( k = 0, j = clone->length - 1; k < j; k++, j-- ) {
            t = s - clone->points[k].x;
            clone->points[k].x = s - clone->points[j].x;
            clone->points[j].x = t;
            t = clone->points[k].y;
            clone->points[k].y = clone->points[j].y;
            clone->points[j].y = t;
        }
        if( ( 2 * ( clone->length / 2 ) ) != clone->length ) clone->points[k].x = s - clone->points[k].x;
        if( ( f2 = ptwXY_domainSlice( smr, clone, x11, x12, 10, 1 ) ) == NULL )
            nfut_printSMRErrorExit2p( smr, "Via." );
        ptwXY_free( clone );

        if( ( mulXY = ptwXY_mul2_ptwXY( smr, d1, f2 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
        if( ( status = ptwXY_integrateDomain( smr, mulXY, &v ) ) != nfu_Okay ) nfut_printSMRError2p( smr, "Via." );

        if( ptwXY_setValueAtX( smr, bfXY, s, v ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
        ptwXY_free( f2 );
    }
    if( ptwXY_setValueAtX( smr, bfXY, sMax, 0. ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    printIfVerbose( bfXY );
}
/*
************************************************************
*/
static void compareDoubles( double d1, double d2, double eps, const char * const funcName ) {

    double s, d, r;

    if( mode != 0 ) return;
    s = 0.5 * ( fabs( d1 ) + fabs( d2 ) );
    d = d2 - d1;
    r = d;
    if( s != 0 ) r /= s;
    if( fabs( r ) > eps ) fprintf( infoF, "ERROR %s: %s compare, %e %e %e %e %e\n", __FILE__, funcName, d1, d2, s, d, r );
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
