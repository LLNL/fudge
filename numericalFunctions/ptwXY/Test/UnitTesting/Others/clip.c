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

#define nSame 6

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";
static double yMin = -0.9, yMax = 0.4;

static int checkClipping( statusMessageReporting *smr, ptwXYPoints *data );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;
    ptwXYPoints *XY, *expXY, *mulXY;
    double x, expXYs[4], domainMin, domainMax;
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

    if( ( XY = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 10, 10, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i = 0; i < nSame; i++ ) {
        if( ptwXY_setValueAtX( &smr, XY, 0.2 * i - .5, yMax + i - .1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkClipping( &smr, XY );
    if( ptwXY_neg( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkClipping( &smr, XY );
    if( ptwXY_neg( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );

    for( ; i < 2 * nSame + 1; i++ ) {
        if( ptwXY_setValueAtX( &smr, XY, 0.2 * i - .5, yMin - i - .1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkClipping( &smr, XY );
    if( ptwXY_neg( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkClipping( &smr, XY );
    if( ptwXY_neg( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );

    if( ptwXY_clear( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );

    for( i = 0; i < 501; i++ ) {
        x = i * M_PI / 50;
        if( ptwXY_setValueAtX( &smr, XY, x, sin( x ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkClipping( &smr, XY );

    if( ptwXY_domainMin( &smr, XY, &domainMin ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_domainMax( &smr, XY, &domainMax ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    expXYs[0] = domainMin;
    expXYs[1] = 0.;
    expXYs[2] = domainMax;
    expXYs[3] = 1.;
    if( ( expXY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 100, 10, 2, expXYs, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    printIfVerbose( expXY );
    if( ptwXY_exp( &smr, expXY, 1. ) != nfu_Okay ) nfut_printSMRError2p( &smr, "Via." );
    printIfVerbose( expXY );
    if( ( mulXY = ptwXY_mul_ptwXY( &smr, XY, expXY ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkClipping( &smr, mulXY );

    ptwXY_free( XY );
    ptwXY_free( expXY );
    ptwXY_free( mulXY );

    exit( errCount );
}
/*
************************************************************
*/
static int checkClipping( statusMessageReporting *smr, ptwXYPoints *data ) {

    int i, errCount = 0;
    ptwXYPoints *clipped, *u;
    ptwXYPoint *point;
    nfu_status status;
    double x1, y1, x2, y2, x, y, s;

    if( verbose ) {
        printf( "# yMin = %.14e\n", yMin );
        printf( "# yMax = %.14e\n", yMax );
    }
    printIfVerbose( data );
    if( ( clipped = ptwXY_clone( smr, data ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_clip( smr, clipped, yMin, yMax ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    printIfVerbose( clipped );

    if( ( u = ptwXY_union( smr, clipped, data, ptwXY_union_fill ) ) == NULL )
        nfut_printSMRErrorExit2p( smr, "Via." );
    for( i = 0; i < u->length; i++ ) {
        point = ptwXY_getPointAtIndex_Unsafely( u, i );
        if( point->y < yMin ) {
            nfu_printMsg( "ERROR %s: at x = %g, point->y = %g < yMin = %g", __FILE__, point->x, point->y, yMin );
            errCount++; }
        else if( point->y > yMax ) {
            nfu_printMsg( "ERROR %s: at x = %g, point->y = %g > yMax = %g", __FILE__, point->x, point->y, yMax );
            errCount++;
        }
        if( ( status = ptwXY_getValueAtX( smr, data, point->x, &y ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
        if( y < yMin ) {
            if( point->y != yMin ) {
                nfu_printMsg( "ERROR %s: at x, y = (%g, %g), point->y = %g != yMin = %g", __FILE__, point->x, y, point->y, yMin );
                errCount++;
            }
        }
        else if( y > yMax ) {
            if( point->y != yMax ) {
                nfu_printMsg( "ERROR %s: at x, y = (%g, %g), point->y = %g != yMax = %g", __FILE__, point->x, y, point->y, yMax );
                errCount++;
            }
        }
    }

    for( i = 0; i < data->length; i++ ) {       /* test insures that points are added at yMin and yMax crossings. */
        point = ptwXY_getPointAtIndex_Unsafely( data, i );
        x2 = point->x;
        y2 = point->y;
        if( i > 0 ) {
            if( ( ( y1 - yMin ) * ( y2 - yMin ) ) < 0 ) {
                x = ( x1 * ( y2 - yMin ) + x2 * ( yMin - y1 ) ) / ( y2 - y1 );
                if( ( status = ptwXY_getValueAtX( smr, clipped, x, &y ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
                s = 0.5 * ( fabs( y ) + fabs( yMin ) );
                if( fabs( y - yMin ) > 1e-14 * s ) {
                    nfu_printMsg( "ERROR %s: at i = %d, at x = %g, y = %g != yMin = %g", __FILE__, i, x, y, yMin );
                    errCount++;
                }
            }
            if( ( ( y1 - yMax ) * ( y2 - yMax ) ) < 0 ) {
                x = ( x1 * ( y2 - yMax ) + x2 * ( yMax - y1 ) ) / ( y2 - y1 );
                if( ( status = ptwXY_getValueAtX( smr, clipped, x, &y ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
                s = 0.5 * ( fabs( y ) + fabs( yMax ) );
                if( fabs( y - yMax ) > 1e-14 * s ) {
                    nfu_printMsg( "ERROR %s: at i = %d, at x = %g, y = %g != yMax = %g", __FILE__, i, x, y, yMax );
                    errCount++;
                }
            }
        }
        x1 = x2;
        y1 = y2;
    }
    ptwXY_free( u );
    ptwXY_free( clipped );

    return( errCount );
}
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    printf( "# length = %d\n", (int) data->length );
    ptwXY_simpleWrite( data, stdout, fmtXY );
    printf( "\n\n" );
}
