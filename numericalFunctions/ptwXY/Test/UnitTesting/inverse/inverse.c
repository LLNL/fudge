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

#define M_SQRT2PI ( 2. * M_SQRT2 / M_2_SQRTPI )

static int verbose = 0;
static char *fmtXY = "%25.17e %25.17e\n";
static FILE *infoF;

static int invert( statusMessageReporting *smr, ptwXYPoints *ptwXY1 );
static void printIfVerbose( ptwXYPoints *data );
static int compare( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errCount = 0;
    ptwXYPoints *ptwXY1, *ptwXY2;
    double accuracy = 1e-3;
    double twoPoints[4] = { 1., 1., 10., 2. };
    int nTwoPoints = sizeof( twoPoints ) / sizeof( twoPoints[0] ) / 2;
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
    if( echo ) printf( "%s\n", __FILE__ );

    if( ( ptwXY1 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, 100, 10, nTwoPoints, twoPoints, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += invert( &smr, ptwXY1 );

    if( ( ptwXY1 = ptwXY_create( &smr, ptwXY_interpolationLogLin, NULL, 10., accuracy, 100, 10, nTwoPoints, twoPoints, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += invert( &smr, ptwXY1 );

    if( ( ptwXY2 = ptwXY_create( &smr, ptwXY_interpolationLogLin, NULL, 10., accuracy, 100, 10, nTwoPoints, twoPoints, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( ptwXY1 = ptwXY_toOtherInterpolation( &smr, ptwXY2, ptwXY_interpolationLinLin, accuracy ) ) == NULL )
        nfut_printSMRError2p( &smr, "Via." );
    ptwXY_free( ptwXY2 );
    errCount += invert( &smr, ptwXY1 );

    if( ( ptwXY1 = ptwXY_create( &smr, ptwXY_interpolationLinLog, NULL, 10., accuracy, 100, 10, nTwoPoints, twoPoints, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += invert( &smr, ptwXY1 );

    if( ( ptwXY2 = ptwXY_create( &smr, ptwXY_interpolationLinLog, NULL, 10., accuracy, 100, 10, nTwoPoints, twoPoints, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( ptwXY1 = ptwXY_toOtherInterpolation( &smr, ptwXY2, ptwXY_interpolationLinLin, accuracy ) ) == NULL )
        nfut_printSMRError2p( &smr, "Via." );
    ptwXY_free( ptwXY2 );
    errCount += invert( &smr, ptwXY1 );

    if( ( ptwXY1 = ptwXY_create( &smr, ptwXY_interpolationLogLog, NULL, 10., accuracy, 100, 10, nTwoPoints, twoPoints, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += invert( &smr, ptwXY1 );

    if( ( ptwXY2 = ptwXY_create( &smr, ptwXY_interpolationLogLog, NULL, 10., accuracy, 100, 10, nTwoPoints, twoPoints, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( ptwXY1 = ptwXY_toOtherInterpolation( &smr, ptwXY2, ptwXY_interpolationLinLin, accuracy ) ) == NULL )
        nfut_printSMRError2p( &smr, "Via." );
    ptwXY_free( ptwXY2 );
    errCount += invert( &smr, ptwXY1 );

    exit( errCount );
}
/*
************************************************************
*/
static int invert( statusMessageReporting *smr, ptwXYPoints *ptwXY1 ) {

    int errCount;
    ptwXYPoints *ptwXY2, *ptwXY3;

    printIfVerbose( ptwXY1 );
    if( ( ptwXY2 = ptwXY_inverse( smr, ptwXY1 ) ) == NULL ) 
        nfut_printSMRError2( smr, __FILE__, __LINE__, __func__, "Via." );
    printIfVerbose( ptwXY2 );
    if( ( ptwXY3 = ptwXY_inverse( smr, ptwXY2 ) ) == NULL ) 
        nfut_printSMRError2( smr, __FILE__, __LINE__, __func__, "Via." );
    printIfVerbose( ptwXY3 );
    errCount = compare( ptwXY1, ptwXY3 );
    ptwXY_free( ptwXY1 );
    ptwXY_free( ptwXY2 );
    ptwXY_free( ptwXY3 );
    return( errCount );
}
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    fprintf( infoF, "# length = %d\n", (int) ptwXY_length( NULL, data ) );
    fprintf( infoF, "# interpolation = %d\n", ptwXY_getInterpolation( data ) );
    fprintf( infoF, "# interpolation string = %s\n", ptwXY_getInterpolationString( data ) );
    ptwXY_simpleWrite( data, infoF, fmtXY );
    fprintf( infoF, "\n\n" );
}
/*
************************************************************
*/
static int compare( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 ) {

    int64_t i1;
    ptwXYPoint *point1, *point2;

    if( ptwXY_length( NULL, ptwXY1 ) != ptwXY_length( NULL, ptwXY2 ) ) return( 1 );

    for( i1 = 0; i1 < ptwXY_length( NULL, ptwXY1 ); ++i1 ) {
        point1 = ptwXY_getPointAtIndex_Unsafely( ptwXY1, i1 );
        point2 = ptwXY_getPointAtIndex_Unsafely( ptwXY2, i1 );
        if( point1->x != point2->x ) return( 1 );
        if( point1->y != point2->y ) return( 1 );
    }
    return( 0 );
}
