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

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_utilities.h>

#define allocatedSize 100

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int compareXYs( statusMessageReporting *smr, ptwXYPoints *XY1, ptwXYPoints *XY2 );
static int compareXYsToCList( statusMessageReporting *smr, ptwXYPoints *XY1, int64_t nPoints, double *xy );
static int compareValues( int64_t i, double x1, double y1, double x2, double y2 );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;
    int64_t returnedPoints;
    double x, xy[2 * allocatedSize];
    ptwXYPoints *XYSrc, *XYDesc;
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

    if( ( XYSrc = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 10, 10, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );

    for( i = 0, x = 1; i < 45; i++, x += 1.1 ) {
        if( ptwXY_setValueAtX( &smr, XYSrc, x, x * x ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }

    if( ( XYDesc = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 10, 10, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );

    if( ptwXY_copy( &smr, XYDesc, XYSrc ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += compareXYs( &smr, XYDesc, XYSrc );
    ptwXY_free( XYDesc );

    if( ptwXY_copyToC_XY( &smr, XYSrc, 0, ptwXY_length( &smr, XYSrc ), allocatedSize, &returnedPoints, xy ) != nfu_Okay )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += compareXYsToCList( &smr, XYSrc, returnedPoints, xy );

    ptwXY_free( XYSrc );
    exit( errCount );
}
/*
************************************************************
*/
static int compareXYs( statusMessageReporting *smr, ptwXYPoints *XY1, ptwXYPoints *XY2 ) {

    int errCount = 0;
    int64_t i, n = ptwXY_length( NULL, XY1 );
    double x1, y1, x2, y2;

    printIfVerbose( XY1 );
    printIfVerbose( XY2 );
    if( n != ptwXY_length( NULL, XY2 ) ) {
        errCount++;
        nfu_printMsg( "ERROR %s: compareXYs, len( XY1 ) = %d != len( XY2 ) = %d", __FILE__, (int) n, (int) ptwXY_length( NULL, XY2 ) ); }
    else {
        for( i = 0; i < n; i++ ) {
            if( ptwXY_getXYPairAtIndex( smr, XY1, i, &x1, &y1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
            if( ptwXY_getXYPairAtIndex( smr, XY2, i, &x2, &y2 ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
            errCount += compareValues( i, x1, y1, x2, y2 );
        }
    }

    return( errCount );
}
/*
************************************************************
*/
static int compareXYsToCList( statusMessageReporting *smr, ptwXYPoints *XY1, int64_t nPoints, double *xy ) {

    int errCount = 0;
    int64_t i, n = ptwXY_length( NULL, XY1 );
    double x1, y1, *p;

    if( n != nPoints ) {
        errCount++;
        nfu_printMsg( "ERROR %s: compareXYsToCList, len( XY1 ) = %d != nPoints = %d", __FILE__, (int) n, (int) nPoints ); }
    else {
        for( i = 0, p = xy; i < n; i++, p += 2 ) {
            if( ptwXY_getXYPairAtIndex( smr, XY1, i, &x1, &y1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
            errCount += compareValues( i, x1, y1, *p, p[1] );
        }
    }

    return( errCount );
}
/*
************************************************************
*/
static int compareValues( int64_t i, double x1, double y1, double x2, double y2 ) {

    if( ( x1 != x2 ) || ( y1 != y2 ) ) {
        nfu_printMsg( "ERROR %s: at index %3d ( x1 = %.17e != x2 = %.17e ) or ( y1 = %.17e != y2 = %.17e )", __FILE__, (int) i, x1, x2, y1, y2 );
        return( 1 );
    }
    return( 0 );
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
