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

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_utilities.h>

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int checkTrim( statusMessageReporting *smr, ptwXYPoints *data );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, n, iarg, echo = 0, errCount = 0;
    ptwXYPoints *XYs;
    double x;
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

    if( ( XYs = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 10, 10, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );

    n = 5;
    for( i = 0; i < n; i++ ) {
        x = .2 * i + .5;
        if( ptwXY_setValueAtX( &smr, XYs, x, 0. ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkTrim( &smr, XYs );

    n += 7;
    for( ; i < n; i++ ) {
        x = .2 * i + .5;
        if( ptwXY_setValueAtX( &smr, XYs, x, 2 * i - 20. ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    n += 3;
    for( ; i < n; i++ ) {
        x = .2 * i + .5;
        if( ptwXY_setValueAtX( &smr, XYs, x, 0. ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkTrim( &smr, XYs );

    if( ptwXY_clear( &smr, XYs ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i = 0; i < 12; i++ ) {
        x = .2 * i + .5;
        if( ptwXY_setValueAtX( &smr, XYs, x, 2 * i - 20. ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkTrim( &smr, XYs );

    ptwXY_free( XYs );

    exit( errCount );
}
/*
************************************************************
*/
static int checkTrim( statusMessageReporting *smr, ptwXYPoints *data ) {

    int allPointsZero = 0;
    int64_t i, i1;
    nfu_status status = nfu_Okay;
    ptwXYPoints *trimmed;
    ptwXYPoint *point1, *point2;

    printIfVerbose( data );
    if( ( trimmed = ptwXY_clone( smr, data ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_trim( smr, trimmed ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    printIfVerbose( trimmed );

    for( i1 = 0; i1 < data->length; i1++ ) {
        point1 = ptwXY_getPointAtIndex_Unsafely( data, i1 );
        if( point1->y != 0. ) break;
    }
    if( i1 == data->length ) {      /* All points are zero. */
        allPointsZero = 1;
        i1 = 0; }
    else if( i1 > 0 ) {
        i1--;
    }
    for( i = 0; ( i < trimmed->length ) && ( i1 < data->length ); i++, i1++ ) {
        point1 = ptwXY_getPointAtIndex_Unsafely( data, i1 );
        point2 = ptwXY_getPointAtIndex_Unsafely( trimmed, i );
        if( ( point1->x != point2->x ) || ( point1->y != point2->y ) ) {
                nfu_printErrorMsg( "ERROR %s: point1 and point2 are not the same at index %d, (x,y) = (%e %e) != (%e %e)", 
                        __FILE__, (int) i, point1->x, point1->y, point2->x, point2->y );
        }
        if( allPointsZero ) i1 = data->length - 2;
    }

    ptwXY_free( trimmed );

    return( status );
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
