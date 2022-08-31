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

static int clone( statusMessageReporting *smr, int n1, int allocatedSize, int overflowAllocatedSize );
static int compareXYs( statusMessageReporting *smr, ptwXYPoints *XY1, ptwXYPoints *XY2 );
static int compareValues( int64_t i, double x1, double y1, double x2, double y2 );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errCount = 0;
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

    errCount += clone( &smr, 5, 10, 10 );
    errCount += clone( &smr, 16, 10, 10 );
    errCount += clone( &smr, 23, 10, 10 );
    errCount += clone( &smr, 56, 10, 10 );
    errCount += clone( &smr, 5, 100, 10 );
    errCount += clone( &smr, 16, 100, 10 );
    errCount += clone( &smr, 123, 100, 10 );
    errCount += clone( &smr, 556, 100, 10 );
    errCount += clone( &smr, 5, 100, 50 );
    errCount += clone( &smr, 16, 100, 50 );
    errCount += clone( &smr, 123, 100, 50 );
    errCount += clone( &smr, 556, 100, 50 );
    return( errCount );
}
/*
************************************************************
*/
static int clone( statusMessageReporting *smr, int n1, int allocatedSize, int overflowAllocatedSize ) {

    int i1, errCount = 0, nonOverflowLength1, nonOverflowLength2;
    double x1;
    ptwXYPoints *XYSrc, *XYClone, *XYClone2;

    XYSrc = ptwXY_new( smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, allocatedSize, overflowAllocatedSize, 0 );
    if( XYSrc == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    for( i1 = 0; i1 < n1; ++i1 ) {
        x1 = 1 - 2 * drand48( );
        if( ptwXY_setValueAtX( smr, XYSrc, x1, x1 * x1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    }

    nonOverflowLength1 = (int) ptwXY_getNonOverflowLength( smr, XYSrc );
    if( verbose ) printf( "nonOverflowLength1 = %d\n", nonOverflowLength1 );

    if( ( XYClone2 = ptwXY_clone2( smr, XYSrc ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    nonOverflowLength2 = (int) ptwXY_getNonOverflowLength( smr, XYSrc );
    if( nonOverflowLength1 != nonOverflowLength2 ) ++errCount;

    if( ( XYClone = ptwXY_clone( smr, XYSrc ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );

    printIfVerbose( XYSrc );
    printIfVerbose( XYClone );
    printIfVerbose( XYClone2 );

    if( ptwXY_simpleCoalescePoints( smr, XYSrc ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_simpleCoalescePoints( smr, XYClone ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );

    errCount += compareXYs( smr, XYSrc, XYClone );
    errCount += compareXYs( smr, XYSrc, XYClone2 );

    ptwXY_free( XYSrc );
    ptwXY_free( XYClone );
    ptwXY_free( XYClone2 );
    return( errCount );
}
/*
************************************************************
*/
static int compareXYs( statusMessageReporting *smr, ptwXYPoints *XY1, ptwXYPoints *XY2 ) {

    int errCount = 0;
    int64_t i, n = ptwXY_length( NULL, XY1 );
    double x1, y1, x2, y2;

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
