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

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int checkMerge( statusMessageReporting *smr, ptwXYPoints *XYs, double epsilon, int n );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i1, i2, iarg, echo = 0, errCount = 0;
    ptwXYPoints *XYs;
    double x, epsilon = 1e-10;
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

    for( i1 = 0, x = 1; i1 < 5; i1++, x += epsilon / 8 ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, ( x - 1 ) / epsilon ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    for( ; i1 < 50; i1++, x += 2 * epsilon ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, ( x - 1 ) / epsilon ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkMerge( &smr, XYs, epsilon, 45 );

    for( i1 = 0, x = 1; i1 < 45; i1++, x += 2 * epsilon ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, ( x - 1 ) / epsilon ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    for( ; i1 < 50; i1++, x += epsilon / 8 ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, ( x - 1 ) / epsilon ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkMerge( &smr, XYs, epsilon, 46 );

    for( i1 = 0, x = 1; i1 < 15; i1++, x += 2 * epsilon ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, ( x - 1 ) / epsilon ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    for( ; i1 < 30; i1++, x += epsilon / 8 ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, ( x - 1 ) / epsilon ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    for( ; i1 < 50; i1++, x += 2 * epsilon ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, ( x - 1 ) / epsilon ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkMerge( &smr, XYs, epsilon, 36 );

    for( i1 = 0, x = 1; i1 < 100; i1++, x += epsilon / 8 ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, i1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkMerge( &smr, XYs, epsilon, 13 );

    for( i1 = 0, x = 1; i1 < 2; i1++, x += epsilon ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, i1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkMerge( &smr, XYs, epsilon / 3, 2 );

    for( i1 = 0, x = 1; i1 < 2; i1++, x += epsilon ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, i1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkMerge( &smr, XYs, 3 * epsilon, 2 );

    for( i2 = 0; i2 < 2; ++i2 ) {
        for( i1 = 0, x = 1; i1 < 3; i1++, x += epsilon ) {
            if( ptwXY_setValueAtX( &smr, XYs, x, i1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
        }
        errCount += checkMerge( &smr, XYs, ( 0.333 + 2.667 * i2 ) * epsilon, 3 - i2);
    }

    ptwXY_free( XYs );
    exit( errCount );
}
/*
************************************************************
*/
static int checkMerge( statusMessageReporting *smr, ptwXYPoints *XYs, double epsilon, int n ) {

    int errCount = 0;
    nfu_status status;

    if( verbose ) printf( "# epsilon = %e\n", epsilon );
    printIfVerbose( XYs );
    if( ( status = ptwXY_mergeClosePoints( smr, XYs, epsilon ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    printIfVerbose( XYs );
    if( ptwXY_length( smr, XYs ) != n ) {
        nfu_printMsg( "ERROR %s: ptwXY_length( smr, XYs ) = %d != n = %d", __FILE__, ptwXY_length( smr, XYs ), n );
        errCount++;
    }
    if( ptwXY_clear( smr, XYs ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
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
