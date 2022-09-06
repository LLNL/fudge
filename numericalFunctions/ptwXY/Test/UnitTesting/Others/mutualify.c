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

static int verbose = 0;
static char *fmtXY = "%17.8e%17.8e\n";

static int checkMutualify( statusMessageReporting *smr, ptwXYPoints *data );
static int checkMutualify2( statusMessageReporting *smr, ptwXYPoints *data, int64_t i1, int64_t i2 );
static int checkMutualify3( statusMessageReporting *smr, ptwXYPoints *d1, ptwXYPoints *d2, int64_t i1, int64_t i2 );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;
    ptwXYPoints *XY;
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
    for( i = 0; i < 10; i++ ) {
        if( ptwXY_setValueAtX( &smr, XY, 0.2 * i - .5, 0.7 + i + .1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkMutualify( &smr, XY );
    if( ptwXY_neg( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkMutualify( &smr, XY );

    ptwXY_free( XY );

    exit( errCount );
}
/*
************************************************************
*/
static int checkMutualify( statusMessageReporting *smr, ptwXYPoints *data ) {

    int errCount = 0;
    int64_t n = data->length - 1;

    errCount += checkMutualify2( smr, data, 2, n );
    errCount += checkMutualify2( smr, data, 0, n - 3 );
    return( errCount );
}
/*
************************************************************
*/
static int checkMutualify2( statusMessageReporting *smr, ptwXYPoints *data, int64_t i1, int64_t i2 ) {

    int errCount = 0;
    ptwXYPoints *clone, *sliced;

    if( ( clone = ptwXY_clone( smr, data ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ( sliced = ptwXY_slice( NULL, data, i1, i2, 0 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    errCount += checkMutualify3( smr, clone, sliced, i1, i2 );

    if( ( clone = ptwXY_clone( smr, data ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ( sliced = ptwXY_slice( smr, data, i1, i2, 0 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    errCount += checkMutualify3( smr, sliced, clone, i1, i2 );

    return( errCount );
}
/*
************************************************************
*/
static int checkMutualify3( statusMessageReporting *smr, ptwXYPoints *d1, ptwXYPoints *d2, int64_t i1, int64_t i2 ) {

    int errCount = 0, positiveXOnly = 1;
    double lowerEps = 1e-6, upperEps = 1e-6;
    nfu_status status;

    if( verbose ) {
        printf( "# i1 = %d\n", (int) i1 );
        printf( "# i2 = %d\n", (int) i2 );
        printf( "# lowerEps = %.14e\n", lowerEps );
        printf( "# upperEps = %.14e\n", upperEps );
        printf( "# positiveXOnly = %d\n", positiveXOnly );
    }
    printIfVerbose( d1 );
    printIfVerbose( d2 );

    if( ( status = ptwXY_mutualifyDomains( smr, d1, lowerEps, upperEps, positiveXOnly, d2, lowerEps, upperEps, positiveXOnly ) ) != nfu_Okay )
        nfut_printSMRErrorExit2p( smr, "Via." );

    printIfVerbose( d1 );
    printIfVerbose( d2 );

    if( ( status = ptwXY_areDomainsMutual( smr, d1, d2 ) ) != nfu_Okay ) {
        errCount++;
        nfu_printMsg( "ERROR %s: ptwXY_MutualifyDomains, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }

    ptwXY_free( d1 );
    ptwXY_free( d2 );

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
