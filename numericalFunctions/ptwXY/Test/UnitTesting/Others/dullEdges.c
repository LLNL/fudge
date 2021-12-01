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
#include <math.h>

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_utilities.h>

#define nSame 6

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int checkDullEdges( statusMessageReporting *smr, ptwXYPoints *data, double lowerEps, double upperEps );
static int checkDullEdges2( statusMessageReporting *smr, ptwXYPoints *data, double lowerEps, double upperEps, int positiveXOnly );
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
    for( i = 0; i < nSame; i++ ) {
        if( ptwXY_setValueAtX( &smr, XY, 0.2 * i - .5, 0.7 + i + .1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkDullEdges( &smr, XY, 1e-10, 1e-10 );
    if( ptwXY_neg( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkDullEdges( &smr, XY, 1e-10, 1e-10 );

    if( ptwXY_clear( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i = 0; i < nSame; i++ ) {
        if( ptwXY_setValueAtX( &smr, XY, 0.2 * i, 0.7 + i + .1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkDullEdges( &smr, XY, 1e-10, 1e-10 );

    if( ptwXY_clear( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i = 0; i < nSame; i++ ) {
        if( ptwXY_setValueAtX( &smr, XY, -0.2 * i, 0.7 + i + .1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkDullEdges( &smr, XY, 1e-10, 1e-10 );

    ptwXY_free( XY );

    exit( errCount );
}
/*
************************************************************
*/
static int checkDullEdges( statusMessageReporting *smr, ptwXYPoints *data, double lowerEps, double upperEps ) {

    int errCount;

    errCount = checkDullEdges2( smr, data, lowerEps, upperEps, 0 );
    errCount += checkDullEdges2( smr, data,  lowerEps, upperEps, 1 );
    errCount += checkDullEdges2( smr, data, -lowerEps, upperEps, 0 );
    errCount += checkDullEdges2( smr, data, -lowerEps, upperEps, 1 );
    errCount += checkDullEdges2( smr, data,  lowerEps, -upperEps, 0 );
    errCount += checkDullEdges2( smr, data,  lowerEps, -upperEps, 1 );
    errCount += checkDullEdges2( smr, data, -lowerEps, -upperEps, 0 );
    errCount += checkDullEdges2( smr, data, -lowerEps, -upperEps, 1 );
    return( errCount );
}
/*
************************************************************
*/
static int checkDullEdges2( statusMessageReporting *smr, ptwXYPoints *data, double lowerEps, double upperEps, int positiveXOnly ) {

    int errCount = 0;
    ptwXYPoints *dullEdges;
    nfu_status status;

    printIfVerbose( data );
    if( verbose ) {
        printf( "# lowerEps = %.14e\n", lowerEps );
        printf( "# upperEps = %.14e\n", upperEps );
        printf( "# positiveXOnly = %d\n", positiveXOnly );
    }
    if( ( dullEdges = ptwXY_clone( smr, data ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ( status = ptwXY_dullEdges( smr, dullEdges, lowerEps, upperEps, positiveXOnly ) ) != nfu_Okay )
        nfut_printSMRErrorExit2p( smr, "Via." );
    printIfVerbose( dullEdges );

    ptwXY_free( dullEdges );

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
