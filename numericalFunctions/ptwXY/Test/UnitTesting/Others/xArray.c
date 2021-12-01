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

#define nXs 99

static int verbose = 0;
static char *fmtX = "%19.12e\n", *fmtXY = "%19.12e %19.12e\n";
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0;
    ptwXYPoints *XYs;
    ptwXPoints *Xs, *xArray;
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

    if( ( Xs = ptwX_new( &smr, nXs ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );

    if( ( XYs = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, nXs / 2, 10, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );

    for( i = 0; i < nXs; i++ ) {
        x = 2.1 * i + 0.1;
        Xs->points[i] = x;
        Xs->length++;
        if( ptwXY_setValueAtX( &smr, XYs, x, i * i ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }

    if( verbose ) {
        printf( "# length = %d\n", (int) XYs->length );
        ptwXY_simpleWrite( XYs, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) Xs->length );
        for( i = 0; i < ptwX_length( &smr, Xs ); i++ ) printf( fmtX, ptwX_getPointAtIndex_Unsafely( Xs, i ) );
        printf( "\n\n" );
    }

    if( ( xArray = ptwXY_getXArray( &smr, XYs ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );

    if( verbose ) {
        printf( "# length = %d\n", (int) nXs );
        for( i = 0; i < nXs; i++ ) printf( fmtX, ptwX_getPointAtIndex_Unsafely( xArray, i ) );
    }

    if( ptwX_length( &smr, Xs ) != ptwX_length( &smr, xArray ) ) nfu_printErrorMsg( "ERROR %s: length not the same", __FILE__ );
    for( i = 0; i < nXs; i++ ) {
        if( Xs->points[i] != xArray->points[i] ) nfu_printErrorMsg( "ERROR %s: at index %d, %e != %e", __FILE__, i, Xs->points[i], xArray->points[i] );
    }

    ptwXY_free( XYs );
    ptwX_free( Xs );
    ptwX_free( xArray );

    exit( EXIT_SUCCESS );
}
