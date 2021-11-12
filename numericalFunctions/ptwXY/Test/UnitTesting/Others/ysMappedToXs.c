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
#include <ptwX.h>
#include <ptwXY.h>
#include <nf_utilities.h>

#define nXYs 8
#define nXs 100

static int verbose = 0;
static char *fmtX = " %19.12e\n";
static char *fmtXY = "%19.12e %19.12e\n";

static int check( statusMessageReporting *smr, ptwXPoints *Xs, ptwXYPoints *XYs, char const *label );
static void printIfVerboseXYs( ptwXYPoints *XYs );
static void printIfVerboseXs( ptwXPoints *Xs );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errCount = 0;
    int64_t i1;
    ptwXYPoints *XYs;
    ptwXYPoint *point;
    ptwXPoints Xs;
    double xValue, domainMin, domainMax, x1, x2;
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

    ptwX_initialize( &smr, &Xs, nXs );

    if( ( XYs = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 10, 10, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i1 = 0; i1 < nXYs; ++i1 ) {
        xValue = 0.2 * i1 - .5;
        if( ptwXY_setValueAtX( &smr, XYs, xValue, xValue * xValue ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    if( ptwXY_domainMin( &smr, XYs, &domainMin ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_domainMax( &smr, XYs, &domainMax ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );

    if( ptwX_setPointAtIndex( &smr, &Xs, 0, domainMin - 5.0 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwX_setPointAtIndex( &smr, &Xs, 1, domainMin - 1.0 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += check( &smr, &Xs, XYs, "0" );

    if( ptwX_setPointAtIndex( &smr, &Xs, 1, domainMin + 0.7  * ( domainMax - domainMin ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += check( &smr, &Xs, XYs, "1" );

    if( ptwX_setPointAtIndex( &smr, &Xs, 0, domainMin + 0.3  * ( domainMax - domainMin ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += check( &smr, &Xs, XYs, "2" );

    if( ptwX_setPointAtIndex( &smr, &Xs, 0, domainMin + 0.03 * ( domainMax - domainMin ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += check( &smr, &Xs, XYs, "3" );

    if( ptwX_setPointAtIndex( &smr, &Xs, 1, domainMin + 0.93 * ( domainMax - domainMin ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += check( &smr, &Xs, XYs, "4" );

    if( ptwX_setPointAtIndex( &smr, &Xs, 1, domainMax + 2.0 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += check( &smr, &Xs, XYs, "5" );

    if( ptwX_setPointAtIndex( &smr, &Xs, 0, domainMax + 1.0 )  != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += check( &smr, &Xs, XYs, "6" );

    if( ptwX_clear( &smr, &Xs ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i1 = 0; i1 < nXYs; ++i1 ) {
        point = ptwXY_getPointAtIndex_Unsafely( XYs, i1 );
        if( ptwX_setPointAtIndex( &smr, &Xs, i1, point->x ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += check( &smr, &Xs, XYs, "7" );

    if( ptwX_clear( &smr, &Xs ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i1 = 0; i1 < nXYs; ++i1 ) {
        point = ptwXY_getPointAtIndex_Unsafely( XYs, i1 );
        x2 = point->x;
        if( i1 > 0 ) {
            if( ptwX_setPointAtIndex( &smr, &Xs, i1 - 1, 0.5 * ( x1 + x2 ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
        }
        x1 = x2;
    }
    errCount += check( &smr, &Xs, XYs, "8" );

    if( ptwX_release( &smr, &Xs ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    ptwXY_free( XYs );

    exit( errCount );
}
/*
************************************************************
*/
static int check( statusMessageReporting *smr, ptwXPoints *Xs, ptwXYPoints *XYs, char const *label ) {

    int errCount = 0;
    int64_t offset;
    ptwXPoints *Ys = NULL;

    if( verbose ) printf( "# New test: %s\n", label );
    printIfVerboseXs( Xs );
    printIfVerboseXYs( XYs );
    if( ( Ys = ptwXY_ysMappedToXs( smr, XYs, Xs, &offset ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( verbose ) printf( "# offset = %d\n", (int) offset );
    printIfVerboseXs( Ys );
    ptwX_free( Ys );

    return( errCount );
}
/*
************************************************************
*/
static void printIfVerboseXYs( ptwXYPoints *XYs ) {

    if( !verbose ) return;
    printf( "# length = %d\n", (int) XYs->length );
    ptwXY_simpleWrite( XYs, stdout, fmtXY );
    printf( "\n\n" );
}
/*
************************************************************
*/
static void printIfVerboseXs( ptwXPoints *Xs ) {

    if( !verbose ) return;
    printf( "# length = %d\n", (int) Xs->length );
    ptwX_simpleWrite( NULL, Xs, stdout, fmtX );
    printf( "\n\n" );
}
