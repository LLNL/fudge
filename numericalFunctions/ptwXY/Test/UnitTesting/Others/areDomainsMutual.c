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

static int verbose = 0;
static char *fmtXY = "%17.8e%17.8e\n";

static int addPointAndCheck( statusMessageReporting *smr, double x, double y, ptwXYPoints *data1, ptwXYPoints *data2, int areMutual );
static int checkAreMutual( statusMessageReporting *smr, ptwXYPoints *data1, ptwXYPoints *data2, int areMutual );
static int checkAreMutual2( statusMessageReporting *smr, ptwXYPoints *data1, ptwXYPoints *data2, int areMutual );
static int checkAreMutual3( statusMessageReporting *smr, ptwXYPoints *data1, ptwXYPoints *data2, int areMutual );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errCount = 0;
    ptwXYPoints *XY1, *XY2;
    double flat1[] = { 1.0, 1.0, 9.0, 1.0 }, triangle1[] = { -1.0, 0.0, 0.0, 1.0, 0.5, 0.0 };
    double xMin = flat1[0], xMax = flat1[2], xMid = 0.5 * ( xMin + xMax );
    double wedge1[] = { 10., 0., 11., 1. }, wedge2[] = { xMid, 0., xMax, 1. }, wedge3[] = { xMin, 1., xMid, 0. };
    double triangle2[] = { xMid - 1., 0., xMid, 1., xMid + 1., 0 };
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

    if( ( XY1 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 20, 10, 2, flat1, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( XY2 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 20, 10, 3, triangle1, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );

    errCount += checkAreMutual( &smr, XY1, XY2, 0 );

    errCount += addPointAndCheck( &smr, flat1[0], flat1[1], XY1, XY2, 0 );
    errCount += addPointAndCheck( &smr, xMid, flat1[3], XY1, XY2, 0 );
    errCount += addPointAndCheck( &smr, flat1[2], flat1[3], XY1, XY2, 0 );
    errCount += addPointAndCheck( &smr, 2 * flat1[2], flat1[3], XY1, XY2, 0 );
    ptwXY_free( XY2 );

    if( ( XY2 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 20, 10, 3, triangle1, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += addPointAndCheck( &smr, flat1[0], 0., XY1, XY2, 0 );
    errCount += addPointAndCheck( &smr, flat1[2], 0., XY1, XY2, 0 );
    ptwXY_free( XY2 );

    if( ( XY2 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 20, 10, 2, wedge1, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkAreMutual( &smr, XY1, XY2, 0 );
    errCount += addPointAndCheck( &smr, flat1[2], 0., XY1, XY2, 0 );
    errCount += addPointAndCheck( &smr, xMid, 0., XY1, XY2, 0 );
    errCount += addPointAndCheck( &smr, flat1[0], 0., XY1, XY2, 0 );
    ptwXY_free( XY2 );

    errCount += checkAreMutual( &smr, XY1, XY1, 1 );

    if( ( XY2 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 20, 10, 2, wedge2, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkAreMutual( &smr, XY1, XY2, 1 );
    ptwXY_free( XY2 );

    if( ( XY2 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 20, 10, 2, wedge3, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkAreMutual( &smr, XY1, XY2, 1 );
    ptwXY_free( XY2 );

    if( ( XY2 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 20, 10, 3, triangle2, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkAreMutual( &smr, XY1, XY2, 1 );
    ptwXY_free( XY2 );

    ptwXY_free( XY1 );

    exit( errCount );
}
/*
************************************************************
*/
static int addPointAndCheck( statusMessageReporting *smr, double x, double y, ptwXYPoints *data1, ptwXYPoints *data2, int areMutual ) {

    if( ptwXY_setValueAtX( smr, data2, x, y ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );

    return( checkAreMutual( smr, data1, data2, areMutual ) );
}
/*
************************************************************
*/
static int checkAreMutual( statusMessageReporting *smr, ptwXYPoints *data1, ptwXYPoints *data2, int areMutual ) {

    ptwXYPoints *n1, *n2;
    int errCount = 0;

    errCount += checkAreMutual2( smr, data1, data2, areMutual );

    if( ( n1 = ptwXY_clone( smr, data1 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_neg( smr, n1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    errCount += checkAreMutual2( smr, n1, data2, areMutual );

    if( ( n2 = ptwXY_clone( smr, data2 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_neg( smr, n2 ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    errCount += checkAreMutual2( smr, data1, n2, areMutual );

    errCount += checkAreMutual2( smr, n1, n2, areMutual );

    ptwXY_free( n1 );
    ptwXY_free( n2 );

    return( errCount );
}
/*
************************************************************
*/
static int checkAreMutual2( statusMessageReporting *smr, ptwXYPoints *d1, ptwXYPoints *d2, int areMutual ) {

    int errCount = checkAreMutual3( smr, d1, d2, areMutual );

    return( errCount + checkAreMutual3( smr, d2, d1, areMutual ) );
}
/*
************************************************************
*/
static int checkAreMutual3( statusMessageReporting *smr, ptwXYPoints *d1, ptwXYPoints *d2, int areMutual ) {

    int errCount = 0;
    nfu_status status;

    if( verbose ) printf( "# areMutual = %d\n", areMutual );
    printIfVerbose( d1 );
    printIfVerbose( d2 );

    if( ( status = ptwXY_areDomainsMutual( smr, d1, d2 ) ) != nfu_Okay ) {
        if( areMutual || ( status != nfu_domainsNotMutual ) ) {
            errCount++;
            nfut_printSMRError2p( smr, "Via." );
        } }
    else {
        if( !areMutual ) {
            errCount++;
            nfu_printMsg( "ERROR %s: ptwXY_areDomainsMutual, is true and it should be false", __FILE__ );
        }
    }

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
