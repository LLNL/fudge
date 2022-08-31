/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/
/*
    This routine 
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_Legendre.h>

#define MAX_ORDER 16

static int verbose = 0, biSectionMax = 12;
static double accuracy = 1e-5;

static int norm( int maxOrder, double scale );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "nf_Legendre: %s\n", __FILE__ );


    for( i = 0; i <= MAX_ORDER; i++ ) errCount += norm( i, 0.412 * i + 0.21 );

    if( errCount ) fprintf( stderr, "%s FAILED\n", __FILE__ );
    return( errCount );
}
/*
************************************************************
*/
static int norm( int maxOrder, double scale ) {

    int i, errCount = 0;
    double integral1, integral2, Cls[MAX_ORDER+1];
    nf_Legendre *nfL;
    ptwXYPoints *XYs;
    nfu_status status;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    Cls[0] = scale;
    for( i = 0; i < MAX_ORDER; i++ ) Cls[i+1] = 0.5 * Cls[i];

    if( ( nfL = nf_Legendre_new( &smr, 10, maxOrder, Cls ) ) == NULL ) {
        nfut_printSMRError2p( &smr, "Via." ); }
    else {
        if( ( XYs = nf_Legendre_to_ptwXY( &smr, nfL, accuracy, biSectionMax, 0 ) ) == NULL ) {
            nfut_printSMRError2p( &smr, "Via." ); }
        else {
            status = ptwXY_integrateDomain( &smr, XYs, &integral1 );
            ptwXY_free( XYs );
            if( status != nfu_Okay ) {
                nfut_printSMRError2p( &smr, "Via." ); }
            else {
                if( nf_Legendre_normalize( &smr, nfL ) != nfu_Okay ) {
                    nfut_printSMRError2p( &smr, "Via." ); }
                else {
                    if( ( XYs = nf_Legendre_to_ptwXY( &smr, nfL, accuracy, biSectionMax, 0 ) ) == NULL ) {
                        nfut_printSMRError2p( &smr, "Via." ); }
                    else {
                        if( ptwXY_integrateDomain( &smr, XYs, &integral2 ) != nfu_Okay ) {
                            nfut_printSMRError2p( &smr, "Via." ); }
                        else {
                            if( verbose ) printf( "integral1 = %.16e,  integral2 = %.16e, %e\n", integral1, integral2, 1.0 - integral2 );
                            if( fabs( 1.0 - integral2 ) > accuracy ) errCount++;
                        }
                        ptwXY_free( XYs );
                    }
                }
            }
        }
        nf_Legendre_free( nfL );
        Cls[0] = 0;
        if( ( nfL = nf_Legendre_new( &smr, 10, maxOrder, Cls ) ) == NULL ) {
            nfut_printSMRError2p( &smr, "Via." ); }
        else {
            if( nf_Legendre_normalize( NULL, nfL ) != nfu_divByZero ) errCount++;
            nf_Legendre_free( nfL );
        }
    }

    return( errCount );
}
