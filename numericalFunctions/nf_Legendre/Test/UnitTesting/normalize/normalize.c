/*
# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
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

    Cls[0] = scale;
    for( i = 0; i < MAX_ORDER; i++ ) Cls[i+1] = 0.5 * Cls[i];

    if( ( nfL = nf_Legendre_new( 10, maxOrder, Cls, &status ) ) == NULL ) {
        nfu_printErrorMsg( "nf_Legendre_new failed with status = %d: '%s'", status, nfu_statusMessage( status ) ); }
    else {
        if( ( XYs = nf_Legendre_to_ptwXY( nfL, accuracy, biSectionMax, 0, &status ) ) == NULL ) {
            nfu_printErrorMsg( "nf_Legendre_to_ptwXY failed with status = %d: '%s'", status, nfu_statusMessage( status ) ); }
        else {
            integral1 = ptwXY_integrateDomain( XYs, &status );
            ptwXY_free( XYs );
            nf_Legendre_normalize( nfL );
            if( ( XYs = nf_Legendre_to_ptwXY( nfL, accuracy, biSectionMax, 0, &status ) ) == NULL ) {
                nfu_printErrorMsg( "nf_Legendre_to_ptwXY failed with status = %d: '%s'", status, nfu_statusMessage( status ) ); }
            else {
                integral2 = ptwXY_integrateDomain( XYs, &status );
                if( verbose ) printf( "integral1 = %.16e,  integral2 = %.16e, %e\n", integral1, integral2, 1.0 - integral2 );
                if( fabs( 1.0 - integral2 ) > accuracy ) errCount++;
                ptwXY_free( XYs );
            }
        }
        nf_Legendre_free( nfL );
        Cls[0] = 0;
        if( ( nfL = nf_Legendre_new( 10, maxOrder, Cls, &status ) ) == NULL ) {
            nfu_printErrorMsg( "nf_Legendre_new failed with status = %d: '%s'", status, nfu_statusMessage( status ) ); }
        else {
            if( nf_Legendre_normalize( nfL ) != nfu_divByZero ) errCount++;
            nf_Legendre_free( nfL );
        }
    }

    return( errCount );
}
