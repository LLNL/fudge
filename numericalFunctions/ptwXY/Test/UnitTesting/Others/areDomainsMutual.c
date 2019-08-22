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

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ptwXY.h>
#include <nf_utilities.h>

static int verbose = 0;
static char *fmtXY = "%17.8e%17.8e\n";

static int addPointAndCheck( double x, double y, ptwXYPoints *data1, ptwXYPoints *data2, int areMutual );
static int checkAreMutual( ptwXYPoints *data1, ptwXYPoints *data2, int areMutual );
static int checkAreMutual2( ptwXYPoints *data1, ptwXYPoints *data2, int areMutual );
static int checkAreMutual3( ptwXYPoints *data1, ptwXYPoints *data2, int areMutual );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errCount = 0;
    nfu_status status;
    ptwXYPoints *XY1, *XY2;
    double flat1[] = { 1.0, 1.0, 9.0, 1.0 }, triangle1[] = { -1.0, 0.0, 0.0, 1.0, 0.5, 0.0 };
    double xMin = flat1[0], xMax = flat1[2], xMid = 0.5 * ( xMin + xMax );
    double wedge1[] = { 10., 0., 11., 1. }, wedge2[] = { xMid, 0., xMax, 1. }, wedge3[] = { xMin, 1., xMid, 0. };
    double triangle2[] = { xMid - 1., 0., xMid, 1., xMid + 1., 0 };

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

    if( ( XY1 = ptwXY_create( ptwXY_interpolationLinLin, 4, 1.e-3, 20, 10, 2, flat1, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: creating XY1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( XY2 = ptwXY_create( ptwXY_interpolationLinLin, 4, 1.e-3, 20, 10, 3, triangle1, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: creating XY2, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    errCount += checkAreMutual( XY1, XY2, 0 );

    errCount += addPointAndCheck( flat1[0], flat1[1], XY1, XY2, 0 );
    errCount += addPointAndCheck( xMid, flat1[3], XY1, XY2, 0 );
    errCount += addPointAndCheck( flat1[2], flat1[3], XY1, XY2, 0 );
    errCount += addPointAndCheck( 2 * flat1[2], flat1[3], XY1, XY2, 0 );
    ptwXY_free( XY2 );

    if( ( XY2 = ptwXY_create( ptwXY_interpolationLinLin, 4, 1.e-3, 20, 10, 3, triangle1, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: creating XY2, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    errCount += addPointAndCheck( flat1[0], 0., XY1, XY2, 0 );
    errCount += addPointAndCheck( flat1[2], 0., XY1, XY2, 0 );
    ptwXY_free( XY2 );

    if( ( XY2 = ptwXY_create( ptwXY_interpolationLinLin, 4, 1.e-3, 20, 10, 2, wedge1, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: creating XY2, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    errCount += checkAreMutual( XY1, XY2, 0 );
    errCount += addPointAndCheck( flat1[2], 0., XY1, XY2, 0 );
    errCount += addPointAndCheck( xMid, 0., XY1, XY2, 0 );
    errCount += addPointAndCheck( flat1[0], 0., XY1, XY2, 0 );
    ptwXY_free( XY2 );

    errCount += checkAreMutual( XY1, XY1, 1 );

    if( ( XY2 = ptwXY_create( ptwXY_interpolationLinLin, 4, 1.e-3, 20, 10, 2, wedge2, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: creating XY2, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    errCount += checkAreMutual( XY1, XY2, 1 );
    ptwXY_free( XY2 );

    if( ( XY2 = ptwXY_create( ptwXY_interpolationLinLin, 4, 1.e-3, 20, 10, 2, wedge3, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: creating XY2, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    errCount += checkAreMutual( XY1, XY2, 1 );
    ptwXY_free( XY2 );

    if( ( XY2 = ptwXY_create( ptwXY_interpolationLinLin, 4, 1.e-3, 20, 10, 3, triangle2, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: creating XY2, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    errCount += checkAreMutual( XY1, XY2, 1 );
    ptwXY_free( XY2 );

    ptwXY_free( XY1 );

    exit( errCount );
}
/*
************************************************************
*/
static int addPointAndCheck( double x, double y, ptwXYPoints *data1, ptwXYPoints *data2, int areMutual ) {

    nfu_status status;

    if( ( status = ptwXY_setValueAtX( data2, x, y ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX on data2, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    return( checkAreMutual( data1, data2, areMutual ) );
}
/*
************************************************************
*/
static int checkAreMutual( ptwXYPoints *data1, ptwXYPoints *data2, int areMutual ) {

    nfu_status status;
    ptwXYPoints *n1, *n2;
    int errCount = 0;

    errCount += checkAreMutual2( data1, data2, areMutual );

    if( ( n1 = ptwXY_clone( data1, &status ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: cloning data1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_neg( n1 ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: negating data1, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    errCount += checkAreMutual2( n1, data2, areMutual );

    if( ( n2 = ptwXY_clone( data2, &status ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: cloning data2, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_neg( n2 ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: negating data2, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    errCount += checkAreMutual2( data1, n2, areMutual );

    errCount += checkAreMutual2( n1, n2, areMutual );

    ptwXY_free( n1 );
    ptwXY_free( n2 );

    return( errCount );
}
/*
************************************************************
*/
static int checkAreMutual2( ptwXYPoints *d1, ptwXYPoints *d2, int areMutual ) {

    int errCount = checkAreMutual3( d1, d2, areMutual );

    return( errCount + checkAreMutual3( d2, d1, areMutual ) );
}
/*
************************************************************
*/
static int checkAreMutual3( ptwXYPoints *d1, ptwXYPoints *d2, int areMutual ) {

    int errCount = 0;
    nfu_status status;

    if( verbose ) printf( "# areMutual = %d\n", areMutual );
    printIfVerbose( d1 );
    printIfVerbose( d2 );

    if( ( status = ptwXY_areDomainsMutual( d1, d2 ) ) != nfu_Okay ) {
        if( areMutual || ( status != nfu_domainsNotMutual ) ) {
            errCount++;
            nfu_printMsg( "ERROR %s: ptwXY_areDomainsMutual, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
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
