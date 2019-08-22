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
#include <time.h>

#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>


static int verbose = 0, timing = 0;
static char *fmtXY = "%25.17e %25.17e\n";
static FILE *infoF;

static int checkThinning( const char * const label, ptwXYPoints *ptwXY1, double accuracy, int deleteXYs );
static void compareDoubles( double d1, double d2, double eps, const char * const label );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0, i, n;
    ptwXYPoints *dataXY;
    nfu_status status;
    double x, xMax;
    double zeros[] = { -1.5, 0., -1.49999999, 0., 0.5, 0., 2., 0., 2.2, 1, 2.20000001, 0 };
    double triangle[] = { -1.0, 0., 0., 1., 1.0, 0. };
    int nZeros = sizeof( zeros ) / ( 2 * sizeof( double ) ), nTriangles = sizeof( triangle ) / ( 2 * sizeof( double ) );

    infoF = stdout;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else if( strcmp( "-t", argv[iarg] ) == 0 ) {
            timing = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    if( ( dataXY = ptwXY_create( ptwXY_interpolationLinLin, 10., 1e-3, 10, 10, nZeros, zeros, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: zeros creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    checkThinning( "zeros", dataXY, 1e-3, 1 );

    if( ( dataXY = ptwXY_create( ptwXY_interpolationLinLin, 10., 1e-3, 10, 10, nTriangles, triangle, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: zeros creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    checkThinning( "triangle", dataXY, 1e-3, 0 );
    ptwXY_thicken( dataXY, 100, 1e-5, 1. + 1e-5 );
    checkThinning( "triangle thickened", dataXY, 1e-3, 1 );

    n = 1001;
    if( ( dataXY = ptwXY_new( ptwXY_interpolationLinLin, 5, 1e-5, n, 40, &status, 0 ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: sin creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    for( i = 0; i < n; i++ ) {
        x = i * M_PI / ( n - 1 );
        dataXY->points[i].x = x;
        dataXY->points[i].y = sin( x );
    }
    dataXY->length = n;
    checkThinning( "sin", dataXY, 1e-3, 1 );

    n = 100001;
    if( ( dataXY = ptwXY_new( ptwXY_interpolationLinLin, 5, 1e-5, n, 40, &status, 0 ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: sin creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    for( i = 0; i < n; i++ ) {
        x = i * M_PI / ( n - 1 );
        dataXY->points[i].x = x;
        dataXY->points[i].y = sin( x );
    }
    dataXY->length = n;
    checkThinning( "sin with many points", dataXY, 1e-3, 1 );

    n = 1000001;
    if( ( dataXY = ptwXY_new( ptwXY_interpolationLinLin, 5, 1e-5, n, 40, &status, 0 ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: sin creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    xMax = 25 * M_PI;
    for( i = 0; i < n; i++ ) {
        x = i  * xMax / ( n - 1 );
        dataXY->points[i].x = x;
        dataXY->points[i].y = exp( 2. * ( x - 0.2 * xMax ) / xMax ) * sin( x );
    }
    dataXY->length = n;
    checkThinning( "exp * sin with many points", dataXY, 1e-2, 1 );

    exit( errCount );
}
/*
************************************************************
*/
static int checkThinning( const char * const label, ptwXYPoints *ptwXY1, double accuracy, int deleteXYs ) {

    int errCount = 0;
    ptwXYPoints *thinned;
    nfu_status status;
    clock_t time0;
    int64_t i;
    double x, y;

    if( verbose ) fprintf( infoF, "# label = %s\n# accuracy = %e\n", label, accuracy );
    printIfVerbose( ptwXY1 );
    time0 = clock( );
    thinned = ptwXY_thin( ptwXY1, accuracy, &status );
    if( timing ) printf( "# %s: time = %.4f sec", label, ( clock( ) - time0 ) / ( (double) CLOCKS_PER_SEC ) );
    if( thinned == NULL ) {
        if( timing ) printf( "\n" );
        errCount++;
        nfu_printMsg( "ERROR %s: %s convolution status = %d: %s", __FILE__, label, status, nfu_statusMessage( status ) ); }
    else {
        if( timing ) printf( "  length = %d\n", (int) thinned->length );
        printIfVerbose( thinned );
        for( i = 0; i < ptwXY1->length; i++ ) {     /* This logic requires that ptwXY_thin coalesed ptwXY1. */
            x = ptwXY1->points[i].x;
            if( ( status = ptwXY_getValueAtX( thinned, x, &y ) ) != nfu_Okay ) {
                nfu_printErrorMsg( "ERROR %s: ptwXY_getValueAtX status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
                break;
            }
            compareDoubles( y, ptwXY1->points[i].y, accuracy, label );
        }
        ptwXY_free( thinned );
    }

    if( deleteXYs ) ptwXY_free( ptwXY1 );
    return( errCount );
}
/*
************************************************************
*/
static void compareDoubles( double d1, double d2, double eps, const char * const label ) {

    double s, d, r;

    s = 0.5 * ( fabs( d1 ) + fabs( d2 ) );
    d = d2 - d1;
    r = d;
    if( s != 0 ) r /= s;
    if( fabs( r ) > eps ) fprintf( infoF, "ERROR %s: %s compare, %e %e %e %e %e\n", __FILE__, label, d1, d2, s, d, r );
}
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    fprintf( infoF, "# length = %d\n", (int) ptwXY_length( data ) );
    ptwXY_simpleWrite( data, infoF, fmtXY );
    fprintf( infoF, "\n\n" );
}
