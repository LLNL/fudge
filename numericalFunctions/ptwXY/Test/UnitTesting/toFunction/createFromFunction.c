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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>

static int verbose = 0, errorCounter, zeroCounter;
static double accuracy = 1e-3, xMax = 20.;

int toPointwise( ptwXPoints *Xs, FILE *f, int biSectionMax );
nfu_status xSinXX_callback( double x, double *y, void *argList );
void printDiff( FILE *f, double x, double y );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, errCount = 0, echo = 0;
    double x;
    ptwXPoints Xs;
    FILE *f;

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

    f = fopen( "e", "w" );

    ptwX_setup( &Xs, 1000 );
    ptwX_setPointAtIndex( &Xs, 0, 1. );
    ptwX_setPointAtIndex( &Xs, 1, 10. );
    errCount += toPointwise( &Xs, f, 16 );
    if( verbose ) printf( "\n\n\n" );

    ptwX_setup( &Xs, 1000 );
    ptwX_setPointAtIndex( &Xs, 0, 1. );
    for( i = 1; ; i++ ) {
        x = sqrt( i * M_PI );
        if( x <= 1 ) continue;
        if( x >= xMax ) break;
        ptwX_setPointAtIndex( &Xs, i, x );
    }
    ptwX_setPointAtIndex( &Xs, i, xMax );
    errCount += toPointwise( &Xs, f, 12 );

    fclose( f );
    exit( errCount );
}
/*
************************************************************
*/
int toPointwise( ptwXPoints *Xs, FILE *f, int biSectionMax ) {

    int i, n, errCount = 0;
    double x, y, s;
    ptwXYPoint *p;
    ptwXYPoints *XYs;
    nfu_status status;

    errorCounter = zeroCounter = 0;
    if( verbose ) {
        printf( "# Xs\n" );
        printf( "# length = %d\n", (int) ptwX_length( Xs ) );
        for( i = 0; i < (int) ptwX_length( Xs ); i++ ) printf( "# Xs[%3d] = %.17e\n", i, ptwX_getPointAtIndex_Unsafely( Xs, i ) );
        printf( "# accuracy = %e\n", accuracy );
        printf( "# biSectionMax = %d\n", biSectionMax );
    }

    XYs = ptwXY_createFromFunction2( Xs, xSinXX_callback, NULL, accuracy, 1, biSectionMax, &status );
    n = (int) ptwXY_length( XYs );
    p = ptwXY_getPointAtIndex_Unsafely( XYs, 0 );
    x = p->x;
    for( i = 1; i < n; i++ ) {
        p = ptwXY_getPointAtIndex_Unsafely( XYs, i );
        s = fabs( p->x ) + fabs( x );
        if( ( p->x - x ) <= ( s * ClosestAllowXFactor * DBL_EPSILON ) ) {
            fprintf( stderr, "Values too close at indices %d and %d: %.17e  %.17e, delta = %e, rel. delta = %e\n", i, i + 1, x, p->x, p->x - x, 2 * ( p->x - x ) / s );
            errCount++;
        }
        x = p->x;
    }

    p = ptwXY_getPointAtIndex_Unsafely( XYs, 0 );
    x = p->x;
    y = p->y;
    fprintf( f, "# Errors\n" );
    printDiff( f, x, y );
    for( i = 1; i < n; i++ ) {
        p = ptwXY_getPointAtIndex_Unsafely( XYs, i );
        printDiff( f, 0.5 * ( x + p->x ), 0.5 * ( y + p->y ) );
        printDiff( f, p->x, p->y );
        x = p->x;
        y = p->y;
    }

    if( errorCounter != zeroCounter ) {
        fprintf( stderr, "errorCounter %d != zeroCounter = %d\n", errorCounter, zeroCounter );
        errCount++;
    }

    if( verbose ) {
        printf( "# length = %d\n", (int) ptwXY_length( XYs ) );
        ptwXY_simplePrint( XYs, "%.12e %.12e\n" );
    }

    ptwX_release( Xs );
    ptwXY_free( XYs );

    return( errCount );
}
/*
************************************************************
*/
nfu_status xSinXX_callback( double x, double *y, void *argList ) {

    *y = x * sin( x * x );
    return( nfu_Okay );
}
/*
************************************************************
*/
void printDiff( FILE *f, double x, double y ) {

    int doPrint = verbose;
    double n, d, r;
    char *s = "";

    xSinXX_callback( x, &n, NULL );
    d = n - y;
    r = 0.5 * ( fabs( y ) + fabs( n ) );
    if( r == 0 ) {
        r = 1.; }
    else {
        r = d / r;
    }
    if( fabs( r ) > accuracy ) {
        doPrint = 1;
        errorCounter++;
        s = "*";
        if( fabs( r ) > 2 * accuracy ) {
            s = "**";
            if( fabs( r ) > 5 * accuracy ) {
                s = "***";
                if( fabs( r ) > 10 * accuracy ) s = "****";
            }
        }
    }
    if( doPrint ) {
        fprintf( f, "%25.17e %25.17e %25.17e, %+e %+e # %s\n", x, y, n, d, r, s );
    }
    if( y == 0. ) zeroCounter++;
}
