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
#include <float.h>

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>

static int verbose = 0, errorCounter, zeroCounter;
static double accuracy = 1e-3, xMax = 20.;

int toPointwise( statusMessageReporting *smr, ptwXPoints *Xs, FILE *f, int biSectionMax );
nfu_status xSinXX_callback( statusMessageReporting *smr, double x, double *y, void *argList );
void printDiff( FILE *f, double x, double y );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, errCount = 0, echo = 0;
    double x;
    ptwXPoints Xs;
    FILE *f;
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

    f = fopen( "e", "w" );

    ptwX_initialize( &smr, &Xs, 1000 );
    ptwX_setPointAtIndex( &smr, &Xs, 0, 1. );
    ptwX_setPointAtIndex( &smr, &Xs, 1, 10. );
    errCount += toPointwise( &smr, &Xs, f, 16 );
    if( verbose ) printf( "\n\n\n" );

    ptwX_initialize( &smr, &Xs, 1000 );
    ptwX_setPointAtIndex( &smr, &Xs, 0, 1. );
    for( i = 1; ; i++ ) {
        x = sqrt( i * M_PI );
        if( x <= 1 ) continue;
        if( x >= xMax ) break;
        ptwX_setPointAtIndex( &smr, &Xs, i, x );
    }
    ptwX_setPointAtIndex( &smr, &Xs, i, xMax );
    errCount += toPointwise( &smr, &Xs, f, 12 );

    fclose( f );
    exit( errCount );
}
/*
************************************************************
*/
int toPointwise( statusMessageReporting *smr, ptwXPoints *Xs, FILE *f, int biSectionMax ) {

    int i, n, errCount = 0;
    double x, y, s;
    ptwXYPoint *p;
    ptwXYPoints *XYs;

    errorCounter = zeroCounter = 0;
    if( verbose ) {
        printf( "# Xs\n" );
        printf( "# length = %d\n", (int) ptwX_length( smr, Xs ) );
        for( i = 0; i < (int) ptwX_length( smr, Xs ); i++ ) printf( "# Xs[%3d] = %.17e\n", i, ptwX_getPointAtIndex_Unsafely( Xs, i ) );
        printf( "# accuracy = %e\n", accuracy );
        printf( "# biSectionMax = %d\n", biSectionMax );
    }

    if( ( XYs = ptwXY_createFromFunction2( smr, Xs, xSinXX_callback, NULL, accuracy, 1, biSectionMax ) ) == NULL )
        nfut_printSMRErrorExit2p( smr, "Via." );
    n = (int) ptwXY_length( smr, XYs );
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
        printf( "# length = %d\n", (int) ptwXY_length( smr, XYs ) );
        ptwXY_simplePrint( XYs, "%.12e %.12e\n" );
    }

    ptwX_release( smr, Xs );
    ptwXY_free( XYs );

    return( errCount );
}
/*
************************************************************
*/
nfu_status xSinXX_callback( statusMessageReporting *smr, double x, double *y, void *argList ) {

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

    xSinXX_callback( NULL, x, &n, NULL );
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
