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


static int verbose = 0;
static char *fmtXY = "%25.17e %25.17e\n";
static FILE *infoF;

static int doubleMath( ptwXYPoints *ptwXY, double d, const char * const label );
static int doubleMath2( ptwXYPoints *ptwXY, double d, int operator, const char * const label );
static int doubleMath3( ptwXYPoints *ptwXY, double d, const char * const label );
static int compareXYs( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, const char * const label );
static int compareDoubles( double d1, double d2, double eps, const char * const label );
static void printIfVerbose( ptwXYPoints *data );
static int binaryMath( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0, i;
    ptwXYPoints *fineXYs, *coarseXYs;
    nfu_status status;
    double *fineYs, *coarseYs, accuracy = 1e-3;
    double fineXs[] = { -2.0, -1.5, -1.4, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 2.2, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
    double coarseXs[] = { -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0 };
    int nFineXs= sizeof( fineXs ) / ( sizeof( double ) ), nCoarseXs = sizeof( coarseXs ) / ( sizeof( double ) );

    infoF = stdout;

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

    if( verbose ) fprintf( infoF, "# accuracy = %e\n", accuracy );

    if( ( fineYs = nfu_malloc( nFineXs * sizeof( double ) ) ) == NULL ) nfu_printErrorMsg( "ERROR %s: nfu_malloc-ing fineXYs", __FILE__ );
    for( i = 0; i < nFineXs; i++ ) fineYs[i] = 2 * i + 3;
    if( ( fineXYs = ptwXY_new( ptwXY_interpolationFlat, 0., accuracy, 10, 10, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: dataXY creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_setXYDataFromXsAndYs( fineXYs, nFineXs, fineXs, fineYs ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: ptwXY_setXYDataFromXsAndYs status = %d for fineXYs: %s", __FILE__, status, nfu_statusMessage( status ) );
    doubleMath( fineXYs, 3.14, "math with pi" );
    doubleMath( fineXYs, -3.14, "math with -pi" );

    if( ( coarseYs = nfu_malloc( nCoarseXs * sizeof( double ) ) ) == NULL ) nfu_printErrorMsg( "ERROR %s: nfu_malloc-ing nCoarseXs", __FILE__ );
    for( i = 0; i < nCoarseXs ; i++ ) coarseYs [i] = -i + 10;
    if( ( coarseXYs = ptwXY_new( ptwXY_interpolationFlat, 0., accuracy, 10, 10, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: coarseXYs creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_setXYDataFromXsAndYs( coarseXYs, nCoarseXs, coarseXs, coarseYs ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: ptwXY_setXYDataFromXsAndYs status = %d for coarseXYs: %s", __FILE__, status, nfu_statusMessage( status ) );

    binaryMath( fineXYs, coarseXYs );
    binaryMath( coarseXYs, fineXYs );

    free( fineYs );
    free( coarseYs );
    ptwXY_free( fineXYs );
    ptwXY_free( coarseXYs );
    exit( errCount );
}
/*
************************************************************
*/
static int doubleMath( ptwXYPoints *ptwXY, double d, const char * const label ) {

    int errs = 0;

    errs += doubleMath2( ptwXY, d, '+', label );
    errs += doubleMath2( ptwXY, d, '-', label );
    errs += doubleMath2( ptwXY, d, '=', label );
    errs += doubleMath2( ptwXY, d, '*', label );
    errs += doubleMath2( ptwXY, d, '/', label );
    errs += doubleMath2( ptwXY, d, '\\', label );
    errs += doubleMath3( ptwXY, d, label );

    return( errs );
}
/*
************************************************************
*/
static int doubleMath2( ptwXYPoints *ptwXY, double d, int operator, const char * const label ) {

    int errs = 0;
    int64_t i;
    ptwXYPoints *XY1, *XY2;
    ptwXYPoint *p1, *p2;
    nfu_status status;
    char Str[128];

    if( verbose ) {
        fprintf( infoF, "# double\n" );
        fprintf( infoF, "# double = %e\n", d );
        fprintf( infoF, "# %c\n", operator );
    }
    printIfVerbose( ptwXY );
    if( ( XY1 = ptwXY_clone( ptwXY, &status ) ) == NULL ) nfu_printErrorMsg( "XY1 clone: status = %d: %s", status, nfu_statusMessage( status ) );
    if( ( XY2 = ptwXY_clone( ptwXY, &status ) ) == NULL ) nfu_printErrorMsg( "XY2 clone: status = %d: %s", status, nfu_statusMessage( status ) );

    if( operator == '+' ) {
        if( ( status = ptwXY_add_double( XY1, d ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: %s: ptwXY_add_double, XY1 + %e with status = %d: %s", __FILE__, label, d, status, nfu_statusMessage( status ) );
        for( i = 0; i < XY2->length; i++ ) {
            p1 = ptwXY_getPointAtIndex_Unsafely( XY2, i );
            p1->y += d;
        } }
    else if( operator == '-' ) {
        if( ( status = ptwXY_sub_doubleFrom( XY1, d ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: %s: ptwXY_sub_doubleFrom, XY1 + %e with status = %d: %s", __FILE__, label, d, status, nfu_statusMessage( status ) );
        for( i = 0; i < XY2->length; i++ ) {
            p1 = ptwXY_getPointAtIndex_Unsafely( XY2, i );
            p1->y -= d;
        } }
    else if( operator == '=' ) {
        if( ( status = ptwXY_sub_fromDouble( XY1, d ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: %s: ptwXY_sub_fromDouble, XY1 + %e with status = %d: %s", __FILE__, label, d, status, nfu_statusMessage( status ) );
        for( i = 0; i < XY2->length; i++ ) {
            p1 = ptwXY_getPointAtIndex_Unsafely( XY2, i );
            p1->y = d - p1->y;
        } }
    else if( operator == '*' ) {
        if( ( status = ptwXY_mul_double( XY1, d ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: %s: ptwXY_mul_double, XY1 + %e with status = %d: %s", __FILE__, label, d, status, nfu_statusMessage( status ) );
        for( i = 0; i < XY2->length; i++ ) {
            p1 = ptwXY_getPointAtIndex_Unsafely( XY2, i );
            p1->y *= d;
        } }
    else if( operator == '/' ) {
        if( ( status = ptwXY_div_doubleFrom( XY1, d ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: %s: ptwXY_div_doubleFrom, XY1 + %e with status = %d: %s", __FILE__, label, d, status, nfu_statusMessage( status ) );
        for( i = 0; i < XY2->length; i++ ) {
            p1 = ptwXY_getPointAtIndex_Unsafely( XY2, i );
            p1->y /= d;
        } }
    else if( operator == '\\' ) {
        if( ( status = ptwXY_div_fromDouble( XY1, d ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: %s: ptwXY_div_fromDouble, XY1 + %e with status = %d: %s", __FILE__, label, d, status, nfu_statusMessage( status ) );
        for( i = 0; i < XY2->length; i++ ) {
            p1 = ptwXY_getPointAtIndex_Unsafely( XY2, i );
            p1->y = d / p1->y;
        }
    }

    printIfVerbose( XY1 );
    if( XY1->length != XY2->length ) 
        nfu_printErrorMsg( "ERROR %s: %s for %c: length = %d != length = %d", __FILE__, label, operator, (int) XY1->length, (int) XY2->length );
    for( i = 0; i < XY2->length; i++ ) {
        p1 = ptwXY_getPointAtIndex_Unsafely( XY1, i );
        p2 = ptwXY_getPointAtIndex_Unsafely( XY2, i );
        sprintf( Str, "%s for x with %c", label, operator );
        errs += compareDoubles( p1->x, p2->x, 1e-14, Str );
        sprintf( Str, "%s for y with %c", label, operator );
        errs += compareDoubles( p1->y, p2->y, 1e-14, Str );
    }

    ptwXY_free( XY1 );
    ptwXY_free( XY2 );
    return( errs );
}
/*
************************************************************
*/
static int doubleMath3( ptwXYPoints *ptwXY, double d, const char * const label ) {

    int errs = 0;
    int64_t i;
    ptwXYPoints *XY1;
    ptwXYPoint *p1, *p2;
    nfu_status status;
    char Str[128];

    if( verbose ) {
        fprintf( infoF, "# all_double\n" );
        fprintf( infoF, "# double = %e\n", d );
    }
    printIfVerbose( ptwXY );
    if( ( XY1 = ptwXY_clone( ptwXY, &status ) ) == NULL ) nfu_printErrorMsg( "XY1 clone: status = %d: %s", status, nfu_statusMessage( status ) );

    if( ( status = ptwXY_add_double( XY1, d ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: %s: ptwXY_add_double, XY1 + %e with status = %d: %s", __FILE__, label, d, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_mul_double( XY1, d ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: %s: ptwXY_mul_double, XY1 + %e with status = %d: %s", __FILE__, label, d, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_sub_doubleFrom( XY1, d ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: %s: ptwXY_sub_doubleFrom, XY1 + %e with status = %d: %s", __FILE__, label, d, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_div_doubleFrom( XY1, d ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: %s: ptwXY_div_doubleFrom, XY1 + %e with status = %d: %s", __FILE__, label, d, status, nfu_statusMessage( status ) );

    if( ( status = ptwXY_mul_double( XY1, d ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: %s: ptwXY_mul_double, XY1 + %e with status = %d: %s", __FILE__, label, d, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_add_double( XY1, d ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: %s: ptwXY_add_double, XY1 + %e with status = %d: %s", __FILE__, label, d, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_div_doubleFrom( XY1, d ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: %s: ptwXY_div_doubleFrom, XY1 + %e with status = %d: %s", __FILE__, label, d, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_sub_doubleFrom( XY1, d ) ) != nfu_Okay )
        nfu_printErrorMsg( "ERROR %s: %s: ptwXY_sub_doubleFrom, XY1 + %e with status = %d: %s", __FILE__, label, d, status, nfu_statusMessage( status ) );

    printIfVerbose( XY1 );
    if( XY1->length != ptwXY->length ) 
        nfu_printErrorMsg( "ERROR %s: %s: length = %d != length = %d", __FILE__, label, (int) XY1->length, (int) ptwXY->length );
    for( i = 0; i < ptwXY->length; i++ ) {
        p1 = ptwXY_getPointAtIndex_Unsafely( XY1, i );
        p2 = ptwXY_getPointAtIndex_Unsafely( ptwXY, i );
        sprintf( Str, "all %s for x with", label );
        errs += compareDoubles( p1->x, p2->x, 1e-14, Str );
        sprintf( Str, "all %s for y", label );
        errs += compareDoubles( p1->y, p2->y, 1e-14, Str );
    }

    ptwXY_free( XY1 );
    return( errs );
}
/*
************************************************************
*/
static int binaryMath( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 ) {

    int errs = 0;
    nfu_status status;
    ptwXYPoints *XY1, *XY2;

    if( verbose ) fprintf( infoF, "# binary_add_sub\n" );
    printIfVerbose( ptwXY1 );
    printIfVerbose( ptwXY2 );
    if( ( XY1 = ptwXY_add_ptwXY( ptwXY1, ptwXY2, &status ) ) == NULL )
        nfu_printErrorMsg( "XY1 ptwXY_add_ptwXY: status = %d: %s", status, nfu_statusMessage( status ) );
    printIfVerbose( XY1 );
    if( ( XY2 = ptwXY_sub_ptwXY( XY1, ptwXY2, &status ) ) == NULL )
        nfu_printErrorMsg( "XY2 ptwXY_sub_ptwXY: status = %d: %s", status, nfu_statusMessage( status ) );
    ptwXY_free( XY1 );
    errs += compareXYs( ptwXY1, XY2, "addition and subtract" );

    if( verbose ) fprintf( infoF, "# binary_mul_div\n" );
    printIfVerbose( ptwXY1 );
    printIfVerbose( ptwXY2 );
    if( ( XY1 = ptwXY_mul_ptwXY( ptwXY1, ptwXY2, &status ) ) == NULL )
        nfu_printErrorMsg( "XY1 ptwXY_mul_ptwXY: status = %d: %s", status, nfu_statusMessage( status ) );
    printIfVerbose( XY1 );
    if( ( XY2 = ptwXY_div_ptwXY( XY1, ptwXY2, &status, 0 ) ) == NULL )
        nfu_printErrorMsg( "XY2 ptwXY_div_ptwXY: status = %d: %s", status, nfu_statusMessage( status ) );
    ptwXY_free( XY1 );
    errs += compareXYs( ptwXY1, XY2, "multiplication and division" );

    return( errs );
}
/*
************************************************************
*/
static int compareXYs( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, const char * const label ) {

    int errs = 0;

    printIfVerbose( ptwXY1 );
    printIfVerbose( ptwXY2 );

    ptwXY_free( ptwXY2 );

    return( errs );
}
/*
************************************************************
*/
static int compareDoubles( double d1, double d2, double eps, const char * const label ) {

    double s, d, r;

    s = 0.5 * ( fabs( d1 ) + fabs( d2 ) );
    d = d2 - d1;
    r = d;
    if( s != 0 ) r /= s;
    if( fabs( r ) > eps ) {
        fprintf( infoF, "ERROR %s: %s compare, %e %e %e %e %e\n", __FILE__, label, d1, d2, s, d, r );
        return( 1 );
    }
    return( 0 );
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
