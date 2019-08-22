/*
# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <ptwXY.h>
#include <nfut_utilities.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>


static int verbose = 0;
static char *fmtXY = "%25.17e %25.17e\n";
static FILE *infoF;

static int doubleMath( statusMessageReporting *smr, ptwXYPoints *ptwXY, double d, const char * const label );
static int doubleMath2( statusMessageReporting *smr, ptwXYPoints *ptwXY, double d, int operator, const char * const label );
static int doubleMath3( statusMessageReporting *smr, ptwXYPoints *ptwXY, double d, const char * const label );
static int compareXYs( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2, const char * const label );
static int compareDoubles( double d1, double d2, double eps, const char * const label );
static void printIfVerbose( ptwXYPoints *data );
static int binaryMath( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0, i;
    ptwXYPoints *fineXYs, *coarseXYs;
    double *fineYs, *coarseYs, accuracy = 1e-3;
    double fineXs[] = { -2.0, -1.5, -1.4, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 2.2, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
    double coarseXs[] = { -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0 };
    int nFineXs= sizeof( fineXs ) / ( sizeof( double ) ), nCoarseXs = sizeof( coarseXs ) / ( sizeof( double ) );
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

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
    if( ( fineXYs = ptwXY_new( &smr, ptwXY_interpolationFlat, NULL, 0., accuracy, 10, 10, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_setXYDataFromXsAndYs( &smr, fineXYs, nFineXs, fineXs, fineYs ) != nfu_Okay )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    doubleMath( &smr, fineXYs, 3.14, "math with pi" );
    doubleMath( &smr, fineXYs, -3.14, "math with -pi" );

    if( ( coarseYs = nfu_malloc( nCoarseXs * sizeof( double ) ) ) == NULL ) nfu_printErrorMsg( "ERROR %s: nfu_malloc-ing nCoarseXs", __FILE__ );
    for( i = 0; i < nCoarseXs ; i++ ) coarseYs [i] = -i + 10;
    if( ( coarseXYs = ptwXY_new( &smr, ptwXY_interpolationFlat, NULL, 0., accuracy, 10, 10, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_setXYDataFromXsAndYs( &smr, coarseXYs, nCoarseXs, coarseXs, coarseYs ) != nfu_Okay )
        nfut_printSMRErrorExit2p( &smr, "Via." );

    binaryMath( &smr, fineXYs, coarseXYs );
    binaryMath( &smr, coarseXYs, fineXYs );

    free( fineYs );
    free( coarseYs );
    ptwXY_free( fineXYs );
    ptwXY_free( coarseXYs );
    exit( errCount );
}
/*
************************************************************
*/
static int doubleMath( statusMessageReporting *smr, ptwXYPoints *ptwXY, double d, const char * const label ) {

    int errs = 0;

    errs += doubleMath2( smr, ptwXY, d, '+', label );
    errs += doubleMath2( smr, ptwXY, d, '-', label );
    errs += doubleMath2( smr, ptwXY, d, '=', label );
    errs += doubleMath2( smr, ptwXY, d, '*', label );
    errs += doubleMath2( smr, ptwXY, d, '/', label );
    errs += doubleMath2( smr, ptwXY, d, '\\', label );
    errs += doubleMath3( smr, ptwXY, d, label );

    return( errs );
}
/*
************************************************************
*/
static int doubleMath2( statusMessageReporting *smr, ptwXYPoints *ptwXY, double d, int operator, const char * const label ) {

    int errs = 0;
    int64_t i;
    ptwXYPoints *XY1, *XY2;
    ptwXYPoint *p1, *p2;
    char Str[128];

    if( verbose ) {
        fprintf( infoF, "# double\n" );
        fprintf( infoF, "# double = %e\n", d );
        fprintf( infoF, "# %c\n", operator );
    }
    printIfVerbose( ptwXY );
    if( ( XY1 = ptwXY_clone( smr, ptwXY ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ( XY2 = ptwXY_clone( smr, ptwXY ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );

    if( operator == '+' ) {
        if( ptwXY_add_double( smr, XY1, d ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
        for( i = 0; i < XY2->length; i++ ) {
            p1 = ptwXY_getPointAtIndex_Unsafely( XY2, i );
            p1->y += d;
        } }
    else if( operator == '-' ) {
        if( ptwXY_sub_doubleFrom( smr, XY1, d ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
        for( i = 0; i < XY2->length; i++ ) {
            p1 = ptwXY_getPointAtIndex_Unsafely( XY2, i );
            p1->y -= d;
        } }
    else if( operator == '=' ) {
        if( ptwXY_sub_fromDouble( smr, XY1, d ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
        for( i = 0; i < XY2->length; i++ ) {
            p1 = ptwXY_getPointAtIndex_Unsafely( XY2, i );
            p1->y = d - p1->y;
        } }
    else if( operator == '*' ) {
        if( ptwXY_mul_double( smr, XY1, d ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
        for( i = 0; i < XY2->length; i++ ) {
            p1 = ptwXY_getPointAtIndex_Unsafely( XY2, i );
            p1->y *= d;
        } }
    else if( operator == '/' ) {
        if( ptwXY_div_doubleFrom( smr, XY1, d ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
        for( i = 0; i < XY2->length; i++ ) {
            p1 = ptwXY_getPointAtIndex_Unsafely( XY2, i );
            p1->y /= d;
        } }
    else if( operator == '\\' ) {
        if( ptwXY_div_fromDouble( smr, XY1, d ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
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
static int doubleMath3( statusMessageReporting *smr, ptwXYPoints *ptwXY, double d, const char * const label ) {

    int errs = 0;
    int64_t i;
    ptwXYPoints *XY1;
    ptwXYPoint *p1, *p2;
    char Str[128];

    if( verbose ) {
        fprintf( infoF, "# all_double\n" );
        fprintf( infoF, "# double = %e\n", d );
    }
    printIfVerbose( ptwXY );
    if( ( XY1 = ptwXY_clone( smr, ptwXY ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );

    if( ptwXY_add_double( smr, XY1, d ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_mul_double( smr, XY1, d ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_sub_doubleFrom( smr, XY1, d ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_div_doubleFrom( smr, XY1, d ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );

    if( ptwXY_mul_double( smr, XY1, d ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_add_double( smr, XY1, d ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_div_doubleFrom( smr, XY1, d ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_sub_doubleFrom( smr, XY1, d ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );

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
static int binaryMath( statusMessageReporting *smr, ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 ) {

    int errs = 0;
    ptwXYPoints *XY1, *XY2;

    if( verbose ) fprintf( infoF, "# binary_add_sub\n" );
    printIfVerbose( ptwXY1 );
    printIfVerbose( ptwXY2 );
    if( ( XY1 = ptwXY_add_ptwXY( smr, ptwXY1, ptwXY2 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    printIfVerbose( XY1 );
    if( ( XY2 = ptwXY_sub_ptwXY( smr, XY1, ptwXY2 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    ptwXY_free( XY1 );
    errs += compareXYs( ptwXY1, XY2, "addition and subtract" );

    if( verbose ) fprintf( infoF, "# binary_mul_div\n" );
    printIfVerbose( ptwXY1 );
    printIfVerbose( ptwXY2 );
    if( ( XY1 = ptwXY_mul_ptwXY( smr, ptwXY1, ptwXY2 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    printIfVerbose( XY1 );
    if( ( XY2 = ptwXY_div_ptwXY( smr, XY1, ptwXY2, 0 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
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
    fprintf( infoF, "# length = %d\n", (int) ptwXY_length( NULL, data ) );
    ptwXY_simpleWrite( data, infoF, fmtXY );
    fprintf( infoF, "\n\n" );
}
