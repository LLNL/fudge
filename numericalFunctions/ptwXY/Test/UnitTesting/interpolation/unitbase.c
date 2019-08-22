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
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include <ptwXY.h>
#include <ptwXY_utilities.h>

static int verbose = 0;
char fmt[] = "%22.14e %22.14e\n";

static void printUnitbasedXY( double w, double wMin, double wMax, ptwXYPoints *p );
void printMsg( const char *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int64_t i;
    int iarg, echo = 0, errCount = 0;
    ptwXYPoints *pXY1, *pXY2, *pl, *pr, *pm1, *pm2, *diff;
    ptwXYPoint *p;
    double y, accuracy = 1e-3, xy1[3*2] = { -1., 0., 0., 1., 1., 0. }, xy2[3*2] = { 8., 0., 10.5, 0.4, 13., 0. };
    nfu_status status;
    ptwXY_interpolation interpolation = ptwXY_interpolationLinLin;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else {
            printMsg( "Error %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );
    
    nfu_setMemoryDebugMode( 0 );

    if( ( pXY1 = ptwXY_create( interpolation, NULL, 5, accuracy, 10, 10,    3, xy1, &status, 0 ) ) == NULL ) 
        printMsg( "pXY1 creation: status = %d: %s", status, nfu_statusMessage( status ) );
    if( ( pXY2 = ptwXY_create( interpolation, NULL, 5, accuracy, 10, 10,    3, xy2, &status, 0 ) ) == NULL ) 
        printMsg( "pXY2 creation: status = %d: %s", status, nfu_statusMessage( status ) );

    if( ( pl = ptwXY_unitbaseInterpolate( 4., 0., pXY1, 20., pXY2, &status ) ) == NULL )
        printMsg( "pl unitbaseInterpolate: status = %d: %s", status, nfu_statusMessage( status ) );

    if( ( pr = ptwXY_unitbaseInterpolate( 12., 0., pXY1, 20., pXY2, &status ) ) == NULL )
        printMsg( "pr unitbaseInterpolate: status = %d: %s", status, nfu_statusMessage( status ) );
    if( ( pm1 = ptwXY_unitbaseInterpolate( 10., 4., pl, 12., pr, &status ) ) == NULL )
        printMsg( "pm1 unitbaseInterpolate: status = %d: %s", status, nfu_statusMessage( status ) );

    if( ( pm2 = ptwXY_unitbaseInterpolate( 10., 0., pXY1, 20., pXY2, &status ) ) == NULL )
        printMsg( "pm2 unitbaseInterpolate: status = %d: %s", status, nfu_statusMessage( status ) );

    if( ( diff = ptwXY_sub_ptwXY( pm1, pm2, &status ) ) == NULL )
        printMsg( "ptwXY_sub_ptwXY: status = %d: %s", status, nfu_statusMessage( status ) );
    for( i = 0; i < ptwXY_length( diff ); i++ ) {
        p = ptwXY_getPointAtIndex_Unsafely( diff, i );
        status = ptwXY_getValueAtX( pm1, p->x, &y );
        switch( status ) {
        case nfu_Okay :
            if( fabs( p->y ) > 1e-12 * fabs( y ) )
                printMsg( "pm1 and pm2 differ at x  = %e by %e:, pm1.y = %e", p->x, p->y, y );
            break;
        case nfu_XOutsideDomain :
            if( ( i == 0 ) || ( i == ( ptwXY_length( diff ) - 1 ) ) ) continue;
            printMsg( "ptwXY_getValueAtX status = %d: %s", status, nfu_statusMessage( status ) );
        default :
            printMsg( "ptwXY_getValueAtX status = %d: %s", status, nfu_statusMessage( status ) );
        }
    }
    ptwXY_free( diff );

    if( verbose ) {
        printf( "\n\n" );
        printf( "# length = %d\n", (int) pXY1->length );
        ptwXY_simpleWrite( pXY1, stdout, fmt );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) pXY2->length );
        ptwXY_simpleWrite( pXY2, stdout, fmt );
        printUnitbasedXY(  4., 0., 20., pl );
        printUnitbasedXY( 12., 0., 20., pr );
        printUnitbasedXY( 10., 4., 12., pm1 );
        printUnitbasedXY( 10., 0., 20., pm2 );
    }

    ptwXY_free( pXY1 );
    ptwXY_free( pXY2 );
    ptwXY_free( pl );
    ptwXY_free( pr );
    ptwXY_free( pm1 );
    ptwXY_free( pm2 );

    exit( errCount ? EXIT_FAILURE : EXIT_SUCCESS );
}
/*
****************************************************************
*/
static void printUnitbasedXY( double w, double wMin, double wMax, ptwXYPoints *p ) {

    printf( "\n\n" );
    printf( "# w = %e\n", w );
    printf( "# wMin = %e\n", wMin );
    printf( "# wMax = %e\n", wMax );
    printf( "# length = %d\n", (int) p->length );
    ptwXY_simpleWrite( p, stdout, fmt );
}
/*
****************************************************************
*/
void printMsg( const char *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    fprintf( stderr, "\n" );
    va_end( args );
    exit( EXIT_FAILURE );
}
