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
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include <nfut_utilities.h>
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
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

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

    if( ( pXY1 = ptwXY_create( &smr, interpolation, NULL, 5, accuracy, 10, 10,    3, xy1, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( pXY2 = ptwXY_create( &smr, interpolation, NULL, 5, accuracy, 10, 10,    3, xy2, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );

    if( ( pl = ptwXY_unitbaseInterpolate( &smr, 4., 0., pXY1, 20., pXY2, 1 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );

    if( ( pr = ptwXY_unitbaseInterpolate( &smr, 12., 0., pXY1, 20., pXY2, 1 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( pm1 = ptwXY_unitbaseInterpolate( &smr, 10., 4., pl, 12., pr, 1 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );

    if( ( pm2 = ptwXY_unitbaseInterpolate( &smr, 10., 0., pXY1, 20., pXY2, 1 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );

    if( ( diff = ptwXY_sub_ptwXY( &smr, pm1, pm2 ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i = 0; i < ptwXY_length( &smr, diff ); i++ ) {
        p = ptwXY_getPointAtIndex_Unsafely( diff, i );
        if( ( status = ptwXY_getValueAtX( &smr, pm1, p->x, &y ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
        switch( status ) {
        case nfu_Okay :
            if( fabs( p->y ) > 1e-12 * fabs( y ) )
                printMsg( "pm1 and pm2 differ at x  = %e by %e:, pm1.y = %e", p->x, p->y, y );
            break;
        case nfu_XOutsideDomain :
            if( ( i == 0 ) || ( i == ( ptwXY_length( &smr, diff ) - 1 ) ) ) continue;
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
