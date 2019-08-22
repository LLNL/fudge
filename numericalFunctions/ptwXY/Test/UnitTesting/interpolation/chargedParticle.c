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

#if 0
nfu_status chargedParticleGetValue( void *argList, double x, double *y, double x1, double y1, double x2, double y2 );
#endif
void printMsg( const char *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errCount = 0, nPoints;
    ptwXYPoints *pSparse, *pDense;
    double accuracy = 1e-2;
    double xys1[] = { /* Following is data from ENDF/B-VII.1 reaction H2(H2,n)He3. */
        1e2, 1.1231e-58,          1e3, 5.5988e-18,          1e4, 8.8216e-6,           2e4, 2.7734e-4,           3e4, 1.1747e-3, 
        4e4, 2.6785e-3,           5e4, 4.6105e-3,           6e4, 6.805e-3,            7e4, 9.1418e-3,           8e4, 0.011541, 
        9e4, 0.01395,             1e5, 0.016339,            1.5e5, 0.027495,          2e5, 0.037113,            2.5e5, 0.045358 };
    nfu_status status = nfu_Okay;
    FILE *ff;
#include "chargedParticle.h"

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

    nPoints = sizeof( xys1 ) / sizeof( xys1[0] ) / 2;
    if( ( pSparse = ptwXY_create( ptwXY_interpolationOther, "charged-particle", 5, accuracy, 10, 10, nPoints, xys1, &status, 0 ) ) == NULL ) 
        printMsg( "pSparse creation: status = %d: %s", status, nfu_statusMessage( status ) );
    ff = fopen( "curve_sparse.dat", "w" );
    fprintf( ff, "# accuracy = %e\n", accuracy );
    fprintf( ff, "# length = %d\n", (int) pSparse->length );
    ptwXY_simpleWrite( pSparse, ff, fmt );
    fclose( ff );

    nPoints = sizeof( xys1_answer ) / sizeof( xys1_answer[0] ) / 2;
    if( ( pDense = ptwXY_create( ptwXY_interpolationLinLin, "charged-particle", 5, accuracy, 10, 10, nPoints, xys1_answer, &status, 0 ) ) == NULL ) 
        printMsg( "pSparse creation: status = %d: %s", status, nfu_statusMessage( status ) );
    ff = fopen( "curve_dense.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) pDense->length );
    ptwXY_simpleWrite( pDense, ff, fmt );
    fclose( ff );

    ptwXY_free( pSparse );
    ptwXY_free( pDense );

    exit( errCount ? EXIT_FAILURE : EXIT_SUCCESS );
}
#if 0
/*
****************************************************************
*/
nfu_status chargedParticleGetValue( void *argList, double x, double *y, double x1, double y1, double x2, double y2 ) {

    double A, B, T = *((double *) argList);

    B = log( x2 * y2 / ( x1 * y1 ) ) / ( 1. / sqrt( x1 - T ) - 1. / sqrt( x2 - T ) );
    A = x1 * y1 * exp( B / sqrt( x1 - T ) );
    *y = A * exp( - B / sqrt( x - T ) ) / x;
    return( nfu_Okay );
}
#endif
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
