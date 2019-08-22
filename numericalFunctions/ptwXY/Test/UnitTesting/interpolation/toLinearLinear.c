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

void printMsg( const char *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, xlog = 0, ylog = 0, errCount;
    ptwXYPoints *pSparse, *pDense, *pLinear;
    double accuracy = 1e-3, xys[2*2] = { 1.0, 2.0, 10.0, 100.0 };
    nfu_status status;
    FILE *ff;
    ptwXY_interpolation interpolation = ptwXY_interpolationLinLin;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-xlog", argv[iarg] ) == 0 ) {
            xlog = 1; }
        else if( strcmp( "-ylog", argv[iarg] ) == 0 ) {
            ylog = 1; }
        else {
            printMsg( "Error %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s%s%s\n", __FILE__, xlog ? " -xlog" : "", ylog ? " -ylog" : "" );
    
    nfu_setMemoryDebugMode( 0 );

    if( xlog ) {
        interpolation = ptwXY_interpolationLogLin;
        if( ylog ) interpolation = ptwXY_interpolationLogLog; }
    else if( ylog ) {
        interpolation = ptwXY_interpolationLinLog;
    }

    if( ( pSparse = ptwXY_create( interpolation, NULL, 5, accuracy, 10, 10,    2, xys, &status, 0 ) ) == NULL ) 
        printMsg( "pSparse creation: status = %d: %s", status, nfu_statusMessage( status ) );
    if( ( pDense = ptwXY_clone( pSparse, &status ) ) == NULL ) printMsg( "pDense creation: status = %d: %s", status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_thicken( pDense, 4000, 0., 1 + 2.e-4 ) ) != nfu_Okay ) printMsg( "thicken: status = %d: %s", status, nfu_statusMessage( status ) );

    ff = fopen( "curve_sparse.dat", "w" );
    fprintf( ff, "# xlog = %d\n", xlog );
    fprintf( ff, "# ylog = %d\n", ylog );
    fprintf( ff, "# accuracy = %e\n", accuracy );
    fprintf( ff, "# length = %d\n", (int) pSparse->length );
    ptwXY_simpleWrite( pSparse, ff, fmt );
    fclose( ff );

    ff = fopen( "curve_dense.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) pDense->length );
    ptwXY_simpleWrite( pDense, ff, fmt );
    fclose( ff );

    if( ( pLinear = ptwXY_toOtherInterpolation( pSparse, ptwXY_interpolationLinLin, accuracy, &status ) ) == NULL ) 
        printMsg( "pLinear creation: status = %d: %s", status, nfu_statusMessage( status ) );

    ff = fopen( "curve_linear.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) pLinear->length );
    ptwXY_simpleWrite( pLinear, ff, fmt );
    fclose( ff );

    if( ( errCount = nfu_ptwXY_cmp( pDense, pLinear, verbose, accuracy ) ) )
        nfu_printMsg( "Error %s: nfu_ptwXY_cmp found %d differences", __FILE__, errCount );

    ptwXY_free( pSparse );
    ptwXY_free( pDense );
    ptwXY_free( pLinear );

    exit( errCount ? EXIT_FAILURE : EXIT_SUCCESS );
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
