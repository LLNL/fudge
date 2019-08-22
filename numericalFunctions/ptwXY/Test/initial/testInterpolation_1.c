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

#define size 1001

static int verbose = 0;
char fmt[] = "%22.14e %22.14e\n";

void printMsg( char const *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0;
    ptwXYPoints *l2, *lLog, *lLin;
    double xy[2*size], xy2[2*2], a, x, y, f;
    nfu_status status;
    FILE *ff;

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

    xy2[0] = 1.0;
    xy2[1] = 1.0;
    xy2[2] = 100.0;
    xy2[3] = 10.0;

    a = log( xy2[3] / xy2[1] ) / log( xy2[2] / xy2[0] );
    f = pow( xy2[2] / xy2[0], 1. / ( size - 1 ) );
    x = xy2[0];

    xy[0] = xy2[0];
    xy[1] = xy2[1];
    for( i = 1; i < size - 1; i++ ) {
        x *= f;
        y = pow( x / xy2[0], a );
        xy[2*i] = x;
        xy[2*i+1] = y;
    }
    xy[2 * size - 2] = xy2[2];
    xy[2 * size - 1] = xy2[3];

    if( ( l2    = ptwXY_create( ptwXY_interpolationLogLog, NULL, 5, 1e-3, 10, 10,    2, xy2, &status, 0 ) ) == NULL ) printMsg( "l2 creation: status = %d: %s", 
        status, nfu_statusMessage( status ) );
    if( ( lLog = ptwXY_create( ptwXY_interpolationLogLog, NULL, 5, 1e-3, 25, 10, size,  xy, &status, 0 ) ) == NULL ) printMsg( "lLog creation: status = %d: %s", 
        status, nfu_statusMessage( status ) );

    ff = fopen( "curve_u_loglog.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) l2->length );
    ptwXY_simpleWrite( l2, ff, fmt );
    fclose( ff );

    ff = fopen( "curve_u_loglog_dense.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) lLog->length );
    ptwXY_simpleWrite( lLog, ff, fmt );
    fclose( ff );

    if( ( lLin = ptwXY_toOtherInterpolation( l2, ptwXY_interpolationLinLin, 1e-3, &status ) ) == NULL ) 
        printMsg( "lLin creation: status = %d: %s", status, nfu_statusMessage( status ) );

    ff = fopen( "curve_u_interpolatedToLinear.dat", "w" );
    fprintf( ff, "# length = %d\n", (int) lLin->length );
    ptwXY_simpleWrite( lLin, ff, fmt );
    fclose( ff );

    ptwXY_free( l2 );
    ptwXY_free( lLog );
    ptwXY_free( lLin );

    exit( EXIT_SUCCESS );
}
/*
****************************************************************
*/
void printMsg( char const *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    fprintf( stderr, "\n" );
    va_end( args );
    exit( EXIT_FAILURE );
}
