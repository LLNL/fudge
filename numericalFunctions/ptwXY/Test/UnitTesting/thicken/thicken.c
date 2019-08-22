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

#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>

static int verbose = 0;
static char *fmt = "%10.5f %12.5f\n";

static int thicken( ptwXYPoints *u, int length, int sectionSubdivideMax, double dx, double fx, char *limitStr );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    double uXY[] = { 1.0, 3.0, 10, 30. };
    ptwXYPoints *u;
    int nuXY = sizeof( uXY ) / ( 2 * sizeof( double ) ), iarg, errCount = 0, echo = 0;
    nfu_status status;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfu_printErrorMsg( "Error %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    if( ( u = ptwXY_create( ptwXY_interpolationLinLin, NULL, 6, 1e-3, 10, 10, nuXY, uXY, &status, 0 ) ) == NULL ) nfu_printErrorMsg( "Error %s: u creation status = %d: %s",
            __FILE__, status, nfu_statusMessage( status ) );
    if( verbose ) {
        fprintf( stdout, "# length = %d\n", (int) u->length );
        ptwXY_simpleWrite( u, stdout, fmt );
        fprintf( stdout, "\n\n" );
    }

    errCount += thicken( u, 11,  10, 0.4, 1.1, "sectionSubdivideMax" );         /* Limit is sectionSubdivideMax. */
    errCount += thicken( u, 24, 100, 0.4, 11., "dx" );                          /* Limit is dx. */
    errCount += thicken( u, 26, 100, 1.4, 1.1, "fx" );                          /* Limit is fx. */
    errCount += thicken( u, 32, 100, 0.4, 1.1, "dx and fx" );                   /* Limit is dx. */

    uXY[2] = 1.0001;
    uXY[3] = 9.;
    if( ( status = ptwXY_setXYData( u, nuXY, uXY ) ) != nfu_Okay )
        nfu_printErrorMsg( "Error %s: ptwXY_setXYData status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    u->interpolation = ptwXY_interpolationLinLog;
    errCount += thicken( u, 101, 100, 0.0000001, 1.0000001, "log-log" );        /* Limit is sectionSubdivideMax. */

    ptwXY_free( u );

    exit( errCount );
}
/*
************************************************************
*/
static int thicken( ptwXYPoints *u, int length, int sectionSubdivideMax, double dx, double fx, char *limitStr ) {

    int i, errs = 0, errCount;
    ptwXYPoints *v;
    ptwXYPoint *point;
    nfu_status status;

    if( verbose ) printf( "# limitStr = '%s': sectionSubdivideMax = %3d  dx = %8.5f  fx = %8.5f\n", limitStr, sectionSubdivideMax, dx, fx );
    if( ( v = ptwXY_clone( u, &status ) ) == NULL ) nfu_printErrorMsg( "Error %s: u clone status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    if( ( status = ptwXY_thicken( v, sectionSubdivideMax, dx, fx ) ) != nfu_Okay ) nfu_printErrorMsg( "Error %s: u thicken status = %d: %s", 
        __FILE__, status, nfu_statusMessage( status ) );
    if( ptwXY_length( v ) != length ) {
        nfu_printMsg( "Error %s: ptwXY_length( v ) = %d != %d for limit '%s'", __FILE__, (int) ptwXY_length( v ), length, limitStr );
        errs++;
    }
    if( verbose ) {
        double x1, x2;

        printf( "# length = %d\n", (int) v->length );
        for( i = 0; i < (int) ptwXY_length( v ); i++ ) {
            point = ptwXY_getPointAtIndex( v, i );
            x2 = point->x;
            printf( "%10.8f %12.8f", x2, point->y );
            if( i != 0 ) printf( "%12.5e %15.8e", x2 - x1, x2 / x1 );
            printf( "\n" );
            x1 = x2;
        }
    }
    if( verbose ) fprintf( stdout, "\n\n" );
    if( ( errCount = nfu_ptwXY_cmp( u, v, verbose, 1e-15 ) ) ) {
        errs += abs( errCount );
        nfu_printMsg( "Error %s: nfu_ptwXY_cmp found %d differences for limit '%s'", __FILE__, errCount, limitStr );
    }

    ptwXY_free( v );

    return( errs );
}
