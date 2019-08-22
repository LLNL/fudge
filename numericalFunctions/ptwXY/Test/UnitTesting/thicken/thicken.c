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

    if( ( u = ptwXY_create( ptwXY_interpolationLinLin, 6, 1e-3, 10, 10, nuXY, uXY, &status, 0 ) ) == NULL ) nfu_printErrorMsg( "Error %s: u creation status = %d: %s",
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
