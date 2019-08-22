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

#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>

static int verbose = 0;

static int integrate( double x1, double y1, double x2, double y2, double nn, double gn, double ng, double gg );
static int integrate2( char *label, ptwXY_interpolation interpolation, double x1, double y1, double x2, double y2, double answer );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0;
    double flat[4] = { 1., 3., 4., 3. }, slope[4] = { 1., 3., 4., 9. };

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

    errCount += integrate( flat[0], flat[1], flat[2], flat[3], 9., 9., 9., 9. );
    errCount += integrate( slope[0], slope[1], slope[2], slope[3], 18., 20.015744631999329, 16.384306079283071, 18.410234412974571 );
    slope[2] = slope[0] * ( 1 + 1e-3 );
    errCount += integrate( slope[0], slope[1], slope[2], slope[3], 0.0059999999999993392, 0.0060004997501579282, 0.0054614353597604226, 0.0054619253462386856 );
    slope[2] = slope[0] * ( 1 + 1e-4 * ( 1. + 1e-8 ) );
    errCount += integrate( slope[0], slope[1], slope[2], slope[3], 6.0000000600046732e-4, 6.000050057508588e-4, 5.4614354143796317e-4, 5.4614844341595715e-4 );
    slope[2] = slope[0] * ( 1 + 1e-4 * ( 1. - 1e-8 ) );
    errCount += integrate( slope[0], slope[1], slope[2], slope[3], 5.9999999399940052e-4, 6.0000499374931635e-4, 5.4614353051412137e-4, 5.4614843249191933e-4 );
    slope[2] = slope[0] * ( 1 + 1e-5 );
    errCount += integrate( slope[0], slope[1], slope[2], slope[3], 6.0000000000393072e-5, 6.0000050000143072e-5, 5.4614353597968027e-5, 5.4614402618738265e-5 );

    exit( errCount );
}
/*
************************************************************
*/
static int integrate( double x1, double y1, double x2, double y2, double nn, double gn, double ng, double gg ) {

    int errCount = 0;

    if( verbose ) printf( "x1 = %.17g\ny1 = %.17g\nx2 = %.17g\ny2 = %.17g\n", x1, y1, x2, y2 );
    errCount += integrate2( "nn", ptwXY_interpolationLinLin, x1, y1, x2, y2, nn );
    errCount += integrate2( "gn", ptwXY_interpolationLogLin, x1, y1, x2, y2, gn );
    errCount += integrate2( "ng", ptwXY_interpolationLinLog, x1, y1, x2, y2, ng );
    errCount += integrate2( "gg", ptwXY_interpolationLogLog, x1, y1, x2, y2, gg );
    return( errCount );
}
/*
************************************************************
*/
static int integrate2( char *label, ptwXY_interpolation interpolation, double x1, double y1, double x2, double y2, double answer ) {

    nfu_status status;
    double integral, s, d, r;
    char str[128], buffer[64], *e;

    if( ( status = ptwXY_f_integrate( interpolation, x1, y1, x2, y2, &integral ) ) != nfu_Okay ) 
        nfu_printErrorMsg( "ERROR %s: ptwXY_f_integrate of '%s' with status %d: %s", __FILE__, label, status, nfu_statusMessage( status ) );
    s = 0.5 * ( fabs( answer ) + fabs( integral ) );
    d = answer - integral;
    r = d;
    if( s != 0 ) r /= s;
    if( fabs( r ) > 1e-12 ) printf( "ERROR %s: %s compare, %e %e %e %e %e\n", __FILE__, label, answer, integral, s, d, r );


    sprintf( str, "%.17g", integral );
    e = strchr( str, 'e' );
    if( e != NULL ) {
        sprintf( buffer, " * 10^(%s)", &(e[1]) );
        strcpy( e, buffer );

    }
    if( verbose ) printf( "Integrate[ %s[x], { x, x1, x2 } ] / ( %s ) - 1\n", label, str );

    return( 0 );
}
