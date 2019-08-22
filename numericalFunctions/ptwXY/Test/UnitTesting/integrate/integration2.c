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
static char *fmtXY = "%19.12e %19.12e\n";

static int integrate( int n1, double *xys1, double xMin, double xMax );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0;
    double XYs1[8] = { 0.00000000e+00, 0.00000000e+00, 1.95312500e+04, 9.59428000e-08, 2.14843750e+05, 2.76211000e-07, 2.00000000e+07, 1.58192490e-12 };

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

    errCount += integrate( 4, XYs1, -1e-5, 1e+3 );
    errCount += integrate( 4, XYs1, 1e-5, 1e+3 );
    errCount += integrate( 4, XYs1, 1e-5, XYs1[2] );
    errCount += integrate( 4, XYs1, 1e-5, 0.5 * ( XYs1[2] + XYs1[4] ) );
    errCount += integrate( 4, XYs1, 0.5 * ( XYs1[2] + XYs1[4] ), 0.6 * ( XYs1[2] + XYs1[4] ) );
    errCount += integrate( 4, XYs1, 0.5 * ( XYs1[2] + XYs1[4] ), 0.6 * ( XYs1[4] + XYs1[6] ) );

    exit( errCount );
}
/*
************************************************************
*/
static int integrate( int n1, double *xys1, double xMin, double xMax ) {

    int errCount = 0;
    nfu_status status;
    double integral;
    ptwXYPoints *ptwXY1;

    if( ( ptwXY1 = ptwXY_create( ptwXY_interpolationLinLin, 8, 1e-5, 10, 10, n1, xys1, &status, 0 ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: from ptwXY_create, status = %d, '%s'", __FILE__, status, nfu_statusMessage( status ) );

    if( verbose ) {
        printf( "# xMin = %.8e\n", xMin );
        printf( "# xMax = %.8e\n", xMax );
    }
    printIfVerbose( ptwXY1 );

    integral = ptwXY_integrate( ptwXY1, xMin, xMax, &status );
    if( status != nfu_Okay ) nfu_printErrorMsg( "ERROR %s: from ptwXY_integrate, status = %d, '%s'", __FILE__, status, nfu_statusMessage( status ) );
    if( verbose ) printf( "# integral = %.8e\n\n\n", integral );
    return( errCount );
}
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    printf( "# length = %d\n", (int) data->length );
    printf( "# interpolation = %d\n", (int) data->interpolation );
    ptwXY_simpleWrite( data, stdout, fmtXY );
    printf( "\n\n" );
}
