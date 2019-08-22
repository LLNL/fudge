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

/*
*   2011 July 11: A bug was discovered in the ptwXY_union in the fill section. This is a test of that the bug is fixed.
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>

struct XYData {
    int points;
    double *XYs;
};

static int verbose = 0;
static double biSectionMax = 5., accuracy = 1e-3;
static char *fmtXY = "%18.11e %18.11e\n";

/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0;
    int64_t i;
    double x, y;
    double XY0[] = { 1.423095e+07, 0.0, 1.440000e+07, 1.745688e-05, 1.460000e+07, 2.285425e-04, 1.480000e+07, 7.085735e-04 };
    double XY1[] = { 1.423095e+07, 0.0, 1.423100e+07, 0.000000e+00, 1.440000e+07, 0.000000e+00, 1.460000e+07, 0.000000e+00, 1.480000e+07, 6.155894e-05 };
    int nXY0 = sizeof( XY0 ) / ( 2 * sizeof( double ) ), nXY1 = sizeof( XY1 ) / ( 2 * sizeof( double ) );
    ptwXYPoints *ptwXY0, *ptwXY1, *add1, *add2, *diff;
    nfu_status status;
                   
    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if(echo ) printf( "%s\n", __FILE__ );

    if( ( ptwXY0 = ptwXY_create( ptwXY_interpolationLinLin, biSectionMax, accuracy, 10, 10, nXY0, XY0, &status, 0 ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: ptwXY0 creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( ptwXY1 = ptwXY_create( ptwXY_interpolationLinLin, biSectionMax, accuracy, 10, 10, nXY1, XY1, &status, 0 ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: ptwXY1 creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    if( ( add1 = ptwXY_add_ptwXY( ptwXY0, ptwXY1, &status ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: add1 status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    if( ( add2 = ptwXY_add_ptwXY( ptwXY1, ptwXY0, &status ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: add2 status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    if( ( diff = ptwXY_sub_ptwXY( add1, add2, &status ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: diff status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    if( verbose ) {
        printf( "# biSectionMax = %.12e\n", biSectionMax );
        printf( "# accuracy = %.12e\n", accuracy );
        printf( "# length = %d\n", (int) ptwXY_length( ptwXY0 ) );
        ptwXY_simpleWrite( ptwXY0, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( ptwXY1 ) );
        ptwXY_simpleWrite( ptwXY1, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( add1 ) );
        ptwXY_simpleWrite( add1, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( add2 ) );
        ptwXY_simpleWrite( add2, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( diff ) );
        ptwXY_simpleWrite( diff, stdout, fmtXY );
        printf( "\n\n" );
    }

    for( i = 0; i < diff->length; i++ ) {
        if( ( status = ptwXY_getXYPairAtIndex( diff, i, &x, &y ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_getXYPairAtIndex status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
        if( y != 0 ) {
            nfu_printMsg( "ERROR %s: value at index %d = %e is not 0.", __FILE__, (int) i, y );
            errCount++;
        }
    }

    ptwXY_free( ptwXY0 );
    ptwXY_free( ptwXY1 );
    ptwXY_free( add1 );
    ptwXY_free( add2 );
    ptwXY_free( diff );

    exit( errCount );
}
