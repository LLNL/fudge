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

struct XYData {
    int points;
    double *XYs;
};

static int verbose = 0;
static double biSectionMax = 5., accuracy = 1e-3;
static char *fmtXY = "%18.11e %18.11e\n";

static int divide( int n1, double *XY1, int n2, double *XY2 );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0, i, j, nDatas;
    double XY1[6] = { 2., 0., 2.2, 3., 2.6, 0. }, XY2[4] = { 2., 0., 2.6, 3. }, XY3[4] = { 2., 2., 2.6, 0. }, XY4[4] = { 2., 2., 2.6, 1. };
    double XY5[10] = { 2., 0., 2.2, 3., 2.3, 0., 2.5, 0., 2.6, 0. }, XY6[10] = { 2., 0., 2.2, 0., 2.3, 0., 2.5, 3., 2.6, 0. };
    double XY7[4] = { 2., 1.0, 2.6, -1.0 }, XY8[6] = { 2., 1., 2.2, 0., 2.6, 2. };
    struct XYData XYDatas[] = { { sizeof( XY1 ) / ( 2 * sizeof( double ) ), XY1 }, { sizeof( XY2 ) / ( 2 * sizeof( double ) ), XY2 },
                                { sizeof( XY3 ) / ( 2 * sizeof( double ) ), XY3 }, { sizeof( XY4 ) / ( 2 * sizeof( double ) ), XY4 },
                                { sizeof( XY5 ) / ( 2 * sizeof( double ) ), XY5 }, { sizeof( XY6 ) / ( 2 * sizeof( double ) ), XY6 },
                                { sizeof( XY7 ) / ( 2 * sizeof( double ) ), XY7 }, { sizeof( XY8 ) / ( 2 * sizeof( double ) ), XY8 } };

    nDatas = sizeof( XYDatas ) / sizeof( struct XYData );

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

    for( i = 0; i < nDatas; i++ ) {
        for( j = 0; j < nDatas; j++ ) {
            errCount += divide( XYDatas[i].points, XYDatas[i].XYs, XYDatas[j].points, XYDatas[j].XYs );
        }
    }

    exit( errCount );
}
/*
************************************************************
*/
static int divide( int n1, double *XY1, int n2, double *XY2 ) {

    int errCount = 0;
    ptwXYPoints *ptwXY1, *ptwXY2, *division;
    nfu_status status;

    if( ( ptwXY1 = ptwXY_create( ptwXY_interpolationLinLin, biSectionMax, accuracy, 10, 10, n1, XY1, &status, 0 ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: ptwXY1 creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( ptwXY2 = ptwXY_create( ptwXY_interpolationLinLin, biSectionMax, accuracy, 10, 10, n2, XY2, &status, 0 ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: ptwXY2 creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    if( verbose ) {
        printf( "# biSectionMax = %.12e\n", biSectionMax );
        printf( "# accuracy = %.12e\n", accuracy );
        printf( "# length = %d\n", (int) ptwXY_length( ptwXY1 ) );
        ptwXY_simpleWrite( ptwXY1, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( ptwXY2 ) );
        ptwXY_simpleWrite( ptwXY2, stdout, fmtXY );
        printf( "\n\n" );
    }

    if( ( division = ptwXY_div_ptwXY( ptwXY1, ptwXY2, &status, 1 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: ptwXY2 division status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    if( verbose ) {
        printf( "# length = %d\n", (int) ptwXY_length( division ) );
        ptwXY_simpleWrite( division, stdout, fmtXY );
        printf( "\n\n" );
    }

    ptwXY_free( ptwXY1 );
    ptwXY_free( ptwXY2 );
    ptwXY_free( division );

    return( errCount );
}
