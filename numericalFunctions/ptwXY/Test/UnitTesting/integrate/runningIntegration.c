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

static int integrate( char *label, int npdf, double *xy_runningSum );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0, npdf;
    double pdf1[] = { 
        -1.00000000e+00, 1.12243000e+00, 0.00000000e+00,
        -9.00000000e-01, 9.09141000e-01, 1.01578550e-01,
        -8.00000000e-01, 6.57468000e-01, 1.79909000e-01,
        -7.00000000e-01, 4.76759000e-01, 2.36620350e-01,
        -6.00000000e-01, 3.76110000e-01, 2.79263800e-01,
        -5.00000000e-01, 3.22837000e-01, 3.14211150e-01,
        -4.00000000e-01, 2.93852000e-01, 3.45045600e-01,
        -3.00000000e-01, 2.72462000e-01, 3.73361300e-01,
        -2.00000000e-01, 2.66565000e-01, 4.00312650e-01,
        -1.00000000e-01, 2.72462000e-01, 4.27264000e-01,
        0.00000000e+00, 2.87255000e-01, 4.55249850e-01,
        1.00000000e-01, 3.13941000e-01, 4.85309650e-01,
        2.00000000e-01, 3.52422000e-01, 5.18627800e-01,
        3.00000000e-01, 3.99798000e-01, 5.56238800e-01,
        4.00000000e-01, 4.56070000e-01, 5.99032200e-01,
        5.00000000e-01, 5.18238000e-01, 6.47747600e-01,
        6.00000000e-01, 5.89302000e-01, 7.03124600e-01,
        7.00000000e-01, 6.63365000e-01, 7.65757950e-01,
        8.00000000e-01, 7.40326000e-01, 8.35942500e-01,
        9.00000000e-01, 8.20286000e-01, 9.13973100e-01,
        1.00000000e+00, 9.00245000e-01, 9.99999650e-01 };

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

    npdf = sizeof( pdf1 ) / sizeof( pdf1[0] ) / 3;
    errCount = integrate( "pdf1", npdf, pdf1 );

    exit( errCount );
}
/*
************************************************************
*/
static int integrate( char *label, int npdf, double *xy_runningSum ) {

    int i1, errCount = 0;
    double d, r, sum1, sum2;
    nfu_status status;
    ptwXYPoints *ptwXY;
    ptwXPoints *sums;

    if( verbose ) {
        printf( "# length = %d\n", npdf );
        for( i1 = 0; i1 < 3 * npdf; i1 += 3 ) printf( " %.12e %.12e %.12e\n", xy_runningSum[i1], xy_runningSum[i1+1], xy_runningSum[i1+2] );
    }
    if( ( ptwXY = ptwXY_new( ptwXY_interpolationLinLin, 6, 1e-3, npdf, 10, &status, 0 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s for '%s': ptwXY new, status = %d: %s", __FILE__, label, status, nfu_statusMessage( status ) );
    for( i1 = 0; i1 < 3 * npdf; i1 += 3 ) {
        status = ptwXY_setValueAtX( ptwXY, xy_runningSum[i1], xy_runningSum[i1+1] );
        if( status != nfu_Okay ) nfu_printErrorMsg( "ERROR %s for '%s': ptwXY_setValueAtX at index %d, (x, y) = (%e, %e): %d %s", __FILE__, label, i1, 
            xy_runningSum[i1], xy_runningSum[i1+1], status, nfu_statusMessage( status ) );
    }
    sums = ptwXY_runningIntegral( ptwXY, &status );
    if( status != nfu_Okay ) nfu_printErrorMsg( "ERROR %s for '%s': ptwXY_runningIntegral: %d %s", __FILE__, label, status, nfu_statusMessage( status ) );
    for( i1 = 0; i1 < npdf; i1++ ) {
        sum1 = xy_runningSum[3*i1+2];
        sum2 = ptwX_getPointAtIndex_Unsafely( sums, i1 );
        d = sum2 - sum1;
        r = d;
        if( sum2 != 0 ) r /= sum2;
        if( fabs( r ) > 1e-10 ) {
            errCount++;
            printf( "ERROR %s: %s compare at index %d, %e %e %e %e\n", __FILE__, label, i1, sum1, sum2, d, r );
        }
    }

    return( errCount );
}
