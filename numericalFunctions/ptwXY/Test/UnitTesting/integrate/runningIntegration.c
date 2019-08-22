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
    if( ( ptwXY = ptwXY_new( ptwXY_interpolationLinLin, NULL, 6, 1e-3, npdf, 10, &status, 0 ) ) == NULL ) 
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
