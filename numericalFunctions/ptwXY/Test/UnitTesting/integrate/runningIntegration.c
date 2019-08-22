/*
# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>

static int verbose = 0;

static int integrate( statusMessageReporting *smr, char *label, int npdf, double *xy_runningSum );
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
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

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
    errCount = integrate( &smr, "pdf1", npdf, pdf1 );

    exit( errCount );
}
/*
************************************************************
*/
static int integrate( statusMessageReporting *smr, char *label, int npdf, double *xy_runningSum ) {

    int i1, errCount = 0;
    double d, r, sum1, sum2;
    ptwXYPoints *ptwXY;
    ptwXPoints *sums;

    if( verbose ) {
        printf( "# length = %d\n", npdf );
        for( i1 = 0; i1 < 3 * npdf; i1 += 3 ) printf( " %.12e %.12e %.12e\n", xy_runningSum[i1], xy_runningSum[i1+1], xy_runningSum[i1+2] );
    }
    if( ( ptwXY = ptwXY_new( smr, ptwXY_interpolationLinLin, NULL, 6, 1e-3, npdf, 10, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( smr, "Via." );
    for( i1 = 0; i1 < 3 * npdf; i1 += 3 ) {
        if( ptwXY_setValueAtX( smr, ptwXY, xy_runningSum[i1], xy_runningSum[i1+1] ) != nfu_Okay )
            nfut_printSMRErrorExit2p( smr, "Via." );
    }
    if( ( sums = ptwXY_runningIntegral( smr, ptwXY ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
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
