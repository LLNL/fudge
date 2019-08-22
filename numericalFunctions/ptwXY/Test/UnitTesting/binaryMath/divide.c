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

struct XYData {
    int points;
    double *XYs;
};

static int verbose = 0;
static double biSectionMax = 5., accuracy = 1e-3;
static char *fmtXY = "%18.11e %18.11e\n";

static int divide( statusMessageReporting *smr, int n1, double *XY1, int n2, double *XY2 );
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

    for( i = 0; i < nDatas; i++ ) {
        for( j = 0; j < nDatas; j++ ) {
            errCount += divide( &smr, XYDatas[i].points, XYDatas[i].XYs, XYDatas[j].points, XYDatas[j].XYs );
        }
    }

    exit( errCount );
}
/*
************************************************************
*/
static int divide( statusMessageReporting *smr, int n1, double *XY1, int n2, double *XY2 ) {

    int errCount = 0;
    ptwXYPoints *ptwXY1, *ptwXY2, *division;

    if( ( ptwXY1 = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, biSectionMax, accuracy, 10, 10, n1, XY1, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( smr, "Via." );
    if( ( ptwXY2 = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, biSectionMax, accuracy, 10, 10, n2, XY2, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( smr, "Via." );

    if( verbose ) {
        printf( "# biSectionMax = %.12e\n", biSectionMax );
        printf( "# accuracy = %.12e\n", accuracy );
        printf( "# length = %d\n", (int) ptwXY_length( smr, ptwXY1 ) );
        ptwXY_simpleWrite( ptwXY1, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( smr, ptwXY2 ) );
        ptwXY_simpleWrite( ptwXY2, stdout, fmtXY );
        printf( "\n\n" );
    }

    if( ( division = ptwXY_div_ptwXY( smr, ptwXY1, ptwXY2, 1 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );

    if( verbose ) {
        printf( "# length = %d\n", (int) ptwXY_length( smr, division ) );
        ptwXY_simpleWrite( division, stdout, fmtXY );
        printf( "\n\n" );
    }

    ptwXY_free( ptwXY1 );
    ptwXY_free( ptwXY2 );
    ptwXY_free( division );

    return( errCount );
}
