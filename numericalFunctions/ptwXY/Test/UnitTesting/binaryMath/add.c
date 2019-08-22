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

/*
*   2011 July 11: A bug was discovered in the ptwXY_union in the fill section. This is a test of that the bug is fixed.
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
    if(echo ) printf( "%s\n", __FILE__ );

    if( ( ptwXY0 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, biSectionMax, accuracy, 10, 10, nXY0, XY0, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( ptwXY1 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, biSectionMax, accuracy, 10, 10, nXY1, XY1, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );

    if( ( add1 = ptwXY_add_ptwXY( &smr, ptwXY0, ptwXY1 ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( add2 = ptwXY_add_ptwXY( &smr, ptwXY1, ptwXY0 ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( diff = ptwXY_sub_ptwXY( &smr, add1, add2 ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );

    if( verbose ) {
        printf( "# biSectionMax = %.12e\n", biSectionMax );
        printf( "# accuracy = %.12e\n", accuracy );
        printf( "# length = %d\n", (int) ptwXY_length( &smr, ptwXY0 ) );
        ptwXY_simpleWrite( ptwXY0, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( &smr, ptwXY1 ) );
        ptwXY_simpleWrite( ptwXY1, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( &smr, add1 ) );
        ptwXY_simpleWrite( add1, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( &smr, add2 ) );
        ptwXY_simpleWrite( add2, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( &smr, diff ) );
        ptwXY_simpleWrite( diff, stdout, fmtXY );
        printf( "\n\n" );
    }

    for( i = 0; i < diff->length; i++ ) {
        if( ( status = ptwXY_getXYPairAtIndex( &smr, diff, i, &x, &y ) ) != nfu_Okay )
            nfut_printSMRErrorExit2p( &smr, "Via." );
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
