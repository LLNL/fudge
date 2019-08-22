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
#include <time.h>

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>

#define M_SQRT2PI ( 2. * M_SQRT2 / M_2_SQRTPI )

static int verbose = 0;
static char *fmtXY = "%25.17e %25.17e\n";
static FILE *infoF;

static int invert( statusMessageReporting *smr, ptwXYPoints *ptwXY1 );
static void printIfVerbose( ptwXYPoints *data );
static int compare( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errCount = 0;
    ptwXYPoints *ptwXY1, *ptwXY2;
    double accuracy = 1e-3;
    double twoPoints[4] = { 1., 1., 10., 2. };
    int nTwoPoints = sizeof( twoPoints ) / sizeof( twoPoints[0] ) / 2;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    infoF = stdout;

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

    if( ( ptwXY1 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., accuracy, 100, 10, nTwoPoints, twoPoints, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += invert( &smr, ptwXY1 );

    if( ( ptwXY1 = ptwXY_create( &smr, ptwXY_interpolationLogLin, NULL, 10., accuracy, 100, 10, nTwoPoints, twoPoints, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += invert( &smr, ptwXY1 );

    if( ( ptwXY2 = ptwXY_create( &smr, ptwXY_interpolationLogLin, NULL, 10., accuracy, 100, 10, nTwoPoints, twoPoints, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( ptwXY1 = ptwXY_toOtherInterpolation( &smr, ptwXY2, ptwXY_interpolationLinLin, accuracy ) ) == NULL )
        nfut_printSMRError2p( &smr, "Via." );
    ptwXY_free( ptwXY2 );
    errCount += invert( &smr, ptwXY1 );

    if( ( ptwXY1 = ptwXY_create( &smr, ptwXY_interpolationLinLog, NULL, 10., accuracy, 100, 10, nTwoPoints, twoPoints, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += invert( &smr, ptwXY1 );

    if( ( ptwXY2 = ptwXY_create( &smr, ptwXY_interpolationLinLog, NULL, 10., accuracy, 100, 10, nTwoPoints, twoPoints, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( ptwXY1 = ptwXY_toOtherInterpolation( &smr, ptwXY2, ptwXY_interpolationLinLin, accuracy ) ) == NULL )
        nfut_printSMRError2p( &smr, "Via." );
    ptwXY_free( ptwXY2 );
    errCount += invert( &smr, ptwXY1 );

    if( ( ptwXY1 = ptwXY_create( &smr, ptwXY_interpolationLogLog, NULL, 10., accuracy, 100, 10, nTwoPoints, twoPoints, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += invert( &smr, ptwXY1 );

    if( ( ptwXY2 = ptwXY_create( &smr, ptwXY_interpolationLogLog, NULL, 10., accuracy, 100, 10, nTwoPoints, twoPoints, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( ptwXY1 = ptwXY_toOtherInterpolation( &smr, ptwXY2, ptwXY_interpolationLinLin, accuracy ) ) == NULL )
        nfut_printSMRError2p( &smr, "Via." );
    ptwXY_free( ptwXY2 );
    errCount += invert( &smr, ptwXY1 );

    exit( errCount );
}
/*
************************************************************
*/
static int invert( statusMessageReporting *smr, ptwXYPoints *ptwXY1 ) {

    int errCount;
    ptwXYPoints *ptwXY2, *ptwXY3;

    printIfVerbose( ptwXY1 );
    if( ( ptwXY2 = ptwXY_inverse( smr, ptwXY1 ) ) == NULL ) 
        nfut_printSMRError2( smr, __FILE__, __LINE__, __func__, "Via." );
    printIfVerbose( ptwXY2 );
    if( ( ptwXY3 = ptwXY_inverse( smr, ptwXY2 ) ) == NULL ) 
        nfut_printSMRError2( smr, __FILE__, __LINE__, __func__, "Via." );
    printIfVerbose( ptwXY3 );
    errCount = compare( ptwXY1, ptwXY3 );
    ptwXY_free( ptwXY1 );
    ptwXY_free( ptwXY2 );
    ptwXY_free( ptwXY3 );
    return( errCount );
}
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    fprintf( infoF, "# length = %d\n", (int) ptwXY_length( NULL, data ) );
    fprintf( infoF, "# interpolation = %d\n", ptwXY_getInterpolation( data ) );
    fprintf( infoF, "# interpolation string = %s\n", ptwXY_getInterpolationString( data ) );
    ptwXY_simpleWrite( data, infoF, fmtXY );
    fprintf( infoF, "\n\n" );
}
/*
************************************************************
*/
static int compare( ptwXYPoints *ptwXY1, ptwXYPoints *ptwXY2 ) {

    int64_t i1;
    ptwXYPoint *point1, *point2;

    if( ptwXY_length( NULL, ptwXY1 ) != ptwXY_length( NULL, ptwXY2 ) ) return( 1 );

    for( i1 = 0; i1 < ptwXY_length( NULL, ptwXY1 ); ++i1 ) {
        point1 = ptwXY_getPointAtIndex_Unsafely( ptwXY1, i1 );
        point2 = ptwXY_getPointAtIndex_Unsafely( ptwXY2, i1 );
        if( point1->x != point2->x ) return( 1 );
        if( point1->y != point2->y ) return( 1 );
    }
    return( 0 );
}
