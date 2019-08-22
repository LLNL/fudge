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
#include <float.h>

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>


static int verbose = 0, timing = 0;
static char *fmtXY = "%25.17e %25.17e\n";
static FILE *infoF;
static double epsilon = 1e-10;

static int checkThinning( statusMessageReporting *smr, const char * const label, double epsilon, int n1, double *input, int n2, double *output );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0;
    statusMessageReporting smr;
#include "thinDomain.h"

    smr_initialize( &smr, smr_status_Ok );

    infoF = stdout;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else if( strcmp( "-t", argv[iarg] ) == 0 ) {
            timing = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    checkThinning( &smr, "test 1.1", epsilon,      n_test1_input, test1_input, n_test1_1output, test1_1output );
    checkThinning( &smr, "test 1.2", epsilon / 3., n_test1_input, test1_input, n_test1_2output, test1_2output );
    checkThinning( &smr, "test 1.3", 3 * epsilon,  n_test1_input, test1_input, n_test1_3output, test1_3output );
/*
*   Test for point close to x[0] but far from x[2].
*/
    checkThinning( &smr, "test 2.1", epsilon,      n_test2_input, test2_input, n_test2_1output, test2_1output );
    checkThinning( &smr, "test 2.2", epsilon / 3., n_test2_input, test2_input, n_test2_2output, test2_2output );
    checkThinning( &smr, "test 2.3", 3 * epsilon,  n_test2_input, test2_input, n_test2_3output, test2_3output );
/*
*   Test for point close to x[2] but far from x[0].
*/
    checkThinning( &smr, "test 3.1", epsilon,      n_test3_input, test3_input, n_test3_1output, test3_1output );
    checkThinning( &smr, "test 3.2", epsilon / 3., n_test3_input, test3_input, n_test3_2output, test3_2output );
    checkThinning( &smr, "test 3.3", 3 * epsilon,  n_test3_input, test3_input, n_test3_3output, test3_3output );

    exit( errCount );
}
/*
************************************************************
*/
static int checkThinning( statusMessageReporting *smr, const char * const label, double epsilon, int n1, double *input, int n2, double *output ) {

    int errCount = 0;
    ptwXYPoint *points1, *points2;
    ptwXYPoints *inputXY, *outputXY, *thinned;
    clock_t time0;
    int64_t i1;
    double x1, x2;

    if( ( inputXY = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, 10., 1e-3, n1, 10, n1, input, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( smr, "Via." );
    if( ( outputXY = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, 10., 1e-3, n2, 10, n2, output, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( smr, "Via." );

    if( verbose ) fprintf( infoF, "# label = %s\n# epsilon = %e\n", label, epsilon );
    printIfVerbose( inputXY );
    if( strcmp( "test 3.3", label ) ) {
        
    }
    printIfVerbose( outputXY );
    time0 = clock( );
    if( ( thinned = ptwXY_thinDomain( smr, inputXY, epsilon ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( timing ) printf( "# %s: time = %.4f sec", label, ( clock( ) - time0 ) / ( (double) CLOCKS_PER_SEC ) );
    if( timing ) printf( "  length = %d\n", (int) thinned->length );
    printIfVerbose( thinned );

    points1 = thinned->points;
    x1 = points1->x;
    for( i1 = 1, ++points1; i1 < thinned->length; ++i1, ++points1 ) {
        x2 = points1->x;
        if( ( x2 - x1 ) < ( 0.5 * epsilon * ( fabs( x1 ) + fabs( x2 ) ) ) ) break;
        x1 = x2;
    }
    if( thinned->length != outputXY->length ) {
        fprintf( infoF, "ERROR %s: %s compare, length not equal: %d vs. %d\n", __FILE__, label, (int) thinned->length, (int) outputXY->length ); }
    else {
        points1 = thinned->points;
        points2 = outputXY->points;
        for( i1 = 0; i1 < thinned->length; ++i1, ++points1, ++points2 ) {
            x1 = points1->x;
            x2 = points2->x;
            if( fabs( x2 - x1 ) > ( 5 * DBL_EPSILON *( fabs( x1 ) + fabs( x2 ) ) ) ) {
                fprintf( infoF, "ERROR %s: %s compare, answer and results differ at index %d: %.17e vs. %.17e\n", 
                        __FILE__, label, (int) i1, x1, x2 );
                break;
            }
        }
    }

    ptwXY_free( thinned );
    ptwXY_free( inputXY );
    ptwXY_free( outputXY );
    return( errCount );
}
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    fprintf( infoF, "# length = %d\n", (int) ptwXY_length( NULL, data ) );
    ptwXY_simpleWrite( data, infoF, fmtXY );
    fprintf( infoF, "\n\n" );
}
