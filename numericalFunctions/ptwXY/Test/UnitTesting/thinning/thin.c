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


static int verbose = 0, timing = 0;
static char *fmtXY = "%25.17e %25.17e\n";
static FILE *infoF;

static int checkThinning( statusMessageReporting *smr, const char * const label, ptwXYPoints *ptwXY1, double accuracy, int deleteXYs );
static void compareDoubles( double d1, double d2, double eps, const char * const label );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0, i, n;
    ptwXYPoints *dataXY;
    double x, xMax;
    double zeros[] = { -1.5, 0., -1.49999999, 0., 0.5, 0., 2., 0., 2.2, 1, 2.20000001, 0 };
    double triangle[] = { -1.0, 0., 0., 1., 1.0, 0. };
    int nZeros = sizeof( zeros ) / ( 2 * sizeof( double ) ), nTriangles = sizeof( triangle ) / ( 2 * sizeof( double ) );
    statusMessageReporting smr;

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

    if( ( dataXY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., 1e-3, 10, 10, nZeros, zeros, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    checkThinning( &smr, "zeros", dataXY, 1e-3, 1 );

    if( ( dataXY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 10., 1e-3, 10, 10, nTriangles, triangle, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    checkThinning( &smr, "triangle", dataXY, 1e-3, 0 );
    if( ptwXY_thicken( &smr, dataXY, 100, 1e-5, 1. + 1e-5 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    checkThinning( &smr, "triangle thickened", dataXY, 1e-3, 1 );

    n = 1001;
    if( ( dataXY = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 5, 1e-5, n, 40, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i = 0; i < n; i++ ) {
        x = i * M_PI / ( n - 1 );
        dataXY->points[i].x = x;
        dataXY->points[i].y = sin( x );
    }
    dataXY->length = n;
    checkThinning( &smr, "sin", dataXY, 1e-3, 1 );

    n = 100001;
    if( ( dataXY = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 5, 1e-5, n, 40, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i = 0; i < n; i++ ) {
        x = i * M_PI / ( n - 1 );
        dataXY->points[i].x = x;
        dataXY->points[i].y = sin( x );
    }
    dataXY->length = n;
    checkThinning( &smr, "sin with many points", dataXY, 1e-3, 1 );

    n = 1000001;
    if( ( dataXY = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 5, 1e-5, n, 40, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    xMax = 25 * M_PI;
    for( i = 0; i < n; i++ ) {
        x = i  * xMax / ( n - 1 );
        dataXY->points[i].x = x;
        dataXY->points[i].y = exp( 2. * ( x - 0.2 * xMax ) / xMax ) * sin( x );
    }
    dataXY->length = n;
    checkThinning( &smr, "exp * sin with many points", dataXY, 1e-2, 1 );

    exit( errCount );
}
/*
************************************************************
*/
static int checkThinning( statusMessageReporting *smr, const char * const label, ptwXYPoints *ptwXY1, double accuracy, int deleteXYs ) {

    int errCount = 0;
    ptwXYPoints *thinned;
    clock_t time0;
    int64_t i;
    double x, y;

    if( verbose ) fprintf( infoF, "# label = %s\n# accuracy = %e\n", label, accuracy );
    printIfVerbose( ptwXY1 );
    time0 = clock( );
    if( ( thinned = ptwXY_thin( smr, ptwXY1, accuracy ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( timing ) printf( "# %s: time = %.4f sec", label, ( clock( ) - time0 ) / ( (double) CLOCKS_PER_SEC ) );
    if( timing ) printf( "  length = %d\n", (int) thinned->length );
    printIfVerbose( thinned );
    for( i = 0; i < ptwXY1->length; i++ ) {     /* This logic requires that ptwXY_thin coalesed ptwXY1. */
        x = ptwXY1->points[i].x;
        if( ptwXY_getValueAtX( smr, thinned, x, &y ) != nfu_Okay )
            nfut_printSMRErrorExit2p( smr, "Via." );
        compareDoubles( y, ptwXY1->points[i].y, accuracy, label );
    }
    ptwXY_free( thinned );

    if( deleteXYs ) ptwXY_free( ptwXY1 );
    return( errCount );
}
/*
************************************************************
*/
static void compareDoubles( double d1, double d2, double eps, const char * const label ) {

    double s, d, r;

    s = 0.5 * ( fabs( d1 ) + fabs( d2 ) );
    d = d2 - d1;
    r = d;
    if( s != 0 ) r /= s;
    if( fabs( r ) > eps ) fprintf( infoF, "ERROR %s: %s compare, %e %e %e %e %e\n", __FILE__, label, d1, d2, s, d, r );
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
