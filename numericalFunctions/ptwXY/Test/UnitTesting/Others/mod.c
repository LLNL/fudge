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

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int checkMod( ptwXYPoints *data, double m, double *Ys );
static void compareDoubles( double x, double m, double d1, double d2, double eps, const char * const funcName );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;
    nfu_status status;
    ptwXYPoints *XYs;
    double Xs[] = {  -8.3381123168917412e+01, -7.8550646919038201e+01, -7.8052706271881675e+01, -7.2676812226386218e+01, -6.7456460315920452e+01, -5.9577593804676845e+01, 
                     -3.9429862990821299e+01, -3.6983850296024201e+01, -1.9070786454251021e+01, -1.4917762730536310e+01,  6.2690762408937957e+00,  1.4018281560416000e+01,  
                      1.8343344831976211e+01,  2.8638310700345272e+01,  3.5726737104437206e+01,  4.9705617757084099e+01,  6.3540253205248973e+01,  6.8112449656470716e+01, 
                      8.6914952896485005e+01,  8.7407476974652553e+01 };
    double YPs[] = {  1.4418784780070020e+00,  3.1307620742964204e+00,  4.8711006786315281e-01,  2.7214114597688166e+00,  1.6585780630549962e+00,  1.1266661352922469e-01,
                      1.4108415058460118e+00,  7.1526154705331635e-01,  2.9203621208775310e+00,  7.9020053741265528e-01,  3.1274835873040026e+00,  1.4519109460568274e+00,  
                      2.6353815640272451e+00,  3.6397681803713411e-01,  1.1692179149494812e+00,  2.5817279532372019e+00,  7.0840013345311093e-01,  2.1390039310850604e+00,
                      2.0919512495605908e+00,  2.5844753277281391e+00 };
    double YNs[] = { -1.699714175582791e+00,  -1.083057929337272e-02,  -2.654482585726640e+00,  -4.201811938209765e-01,  -1.483014590534797e+00,  -3.028926040060568e+00, 
                     -1.730751147743781e+00,  -2.426331106536477e+00,  -2.212305327122621e-01,  -2.351392116177138e+00,  -1.410906628579056e-02,  -1.689681707532966e+00,  
                     -5.062110895625480e-01,  -2.777615835552659e+00,  -1.972374738640312e+00,  -5.598647003525912e-01,  -2.433192520136682e+00,  -1.002588722504733e+00, 
                     -1.049641404029202e+00,  -5.571173258616540e-01 };
    int nYs = sizeof( Xs ) / sizeof( double );

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

    if( ( XYs = ptwXY_new( ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 10, 10, &status, 0 ) ) == NULL ) 
            nfu_printErrorMsg( "ERROR %s: XYs creation, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    for( i = 0; i < nYs; i++ ) {
        if( ( status = ptwXY_setValueAtX( XYs, Xs[i], Xs[i] ) ) != nfu_Okay )
                    nfu_printErrorMsg( "ERROR %s: ptwXY_setValueAtX, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    }
    errCount += checkMod( XYs,  M_PI, YPs );
    errCount += checkMod( XYs, -M_PI, YNs );

    ptwXY_free( XYs );

    exit( errCount );
}
/*
************************************************************
*/
static int checkMod( ptwXYPoints *data, double m, double *Ys ) {

    int64_t i, errCount = 0;
    nfu_status status = nfu_Okay;
    ptwXYPoints *moded;
    double x, y, yc;

    if( verbose ) printf( "# mod = %.15e\n", m );
    printIfVerbose( data );
    if( ( moded = ptwXY_clone( data, &status ) ) == NULL )
            nfu_printErrorMsg( "ERROR %s: data clone, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_mod( moded, m, 0 ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_mod, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    printIfVerbose( moded );

    if( ( status = ptwXY_simpleCoalescePoints( moded ) ) != nfu_Okay ) 
            nfu_printErrorMsg( "ERROR %s: coalescing, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    for( i = 0; i < moded->length; i++ ) {
        x = moded->points[i].x;
        y = moded->points[i].y;
        yc = fmod( x, m );
        compareDoubles( x, m, y, yc, 1e-12, "checkMod" );
    }
    ptwXY_free( moded );

    if( ( moded = ptwXY_clone( data, &status ) ) == NULL )
            nfu_printErrorMsg( "ERROR %s: data clone, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_mod( moded, m, 1 ) ) != nfu_Okay )
            nfu_printErrorMsg( "ERROR %s: ptwXY_mod, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( status = ptwXY_simpleCoalescePoints( moded ) ) != nfu_Okay ) 
            nfu_printErrorMsg( "ERROR %s: coalescing, status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    for( i = 0; i < moded->length; i++ ) {
        x = moded->points[i].x;
        y = moded->points[i].y;
        yc = Ys[i];
        compareDoubles( x, m, y, yc, 1e-12, "checkMod python" );
    }
    ptwXY_free( moded );

    return( errCount );
}
/*
************************************************************
*/
static void compareDoubles( double x, double m, double d1, double d2, double eps, const char * const funcName ) {

    double s, d, r;

    s = 0.5 * ( fabs( d1 ) + fabs( d2 ) );
    d = d2 - d1;
    r = d;
    if( s != 0 ) r /= s;
    if( fabs( r ) > eps ) fprintf( stdout, "ERROR %s: %s compare x mod m = %e mod %e, %e %e %e %e %e\n", __FILE__, funcName, x, m, d1, d2, s, d, r );
}
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    printf( "# length = %d\n", (int) data->length );
    ptwXY_simpleWrite( data, stdout, fmtXY );
    printf( "\n\n" );
}
