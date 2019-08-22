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

#define nSame 6

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";
static double yMin = -0.9, yMax = 0.4;

static int checkClipping( statusMessageReporting *smr, ptwXYPoints *data );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;
    ptwXYPoints *XY, *expXY, *mulXY;
    double x, expXYs[4], domainMin, domainMax;
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

    if( ( XY = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 10, 10, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i = 0; i < nSame; i++ ) {
        if( ptwXY_setValueAtX( &smr, XY, 0.2 * i - .5, yMax + i - .1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkClipping( &smr, XY );
    if( ptwXY_neg( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkClipping( &smr, XY );
    if( ptwXY_neg( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );

    for( ; i < 2 * nSame + 1; i++ ) {
        if( ptwXY_setValueAtX( &smr, XY, 0.2 * i - .5, yMin - i - .1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkClipping( &smr, XY );
    if( ptwXY_neg( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkClipping( &smr, XY );
    if( ptwXY_neg( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );

    if( ptwXY_clear( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );

    for( i = 0; i < 501; i++ ) {
        x = i * M_PI / 50;
        if( ptwXY_setValueAtX( &smr, XY, x, sin( x ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkClipping( &smr, XY );

    if( ptwXY_domainMin( &smr, XY, &domainMin ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_domainMax( &smr, XY, &domainMax ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    expXYs[0] = domainMin;
    expXYs[1] = 0.;
    expXYs[2] = domainMax;
    expXYs[3] = 1.;
    if( ( expXY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 100, 10, 2, expXYs, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    printIfVerbose( expXY );
    if( ptwXY_exp( &smr, expXY, 1. ) != nfu_Okay ) nfut_printSMRError2p( &smr, "Via." );
    printIfVerbose( expXY );
    if( ( mulXY = ptwXY_mul_ptwXY( &smr, XY, expXY ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkClipping( &smr, mulXY );

    ptwXY_free( XY );
    ptwXY_free( expXY );
    ptwXY_free( mulXY );

    exit( errCount );
}
/*
************************************************************
*/
static int checkClipping( statusMessageReporting *smr, ptwXYPoints *data ) {

    int i, errCount = 0;
    ptwXYPoints *clipped, *u;
    ptwXYPoint *point;
    nfu_status status;
    double x1, y1, x2, y2, x, y, s;

    if( verbose ) {
        printf( "# yMin = %.14e\n", yMin );
        printf( "# yMax = %.14e\n", yMax );
    }
    printIfVerbose( data );
    if( ( clipped = ptwXY_clone( smr, data ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_clip( smr, clipped, yMin, yMax ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    printIfVerbose( clipped );

    if( ( u = ptwXY_union( smr, clipped, data, ptwXY_union_fill ) ) == NULL )
        nfut_printSMRErrorExit2p( smr, "Via." );
    for( i = 0; i < u->length; i++ ) {
        point = ptwXY_getPointAtIndex_Unsafely( u, i );
        if( point->y < yMin ) {
            nfu_printMsg( "ERROR %s: at x = %g, point->y = %g < yMin = %g", __FILE__, point->x, point->y, yMin );
            errCount++; }
        else if( point->y > yMax ) {
            nfu_printMsg( "ERROR %s: at x = %g, point->y = %g > yMax = %g", __FILE__, point->x, point->y, yMax );
            errCount++;
        }
        if( ( status = ptwXY_getValueAtX( smr, data, point->x, &y ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
        if( y < yMin ) {
            if( point->y != yMin ) {
                nfu_printMsg( "ERROR %s: at x, y = (%g, %g), point->y = %g != yMin = %g", __FILE__, point->x, y, point->y, yMin );
                errCount++;
            }
        }
        else if( y > yMax ) {
            if( point->y != yMax ) {
                nfu_printMsg( "ERROR %s: at x, y = (%g, %g), point->y = %g != yMax = %g", __FILE__, point->x, y, point->y, yMax );
                errCount++;
            }
        }
    }

    for( i = 0; i < data->length; i++ ) {       /* test insures that points are added at yMin and yMax crossings. */
        point = ptwXY_getPointAtIndex_Unsafely( data, i );
        x2 = point->x;
        y2 = point->y;
        if( i > 0 ) {
            if( ( ( y1 - yMin ) * ( y2 - yMin ) ) < 0 ) {
                x = ( x1 * ( y2 - yMin ) + x2 * ( yMin - y1 ) ) / ( y2 - y1 );
                if( ( status = ptwXY_getValueAtX( smr, clipped, x, &y ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
                s = 0.5 * ( fabs( y ) + fabs( yMin ) );
                if( fabs( y - yMin ) > 1e-14 * s ) {
                    nfu_printMsg( "ERROR %s: at i = %d, at x = %g, y = %g != yMin = %g", __FILE__, i, x, y, yMin );
                    errCount++;
                }
            }
            if( ( ( y1 - yMax ) * ( y2 - yMax ) ) < 0 ) {
                x = ( x1 * ( y2 - yMax ) + x2 * ( yMax - y1 ) ) / ( y2 - y1 );
                if( ( status = ptwXY_getValueAtX( smr, clipped, x, &y ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
                s = 0.5 * ( fabs( y ) + fabs( yMax ) );
                if( fabs( y - yMax ) > 1e-14 * s ) {
                    nfu_printMsg( "ERROR %s: at i = %d, at x = %g, y = %g != yMax = %g", __FILE__, i, x, y, yMax );
                    errCount++;
                }
            }
        }
        x1 = x2;
        y1 = y2;
    }
    ptwXY_free( u );
    ptwXY_free( clipped );

    return( errCount );
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
