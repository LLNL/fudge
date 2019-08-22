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

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_utilities.h>

static int verbose = 0;
static char *fmtXY = "%19.12e %19.12e\n";

static int clone( statusMessageReporting *smr, int n1, int allocatedSize, int overflowAllocatedSize );
static int compareXYs( statusMessageReporting *smr, ptwXYPoints *XY1, ptwXYPoints *XY2 );
static int compareValues( int64_t i, double x1, double y1, double x2, double y2 );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errCount = 0;
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

    errCount += clone( &smr, 5, 10, 10 );
    errCount += clone( &smr, 16, 10, 10 );
    errCount += clone( &smr, 23, 10, 10 );
    errCount += clone( &smr, 56, 10, 10 );
    errCount += clone( &smr, 5, 100, 10 );
    errCount += clone( &smr, 16, 100, 10 );
    errCount += clone( &smr, 123, 100, 10 );
    errCount += clone( &smr, 556, 100, 10 );
    errCount += clone( &smr, 5, 100, 50 );
    errCount += clone( &smr, 16, 100, 50 );
    errCount += clone( &smr, 123, 100, 50 );
    errCount += clone( &smr, 556, 100, 50 );
    return( errCount );
}
/*
************************************************************
*/
static int clone( statusMessageReporting *smr, int n1, int allocatedSize, int overflowAllocatedSize ) {

    int i1, errCount = 0, nonOverflowLength1, nonOverflowLength2;
    double x1;
    ptwXYPoints *XYSrc, *XYClone, *XYClone2;

    XYSrc = ptwXY_new( smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, allocatedSize, overflowAllocatedSize, 0 );
    if( XYSrc == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    for( i1 = 0; i1 < n1; ++i1 ) {
        x1 = 1 - 2 * drand48( );
        if( ptwXY_setValueAtX( smr, XYSrc, x1, x1 * x1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    }

    nonOverflowLength1 = (int) ptwXY_getNonOverflowLength( smr, XYSrc );
    if( verbose ) printf( "nonOverflowLength1 = %d\n", nonOverflowLength1 );

    if( ( XYClone2 = ptwXY_clone2( smr, XYSrc ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    nonOverflowLength2 = (int) ptwXY_getNonOverflowLength( smr, XYSrc );
    if( nonOverflowLength1 != nonOverflowLength2 ) ++errCount;

    if( ( XYClone = ptwXY_clone( smr, XYSrc ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );

    printIfVerbose( XYSrc );
    printIfVerbose( XYClone );
    printIfVerbose( XYClone2 );

    if( ptwXY_simpleCoalescePoints( smr, XYSrc ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_simpleCoalescePoints( smr, XYClone ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );

    errCount += compareXYs( smr, XYSrc, XYClone );
    errCount += compareXYs( smr, XYSrc, XYClone2 );

    ptwXY_free( XYSrc );
    ptwXY_free( XYClone );
    ptwXY_free( XYClone2 );
    return( errCount );
}
/*
************************************************************
*/
static int compareXYs( statusMessageReporting *smr, ptwXYPoints *XY1, ptwXYPoints *XY2 ) {

    int errCount = 0;
    int64_t i, n = ptwXY_length( NULL, XY1 );
    double x1, y1, x2, y2;

    if( n != ptwXY_length( NULL, XY2 ) ) {
        errCount++;
        nfu_printMsg( "ERROR %s: compareXYs, len( XY1 ) = %d != len( XY2 ) = %d", __FILE__, (int) n, (int) ptwXY_length( NULL, XY2 ) ); }
    else {
        for( i = 0; i < n; i++ ) {
            if( ptwXY_getXYPairAtIndex( smr, XY1, i, &x1, &y1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
            if( ptwXY_getXYPairAtIndex( smr, XY2, i, &x2, &y2 ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
            errCount += compareValues( i, x1, y1, x2, y2 );
        }
    }

    return( errCount );
}
/*
************************************************************
*/
static int compareValues( int64_t i, double x1, double y1, double x2, double y2 ) {

    if( ( x1 != x2 ) || ( y1 != y2 ) ) {
        nfu_printMsg( "ERROR %s: at index %3d ( x1 = %.17e != x2 = %.17e ) or ( y1 = %.17e != y2 = %.17e )", __FILE__, (int) i, x1, x2, y1, y2 );
        return( 1 );
    }
    return( 0 );
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
