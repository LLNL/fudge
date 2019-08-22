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

static int checkDullEdges( statusMessageReporting *smr, ptwXYPoints *data, double lowerEps, double upperEps );
static int checkDullEdges2( statusMessageReporting *smr, ptwXYPoints *data, double lowerEps, double upperEps, int positiveXOnly );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;
    ptwXYPoints *XY;
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
        if( ptwXY_setValueAtX( &smr, XY, 0.2 * i - .5, 0.7 + i + .1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkDullEdges( &smr, XY, 1e-10, 1e-10 );
    if( ptwXY_neg( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkDullEdges( &smr, XY, 1e-10, 1e-10 );

    if( ptwXY_clear( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i = 0; i < nSame; i++ ) {
        if( ptwXY_setValueAtX( &smr, XY, 0.2 * i, 0.7 + i + .1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkDullEdges( &smr, XY, 1e-10, 1e-10 );

    if( ptwXY_clear( &smr, XY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i = 0; i < nSame; i++ ) {
        if( ptwXY_setValueAtX( &smr, XY, -0.2 * i, 0.7 + i + .1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkDullEdges( &smr, XY, 1e-10, 1e-10 );

    ptwXY_free( XY );

    exit( errCount );
}
/*
************************************************************
*/
static int checkDullEdges( statusMessageReporting *smr, ptwXYPoints *data, double lowerEps, double upperEps ) {

    int errCount;

    errCount = checkDullEdges2( smr, data, lowerEps, upperEps, 0 );
    errCount += checkDullEdges2( smr, data,  lowerEps, upperEps, 1 );
    errCount += checkDullEdges2( smr, data, -lowerEps, upperEps, 0 );
    errCount += checkDullEdges2( smr, data, -lowerEps, upperEps, 1 );
    errCount += checkDullEdges2( smr, data,  lowerEps, -upperEps, 0 );
    errCount += checkDullEdges2( smr, data,  lowerEps, -upperEps, 1 );
    errCount += checkDullEdges2( smr, data, -lowerEps, -upperEps, 0 );
    errCount += checkDullEdges2( smr, data, -lowerEps, -upperEps, 1 );
    return( errCount );
}
/*
************************************************************
*/
static int checkDullEdges2( statusMessageReporting *smr, ptwXYPoints *data, double lowerEps, double upperEps, int positiveXOnly ) {

    int errCount = 0;
    ptwXYPoints *dullEdges;
    nfu_status status;

    printIfVerbose( data );
    if( verbose ) {
        printf( "# lowerEps = %.14e\n", lowerEps );
        printf( "# upperEps = %.14e\n", upperEps );
        printf( "# positiveXOnly = %d\n", positiveXOnly );
    }
    if( ( dullEdges = ptwXY_clone( smr, data ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ( status = ptwXY_dullEdges( smr, dullEdges, lowerEps, upperEps, positiveXOnly ) ) != nfu_Okay )
        nfut_printSMRErrorExit2p( smr, "Via." );
    printIfVerbose( dullEdges );

    ptwXY_free( dullEdges );

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
