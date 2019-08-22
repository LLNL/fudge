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

static int checkMerge( statusMessageReporting *smr, ptwXYPoints *XYs, double epsilon, int n );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i1, i2, iarg, echo = 0, errCount = 0;
    ptwXYPoints *XYs;
    double x, epsilon = 1e-10;
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

    if( ( XYs = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 10, 10, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );

    for( i1 = 0, x = 1; i1 < 5; i1++, x += epsilon / 8 ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, ( x - 1 ) / epsilon ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    for( ; i1 < 50; i1++, x += 2 * epsilon ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, ( x - 1 ) / epsilon ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkMerge( &smr, XYs, epsilon, 45 );

    for( i1 = 0, x = 1; i1 < 45; i1++, x += 2 * epsilon ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, ( x - 1 ) / epsilon ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    for( ; i1 < 50; i1++, x += epsilon / 8 ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, ( x - 1 ) / epsilon ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkMerge( &smr, XYs, epsilon, 46 );

    for( i1 = 0, x = 1; i1 < 15; i1++, x += 2 * epsilon ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, ( x - 1 ) / epsilon ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    for( ; i1 < 30; i1++, x += epsilon / 8 ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, ( x - 1 ) / epsilon ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    for( ; i1 < 50; i1++, x += 2 * epsilon ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, ( x - 1 ) / epsilon ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkMerge( &smr, XYs, epsilon, 36 );

    for( i1 = 0, x = 1; i1 < 100; i1++, x += epsilon / 8 ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, i1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkMerge( &smr, XYs, epsilon, 13 );

    for( i1 = 0, x = 1; i1 < 2; i1++, x += epsilon ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, i1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkMerge( &smr, XYs, epsilon / 3, 2 );

    for( i1 = 0, x = 1; i1 < 2; i1++, x += epsilon ) {
        if( ptwXY_setValueAtX( &smr, XYs, x, i1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += checkMerge( &smr, XYs, 3 * epsilon, 2 );

    for( i2 = 0; i2 < 2; ++i2 ) {
        for( i1 = 0, x = 1; i1 < 3; i1++, x += epsilon ) {
            if( ptwXY_setValueAtX( &smr, XYs, x, i1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
        }
        errCount += checkMerge( &smr, XYs, ( 0.333 + 2.667 * i2 ) * epsilon, 3 - i2);
    }

    ptwXY_free( XYs );
    exit( errCount );
}
/*
************************************************************
*/
static int checkMerge( statusMessageReporting *smr, ptwXYPoints *XYs, double epsilon, int n ) {

    int errCount = 0;
    nfu_status status;

    if( verbose ) printf( "# epsilon = %e\n", epsilon );
    printIfVerbose( XYs );
    if( ( status = ptwXY_mergeClosePoints( smr, XYs, epsilon ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    printIfVerbose( XYs );
    if( ptwXY_length( smr, XYs ) != n ) {
        nfu_printMsg( "ERROR %s: ptwXY_length( smr, XYs ) = %d != n = %d", __FILE__, ptwXY_length( smr, XYs ), n );
        errCount++;
    }
    if( ptwXY_clear( smr, XYs ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
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
