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
#include <ptwX.h>
#include <ptwXY.h>
#include <nf_utilities.h>

#define nXYs 8
#define nXs 100

static int verbose = 0;
static char *fmtX = " %19.12e\n";
static char *fmtXY = "%19.12e %19.12e\n";

static int check( statusMessageReporting *smr, ptwXPoints *Xs, ptwXYPoints *XYs, char const *label );
static void printIfVerboseXYs( ptwXYPoints *XYs );
static void printIfVerboseXs( ptwXPoints *Xs );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errCount = 0;
    int64_t i1;
    ptwXYPoints *XYs;
    ptwXYPoint *point;
    ptwXPoints Xs;
    double xValue, domainMin, domainMax, x1, x2;
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

    ptwX_initialize( &smr, &Xs, nXs );

    if( ( XYs = ptwXY_new( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 10, 10, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i1 = 0; i1 < nXYs; ++i1 ) {
        xValue = 0.2 * i1 - .5;
        if( ptwXY_setValueAtX( &smr, XYs, xValue, xValue * xValue ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    if( ptwXY_domainMin( &smr, XYs, &domainMin ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_domainMax( &smr, XYs, &domainMax ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );

    if( ptwX_setPointAtIndex( &smr, &Xs, 0, domainMin - 5.0 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwX_setPointAtIndex( &smr, &Xs, 1, domainMin - 1.0 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += check( &smr, &Xs, XYs, "0" );

    if( ptwX_setPointAtIndex( &smr, &Xs, 1, domainMin + 0.7  * ( domainMax - domainMin ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += check( &smr, &Xs, XYs, "1" );

    if( ptwX_setPointAtIndex( &smr, &Xs, 0, domainMin + 0.3  * ( domainMax - domainMin ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += check( &smr, &Xs, XYs, "2" );

    if( ptwX_setPointAtIndex( &smr, &Xs, 0, domainMin + 0.03 * ( domainMax - domainMin ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += check( &smr, &Xs, XYs, "3" );

    if( ptwX_setPointAtIndex( &smr, &Xs, 1, domainMin + 0.93 * ( domainMax - domainMin ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += check( &smr, &Xs, XYs, "4" );

    if( ptwX_setPointAtIndex( &smr, &Xs, 1, domainMax + 2.0 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += check( &smr, &Xs, XYs, "5" );

    if( ptwX_setPointAtIndex( &smr, &Xs, 0, domainMax + 1.0 )  != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += check( &smr, &Xs, XYs, "6" );

    if( ptwX_clear( &smr, &Xs ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i1 = 0; i1 < nXYs; ++i1 ) {
        point = ptwXY_getPointAtIndex_Unsafely( XYs, i1 );
        if( ptwX_setPointAtIndex( &smr, &Xs, i1, point->x ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    errCount += check( &smr, &Xs, XYs, "7" );

    if( ptwX_clear( &smr, &Xs ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i1 = 0; i1 < nXYs; ++i1 ) {
        point = ptwXY_getPointAtIndex_Unsafely( XYs, i1 );
        x2 = point->x;
        if( i1 > 0 ) {
            if( ptwX_setPointAtIndex( &smr, &Xs, i1 - 1, 0.5 * ( x1 + x2 ) ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
        }
        x1 = x2;
    }
    errCount += check( &smr, &Xs, XYs, "8" );

    if( ptwX_release( &smr, &Xs ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    ptwXY_free( XYs );

    exit( errCount );
}
/*
************************************************************
*/
static int check( statusMessageReporting *smr, ptwXPoints *Xs, ptwXYPoints *XYs, char const *label ) {

    int errCount = 0;
    int64_t offset;
    ptwXPoints *Ys = NULL;

    if( verbose ) printf( "# New test: %s\n", label );
    printIfVerboseXs( Xs );
    printIfVerboseXYs( XYs );
    if( ( Ys = ptwXY_ysMappedToXs( smr, XYs, Xs, &offset ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( verbose ) printf( "# offset = %d\n", (int) offset );
    printIfVerboseXs( Ys );
    ptwX_free( Ys );

    return( errCount );
}
/*
************************************************************
*/
static void printIfVerboseXYs( ptwXYPoints *XYs ) {

    if( !verbose ) return;
    printf( "# length = %d\n", (int) XYs->length );
    ptwXY_simpleWrite( XYs, stdout, fmtXY );
    printf( "\n\n" );
}
/*
************************************************************
*/
static void printIfVerboseXs( ptwXPoints *Xs ) {

    if( !verbose ) return;
    printf( "# length = %d\n", (int) Xs->length );
    ptwX_simpleWrite( NULL, Xs, stdout, fmtX );
    printf( "\n\n" );
}
