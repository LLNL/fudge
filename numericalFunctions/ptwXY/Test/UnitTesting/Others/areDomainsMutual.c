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

static int verbose = 0;
static char *fmtXY = "%17.8e%17.8e\n";

static int addPointAndCheck( statusMessageReporting *smr, double x, double y, ptwXYPoints *data1, ptwXYPoints *data2, int areMutual );
static int checkAreMutual( statusMessageReporting *smr, ptwXYPoints *data1, ptwXYPoints *data2, int areMutual );
static int checkAreMutual2( statusMessageReporting *smr, ptwXYPoints *data1, ptwXYPoints *data2, int areMutual );
static int checkAreMutual3( statusMessageReporting *smr, ptwXYPoints *data1, ptwXYPoints *data2, int areMutual );
static void printIfVerbose( ptwXYPoints *data );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errCount = 0;
    ptwXYPoints *XY1, *XY2;
    double flat1[] = { 1.0, 1.0, 9.0, 1.0 }, triangle1[] = { -1.0, 0.0, 0.0, 1.0, 0.5, 0.0 };
    double xMin = flat1[0], xMax = flat1[2], xMid = 0.5 * ( xMin + xMax );
    double wedge1[] = { 10., 0., 11., 1. }, wedge2[] = { xMid, 0., xMax, 1. }, wedge3[] = { xMin, 1., xMid, 0. };
    double triangle2[] = { xMid - 1., 0., xMid, 1., xMid + 1., 0 };
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

    if( ( XY1 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 20, 10, 2, flat1, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( XY2 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 20, 10, 3, triangle1, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );

    errCount += checkAreMutual( &smr, XY1, XY2, 0 );

    errCount += addPointAndCheck( &smr, flat1[0], flat1[1], XY1, XY2, 0 );
    errCount += addPointAndCheck( &smr, xMid, flat1[3], XY1, XY2, 0 );
    errCount += addPointAndCheck( &smr, flat1[2], flat1[3], XY1, XY2, 0 );
    errCount += addPointAndCheck( &smr, 2 * flat1[2], flat1[3], XY1, XY2, 0 );
    ptwXY_free( XY2 );

    if( ( XY2 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 20, 10, 3, triangle1, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += addPointAndCheck( &smr, flat1[0], 0., XY1, XY2, 0 );
    errCount += addPointAndCheck( &smr, flat1[2], 0., XY1, XY2, 0 );
    ptwXY_free( XY2 );

    if( ( XY2 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 20, 10, 2, wedge1, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkAreMutual( &smr, XY1, XY2, 0 );
    errCount += addPointAndCheck( &smr, flat1[2], 0., XY1, XY2, 0 );
    errCount += addPointAndCheck( &smr, xMid, 0., XY1, XY2, 0 );
    errCount += addPointAndCheck( &smr, flat1[0], 0., XY1, XY2, 0 );
    ptwXY_free( XY2 );

    errCount += checkAreMutual( &smr, XY1, XY1, 1 );

    if( ( XY2 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 20, 10, 2, wedge2, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkAreMutual( &smr, XY1, XY2, 1 );
    ptwXY_free( XY2 );

    if( ( XY2 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 20, 10, 2, wedge3, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkAreMutual( &smr, XY1, XY2, 1 );
    ptwXY_free( XY2 );

    if( ( XY2 = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 4, 1.e-3, 20, 10, 3, triangle2, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    errCount += checkAreMutual( &smr, XY1, XY2, 1 );
    ptwXY_free( XY2 );

    ptwXY_free( XY1 );

    exit( errCount );
}
/*
************************************************************
*/
static int addPointAndCheck( statusMessageReporting *smr, double x, double y, ptwXYPoints *data1, ptwXYPoints *data2, int areMutual ) {

    if( ptwXY_setValueAtX( smr, data2, x, y ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );

    return( checkAreMutual( smr, data1, data2, areMutual ) );
}
/*
************************************************************
*/
static int checkAreMutual( statusMessageReporting *smr, ptwXYPoints *data1, ptwXYPoints *data2, int areMutual ) {

    ptwXYPoints *n1, *n2;
    int errCount = 0;

    errCount += checkAreMutual2( smr, data1, data2, areMutual );

    if( ( n1 = ptwXY_clone( smr, data1 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_neg( smr, n1 ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    errCount += checkAreMutual2( smr, n1, data2, areMutual );

    if( ( n2 = ptwXY_clone( smr, data2 ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_neg( smr, n2 ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    errCount += checkAreMutual2( smr, data1, n2, areMutual );

    errCount += checkAreMutual2( smr, n1, n2, areMutual );

    ptwXY_free( n1 );
    ptwXY_free( n2 );

    return( errCount );
}
/*
************************************************************
*/
static int checkAreMutual2( statusMessageReporting *smr, ptwXYPoints *d1, ptwXYPoints *d2, int areMutual ) {

    int errCount = checkAreMutual3( smr, d1, d2, areMutual );

    return( errCount + checkAreMutual3( smr, d2, d1, areMutual ) );
}
/*
************************************************************
*/
static int checkAreMutual3( statusMessageReporting *smr, ptwXYPoints *d1, ptwXYPoints *d2, int areMutual ) {

    int errCount = 0;
    nfu_status status;

    if( verbose ) printf( "# areMutual = %d\n", areMutual );
    printIfVerbose( d1 );
    printIfVerbose( d2 );

    if( ( status = ptwXY_areDomainsMutual( smr, d1, d2 ) ) != nfu_Okay ) {
        if( areMutual || ( status != nfu_domainsNotMutual ) ) {
            errCount++;
            nfut_printSMRError2p( smr, "Via." );
        } }
    else {
        if( !areMutual ) {
            errCount++;
            nfu_printMsg( "ERROR %s: ptwXY_areDomainsMutual, is true and it should be false", __FILE__ );
        }
    }

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
