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

#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>
#include <nfut_utilities.h>

static int verbose = 0;
static int Bins[] = { 2, 3, 4, 5, 7, 10, 20, 100 };
static int nBins = sizeof( Bins ) / sizeof( Bins[0] );

static int epb( statusMessageReporting *smr, ptwXY_interpolation interpolation, int nXYs, double *XYs );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0;
    double XYs1[] = { 0, 2.643700e-2, 0.4, 9.808100e-1, 0.8, 1.435800e+0, 1.2, 5.697500e-2, 1.6, 0 };
    double XYs2[] = { 0, 9.164300e-4, 0.4, 3.272800e-2, 0.8, 9.592600e-1, 1.2, 1.250000e+0, 1.6, 2.571200e-1, 2, 0 };
    statusMessageReporting smr;
    int nXYs1 = sizeof( XYs1 ) / sizeof( XYs1[0] ) / 2;
    int nXYs2 = sizeof( XYs2 ) / sizeof( XYs2[0] ) / 2;

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

    errCount += epb( &smr, ptwXY_interpolationLinLin, nXYs1, XYs1 );
    errCount += epb( &smr, ptwXY_interpolationFlat, nXYs1, XYs1 );
    errCount += epb( &smr, ptwXY_interpolationLinLin, nXYs2, XYs2 );
    errCount += epb( &smr, ptwXY_interpolationFlat, nXYs2, XYs2 );

    exit( errCount );
}
/*
************************************************************
*/
static int epb( statusMessageReporting *smr, ptwXY_interpolation interpolation, int nXYs, double *XYs ) {

    int i1, i2, errCount = 0;
    ptwXYPoints *ptwXY;
    ptwXPoints *ptwX;
    double integral;

    if( ( ptwXY = ptwXY_create( smr, interpolation, NULL, 4, 1.e-3, 10, 10, nXYs, XYs, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( smr, "Via." );

    if( ptwXY_integrateDomain( smr, ptwXY, &integral ) != nfu_Okay )
        nfut_printSMRErrorExit2p( smr, "Via." );

    nfu_printXYDataOnVerbosity( verbose, ptwXY );
    for( i1 = 0; i1 < nBins; ++i1 ) {
        if( verbose ) printf( "# numberOfBins = %d\n", Bins[i1] );
        if( ( ptwX = ptwXY_equalProbableBins( smr, ptwXY, Bins[i1] ) ) == NULL ) {
            errCount++; }
        else {
            nfu_printXDataOnVerbosity( verbose, ptwX );
            for( i2 = 0; i2 < Bins[i1]; ++i2 ) {
                double x1 = ptwX_getPointAtIndex_Unsafely( ptwX, i2 );
                double x2 = ptwX_getPointAtIndex_Unsafely( ptwX, i2 + 1 );
                double subIntegral;
                nfu_status status = ptwXY_integrate( smr, ptwXY, x1, x2, &subIntegral );

                if( status != nfu_Okay ) {
                    nfu_printSMRError2p( smr, "Via" );
                    break;
                }
                if( nfu_cmpDoubles( subIntegral, integral / Bins[i1], 1e-13 ) ) {
                    fprintf( stderr, "Boundary error: nXYs = %d  i1 = %d  number of bins = %d\n", nXYs, i1, Bins[i1] );
                    ++errCount;
                    break;
                }
            }
            ptwX_free( ptwX );
        }
    }

    ptwXY_free( ptwXY );
    return( errCount );
}
