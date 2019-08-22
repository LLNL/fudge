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
/*
    This routine 
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_Legendre.h>

#define MAX_ORDER 16

static int verbose = 0, biSectionMax = 12;
static double accuracy = 1e-5;

static int norm( int maxOrder, double scale );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "nf_Legendre: %s\n", __FILE__ );


    for( i = 0; i <= MAX_ORDER; i++ ) errCount += norm( i, 0.412 * i + 0.21 );

    if( errCount ) fprintf( stderr, "%s FAILED\n", __FILE__ );
    return( errCount );
}
/*
************************************************************
*/
static int norm( int maxOrder, double scale ) {

    int i, errCount = 0;
    double integral1, integral2, Cls[MAX_ORDER+1];
    nf_Legendre *nfL;
    ptwXYPoints *XYs;
    nfu_status status;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    Cls[0] = scale;
    for( i = 0; i < MAX_ORDER; i++ ) Cls[i+1] = 0.5 * Cls[i];

    if( ( nfL = nf_Legendre_new( &smr, 10, maxOrder, Cls ) ) == NULL ) {
        nfut_printSMRError2p( &smr, "Via." ); }
    else {
        if( ( XYs = nf_Legendre_to_ptwXY( &smr, nfL, accuracy, biSectionMax, 0 ) ) == NULL ) {
            nfut_printSMRError2p( &smr, "Via." ); }
        else {
            status = ptwXY_integrateDomain( &smr, XYs, &integral1 );
            ptwXY_free( XYs );
            if( status != nfu_Okay ) {
                nfut_printSMRError2p( &smr, "Via." ); }
            else {
                if( nf_Legendre_normalize( &smr, nfL ) != nfu_Okay ) {
                    nfut_printSMRError2p( &smr, "Via." ); }
                else {
                    if( ( XYs = nf_Legendre_to_ptwXY( &smr, nfL, accuracy, biSectionMax, 0 ) ) == NULL ) {
                        nfut_printSMRError2p( &smr, "Via." ); }
                    else {
                        if( ptwXY_integrateDomain( &smr, XYs, &integral2 ) != nfu_Okay ) {
                            nfut_printSMRError2p( &smr, "Via." ); }
                        else {
                            if( verbose ) printf( "integral1 = %.16e,  integral2 = %.16e, %e\n", integral1, integral2, 1.0 - integral2 );
                            if( fabs( 1.0 - integral2 ) > accuracy ) errCount++;
                        }
                        ptwXY_free( XYs );
                    }
                }
            }
        }
        nf_Legendre_free( nfL );
        Cls[0] = 0;
        if( ( nfL = nf_Legendre_new( &smr, 10, maxOrder, Cls ) ) == NULL ) {
            nfut_printSMRError2p( &smr, "Via." ); }
        else {
            if( nf_Legendre_normalize( NULL, nfL ) != nfu_divByZero ) errCount++;
            nf_Legendre_free( nfL );
        }
    }

    return( errCount );
}
