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
#include <ptwXY_utilities.h>

static int verbose = 0;
static char *fmtXY = "%.12e %.12e\n";

static void printIfVerbose( const char * const label, ptwXYPoints *data, double accuracy, double offset, double sigma, double amplitude, 
    double xMin, double xMax, double dullEps );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0;
    double accuracy = 1e-5, area, offset, amplitude = 5., sigma = 10., xMin, xMax;
    ptwXYPoints *gaussian;
    nfu_status status;

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

    if( ( gaussian = ptwXY_createGaussianCenteredSigma1( accuracy, &status ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: ptwXY_createGaussianCenteredSigma1, status = %d, '%s'", __FILE__, status, nfu_statusMessage( status ) );
    printIfVerbose( "ptwXY_createGaussianCenteredSigma1", gaussian, accuracy, 0., 1., 1., 0., 0., 0. );

    area = ptwXY_integrateDomain( gaussian, &status );
    if( status != nfu_Okay ) nfu_printErrorMsg( "ERROR %s: ptwXY_integrateDomain, status = %d, '%s'", __FILE__, status, nfu_statusMessage( status ) );
    if( verbose ) printf( "# area = %.12e  area err = %e\n", area, area / 2.50662827463 - 1 );
    if( fabs( area / 2.50662827463 - 1 ) > 0.2 * accuracy )
        nfu_printErrorMsg( "ERROR %s: area err = %e, area = %.12e", __FILE__, area / 2.50662827463 - 1, area );
    ptwXY_free( gaussian );

    offset = 110.;
    xMin = offset - 66.;
    xMax = offset + 66.;
    if( ( gaussian = ptwXY_createGaussian( accuracy, offset, sigma, amplitude, xMin, xMax, 0., &status ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: ptwXY_createGaussian, status = %d, '%s'", __FILE__, status, nfu_statusMessage( status ) );
    printIfVerbose( "ptwXY_createGaussian", gaussian, accuracy, offset, sigma, amplitude, xMin, xMax, 0. );
    area = ptwXY_integrateDomain( gaussian, &status );
    if( status != nfu_Okay ) nfu_printErrorMsg( "ERROR %s: ptwXY_integrateDomain, status = %d, '%s'", __FILE__, status, nfu_statusMessage( status ) );
    if( verbose ) printf( "# area = %.12e  area err = %e\n", area, area / ( sigma * amplitude * 2.50662827463 ) - 1 );
    if( fabs( area / ( sigma * amplitude * 2.50662827463 ) - 1 ) > 0.2 * accuracy )
        nfu_printErrorMsg( "ERROR %s: area err = %e, area = %.12e", __FILE__, area / ( sigma * amplitude * 2.50662827463 ) - 1, area );
    ptwXY_free( gaussian );

    return( 0 );
}
/*
************************************************************
*/
static void printIfVerbose( const char * const label, ptwXYPoints *data, double accuracy, double offset, double sigma, double amplitude, 
    double xMin, double xMax, double dullEps ) {

    if( !verbose ) return;
    printf( "# %s\n", label );
    printf( "# accuracy = %g\n", accuracy );
    printf( "# offset = %g\n", offset );
    printf( "# sigma = %g\n", sigma );
    printf( "# amplitude = %g\n", amplitude );
    printf( "# xMin = %.17e\n", xMin );
    printf( "# xMax = %.17e\n", xMax );
    printf( "# dullEps = %g\n", dullEps );
    printf( "# length = %d\n", (int) data->length );
    ptwXY_simpleWrite( data, stdout, fmtXY );
    printf( "\n\n" );
}
