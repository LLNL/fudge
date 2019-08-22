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
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
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
