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
/*
    This routine tests the nf_Legendre_GaussianQuadrature integrator by comparing the integral of x^n (n = 0 to MAX_DEGREE)
between x1 and x2 to the soluation ( x2^(n+1) - x1^(n+1) ) / ( n + 1 ).
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <ptwXY.h>
#include <nf_Legendre.h>

#define MAX_DEGREE 80
#define TOL 2e-13

static int verbose = 0;

static int GaussianQuadrature( int degree, double x1, double x2 );
static int GaussianQuadrature2( int degree, double x1, double x2 );
static nfu_status GaussianQuadrature_callback( double x, double *f, void *argList );
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


    for( i = 0; i <= MAX_DEGREE; i++ ) errCount += GaussianQuadrature( i, 0, 1 );
    for( i = 0; i <= MAX_DEGREE; i++ ) errCount += GaussianQuadrature( i, 4.5, 6.2 );
    for( i = 0; i <= MAX_DEGREE; i++ ) errCount += GaussianQuadrature( i, -1.2, .7 );

    if( errCount ) fprintf( stderr, "%s FAILED\n", __FILE__ );
    return( errCount );
}
/*
************************************************************
*/
static int GaussianQuadrature( int degree, double x1, double x2 ) {

    return( GaussianQuadrature2( degree, x1, x2 ) + GaussianQuadrature2( degree, x2, x1 ) );
}
/*
************************************************************
*/
static int GaussianQuadrature2( int degree, double x1, double x2 ) {

    int errCount = 0;
    double integral, exact, r;
    nfu_status status;

    status = nf_Legendre_GaussianQuadrature( degree, x1, x2, GaussianQuadrature_callback, (void *) &degree, &integral );
    if( status == nfu_Okay ) {
        exact = ( pow( x2, degree + 1 ) - pow( x1, degree + 1 ) ) / ( degree + 1 );
        r = integral / exact - 1.0;
        if( fabs( r ) > TOL ) errCount++;
        if( verbose ) printf( "%3d %25.16e %25.16e %20.12e %20.12e %12.5e\n", degree, x1, x2, exact, integral, r ); }
    else {
        errCount++;
    }
    return( errCount );
}
/*
************************************************************
*/
static nfu_status GaussianQuadrature_callback( double x, double *f, void *argList ) {

    int degree = *((int *) argList);

    *f = pow( x, degree );
    return( nfu_Okay );
}
