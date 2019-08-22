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

    int i, iarg, echo = 0, errCount;

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
