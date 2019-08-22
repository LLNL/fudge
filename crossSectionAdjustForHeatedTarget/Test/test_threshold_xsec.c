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
*    This routine test the ???.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../Src/crossSectionAdjustForHeatedTarget.h"

static long nPoints = 200;
static double Temperature = 1e-3, EMin = 1., EMax = 20., massRatio = 8.;

static double Theory( double E );
static void PrintMsg( char *s, int Usage );
/*
***************************************
*/
int main( int argc, char **argv ) {

    int err;
    long i;
    double f, *E_cs, E, *p, t;
    char *e;
    crossSectionAdjustForHeatedTarget_info info = { 0, 0 };

    if( argc > 6 ) PrintMsg( "too many parameters", 1 );
    if( argc > 1 ) {
        if( argv[1][0] == '?' ) PrintMsg( NULL, 1 );
        f = strtod( argv[1], &e );
        if( ( e == argv[1] ) || ( *e != 0 ) ) PrintMsg( "cannot convert temperature parameter to double", 1 );
        Temperature = f;
        if( argc > 2 ) {
            f = strtod( argv[2], &e );
            if( ( e == argv[2] ) || ( *e != 0 ) ) PrintMsg( "cannot convert EMin parameter to double", 1 );
            EMin = f;
            if( argc > 3 ) {
                f = strtod( argv[3], &e );
                if( ( e == argv[3] ) || ( *e != 0 ) ) PrintMsg( "cannot convert EMax parameter to double", 1 );
                EMax = f;
                if( argc > 4 ) {
                    f = strtod( argv[4], &e );
                    if( ( e == argv[4] ) || ( *e != 0 ) ) PrintMsg( "cannot convert massRatio parameter to double", 1 );
                    massRatio = f;
                    if( argc > 5 ) {
                        i = strtol( argv[5], &e, 10 );
                        if( ( e == argv[5] ) || ( *e != 0 ) ) PrintMsg( "cannot convert nPoints parameter to long", 1 );
                        nPoints = i;
                    }
                }
            }
        }
    }
    if( ( EMin >= EMax ) || ( EMin <= 0. ) ) PrintMsg( "EMax must be greater than EMin and EMin must be greater than 0.", 0 );
    if( nPoints < 2 ) PrintMsg( "nPoints must be greater than 2", 0 );
    printf( "# T = %e  EMin = %e  EMax = %e  massRatio = %e  nPoints = %d\n", Temperature, EMin, EMax, massRatio, (int) nPoints );
    f = pow( EMax / EMin, 1. / ( nPoints - 1 ) );
    E_cs = malloc( 2 * nPoints * sizeof( double ) );
    if( E_cs == NULL ) PrintMsg( "cannot allocate memory for E_cs", 0 );

    E = EMin;
    p = E_cs;
    for( i = 0; i < nPoints; i++, E *= f ) {
        *(p++) = E;
        *(p++) = 1.;
    }
    p[-2] = EMax;

    err = crossSectionAdjustForHeatedTarget( crossSectionAdjustForHeatedTarget_limit_threshold, crossSectionAdjustForHeatedTarget_limit_constant, 
        &info, 1e-10, massRatio, Temperature, 5e-3, nPoints, E_cs, &p );
    if( err < 0 ) {
        fprintf( stderr, "\ncrossSectionAdjustForHeatedTarget returned err = %d\n", err );
        PrintMsg( "bad crossSectionAdjustForHeatedTarget return value", 0 );
    }
    for( i = 0; i < err; i++, p += 2 ) {
        t = Theory( p[0] );
        printf( "%16.8e %16.8e %16.8e %12.4e\n", p[0], p[1], t, p[1] / t - 1. );
    }
    exit( EXIT_SUCCESS );
}
/*
***************************************
*/
static double Theory( double E ) {

    double x, x2 = massRatio * E / Temperature, sqrtpi = 1.7724538509055160273, cs;

    x = sqrt( x2 );
    cs =  erf( x ) * ( 1. + 0.5 / x2 ) + exp( -x2 ) / ( sqrtpi * x );
    return( cs );
}
/*
***************************************
*/
static void PrintMsg( char *s, int Usage ) {

    if( s != NULL ) {
        fprintf( stderr, "\nError - %s\n\n", s );
    }

    if( Usage ) {
        printf( "\nUSAGE:\n" );
        printf( "    test_threshold_xsec Temperature, EMin, EMax, massRatio, nPoints\n" );
        printf( "\n" );
        printf( "PARAMETERS\n" );
        printf( "    Temperature    Target's temperature (default = %e).\n", Temperature );
        printf( "    EMin           Lowest cross section energy point (default = %e).\n", EMin );
        printf( "    EMax           Highest cross section energy point (default = %e).\n", EMax );
        printf( "    massRatio      Target to incident particle mass ratio (default = %e).\n", massRatio );
        printf( "    nPoints        Number of pairs of (E, xsec) points (default = %d).\n", (int) nPoints );
        printf( "\n" );
    }
    exit( EXIT_FAILURE );
}
