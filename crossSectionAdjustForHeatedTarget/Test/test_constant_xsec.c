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
*    This routine test the results for a constant cross-section to theory.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../Src/crossSectionAdjustForHeatedTarget.h"

static long nPoints = 200;
static double Temperature = 1e-3, EMin = 1e-10, EMax = 20., massRatio = 8.;

static double Theory( double E );
static void PrintMsg( char *s, int Usage );
/*
***************************************
*/
int main( int argc, char **argv ) {

    int err, checker = 0;
    long i;
    double f, *E_cs, E, *p, t, rErr, absrErr, rErrMax = 0.;
    char *e;
    crossSectionAdjustForHeatedTarget_info info;
    FILE *fOut;

    info.mode = 0;
    info.verbose = 0;

    if( argc > 6 ) PrintMsg( "too many parameters", 1 );
    if( argc > 1 ) {
        if( strcmp( argv[1], "-check" ) == 0 ) {
            checker = 1; }
        else {
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
    }
    if( ( EMin >= EMax ) || ( EMin <= 0. ) ) PrintMsg( "EMax must be greater than EMin and EMin must be greater than 0.", 0 );
    if( nPoints < 2 ) PrintMsg( "nPoints must be greater than 2", 0 );
    if( !checker ) printf( "# T = %e  EMin = %e  EMax = %e  massRatio = %e  nPoints = %ld\n", Temperature, EMin, EMax, massRatio, nPoints );
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

    if( !checker ) {
        fOut = fopen( "test_constant_xsec.data", "w" );
        for( i = 0, p = E_cs; i < nPoints; i++, p += 2 ) fprintf( fOut, "%13.5e %13.5e\n", *p, *(p+1) );
        fclose( fOut );
    }

    err = crossSectionAdjustForHeatedTarget( crossSectionAdjustForHeatedTarget_limit_constant, crossSectionAdjustForHeatedTarget_limit_constant,
        &info, 2e-10, massRatio, Temperature, 1e-4, nPoints, E_cs, &p );
    if( err < 0 ) {
        fprintf( stderr, "\ncrossSectionAdjustForHeatedTarget returned err = %d, bytes = %d\n", err, info.bytes );
        PrintMsg( "bad crossSectionAdjustForHeatedTarget return value", 0 );
    }
    if( !checker ) printf( "# WarningStats = %d evaluations = %d, pairs = %d\n", info.WarningStats, info.evaluations, err );
    for( i = 0; i < err; i++, p += 2 ) {
        t = Theory( p[0] );
        rErr = p[1] / t - 1.;
        absrErr = fabs( rErr );
        if( absrErr > rErrMax ) rErrMax = absrErr;
        if( !checker ) printf( "%16.8e %16.8e %16.8e %12.4e %s\n", p[0], p[1], t, rErr, ( absrErr > 1e-6 ? " ****" : "" ) );
    }
    if( rErrMax > 1e-8 ) exit( EXIT_FAILURE );
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
        printf( "    test_constant_xsec Temperature, EMin, EMax, massRatio, nPoints\n" );
        printf( "\n" );
        printf( "PARAMETERS\n" );
        printf( "    Temperature    Target's temperature (default = %e).\n", Temperature );
        printf( "    EMin           Lowest cross section energy point (default = %e).\n", EMin );
        printf( "    EMax           Highest cross section energy point (default = %e).\n", EMax );
        printf( "    massRatio      Target to incident particle mass ratio (default = %e).\n", massRatio );
        printf( "    nPoints        Number of pairs of (E, xsec) points (default = %ld).\n", nPoints );
        printf( "\n" );
    }
    exit( EXIT_FAILURE );
}
