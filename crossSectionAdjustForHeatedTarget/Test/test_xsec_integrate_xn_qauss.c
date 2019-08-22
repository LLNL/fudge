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
*    This routine test the results from crossSectionAdjustForHeatedTarget_integrate_xn_qauss.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <math.h>

#include "../Src/crossSectionAdjustForHeatedTarget.h"

char *toMathematicaDouble( double d );
/*
***************************************
*/
int main( int argc, char **argv ) {

    long long j, i, k, ni = 6, nj = 6, in, id;
    double a, b, d, dxnerf_e[5], dxnerf_t[5];
    int o = -1, t1;

    if( argc == 2 ) o = atoi( argv[1] );

    if(      o == 0 ) {
        printf( "f[x_] = Erf[ x ] / 2\n" ); }
    else if( o == 1 ) {
        printf( "f[x_] = - Exp[ -x^2 ] / ( 2 * Sqrt[ Pi ] )\n" ); }
    else if( o == 2 ) {
        printf( "f[x_] = Erf[ x ] / 4 - x * Exp[ -x^2 ] / ( 2 * Sqrt[ Pi ] )\n" ); }
    else if( o == 3 ) {
        printf( "f[x_] = - ( x^2 + 1 ) * Exp[ -x^2 ] / ( 2 * Sqrt[ Pi ] )\n" ); }
    else if( o == 4 ) {
        printf( "f[x_] = 3 * Erf[ x ] / 8 - x * ( 2 * x^2 + 3 ) * Exp[ -x^2 ] / ( 4 * Sqrt[ Pi ] )\n" ); }
    else {
        o = -1;
    }

    for( j = 1; j <= nj; j++ ) {
        in = 1;
        id = 10;
        for( i = 0; i < ni; i++, id *= 10 ) {
            b = j / 2.;
            a = b  - ( (double) in ) / id;
            crossSectionAdjustForHeatedTarget_integrate_xn_qauss( crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode_exact, a, b, dxnerf_e );
            if( o >= 0 ) printf( "Print[ N[ f[%d / 2] - f[%d / 2 - %d / %d], 20 ]  / ( %s ) - 1.]\n", 
                (int) j, (int) j, (int) in, (int) id, toMathematicaDouble( dxnerf_e[o] ) );
            crossSectionAdjustForHeatedTarget_integrate_xn_qauss( crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode_taylor, a, b, dxnerf_t );
            if( o >= 0 ) printf( "Print[ N[ f[%d / 2] - f[%d / 2 - %d / %d], 20 ] / ( %s ) - 1.]\n",
                (int) j, (int) j, (int) in, (int) id, toMathematicaDouble( dxnerf_t[o] ) );
            printf( "1e-%d  ", (int) (i + 1) );
            for( k = 0; k < 5; k++ ) {
                d = dxnerf_e[k] / dxnerf_t[k] - 1.;
                printf( "%14.6e%s", d, fabs( d ) > 1e-7 ? "***  " : "     " );
            }
            printf( "\n" );
            if( o >= 0 ) printf( "\n" );
        }
        printf( "\n" );
    }
    if( o < 0 ) {
        crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode Mode = crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode_taylor;
        in = 100000;
        for( k = 0; k < 2; k++ ) {
            id = 0;
            t1 = clock( );
            for( i = 0; i < 100; i++ ) {
                a = 0.03 * i;
                b = a + .001;
                for( j = 0; j < in; j++ ) crossSectionAdjustForHeatedTarget_integrate_xn_qauss( Mode, a, b, dxnerf_t );
                id += in;
            }
            t1 = clock( ) - t1;
            printf( "Timing for %d calls = %.3f sec for %s\n", (int) id, ( (double) t1 ) / ((double) CLOCKS_PER_SEC),
                ( k == 0 ? "Taylor series" : "exact" ) );
            Mode = crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode_exact;
        }
    }
    exit( EXIT_SUCCESS );
}
/*
***************************************
*/
char *toMathematicaDouble( double d ) {

    static char *c, Str[132], PowerTen[32];

    sprintf( Str, "%25.16e", d );
    c = strchr( Str, 'e' );
    if( c != NULL ) {
        *(c++) = ' ';
        strcpy( PowerTen, c );
        strcpy( c, "10^" );
        c += 3;
        strcpy( c, PowerTen );
    }
    return( Str );
}
