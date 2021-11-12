/*
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void printMathCommand( double a, double x );
static char const *value2MathValue( double value, char str[128] );
/*
************************************************************
*/
char const *value2MathValue( double value, char str[128] );
int main( int argc, char **argv ) {

    double a, x;

    for( a = 1e-12; a < 20; a *= 2 ) {
        printMathCommand( a, 0. );
        for( x = 1e-10; x < 500; x *= 1.2 ) printMathCommand( a, x );
    }

    exit( EXIT_SUCCESS );
}
/*
************************************************************
*/
static void printMathCommand( double a, double x ) {

    char str1[128], str2[128];

    value2MathValue( a, str1 );
    value2MathValue( x, str2 );
    printf( "Print[ \"{ %+25.17e, %+25.17e, \", CForm[ N[ Gamma[ %s, 0, %s ], 20 ] ], \" },\" ]\n", a, x, str1, str2 );
}
/*
************************************************************
*/
static char const *value2MathValue( double value, char str[128] ) {

    int i;

    value *= 1e27;
    sprintf( str, "%.0f / 1", value );
    for( i = 0; i < 27; i++ ) strcat( str, "0" );
    return( str );
}
