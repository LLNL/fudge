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
#include <stdlib.h>
#include <string.h>

static char const *value2MathValue( double value, char str[128] );
/*
************************************************************
*/
char const *value2MathValue( double value, char str[128] );
int main( int argc, char **argv ) {

    double a;
    char str1[128];

    for( a = 1e-12; a < 20; a *= 1.1 ) {
        value2MathValue( a, str1 );
        printf( "Print[ \"{ %+25.17e, \", CForm[ N[ Gamma[ %s, 0 ], 20 ] ], \" },\" ]\n", a, str1 );
    }

    exit( EXIT_SUCCESS );
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
