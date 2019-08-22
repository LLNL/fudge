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
#include <math.h>

void printData( double E, double cs );
/*
******************************************************
*/
int main( ) {

    int i, n = 41;
    double x, w = 1e-3, E, cs, fE, fw = 3., y, offset = 0.;

    fE = pow( ( 1. + fw * w ) / ( 1. - fw * w ), 1. / ( n - 1 ) );
    printData( 1e-10, offset );
    for( x = 1e-8; x < 11; x *= 10. ) {
        E = x * ( 1. - fw * w );
        for( i = 0; i < n ; i++ ) {
            y = ( E / x - 1. ) / w;
            cs = offset + 100. * exp( - y * y );
            printData( E, cs );
            E *= fE;
        }
    }
    printData( 20., offset );
    exit( EXIT_SUCCESS );
}
/*
******************************************************
*/
void printData( double E, double cs ) {

    printf( "%18.10e %15.6e\n", E, cs );
}
