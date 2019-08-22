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

#include <ptwXY.h>
#define nPairs 7

/*
********************************************************
*/
int main( int argv, char **agrv ) {

    double c;
    double maleData[2 * ( nPairs - 1 )] =   { 1871, 1212, 1883, 1215, 1889,  51, 1895, 11, 1905,  9, 1915,  9 };
    double femaleData[2 * nPairs] = { 1871, 1231, 1883, 1241, 1885, 621, 1889, 229, 1895, 31, 1905, 23, 1915, 21 };
    ptwXYPoints *males, *females, *total;
    nfu_status status;

    males = ptwXY_new( ptwXY_interpolationLinLin, 5, 1e-3, nPairs, 4, &status, 0 );
    ptwXY_setXYData( males, nPairs - 1, maleData );

    females = ptwXY_new( ptwXY_interpolationLinLin, 5, 1e-3, nPairs, 4, &status, 0 );
    ptwXY_setXYData( females, nPairs, femaleData );

    total = ptwXY_add_ptwXY( males, females, &status );

    printf( "\nMale population\n" );
    printf( "  Year | Count\n" );
    printf( " -------+-----\n" );
    ptwXY_simplePrint( males, " %5.0f | %5.0f\n" );

    printf( "\nFemale population\n" );
    printf( "  Year | Count\n" );
    printf( " ------+------\n" );
    ptwXY_simplePrint( females, " %5.0f | %5.0f\n" );

    printf( "\nTotal population\n" );
    printf( "  Year | Count\n" );
    printf( " ------+------\n" );
    ptwXY_simplePrint( total, " %5.0f | %5.0f\n" );

    ptwXY_getValueAtX( males, 1885, &c );
    printf( "\nInterpolated male population at 1885 is %.0f\n", c );

    exit( EXIT_SUCCESS );
}
