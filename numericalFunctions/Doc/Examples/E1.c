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

    males = ptwXY_new( ptwXY_interpolationLinLin, NULL, 5, 1e-3, nPairs, 4, &status, 0 );
    ptwXY_setXYData( males, nPairs - 1, maleData );

    females = ptwXY_new( ptwXY_interpolationLinLin, NULL, 5, 1e-3, nPairs, 4, &status, 0 );
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
