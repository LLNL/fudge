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
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ptwXY.h>
#include <nf_utilities.h>

static char *fmtXY = "%18.11e %18.11e\n";

void doProduct( ptwXYPoints *ptwXY3, ptwXYPoints *ptwXY4, double biSection );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    double f3[4] = { 0., 0., 1., 1. }, f4[4] = { 0., 1., 1., 0. }, accuracy = 1e-3;
    ptwXYPoints *ptwXY3, *ptwXY4;
    nfu_status status;

    if( ( ptwXY3 = ptwXY_create( ptwXY_interpolationLinLin, NULL, 0, accuracy, 10, 10, 2, f3, &status, 0 ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: ptwXY3 creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( ptwXY4 = ptwXY_create( ptwXY_interpolationLinLin, NULL, 0, accuracy, 10, 10, 2, f4, &status, 0 ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: ptwXY4 creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    printf( "# f3\n" );
    ptwXY_simpleWrite( ptwXY3, stdout, fmtXY );
    printf( "\n\n# f4\n" );
    ptwXY_simpleWrite( ptwXY4, stdout, fmtXY );

    doProduct( ptwXY3, ptwXY4, 0 );
    doProduct( ptwXY3, ptwXY4, 1 );
    doProduct( ptwXY3, ptwXY4, 2 );
    doProduct( ptwXY3, ptwXY4, 3 );
    doProduct( ptwXY3, ptwXY4, 4 );
    doProduct( ptwXY3, ptwXY4, 10 );

    ptwXY_free( ptwXY3 );
    ptwXY_free( ptwXY4 );
    exit( EXIT_SUCCESS );

}
/*
************************************************************
*/
void doProduct( ptwXYPoints *ptwXY3, ptwXYPoints *ptwXY4, double biSection ) {

    ptwXYPoints *product;
    nfu_status status;

    ptwXY_setBiSectionMax( ptwXY3, biSection );
    if( ( product = ptwXY_mul2_ptwXY( ptwXY3, ptwXY4, &status ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: multiplication status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    printf( "\n\n# biSection = %.1f\n", biSection );
    printf( "# length = %d\n", (int) ptwXY_length( product ) );
    ptwXY_simpleWrite( product, stdout, fmtXY );
    ptwXY_free( product );
}
