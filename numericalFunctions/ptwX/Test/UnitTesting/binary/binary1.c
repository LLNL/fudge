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

#include <nf_utilities.h>
#include <ptwX.h>

#define nXs 23

static int verbose = 0;

int check( ptwXPoints *ptwX1, ptwXPoints *ptwX2, char const *msg );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errors = 0;
    nfu_status status;
    ptwXPoints *ptwX1, *ptwX2;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfu_printErrorMsg( "Error %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    ptwX1 = ptwX_createLine( 10, nXs, -2, 100, &status );

    ptwX2 = ptwX_createLine( 10, nXs, 1, 0, &status );
    if( ( status = ptwX_slopeOffset( ptwX2, -2, 100 ) ) != nfu_Okay ) 
        nfu_printMsg( "Error %s: status = %d from ptwX_slopeOffset( ptwX2, -2, 100 )", __FILE__, status );
    errors += check( ptwX1, ptwX2, "ptwX_slopeOffset( ptwX2, -2, 100 )" );

    ptwX2 = ptwX_createLine( 10, nXs, 1, 0, &status );
    if( ( status = ptwX_mul_double( ptwX2, -2 ) ) != nfu_Okay ) 
        nfu_printMsg( "Error %s: status = %d from ptwX_mul_double( ptwX2, -2 )", __FILE__, status );
    if( ( status = ptwX_add_double( ptwX2, 100 ) ) != nfu_Okay ) 
        nfu_printMsg( "Error %s: status = %d from ptwX_add_double( ptwX2, 100 )", __FILE__, status );
    errors += check( ptwX1, ptwX2, "ptwX_mul_double then ptwX_add_double" );

    ptwX_free( ptwX1 );
    exit( errors );
}
/*
****************************************************************
*/
int check( ptwXPoints *ptwX1, ptwXPoints *ptwX2, char const *msg ) {

    nfu_status status;
    int close = ptwX_close( ptwX1, ptwX2, 3, 0., &status );

    if( close > 0 ) {
        nfu_printMsg( "Error %s: for %s values differ at index %d", __FILE__, msg, close - 1 ); }
    else if( close < 0 ) {
        nfu_printMsg( "Error %s: for %s status = %d: %s", __FILE__, msg, status, nfu_statusMessage( status ) );
    }
    ptwX_free( ptwX2 );
    return( close );
}
