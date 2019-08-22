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

#include <stdlib.h>
#include <math.h>

#include "nf_utilities.h"

#define numberOfStaticDoubles ( 100 * 1000 )

static nfu_status nfu_stringToListOfDoubles2( char const *str, int64_t *numberConverted, double **doublePtr, char **endCharacter );
/*
========================================================================
*/
nfu_status nfu_stringToListOfDoubles( char const *str, int64_t *numberConverted, double **doublePtr, char **endCharacter ) {

    *numberConverted = 0;
    *doublePtr = NULL;
    return( nfu_stringToListOfDoubles2( str, numberConverted, doublePtr, endCharacter ) );
}
/*
========================================================================
*/
static nfu_status nfu_stringToListOfDoubles2( char const *str, int64_t *numberConverted, double **doublePtr, char **endCharacter ) {

    int64_t i1, i2, numberConverted_initial = *numberConverted;
    double staticDoubles[numberOfStaticDoubles];
    nfu_status status = nfu_Okay;

    for( i1 = 0; i1 < numberOfStaticDoubles; i1++, (*numberConverted)++ ) {
        staticDoubles[i1] = strtod( str, endCharacter );
        if( str == (char const *) *endCharacter ) {
            if( *numberConverted > 0 ) {
                if( ( *doublePtr = (double *) nfu_malloc( (size_t) *numberConverted * sizeof( double ) ) ) == NULL ) status = nfu_mallocError;
            }
            break;
        }
        str = (char const *) *endCharacter;
    }

    if( ( status == nfu_Okay ) && ( *doublePtr == NULL ) ) status = nfu_stringToListOfDoubles2( str, numberConverted, doublePtr, endCharacter );
    if( *doublePtr != NULL ) {
        double *doublePtr2 = &((*doublePtr)[numberConverted_initial]);

        for( i2 = 0; i2 < i1; i2++, doublePtr2++ ) *doublePtr2 = staticDoubles[i2];
    }
    return( status );
}
/*
============================================================
*/
char *nf_floatToShortestString( double value, int significantDigits, int favorEFormBy, int flags ) {

    int n1, ne, nf, digitsRightOfPeriod_f, exponent;
    char Str_e[512], Str_f[512], *Str_r = Str_e, Fmt[32], *e1, *e2, *sign = "";

    if( flags & nf_floatToShortestString_includeSign ) sign = "+";

    if( !isfinite( value ) ) {
        sprintf( Fmt, "%%%sf", sign );
        sprintf( Str_e, Fmt, value );
        return( strdup( Str_e ) );
    }

    significantDigits--;
    if( significantDigits < 0 ) significantDigits = 0;
    if( significantDigits > 24 ) significantDigits = 24;

    sprintf( Fmt, "%%%s.%de", sign, significantDigits );
    sprintf( Str_e, Fmt, value );

    e1 = strchr( Str_e, 'e' );
    if( significantDigits == 0 ) {
        if( *(e1 - 1) != '.' ) {
            char *e3;

            e2 = strchr( e1, 0 );
            e3 = e2 + 1;
            for( ; e2 != e1; e2--, e3-- ) *e3 = *e2;
            *(e1++) = '.';
        }
    }
    *e1 = 0;
    n1 = (int) strlen( Str_e ) - 1;
    if( flags & nf_floatToShortestString_trimZeros ) while( Str_e[n1] == '0' ) n1--;
    ne = flags & nf_floatToShortestString_keepPeriod;
    if( !( flags & nf_floatToShortestString_keepPeriod ) ) if( Str_e[n1] == '.' ) n1--;
    n1++;
    Str_e[n1] = 0;

    e1++;
    exponent = (int) strtol( e1, &e2, 10 );
    if( exponent != 0 ) {               /* If 0, the exponent was "e+00". */
        for( e1 = Str_e; *e1 != 0; e1++ ) ;
        sprintf( e1, "e%d", exponent );

        digitsRightOfPeriod_f = significantDigits - exponent;
        if( ( digitsRightOfPeriod_f > 25 ) || ( exponent > 50 ) ) return( strdup( Str_r ) );
        if( digitsRightOfPeriod_f < 0 ) digitsRightOfPeriod_f = 0;

        sprintf( Fmt, "%%%s.%df", sign, digitsRightOfPeriod_f );
        sprintf( Str_f, Fmt, value );

        ne = (int) strlen( Str_e );
        nf = (int) strlen( Str_f );
        if( strchr( Str_f, '.' ) != NULL ) {        /* '.' in string. */
            if( flags & nf_floatToShortestString_trimZeros ) while( Str_f[nf-1] == '0' ) nf--;
            if( Str_f[nf-1] == '.' ) {
                if( !( flags & nf_floatToShortestString_keepPeriod ) ) nf--;
            } }
        else {      /* Maybe we want a '.' else it looks like an integer, "12345." vs "12345". */
            if( flags & nf_floatToShortestString_keepPeriod ) {
                Str_f[nf] = '.';
                nf++;
            }
        }
        Str_f[nf] = 0;

        if( ( nf + favorEFormBy ) < ne ) Str_r = Str_f;
    }
    return( strdup( Str_r ) );
}
