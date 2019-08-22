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

#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "nf_utilities.h"

#ifdef WIN32
#include <float.h>
#define isfinite _finite
#endif

#define numberOfStaticDoubles ( 100 * 1000 )

static double *nfu_stringToListOfDoubles2( statusMessageReporting *smr, char const *str, char sep, int64_t *numberConverted, 
        char **endCharacter );
/*
========================================================================
*/
double *nfu_stringToListOfDoubles( statusMessageReporting *smr, char const *str, char sep, int64_t *numberConverted, 
        char **endCharacter ) {

    if( strchr( "0123456789.+-eE", sep ) != NULL ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badInput, "Invalid sep ='%c'.", sep );
        return( NULL );
    }

    *numberConverted = 0;
    *endCharacter = (char *) str;
    if( isspace( sep ) ) sep = ' ';             /* Make it the space character if any white space as it simplifies logic below. */
    return( nfu_stringToListOfDoubles2( smr, str, sep, numberConverted, endCharacter ) );
}
/*
========================================================================
*/
static double *nfu_stringToListOfDoubles2( statusMessageReporting *smr, char const *str, char sep, int64_t *numberConverted, 
        char **endCharacter ) {

    int64_t i1, i2, numberConverted_initial = *numberConverted;
    double *doublePtr = NULL, staticDoubles[numberOfStaticDoubles];
    nfu_status status = nfu_Okay;

    for( i1 = 0; i1 < numberOfStaticDoubles; i1++, (*numberConverted)++ ) {
        if(  *numberConverted == 0 ) {
            staticDoubles[i1] = strtod( str, endCharacter ); }
        else {                                  /* Check that there is one sep character and allow for arbitrary number of white spaces. */
            char const *str2 = str;

            while( isspace( *str2 ) ) ++str2;   /* Only need to check for white spaces before sep character as strtod will ignore after. */
            if( sep != ' ' ) {
                if( *str2 == sep ) {
                    ++str2; }
                else {
                    str2 = str;
                }
            }
            if( str < str2 ) staticDoubles[i1] = strtod( str2, endCharacter );
            if( str2 == (char const *) *endCharacter ) *endCharacter = (char *) str;
        }
        if( str == (char const *) *endCharacter ) {
            int64_t number = *numberConverted;
            if( *numberConverted == 0 ) number = 1;
            if( ( doublePtr = (double *) smr_malloc2( smr, (size_t) number * sizeof( double ), 0, "doublePtr" ) ) == NULL ) {
                smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                return( NULL );
            }
            break;
        }
        str = (char const *) *endCharacter;
    }

    if( ( status == nfu_Okay ) && ( doublePtr == NULL ) )
        doublePtr = nfu_stringToListOfDoubles2( smr, str, sep, numberConverted, endCharacter );
    if( doublePtr != NULL ) {
        double *doublePtr2 = &(doublePtr[numberConverted_initial]);
        char *end = *endCharacter;

        for( i2 = 0; i2 < i1; i2++, doublePtr2++ ) *doublePtr2 = staticDoubles[i2];
        while( isspace( *end ) ) ++end;
        if( *end == 0 ) *endCharacter = end;
    }
    return( doublePtr );
}
/*
============================================================
*/
char *nf_floatToShortestString( double value, int significantDigits, int favorEFormBy, int flags ) {

    int n1, ne, nf, digitsRightOfPeriod_f, exponent;
    char Str_e[512], Str_f[512], *Str_r = Str_e, Fmt[32], *e1, *e2;
    const char *sign = "";

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
