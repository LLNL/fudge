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
#include <stdint.h>
#include <math.h>
#include <ctype.h>

#include "nf_utilities.h"

#define numberOfStaticInt32s ( 100 * 1000 )

#ifndef INT32_MIN
#define INT32_MIN -2147483648
#define INT32_MAX 2147483647
#endif

static int32_t *nfu_stringToListOfInt32s_2( statusMessageReporting *smr, char const *str, char sep, int64_t *numberConverted, 
        char **endCharacter );
static int nfu_stringToInt32( statusMessageReporting *smr, char const *str, char **endCharacter, int32_t *value );
/*
========================================================================
*/
int32_t *nfu_stringToListOfInt32s( statusMessageReporting *smr, char const *str, char sep, int64_t *numberConverted, 
        char **endCharacter ) {

    if( strchr( "0123456789.+-eE", sep ) != NULL ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badInput, "Invalid sep ='%c'.", sep );
        return( NULL );
    }

    *numberConverted = 0;
    *endCharacter = (char *) str;
    if( isspace( sep ) ) sep = ' ';             /* Make it the space character if any white space as it simplifies logic below. */
    return( nfu_stringToListOfInt32s_2( smr, str, sep, numberConverted, endCharacter ) );
}
/*
========================================================================
*/
static int32_t *nfu_stringToListOfInt32s_2( statusMessageReporting *smr, char const *str, char sep, int64_t *numberConverted, 
        char **endCharacter ) {

    int64_t i1, i2, numberConverted_initial = *numberConverted;
    int32_t *Int32Ptr = NULL, staticInt32s[numberOfStaticInt32s];

    for( i1 = 0; i1 < numberOfStaticInt32s; i1++, (*numberConverted)++ ) {
        if(  *numberConverted == 0 ) {
            if( nfu_stringToInt32( smr, str, endCharacter, &staticInt32s[i1] ) != 0 ) {
                *endCharacter = (char *) str;
                smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                return( NULL );
            } }
        else {                                  /* Check that there is one sep character and allow for arbitrary number of white spaces. */
            char const *str2 = str;

            while( isspace( *str2 ) ) ++str2;   /* Only need to check for white spaces before sep character as strtol will ignore after. */
            if( sep != ' ' ) {
                if( *str2 == sep ) {
                    ++str2; }
                else {
                    str2 = str;
                }
            }
            if( str < str2 ) {
                if( nfu_stringToInt32( smr, str2, endCharacter, &staticInt32s[i1] ) != 0 ) {
                    *endCharacter = (char *) str;
                    smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                    return( NULL );
                }
            }
            if( str2 == (char const *) *endCharacter ) *endCharacter = (char *) str;
        }
        if( str == (char const *) *endCharacter ) {
            int64_t number = *numberConverted;

            if( *numberConverted == 0 ) number = 1;
            if( ( Int32Ptr = (int32_t *) smr_malloc2( smr, (size_t) number * sizeof( int32_t ), 0, "Int32Ptr" ) ) == NULL ) {
                smr_setReportError2p( smr, nfu_SMR_libraryID, nfu_Error, "Via." );
                return( NULL );
            }
            break;
        }
        str = (char const *) *endCharacter;
    }

    if( Int32Ptr == NULL ) Int32Ptr = nfu_stringToListOfInt32s_2( smr, str, sep, numberConverted, endCharacter );
    if( Int32Ptr != NULL ) {
        int32_t *Int32Ptr2 = &(Int32Ptr[numberConverted_initial]);
        char *end = *endCharacter;

        for( i2 = 0; i2 < i1; i2++, Int32Ptr2++ ) *Int32Ptr2 = staticInt32s[i2];
        while( isspace( *end ) ) ++end;
        if( *end == 0 ) *endCharacter = end;
    }
    return( Int32Ptr );
}
/*
========================================================================
*/
static int nfu_stringToInt32( statusMessageReporting *smr, char const *str, char **endCharacter, int32_t *value ) {

    long lValue = strtol( str, endCharacter, 10 );

    if( lValue < INT32_MIN ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badInput, "int32_t underflow: %l", lValue );
        return( -1 ); }
    else if( lValue > INT32_MAX ) {
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badInput, "int32_t overflow: %l", lValue );
        return( 1 );
    }
    *value = (int) lValue;
    return( 0 );
}
