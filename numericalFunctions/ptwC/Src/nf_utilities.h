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

#ifndef nf_utilities_h_included
#define nf_utilities_h_included

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdarg.h>

#define NUMERICALFUNCTIONS_SVN_VERSION 110+

#define nf_floatToShortestString_trimZeros   ( 1 << 0 )
#define nf_floatToShortestString_keepPeriod  ( 1 << 1 )
#define nf_floatToShortestString_includeSign ( 1 << 2 )

#if defined __cplusplus
    extern "C" {
#endif

typedef enum nfu_status_e {         nfu_Okay,               nfu_mallocError,            nfu_insufficientMemory,     
    nfu_badIndex,                   nfu_XNotAscending,      nfu_badIndexForX,           nfu_XOutsideDomain,             
    nfu_invalidInterpolation,       nfu_badSelf,            nfu_divByZero,              nfu_unsupportedInterpolationConversion, 
    nfu_unsupportedInterpolation,   nfu_empty,              nfu_tooFewPoints,           nfu_domainsNotMutual,                   
    nfu_badInput,                   nfu_badNorm,            nfu_badIntegrationInput,    nfu_otherInterpolation,
    nfu_failedToConverge,           nfu_oddNumberOfValues,  nfu_badLogValue } nfu_status;

/*
* Functions in nf_utilities.c
*/
double nfu_getNAN( void );
int nfu_isNAN( double d );
double nfu_getInfinity( double sign );
const char *nfu_statusMessage( nfu_status status );
void nfu_setMemoryDebugMode( int mode );
void *nfu_malloc( size_t size );
void *nfu_calloc( size_t size, size_t n );
void *nfu_realloc( size_t size, void *old );
void *nfu_free( void *p );
void nfu_printMsg( char *fmt, ... );
void nfu_printErrorMsg( char *fmt, ... );

/*
* Functions in nf_stringToDoubles.c
*/
double *nfu_stringToListOfDoubles( char const *str, char sep, int64_t *numberConverted, char **endCharacter, nfu_status *status );
char *nf_floatToShortestString( double value, int significantDigits, int favorEFormBy, int flags );

#if defined __cplusplus
    }
#endif

#endif          /* End of nf_utilities_h_included. */
