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

typedef enum nfu_status_e {         nfu_Okay,           nfu_mallocError,            nfu_insufficientMemory,     
    nfu_badIndex,                   nfu_XNotAscending,  nfu_badIndexForX,           nfu_XOutsideDomain,             
    nfu_invalidInterpolation,       nfu_badSelf,        nfu_divByZero,              nfu_unsupportedInterpolationConversion, 
    nfu_unsupportedInterpolation,   nfu_empty,          nfu_tooFewPoints,           nfu_domainsNotMutual,                   
    nfu_badInput,                   nfu_badNorm,        nfu_badIntegrationInput,    nfu_otherInterpolation,
    nfu_failedToConverge,           nfu_oddNumberOfValues } nfu_status;

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
nfu_status nfu_stringToListOfDoubles( char const *str, int64_t *numberConverted, double **doublePtr, char **endCharacter );
char *nf_floatToShortestString( double value, int significantDigits, int favorEFormBy, int flags );

#if defined __cplusplus
    }
#endif

#endif          /* End of nf_utilities_h_included. */
