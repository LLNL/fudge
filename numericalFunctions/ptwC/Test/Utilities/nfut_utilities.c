/*
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include <statusMessageReporting.h>

#include "nfut_utilities.h"

/*
********************************************************
*/
long nfut_charToLong( statusMessageReporting *smr, char const *msg, char const *stringValue ) {

    char *endCharacter;
    long value = strtol( stringValue, &endCharacter, 10 );

    if( *endCharacter != 0 ) nfut_printSMRErrorExit2( smr, "ERROR %s: <%s> not a valid long", msg, stringValue );
    return( value );
}
/*
********************************************************
*/
int nfut_cmpDoubles( double d1, double d2, double espilon ) {

    double diff = d2 - d1, dMax = fabs( d1 );

    if( diff == 0 ) return( 0 );
    if( fabs( d2 ) > dMax ) dMax = fabs( d2 );
    if( fabs( diff ) > espilon * dMax ) return( 1 );
    return( 0 );
}
/*
********************************************************
*/
void nfut_printSMRError( statusMessageReporting *smr, char const *file, int line, char const *function, char const *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    smr_vsetReportError( smr, NULL, file, line, function, nfu_SMR_libraryID, 0, fmt, &args ); /* FIXME 0 needs to be something else. */
    va_end( args );

    smr_write( smr, stderr, 1 );
}
/*
********************************************************
*/
void nfut_printSMRErrorExit( statusMessageReporting *smr, char const *file, int line, char const *function, char const *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    smr_vsetReportError( smr, NULL, file, line, function, nfu_SMR_libraryID, 0, fmt, &args ); /* FIXME 0 needs to be something else. */
    va_end( args );

    smr_write( smr, stderr, 1 );
    exit( EXIT_FAILURE );
}
