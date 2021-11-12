/*
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdlib.h>
#include <stdarg.h>
#include "statusMessageReporting.h"

static int verbose = 0, ID;

void _smr_setReportInfo2( statusMessageReporting *smr, char const *fmt, ... );
void _smr_setReportWarning2( statusMessageReporting *smr, char const *fmt, ... );
void _smr_setReportError2( statusMessageReporting *smr, char const *fmt, ... );
/*
============================================================
*/
int main( int argc, char **argv ) {

    int i;
    statusMessageReporting static_smr, *smr1 = &static_smr;

    if( argc > 1 ) verbose = 1;
    smr_setup( );

    for( i = 0; i < smr_numberOfRegisteredLibraries( ); ++i ) printf( "%3d <%s>\n", i, smr_getRegisteredLibrarysName( i ) );
    printf( "\n" );

    ID = smr_registerLibrary( "check3");
    for( i = 0; i < smr_numberOfRegisteredLibraries( ); ++i ) printf( "%3d <%s>\n", i, smr_getRegisteredLibrarysName( i ) );
    printf( "\n" );

    if( smr_initialize( smr1, smr_status_Info ) != 0 ) {
        fprintf( stderr, "smr_initialize failed for smr1\n" );
        exit( EXIT_FAILURE );
    }


    smr_setReportInfo( smr1, NULL, __FILE__, __LINE__, __func__, ID, 7, "smr_setReportInfo, MACRO = %d", 0 );
    smr_print( smr1, 1 );
    smr_setReportInfo2( smr1, ID, 7, "smr_setReportInfo2, MACRO = %d", 1 );
    smr_print( smr1, 1 );
    _smr_setReportInfo2( smr1, "smr_vsetReportInfo2, MACRO = %d", 1 );
    printf( "\n" );

    smr_setReportWarning( smr1, NULL, __FILE__, __LINE__, __func__, ID, 7, "smr_setReportWarning, MACRO = %d", 0 );
    smr_print( smr1, 1 );
    smr_setReportWarning2( smr1, ID, 7, "smr_setReportWarning2, MACRO = %d", 1 );
    smr_print( smr1, 1 );
    _smr_setReportWarning2( smr1, "smr_vsetReportWarning2, MACRO = %d", 1 );
    printf( "\n" );

    smr_setReportError( smr1, NULL, __FILE__, __LINE__, __func__, ID, 7, "smr_setReportError, MACRO = %d", 0 );
    smr_print( smr1, 1 );
    smr_setReportError2( smr1, ID, 7, "smr_setReportError2, MACRO = %d", 1 );
    smr_print( smr1, 1 );
    _smr_setReportError2( smr1, "smr_vsetReportError2, MACRO = %d", 1 );
    printf( "\n" );

    smr_cleanup( );
    exit( EXIT_SUCCESS );
}
/*
============================================================
*/
void _smr_setReportInfo2( statusMessageReporting *smr, char const *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    smr_vsetReportInfo2( smr, ID, 7, fmt, &args );
    smr_print( smr, 1 );
    va_end( args );
}
/*
============================================================
*/
void _smr_setReportWarning2( statusMessageReporting *smr, char const *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    smr_vsetReportWarning2( smr, ID, 7, fmt, &args );
    smr_print( smr, 1 );
    va_end( args );
}
/*
============================================================
*/
void _smr_setReportError2( statusMessageReporting *smr, char const *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    smr_vsetReportError2( smr, ID, 7, fmt, &args );
    smr_print( smr, 1 );
    va_end( args );
}
