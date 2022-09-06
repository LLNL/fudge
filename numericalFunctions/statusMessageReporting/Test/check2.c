/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include "statusMessageReporting.h"
#include <stdlib.h>

static int verbose = 0;

/*
============================================================
*/
int main( int argc, char **argv ) {

    int i, ID, check2ID, check8ID;
    statusMessageReporting static_smr, *smr1 = &static_smr, *smr2;
    char lName[32];

    if( argc > 1 ) verbose = 1;
    smr_setup( );

    for( i = 0; i < smr_numberOfRegisteredLibraries( ); ++i ) printf( "%3d <%s>\n", i, smr_getRegisteredLibrarysName( i ) );
    printf( "\n" );

    for( i = 0; i < 10; ++i ) {
        sprintf( lName, "lib%3.3d", i );
        ID = smr_registerLibrary( lName );
        if( i == 2 ) check2ID = ID;
        if( i == 8 ) check8ID = ID;
    }
    for( i = 0; i < smr_numberOfRegisteredLibraries( ); ++i ) printf( "%3d <%s>\n", i, smr_getRegisteredLibrarysName( i ) );

    if( smr_initialize( smr1, smr_status_Warning ) != 0 ) {
        fprintf( stderr, "smr_initialize failed for smr1\n" );
        exit( EXIT_FAILURE );
    }
    if( ( smr2 = smr_new( smr1, smr_status_Ok ) ) == NULL ) {
        smr_print( smr1, 1 );
        exit( EXIT_FAILURE );
    }

    smr_setReportInfo( smr1, NULL, __FILE__, __LINE__, __func__, check2ID, 5, "smr1" );
    smr_setReportInfo( smr2, NULL, __FILE__, __LINE__, __func__, check8ID, 40, "smr2" );

    printf( "\n" );
    i = smr_getLibraryID( smr_firstReport( smr1 ) );
    printf( "smr1 ID = %3d libararyName = %s\n", i, smr_getRegisteredLibrarysName( i ) );
    i = smr_getLibraryID( smr_firstReport( smr2 ) );
    printf( "smr2 ID = %3d libararyName = %s\n", i, smr_getRegisteredLibrarysName( i ) );

    smr_free( &smr2 );
    smr_cleanup( );

    exit( EXIT_SUCCESS );
}
