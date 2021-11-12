/*
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include "statusMessageReporting.h"
#include <stdlib.h>

static int verbose = 0;
static int check1ID = smr_unknownID;

static char const text1[] = "Hi";

void addMessage( statusMessageReporting *smr1, statusMessageReporting *smr2, enum smr_status status, int count, char const *message );
void printMessages( statusMessageReporting *smr1, statusMessageReporting *smr2, int clear );
void releaseMessages( statusMessageReporting *smr1, statusMessageReporting *smr2 );
enum smr_status getStatus( void );
void timeToAbort( int i );
/*
============================================================
*/
int main( int argc, char **argv ) {

    int i, counter;
    enum smr_status status;
    statusMessageReporting static_smr, *smr1 = &static_smr, *smr2;
    char msg[128];

    if( argc > 1 ) verbose = 1;
    smr_setup( );
    check1ID = smr_registerLibrary( "check1" );

    smr_initialize( smr1, smr_status_Warning );
    if( ( smr2 = smr_new( smr1, smr_status_Ok ) ) == NULL ) {
        smr_print( smr1, 1 );
        exit( EXIT_FAILURE );
    }

    addMessage( smr1, smr2, smr_status_Info, 1, text1 );
    printMessages( smr1, smr2, 0 );
    releaseMessages( smr1, smr2 );

    for( i = 0, counter = 0; i < 400000; i++ ) {
        counter++;
        status = getStatus( );
        if( verbose ) printf( "\ni = %d, status = %d, <%s>\n", i, status, smr_statusToString( status ) );
        sprintf( msg, "This is message %d with status = %d, string = <%s>\n", i, status, smr_statusToString( status ) );
        addMessage( smr1, smr2, status, 1, msg );
        if( drand48( ) < .1 ) {
            if( smr_numberOfReports( smr2 ) != counter ) fprintf( stderr, "smr_numberOfReports( smr2 ) = %d != counter = %d at index = %d\n", 
                smr_numberOfReports( smr2 ), counter, i );
            printMessages( smr1, smr2, 1 );
            counter = 0;
        }
    }
    printMessages( smr1, smr2, 1 );
    releaseMessages( smr1, smr2 );
    smr_free( &smr2 );
    smr_cleanup( );

    exit( EXIT_SUCCESS );
}
/*
============================================================
*/
void addMessage( statusMessageReporting *smr1, statusMessageReporting *smr2, enum smr_status status, int count, char const *message ) {

    int i;

    switch( status ) {
    case smr_status_Ok :
        smr_release( smr1 );
        smr_release( smr2 );
        break;
    case smr_status_Info :
        if( ( i = smr_setReportInfo( smr1, NULL, __FILE__, __LINE__, __func__, check1ID, 30, "count = %3d, message = %s", count, message ) ) ) timeToAbort( i );
        if( ( i = smr_setReportInfo( smr2, NULL, __FILE__, __LINE__, __func__, check1ID, 30, "count = %3d, message = %s", count, message ) ) ) timeToAbort( i );
        break;
    case smr_status_Warning :
        if( ( i = smr_setReportWarning( smr1, NULL, __FILE__, __LINE__, __func__, check1ID, 110, "count = %3d, message = %s", count, message ) ) ) timeToAbort( i );
        if( ( i = smr_setReportWarning( smr2, NULL, __FILE__, __LINE__, __func__, check1ID, 110, "count = %3d, message = %s", count, message ) ) ) timeToAbort( i );
        break;
    case smr_status_Error :
        if( ( i = smr_setReportError( smr1, NULL, __FILE__, __LINE__, __func__, check1ID, 510, "count = %3d, message = %s", count, message ) ) ) timeToAbort( i );
        if( ( i = smr_setReportError( smr2, NULL, __FILE__, __LINE__, __func__, check1ID, 510, "count = %3d, message = %s", count, message ) ) ) timeToAbort( i );
        break;
    default :
        fprintf( stderr, "Invalid status = %d\n", status );
        timeToAbort( 2 );
    }
}
/*
============================================================
*/
void printMessages( statusMessageReporting *smr1, statusMessageReporting *smr2, int clear ) {

    if( verbose ) {
        printf( "Printing smr1 (clear = %d) that has %d messages\n", clear, smr_numberOfReports( smr1 ) );
        smr_print( smr1, clear );
        printf( "Printing smr2 (clear = %d) that has %d messages\n", clear, smr_numberOfReports( smr2 ) );
        smr_print( smr2, clear ); }
    else {
        if( clear ) releaseMessages( smr1, smr2 );
    }
}
/*
============================================================
*/
void releaseMessages( statusMessageReporting *smr1, statusMessageReporting *smr2 ) {

    if( verbose ) printf( "Releasing smr1\n" );
    smr_release( smr1 );
    if( verbose ) printf( "Releasing smr2\n" );
    smr_release( smr2 );
}
/*
============================================================
*/
enum smr_status getStatus( void ) {

    int i = (int) 3 * drand48( );

    switch( i ) {
    case 0 : return( smr_status_Info );
    case 1 : return( smr_status_Warning );
    case 2 : return( smr_status_Error );
    }
    fprintf( stderr, "getStatus i value = %d is invalid\n", i );
    timeToAbort( 0 );
    return( smr_status_Ok );        /* Will never get here but need to make the compilers happy. */
}
/*
============================================================
*/
void timeToAbort( int i ) {

    fprintf( stderr, "i = %d\n", i );
    exit( EXIT_FAILURE );
}
