/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "statusMessageReporting.h"

#define SMR_InitialMessageSize 1024
#define SMR_IncrementMessageSize 1024

static int smrIsSetup = 0;
static char smr_mallocFailed[] = "statusMessageReporting could not allocate memory for message";
static char statusStringOk[] = "Ok",  statusStringInfo[] = "Info", 
    statusStringWarning[] = "Warning", statusStringError[] = "Error", statusStringInvalid[] = "Invalid";

static int numberOfRegisteredLibraries = 0;
static char unknownLibrary[] = "unknownID";
static char tooManyLibrary[] = "tooManyIDs";
static char invalidLibrary[] = "invalidID";
static char errnoLibrary[] = "errnoID";
static char smrLibrary[] = "statusMessageReporting";
static char *registeredLibraries[smr_maximumNumberOfRegisteredLibraries];

static statusMessageReport *smr_reportNew( void );
static int smr_reportInitialize( statusMessageReport *report );
static void smr_reportRelease( statusMessageReport *report );
static int smr_setReport( statusMessageReporting *smr, void *userInterface, char const *file, int line, char const *function, int libraryID, int code, 
    enum smr_status status, char const *fmt, va_list *args );
static int smr_setAllocationFailure( statusMessageReport *report, char const *file, int line, char const *function, char const *fmt, va_list *args );
static statusMessageReport *smr_firstReport2( statusMessageReporting const *smr );
statusMessageReport *smr_nextReport2( statusMessageReport const *report );
static void smr_write2( statusMessageReport const *report, FILE *f );
/*
============================================================
*/
int smr_setup( void ) {

    int i;

    if( smrIsSetup ) return( 0 );
    smrIsSetup = 1;
    for( i = 0; i < smr_maximumNumberOfRegisteredLibraries; ++i ) registeredLibraries[i] = NULL;
    registeredLibraries[smr_unknownID] = unknownLibrary;
    ++numberOfRegisteredLibraries;
    registeredLibraries[smr_tooManyIDs] = tooManyLibrary;
    ++numberOfRegisteredLibraries;
    registeredLibraries[smr_invalidID] = invalidLibrary;
    ++numberOfRegisteredLibraries;
    registeredLibraries[smr_errnoID] = errnoLibrary;
    ++numberOfRegisteredLibraries;
    registeredLibraries[smr_smrID] = smrLibrary;
    ++numberOfRegisteredLibraries;
    return( 1 );
}
/*
============================================================
*/
int smr_cleanup( void ) {

    int i;

    if( smrIsSetup == 0 ) return( 0 );
    for( i = smr_smrID + 1; i < numberOfRegisteredLibraries; ++i ) smr_freeMemory( (void **) &(registeredLibraries[i]) );
    numberOfRegisteredLibraries = 0;
    smrIsSetup = 0;

    return( 0 );
}
/*
============================================================
*/
int smr_registerLibrary( char const *libraryName ) {

    int i;

    if( smrIsSetup == 0 ) smr_setup( );
    for( i = 0; i < numberOfRegisteredLibraries; ++i ) {             /* Check if name is already registered. */
        if( strcmp( libraryName, registeredLibraries[i] ) == 0 ) return( i );
    }
    if( numberOfRegisteredLibraries == smr_maximumNumberOfRegisteredLibraries ) return( smr_tooManyIDs );
    if( ( registeredLibraries[numberOfRegisteredLibraries] = strdup( libraryName ) ) == NULL ) return( -1 );
    ++numberOfRegisteredLibraries;
    return( numberOfRegisteredLibraries - 1 );
}
/*
============================================================
*/
int smr_numberOfRegisteredLibraries( void ) {

    return( numberOfRegisteredLibraries );
}
/*
============================================================
*/
char const *smr_getRegisteredLibrarysName( int ID ) {

    if( ( ID < 0 ) || ( ID >= smr_maximumNumberOfRegisteredLibraries ) ) return( NULL );
    return( registeredLibraries[ID] );
}
/*
============================================================
*/
statusMessageReporting *smr_new( statusMessageReporting *smr, enum smr_status verbosity ) {

    statusMessageReporting *new_SMR;

    if( ( new_SMR = (statusMessageReporting *) smr_malloc2( smr, sizeof( statusMessageReporting ), 0, "new_SMR" ) ) == NULL ) return( NULL );
    smr_initialize( new_SMR, verbosity );
    return( new_SMR );
}
/*
============================================================
*/
int smr_initialize( statusMessageReporting *smr, enum smr_status verbosity ) {

    if( smr == NULL ) return( 0 );
    smr->verbosity = verbosity;
    smr_reportInitialize( &(smr->report) );
    return( 0 );
}
/*
============================================================
*/
statusMessageReporting *smr_clone( statusMessageReporting const *smr ) {

    if( smr == NULL ) return( NULL );
    return( smr_new( NULL, smr->verbosity ) );
}
/*
============================================================
*/
void smr_release( statusMessageReporting *smr ) {

    statusMessageReport *current, *next, *first = smr_firstReport2( smr );

    if( smr == NULL ) return;
    for( current = first; current != NULL; current = next ) {
        next = smr_nextReport2( current );
        smr_reportRelease( current );
        if( current != first ) smr_freeMemory( (void **) &current );
    }
    smr_initialize( smr, smr->verbosity );
}
/*
============================================================
*/
void *smr_free( statusMessageReporting **smr ) {

    if( smr == NULL ) return( NULL );
    if( *smr != NULL ) {
        smr_release( *smr );
        smr_freeMemory( (void **) smr );
    }
    return( *smr );
}
/*
============================================================
*/
static statusMessageReport *smr_reportNew( void ) {

    statusMessageReport *report;

    if( ( report = (statusMessageReport *) smr_malloc2( NULL, sizeof( statusMessageReport ), 0, "report" ) ) == NULL ) return( NULL );
    smr_reportInitialize( report );
    return( report );
}
/*
============================================================
*/
static int smr_reportInitialize( statusMessageReport *report ) {

    report->next = NULL;
    report->status = smr_status_Ok;
    report->libraryID = smr_unknownID;
    report->code = smr_codeNULL;
    report->line = -1;
    report->fileName[0] = 0;
    report->function[0] = 0;
    report->message = NULL;
    return( 0 );
}
/*
============================================================
*/
static void smr_reportRelease( statusMessageReport *report ) {

    if( report->message != NULL ) {
        if( report->message != smr_mallocFailed ) smr_freeMemory( (void **) &(report->message) );
    }
    smr_reportInitialize( report );
}
/*
============================================================
*/
static int smr_setReport( statusMessageReporting *smr, void *userInterface, char const *file, int line, char const *function, int libraryID, int code, 
    enum smr_status status, char const *fmt, va_list *args ) {

    char *userMsg;
    statusMessageReport *report, *next;

    if( smr == NULL ) return( 0 );
    if( (int) status < (int) smr->verbosity ) return( 0 );
    if( status == smr_status_Ok ) return( 0 );
	if( smr->report.status == smr_status_Ok ) {
		report = &smr->report; }
	else {
    	if( ( report = smr_reportNew( ) ) == NULL ) return( smr_setAllocationFailure( NULL, file, line, function, fmt, args ) );
    	for( next = smr_firstReport2( smr ); next->next != NULL; next = next->next );
    	next->next = report;
	}
    report->status = status;
    if( ( libraryID < 0 ) || ( libraryID >= numberOfRegisteredLibraries ) ) libraryID = smr_invalidID;
    report->libraryID = libraryID;
    report->code = code;
    report->line = line;
    if( file != NULL ) strncpy( report->fileName, file, smr_maximumFileNameSize );
    report->fileName[smr_maximumFileNameSize] = 0;
    if( function != NULL ) strncpy( report->function, function, smr_maximumFileNameSize );
    report->function[smr_maximumFileNameSize] = 0;

    if( ( report->message = smr_vallocateFormatMessage( fmt, args ) ) == NULL ) return( smr_setAllocationFailure( report, file, line, function, fmt, args ) );
    if( userInterface != NULL ) {
        if( ( userMsg = (*(smr_userInterface *) userInterface)( (void *) userInterface ) ) != NULL ) {
            int userSize = (int) strlen( userMsg );
            if( ( report->message = (char *) smr_realloc2( NULL, report->message, strlen( report->message ) + userSize + 2, "report->message" ) ) == NULL ) {
                free( userMsg );
                return( smr_setAllocationFailure( report, file, line, function, fmt, args ) );
            }
            strcat( report->message, userMsg );
            free( userMsg );
        }
    }
    return( 0 );
}
/*
============================================================
*/
static int smr_setAllocationFailure( statusMessageReport *report, char const *file, int line, char const *function, char const *fmt, va_list *args ) {

    vfprintf( stderr, fmt, *args );
    va_end( *args );
    fprintf( stderr, "\n    At line %d of %s in function %s\n", line, file, function );
    if( report != NULL ) {
        report->status = smr_status_Error;
        report->message = (char *) smr_mallocFailed;
        return( 1 );
    }
    return( -1 );
}
/*
============================================================
*/
int smr_setReportInfo( statusMessageReporting *smr, void *userInterface, char const *file, int line, char const *function, int libraryID, int code, char const *fmt, ... ) {

    int status;
    va_list args;

    va_start( args, fmt );
    status = smr_setReport( smr, userInterface, file, line, function, libraryID, code, smr_status_Info, fmt, &args );
    va_end( args );
    return( status );
}
/*
============================================================
*/
int smr_vsetReportInfo( statusMessageReporting *smr, void *userInterface, char const *file, int line, char const *function, int libraryID, int code, char const *fmt, va_list *args ) {

    return( smr_setReport( smr, userInterface, file, line, function, libraryID, code, smr_status_Info, fmt, args ) );
}
/*
============================================================
*/
int smr_setReportWarning( statusMessageReporting *smr, void *userInterface, char const *file, int line, char const *function, int libraryID, int code, char const *fmt, ... ) {

    int status;
    va_list args;

    va_start( args, fmt );
    status = smr_setReport( smr, userInterface, file, line, function, libraryID, code, smr_status_Warning, fmt, &args );
    va_end( args );
    return( status );
}
/*
============================================================
*/
int smr_vsetReportWarning( statusMessageReporting *smr, void *userInterface, char const *file, int line, char const *function, int libraryID, int code, char const *fmt, va_list *args ) {

    return( smr_setReport( smr, userInterface, file, line, function, libraryID, code, smr_status_Warning, fmt, args ) );
}
/*
============================================================
*/
int smr_setReportError( statusMessageReporting *smr, void *userInterface, char const *file, int line, char const *function, int libraryID, int code, char const *fmt, ... ) {

    int status;
    va_list args;

    va_start( args, fmt );
    status = smr_setReport( smr, userInterface, file, line, function, libraryID, code, smr_status_Error, fmt, &args );
    va_end( args );
    return( status );
}
/*
============================================================
*/
int smr_vsetReportError( statusMessageReporting *smr, void *userInterface, char const *file, int line, char const *function, int libraryID, int code, char const *fmt, va_list *args ) {

    return( smr_setReport( smr, userInterface, file, line, function, libraryID, code, smr_status_Error, fmt, args ) );
}
/*
============================================================
*/
enum smr_status smr_highestStatus( statusMessageReporting const *smr ) {

    enum smr_status status = smr_status_Ok;
    statusMessageReport const *report;

    if( smr == NULL ) return( smr_status_Ok );
    for( report = smr_firstReport( smr ); report != NULL; report = smr_nextReport( report ) ) if( report->status > status ) status = report->status;
    return( status );
}
/*
============================================================
*/
int smr_isOk( statusMessageReporting const *smr ) { 

    return( smr_highestStatus( smr ) == smr_status_Ok );
}
/*
============================================================
*/
int smr_isInfo( statusMessageReporting const *smr ) { 

    return( smr_highestStatus( smr ) == smr_status_Info );
}
/*
============================================================
*/
int smr_isWarning( statusMessageReporting const *smr ) { 

    return( smr_highestStatus( smr ) == smr_status_Warning );
}
/*
============================================================
*/
int smr_isError( statusMessageReporting const *smr ) { 

    return( smr_highestStatus( smr ) == smr_status_Error );
}
/*
============================================================
*/
int smr_isWarningOrError( statusMessageReporting const *smr ) { 

    enum smr_status status = smr_highestStatus( smr );

    return( ( status == smr_status_Warning ) || ( status == smr_status_Error ) );
}
/*
============================================================
*/
int smr_isReportOk( statusMessageReport const *report ) { 

    if( report == NULL ) return( 0 );
    return( report->status == smr_status_Ok );
}
/*
============================================================
*/
int smr_isReportInfo( statusMessageReport const *report ) { 

    if( report == NULL ) return( 0 );
    return( report->status == smr_status_Info );
}
/*
============================================================
*/
int smr_isReportWarning( statusMessageReport const *report ) { 

    if( report == NULL ) return( 0 );
    return( report->status == smr_status_Warning );
}
/*
============================================================
*/
int smr_isReportError( statusMessageReport const *report ) { 

    if( report == NULL ) return( 0 );
    return( report->status == smr_status_Error );
}
/*
============================================================
*/
int smr_isReportWarningOrError( statusMessageReport const *report ) { 

    if( report == NULL ) return( 0 );
    return( ( report->status == smr_status_Warning ) || ( report->status == smr_status_Error ) );
}
/*
============================================================
*/
int smr_numberOfReports( statusMessageReporting const *smr ) {

    int n = 0;
    statusMessageReport const *report;

    if( smr == NULL ) return( 0 );
    if( smr->report.status == smr_status_Ok ) return( 0 );
    for( report = smr_firstReport( smr ); report != NULL; report = smr_nextReport( report ) ) ++n;
    return( n );
}
/*
============================================================
*/
statusMessageReport const *smr_firstReport( statusMessageReporting const *smr ) {

    return( (statusMessageReport const *) smr_firstReport2( smr ) );
}
/*
============================================================
*/
static statusMessageReport *smr_firstReport2( statusMessageReporting const *smr ) {

    if( smr == NULL ) return( NULL );
    if( smr->report.status == smr_status_Ok ) return( NULL );
    return( &(((statusMessageReporting *) smr)->report) );
}
/*
============================================================
*/
statusMessageReport const *smr_nextReport( statusMessageReport const *report ) {

    return( smr_nextReport2( report ) );
}
/*
============================================================
*/
statusMessageReport *smr_nextReport2( statusMessageReport const *report ) {

    if( report == NULL ) return( NULL );
    return( report->next );
}
/*
============================================================
*/
enum smr_status smr_getVerbosity( statusMessageReporting const *smr ) {

    if( smr == NULL ) return( smr_status_Ok );
    return( smr->verbosity );
}
/*
============================================================
*/
int smr_getLibraryID( statusMessageReport const *report ) {

    if( report == NULL ) return( 0 );
    return( report->libraryID );
}
/*
============================================================
*/
int smr_getCode( statusMessageReport const *report ) {

    if( report == NULL ) return( -1 );
    return( report->code );
}
/*
============================================================
*/
int smr_getLine( statusMessageReport const *report ) {

    if( report == NULL ) return( -1 );
    return( report->line );
}
/*
============================================================
*/
char const *smr_getFile( statusMessageReport const *report ) {

    if( report == NULL ) return( NULL );
    return( report->fileName );
}
/*
============================================================
*/
char const *smr_getFunction( statusMessageReport const *report ) {

    if( report == NULL ) return( NULL );
    return( report->function );
}
/*
============================================================
*/
char const *smr_getMessage( statusMessageReport const *report ) {

    if( report == NULL ) return( NULL );
    return( report->message );
}
/*
============================================================
*/
char *smr_copyMessage( statusMessageReport const *report ) {

    if( report == NULL ) return( NULL );
    if( report->status == smr_status_Ok ) return( NULL );
    return( smr_allocateFormatMessage( report->message ) );
}
/*
============================================================
*/
char *smr_copyFullMessage( statusMessageReport const *report ) {

    if( report == NULL ) return( NULL );
    if( report->status == smr_status_Ok ) return( NULL );
    return( smr_allocateFormatMessage( "%s\n    At line %d of %s in function %s", report->message, report->line, report->fileName, report->function ) );
}
/*
============================================================
*/
void smr_print( statusMessageReporting *smr, int clear ) {

    smr_write( smr, stdout, clear );
}
/*
============================================================
*/
void smr_write( statusMessageReporting *smr, FILE *f, int clear ) {

    if( smr == NULL ) return;
    fprintf( f, "======= Status message reports =======\n" );
    smr_write2( smr_firstReport( smr ), f );
    if( clear ) smr_release( smr );
}
/*
============================================================
*/
static void smr_write2( statusMessageReport const *report, FILE *f ) {

    if( report == NULL ) return;

    smr_write2( smr_nextReport( report ), f );
    smr_reportWrite( report, f );
}
/*
============================================================
*/
void smr_reportPrint( statusMessageReport const *report ) {

    smr_reportWrite( report, stdout );
}
/*
============================================================
*/
void smr_reportWrite( statusMessageReport const *report, FILE *f ) {

    if( report->message != NULL ) fprintf( f, "%s\n    At line %d of %s in function %s\n", report->message, report->line, report->fileName, report->function );
}
/*
============================================================
*/
char const *smr_statusToString( enum smr_status status ) {

    switch( status ) {
    case smr_status_Ok : return( statusStringOk );
    case smr_status_Info : return( statusStringInfo );
    case smr_status_Warning : return( statusStringWarning );
    case smr_status_Error : return( statusStringError );
    }
    return( statusStringInvalid );
}
/*
============================================================
*/
char *smr_allocateFormatMessage( char const *fmt, ... ) {

    char *s;
    va_list args;

    va_start( args, fmt );
    s = smr_vallocateFormatMessage( fmt, &args );
    va_end( args );
    return( s );
}
/*
============================================================
*/
char *smr_vallocateFormatMessage( char const *fmt, va_list *args ) {

    int n, size = SMR_InitialMessageSize;
    char buffer[SMR_InitialMessageSize], *message = buffer;
    va_list args_;

    while( 1 ) {
        va_copy( args_, *args );
        n = vsnprintf( message, size, fmt, args_ );
        va_end( args_ );
        if( ( n > -1 ) && ( n < size ) ) break;
        if( n > -1 ) {      /* glibc 2.1 */
            size = n + 3; }
        else {              /* glibc 2.0 */
            size *= 2;
        }
        if( message == buffer ) message = NULL;
        if( ( message = (char *) realloc( message, size ) ) == NULL ) return( NULL );
    }
    if( message == buffer ) {
        if( ( message = (char *) malloc( n + 1 ) ) == NULL ) return( NULL );
        strcpy( message, buffer ); }
    else {
        char *old = message;
        if( ( message = (char *) realloc( message, n + 1 ) ) == NULL ) {
            if( old != NULL ) free( old );
        }
    }
    return( message );
}
/*
============================================================
*/
void *smr_malloc( statusMessageReporting *smr, size_t size, int zero, char const *forItem, char const *file, int line, char const *function ) {

    void *p = smr_realloc( smr, NULL, size, forItem, file, line, function );
    size_t i;
    char *c;
    long long *l;

    if( ( p != NULL ) && zero ) {
        for( i = 0, l = (long long *) p; i < size / sizeof( long long ); i++, l++ ) *l = 0;
        for( i *= sizeof( long long ), c = (char *) l; i < size; i++, c++ ) *c = 0;
    }

    return( p );
}
/*
============================================================
*/
void *smr_realloc( statusMessageReporting *smr, void *pOld, size_t size, char const *forItem, char const *file, int line, char const *function ) {

    void *p = realloc( pOld, size );

    if( ( p == NULL ) && ( smr != NULL ) ) {
        smr_setReportError( smr, NULL, file, line, function, smr_smrID, smr_codeMemoryAllocating, " smr_realloc: failed to realloc size = %z for variable %s\n", size, forItem );
    }
    return( p );
}
/*
============================================================
*/
void *smr_freeMemory( void **p ) {

    if( p == NULL ) return( NULL );
    if( *p != NULL ) {
        free( *p );
        *p = NULL;
    }
    return( *p );
}
/*
============================================================
*/
char *smr_allocateCopyString( statusMessageReporting *smr, char const *s, char const *forItem, char const *file, int line, char const *function ) {
/*
*   User must free returned string.
*/
    char *c = strdup( s );

    if( c == NULL ) smr_setReportError( smr, NULL, file, line, function, smr_smrID, smr_codeMemoryAllocating, " smr_allocateCopyString: strdup failed for strlen( s ) = %z for variable %s",
    		strlen( s ), forItem );
    return( c );
}
/*
============================================================
*/
char *smr_allocateCopyStringN( statusMessageReporting *smr, char const *s, size_t n, char const *forItem, char const *file, int line, char const *function ) {
/*
*   User must free returned string.
*/
    size_t l = strlen( s );
    char *c;

    if( l > n ) l = n;
    if( ( c = (char *) smr_malloc( smr, l + 1, 0, forItem, file, line, function ) ) != NULL ) {
        strncpy( c, s, n );
        c[l] = 0;
    }
/*
    c = strndup( s, l );        # Not standard on enough systems.
    if( c != NULL ) {
        c[l] = 0; }
    else {
         smr_setReportError( smr, NULL, file, line, function, smr_smrID, smr_codeMemoryAllocating, " smr_allocateCopyStringN: strndup failed for strlen( s ) = %z for variable %s",
            strlen( s ), forItem );
    }
*/
    return( c );
}
