/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <nf_utilities.h>
#include <nfut_utilities.h>

#define nBins 17
#define nStringSize 16

static int verbose = 0;
static long numberOfCompareSamples = 1000 * 1000;

int compare( statusMessageReporting *smr, long numberOfCompareSamples );
void parseSingle( char const *single );
void timing( statusMessageReporting *smr );
void timing2( statusMessageReporting *smr, long numberOfValuesToGenerate );
char *generateValues( statusMessageReporting *smr, long numberOfValuesToGenerate, int32_t **values );
void addToString( char **dest, char *src );
void printUsage( );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, doTiming = 0, errCount = 0;
    statusMessageReporting smr;
    char *single = NULL;

    smr_initialize( &smr, smr_status_Ok );

    for( iarg = 1; iarg < argc; ++iarg ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else if( strcmp( "-h", argv[iarg] ) == 0 ) {
            printUsage( ); }
        else if( strcmp( "-c", argv[iarg] ) == 0 ) {
            ++iarg;
            if( iarg == argc ) nfu_printErrorMsg( "ERROR %s: -c options needs float string", __FILE__ );
            numberOfCompareSamples = nfut_charToLong( &smr, "numberOfCompareSamples", argv[iarg] ); }
        else if( strcmp( "-s", argv[iarg] ) == 0 ) {
            ++iarg;
            if( iarg == argc ) nfu_printErrorMsg( "ERROR %s: -s options needs float string", __FILE__ );
            single = argv[iarg]; }
        else if( strcmp( "-t", argv[iarg] ) == 0 ) {
            doTiming = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    if( single != NULL ) parseSingle( single );
    if( doTiming > 0 ) timing( &smr );

    errCount = compare( &smr, numberOfCompareSamples );

    exit( errCount );
}
/*
************************************************************
*/
int compare( statusMessageReporting *smr, long numberOfCompareSamples ) {

    int errCount = 0;
    char *stringOfValues;
    char *endCharacter;
    long i1;
    int64_t numberConverted;
    int32_t *values1, *values = NULL;

    stringOfValues = generateValues( smr, numberOfCompareSamples, &values );

    values1 = nfu_stringToListOfInt32s( smr, stringOfValues, ' ', &numberConverted, &endCharacter );
    if( values1 == NULL ) nfut_printSMRErrorExit2p( smr, "compare" );
    if( verbose ) {
        long length = endCharacter - stringOfValues;
        long size = strlen( endCharacter );

        printf( "(1)\n" );
        printf( "   length = %ld  size = %ld\n", length, size );
        printf( "   endCharacter = <%s>\n", endCharacter );
    }
    if( (long) numberConverted != numberOfCompareSamples )
        nfut_printSMRErrorExit2( smr, "nfu_stringToListOfInt32s: numberConverted = %ld != numberOfCompareSamples = %ld", numberConverted, numberOfCompareSamples );

    for( i1 = 0; i1 < numberOfCompareSamples; ++i1 ) {
        if( values[i1] != values1[i1] ) {
            ++errCount;
            if( verbose ) {
                if( values[i1] != values1[i1] ) printf( "%9ld  %9d %9d\n", i1, values[i1], values1[i1] );
            }
        }
    }
    if( errCount != 0 ) printf( "errCount = %9d\n", errCount );

    free( stringOfValues );
    free( values );
    free( values1 );

    return( errCount );
}
/*
************************************************************
*/
void parseSingle( char const *single ) {

    char *endCharacter;
    int32_t value;

    nfu_stringToInt32( NULL, single, &endCharacter, &value );

    printf( "string = <%s>   value = <%d>   endCharacters = <%s>\n", single, value, endCharacter );

    exit( 0 );
}
/*
************************************************************
*/
void timing( statusMessageReporting *smr ) {

    timing2( smr,      10 );
    timing2( smr,     100 );
    timing2( smr,    1000 );
    timing2( smr,   10000 );
    timing2( smr,  100000 );
    timing2( smr, 1000000 );
    exit( 0 );
}
/*
************************************************************
*/
void timing2( statusMessageReporting *smr, long numberOfValuesToGenerate ) {

    long i1, loops = 1 * 1000 * 1000 / numberOfValuesToGenerate;
    int64_t numberConverted;
    char *stringOfValues;
    char *endCharacter;
    int32_t *values;
    clock_t time0;
    double cpuTime, cpuGenerationTime, rate;

    if( loops == 0 ) loops = 1;

    time0 = clock( );
    stringOfValues = generateValues( smr, numberOfValuesToGenerate, NULL );
    cpuGenerationTime = ( clock( ) - time0 ) / ( (double) CLOCKS_PER_SEC );

    time0 = clock( );
    for( i1 = 0; i1 < loops; ++i1 ) {
        values = nfu_stringToListOfInt32s( smr, stringOfValues, ' ', &numberConverted, &endCharacter );
        if( values == NULL ) nfut_printSMRErrorExit2p( smr, "" );
        free( values );
    }
    cpuTime = ( clock( ) - time0 ) / ( (double) CLOCKS_PER_SEC ) / loops;
    rate = numberOfValuesToGenerate / cpuTime;
    printf( "cpu time for %8ld strings to int32s = %9.3e (%.2e/s) --- values generation time = %9.3e\n", numberOfValuesToGenerate, cpuTime, rate, cpuGenerationTime );

    free( stringOfValues );
}
/*
************************************************************
*/
char *generateValues( statusMessageReporting *smr, long numberOfValuesToGenerate, int32_t **values ) {

    long i1;
    long stringSize = numberOfValuesToGenerate;
    char *stringOfValues, *p1;

    stringSize *= nStringSize;
    if( ( stringOfValues = malloc( stringSize ) ) == NULL ) nfut_printSMRErrorExit2( smr, "malloc of %ld bytes failed for stringOfValues", stringSize );
    if( values != NULL ) {
        long size = numberOfValuesToGenerate * sizeof( int32_t );

        *values = malloc( size );
        if( *values == NULL ) nfut_printSMRErrorExit2( smr, "malloc of %ld bytes failed for values", size );
    }

    p1 = stringOfValues;
    addToString( &p1, "   " );

    for( i1 = 0; i1 < numberOfValuesToGenerate; ++i1 ) {
        char str[nStringSize+1];
        int32_t value = (int32_t) mrand48( );

        if( values != NULL ) (*values)[i1] = value;
        sprintf( str, " %d\n", value );
        addToString( &p1, str );
    }

    return( stringOfValues );
}
/*
************************************************************
*/
void addToString( char **dest, char *src ) {

    for( ; *src != 0; ++src, ++(*dest) ) **dest = *src;
    **dest = 0;
}
/*
************************************************************
*/
void printUsage( ) {

    printf( "\nUSAGE:\n" );
    printf( "    stringToListOfInt32s [-h] [-e] [-v] [-s VALUE] [-c SAMPLES] [-t]\n" );
    printf( "\n" );
    printf( "    -h             print this message.\n" );
    printf( "    -e             print file name.\n" );
    printf( "    -v             verbose flag.\n" );
    printf( "    -s VALUE       print the results of converting VALUE with system stdtod and nf_stdtod.\n" );
    printf( "    -c SAMPLES     compare system stdtod to nf_stdtod with SAMPLES samples. Default option with %ld samples.\n", numberOfCompareSamples );
    printf( "    -t             run timing.\n" );

    exit( EXIT_SUCCESS );
}
