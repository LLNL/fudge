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
#define nStringSize 64

static int verbose = 0;
static long numberOfCompareSamples = 1000 * 1000;

int compare( statusMessageReporting *smr, long numberOfCompareSamples );
void parseSingle( char const *single );
void timing( statusMessageReporting *smr );
void timing2( statusMessageReporting *smr, int useSystem_strtod );
void timing3( statusMessageReporting *smr, long numberOfValuesToGenerate, int useSystem_strtod );
char *generateValues( statusMessageReporting *smr, long numberOfValuesToGenerate, double **values );
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
    int64_t numberConverted;
    long i1, bins[nBins+1], bins2[nBins+1];
    double *values1, *values2, *values = NULL;

    for( i1 = 0; i1 <= nBins; ++i1 ) bins[i1] = bins2[i1] = 0;

    stringOfValues = generateValues( smr, numberOfCompareSamples, &values );

    values1 = nfu_stringToListOfDoubles( smr, stringOfValues, ' ', &numberConverted, &endCharacter, 1 );
    if( values1 == NULL ) nfut_printSMRErrorExit2p( smr, "system strtod" );
    if( verbose ) {
        long length = endCharacter - stringOfValues;
        long size = strlen( endCharacter );

        printf( "(1)\n" );
        printf( "   length = %ld  size = %ld\n", length, size );
        printf( "   endCharacter = <%s>\n", endCharacter );
    }
    if( (long) numberConverted != numberOfCompareSamples )
        nfut_printSMRErrorExit2( smr, "stdtod: numberConverted = %ld != numberOfCompareSamples = %ld", numberConverted, numberOfCompareSamples );

    values2 = nfu_stringToListOfDoubles( smr, stringOfValues, ' ', &numberConverted, &endCharacter, 0 );
    if( values2 == NULL ) nfut_printSMRErrorExit2p( smr, "system strtod" );
    if( verbose ) {
        long length = endCharacter - stringOfValues;
        long size = strlen( endCharacter );

        printf( "(2)\n" );
        printf( "   length = %ld  size = %ld\n", length, size );
        printf( "   endCharacter = <%s>\n", endCharacter );
    }
    if( (long) numberConverted != numberOfCompareSamples ) {
        for( i1 = 0; i1 < numberOfCompareSamples; ++i1 ) {
            if( fabs( values[i1] - values2[i1] ) > 5e-3 * fabs( values[i1] ) ) {
                printf( " %10ld %25.17e %25.17e\n", i1, values[i1], values2[i1] );
                break;
            }
        }
        printf( "\n" );
        for( i1 = numberConverted - 3; i1 < numberOfCompareSamples; ++i1 ) printf( " %8ld %25.17e %25.17e %25.17e\n", i1, values[i1], values1[i1], values2[i1] );
        nfut_printSMRErrorExit2( smr, "nf_strtod: numberConverted = %ld != numberOfCompareSamples = %ld", numberConverted, numberOfCompareSamples );
    }

    for( i1 = 0; i1 < numberOfCompareSamples; ++i1 ) {
        int bin;
        double value1 = values1[i1], value2 = values2[i1];
        double diff = value1 - value2, relDiff;
        double max = fabs( value1 );

        if( diff == 0.0 ) continue;
        if( fabs( value2 ) > max ) max = fabs( value2 );

        relDiff = fabs( diff / max );
        if( relDiff > 2e-15 ) ++errCount;

        bin = -log10( relDiff ) + 1;
        if( bin < 0 ) bin = 0;
        if( bin > nBins ) bin = nBins;
        ++bins[bin];
        if( verbose ) {
            if( relDiff > 2e-15 ) printf( "bin = %d  %25.17e %25.17e %11.3e %11.3e %11.3e\n", bin, value1, value2, diff, relDiff, -log10( relDiff ) );
        }

        for( bin = 0; bin < nBins; ++bin ) {
            if( relDiff < 2.2e-16 ) break;
            relDiff /= 2.0;
        }
        ++bins2[bin];
    }

    if( verbose > 0 ) {
        printf( " bin     count  ~bit err\n" );
        for( i1 = 0; i1 <= nBins; ++i1 ) printf( " %3ld  %8ld  %8ld\n", i1, bins[i1], bins2[i1] );
    }

    free( stringOfValues );
    free( values );
    free( values1 );
    free( values2 );

    return( errCount );
}
/*
************************************************************
*/
void parseSingle( char const *single ) {

    char *endCharacter1;
    char *endCharacter2;
    double value1 = strtod( single, &endCharacter1 );

    printf( "string = <%s>   value = <%.17e>   endCharacters = <%s>\n", single, value1, endCharacter1 );

    double value2 = nf_strtod( single, &endCharacter2 );
    printf( "string = <%s>   value = <%.17e>   endCharacters = <%s>\n", single, value2, endCharacter2 );

    if( value1 != value2 ) {
        double diff = value1 - value2;
        double max = fabs( value1 );

        if( value1 == 0.0 ) max = fabs( value2 );
        double relDiff = diff / max;

        printf( "    diff = %.3e  relDiff = %.3e\n", diff, relDiff );
    }

    exit( 0 );
}
/*
************************************************************
*/
void timing( statusMessageReporting *smr ) {

    timing2( smr, 0 );
    timing2( smr, 1 );
    exit( 0 );
}
/*
************************************************************
*/
void timing2( statusMessageReporting *smr, int useSystem_strtod ) {

    if( useSystem_strtod ) {
        printf( "Using system strtod\n" ); }
    else {
        printf( "Using nf_strtod\n" );
    }

    timing3( smr,      10, useSystem_strtod );
    timing3( smr,     100, useSystem_strtod );
    timing3( smr,    1000, useSystem_strtod );
    timing3( smr,   10000, useSystem_strtod );
    timing3( smr,  100000, useSystem_strtod );
    timing3( smr, 1000000, useSystem_strtod );
}
/*
************************************************************
*/
void timing3( statusMessageReporting *smr, long numberOfValuesToGenerate, int useSystem_strtod ) {

    long i1, loops = 1 * 1000 * 1000 / numberOfValuesToGenerate;
    int64_t numberConverted;
    char *stringOfValues;
    char *endCharacter;
    double *values;
    clock_t time0;
    double cpuTime, cpuGenerationTime, rate;

    if( loops == 0 ) loops = 1;

    time0 = clock( );
    stringOfValues = generateValues( smr, numberOfValuesToGenerate, NULL );
    cpuGenerationTime = ( clock( ) - time0 ) / ( (double) CLOCKS_PER_SEC );

    time0 = clock( );
    for( i1 = 0; i1 < loops; ++i1 ) {
        values = nfu_stringToListOfDoubles( smr, stringOfValues, ' ', &numberConverted, &endCharacter, useSystem_strtod );
        if( values == NULL ) nfut_printSMRErrorExit2p( smr, "" );
        free( values );
    }
    cpuTime = ( clock( ) - time0 ) / ( (double) CLOCKS_PER_SEC ) / loops;
    rate = numberOfValuesToGenerate / cpuTime;
    printf( "    cpu time for %8ld strings to doubles = %9.3e (%.2e/s) --- values generation time = %9.3e\n", numberOfValuesToGenerate, cpuTime, rate, cpuGenerationTime );

    free( stringOfValues );
}
/*
************************************************************
*/
char *generateValues( statusMessageReporting *smr, long numberOfValuesToGenerate, double **values ) {

    long i1;
    long stringSize = numberOfValuesToGenerate;
    char *stringOfValues, *p1;

    stringSize *= nStringSize;
    if( ( stringOfValues = malloc( stringSize ) ) == NULL ) 
        nfut_printSMRErrorExit2( smr, "malloc of %ld bytes failed for stringOfValues", stringSize );
    if( values != NULL ) {
        long size = numberOfValuesToGenerate * sizeof( double );

        *values = malloc( size );
        if( *values == NULL ) nfut_printSMRErrorExit2( smr, "malloc of %ld bytes failed for values", size );
    }
    p1 = stringOfValues;
    addToString( &p1, "   " );

    for( i1 = 0; i1 < numberOfValuesToGenerate; ++i1 ) {
        int eFormat = ( drand48( ) < 0.5 ) ? 0 : 1;
        int precision = (int) ( 3 + 15 * drand48( ) );
        char fmt[12], str[nStringSize+1];

        double value;
        double power = -25.0 * log( drand48( ) );       /* For -25.0, will product a power > 300 about 4 in every million samples. */

        while( fabs( power ) > 308 ) power = 10.0 * drand48( );
        if( drand48( ) < 0.5 ) power = -power;
        value = drand48( ) * pow( 10.0, power );
        if( drand48( ) < 0.5 ) value = -value;

        if( values != NULL ) (*values)[i1] = value;
        sprintf( fmt, " %%.%d%c", precision, ( eFormat ? 'e' : 'g' ) );
        sprintf( str, fmt, value );
        addToString( &p1, str );
        if( drand48( ) < 0.5 ) addToString( &p1, " " );
        if( drand48( ) < 0.5 ) addToString( &p1, " " );
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
    printf( "    stringToListOfDoubles [-h] [-e] [-v] [-s VALUE] [-c SAMPLES] [-t]\n" );
    printf( "\n" );
    printf( "    -h             print this message.\n" );
    printf( "    -e             print file name.\n" );
    printf( "    -v             verbose flag.\n" );
    printf( "    -s VALUE       print the results of converting VALUE with system stdtod and nf_stdtod.\n" );
    printf( "    -c SAMPLES     compare system stdtod to nf_stdtod with SAMPLES samples. Default option with %ld samples.\n", numberOfCompareSamples );
    printf( "    -t             run timing.\n" );

    exit( EXIT_SUCCESS );
}
