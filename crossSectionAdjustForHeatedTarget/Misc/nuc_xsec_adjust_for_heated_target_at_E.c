/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

/*
*   This routine heats cross section data in a file for a single incident energy E.
*   USAGE:
*       nuc_xsec_adjust_for_heated_target_at_E file Temperature massRatio Energy
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>

#include "crossSectionAdjustForHeatedTarget.h"

#define bufferSize 1024

static double getDouble( const char *str, const char *name );
void PrintMsg( char *fmt, ... );
/*
***************************************************
*/
int main( int argc, char **argv ) {

    int i, n = 0, err;
    char buffer[bufferSize];
    FILE *f;
    double E, cs, *ECs, T, massRatio, energy, f_interpolation = 2e-3, t1, t2, t3, t4;
    crossSectionAdjustForHeatedTarget_info info;
    crossSectionAdjustForHeatedTarget_limit lowerlimit = crossSectionAdjustForHeatedTarget_limit_constant;
    E_cs_heated_point_Info E_cs_Info;
    E_cs_heated_point E_cs_point;

    info.mode = 0;
    info.verbose = 0;

    if( argc != 5 ) PrintMsg( "Need dataFile, Temperature, massRatio and Energy" );
    T = getDouble( argv[2], "temperature" );
    massRatio = getDouble( argv[3], "massRatio" );
    energy = getDouble( argv[4], "energy" );

    if( ( f = fopen( argv[1], "r" ) ) == NULL ) PrintMsg( "Cannot open data file %s", argv[1] );
    while( 1 ) {
        if( fgets( buffer, bufferSize, f ) == NULL ) break;
        n++;
    }
    if( n == 0 ) PrintMsg( "data file is empty" );
    fseek( f, 0, SEEK_SET );
    ECs = (double *) malloc( 2 * n * sizeof( double ) );
    t1 = clock( ) / (double) CLOCKS_PER_SEC;
    for( i = 0; i < n; i++ ) {
        fgets( buffer, bufferSize, f );
#if 0
        if( sscanf( buffer, "%le %le", &E, &cs ) != 2 ) PrintMsg( "converting E and cs to doubles at line %d\n%s", i + 1, buffer );
#else       /* This way is 1.6 times faster. */
        {
        char *e;
        E = strtod( buffer, &e );
        if( ( *e != ' ' ) && ( *e != '\t' ) && ( *e != ',' ) ) PrintMsg( "converting E to double at line %d\n%s%s", i + 1, buffer, e );
        cs = strtod( e, &e );
        if( ( *e != ' ' ) && ( *e != '\t' ) && ( *e != ',' ) && ( *e != '\n' ) ) PrintMsg( "converting cs to double at line %d\n%s%s", i + 1, buffer, e );
        }
#endif
        ECs[2*i] = E;
        ECs[2*i+1] = cs;
    }
    t2 = clock( ) / (double) CLOCKS_PER_SEC;
    fclose( f );

    if( ECs[0] > 1.1e-9 ) lowerlimit = crossSectionAdjustForHeatedTarget_limit_threshold;
    err = crossSectionAdjustForHeatedTarget_init( lowerlimit, crossSectionAdjustForHeatedTarget_limit_constant, &info, massRatio, T, f_interpolation, 
        n, ECs, &E_cs_Info );
    t3 = clock( ) / (double) CLOCKS_PER_SEC;
    if( err < 0 ) PrintMsg( "Error - crossSectionAdjustForHeatedTarget returned %d", err );

    crossSectionAdjustForHeatedTarget_heat_at_E( energy, &E_cs_Info, &E_cs_point );
    t4 = clock( ) / (double) CLOCKS_PER_SEC;
    printf( "%20.12e %15.7e   # Read time = %.3f sec.  Heating setup time = %.3f sec. Heating time = %.3f sec.\n", 
        E_cs_point.E, E_cs_point.cs, t2 - t1, t3 - t2, t4 - t3 );
    free( ECs );
    exit( EXIT_SUCCESS );
}
/*
***************************************************
*/
static double getDouble( const char *str, const char *name ) {

    char *e;
    double d = strtod( str, &e );

    if( *e != 0 ) PrintMsg( "could not convert '%s' to double for %s", str, name );
    return( d );
}
/*
***************************************************
*/
void PrintMsg( char *fmt, ... ) {

    char Str[64 * 1024];
    va_list args;

    va_start( args, fmt );
    vsnprintf( Str, sizeof( Str ), fmt, args );
    fprintf( stderr, "\nError: %s\n", Str );
    va_end( args );
    exit( EXIT_FAILURE );
}
