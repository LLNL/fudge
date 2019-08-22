/*
# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>
*/

/*
*   This routine heats cross section data in a file and writes the output to a file
*   with the same name and the added extension ".heated".
*   USAGE:
*       nuc_xsec_adjust_for_heated_target file Temperature massRatio
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
    double E, cs, *ECs, *E_cs_out, T, massRatio, f_interpolation = 2e-3, t1, t2, t3;
    crossSectionAdjustForHeatedTarget_info info;
    crossSectionAdjustForHeatedTarget_limit lowerlimit = crossSectionAdjustForHeatedTarget_limit_constant;

    info.mode = 0;
    info.verbose = 0;

    if( argc != 4 ) PrintMsg( "Need dataFile, Temperature and massRatio" );
    T = getDouble( argv[2], "temperature" );
    massRatio = getDouble( argv[3], "massRatio" );

    if( ( f = fopen( argv[1], "r" ) ) == NULL ) PrintMsg( "Cannot open data file %s", argv[1] );
    while( 1 ) {
        if( fgets( buffer, bufferSize, f ) == NULL ) break;
        if( buffer[0] == '#' ) continue;
        n++;
    }
    if( n == 0 ) PrintMsg( "data file is empty" );
    fseek( f, 0, SEEK_SET );
    ECs = (double *) malloc( 2 * n * sizeof( double ) );
    t1 = clock( ) / (double) CLOCKS_PER_SEC;
    for( i = 0; i < n; ) {
        fgets( buffer, bufferSize, f );
        if( buffer[0] == '#' ) continue;
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
        i++;
    }
    t2 = clock( ) / (double) CLOCKS_PER_SEC;
    fclose( f );

    if( ECs[0] > 1.1e-9 ) lowerlimit = crossSectionAdjustForHeatedTarget_limit_threshold;
    err = crossSectionAdjustForHeatedTarget( lowerlimit, crossSectionAdjustForHeatedTarget_limit_constant, &info, 1e-11, massRatio, T, f_interpolation, 
        n, ECs, &E_cs_out );
    t3 = clock( ) / (double) CLOCKS_PER_SEC;
    if( err < 0 ) PrintMsg( "Error - crossSectionAdjustForHeatedTarget returned %d", err );

    strcpy( buffer, argv[1] );
    strcat( buffer, ".heated" );
    if( ( f = fopen( buffer, "w" ) ) == NULL ) PrintMsg( "Cannot open output file %s", buffer );
    fprintf( f, "# Temperature = %g  massRatio = %g  Read time = %.3f sec.  Heating time = %.3f sec.\n", T, massRatio, t2 - t1, t3 - t2 );
    for( i = 0; i < err; i++ ) fprintf( f, " %13.7E %12.5E\n", E_cs_out[2 * i], E_cs_out[2 * i + 1] );
    fclose( f );

    free( ECs );
    free( E_cs_out );
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
