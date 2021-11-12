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
#include <math.h>
#include <stdarg.h>

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <ptwXY_utilities.h>

static int verbose = 0;

int interpolationCheck( ptwXY_interpolation interpolation );
void printMsg( const char *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errs = 0;
    statusMessageReporting smr;
    char const *other;
    ptwXY_interpolation interpolation;

    smr_initialize( &smr, smr_status_Ok );

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else {
            printMsg( "Error %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) fprintf( stderr, "%s\n", __FILE__ );

    errs += interpolationCheck( ptwXY_interpolationLinLin );
    errs += interpolationCheck( ptwXY_interpolationLogLin );
    errs += interpolationCheck( ptwXY_interpolationLinLog );
    errs += interpolationCheck( ptwXY_interpolationLogLog );
    errs += interpolationCheck( ptwXY_interpolationFlat );

    other = ptwXY_interpolationToString( ptwXY_interpolationOther );
    if( other != NULL ) errs += 1;
    if( verbose ) printf( "other = %p\n", other );

    interpolation = ptwXY_stringToInterpolation( "hi" );
    if( interpolation != ptwXY_interpolationOther ) errs += 1;
    if( verbose ) printf( "'hi' = %d\n", interpolation );

    exit( errs );
}
/*
****************************************************************
*/
int interpolationCheck( ptwXY_interpolation interpolation ) {

    int errs = 0;
    char const *interpolationString = ptwXY_interpolationToString( interpolation );

    ptwXY_interpolation interpolation2 = ptwXY_stringToInterpolation( interpolationString );
    if( interpolation != interpolation2 ) errs += 1;

    if( verbose ) printf( "%d --> '%s' --> %d\n", interpolation, interpolationString, interpolation2 );
    return( errs );
}
/*
****************************************************************
*/
void printMsg( const char *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    fprintf( stderr, "\n" );
    va_end( args );
    exit( EXIT_FAILURE );
}
