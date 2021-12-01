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

nfu_status chargedParticleGetValue( void *argList, double x, double *y, double x1, double y1, double x2, double y2 );
void printMsg( const char *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, nPoints;
    ptwXYPoints *ptwXY, *ptwXY2, *ptwXY3;
    double accuracy = 1e-3;
    double xys[] = { 1, 1, 10, 10 };
    statusMessageReporting smr;

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
    
    nPoints = sizeof( xys ) / sizeof( xys[0] ) / 2;

    if( ( ptwXY2 = ptwXY_create( &smr, ptwXY_interpolationOther, "charged-particle", 5, accuracy, 10, 10, nPoints, xys, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );

    if( ( ptwXY = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 5, accuracy, 10, 10, nPoints, xys, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "linear-linear string = <%s>\n", ptwXY_getInterpolationString( ptwXY ) );
    if( ptwXY_copy( &smr, ptwXY2, ptwXY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "string = <%s>\n", ptwXY_getInterpolationString( ptwXY2 ) );
    if( ( ptwXY3 = ptwXY_clone( &smr, ptwXY2 ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "string = <%s>\n", ptwXY_getInterpolationString( ptwXY3 ) );
    ptwXY_free( ptwXY );
    ptwXY_free( ptwXY3 );

    if( verbose ) printf( "\n" );
    if( ( ptwXY = ptwXY_create( &smr, ptwXY_interpolationLogLin, NULL, 5, accuracy, 10, 10, nPoints, xys, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "linear-log string = <%s>\n", ptwXY_getInterpolationString( ptwXY ) );
    if( ptwXY_copy( &smr, ptwXY2, ptwXY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "string = <%s>\n", ptwXY_getInterpolationString( ptwXY2 ) );
    if( ( ptwXY3 = ptwXY_clone( &smr, ptwXY2 ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "string = <%s>\n", ptwXY_getInterpolationString( ptwXY3 ) );
    ptwXY_free( ptwXY );
    ptwXY_free( ptwXY3 );

    if( verbose ) printf( "\n" );
    if( ( ptwXY = ptwXY_create( &smr, ptwXY_interpolationLinLog, NULL, 5, accuracy, 10, 10, nPoints, xys, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "log-linear string = <%s>\n", ptwXY_getInterpolationString( ptwXY ) );
    if( ptwXY_copy( &smr, ptwXY2, ptwXY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "string = <%s>\n", ptwXY_getInterpolationString( ptwXY2 ) );
    if( ( ptwXY3 = ptwXY_clone( &smr, ptwXY2 ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "string = <%s>\n", ptwXY_getInterpolationString( ptwXY3 ) );
    ptwXY_free( ptwXY );
    ptwXY_free( ptwXY3 );

    if( verbose ) printf( "\n" );
    if( ( ptwXY = ptwXY_create( &smr, ptwXY_interpolationLogLog, NULL, 5, accuracy, 10, 10, nPoints, xys, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "log-log string = <%s>\n", ptwXY_getInterpolationString( ptwXY ) );
    if( ptwXY_copy( &smr, ptwXY2, ptwXY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "string = <%s>\n", ptwXY_getInterpolationString( ptwXY2 ) );
    if( ( ptwXY3 = ptwXY_clone( &smr, ptwXY2 ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "string = <%s>\n", ptwXY_getInterpolationString( ptwXY3 ) );
    ptwXY_free( ptwXY );
    ptwXY_free( ptwXY3 );

    if( verbose ) printf( "\n" );
    if( ( ptwXY = ptwXY_create( &smr, ptwXY_interpolationFlat, NULL, 5, accuracy, 10, 10, nPoints, xys, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "flat string = <%s>\n", ptwXY_getInterpolationString( ptwXY ) );
    if( ptwXY_copy( &smr, ptwXY2, ptwXY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "string = <%s>\n", ptwXY_getInterpolationString( ptwXY2 ) );
    if( ( ptwXY3 = ptwXY_clone( &smr, ptwXY2 ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "string = <%s>\n", ptwXY_getInterpolationString( ptwXY3 ) );
    ptwXY_free( ptwXY );
    ptwXY_free( ptwXY3 );

    if( verbose ) printf( "\n" );
    if( ( ptwXY = ptwXY_create( &smr, ptwXY_interpolationOther, "charged-particle", 5, accuracy, 10, 10, nPoints, xys, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "other string = <%s>\n", ptwXY_getInterpolationString( ptwXY ) );
    if( ptwXY_copy( &smr, ptwXY2, ptwXY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "string = <%s>\n", ptwXY_getInterpolationString( ptwXY2 ) );
    if( ( ptwXY3 = ptwXY_clone( &smr, ptwXY2 ) ) == NULL ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) printf( "string = <%s>\n", ptwXY_getInterpolationString( ptwXY3 ) );
    ptwXY_free( ptwXY );
    ptwXY_free( ptwXY3 );

    ptwXY_free( ptwXY2 );

    exit( EXIT_SUCCESS );
}
/*
****************************************************************
*/
nfu_status chargedParticleGetValue( void *argList, double x, double *y, double x1, double y1, double x2, double y2 ) {

    double A, B, T = *((double *) argList);

    B = log( x2 * y2 / ( x1 * y1 ) ) / ( 1. / sqrt( x1 - T ) - 1. / sqrt( x2 - T ) );
    A = x1 * y1 * exp( B / sqrt( x1 - T ) );
    *y = A * exp( - B / sqrt( x - T ) ) / x;
    return( nfu_Okay );
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
