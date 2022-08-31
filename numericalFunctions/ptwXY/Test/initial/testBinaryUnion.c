/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

/*
    This routine tests the not mutual domain logic. It multiplies or divides two ptwXYPoints (u and v) and checks the return status 
    for the following cases.
      Multiplication tests.
        1) Domain of u greater than that of v at both ends and end points of v contain non-zero y-values.
        2) Same as test 1 except lower v point has zero y-value.
        3) Same as test 1 except upper v point has zero y-value.
        4) Same as test 1 except lower and upper v points have zero y-values.
      Division tests with same domain as 4.
        5) v / u, no safe divide should fail with nfu_divByZero
        6) v / u, safe divide should pass.
        7) u / v, safe divide should fail with nfu_XOutsideDomain.
        8) v / u, same as 6 execpt lower v point has non-zero y-value.
*/
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>

#include <nfut_utilities.h>
#include <ptwXY.h>

static int verbose = 0;

nfu_status getStatus( statusMessageReporting *smr );
void printMsg( char const *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int Test, iarg, echo = 0;
    ptwXYPoints *u, *v, *m;
    double uXY[] = { 1.0, 11.0, 2.66, 1.7, 3.14, 3.12, 3.66, -1.6, 7.66, 1.5 };
    double vXY[] = { 3.732, 1.760, 5.598, 4.349, 5.679, 7.934, 6.459, 4.288, 6.557, 6.787, 6.838, -0.92423, 7.1138, 6.5695 };
    int nuXY = sizeof( uXY ) / ( 2 * sizeof( double ) ), nvXY = sizeof( vXY ) / ( 2 * sizeof( double ) );
    FILE *ff;
    char *fmt = "%10.5f %12.5f\n";
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    nfu_setMemoryDebugMode( 0 );

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            printMsg( "Error %s: unsupported option = '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    if( ( u = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 3, 1e-3, 10, 10, nuXY, uXY, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( v = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 3, 1e-3, 10, 10, nvXY, vXY, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );

    ff = fopen( "curve_u.dat", "w" );
    ptwXY_simpleWrite( u, ff, fmt );
    if( verbose ) ptwXY_simpleWrite( u, stdout, fmt );
    if( verbose ) printf( "\n" );
    fclose( ff );

    ff = fopen( "curve_v.dat", "w" );
    ptwXY_simpleWrite( v, ff, fmt );
    if( verbose ) ptwXY_simpleWrite( v, stdout, fmt );
    if( verbose ) printf( "\n" );
    fclose( ff );

    ff = fopen( "curve_o.dat", "w" );

Test = 1;                               /* This must produce a nfu_domainsNotMutual error. */
#ifdef DEBUG
printf( "Test %d\n", Test );
#endif
    if( ( m = ptwXY_mul2_ptwXY( &smr, u, v ) ) == NULL ) {
        if( getStatus( &smr ) != nfu_domainsNotMutual ) nfut_printSMRError2p( &smr, "Via." );
        smr_release( &smr ); }
    else {
        printMsg( "Test %d multiplication: should have failed with '%s'", Test, nfu_statusMessage( nfu_domainsNotMutual ) );
        ptwXY_free( m );
    }

Test = 2;                               /* This must produce a nfu_domainsNotMutual error. */
#ifdef DEBUG
printf( "Test %d\n", Test );
#endif
    if( ptwXY_setValueAtX( &smr, v, vXY[0], 0. ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) ptwXY_simpleWrite( v, stdout, fmt );
    if( verbose ) printf( "\n" );
    if( ( m = ptwXY_mul2_ptwXY( &smr, u, v ) ) == NULL ) {
        if( getStatus( &smr ) != nfu_domainsNotMutual ) nfut_printSMRError2p( &smr, "Via." );
        smr_release( &smr ); }
    else {
        printMsg( "Test %d multiplication: should have failed with: %s", Test, nfu_statusMessage( nfu_domainsNotMutual ) );
        ptwXY_free( m );
    }

Test = 3;                               /* This must produce a nfu_domainsNotMutual error. */
#ifdef DEBUG
printf( "Test %d\n", Test );
#endif
    if( ptwXY_setValueAtX( &smr, v, vXY[0], 1. ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_setValueAtX( &smr, v, vXY[2 * ( nvXY - 1 )], 0. ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) ptwXY_simpleWrite( v, stdout, fmt );
    if( verbose ) printf( "\n" );
    if( ( m = ptwXY_mul2_ptwXY( &smr, u, v ) ) == NULL ) {
        if( getStatus( &smr ) != nfu_domainsNotMutual ) nfut_printSMRError2p( &smr, "Via." );
        smr_release( &smr ); }
    else {
        printMsg( "Test %d multiplication: should have failed with: %s", Test, nfu_statusMessage( nfu_domainsNotMutual ) );
        ptwXY_free( m );
    }

Test = 4;                               /* This must work. */
#ifdef DEBUG
printf( "Test %d\n", Test );
#endif
    if( ptwXY_setValueAtX( &smr, v, vXY[0], 0. ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) ptwXY_simpleWrite( v, stdout, fmt );
    if( verbose ) printf( "\n" );
    if( ( m = ptwXY_mul2_ptwXY( &smr, u, v ) ) == NULL ) {
        nfut_printSMRError2p( &smr, "Via." ); }
    else {
        ptwXY_simpleWrite( m, ff, fmt );
        fprintf( ff, "\n" );
        if( verbose ) ptwXY_simpleWrite( m, stdout, fmt );
        if( verbose ) printf( "\n" );
        ptwXY_free( m );
    }

Test = 5;                               /* This must cause a nfu_divByZero error. */
#ifdef DEBUG
printf( "Test %d\n", Test );
#endif
    if( ( m = ptwXY_div_ptwXY( &smr, v, u, 0 ) ) == NULL ) {
        if( getStatus( &smr ) != nfu_divByZero ) nfut_printSMRError2p( &smr, "Via." );
        smr_release( &smr ); }
    else {
        printMsg( "Test %d division: should have failed with: %s", Test, nfu_statusMessage( nfu_divByZero ) );
        ptwXY_free( m );
    }

Test = 6;                               /* This must work. */
#ifdef DEBUG
printf( "Test %d\n", Test );
#endif
    if( ( m = ptwXY_div_ptwXY( &smr, v, u, 1 ) ) == NULL ) {
        nfut_printSMRError2p( &smr, "Via." ); }
    else {
        ptwXY_simpleWrite( m, ff, fmt );
        fprintf( ff, "\n" );
        if( verbose ) ptwXY_simpleWrite( m, stdout, fmt );
        if( verbose ) printf( "\n" );
        ptwXY_free( m );
    }

Test = 7;                               /* This must cause a nfu_XOutsideDomain error. */
#ifdef DEBUG
printf( "Test %d\n", Test );
#endif
    if( ( m = ptwXY_div_ptwXY( &smr, u, v, 1 ) ) == NULL ) {
        if( getStatus( &smr ) != nfu_XOutsideDomain ) nfut_printSMRError2p( &smr, "Via." );
        smr_release( &smr ); }
    else {
        printMsg( "Test %d division: should have failed with: %s", Test, nfu_statusMessage( nfu_XOutsideDomain ) );
        ptwXY_free( m );
    }

Test = 8;                               /* This must cause a nfu_domainsNotMutual error. */
#ifdef DEBUG
printf( "Test %d\n", Test );
#endif
    if( ptwXY_setValueAtX( &smr, v, vXY[0], 1. ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ( m = ptwXY_div_ptwXY( &smr, v, u, 1 ) ) == NULL ) {
        if( getStatus( &smr ) != nfu_domainsNotMutual ) nfut_printSMRError2p( &smr, "Via." );
        smr_release( &smr ); }
    else {
        printMsg( "Test %d division: should have failed with: %s", Test, nfu_statusMessage( nfu_domainsNotMutual ) );
        ptwXY_free( m );
    }

    fclose( ff );
    ptwXY_free( u );
    ptwXY_free( v );

    exit( EXIT_SUCCESS );
}
/*
****************************************************************
*/
nfu_status getStatus( statusMessageReporting *smr ) {

    statusMessageReport const *report = smr_firstReport( smr );

    if( report == NULL ) return( nfu_Okay );
    return( smr_getCode( report ) );
}
/*
****************************************************************
*/
void printMsg( char const *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    fprintf( stderr, "\n" );
    va_end( args );
    exit( EXIT_FAILURE );
}
