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

#include <ptwXY.h>
#include <nfut_utilities.h>
#include <ptwXY_utilities.h>

static int verbose = 0;
static char *fmtXY = "%25.15e %25.15e\n";
static FILE *infoF;

void printMsg( const char *fmt, ... );
static void printIfVerbose( ptwXYPoints *data );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, nXYs;
    ptwXYPoints *pFlat, *pLinear;
    double accuracy = 1e-3, lowerEps = 1e-6, upperEps = 1e-6, xys[] = { 
        0.00000000e+00, 2.44941900e-08, 3.00000000e+05, 2.95502600e-08, 5.00000000e+05, 1.60509700e-07, 7.00000000e+05, 7.29768400e-07,
        9.00000000e+05, 7.16263300e-08, 1.10000000e+06, 3.13305400e-07, 1.30000000e+06, 9.24667600e-08, 1.50000000e+06, 1.33916500e-07,
        1.70000000e+06, 2.15220800e-07, 1.90000000e+06, 4.09025600e-07, 2.10000000e+06, 7.94523400e-08, 2.30000000e+06, 1.26046700e-07,
        2.50000000e+06, 9.74776700e-08, 2.70000000e+06, 1.82226500e-07, 2.90000000e+06, 1.40790400e-07, 3.10000000e+06, 1.26457400e-07,
        3.30000000e+06, 1.11357400e-07, 3.50000000e+06, 1.41709800e-07, 3.70000000e+06, 1.34958900e-07, 3.90000000e+06, 1.28301200e-07,
        4.10000000e+06, 1.21774200e-07, 4.30000000e+06, 1.15381700e-07, 4.50000000e+06, 1.09135600e-07, 4.70000000e+06, 1.03051700e-07,
        4.90000000e+06, 9.71452600e-08, 5.10000000e+06, 9.14627000e-08, 5.30000000e+06, 8.59984900e-08, 5.50000000e+06, 8.06598500e-08,
        5.70000000e+06, 7.55374200e-08, 5.90000000e+06, 7.06618500e-08, 6.10000000e+06, 2.96670100e-08, 6.30000000e+06, 5.85612600e-08,
        6.50000000e+06, 7.42873500e-08, 6.70000000e+06, 5.47748800e-08, 6.90000000e+06, 1.03412200e-07, 7.10000000e+06, 5.89404900e-09,
        7.30000000e+06, 4.59737500e-09, 7.50000000e+06, 7.45342100e-08, 7.70000000e+06, 2.55224600e-09, 7.90000000e+06, 2.90521500e-08,
        8.10000000e+06, 1.19440100e-09, 8.30000000e+06, 3.76551100e-08, 8.50000000e+06, 3.44437000e-10, 8.70000000e+06, 1.11754500e-07,
        8.90000000e+06, 5.53869000e-13, 9.10000000e+06, 8.23249000e-14, 9.30000000e+06, 7.15863000e-15, 9.49999990e+06, 1.47460800e-07,
        9.50000000e+06, 0.00000000e+00 };
    ptwXY_interpolation interpolation = ptwXY_interpolationFlat;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    infoF = stdout;

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

    nfu_setMemoryDebugMode( 0 );
    nXYs = sizeof( xys ) / ( 2 * sizeof( xys[0] ) );
    if( ( pFlat = ptwXY_create( &smr, interpolation, NULL, 5, accuracy, 10, 10, nXYs, xys, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    printIfVerbose( pFlat );

    if( ( pLinear = ptwXY_flatInterpolationToLinear( &smr, pFlat, lowerEps, upperEps ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    printIfVerbose( pLinear );

    ptwXY_free( pFlat );
    ptwXY_free( pLinear );

    exit( EXIT_SUCCESS );
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
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    fprintf( infoF, "# length = %d\n", (int) ptwXY_length( NULL, data ) );
    ptwXY_simpleWrite( data, infoF, fmtXY );
    fprintf( infoF, "\n\n" );
}
