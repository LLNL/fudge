/*
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>

static int verbose = 0;
static char *fmt = "%10.5f %12.5f\n";

static int thicken( statusMessageReporting *smr, ptwXYPoints *u, int length, int sectionSubdivideMax, double dx, double fx, char *limitStr );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    double uXY[] = { 1.0, 3.0, 10, 30. };
    ptwXYPoints *u;
    int nuXY = sizeof( uXY ) / ( 2 * sizeof( double ) ), iarg, errCount = 0, echo = 0;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfu_printErrorMsg( "Error %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    if( ( u = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 6, 1e-3, 10, 10, nuXY, uXY, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    if( verbose ) {
        fprintf( stdout, "# length = %d\n", (int) u->length );
        ptwXY_simpleWrite( u, stdout, fmt );
        fprintf( stdout, "\n\n" );
    }

    errCount += thicken( &smr, u, 11,  10, 0.4, 1.1, "sectionSubdivideMax" );         /* Limit is sectionSubdivideMax. */
    errCount += thicken( &smr, u, 24, 100, 0.4, 11., "dx" );                          /* Limit is dx. */
    errCount += thicken( &smr, u, 26, 100, 1.4, 1.1, "fx" );                          /* Limit is fx. */
    errCount += thicken( &smr, u, 32, 100, 0.4, 1.1, "dx and fx" );                   /* Limit is dx. */

    uXY[2] = 1.0001;
    uXY[3] = 9.;
    if( ptwXY_setXYData( &smr, u, nuXY, uXY ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    u->interpolation = ptwXY_interpolationLogLin;
    errCount += thicken( &smr, u, 101, 100, 0.0000001, 1.0000001, "log-log" );        /* Limit is sectionSubdivideMax. */

    ptwXY_free( u );

    exit( errCount );
}
/*
************************************************************
*/
static int thicken( statusMessageReporting *smr, ptwXYPoints *u, int length, int sectionSubdivideMax, double dx, double fx, char *limitStr ) {

    int i, errs = 0, errCount;
    ptwXYPoints *v;
    ptwXYPoint *point;

    if( verbose ) printf( "# limitStr = '%s': sectionSubdivideMax = %3d  dx = %8.5f  fx = %8.5f\n", limitStr, sectionSubdivideMax, dx, fx );
    if( ( v = ptwXY_clone( smr, u ) ) == NULL ) nfut_printSMRErrorExit2p( smr, "Via." );

    if( ptwXY_thicken( smr, v, sectionSubdivideMax, dx, fx ) != nfu_Okay ) nfut_printSMRErrorExit2p( smr, "Via." );
    if( ptwXY_length( smr, v ) != length ) {
        nfu_printMsg( "Error %s: ptwXY_length( smr, v ) = %d != %d for limit '%s'", __FILE__, (int) ptwXY_length( smr, v ), length, limitStr );
        errs++;
    }
    if( verbose ) {
        double x1, x2;

        printf( "# length = %d\n", (int) v->length );
        for( i = 0; i < (int) ptwXY_length( smr, v ); i++ ) {
            point = ptwXY_getPointAtIndex_Unsafely( v, i );
            x2 = point->x;
            printf( "%10.8f %12.8f", x2, point->y );
            if( i != 0 ) printf( "%12.5e %15.8e", x2 - x1, x2 / x1 );
            printf( "\n" );
            x1 = x2;
        }
    }
    if( verbose ) fprintf( stdout, "\n\n" );
    if( ( errCount = nfu_ptwXY_cmp( u, v, verbose, 1e-15 ) ) ) {
        errs += abs( errCount );
        nfu_printMsg( "Error %s: nfu_ptwXY_cmp found %d differences for limit '%s'", __FILE__, errCount, limitStr );
    }

    ptwXY_free( v );

    return( errs );
}
