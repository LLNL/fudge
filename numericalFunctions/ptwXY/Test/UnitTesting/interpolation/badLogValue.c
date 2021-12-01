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

static int verbose = 0;

static int loop( statusMessageReporting *smr, int n1, double *xys, nfu_status *statuses );
void printMsg( const char *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errCounts = 0;
    double xys1[4] = { 1.0, 1.0, 10.0, 100.0 };
    nfu_status status1[4] = { nfu_Okay, nfu_Okay, nfu_Okay, nfu_Okay };
    double xys2[4] = { 1.0, 0.0, 10.0, 1.0 };
    nfu_status status2[4] = { nfu_Okay, nfu_badLogValue, nfu_Okay, nfu_badLogValue };
    double xys3[4] = { -1.0, 1.0, 10.0, 1.0 };
    nfu_status status3[4] = { nfu_Okay, nfu_Okay, nfu_badLogValue, nfu_badLogValue };
    double xys4[4] = { -1.0, 1.0, 10.0, -1.0 };
    nfu_status status4[4] = { nfu_Okay, nfu_badLogValue, nfu_badLogValue, nfu_badLogValue };
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

    errCounts += loop( &smr, 2, xys1, status1 );
    errCounts += loop( &smr, 2, xys2, status2 );
    errCounts += loop( &smr, 2, xys3, status3 );
    errCounts += loop( &smr, 2, xys4, status4 );
    exit( errCounts );
}
/*
****************************************************************
*/
static int loop( statusMessageReporting *smr, int n1, double *xys, nfu_status *statuses ) {

    int i1, errCounts = 0;
    double accuracy = 1e-3;
    ptwXYPoints *p1, *p2;
    nfu_status status;
    statusMessageReport const *report;

    if( ( p1 = ptwXY_create( smr, ptwXY_interpolationLogLog, NULL, 5, accuracy, 10, 10,    2, xys, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( smr, "Via." );

    for( i1 = ptwXY_interpolationLinLin; i1 <= ptwXY_interpolationLogLog; ++i1 ) {
        p1->interpolation = (ptwXY_interpolation) i1;
        p2 = ptwXY_toOtherInterpolation( smr, p1, ptwXY_interpolationLinLin, accuracy );
        ptwXY_free( p2 );
        report = smr_firstReport( smr );
        if( report == NULL ) {
            status = nfu_Okay; }
        else {
            status = smr_getCode( report );
        }
        if( status != statuses[i1] ) ++errCounts;
        smr_release( smr );
    }
    ptwXY_free( p1 );
    return( errCounts );
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
