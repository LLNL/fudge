#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include <ptwXY.h>

static int verbose = 0;

static int loop( int n1, double *xys, nfu_status *statuses );
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

    errCounts += loop( 2, xys1, status1 );
    errCounts += loop( 2, xys2, status2 );
    errCounts += loop( 2, xys3, status3 );
    errCounts += loop( 2, xys4, status4 );
    exit( errCounts );
}
/*
****************************************************************
*/
static int loop( int n1, double *xys, nfu_status *statuses ) {

    int i1, errCounts = 0;
    double accuracy = 1e-3;
    ptwXYPoints *p1, *p2;
    nfu_status status;

    p1 = ptwXY_create( ptwXY_interpolationLogLog, NULL, 5, accuracy, 10, 10,    2, xys, &status, 0 );

    for( i1 = ptwXY_interpolationLinLin; i1 <= ptwXY_interpolationLogLog; ++i1 ) {
        p1->interpolation = (ptwXY_interpolation) i1;
        p2 = ptwXY_toOtherInterpolation( p1, ptwXY_interpolationLinLin, accuracy, &status );
        if( status != statuses[i1] ) ++errCounts;
    }
    ptwXY_free( p1 );
    exit( errCounts ? EXIT_FAILURE : EXIT_SUCCESS );
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
