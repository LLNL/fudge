/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
/*
    This routine creates a nf_Legendre using nf_Legendre_new, converts it to a ptwXYPoints using nf_Legendre_to_ptwXY
and then converts it back to an nf_Legendre using nf_Legendre_from_ptwXY. The Legendre coefficients between the
orginal and re-converts nf_Legendre are compared.
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <ptwXY.h>
#include <nf_Legendre.h>

#define MAX_ORDER 16

static int verbose = 0, biSectionMax = 12;
static double accuracy = 1e-5;

static int to_ptwXY( int maxOrder, double *Cls );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount = 0;
    double Cls[MAX_ORDER+1];

    Cls[0] = 1.;
    for( i = 0; i < MAX_ORDER; i++ ) Cls[i+1] = 0.5 * Cls[i];

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "nf_Legendre: %s\n", __FILE__ );


    for( i = 0; i <= MAX_ORDER; i++ ) errCount += to_ptwXY( i, Cls );

    if( errCount ) fprintf( stderr, "%s FAILED\n", __FILE__ );
    return( errCount );
}
/*
************************************************************
*/
static int to_ptwXY( int maxOrder, double *Cls ) {

    int i, errCount, errCounts = 1;
    double d, r;
    nf_Legendre *nfL1, *nfL2;
    ptwXYPoints *XYs;
    nfu_status status;

    if( ( nfL1 = nf_Legendre_new( 10, maxOrder, Cls, &status ) ) == NULL ) {
        nfu_printErrorMsg( "nf_Legendre_new failed with status = %d: '%s'", status, nfu_statusMessage( status ) ); }
    else {
        if( ( XYs = nf_Legendre_to_ptwXY( nfL1, accuracy, biSectionMax, 0, &status ) ) == NULL ) {
            nfu_printErrorMsg( "nf_Legendre_to_ptwXY failed with status = %d: '%s'", status, nfu_statusMessage( status ) ); }
        else {
            if( ( nfL2 = nf_Legendre_from_ptwXY( XYs, maxOrder, &status ) ) != NULL ) {
                if( verbose ) {
                    printf( "%d %d\n", nf_Legendre_maxOrder( nfL1 ), nf_Legendre_maxOrder( nfL2 ) );
                    for( i = 0; i <= nf_Legendre_maxOrder( nfL1 ); i++ ) printf( "%3d %.12e\n", i, nf_Legendre_getCl( nfL1, i, &status ) );
                    printf( "\n" );
                }
                errCounts = 0;
                for( i = 0; i <= nf_Legendre_maxOrder( nfL2 ); i++ ) {
                    r = nf_Legendre_getCl( nfL1, i, &status );
                    d = r - nf_Legendre_getCl( nfL2, i, &status );
                    if( r != 0 ) {
                        r = d / r; }
                    else {
                        if( d != 0 ) r = 1;
                    }
                    errCount = 0;
                    if( i < 7 ) {               /* These values are all empirical. */
                        if( fabs( r ) > accuracy ) errCount = 1; }
                    else if( i < 10 ) {
                        if( fabs( r ) > 10 * accuracy ) errCount = 1; }
                    else {
                        if( fabs( r ) > 200 * accuracy ) errCount = 1;
                    }
                    if( verbose ) printf( "%3d %.12e  %12.4e %12.4e%s\n", i, nf_Legendre_getCl( nfL2, i, &status ), d, r, ( errCount ? " ====" : "" ) );
                    errCounts += errCount;
                }
                if( verbose ) printf( "\n" );
                nf_Legendre_free( nfL2 );
            }
            ptwXY_free( XYs );
        }
        nf_Legendre_free( nfL1 );
    }

    return( errCounts );
}
