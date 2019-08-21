/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/
/*
    This routine tests the nf_Legendre_GaussianQuadrature integrator by comparing the integral of x^n (n = 0 to MAX_DEGREE)
between x1 and x2 to the soluation ( x2^(n+1) - x1^(n+1) ) / ( n + 1 ).
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <ptwXY.h>
#include <nf_Legendre.h>

#define MAX_DEGREE 80
#define TOL 2e-13

static int verbose = 0;

static int GaussianQuadrature( int degree, double x1, double x2 );
static int GaussianQuadrature2( int degree, double x1, double x2 );
static nfu_status GaussianQuadrature_callback( double x, double *f, void *argList );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, errCount;

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


    for( i = 0; i <= MAX_DEGREE; i++ ) errCount += GaussianQuadrature( i, 0, 1 );
    for( i = 0; i <= MAX_DEGREE; i++ ) errCount += GaussianQuadrature( i, 4.5, 6.2 );
    for( i = 0; i <= MAX_DEGREE; i++ ) errCount += GaussianQuadrature( i, -1.2, .7 );

    if( errCount ) fprintf( stderr, "%s FAILED\n", __FILE__ );
    return( errCount );
}
/*
************************************************************
*/
static int GaussianQuadrature( int degree, double x1, double x2 ) {

    return( GaussianQuadrature2( degree, x1, x2 ) + GaussianQuadrature2( degree, x2, x1 ) );
}
/*
************************************************************
*/
static int GaussianQuadrature2( int degree, double x1, double x2 ) {

    int errCount = 0;
    double integral, exact, r;
    nfu_status status;

    status = nf_Legendre_GaussianQuadrature( degree, x1, x2, GaussianQuadrature_callback, (void *) &degree, &integral );
    if( status == nfu_Okay ) {
        exact = ( pow( x2, degree + 1 ) - pow( x1, degree + 1 ) ) / ( degree + 1 );
        r = integral / exact - 1.0;
        if( fabs( r ) > TOL ) errCount++;
        if( verbose ) printf( "%3d %25.16e %25.16e %20.12e %20.12e %12.5e\n", degree, x1, x2, exact, integral, r ); }
    else {
        errCount++;
    }
    return( errCount );
}
/*
************************************************************
*/
static nfu_status GaussianQuadrature_callback( double x, double *f, void *argList ) {

    int degree = *((int *) argList);

    *f = pow( x, degree );
    return( nfu_Okay );
}
