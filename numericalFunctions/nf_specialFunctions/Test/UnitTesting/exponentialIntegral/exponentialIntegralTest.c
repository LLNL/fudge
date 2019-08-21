/*
# <<BEGIN-copyright>>   
# <<END-copyright>> 
*/

#include <stdio.h>
#include <stdlib.h>
#include <nf_specialFunctions.h>

struct nXEnOfX {
    int n;
    double x, En;
};

#include "exponentialIntegralTest.dat"

#define nPowerErrors 10
static int verbose = 0;

/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i1, i2, nData = sizeof( data ) / sizeof( data[0] ), counts = 0, iarg, echo = 0, powerErrors[nPowerErrors+1], info = 0;
    double En, r, r2;
    nfu_status status;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else if( strcmp( "-i", argv[iarg] ) == 0 ) {
            info = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    for( i2 = 0; i2 <= nPowerErrors; i2++ ) powerErrors[i2] = 0;

    for( i1 = 0; i1 < nData; i1++ ) {
        En = nf_exponentialIntegral( data[i1].n, data[i1].x, &status );
        r = 1;
        if( data[i1].En != 0. ) r = En / data[i1].En - 1;
        if( ( r > 3e-14 ) || status ) {
            printf( "%3d %.17e %.17e %.17e %+.3e  %d\n", data[i1].n, data[i1].x, data[i1].En, En, r, status );
            counts++;
        }
        for( i2 = 0, r2 = 1e-16; i2 < nPowerErrors; i2++, r2 *= 10. ) {
            if( fabs( r ) < r2 ) break;
        }
        powerErrors[i2]++;
    }
    if( info ) {
        printf( "relative" );
        for( i2 = 0; i2 <= nPowerErrors; i2++ ) printf( " %7d", -i2 - 6 );
        printf( "\n" );
        printf( "error:  " );
        for( i2 = 0; i2 <= nPowerErrors; i2++ ) printf( "%s%5d  %s", ( ( i2 < 4 ) ? " " : "" ), 10, ( ( i2 >= 4 ) ? " " : "" ) );
        printf( "\n" );

        printf( "--------" );
        for( i2 = 0; i2 <= nPowerErrors; i2++ ) printf( "--------" );
        printf( "\n" );

        printf( "counts: " );
        for( i2 = nPowerErrors; i2 >= 0; i2-- ) printf( " %7d", powerErrors[i2] );
        printf( "\n" );
    }
    exit( counts );
}
