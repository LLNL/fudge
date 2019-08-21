/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void printData( double E, double cs );
/*
******************************************************
*/
int main( ) {

	int i, n = 41;
	double x, w = 1e-3, E, cs, fE, fw = 3., y, offset = 0.;

	fE = pow( ( 1. + fw * w ) / ( 1. - fw * w ), 1. / ( n - 1 ) );
	printData( 1e-10, offset );
	for( x = 1e-8; x < 11; x *= 10. ) {
		E = x * ( 1. - fw * w );
		for( i = 0; i < n ; i++ ) {
			y = ( E / x - 1. ) / w;
			cs = offset + 100. * exp( - y * y );
			printData( E, cs );
			E *= fE;
		}
	}
	printData( 20., offset );
	exit( EXIT_SUCCESS );
}
/*
******************************************************
*/
void printData( double E, double cs ) {

	printf( "%18.10e %15.6e\n", E, cs );
}
