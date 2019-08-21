/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char const *value2MathValue( double value, char str[128] );
/*
************************************************************
*/
char const *value2MathValue( double value, char str[128] );
int main( int argc, char **argv ) {

    int n = 1;
    double x;
    char str[128];

    for( n = 0; n < 20; n++ ) {
        for( x = 1e-10; x < 500; x *= 1.2 ) {
            value2MathValue( x, str );
            printf( "Print[ \"{ \", %3d, \", %+25.17e, \", CForm[ N[ ExpIntegralE[ %d, %s ], 20 ] ], \" },\" ]\n", n, x, n, str );
        }
    }

    exit( EXIT_SUCCESS );
}
/*
************************************************************
*/
static char const *value2MathValue( double value, char str[128] ) {

    int i;

    value *= 1e27;
    sprintf( str, "%.0f / 1", value );
    for( i = 0; i < 27; i++ ) strcat( str, "0" );
    return( str );
}
