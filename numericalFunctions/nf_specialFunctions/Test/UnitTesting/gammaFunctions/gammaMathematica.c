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

    double a;
    char str1[128];

    for( a = 1e-12; a < 20; a *= 1.1 ) {
        value2MathValue( a, str1 );
        printf( "Print[ \"{ %+25.17e, \", CForm[ N[ Gamma[ %s, 0 ], 20 ] ], \" },\" ]\n", a, str1 );
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
