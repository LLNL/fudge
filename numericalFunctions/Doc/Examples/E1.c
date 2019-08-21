/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>

#include <ptwXY.h>
#define nPairs 7

/*
********************************************************
*/
int main( int argv, char **agrv ) {

    double c;
    double maleData[2 * ( nPairs - 1 )] =   { 1871, 1212, 1883, 1215, 1889,  51, 1895, 11, 1905,  9, 1915,  9 };
    double femaleData[2 * nPairs] = { 1871, 1231, 1883, 1241, 1885, 621, 1889, 229, 1895, 31, 1905, 23, 1915, 21 };
    ptwXYPoints *males, *females, *total;
    nfu_status status;

    males = ptwXY_new( ptwXY_interpolationLinLin, 5, 1e-3, nPairs, 4, &status, 0 );
    ptwXY_setXYData( males, nPairs - 1, maleData );

    females = ptwXY_new( ptwXY_interpolationLinLin, 5, 1e-3, nPairs, 4, &status, 0 );
    ptwXY_setXYData( females, nPairs, femaleData );

    total = ptwXY_add_ptwXY( males, females, &status );

    printf( "\nMale population\n" );
    printf( "  Year | Count\n" );
    printf( " -------+-----\n" );
    ptwXY_simplePrint( males, " %5.0f | %5.0f\n" );

    printf( "\nFemale population\n" );
    printf( "  Year | Count\n" );
    printf( " ------+------\n" );
    ptwXY_simplePrint( females, " %5.0f | %5.0f\n" );

    printf( "\nTotal population\n" );
    printf( "  Year | Count\n" );
    printf( " ------+------\n" );
    ptwXY_simplePrint( total, " %5.0f | %5.0f\n" );

    ptwXY_getValueAtX( males, 1885, &c );
    printf( "\nInterpolated male population at 1885 is %.0f\n", c );

    exit( EXIT_SUCCESS );
}
