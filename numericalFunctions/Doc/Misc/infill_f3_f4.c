/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ptwXY.h>
#include <nf_utilities.h>

static char *fmtXY = "%18.11e %18.11e\n";

void doProduct( ptwXYPoints *ptwXY3, ptwXYPoints *ptwXY4, double biSection );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    double f3[4] = { 0., 0., 1., 1. }, f4[4] = { 0., 1., 1., 0. }, accuracy = 1e-3;
    ptwXYPoints *ptwXY3, *ptwXY4;
    nfu_status status;

    if( ( ptwXY3 = ptwXY_create( ptwXY_interpolationLinLin, 0, accuracy, 10, 10, 2, f3, &status, 0 ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: ptwXY3 creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );
    if( ( ptwXY4 = ptwXY_create( ptwXY_interpolationLinLin, 0, accuracy, 10, 10, 2, f4, &status, 0 ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: ptwXY4 creation status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    printf( "# f3\n" );
    ptwXY_simpleWrite( ptwXY3, stdout, fmtXY );
    printf( "\n\n# f4\n" );
    ptwXY_simpleWrite( ptwXY4, stdout, fmtXY );

    doProduct( ptwXY3, ptwXY4, 0 );
    doProduct( ptwXY3, ptwXY4, 1 );
    doProduct( ptwXY3, ptwXY4, 2 );
    doProduct( ptwXY3, ptwXY4, 3 );
    doProduct( ptwXY3, ptwXY4, 4 );
    doProduct( ptwXY3, ptwXY4, 10 );

    ptwXY_free( ptwXY3 );
    ptwXY_free( ptwXY4 );
    exit( EXIT_SUCCESS );

}
/*
************************************************************
*/
void doProduct( ptwXYPoints *ptwXY3, ptwXYPoints *ptwXY4, double biSection ) {

    ptwXYPoints *product;
    nfu_status status;

    ptwXY_setBiSectionMax( ptwXY3, biSection );
    if( ( product = ptwXY_mul2_ptwXY( ptwXY3, ptwXY4, &status ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: multiplication status = %d: %s", __FILE__, status, nfu_statusMessage( status ) );

    printf( "\n\n# biSection = %.1f\n", biSection );
    printf( "# length = %d\n", (int) ptwXY_length( product ) );
    ptwXY_simpleWrite( product, stdout, fmtXY );
    ptwXY_free( product );
}
