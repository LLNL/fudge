/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>

#include <ptwXY.h>

void doProductAccuracy( ptwXYPoints *ptwXY3, ptwXYPoints *ptwXY4, int biSectionMax, double accuracy );
/*
********************************************************
*/
int main( int argc, char **argv ) {
    
    double f3[4] = { 0.0, 1.0, 1.0, 2.0 }, f4[4] = { 0.0, 1.0, 1.0, 2.0 };
    ptwXYPoints *ptwXY3, *ptwXY4;
    ptwXY_interpolation linlin = ptwXY_interpolationLinLin;
    
    ptwXY3 = ptwXY_create2( NULL, linlin, 10, 10, 2, f3, 0 );
    ptwXY4 = ptwXY_create2( NULL, linlin, 10, 10, 2, f4, 0 );
    
    doProductAccuracy( ptwXY3, ptwXY4, 12, 1e-3 );
    doProductAccuracy( ptwXY3, ptwXY4, 12, 1e-2 );
    doProductAccuracy( ptwXY3, ptwXY4, 12, 1e-1 );
    doProductAccuracy( ptwXY3, ptwXY4, 2, 1e-3 );

    ptwXY_free( ptwXY3 );
    ptwXY_free( ptwXY4 );
}
/*
********************************************************
*/
void doProductAccuracy( ptwXYPoints *ptwXY3, ptwXYPoints *ptwXY4, int biSectionMax, double accuracy ) {

    ptwXYPoints *product;

    printf( "\n\n# accuracy = %.2g\n", accuracy );
    printf( "# biSection = %d\n", biSectionMax );

    ptwXY_setBiSectionMax( ptwXY3, biSectionMax );
    ptwXY_setBiSectionMax( ptwXY4, biSectionMax );
    ptwXY_setAccuracy( ptwXY3, accuracy );
    ptwXY_setAccuracy( ptwXY4, accuracy );

    product = ptwXY_mul2_ptwXY( NULL, ptwXY3, ptwXY4 );

    ptwXY_simplePrint( product, " %.12f  %.12f\n" );

    ptwXY_free( product );
}
