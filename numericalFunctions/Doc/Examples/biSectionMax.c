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

void doProductBiSection( ptwXYPoints *ptwXY3, ptwXYPoints *ptwXY4, int biSection );
/*
********************************************************
*/
int main( int argc, char **argv ) {
    
    double f3[4] = { 0., 0., 1., 1. }, f4[4] = { 0., 1., 1., 0. };
    ptwXYPoints *ptwXY3, *ptwXY4;
    ptwXY_interpolation linlin = ptwXY_interpolationLinLin;
    
    ptwXY3 = ptwXY_create2( NULL, linlin, 10, 10, 2, f3, 0 );
    ptwXY4 = ptwXY_create2( NULL, linlin, 10, 10, 2, f4, 0 );
    
    doProductBiSection( ptwXY3, ptwXY4, 0 );
    doProductBiSection( ptwXY3, ptwXY4, 1 );
    doProductBiSection( ptwXY3, ptwXY4, 2 );
    doProductBiSection( ptwXY3, ptwXY4, 3 );

    ptwXY_free( ptwXY3 );
    ptwXY_free( ptwXY4 );
}

/*
********************************************************
*/
void doProductBiSection( ptwXYPoints *ptwXY3, ptwXYPoints *ptwXY4, int biSection ) {

    ptwXYPoints *product;

    printf( "\n\n# biSection = %d\n", biSection );

    ptwXY_setBiSectionMax( ptwXY3, biSection );
    ptwXY_setBiSectionMax( ptwXY4, biSection );
    product = ptwXY_mul2_ptwXY( NULL, ptwXY3, ptwXY4 );

    ptwXY_simplePrint( product, " %.12f | %.12f\n" );

    ptwXY_free( product );
}
