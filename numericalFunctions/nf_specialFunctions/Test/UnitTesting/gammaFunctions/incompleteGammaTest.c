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

#include <nfut_utilities.h>
#include <nf_specialFunctions.h>

struct nXFOfX {
    double a, x, f;
};

#include "incompleteGammaTest.dat"

#define nPowerErrors 10
static int verbose = 0;

/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i1, i2, nData = sizeof( data ) / sizeof( data[0] ), counts = 0, iarg, echo = 0, powerErrors[nPowerErrors+1], info = 0;
    double f, r, r2;
    nfu_status status;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

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
    if( echo ) printf( "nf_specialFunctions: %s\n", __FILE__ );

    for( i2 = 0; i2 <= nPowerErrors; i2++ ) powerErrors[i2] = 0;

    for( i1 = 0; i1 < nData; i1++ ) {
        status = nf_incompleteGammaFunction( &smr, data[i1].a, data[i1].x, &f );
        r = 1;
        if( data[i1].f != 0. ) {
            r = f / data[i1].f - 1; }
        else {
            if( f == 0. ) r = 0.;
        }
        if( ( fabs( r ) > 1e-13 ) || status ) {
            printf( "%.17e %.17e %.17e %.17e %+.3e  %d\n", data[i1].a, data[i1].x, data[i1].f, f, r, status );
            nfut_printSMRError2p( &smr, "Via." );
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
