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

#include <nf_utilities.h>
#include <nfut_utilities.h>
#include <ptwX.h>

#define nXs 23

static int verbose = 0;

int check( statusMessageReporting *smr, ptwXPoints *ptwX1, ptwXPoints *ptwX2, char const *msg );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, errors = 0;
    ptwXPoints *ptwX1, *ptwX2;
    statusMessageReporting smr;

    nfu_setup( );
    smr_initialize( &smr, smr_status_Ok );

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfut_printSMRError2( &smr, "Error: invalid input option '%s'", argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    if( ( ptwX1 = ptwX_createLine( &smr, 10, nXs, -2, 100 ) ) == NULL ) nfut_printSMRError2p( &smr, "Via." );

    if( ( ptwX2 = ptwX_createLine( &smr, 10, nXs, 1, 0 ) ) == NULL ) nfut_printSMRError2p( &smr, "Via." );
    if( ptwX_slopeOffset( &smr, ptwX2, -2, 100 ) != nfu_Okay ) nfut_printSMRError2p( &smr, "via" );
    errors += check( &smr, ptwX1, ptwX2, "ptwX_slopeOffset( ptwX2, -2, 100 )" );

    if( ( ptwX2 = ptwX_createLine( &smr, 10, nXs, 1, 0 ) ) == NULL ) nfut_printSMRError2p( &smr, "via" );
    if( ptwX_mul_double( &smr, ptwX2,  -2 ) != nfu_Okay ) nfut_printSMRError2p( &smr, "via" );
    if( ptwX_add_double( &smr, ptwX2, 100 ) != nfu_Okay ) nfut_printSMRError2p( &smr, "via" );
    errors += check( &smr, ptwX1, ptwX2, "ptwX_mul_double then ptwX_add_double" );

    ptwX_free( ptwX1 );
    exit( errors );
}
/*
****************************************************************
*/
int check( statusMessageReporting *smr, ptwXPoints *ptwX1, ptwXPoints *ptwX2, char const *msg ) {

    int close;
    nfu_status status = ptwX_close( smr, ptwX1, ptwX2, 3, 0., &close );

    if( status != nfu_Okay ) nfut_printSMRError2p( smr, "via" );
    close -= (int) ptwX_length( smr, ptwX1 );
    ptwX_free( ptwX2 );
    return( close );
}
