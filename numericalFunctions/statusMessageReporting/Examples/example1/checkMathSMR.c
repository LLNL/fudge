/*
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include "checkMathSMR.h"

static int checkMathSMRID = smr_unknownID;
/*
=================================================
*/
int checkMathSMR_setup( void ) {

    checkMathSMRID = smr_registerLibrary( "checkMathSMRLibrary" );
    return( checkMathSMRID );
}
/*
=================================================
*/
int checkMathSMR_getLibrarysID( void ) {

    return( checkMathSMRID );
}
/*
=================================================
*/
int checkMathSMR_test( statusMessageReporting *smr, int verbose ) {

    double x1, xMin = 1e-5, xMax = 1e5, value;

    for( x1 = xMin; x1 < xMax; x1 *= 10. ) {
        value = mathSMR_sin( smr, x1 );
        if( smr_isOk( smr ) == 0 ) {
            if( verbose )
                smr_setReportError2p( smr, checkMathSMRID, checkMathSMR_code_sin, "Checking math failed for mathSMR_sin." );
            return( 1 );
        }
    }

    for( x1 = xMin; x1 < xMax; x1 *= 10. ) {
        value = mathSMR_exp( smr, x1 );
        if( smr_isOk( smr ) == 0 ) {
            if( verbose )
                smr_setReportError2p( smr, checkMathSMRID, checkMathSMR_code_exp, "Checking math failed for mathSMR_exp." );
            return( 1 );
        }
    }

    return( 0 );
}
