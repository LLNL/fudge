/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <math.h>
#include "mathSMR.h"

static int mathSMRID = smr_unknownID;
/*
=================================================
*/
int mathSMR_setup( void ) {

	mathSMRID = smr_registerLibrary( "mathSMRLibrary" );
	return( mathSMRID );
}
/*
=================================================
*/
int mathSMR_getLibrarysID( void ) {

	return( mathSMRID );
}
/*
=================================================
*/
double mathSMR_sin( statusMessageReporting *smr, double x1 ) {

	double returnValue = sin( x1 );

	if( isfinite( returnValue ) == 0 ) {
    	smr_setReportError2( smr, mathSMRID, mathSMR_code_sin, "returnValue of <%25.17e> for input of <%25.17e>.", returnValue, x1 );
	}
	return( returnValue );
}
/*
=================================================
*/
double mathSMR_exp( statusMessageReporting *smr, double x1 ) {

	double returnValue = exp( x1 );

	if( isfinite( returnValue ) == 0 ) {
    	smr_setReportError2( smr, mathSMRID, mathSMR_code_exp, "returnValue of <%25.17e> for input of <%25.17e>.", returnValue, x1 );
	}
	return( returnValue );
}
