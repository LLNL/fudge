/*
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdlib.h>
#include <statusMessageReporting.h>
#include <checkMathSMR.h>

int main( int argc, char **argv ) {

	statusMessageReporting smr;
   	int libID = smr_registerLibrary( "MyMainFunction1111" );
    int verbose = argc == 1;

   	smr_initialize( &smr, smr_status_Ok );
    checkMathSMR_setup( );
    mathSMR_setup( );

	if( checkMathSMR_test( &smr, verbose ) && verbose ) {
   		smr_setReportError2p( &smr, libID, 1, "Let's add one Error report at the top level for fun." );
	}
  
   	printf( ">> number of reports = %d\n", smr_numberOfReports( &smr ) );
   	smr_print( &smr, 1 );

   	smr_cleanup( );
	exit( EXIT_SUCCESS );
}
