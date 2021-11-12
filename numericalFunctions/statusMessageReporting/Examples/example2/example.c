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

int main( int argc, char **argv ) {

	statusMessageReporting smr;
   	int libID = smr_registerLibrary( "anotherLibrary2222" );

   	smr_initialize( &smr, smr_status_Ok );
   	smr_setReportInfo2p( &smr, libID, 1, "INFO report with code 1" );
   	smr_setReportWarning2p( &smr, libID, 22, "WARNING report with code 22" );
   	smr_setReportError2p( &smr, libID, 3, "ERROR report with code 3" );
   	smr_setReportError2p( &smr, libID, 33, "second ERROR report with code 33" );
  
   	printf( ">> number of reports = %d\n", smr_numberOfReports( &smr ) );
   	smr_print( &smr, 1 );

   	smr_cleanup( );
	exit( EXIT_SUCCESS );
}
