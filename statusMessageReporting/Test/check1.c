/*
# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
# <<END-copyright>>
*/

#include "statusMessageReporting.h"
#include <stdlib.h>

static int verbose = 0;
static int check1ID = smr_unknownID;

static char const text1[] = "Hi";

void addMessage( statusMessageReporting *smr1, statusMessageReporting *smr2, enum smr_status status, int count, char const *message );
void printMessages( statusMessageReporting *smr1, statusMessageReporting *smr2, int clear );
void releaseMessages( statusMessageReporting *smr1, statusMessageReporting *smr2 );
enum smr_status getStatus( void );
void timeToAbort( int i );
/*
============================================================
*/
int main( int argc, char **argv ) {

    int i, counter;
    enum smr_status status;
    statusMessageReporting static_smr, *smr1 = &static_smr, *smr2;
    char msg[128];

    if( argc > 1 ) verbose = 1;
    smr_setup( );
    check1ID = smr_registerLibrary( "check1" );

    smr_initialize( smr1, smr_status_Warning );
    if( ( smr2 = smr_new( smr1, smr_status_Ok ) ) == NULL ) {
        smr_print( smr1, 1 );
        exit( EXIT_FAILURE );
    }

    addMessage( smr1, smr2, smr_status_Info, 1, text1 );
    printMessages( smr1, smr2, 0 );
    releaseMessages( smr1, smr2 );

    for( i = 0, counter = 0; i < 400000; i++ ) {
        counter++;
        status = getStatus( );
        if( verbose ) printf( "\ni = %d, status = %d, <%s>\n", i, status, smr_statusToString( status ) );
        sprintf( msg, "This is message %d with status = %d, string = <%s>\n", i, status, smr_statusToString( status ) );
        addMessage( smr1, smr2, status, 1, msg );
        if( drand48( ) < .1 ) {
            if( smr_numberOfReports( smr2 ) != counter ) fprintf( stderr, "smr_numberOfReports( smr2 ) = %d != counter = %d at index = %d\n", 
                smr_numberOfReports( smr2 ), counter, i );
            printMessages( smr1, smr2, 1 );
            counter = 0;
        }
    }
    printMessages( smr1, smr2, 1 );
    releaseMessages( smr1, smr2 );
    smr_free( &smr2 );
    smr_cleanup( );

    exit( EXIT_SUCCESS );
}
/*
============================================================
*/
void addMessage( statusMessageReporting *smr1, statusMessageReporting *smr2, enum smr_status status, int count, char const *message ) {

    int i;

    switch( status ) {
    case smr_status_Ok :
        smr_release( smr1 );
        smr_release( smr2 );
        break;
    case smr_status_Info :
        if( ( i = smr_setReportInfo( smr1, NULL, __FILE__, __LINE__, __func__, check1ID, 30, "count = %3d, message = %s", count, message ) ) ) timeToAbort( i );
        if( ( i = smr_setReportInfo( smr2, NULL, __FILE__, __LINE__, __func__, check1ID, 30, "count = %3d, message = %s", count, message ) ) ) timeToAbort( i );
        break;
    case smr_status_Warning :
        if( ( i = smr_setReportWarning( smr1, NULL, __FILE__, __LINE__, __func__, check1ID, 110, "count = %3d, message = %s", count, message ) ) ) timeToAbort( i );
        if( ( i = smr_setReportWarning( smr2, NULL, __FILE__, __LINE__, __func__, check1ID, 110, "count = %3d, message = %s", count, message ) ) ) timeToAbort( i );
        break;
    case smr_status_Error :
        if( ( i = smr_setReportError( smr1, NULL, __FILE__, __LINE__, __func__, check1ID, 510, "count = %3d, message = %s", count, message ) ) ) timeToAbort( i );
        if( ( i = smr_setReportError( smr2, NULL, __FILE__, __LINE__, __func__, check1ID, 510, "count = %3d, message = %s", count, message ) ) ) timeToAbort( i );
        break;
    default :
        fprintf( stderr, "Invalid status = %d\n", status );
        timeToAbort( 2 );
    }
}
/*
============================================================
*/
void printMessages( statusMessageReporting *smr1, statusMessageReporting *smr2, int clear ) {

    if( verbose ) {
        printf( "Printing smr1 (clear = %d) that has %d messages\n", clear, smr_numberOfReports( smr1 ) );
        smr_print( smr1, clear );
        printf( "Printing smr2 (clear = %d) that has %d messages\n", clear, smr_numberOfReports( smr2 ) );
        smr_print( smr2, clear ); }
    else {
        if( clear ) releaseMessages( smr1, smr2 );
    }
}
/*
============================================================
*/
void releaseMessages( statusMessageReporting *smr1, statusMessageReporting *smr2 ) {

    if( verbose ) printf( "Releasing smr1\n" );
    smr_release( smr1 );
    if( verbose ) printf( "Releasing smr2\n" );
    smr_release( smr2 );
}
/*
============================================================
*/
enum smr_status getStatus( void ) {

    int i = (int) 3 * drand48( );

    switch( i ) {
    case 0 : return( smr_status_Info );
    case 1 : return( smr_status_Warning );
    case 2 : return( smr_status_Error );
    }
    fprintf( stderr, "getStatus i value = %d is invalid\n", i );
    timeToAbort( 0 );
    return( smr_status_Ok );        /* Will never get here but need to make the compilers happy. */
}
/*
============================================================
*/
void timeToAbort( int i ) {

    fprintf( stderr, "i = %d\n", i );
    exit( EXIT_FAILURE );
}
