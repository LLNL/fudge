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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include <ptwXY.h>
#include <nfut_utilities.h>
#include <ptwXY_utilities.h>

static int verbose = 0;
static char *fmtXY = "%25.15e %25.15e\n";
static FILE *infoF;

void printMsg( const char *fmt, ... );
static void printIfVerbose( ptwXYPoints *data );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, echo = 0, nXYs;
    ptwXYPoints *pFlat, *pLinear;
    double accuracy = 1e-3, lowerEps = 1e-6, upperEps = 1e-6, xys[] = { 
        0.00000000e+00, 2.44941900e-08, 3.00000000e+05, 2.95502600e-08, 5.00000000e+05, 1.60509700e-07, 7.00000000e+05, 7.29768400e-07,
        9.00000000e+05, 7.16263300e-08, 1.10000000e+06, 3.13305400e-07, 1.30000000e+06, 9.24667600e-08, 1.50000000e+06, 1.33916500e-07,
        1.70000000e+06, 2.15220800e-07, 1.90000000e+06, 4.09025600e-07, 2.10000000e+06, 7.94523400e-08, 2.30000000e+06, 1.26046700e-07,
        2.50000000e+06, 9.74776700e-08, 2.70000000e+06, 1.82226500e-07, 2.90000000e+06, 1.40790400e-07, 3.10000000e+06, 1.26457400e-07,
        3.30000000e+06, 1.11357400e-07, 3.50000000e+06, 1.41709800e-07, 3.70000000e+06, 1.34958900e-07, 3.90000000e+06, 1.28301200e-07,
        4.10000000e+06, 1.21774200e-07, 4.30000000e+06, 1.15381700e-07, 4.50000000e+06, 1.09135600e-07, 4.70000000e+06, 1.03051700e-07,
        4.90000000e+06, 9.71452600e-08, 5.10000000e+06, 9.14627000e-08, 5.30000000e+06, 8.59984900e-08, 5.50000000e+06, 8.06598500e-08,
        5.70000000e+06, 7.55374200e-08, 5.90000000e+06, 7.06618500e-08, 6.10000000e+06, 2.96670100e-08, 6.30000000e+06, 5.85612600e-08,
        6.50000000e+06, 7.42873500e-08, 6.70000000e+06, 5.47748800e-08, 6.90000000e+06, 1.03412200e-07, 7.10000000e+06, 5.89404900e-09,
        7.30000000e+06, 4.59737500e-09, 7.50000000e+06, 7.45342100e-08, 7.70000000e+06, 2.55224600e-09, 7.90000000e+06, 2.90521500e-08,
        8.10000000e+06, 1.19440100e-09, 8.30000000e+06, 3.76551100e-08, 8.50000000e+06, 3.44437000e-10, 8.70000000e+06, 1.11754500e-07,
        8.90000000e+06, 5.53869000e-13, 9.10000000e+06, 8.23249000e-14, 9.30000000e+06, 7.15863000e-15, 9.49999990e+06, 1.47460800e-07,
        9.50000000e+06, 0.00000000e+00 };
    ptwXY_interpolation interpolation = ptwXY_interpolationFlat;
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    infoF = stdout;

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else {
            printMsg( "Error %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) fprintf( stderr, "%s\n", __FILE__ );

    nfu_setMemoryDebugMode( 0 );
    nXYs = sizeof( xys ) / ( 2 * sizeof( xys[0] ) );
    if( ( pFlat = ptwXY_create( &smr, interpolation, NULL, 5, accuracy, 10, 10, nXYs, xys, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    printIfVerbose( pFlat );

    if( ( pLinear = ptwXY_flatInterpolationToLinear( &smr, pFlat, lowerEps, upperEps ) ) == NULL ) 
        nfut_printSMRErrorExit2p( &smr, "Via." );
    printIfVerbose( pLinear );

    ptwXY_free( pFlat );
    ptwXY_free( pLinear );

    exit( EXIT_SUCCESS );
}
/*
****************************************************************
*/
void printMsg( const char *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    fprintf( stderr, "\n" );
    va_end( args );
    exit( EXIT_FAILURE );
}
/*
************************************************************
*/
static void printIfVerbose( ptwXYPoints *data ) {

    if( !verbose ) return;
    fprintf( infoF, "# length = %d\n", (int) ptwXY_length( NULL, data ) );
    ptwXY_simpleWrite( data, infoF, fmtXY );
    fprintf( infoF, "\n\n" );
}
