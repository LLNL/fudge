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

/*
    This routine tests the insertion of points where a point already exists.
    The number of initial points, currently 7, must be a little smaller than the space allocated for 
    points (currently 10) and overflowPoints (currently 10).
*/
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>

#include <nfut_utilities.h>
#include <ptwXY.h>

static int verbose = 0;

void printMsg( char const *fmt, ... );
/*
****************************************************************
*/
int main( int argc, char **argv ) {

    ptwXYPoint *point;
    ptwXYPoints *u;
    double uXY[] = { 3.732, 0, 5.598, 0, 5.679, 0, 6.459, 0, 6.557, 0, 6.838, 0, 7.1138, 0 };
    int i, j, n, nuXY = sizeof( uXY ) / ( 2 * sizeof( double ) ), iarg, echo = 0;
    FILE *ff;
    char *fmt = "%7.4f, %5.1f\n";
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            printMsg( "Error %s: unsupported option = '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    nfu_setMemoryDebugMode( 0 );

    if( ( u = ptwXY_create( &smr, ptwXY_interpolationLinLin, NULL, 5, 1e-3, 10, 10, nuXY, uXY, 0 ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );

    if( ( ff = fopen( "curve_u.dat", "w" ) ) == NULL ) printMsg( "Error %s: could not open output file", __FILE__ );
    ptwXY_simpleWrite( u, ff, fmt );
    fprintf( ff, "\n" );
    if( verbose ) ptwXY_simpleWrite( u, stdout, fmt );
    if( verbose ) printf( "\n" );

/*
 Overwrite existing y values for each point. This should not change the number of points.
*/
    for( i = 0; i < nuXY; i++ ) {
        if( ptwXY_setValueAtX( &smr, u, uXY[2 * i], -i - 1. ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    if( ptwXY_length( &smr, u ) != nuXY )
        printMsg( "Error %s: ptwXY_length( &smr, u ) = %d != nuXY = %d", __FILE__, (int) ptwXY_length( &smr, u ), nuXY );
    ptwXY_simpleWrite( u, ff, fmt );
    fprintf( ff, "\n" );
    if( verbose ) ptwXY_simpleWrite( u, stdout, fmt );
    if( verbose ) printf( "\n" );

/*
 Add a point before xMin, between each point and after xMax. This should change the number of points to 2 * nuXY + 1.
*/
    if( ptwXY_setValueAtX( &smr, u, uXY[0] - .1, 4 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    for( i = 0; i < nuXY - 1; i++ ) {
        if( ptwXY_setValueAtX( &smr, u, 0.5 * ( uXY[2 * i] + uXY[2 * i + 2] ), i + 5. ) != nfu_Okay )
            nfut_printSMRErrorExit2p( &smr, "Via." );
    }
    if( ptwXY_setValueAtX( &smr, u, uXY[2 * nuXY - 2] + .1, 11 ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
    if( ptwXY_length( &smr, u ) != 2 * nuXY + 1 )
        printMsg( "Error %s: ptwXY_length( &smr, u ) = %d != nuXY = %d", __FILE__, (int) ptwXY_length( &smr, u ), nuXY );
    ptwXY_simpleWrite( u, ff, fmt );
    fprintf( ff, "\n" );
    if( verbose ) ptwXY_simpleWrite( u, stdout, fmt );
    if( verbose ) printf( "\n" );
    if( verbose ) ptwXY_showInteralStructure( u, stdout, 1 );

/*
 Adds some points and after each add overwrite all existing y values for each point. The overwrite should not change the number of points.
*/
    for( i = 0; i < nuXY - 1; i++ ) {
        if( ptwXY_setValueAtX( &smr, u, 0.75 * uXY[2 * i] + 0.25 * uXY[2 * i + 2], i + 55. ) != nfu_Okay )
            nfut_printSMRErrorExit2p( &smr, "Via." );
        n = ptwXY_length( &smr, u );
        for( j = 0; j < n; j++ ) {
            point = ptwXY_getPointAtIndex_Unsafely( u, j );
            if( ptwXY_setValueAtX( &smr, u, point->x, i + j ) != nfu_Okay ) nfut_printSMRErrorExit2p( &smr, "Via." );
            if( ptwXY_length( &smr, u ) != n )
                printMsg( "Error %s: ptwXY_length( &smr, u ) = %d != n = %d", __FILE__, (int) ptwXY_length( &smr, u ), n );
        }
    }
    ptwXY_simpleWrite( u, ff, fmt );
    if( verbose ) ptwXY_simpleWrite( u, stdout, fmt );
    if( verbose ) printf( "\n" );
    if( verbose ) ptwXY_showInteralStructure( u, stdout, 1 );

    ptwXY_free( u );
    fclose( ff );

    exit( EXIT_SUCCESS );
}
/*
****************************************************************
*/
void printMsg( char const *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    vfprintf( stderr, fmt, args );
    fprintf( stderr, "\n" );
    va_end( args );
    exit( EXIT_FAILURE );
}
