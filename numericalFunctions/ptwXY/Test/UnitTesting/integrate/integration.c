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
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>

static int verbose = 0;

static int integrate( statusMessageReporting *smr, double x1, double y1, double x2, double y2, double nn, double gn, double ng, double gg );
static int integrate2( statusMessageReporting *smr, char *label, ptwXY_interpolation interpolation, 
        double x1, double y1, double x2, double y2, double answer );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0;
    double flat[4] = { 1., 3., 4., 3. }, slope[4] = { 1., 3., 4., 9. };
    statusMessageReporting smr;

    smr_initialize( &smr, smr_status_Ok );

    for( iarg = 1; iarg < argc; iarg++ ) {
        if( strcmp( "-v", argv[iarg] ) == 0 ) {
            verbose = 1; }
        else if( strcmp( "-e", argv[iarg] ) == 0 ) {
            echo = 1; }
        else {
            nfu_printErrorMsg( "ERROR %s: invalid input option '%s'", __FILE__, argv[iarg] );
        }
    }
    if( echo ) printf( "%s\n", __FILE__ );

    errCount += integrate( &smr, flat[0], flat[1], flat[2], flat[3], 9., 9., 9., 9. );
    errCount += integrate( &smr, slope[0], slope[1], slope[2], slope[3], 18., 20.015744631999329, 16.384306079283071, 1.84102344129745710e+01 );
    slope[2] = slope[0] * ( 1 + 1e-3 );
    errCount += integrate( &smr, slope[0], slope[1], slope[2], slope[3], 0.0059999999999993392, 0.0060004997501579282, 0.0054614353597604226, 5.46192534635673351e-03 );
    slope[2] = slope[0] * ( 1 + 1e-4 * ( 1. + 1e-8 ) );
    errCount += integrate( &smr, slope[0], slope[1], slope[2], slope[3], 6.0000000600046732e-4, 6.000050057508588e-4, 5.4614354143796317e-4, 5.46148443428480447e-04 );
    slope[2] = slope[0] * ( 1 + 1e-4 * ( 1. - 1e-8 ) );
    errCount += integrate( &smr, slope[0], slope[1], slope[2], slope[3], 5.9999999399940052e-4, 6.0000499374931635e-4, 5.4614353051412137e-4, 5.46148432504442623e-04 );
    slope[2] = slope[0] * ( 1 + 1e-5 );
    errCount += integrate( &smr, slope[0], slope[1], slope[2], slope[3], 6.0000000000393072e-5, 6.0000050000143072e-5, 5.4614353597968027e-5, 5.46144026199979856e-05 );

    exit( errCount );
}
/*
************************************************************
*/
static int integrate( statusMessageReporting *smr, double x1, double y1, double x2, double y2, double nn, double gn, double ng, double gg ) {

    int errCount = 0;

    if( verbose ) printf( "x1 = %.17g\ny1 = %.17g\nx2 = %.17g\ny2 = %.17g\n", x1, y1, x2, y2 );
    errCount += integrate2( smr, "nn", ptwXY_interpolationLinLin, x1, y1, x2, y2, nn );
    errCount += integrate2( smr, "gn", ptwXY_interpolationLinLog, x1, y1, x2, y2, gn );
    errCount += integrate2( smr, "ng", ptwXY_interpolationLogLin, x1, y1, x2, y2, ng );
    errCount += integrate2( smr, "gg", ptwXY_interpolationLogLog, x1, y1, x2, y2, gg );
    return( errCount );
}
/*
************************************************************
*/
static int integrate2( statusMessageReporting *smr, char *label, ptwXY_interpolation interpolation, double x1, double y1, double x2, double y2, double answer ) {

    nfu_status status;
    double integral, s, d, r;
    char str[128], buffer[64], *e;

    if( ( status = ptwXY_f_integrate( smr, interpolation, x1, y1, x2, y2, &integral ) ) != nfu_Okay ) 
        smr_setReportError2( smr, nfu_SMR_libraryID, nfu_badSelf, "Via for label = '%s'", label );
    s = 0.5 * ( fabs( answer ) + fabs( integral ) );
    d = answer - integral;
    r = d;
    if( s != 0 ) r /= s;
    if( fabs( r ) > 1e-12 ) printf( "ERROR %s: %s compare, %e %e %e %e %e\n", __FILE__, label, answer, integral, s, d, r );

    sprintf( str, "%.17g", integral );
    e = strchr( str, 'e' );
    if( e != NULL ) {
        sprintf( buffer, " * 10^(%s)", &(e[1]) );
        strcpy( e, buffer );

    }
    if( verbose ) printf( "Integrate[ %s[x], { x, x1, x2 } ] / ( %s ) - 1\n", label, str );

    return( 0 );
}
