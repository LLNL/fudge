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

#include <nfut_utilities.h>
#include <nf_specialFunctions.h>

struct aFOfA {
    double a, f;
};

#include "gammaTest.dat"

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
        status = nf_gammaFunction( &smr, data[i1].a, &f );
        r = 1;
        if( data[i1].f != 0. ) {
            r = f / data[i1].f - 1; }
        else {
            if( f == 0. ) r = 0.;
        }
        if( ( fabs( r ) > 3e-15 ) || ( status != nfu_Okay ) ) {
            printf( "%.17e %.17e %.17e %+.3e  %d\n", data[i1].a, data[i1].f, f, r, status );
            if( status != nfu_Okay ) nfut_printSMRError2p( &smr, "Via." );
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
