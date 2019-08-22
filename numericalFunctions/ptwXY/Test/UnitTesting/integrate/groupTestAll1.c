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

#include <nfut_utilities.h>
#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>

static int verbose = 0;
static char *fmtX = "%18.11e\n", *fmtXY = "%18.11e %18.11e\n";

static void compareResults( statusMessageReporting *smr, ptwXPoints *groupedData, const char * const fileName );
static void compareNone2dx( statusMessageReporting *smr, ptwXPoints *groupBoundaries, ptwXPoints *fluxGrouped_None, ptwXPoints *groupedData, char *sData );
static void compareNone2Normed( statusMessageReporting *smr, ptwXPoints *groupNorm, ptwXPoints *groupedData_None, ptwXPoints *groupedData_Norm, char *sData );
static void compareDoubles( int index, double d1, double d2, double eps, char *sData, char *s );
static ptwXPoints *getGroupBoundaries( statusMessageReporting *smr );
static ptwXYPoints *getFluxData( statusMessageReporting *smr );
static ptwXYPoints *getCrossSectionData( statusMessageReporting *smr );
static ptwXYPoints *getMultiplicityData( statusMessageReporting *smr );
static void writeXYDataOnVerbosity( ptwXYPoints *data, const char * const fileName );
static void writeXDataOnVerbosity( statusMessageReporting *smr, ptwXPoints *data, const char * const fileName );
/*
************************************************************
*/
int main( int argc, char **argv ) {

    int iarg, errCount = 0, echo = 0;
    ptwXYPoints *flux, *crossSection, *multiplicity;
    ptwXPoints *groupBoundaries, *groupedData, *fluxGrouped_None, *crossSectionGrouped_None, *multiplicityGrouped_None;
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

    groupBoundaries = getGroupBoundaries( &smr );
    writeXDataOnVerbosity( &smr, groupBoundaries, "groupBoundaries.dat" );

    flux = getFluxData( &smr );
    writeXYDataOnVerbosity( flux, "flux.dat" );

    crossSection = getCrossSectionData( &smr );
    writeXYDataOnVerbosity( crossSection, "crossSection.dat" );

    multiplicity = getMultiplicityData( &smr );
    writeXYDataOnVerbosity( multiplicity, "multiplicity.dat" );

/*
*   ptwXY_groupOneFunction testing.
*/
    if( ( fluxGrouped_None = ptwXY_groupOneFunction( &smr, flux, groupBoundaries, ptwXY_group_normType_none, NULL ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    writeXDataOnVerbosity( &smr, fluxGrouped_None, "flux_None.dat" );
    compareResults( &smr, fluxGrouped_None, "Data/flux_None.dat" );

    if( ( groupedData = ptwXY_groupOneFunction( &smr, flux, groupBoundaries, ptwXY_group_normType_dx, NULL ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    compareNone2dx( &smr, groupBoundaries, fluxGrouped_None, groupedData, "flux" );
    writeXDataOnVerbosity( &smr, groupedData, "flux_dx.dat" );
    ptwX_free( groupedData );

    if( ( groupedData = ptwXY_groupOneFunction( &smr, flux, groupBoundaries, ptwXY_group_normType_norm, fluxGrouped_None ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    compareNone2Normed( &smr, fluxGrouped_None, fluxGrouped_None, groupedData, "flux" );
    writeXDataOnVerbosity( &smr, groupedData, "flux_norm.dat" );
    ptwX_free( groupedData );

/*
*   ptwXY_groupTwoFunctions testing.
*/
    if( ( crossSectionGrouped_None = ptwXY_groupTwoFunctions( &smr, crossSection, flux, groupBoundaries, ptwXY_group_normType_none, NULL ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    writeXDataOnVerbosity( &smr, crossSectionGrouped_None, "crossSection_None.dat" );
    compareResults( &smr, crossSectionGrouped_None, "Data/crossSection_None.dat" );

    if( ( groupedData = ptwXY_groupTwoFunctions( &smr, crossSection, flux, groupBoundaries, ptwXY_group_normType_dx, NULL ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    compareNone2dx( &smr, groupBoundaries, crossSectionGrouped_None, groupedData, "cross section" );
    writeXDataOnVerbosity( &smr, groupedData, "crossSection_dx.dat" );
    ptwX_free( groupedData );

    if( ( groupedData = ptwXY_groupTwoFunctions( &smr, crossSection, flux, groupBoundaries, ptwXY_group_normType_norm, fluxGrouped_None ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    compareNone2Normed( &smr, fluxGrouped_None, crossSectionGrouped_None, groupedData, "cross section" );
    writeXDataOnVerbosity( &smr, groupedData, "crossSection_norm.dat" );
    ptwX_free( groupedData );

/*
*   ptwXY_groupThreeFunctions testing.
*/
    if( ( multiplicityGrouped_None = ptwXY_groupThreeFunctions( &smr, crossSection, multiplicity, flux, groupBoundaries, ptwXY_group_normType_none, NULL ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    writeXDataOnVerbosity( &smr, multiplicityGrouped_None, "multiplicity_None.dat" );
    compareResults( &smr, multiplicityGrouped_None, "Data/multiplicity_None.dat" );

    if( ( groupedData = ptwXY_groupThreeFunctions( &smr, crossSection, multiplicity, flux, groupBoundaries, ptwXY_group_normType_dx, NULL ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    compareNone2dx( &smr, groupBoundaries, multiplicityGrouped_None, groupedData, "multiplcity" );
    writeXDataOnVerbosity( &smr, groupedData, "multiplicity_dx.dat" );
    ptwX_free( groupedData );

    if( ( groupedData = ptwXY_groupThreeFunctions( &smr, crossSection, multiplicity, flux, groupBoundaries, ptwXY_group_normType_norm, fluxGrouped_None ) ) == NULL )
        nfut_printSMRErrorExit2p( &smr, "Via." );
    compareNone2Normed( &smr, fluxGrouped_None, multiplicityGrouped_None, groupedData, "multiplicity" );
    writeXDataOnVerbosity( &smr, groupedData, "multiplicity_norm.dat" );
    ptwX_free( groupedData );

    ptwX_free( groupBoundaries );
    ptwXY_free( flux );
    ptwX_free( fluxGrouped_None );
    ptwXY_free( crossSection );
    ptwX_free( crossSectionGrouped_None );
    ptwXY_free( multiplicity );
    ptwX_free( multiplicityGrouped_None );

    exit( errCount );
}
/*
************************************************************
*/
static void compareResults( statusMessageReporting *smr, ptwXPoints *groupedData, const char * const fileName ) {

    int i, j, n = (int) ptwX_length( smr, groupedData );
    ptwXPoints *old;
    FILE *f;
    double x1, x2, s, d, r;

    if( ( old = ptwX_new( smr, n ) ) == NULL ) 
        nfut_printSMRErrorExit2p( smr, "Via." );
    if( ( f = fopen( fileName, "r" ) ) == NULL ) nfu_printErrorMsg( "ERROR %s: could not open file '%s' for comparison", __FILE__, fileName );
    for( i = 0; i < n; i++ ) {
        if( ( j = fscanf( f, "%le", &d ) ) != 1 ) nfu_printErrorMsg( "ERROR %s: reading double at line %d of file '%s' for comparison", 
                __FILE__, i, fileName );
        old->points[i] = d;
    }
    fclose( f );
    for( i = 0; i < n; i++ ) {
        x1 = groupedData->points[i];
        x2 = old->points[i];
        s = fabs( x1 ) + fabs( x2 );
        d = x1 - x2;
        r = d;
        if( s != 0 ) r /= s;
        if( fabs( r ) > 1e-10 ) printf( "ERROR %s: compare for '%s' at line %3d, %e %e %e %e %e\n", __FILE__, fileName, i, x1, x2, s, d, r );
    }
    ptwX_free( old );
}
/*
************************************************************
*/
static void compareNone2dx( statusMessageReporting *smr, ptwXPoints *groupBoundaries, ptwXPoints *groupedData_None, ptwXPoints *groupedData, char *sData ) {

    int i, n = ptwX_length( smr, groupedData_None );
    double dx, d1, d2;

    for( i = 0; i < n; i++ ) {
        dx = ptwX_getPointAtIndex_Unsafely( groupBoundaries, i + 1 ) - ptwX_getPointAtIndex_Unsafely( groupBoundaries, i );
        d1 = ptwX_getPointAtIndex_Unsafely( groupedData_None, i );
        d2 = ptwX_getPointAtIndex_Unsafely( groupedData, i );
        compareDoubles( i, d1 / dx, d2, 1e-15, sData, "compareNone2dx" );
    }
}
/*
************************************************************
*/
static void compareNone2Normed( statusMessageReporting *smr, ptwXPoints *groupNorm, ptwXPoints *groupedData_None, ptwXPoints *groupedData_Norm, char *sData ) {

    int i, n = ptwX_length( smr, groupedData_None );
    double norm, d1, d2;

    for( i = 0; i < n; i++ ) {
        norm = ptwX_getPointAtIndex_Unsafely( groupNorm, i );
        d1 = ptwX_getPointAtIndex_Unsafely( groupedData_None, i );
        d2 = ptwX_getPointAtIndex_Unsafely( groupedData_Norm, i );
        compareDoubles( i, d1 / norm, d2, 1e-15, sData, "compareNone2Normed" );
    }
}
/*
************************************************************
*/
static void compareDoubles( int index, double d1, double d2, double eps, char *sData, char *funcName ) {

    double s, d, r;

    s = 0.5 * ( fabs( d1 ) + fabs( d2 ) );
    d = d2 - d1;
    r = d;
    if( s != 0 ) r /= s;
    if( fabs( r ) > eps ) printf( "ERROR %s: %s compare for '%s' at index %3d, %e %e %e %e %e\n", __FILE__, sData, funcName, index, d1, d2, s, d, r );
}
/*
************************************************************
*/
static ptwXPoints *getGroupBoundaries( statusMessageReporting *smr ) {

    double data[] = {
        1.30680e-09, 2.09080e-08, 1.30680e-07, 3.34530e-07, 1.17610e-06, 2.09080e-06, 5.65780e-06, 1.30680e-05, 2.07460e-05, 5.12300e-05,
        1.02450e-04, 2.09080e-04, 3.81050e-04, 5.65780e-04, 7.15580e-04, 1.05850e-03, 1.30680e-03, 1.88170e-03, 2.94020e-03, 3.34530e-03,
        4.23390e-03, 5.76280e-03, 7.52700e-03, 1.02450e-02, 1.51060e-02, 2.09080e-02, 2.64620e-02, 3.26690e-02, 3.95300e-02, 7.00200e-02,
        9.89090e-02, 1.30680e-01, 1.81950e-01, 2.07460e-01, 2.41700e-01, 2.70970e-01, 2.94020e-01, 3.34530e-01, 3.77650e-01, 5.12300e-01,
        6.32470e-01, 7.52700e-01, 8.83370e-01, 1.02450e+00, 1.17610e+00, 1.33810e+00, 1.51060e+00, 1.69360e+00, 2.09080e+00, 2.30510e+00,
        2.52990e+00, 2.74110e+00, 3.01080e+00, 3.53350e+00, 4.06880e+00, 4.39600e+00, 4.70440e+00, 4.99080e+00, 5.35250e+00, 5.65780e+00,
        6.04250e+00, 6.36660e+00, 6.73670e+00, 7.15580e+00, 7.54790e+00, 7.90960e+00, 8.32150e+00, 8.78670e+00, 9.17670e+00, 9.66480e+00,
        1.01200e+01, 1.05850e+01, 1.10120e+01, 1.15470e+01, 1.19930e+01, 1.24990e+01, 1.30680e+01, 1.35420e+01, 1.38630e+01, 1.41340e+01,
        1.44070e+01, 1.46830e+01, 1.51860e+01, 1.57540e+01, 1.63340e+01, 1.69230e+01, 1.81340e+01, 2.00000e+01 };
    int nData = sizeof( data ) / sizeof( double );
    ptwXPoints *groupBoundaries;

    if( ( groupBoundaries = ptwX_create( smr, nData, nData, data ) ) == NULL ) 
        nfut_printSMRErrorExit2p( smr, "Via." );
    return( groupBoundaries );
}
/*
************************************************************
*/
static ptwXYPoints *getFluxData( statusMessageReporting *smr ) {

    double data[] = {  
        0.0000000e+00, 5.4598150e+01, 5.9634465e+00, 5.6305711e+00, 8.2124320e+00, 2.3903851e+00, 9.6171596e+00, 1.3997871e+00, 1.0649959e+01, 9.4447384e-01,
        1.1453293e+01, 6.9547558e-01, 1.2120111e+01, 5.3946075e-01, 1.2689944e+01, 4.3419372e-01, 1.3187622e+01, 3.5920720e-01, 1.3620180e+01, 3.0463531e-01,
        1.4009880e+01, 2.6260684e-01, 1.4364428e+01, 2.2942870e-01, 1.4689668e+01, 2.0269278e-01, 1.4989942e+01, 1.8078369e-01, 1.5268929e+01, 1.6255565e-01,
        1.5529441e+01, 1.4719790e-01, 1.5773780e+01, 1.3411483e-01, 1.5998524e+01, 1.2311016e-01, 1.6211129e+01, 1.1353226e-01, 1.6412837e+01, 1.0513506e-01,
        1.6604714e+01, 9.7724282e-02, 1.6787655e+01, 9.1145600e-02, 1.6962472e+01, 8.5273285e-02, 1.7129855e+01, 8.0005574e-02, 1.7290411e+01, 7.5258733e-02,
        1.7444603e+01, 7.0965397e-02, 1.7592983e+01, 6.7065295e-02, 1.7735973e+01, 6.3509804e-02, 1.7873951e+01, 6.0257736e-02, 1.8007253e+01, 5.7274143e-02,
        1.8136190e+01, 5.4528889e-02, 1.8261037e+01, 5.1996154e-02, 1.8382047e+01, 4.9653589e-02, 1.8607903e+01, 4.5560014e-02, 1.8821502e+01, 4.1999560e-02,
        1.9024105e+01, 3.8879888e-02, 1.9216791e+01, 3.6128174e-02, 1.9400467e+01, 3.3686629e-02, 1.9575956e+01, 3.1508211e-02, 1.9743954e+01, 2.9554878e-02,
        1.9905077e+01, 2.7795351e-02, 2.0059790e+01, 2.6204483e-02, 2.0208653e+01, 2.4759787e-02, 2.0352091e+01, 2.3443133e-02, 2.0490488e+01, 2.2239172e-02,
        2.0624179e+01, 2.1134888e-02, 2.0753480e+01, 2.0119061e-02, 2.0878669e+01, 1.9182082e-02, 2.1000000e+01, 1.8315639e-02 };
    int nData = sizeof( data ) / ( 2 * sizeof( double ) );
    ptwXYPoints *flux;

    if( ( flux = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, 6, 1e-3, 10, 10, nData, data, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( smr, "Via." );
    return( flux );
}
/*
************************************************************
*/
static ptwXYPoints *getCrossSectionData( statusMessageReporting *smr ) {

#include "Data/crossSection.dat.h"     /* This lines include data. */
    int nData = sizeof( data ) / ( 2 * sizeof( double ) );
    ptwXYPoints *crossSection;

    if( ( crossSection = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, 6, 1e-3, 10, 10, nData, data, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( smr, "Via." );
    return( crossSection);
}
/*
************************************************************
*/
static ptwXYPoints *getMultiplicityData( statusMessageReporting *smr ) {

    double data[] = {  
            1.00000000e-11,  3.95456000e+00, 2.00000000e-02,  3.95456000e+00, 3.00000000e-02,  3.95403000e+00, 4.00000000e-02,  3.95466000e+00,
            5.00000000e-02,  3.95591000e+00, 6.00000000e-02,  3.95750000e+00, 8.00000000e-02,  3.96126000e+00, 9.00000000e-02,  3.96508000e+00,
            1.00000000e-01,  3.96818000e+00, 1.20000000e-01,  3.97448000e+00, 1.50000000e-01,  3.98320000e+00, 2.00000000e-01,  3.99771000e+00,
            3.00000000e-01,  4.02781000e+00, 4.00000000e-01,  4.05829000e+00, 5.00000000e-01,  4.08896000e+00, 6.00000000e-01,  4.11958000e+00,
            7.00000000e-01,  4.15106000e+00, 8.00000000e-01,  4.17956000e+00, 9.00000000e-01,  4.20731000e+00, 1.00000000e+00,  4.23581000e+00,
            1.10000000e+00,  4.26337000e+00, 1.20000000e+00,  4.28905000e+00, 1.30000000e+00,  4.31651000e+00, 1.40000000e+00,  4.34344000e+00,
            1.50000000e+00,  4.36066000e+00, 1.60000000e+00,  4.38702000e+00, 1.70000000e+00,  4.41008000e+00, 1.80000000e+00,  4.43265000e+00,
            1.90000000e+00,  4.45567000e+00, 2.00000000e+00,  4.47668000e+00, 2.50000000e+00,  4.58012000e+00, 3.00000000e+00,  4.67677000e+00,
            3.50000000e+00,  4.76295000e+00, 4.00000000e+00,  4.84127000e+00, 5.00000000e+00,  4.97391000e+00, 6.00000000e+00,  5.03102000e+00,
            7.00000000e+00,  4.96053000e+00, 8.00000000e+00,  4.75777000e+00, 9.00000000e+00,  4.50265000e+00, 1.00000000e+01,  4.32857000e+00,
            1.10000000e+01,  4.34556000e+00, 1.20000000e+01,  4.45396000e+00, 1.30000000e+01,  4.59275000e+00, 1.40000000e+01,  4.71241000e+00,
            1.45000000e+01,  4.76865000e+00, 1.50000000e+01,  4.82063000e+00, 1.60000000e+01,  4.89104000e+00, 1.70000000e+01,  4.92432000e+00,
            1.80000000e+01,  4.90964000e+00, 1.90000000e+01,  4.87693000e+00, 2.00000000e+01,  4.86958000e+00 };
    int nData = sizeof( data ) / ( 2 * sizeof( double ) );
    ptwXYPoints *multiplicity;

    if( ( multiplicity = ptwXY_create( smr, ptwXY_interpolationLinLin, NULL, 6, 1e-3, 10, 10, nData, data, 0 ) ) == NULL ) 
        nfut_printSMRErrorExit2p( smr, "Via." );
    return( multiplicity );
}
/*
************************************************************
*/
static void writeXYDataOnVerbosity( ptwXYPoints *data, const char * const fileName ) {

    FILE *f;

    if( !verbose ) return;
    if( ( f = fopen( fileName, "w" ) ) == NULL ) nfu_printErrorMsg( "ERROR %s: could not open file %s\n", __FILE__, fileName );
    fprintf( f, "# length = %d\n", (int) ptwXY_length( NULL, data ) );
    ptwXY_simpleWrite( data, f, fmtXY );
    fclose( f );
}
/*
************************************************************
*/
static void writeXDataOnVerbosity( statusMessageReporting *smr, ptwXPoints *data, const char * const fileName ) {

    int64_t i;
    FILE *f;

    if( !verbose ) return;
    if( ( f = fopen( fileName, "w" ) ) == NULL ) nfu_printErrorMsg( "ERROR %s: could not open file %s\n", __FILE__, fileName );
    fprintf( f, "# length = %d\n", (int) ptwX_length( smr, data ) );
    for( i = 0; i < ptwX_length( smr, data ); i++ ) fprintf( f, fmtX, ptwX_getPointAtIndex_Unsafely( data, i ) );
    fclose( f );
}
