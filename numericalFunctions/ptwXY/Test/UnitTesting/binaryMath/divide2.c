/*
# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>
*/
/*
    This test was the reason ClosestAllowXFactor as implemented. The implementation of ClosestAllowXFactor sures that no
    point will be added by some of the bi-secting routines that is too close to its neighboring points.
*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ptwXY.h>
#include <nf_utilities.h>
#include <ptwXY_utilities.h>

struct XYData {
    int points;
    double *XYs;
};

static int verbose = 0;
static double biSectionMax = 5., accuracy = 1e-3;
static char *fmtXY = "%18.11e %18.11e\n";

/*
************************************************************
*/
int main( int argc, char **argv ) {

    int i, iarg, echo = 0, status = EXIT_SUCCESS;
    double x, y;
    double numerator[] = {  1.99472656250000013e+00, -8.50773902802757220e-01, 1.99650065104166696e+00, -5.62986358469515835e-01,
                            1.99827473958333357e+00, -2.76792310494329286e-01, 1.99999999999999911e+00,  0.00000000000000000e+00,
                            2.00004882812500018e+00,  7.81189919780445052e-03, 2.00182291666666679e+00,  2.90829924402260076e-01,
                            2.00359700520833339e+00,  5.72265414633420733e-01 };
    double denominator[] = { 1.99472656250000013e+00, -8.50773902802757220e-01, 1.99517008463541679e+00, -7.78677426224021474e-01,
                             1.99561360677083366e+00, -7.06680714794401865e-01, 1.99605712890625031e+00, -6.34783711281670548e-01,
                             1.99650065104166696e+00, -5.62986358469515835e-01, 1.99694417317708361e+00, -4.91288599159361183e-01,
                             1.99738769531250027e+00, -4.19690376167636714e-01, 1.99783121744791692e+00, -3.48191632329417189e-01,
                             1.99827473958333357e+00, -2.76792310494329286e-01, 1.99871826171875022e+00, -2.05492353531553817e-01,
                             1.99916178385416687e+00, -1.34291704325278261e-01, 1.99960530598958353e+00, -6.31903057765157428e-02,
                             2.00000000000000000e+00,  0.00000000000000000e+00, 2.00004882812500018e+00,  7.81189919780445052e-03,
                             2.00049235026041661e+00,  7.87149676625631400e-02, 2.00093587239583348e+00,  1.49518956666270242e-01,
                             2.00137939453125036e+00,  2.20223923241064767e-01, 2.00182291666666679e+00,  2.90829924402260076e-01,
                             2.00226643880208321e+00,  3.61337017148343875e-01, 2.00270996093750009e+00,  4.31745258461432968e-01,
                             2.00315348307291696e+00,  5.02054705307273252e-01, 2.00359700520833339e+00,  5.72265414633420733e-01 };
    int nN = sizeof( numerator ) / 2 / sizeof( numerator[0] );
    int nD = sizeof( denominator ) / 2 / sizeof( denominator[0] );
    ptwXYPoints *ptwXY1, *ptwXY2, *division;
    nfu_status status_nf;

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

    if( ( ptwXY1 = ptwXY_create( ptwXY_interpolationLinLin, biSectionMax, accuracy, 100, 10, nN, numerator, &status_nf, 0 ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: ptwXY1 creation status = %d: %s", __FILE__, status_nf, nfu_statusMessage( status_nf ) );
    if( ( ptwXY2 = ptwXY_create( ptwXY_interpolationLinLin, biSectionMax, accuracy, 100, 10, nD, denominator, &status_nf, 0 ) ) == NULL )
        nfu_printErrorMsg( "ERROR %s: ptwXY2 creation status = %d: %s", __FILE__, status_nf, nfu_statusMessage( status_nf ) );

    if( verbose ) {
        printf( "# biSectionMax = %.12e\n", biSectionMax );
        printf( "# accuracy = %.12e\n", accuracy );
        printf( "# length = %d\n", (int) ptwXY_length( ptwXY1 ) );
        ptwXY_simpleWrite( ptwXY1, stdout, fmtXY );
        printf( "\n\n" );
        printf( "# length = %d\n", (int) ptwXY_length( ptwXY2 ) );
        ptwXY_simpleWrite( ptwXY2, stdout, fmtXY );
        printf( "\n\n" );
    }

    if( ( division = ptwXY_div_ptwXY( ptwXY1, ptwXY2, &status_nf, 1 ) ) == NULL ) 
        nfu_printErrorMsg( "ERROR %s: ptwXY2 division status = %d: %s", __FILE__, status_nf, nfu_statusMessage( status_nf ) );

    if( verbose ) {
        printf( "# length = %d\n", (int) ptwXY_length( division ) );
        ptwXY_simpleWrite( division, stdout, fmtXY );
        printf( "\n\n" );
    }

    for( i = 0; i < division->length; i++ ) {
        ptwXY_getXYPairAtIndex( division, i, &x, &y );
        if( !isfinite( y ) ) status++;
    }

    ptwXY_free( ptwXY1 );
    ptwXY_free( ptwXY2 );
    ptwXY_free( division );

    return( status );
}
