/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ptwXY_utilities.h"

static char const *fmtX = " %20.12e";
static char const *fmtXY = " %20.12e %20.12e";

/*
********************************************************
*/
int nfu_cmpDoubles( double d1, double d2, double espilon ) {

    double diff = d2 - d1, dMax = fabs( d1 );

    if( diff == 0 ) return( 0 );
    if( fabs( d2 ) > dMax ) dMax = fabs( d2 );
    if( fabs( diff ) > espilon * dMax ) {
#if 0
        double rDiff = diff / dMax;                 /* For debugging. */
#endif
        return( 1 );
    }
    return( 0 );
}
/*
********************************************************
*/
int nfu_ptwXY_cmp( ptwXYPoints *p1, ptwXYPoints *p2, int verbose, double frac ) {
/*
c   This rotuine compares the y-values for p1 and p2 for each x-value in their union.
*/
    int64_t i;
    int errCount = 0;
    nfu_status status;
    ptwXYPoints *u;
    ptwXYPoint *point;
    double y1,  y2, ratio;

    if( ( u = ptwXY_union( NULL, p1, p2, 0 ) ) == NULL ) return( -2 );

    for( i = 0; i < ptwXY_length( NULL, u ); i++ ) {
        point = ptwXY_getPointAtIndex_Unsafely( u, i );
        if( ( status = ptwXY_getValueAtX( NULL, p1, point->x, &y1 ) ) != nfu_Okay ) {
            nfu_printMsg( "Error nfu_ptwXY_cmp: for p1 ptwXY_getValueAtX return status = %d: %s", status, nfu_statusMessage( status ) );
            ptwXY_free( u );
            return( -1 );
        }
        if( ( status = ptwXY_getValueAtX( NULL, p2, point->x, &y2 ) ) != nfu_Okay ) {
            nfu_printMsg( "Error nfu_ptwXY_cmp: for p2 ptwXY_getValueAtX return status = %d: %s", status, nfu_statusMessage( status ) );
            ptwXY_free( u );
            return( -1 );
        }
        if( fabs( y2 - y1 ) > frac * ( fabs( y1 ) + fabs( y2 ) ) ) {
            errCount++;
            ratio = fabs( y1 );
            if( fabs( y2 ) > ratio ) ratio = fabs( y2 );
            ratio = ( y2 - y1 ) / ratio;
            if( verbose ) nfu_printMsg( "Error nfu_ptwXY_cmp: difference at x = %16.8e, y1 = %16.8e, y2 = %16.8e, y2 - y1 = %16.8e, ratio = %.2e", 
                point->x, y1, y2, y2 - y1, ratio );
        }
    }

    ptwXY_free( u );
    return( errCount );
}
/*
********************************************************
*/
void nfu_printSMRError( statusMessageReporting *smr, char const *file, int line, char const *function, char const *fmt, ... ) {

    va_list args;

    va_start( args, fmt );
    smr_vsetReportError( smr, NULL, file, line, function, nfu_SMR_libraryID, 0, fmt, &args ); /* FIXME 0 needs to be something else. */
    va_end( args );

    smr_write( smr, stderr, 1 );
    exit( EXIT_FAILURE );
}
/*
************************************************************
*/
void nfu_printXYDataOnVerbosity( int verbose, ptwXYPoints *data ) {

    if( !verbose ) return;
    printf( "# length = %d\n", (int) ptwXY_length( NULL, data ) );
    printf( "# interpolation = %s\n", ptwXY_getInterpolationString( data ) );
    ptwXY_simplePrint( data, fmtXY );
}
/*
************************************************************
*/
void nfu_printXDataOnVerbosity( int verbose, ptwXPoints *data ) {

    if( !verbose ) return;
    printf( "# length = %d\n", (int) ptwX_length( NULL, data ) );
    ptwX_simplePrint( NULL, data, fmtX );
}
