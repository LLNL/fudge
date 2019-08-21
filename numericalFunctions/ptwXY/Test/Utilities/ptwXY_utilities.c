/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ptwXY_utilities.h"

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
    double y1,  y2;

    if( ( u = ptwXY_union( p1, p2, &status, 0 ) ) == NULL ) return( -2 );

    for( i = 0; i < ptwXY_length( u ); i++ ) {
        point = ptwXY_getPointAtIndex_Unsafely( u, i );
        if( ( status = ptwXY_getValueAtX( p1, point->x, &y1 ) ) != nfu_Okay ) {
            nfu_printMsg( "Error nfu_ptwXY_cmp: for p1 ptwXY_getValueAtX return status = %d: %s", status, nfu_statusMessage( status ) );
            ptwXY_free( u );
            return( -1 );
        }
        if( ( status = ptwXY_getValueAtX( p2, point->x, &y2 ) ) != nfu_Okay ) {
            nfu_printMsg( "Error nfu_ptwXY_cmp: for p2 ptwXY_getValueAtX return status = %d: %s", status, nfu_statusMessage( status ) );
            ptwXY_free( u );
            return( -1 );
        }
        if( fabs( y2 - y1 ) > frac * ( fabs( y1 ) + fabs( y2 ) ) ) {
            errCount++;
            if( verbose ) nfu_printMsg( "Error nfu_ptwXY_cmp: difference at x = %g, y1 = %16g, y2 = %16g, y2 - y1 = %g", point->x, y1, y2, y2 - y1 );
        }
    }

    ptwXY_free( u );
    return( errCount );
}
