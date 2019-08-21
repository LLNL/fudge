/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <math.h>
#include <float.h>

#include "ptwXY.h"

/*
************************************************************
*/
nfu_status ptwXY_abs( ptwXYPoints *ptwXY ) {

    int64_t i, nonOverflowLength = ptwXY_getNonOverflowLength( ptwXY );
    ptwXYPoint *p;
    ptwXYOverflowPoint *o, *overflowHeader = &(ptwXY->overflowHeader);

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );

    for( i = 0, p = ptwXY->points; i < nonOverflowLength; i++, p++ ) p->y = fabs( p->y );
    for( o = overflowHeader->next; o != overflowHeader; o = o->next ) o->point.y = fabs( o->point.y );
    return( ptwXY->status );
}
/*
************************************************************
*/
nfu_status ptwXY_neg( ptwXYPoints *ptwXY ) {

    int64_t i, nonOverflowLength = ptwXY_getNonOverflowLength( ptwXY );
    ptwXYPoint *p;
    ptwXYOverflowPoint *o, *overflowHeader = &(ptwXY->overflowHeader);

    if( ptwXY->status != nfu_Okay ) return( ptwXY->status );

    for( i = 0, p = ptwXY->points; i < nonOverflowLength; i++, p++ ) p->y = -p->y;
    for( o = overflowHeader->next; o != overflowHeader; o = o->next ) o->point.y = -o->point.y;
    return( ptwXY->status );
}
