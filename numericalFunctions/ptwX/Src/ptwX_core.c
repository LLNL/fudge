/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ptwX.h"

/*
************************************************************
*/
ptwXPoints *ptwX_new( int64_t size, nfu_status *status ) {

    ptwXPoints *ptwX = (ptwXPoints *) nfu_calloc( sizeof( ptwXPoints ), 1 );

    *status = nfu_mallocError;
    if( ptwX == NULL ) return( NULL );
    ptwX_setup( ptwX, size );
    if( ( *status = ptwX->status ) != nfu_Okay ) ptwX = nfu_free( ptwX );
    return( ptwX );
}
/*
************************************************************
*/
nfu_status ptwX_setup( ptwXPoints *ptwX, int64_t size ) {

    ptwX->status = nfu_Okay;
    ptwX->length = 0;
    ptwX->allocatedSize = 0;
    ptwX->mallocFailedSize = 0;
    ptwX->points = NULL;
    ptwX_reallocatePoints( ptwX, size, 0 );
    return( ptwX->status );
}
/*
************************************************************
*/
ptwXPoints *ptwX_create( int64_t size, int64_t length, double *xs, nfu_status *status ) {

    ptwXPoints *ptwX = ptwX_new( size, status );

    if( ( *status = ptwX_setData( ptwX, length, xs ) ) != nfu_Okay ) ptwX = ptwX_free( ptwX );
    return( ptwX );
}
/*
************************************************************
*/
nfu_status ptwX_copy( ptwXPoints *dest, ptwXPoints *src ) {

    if( dest->status == nfu_Okay ) return( dest->status );
    if( src->status == nfu_Okay ) return( src->status );
    ptwX_clear( dest );
    return( ptwX_setData( dest, src->length, src->points ) );
}
/*
************************************************************
*/
ptwXPoints *ptwX_clone( ptwXPoints *ptwX, nfu_status *status ) {

    return( ptwX_slice( ptwX, 0, ptwX->length, status ) );
}
/*
************************************************************
*/
ptwXPoints *ptwX_slice( ptwXPoints *ptwX, int64_t index1, int64_t index2, nfu_status *status ) {

    int64_t i, j, length;
    ptwXPoints *n;

    *status = nfu_badSelf;
    if( ptwX->status != nfu_Okay ) return( NULL );
    *status = nfu_badIndex;
    if( index1 < 0 ) return( NULL );
    if( index2 < index1 ) return( NULL );
    if( index2 > ptwX->length ) return( NULL );
    length = ( index2 - index1 );
    if( ( n = ptwX_new( length, status ) ) == NULL ) return( n );
    *status = n->status;
    for( j = 0, i = index1; i < index2; i++ ) n->points[j] = ptwX->points[i];
    n->length = length;
    return( n );
}
/*
************************************************************
*/
nfu_status ptwX_reallocatePoints( ptwXPoints *ptwX, int64_t size, int forceSmallerResize ) {

    if( size < ptwX_minimumSize ) size = ptwX_minimumSize;                        /* ptwX_minimumSize must be > 0 for other routines to work properly. */
    if( size < ptwX->length ) size = ptwX->length;
    if( size != ptwX->allocatedSize ) {
        if( size > ptwX->allocatedSize ) {                                         /* Increase size of allocated points. */
             ptwX->points = nfu_realloc( size * sizeof( double ), ptwX->points ); }
        else if( ( ptwX->allocatedSize > 2 * size ) || forceSmallerResize ) {      /* Decrease size, if at least 1/2 size reduction or if forced to. */
            ptwX->points = nfu_realloc( size * sizeof( double ), ptwX->points );
        }
        if( ptwX->points == NULL ) {
            ptwX->mallocFailedSize = size;
            size = 0;
            ptwX->status = nfu_mallocError;
        }
        ptwX->allocatedSize = size;
    }

    return( ptwX->status );
}
/*
************************************************************
*/
nfu_status ptwX_clear( ptwXPoints *ptwX ) {

    ptwX->length = 0;
    if( ptwX->status != nfu_Okay ) return( ptwX->status );

    return( nfu_Okay );
}
/*
************************************************************
*/
nfu_status ptwX_release( ptwXPoints *ptwX ) {

    ptwX->length = 0;
    ptwX->allocatedSize = 0;
    ptwX->points = nfu_free( ptwX->points );

    return( nfu_Okay );
}
/*
************************************************************
*/
ptwXPoints *ptwX_free( ptwXPoints *ptwX ) {

    ptwX_release( ptwX );
    return( (ptwXPoints *) nfu_free( ptwX ) );
}
/*
************************************************************
*/
int64_t ptwX_length( ptwXPoints *ptwX ) {

    return( ptwX->length );
}
/*
************************************************************
*/
nfu_status ptwX_setData( ptwXPoints *ptwX, int64_t length, double *xs ) {

    int64_t  i;

    if( ptwX->status != nfu_Okay ) return( ptwX->status );

    if( length > ptwX->allocatedSize ) {
        ptwX_reallocatePoints( ptwX, length, 0 );
        if( ptwX->status != nfu_Okay ) return( ptwX->status );
    }
    for( i = 0; i < length; i++ ) ptwX->points[i] = xs[i];
    ptwX->length = length;

    return( ptwX->status );
}
/*
************************************************************
*/
double *ptwX_getPointAtIndex( ptwXPoints *ptwX, int64_t index ) {

    if( ptwX->status != nfu_Okay ) return( NULL );
    if( ( index < 0 ) || ( index >= ptwX->length ) ) return( NULL );
    return( &(ptwX->points[index]) );
}
/*
************************************************************
*/
double ptwX_getPointAtIndex_Unsafely( ptwXPoints *ptwX, int64_t index ) {

    return( ptwX->points[index] );
}
/*
************************************************************
*/
nfu_status ptwX_setPointAtIndex( ptwXPoints *ptwX, int64_t index, double x ) {

    nfu_status status;

    if( ptwX->status != nfu_Okay ) return( ptwX->status);
    if( ( index < 0 ) || ( index > ptwX->length ) ) return( nfu_badIndex );
    if( index == ptwX->allocatedSize ) {
        if( ( status = ptwX_reallocatePoints( ptwX, ptwX->allocatedSize + 10, 0 ) ) != nfu_Okay ) return( status );
    }
    ptwX->points[index] = x;
    if( index == ptwX->length ) ptwX->length++;
    return( nfu_Okay );
}
/*
************************************************************
*/
int ptwX_ascendingOrder( ptwXPoints *ptwX ) {

    int order = 1;
    int64_t i;
    double x1, x2;

    if( ptwX->length < 2 ) return( 0 );

    if( ( x1 = ptwX->points[0] ) < ( x2 = ptwX->points[1] ) ) {     /* Check for ascending order. */
        for( i = 2; i < ptwX->length; i++ ) {
            x1 = x2;
            x2 = ptwX->points[i];
            if( x2 <= x1 ) return( 0 );
        } }
    else {
        if( x1 == x2 ) return( 0 );
        order = -1;                                                 /* Check for descending order. */
        for( i = 2; i < ptwX->length; i++ ) {
            x1 = x2;
            x2 = ptwX->points[i];
            if( x1 <= x2 ) return( 0 );
        }
    }
    return( order );
}
