/*
# <<BEGIN-copyright>>
# <<END-copyright>>
*/

#ifndef ptwX_h_included
#define ptwX_h_included

#include <stdio.h>
#include <stdint.h>

#include <nf_utilities.h>

#if defined __cplusplus
    extern "C" {
#endif

#define ptwX_minimumSize 10

typedef
    struct ptwXPoints_s {
        nfu_status status;
        int64_t length;
        int64_t allocatedSize;
        int64_t mallocFailedSize;
        double *points;
    } ptwXPoints;

/*
* Routines in ptwX_core.c
*/
ptwXPoints *ptwX_new( int64_t size, nfu_status *status );
nfu_status ptwX_setup( ptwXPoints *ptwX, int64_t size );
ptwXPoints *ptwX_create( int64_t size, int64_t length, double *xs, nfu_status *status );
nfu_status ptwX_copy( ptwXPoints *dest, ptwXPoints *src );
ptwXPoints *ptwX_clone( ptwXPoints *ptwX, nfu_status *status );
ptwXPoints *ptwX_slice( ptwXPoints *ptwX, int64_t index1, int64_t index2, nfu_status *status );
nfu_status ptwX_reallocatePoints( ptwXPoints *ptwX, int64_t size, int forceSmallerResize );
nfu_status ptwX_clear( ptwXPoints *ptwX );
nfu_status ptwX_release( ptwXPoints *ptwX );
ptwXPoints *ptwX_free( ptwXPoints *ptwX );

int64_t ptwX_length( ptwXPoints *ptwX );
nfu_status ptwX_setData( ptwXPoints *ptwX, int64_t length, double *xs );
double *ptwX_getPointAtIndex( ptwXPoints *ptwX, int64_t index );
double ptwX_getPointAtIndex_Unsafely( ptwXPoints *ptwX, int64_t index );
nfu_status ptwX_setPointAtIndex( ptwXPoints *ptwX, int64_t index, double x );
int ptwX_ascendingOrder( ptwXPoints *ptwX );

#if defined __cplusplus
    }
#endif

#endif          /* End of ptwX_h_included. */
