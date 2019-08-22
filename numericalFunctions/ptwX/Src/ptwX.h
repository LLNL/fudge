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

#ifndef ptwX_h_included
#define ptwX_h_included

#include <stdio.h>
#include <stdint.h>

#include <nf_utilities.h>

#if defined __cplusplus
    extern "C" {
#endif

#define ptwX_minimumSize 10

enum ptwX_sort_order { ptwX_sort_order_descending, ptwX_sort_order_ascending };

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
nfu_status ptwX_deletePoints( ptwXPoints *ptwX, int64_t i1, int64_t i2 );
double *ptwX_getPointAtIndex( ptwXPoints *ptwX, int64_t index );
double ptwX_getPointAtIndex_Unsafely( ptwXPoints *ptwX, int64_t index );
nfu_status ptwX_setPointAtIndex( ptwXPoints *ptwX, int64_t index, double x );
nfu_status ptwX_insertPointsAtIndex( ptwXPoints *ptwX, int64_t index, int64_t n1, double *xs );
int ptwX_ascendingOrder( ptwXPoints *ptwX );
ptwXPoints *ptwX_fromString( char const *str, char **endCharacter, nfu_status *status );
nfu_status ptwX_countOccurrences( ptwXPoints *ptwX, double value, int *count );
nfu_status ptwX_reverse( ptwXPoints *ptwX );
nfu_status ptwX_sort( ptwXPoints *ptwX, enum ptwX_sort_order order );
nfu_status ptwX_closesDifference( ptwXPoints *ptwX, double value, int64_t *index, double *difference );
nfu_status ptwX_closesDifferenceInRange( ptwXPoints *ptwX, int64_t i1, int64_t i2, double value, int64_t *index, double *difference );
ptwXPoints *ptwX_unique( ptwXPoints *ptwX, int order, nfu_status *status );
nfu_status ptwX_abs( ptwXPoints *ptwX );
nfu_status ptwX_neg( ptwXPoints *ptwX );
nfu_status ptwX_compare( ptwXPoints *ptwX1, ptwXPoints *ptwX2, int *comparison );

#if defined __cplusplus
    }
#endif

#endif          /* End of ptwX_h_included. */
