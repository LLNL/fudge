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
