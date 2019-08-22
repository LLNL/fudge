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
