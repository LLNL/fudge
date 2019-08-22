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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    double y1,  y2, ratio;

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
