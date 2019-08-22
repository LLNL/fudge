/*
# <<BEGIN-copyright>>
# Copyright (c) 2016, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Nuclear Data and Theory group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-683960.
# All rights reserved.
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
#       notice, this list of conditions and the disclaimer below.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the disclaimer (as noted below) in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of LLNS/LLNL nor the names of its contributors may be used
#       to endorse or promote products derived from this software without specific
#       prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
# THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# 
# Additional BSD Notice
# 
# 1. This notice is required to be provided under our contract with the U.S.
# Department of Energy (DOE). This work was produced at Lawrence Livermore
# National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
# 
# 2. Neither the United States Government nor Lawrence Livermore National Security,
# LLC nor any of their employees, makes any warranty, express or implied, or assumes
# any liability or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that its use
# would not infringe privately-owned rights.
# 
# 3. Also, reference herein to any specific commercial products, process, or services
# by trade name, trademark, manufacturer or otherwise does not necessarily constitute
# or imply its endorsement, recommendation, or favoring by the United States Government
# or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed
# herein do not necessarily state or reflect those of the United States Government or
# Lawrence Livermore National Security, LLC, and shall not be used for advertising or
# product endorsement purposes.
# 
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
ptwXPoints *ptwX_new( statusMessageReporting *smr, int64_t size );
nfu_status ptwX_initialize( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t size );
ptwXPoints *ptwX_create( statusMessageReporting *smr, int64_t size, int64_t length, double const *xs );
ptwXPoints *ptwX_createLine( statusMessageReporting *smr, int64_t size, int64_t length, double slope, double offset );
nfu_status ptwX_copy( statusMessageReporting *smr, ptwXPoints *dest, ptwXPoints *src );
ptwXPoints *ptwX_clone( statusMessageReporting *smr, ptwXPoints *ptwX );
ptwXPoints *ptwX_slice( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t index1, int64_t index2 );
nfu_status ptwX_reallocatePoints( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t size, int forceSmallerResize );
nfu_status ptwX_clear( statusMessageReporting *smr, ptwXPoints *ptwX );
nfu_status ptwX_release( statusMessageReporting *smr, ptwXPoints *ptwX );
ptwXPoints *ptwX_free( ptwXPoints *ptwX );

int64_t ptwX_length( statusMessageReporting *smr, ptwXPoints *ptwX );
nfu_status ptwX_setData( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t length, double const *xs );
nfu_status ptwX_deletePoints( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t i1, int64_t i2 );
double *ptwX_getPointAtIndex( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t index );
double ptwX_getPointAtIndex_Unsafely( ptwXPoints *ptwX, int64_t index );
nfu_status ptwX_setPointAtIndex( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t index, double x );
nfu_status ptwX_insertPointsAtIndex( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t index, int64_t n1, double const *xs );
nfu_status ptwX_ascendingOrder( statusMessageReporting *smr, ptwXPoints *ptwX, int *order );
ptwXPoints *ptwX_fromString( statusMessageReporting *smr, char const *str, char sep, char **endCharacter );
int ptwX_countOccurrences( statusMessageReporting *smr, ptwXPoints *ptwX, double value );
nfu_status ptwX_reverse( statusMessageReporting *smr, ptwXPoints *ptwX );
nfu_status ptwX_sort( statusMessageReporting *smr, ptwXPoints *ptwX, enum ptwX_sort_order order );
nfu_status ptwX_closesDifference( statusMessageReporting *smr, ptwXPoints *ptwX, double value, int64_t *index, double *difference );
nfu_status ptwX_closesDifferenceInRange( statusMessageReporting *smr, ptwXPoints *ptwX, int64_t i1, int64_t i2, 
        double value, int64_t *index, double *difference );
ptwXPoints *ptwX_unique( statusMessageReporting *smr, ptwXPoints *ptwX, int order );

nfu_status ptwX_abs( statusMessageReporting *smr, ptwXPoints *ptwX );
nfu_status ptwX_neg( statusMessageReporting *smr, ptwXPoints *ptwX );
nfu_status ptwX_add_double( statusMessageReporting *smr, ptwXPoints *ptwX, double value );
nfu_status ptwX_mul_double( statusMessageReporting *smr, ptwXPoints *ptwX, double value );
nfu_status ptwX_slopeOffset( statusMessageReporting *smr, ptwXPoints *ptwX, double slope, double offset );
nfu_status ptwX_add_ptwX( statusMessageReporting *smr, ptwXPoints *ptwX1, ptwXPoints *ptwX2 );
nfu_status ptwX_sub_ptwX( statusMessageReporting *smr, ptwXPoints *ptwX1, ptwXPoints *ptwX2 );

nfu_status ptwX_range( statusMessageReporting *smr, ptwXPoints *ptwX, double *rangeMin, double *rangeMax );

nfu_status ptwX_compare( statusMessageReporting *smr, ptwXPoints *ptwX1, ptwXPoints *ptwX2, int *comparison );
nfu_status ptwX_close( statusMessageReporting *smr, ptwXPoints *ptwX1, ptwXPoints *ptwX2, int epsilonFactor, double epsilon,
        int *index );

/*
* Routines in ptwX_misc.c
*/
nfu_status ptwX_simpleWrite( statusMessageReporting *smr, ptwXPoints const *ptwX, FILE *f, char const *format );
nfu_status ptwX_simplePrint( statusMessageReporting *smr, ptwXPoints const *ptwX, char const *format );

#if defined __cplusplus
    }
#endif

#endif          /* End of ptwX_h_included. */
