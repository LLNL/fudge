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

#ifndef nf_Legendre_h_included
#define nf_Legendre_h_included

#include <nf_utilities.h>
#include <ptwXY.h>

#if defined __cplusplus
    extern "C" {
#endif

#define nf_Legendre_minMaxOrder 4
#define nf_Legendre_maxMaxOrder 64
#define nf_Legendre_sizeIncrement 8

typedef struct nf_Legendre_s nf_Legendre;

struct nf_Legendre_s {
    int maxOrder;
    int allocated;          /* Will never be less than nf_Legendre_minMaxOrder. */
    double *Cls;
};

typedef nfu_status (*nf_Legendre_GaussianQuadrature_callback)( double x, double *y, void *argList );

/*
* Methods in nf_Legendre.c
*/
nf_Legendre *nf_Legendre_new( int initialSize, int maxOrder, double *Cls, nfu_status *status );
nfu_status nf_Legendre_setup( nf_Legendre *nfL, int initialSize, int maxOrder );
nfu_status nf_Legendre_release( nf_Legendre *nfL );
nf_Legendre *nf_Legendre_free( nf_Legendre *nfL );
nf_Legendre *nf_Legendre_clone( nf_Legendre *nfL, nfu_status *status );
nfu_status nf_Legendre_reallocateCls( nf_Legendre *Legendre, int size, int forceSmallerResize );
int nf_Legendre_maxOrder( nf_Legendre *Legendre );
int nf_Legendre_allocated( nf_Legendre *Legendre );
double nf_Legendre_getCl( nf_Legendre *Legendre, int l, nfu_status *status );
nfu_status nf_Legendre_setCl( nf_Legendre *Legendre, int l, double Cl );
nfu_status nf_Legendre_normalize( nf_Legendre *Legendre );
double nf_Legendre_evauluateAtMu( nf_Legendre *nfL, double mu, nfu_status *status );
double nf_Legendre_PofL_atMu( int l, double mu );
ptwXYPoints *nf_Legendre_to_ptwXY( nf_Legendre *nfL, double accuracy, int biSectionMax, int checkForRoots, nfu_status *status );
nf_Legendre *nf_Legendre_from_ptwXY( ptwXYPoints *ptwXY, int maxOrder, nfu_status *status );

/*
* Methods in nf_Legendre_GaussianQuadrature.c
*/
nfu_status nf_Legendre_GaussianQuadrature( int degree, double x1, double x2, nf_Legendre_GaussianQuadrature_callback func, void *argList, double *integral );

#if defined __cplusplus
    }
#endif

#endif          /* End of nf_Legendre_h_included. */
