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
