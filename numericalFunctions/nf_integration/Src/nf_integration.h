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

#ifndef nf_integration_h_included
#define nf_integration_h_included

#include <nf_utilities.h>
#include <nf_Legendre.h>

#if defined __cplusplus
    extern "C" {
#endif

#define nf_GnG_adaptiveQuadrature_MaxMaxDepth 20

typedef nfu_status (*nf_GnG_adaptiveQuadrature_callback)( nf_Legendre_GaussianQuadrature_callback integrandFunction, void *argList, double x1, 
    double x2, double *integral );

nfu_status nf_GnG_adaptiveQuadrature( nf_GnG_adaptiveQuadrature_callback quadratureFunction, nf_Legendre_GaussianQuadrature_callback integrandFunction, 
    void *argList, double x1, double x2, int maxDepth, double tolerance, double *integral, long *evaluations );

#if defined __cplusplus
    }
#endif

#endif          /* End of nf_integration_h_included. */

