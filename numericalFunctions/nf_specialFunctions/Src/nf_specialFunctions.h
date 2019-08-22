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

#ifndef specialFunctions_h_included
#define specialFunctions_h_included

#include <math.h>
#include <float.h>
#include <nf_utilities.h>

#ifdef WIN32
#define isfinite _finite
#define M_PI 3.141592653589793238463
#define INFINITY (DBL_MAX+DBL_MAX)
#endif

#if defined __cplusplus
    extern "C" {
#endif

double nf_polevl( double x, double coef[], int N );
double nf_p1evl( double x, double coef[], int N );
double nf_exponentialIntegral( int n, double x, nfu_status *status );
double nf_gammaFunction( double x, nfu_status *status );
double nf_logGammaFunction( double x, nfu_status *status );
double nf_incompleteGammaFunction( double a, double x, nfu_status *status );
double nf_incompleteGammaFunctionComplementary( double a, double x, nfu_status *status );

double  nf_amc_log_factorial( int );
double  nf_amc_factorial( int );
double  nf_amc_wigner_3j( int, int, int, int, int, int );
double  nf_amc_wigner_6j( int, int, int, int, int, int );
double  nf_amc_wigner_9j( int, int, int, int, int, int, int, int, int );
double  nf_amc_racah( int, int, int, int, int, int );
double  nf_amc_clebsh_gordan( int, int, int, int, int );
double  nf_amc_z_coefficient( int, int, int, int, int, int );
double  nf_amc_zbar_coefficient( int, int, int, int, int, int );
double  nf_amc_reduced_matrix_element( int, int, int, int, int, int, int );

#if defined __cplusplus
    }
#endif

#endif          /* End of ptwXY_h_included. */
