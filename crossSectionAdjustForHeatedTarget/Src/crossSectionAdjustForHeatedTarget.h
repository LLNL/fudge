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

#ifndef crossSectionAdjustForHeatedTarget_h_included
#define crossSectionAdjustForHeatedTarget_h_included

#if defined __cplusplus
    extern "C" {
#endif

typedef enum crossSectionAdjustForHeatedTarget_limit { crossSectionAdjustForHeatedTarget_limit_one_over_v, crossSectionAdjustForHeatedTarget_limit_constant,
    crossSectionAdjustForHeatedTarget_limit_threshold } crossSectionAdjustForHeatedTarget_limit;

#define crossSectionAdjustForHeatedTarget_mode_all 1
#define crossSectionAdjustForHeatedTarget_mode_do_not_thin 2
#define crossSectionAdjustForHeatedTarget_mode_allEDomain 4

typedef
    struct crossSectionAdjustForHeatedTarget_info {
        int mode;
        int verbose;
        int InfoStats;
        int WarningStats;
        int ErrorStats;
        int resolutionStats;
        int bytes;
        int evaluations;
        double fInterpolation;
    } crossSectionAdjustForHeatedTarget_info;

typedef
    struct {
        int index, prior, next, thinnable;
        double E, cs;
    } E_cs_heated_point;

typedef
    struct {
        double mass_Ratio, T, csMax;
        crossSectionAdjustForHeatedTarget_limit lowerlimit, upperlimit;
        crossSectionAdjustForHeatedTarget_info *info;
        int n_pairs;
        int iStart_orig, iEnd_orig, iStart, iEnd;
        double *E_cs_in_orig, *E_cs_in;
        int n_done;
        int n_allocated;
        double dEMin;           /* Estimate of minimum dE required to meet fInterpolation. Set for each point of the inputted cross section data. */
        double fInterpolationSlope, fInterpolationOffset;
        E_cs_heated_point *E_cs;
        double *sqrtEi;
        double *SEi;
    } E_cs_heated_point_Info;


typedef enum crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode { crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode_default, 
    crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode_exact, crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode_taylor } 
    crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode;

int crossSectionAdjustForHeatedTarget_init( crossSectionAdjustForHeatedTarget_limit lowerlimit, crossSectionAdjustForHeatedTarget_limit upperlimit,
    crossSectionAdjustForHeatedTarget_info *info, double mass_Ratio, double T, double fInterpolation, int n_pairs, double *E_cs_in,
    E_cs_heated_point_Info *E_cs_Info );

int crossSectionAdjustForHeatedTarget( crossSectionAdjustForHeatedTarget_limit lowerlimit, crossSectionAdjustForHeatedTarget_limit upperlimit, 
    crossSectionAdjustForHeatedTarget_info *info, double EMin, double mass_Ratio, double T, double fInterpolation, int n_pairs, 
    double *E_cs_in, double **E_cs_out );
void crossSectionAdjustForHeatedTarget_heat_at_E( double E, E_cs_heated_point_Info *E_cs_Info, E_cs_heated_point *E_cs_point );
int crossSectionAdjustForHeatedTarget_integrate_xn_qauss( crossSectionAdjustForHeatedTarget_integrate_xn_qauss_mode Mode, double a, double b, double dxnerf[5] );

#if defined __cplusplus
    }
#endif

#endif /* End of crossSectionAdjustForHeatedTarget_h_included */
