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
