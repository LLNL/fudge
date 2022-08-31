/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
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
