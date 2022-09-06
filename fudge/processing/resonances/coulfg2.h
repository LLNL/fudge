/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#ifndef coulfg2_h_included
#define coulfg2_h_included

#if defined __cplusplus
    extern "C" {
#endif

int Coulomb_FG( double XX, double ETA, double XLMIN, double XLMAX, double *FC, double *GC, int *M1 );

#if defined __cplusplus
    }
#endif

#endif              /* End of coulfg2_h_included. */
