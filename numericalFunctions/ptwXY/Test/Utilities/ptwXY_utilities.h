/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#ifndef nf_utility_h_included
#define nf_utility_h_included

#include <nf_utilities.h>
#include <ptwX.h>
#include <ptwXY.h>

#if defined __cplusplus
    extern "C" {
#endif

int nfu_cmpDoubles( double d1, double d2, double espilon );
int nfu_ptwXY_cmp( ptwXYPoints *p1, ptwXYPoints *p2, int verbose, double frac );
void nfu_printSMRError( statusMessageReporting *smr, char const *file, int line, char const *function, char const *fmt, ... );

#define nfu_printSMRError2( smr, fmt, ... ) nfu_printSMRError( smr, __FILE__, __LINE__, __func__, fmt, __VA_ARGS__ )
#define nfu_printSMRError2p( smr, fmt )     nfu_printSMRError( smr, __FILE__, __LINE__, __func__, fmt )

void nfu_printXYDataOnVerbosity( int verbose, ptwXYPoints *data );
void nfu_printXDataOnVerbosity( int verbose, ptwXPoints *data );

#if defined __cplusplus
    }
#endif

#endif              /* End of nf_utility_h_included. */
