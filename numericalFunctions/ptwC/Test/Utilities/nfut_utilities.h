/*
# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#ifndef nfut_utility_h_included
#define nfut_utility_h_included

#include <nf_utilities.h>

#if defined __cplusplus
    extern "C" {
#endif

long nfut_charToLong( statusMessageReporting *smr, char const *msg, char const *stringValue );

int nfut_cmpDoubles( double d1, double d2, double espilon );

void nfut_printSMRError( statusMessageReporting *smr, char const *file, int line, char const *function, char const *fmt, ... );
void nfut_printSMRErrorExit( statusMessageReporting *smr, char const *file, int line, char const *function, char const *fmt, ... );

#define nfut_printSMRError2( smr, fmt, ... ) nfut_printSMRError( smr, __FILE__, __LINE__, __func__, fmt, __VA_ARGS__ )
#define nfut_printSMRError2p( smr, fmt )     nfut_printSMRError( smr, __FILE__, __LINE__, __func__, fmt )
#define nfut_printSMRErrorExit2( smr, fmt, ... ) nfut_printSMRErrorExit( smr, __FILE__, __LINE__, __func__, fmt, __VA_ARGS__ )
#define nfut_printSMRErrorExit2p( smr, fmt )     nfut_printSMRErrorExit( smr, __FILE__, __LINE__, __func__, fmt )

#if defined __cplusplus
    }
#endif

#endif              /* End of nfut_utility_h_included. */
