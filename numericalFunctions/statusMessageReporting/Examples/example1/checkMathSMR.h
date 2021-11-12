/*
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <statusMessageReporting.h>
#include <mathSMR.h>

#define checkMathSMR_code_sin 1
#define checkMathSMR_code_exp 2

int checkMathSMR_setup( void );
int checkMathSMR_getLibrarysID( void );
int checkMathSMR_test( statusMessageReporting *smr, int verbose );
