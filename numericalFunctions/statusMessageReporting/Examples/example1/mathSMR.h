/*
# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>
*/

#include <statusMessageReporting.h>

#define mathSMR_code_sin 1
#define mathSMR_code_exp 2

int mathSMR_setup( void );
int mathSMR_getLibrarysID( void );
double mathSMR_sin( statusMessageReporting *smr, double x1 );
double mathSMR_exp( statusMessageReporting *smr, double x1 );
