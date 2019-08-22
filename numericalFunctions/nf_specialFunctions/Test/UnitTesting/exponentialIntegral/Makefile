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

TARGETS = exponentialIntegralTest

NUMERICALFUNCTIONS_PATH = ../../../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include
NUMERICALFUNCTIONS_LIB_PATH = $(NUMERICALFUNCTIONS_PATH)/lib


CFLAGS = -g -O0 -Wall
localCFLAGS = $(CFLAGS) $(NFUCFLAGS) -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -I./

.PHONY: default clean realclean

default: $(TARGETS)

check: default
	./exponentialIntegralTest $(CHECKOPTIONS)

exponentialIntegralTest: exponentialIntegralTest.c
	$(CC) $(localCFLAGS) -o exponentialIntegralTest exponentialIntegralTest.c -L$(NUMERICALFUNCTIONS_LIB_PATH) -lnumericalfunctions -lm

exponentialIntegralMathematica:
	$(CC) $(CFLAGS) -o exponentialIntegralMathematica exponentialIntegralMathematica.c
	./exponentialIntegralMathematica > exponentialIntegralMathematica.m

clean:
	rm -f exponentialIntegralTest exponentialIntegralMathematica exponentialIntegralMathematica.m

realclean: clean
