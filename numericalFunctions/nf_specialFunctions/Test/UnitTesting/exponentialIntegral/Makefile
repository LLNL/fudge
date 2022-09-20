# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

TARGETS = exponentialIntegralTest

SMR_PATH = ../../../../statusMessageReporting
NUMERICALFUNCTIONS_PATH = ../../../..

ptwXYCFLAGS = -I$(NUMERICALFUNCTIONS_PATH)/include -L$(NUMERICALFUNCTIONS_PATH)/lib -lnumericalfunctions \
              -I$(SMR_PATH)/include -L$(SMR_PATH)/lib -lstatusMessageReporting \
              -lptwXY_utilities -lnfut_utilities
localCFLAGS = $(CFLAGS) $(ptwXYCFLAGS)

.PHONY: default clean realclean

default: $(TARGETS)

check: default
	./exponentialIntegralTest $(CHECKOPTIONS)

exponentialIntegralTest: exponentialIntegralTest.c
	$(CC) -o exponentialIntegralTest exponentialIntegralTest.c $(localCFLAGS) -lm

exponentialIntegralMathematica:
	$(CC) $(CFLAGS) -o exponentialIntegralMathematica exponentialIntegralMathematica.c
	./exponentialIntegralMathematica > exponentialIntegralMathematica.m

clean:
	rm -f exponentialIntegralTest exponentialIntegralMathematica exponentialIntegralMathematica.m

realclean: clean
	rm -rf *.dSYM
