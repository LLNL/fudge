# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

SMR_PATH = ../../../statusMessageReporting
NUMERICALFUNCTIONS_PATH = ../../..

ptwXYCFLAGS = -I$(NUMERICALFUNCTIONS_PATH)/include -L$(NUMERICALFUNCTIONS_PATH)/lib -lnumericalfunctions \
              -I$(SMR_PATH)/include -L$(SMR_PATH)/lib -lstatusMessageReporting \
              -lptwXY_utilities -lnfut_utilities
localCFLAGS = $(CFLAGS) $(ptwXYCFLAGS)

TARGETS = testBinaryUnion testOverwrite testAddRandomPoints testAddRandomPoints2 testInterpolation_1 \
          testMultiPly2_1 testDivide2_1 testExp_1

PYTHON = python

.PHONY: default check clean realclean

default:
	$(MAKE) $(TARGETS)

check: default
	for f in $(TARGETS); do $(PYTHON) check.py $$f; done

testBinaryUnion: testBinaryUnion.c
	$(CC) -o testBinaryUnion testBinaryUnion.c $(localCFLAGS) -lm

testOverwrite: testOverwrite.c
	$(CC) -o testOverwrite testOverwrite.c $(localCFLAGS) -lm

testAddRandomPoints2: testAddRandomPoints2.c
	$(CC) -o testAddRandomPoints2 testAddRandomPoints2.c $(localCFLAGS) -lm

testAddRandomPoints: testAddRandomPoints.c
	$(CC) -o testAddRandomPoints testAddRandomPoints.c $(localCFLAGS) -lm

testInterpolation_1: testInterpolation_1.c
	$(CC) -o testInterpolation_1 testInterpolation_1.c $(localCFLAGS) -lm

testMultiPly2_1: testMultiPly2_1.c
	$(CC) -o testMultiPly2_1 testMultiPly2_1.c $(localCFLAGS) -lm

testDivide2_1: testDivide2_1.c
	$(CC) -L. -o testDivide2_1 testDivide2_1.c $(localCFLAGS) -lm

testExp_1: testExp_1.c
	$(CC) -L. -o testExp_1 testExp_1.c $(localCFLAGS) -lm

clean:
	rm -f *.dat t_*

realclean: clean
	rm -fr $(TARGETS) *.dSYM
