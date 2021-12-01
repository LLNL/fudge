# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

SMR_PATH = ../../../statusMessageReporting
NUMERICALFUNCTIONS_PATH = ../..

ptwXYCFLAGS = -I$(NUMERICALFUNCTIONS_PATH)/include -L$(NUMERICALFUNCTIONS_PATH)/lib -lnumericalfunctions \
              -I$(SMR_PATH)/include -L$(SMR_PATH)/lib -lstatusMessageReporting
localCFLAGS = $(ptwXYCFLAGS)

.PHONY: clean realclean

default: E1 biSectionMax accuracy

E1: E1.c
	$(CC) $(CFLAGS) $(localCFLAGS) -o E1 E1.c $(Libs) -lm

biSectionMax: biSectionMax.c
	$(CC) $(CFLAGS) $(localCFLAGS) -o biSectionMax biSectionMax.c $(Libs) -lm

accuracy: accuracy.c
	$(CC) $(CFLAGS) $(localCFLAGS) -o accuracy accuracy.c $(Libs) -lm

clean:
	rm -f E1 biSectionMax accuracy

realclean: clean
