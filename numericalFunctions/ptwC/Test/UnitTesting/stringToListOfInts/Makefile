# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

SMR_PATH = ../../../../statusMessageReporting
NUMERICALFUNCTIONS_PATH = ../../../..

ptwXYCFLAGS = -lnfut_utilities \
			  -I$(NUMERICALFUNCTIONS_PATH)/include -L$(NUMERICALFUNCTIONS_PATH)/lib -lnumericalfunctions \
              -I$(SMR_PATH)/include -L$(SMR_PATH)/lib -lstatusMessageReporting
localCFLAGS = $(CFLAGS) $(ptwXYCFLAGS)

TARGETS = stringToListOfInts

.PHONY: default check clean realclean $(TARGETS)

default: $(TARGETS)

stringToListOfInts: stringToListOfInts.c
	$(CC) -o stringToListOfInts stringToListOfInts.c $(localCFLAGS) -lm

check: default
	./stringToListOfInts $(CHECKOPTIONS)

clean:
	rm -f *.dat v

realclean: clean
	rm -fr $(TARGETS) *.dSYM
