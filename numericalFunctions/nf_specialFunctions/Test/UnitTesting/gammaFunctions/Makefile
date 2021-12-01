# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

TARGETS = gammaTest incompleteGammaTest

SMR_PATH = ../../../../statusMessageReporting
NUMERICALFUNCTIONS_PATH = ../../../..

ptwXYCFLAGS = -I$(NUMERICALFUNCTIONS_PATH)/include -L$(NUMERICALFUNCTIONS_PATH)/lib -lnumericalfunctions \
              -I$(SMR_PATH)/include -L$(SMR_PATH)/lib -lstatusMessageReporting \
              -lptwXY_utilities -lnfut_utilities
localCFLAGS = $(CFLAGS) $(ptwXYCFLAGS)

.PHONY: default clean realclean m gammaMathematica.m incompleteGammaMathematica.m

default: $(TARGETS)

check: default
	./gammaTest $(CHECKOPTIONS)
	./incompleteGammaTest $(CHECKOPTIONS)

m: gammaMathematica.m incompleteGammaMathematica.m

gammaTest: gammaTest.c
	$(CC) -o gammaTest gammaTest.c $(localCFLAGS) -lm

incompleteGammaTest: incompleteGammaTest.c
	$(CC) -o incompleteGammaTest incompleteGammaTest.c $(localCFLAGS) -lm

gammaMathematica.m:
	$(CC) $(CFLAGS) -o gammaMathematica gammaMathematica.c
	./gammaMathematica > gammaMathematica.m

incompleteGammaMathematica.m:
	$(CC) $(CFLAGS) -o incompleteGammaMathematica incompleteGammaMathematica.c
	./incompleteGammaMathematica > incompleteGammaMathematica.m

clean:
	rm -f gammaTest gammaMathematica incompleteGammaTest incompleteGammaMathematica *.m

realclean: clean
	rm -fr *.dSYM
