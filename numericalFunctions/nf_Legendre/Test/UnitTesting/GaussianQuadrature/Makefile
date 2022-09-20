# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

SMR_PATH = ../../../../statusMessageReporting
NUMERICALFUNCTIONS_PATH = ../../../..

ptwXYCFLAGS = -I$(NUMERICALFUNCTIONS_PATH)/include -L$(NUMERICALFUNCTIONS_PATH)/lib -lnumericalfunctions \
              -I$(SMR_PATH)/include -L$(SMR_PATH)/lib -lstatusMessageReporting \
              -lptwXY_utilities -lnfut_utilities
localCFLAGS = $(CFLAGS) $(ptwXYCFLAGS)

TARGETS = GaussianQuadrature

.PHONY: default check clean realclean $(TARGETS)

default: $(TARGETS)

GaussianQuadrature: GaussianQuadrature.c
	$(CC) -o GaussianQuadrature GaussianQuadrature.c -lnumericalfunctions $(localCFLAGS) -lm

check: default
	for f in $(TARGETS); do ./$$f $(CHECKOPTIONS); done

clean:
	rm -f *.dat v e

realclean: clean
	rm -fr $(TARGETS) *.dSYM
