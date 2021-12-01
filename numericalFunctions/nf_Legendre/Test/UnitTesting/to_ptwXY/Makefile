# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
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

TARGETS = to_ptwXY

.PHONY: default check clean realclean $(TARGETS)

default: $(TARGETS)

to_ptwXY: to_ptwXY.c
	$(CC) -o to_ptwXY to_ptwXY.c -lnumericalfunctions $(localCFLAGS) -lm

check: default
	for f in $(TARGETS); do ./$$f $(CHECKOPTIONS); done

clean:
	rm -f *.dat v e

realclean: clean
	rm -fr $(TARGETS) *.dSYM
