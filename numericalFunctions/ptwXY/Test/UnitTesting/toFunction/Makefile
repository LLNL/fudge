# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

SMR_PATH = ../../../../statusMessageReporting
NUMERICALFUNCTIONS_PATH = ../../../..

ptwXYCFLAGS = -lnfut_utilities -lptwXY_utilities \
			  -I$(NUMERICALFUNCTIONS_PATH)/include -L$(NUMERICALFUNCTIONS_PATH)/lib -lnumericalfunctions \
              -I$(SMR_PATH)/include -L$(SMR_PATH)/lib -lstatusMessageReporting
localCFLAGS = $(CFLAGS) $(ptwXYCFLAGS)

TARGETS = createFromFunction

.PHONY: default check clean realclean $(TARGETS)

default: $(TARGETS)

createFromFunction: createFromFunction.c
	$(CC) -o createFromFunction createFromFunction.c $(localCFLAGS) -lm

check: default
	for f in $(TARGETS); do ./$$f $(CHECKOPTIONS); done

clean:
	rm -f *.dat v e

realclean: clean
	rm -fr $(TARGETS) *.dSYM
