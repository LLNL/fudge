# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
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

TARGETS = toLinearLinear unitbase unitbase2 chargedParticle strings badLogValue strings2

.PHONY: default check clean realclean $(TARGETS)

default: $(TARGETS)

toLinearLinear: toLinearLinear.c
	$(CC) -o toLinearLinear toLinearLinear.c $(localCFLAGS) -lm

unitbase: unitbase.c
	$(CC) -o unitbase unitbase.c $(localCFLAGS) -lm

unitbase2: unitbase2.c
	$(CC) -o unitbase2 unitbase2.c $(localCFLAGS) -lm

chargedParticle: chargedParticle.c
	$(CC) -o chargedParticle chargedParticle.c $(localCFLAGS) -lm

strings: strings.c
	$(CC) -o strings strings.c $(localCFLAGS) -lm

strings2: strings2.c
	$(CC) -o strings2 strings2.c $(localCFLAGS) -lm

badLogValue: badLogValue.c
	$(CC) -o badLogValue badLogValue.c $(localCFLAGS) -lm

check: default
	./test.com toLinearLinear $(CHECKOPTIONS)
	./test.com toLinearLinear $(CHECKOPTIONS) -xlog
	./test.com toLinearLinear $(CHECKOPTIONS) -ylog
	./test.com toLinearLinear $(CHECKOPTIONS) -xlog -ylog
	./test.com unitbase $(CHECKOPTIONS)
	./test.com unitbase2 $(CHECKOPTIONS)
	./test.com chargedParticle $(CHECKOPTIONS)
	./test.com strings.com $(CHECKOPTIONS)
	./test.com badLogValue $(CHECKOPTIONS)
	./test.com strings2 $(CHECKOPTIONS)

clean:
	rm -f *.dat v

realclean: clean
	rm -fr $(TARGETS) *.dSYM
