# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: default check clean realclean

NUMERICALFUNCTIONS_PATH = ../../../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include
NUMERICALFUNCTIONS_LIB_PATH = $(NUMERICALFUNCTIONS_PATH)/lib

ptwXCFLAGS =  -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -L$(NUMERICALFUNCTIONS_LIB_PATH) -lnumericalfunctions    \
              -lnfut_utilities                                                                              \
              -I$(SMR_PATH)/include -L$(SMR_PATH)/lib -lstatusMessageReporting
localCFLAGS = $(CFLAGS) $(ptwXCFLAGS)

default: testCreate

testCreate: testCreate.c $(PTWXY_ARCHIVE_FILE)
	$(CC) -o testCreate testCreate.c $(localCFLAGS)

check: default
	./testCreate $(CHECKOPTIONS)

clean:
	rm -f testCreate

realclean: clean
	rm -fr *.dSYM
