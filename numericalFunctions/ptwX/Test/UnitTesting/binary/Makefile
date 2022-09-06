# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: default check clean realclean doTarget

NUMERICALFUNCTIONS_PATH = ../../../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include
NUMERICALFUNCTIONS_LIB_PATH = $(NUMERICALFUNCTIONS_PATH)/lib

PTW1D_PATH = ../../..
PTW1D_INCLUDE_PATH = $(PTW1D_PATH)/include
PTW1D_LIB_PATH = $(PTW1D_PATH)/lib
PTW1D_ARCHIVE_FILE = $(PTW1D_LIB_PATH)/libptwX.a

ptwXCFLAGS = -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -L$(NUMERICALFUNCTIONS_LIB_PATH) -lnumericalfunctions    \
              -lnfut_utilities                                                                              \
              -I$(SMR_PATH)/include -L$(SMR_PATH)/lib -lstatusMessageReporting
localCFLAGS = $(CFLAGS) $(ptwXCFLAGS)

default: binary1

binary1: binary1.c $(PTWXY_ARCHIVE_FILE)
	$(CC) -o binary1 binary1.c $(localCFLAGS) -lm

check: default
	./binary1 $(CHECKOPTIONS)

clean:
	rm -f binary1

realclean: clean
	rm -fr *.dSYM
