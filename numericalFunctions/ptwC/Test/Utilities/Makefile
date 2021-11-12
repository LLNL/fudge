# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: default libnfut_utilities.a nfut_utilities.o utilities check clean realclean

ARFLAGS = rc

NUMERICALFUNCTIONS_PATH = ../../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include

localCFLAGS = $(CFLAGS) -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -I$(SMR_PATH)/include

default: libnfut_utilities.a 
	cp nfut_utilities.h ../../../include
	cp libnfut_utilities.a ../../../lib

libnfut_utilities.a: nfut_utilities.o
	$(AR) $(ARFLAGS) $@ $^
	ranlib $@

nfut_utilities.o: nfut_utilities.c nfut_utilities.h
	$(CC) -c -o nfut_utilities.o nfut_utilities.c $(localCFLAGS)

utilities: default

check: default

clean:
	rm -f *.a *.o

realclean: clean
