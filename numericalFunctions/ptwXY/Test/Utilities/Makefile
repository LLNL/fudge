# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

.PHONY: default libptwXY_utilities.a ptwXY_utilities.o utilities check clean realclean

ARFLAGS = rc

SMR_PATH = ../../../statusMessageReporting
NUMERICALFUNCTIONS_PATH = ../../..

localCFLAGS = $(CFLAGS) -I$(NUMERICALFUNCTIONS_PATH)/include -I$(SMR_PATH)/include

default: libptwXY_utilities.a
	cp ptwXY_utilities.h ../../../include
	cp libptwXY_utilities.a ../../../lib

libptwXY_utilities.a: ptwXY_utilities.o
	$(AR) $(ARFLAGS) $@ $^
	ranlib $@

ptwXY_utilities.o: ptwXY_utilities.c ptwXY_utilities.h
	$(CC) -c -o ptwXY_utilities.o ptwXY_utilities.c $(localCFLAGS)

utilities: default

check: default

clean:
	rm -f *.a *.o

realclean: clean
