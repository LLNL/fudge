# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

PTWC_PATH = ../../ptwC
PTWC_INCLUDE_PATH = $(PTWC_PATH)/include

localCFLAGS = $(CFLAGS) $(ptwXCFLAGS) -I$(PTWC_INCLUDE_PATH) -I$(SMR_PATH)/include
objects = ptwX_core.o ptwX_misc.o

.PHONY: default check clean realclean

default: libptwX.a

check:

libptwX.a: $(objects)
	$(AR) -cr libptwX.a *.o

ptwX_core.o: ptwX_core.c ptwX.h
	$(CC) $(localCFLAGS) -c -o ptwX_core.o ptwX_core.c

ptwX_misc.o: ptwX_misc.c ptwX.h
	$(CC) $(localCFLAGS) -c -o ptwX_misc.o ptwX_misc.c

clean:
	rm -f *.o *.a

realclean: clean
