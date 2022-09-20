# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

NUMERICALFUNCTIONS_PATH = ../..

localCFLAGS = $(CFLAGS) -I$(SMR_PATH)/include
objects = nf_utilities.o nf_stringToInt32s.o nf_stringToDoubles.o

.PHONY: default clean realclean

default: libnf_utilities.a

libnf_utilities.a: $(objects)
	$(AR) -cr libnf_utilities.a *.o

nf_utilities.o: nf_utilities.c nf_utilities.h
	$(CC) $(localCFLAGS) -c -o nf_utilities.o nf_utilities.c

nf_stringToInt32s.o: nf_stringToInt32s.c nf_utilities.h
	$(CC) $(localCFLAGS) -c -o nf_stringToInt32s.o nf_stringToInt32s.c

nf_stringToDoubles.o: nf_stringToDoubles.c nf_utilities.h
	$(CC) $(localCFLAGS) -c -o nf_stringToDoubles.o nf_stringToDoubles.c

nf_stringToDoubles_main: nf_stringToDoubles_main.c nf_stringToDoubles.o nf_utilities.o nf_utilities.h
	$(CC) $(localCFLAGS) -o nf_stringToDoubles_main nf_stringToDoubles_main.c nf_stringToDoubles.o nf_utilities.o

clean:
	rm -f *.o *.a

realclean: clean
