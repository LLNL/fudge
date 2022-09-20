# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

NUMERICALFUNCTIONS_PATH = ../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include

ptwXYCFLAGS = -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -I$(SMR_PATH)/include
localCFLAGS = $(CFLAGS) $(ptwXYCFLAGS)
objects = ptwXY_core.o ptwXY_methods.o ptwXY_misc.o ptwXY_functions.o ptwXY_convenient.o ptwXY_integration.o \
    ptwXY_binaryOperators.o ptwXY_unitaryOperators.o ptwXY_interpolation.o

.PHONY: default clean realclean

default: libptwXY.a

libptwXY.a: $(objects)
	$(AR) -cr libptwXY.a *.o

ptwXY_core.o: ptwXY_core.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_core.o ptwXY_core.c

ptwXY_methods.o: ptwXY_methods.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_methods.o ptwXY_methods.c

ptwXY_functions.o: ptwXY_functions.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_functions.o ptwXY_functions.c

ptwXY_misc.o: ptwXY_misc.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_misc.o ptwXY_misc.c

ptwXY_convenient.o: ptwXY_convenient.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_convenient.o ptwXY_convenient.c

ptwXY_integration.o: ptwXY_integration.c ptwXY.h
	$(CC) $(localCFLAGS) -I./ -I../../nf_Legendre/Src -I../../nf_integration/Src -c -o ptwXY_integration.o ptwXY_integration.c

ptwXY_binaryOperators.o: ptwXY_binaryOperators.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_binaryOperators.o ptwXY_binaryOperators.c

ptwXY_unitaryOperators.o: ptwXY_unitaryOperators.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_unitaryOperators.o ptwXY_unitaryOperators.c

ptwXY_interpolation.o: ptwXY_interpolation.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_interpolation.o ptwXY_interpolation.c

clean:
	rm -f *.o *.a

realclean: clean
