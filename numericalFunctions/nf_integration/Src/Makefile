# <<BEGIN-copyright>>
# Copyright 2021, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

NUMERICALFUNCTIONS_PATH = ../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include

localCFLAGS = $(CFLAGS) -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -I$(SMR_PATH)/include
objects = nf_GnG_adaptiveQuadrature.o

.PHONY: default clean realclean

default: libnf_integration.a

libnf_integration.a: $(objects)
	$(AR) -cr libnf_integration.a *.o

nf_GnG_adaptiveQuadrature.o: nf_GnG_adaptiveQuadrature.c nf_integration.h
	$(CC) $(localCFLAGS) -c -o nf_GnG_adaptiveQuadrature.o nf_GnG_adaptiveQuadrature.c

clean:
	rm -f *.o *.a

realclean: clean
