# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

SHELL = /bin/sh

NUMERICALFUNCTIONS_PATH = ../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include

localCFLAGS = $(CFLAGS) -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -I$(SMR_PATH)/include

objects = nf_Legendre.o nf_Legendre_GaussianQuadrature.o

.PHONY: default clean realclean

default: libnf_Legendre.a

libnf_Legendre.a: $(objects)
	$(AR) -cr libnf_Legendre.a *.o

nf_Legendre.o: nf_Legendre.c nf_Legendre.h
	$(CC) $(localCFLAGS) -c -o nf_Legendre.o nf_Legendre.c

nf_Legendre_GaussianQuadrature.o: nf_Legendre_GaussianQuadrature.c nf_Legendre.h
	$(CC) $(localCFLAGS) -c -o nf_Legendre_GaussianQuadrature.o nf_Legendre_GaussianQuadrature.c

clean:
	rm -f *.o *.a

realclean: clean
