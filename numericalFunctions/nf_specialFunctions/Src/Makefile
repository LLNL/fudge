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

objects = nf_exponentialIntegral.o nf_gammaFunctions.o nf_incompleteGammaFunctions.o nf_polevl.o nf_angularMomentumCoupling.o

.PHONY: default clean realclean

default: $(objects)

nf_exponentialIntegral.o: nf_exponentialIntegral.c
	$(CC) $(localCFLAGS) -c -o nf_exponentialIntegral.o nf_exponentialIntegral.c

nf_gammaFunctions.o: nf_gammaFunctions.c
	$(CC) $(localCFLAGS) -c -o nf_gammaFunctions.o nf_gammaFunctions.c

nf_incompleteGammaFunctions.o: nf_incompleteGammaFunctions.c
	$(CC) $(localCFLAGS) -c -o nf_incompleteGammaFunctions.o nf_incompleteGammaFunctions.c

nf_polevl.o: nf_polevl.c
	$(CC) $(localCFLAGS) -c -o nf_polevl.o nf_polevl.c

nf_angularMomentumCoupling.o: nf_angularMomentumCoupling.c
	$(CC) $(localCFLAGS) -c -o nf_angularMomentumCoupling.o nf_angularMomentumCoupling.c

clean:
	rm -f *.o *.a

realclean: clean

test_amc: test_amc.o nf_angularMomentumCoupling.o
	$(CC) $(localCFLAGS) -o test_amc nf_angularMomentumCoupling.o test_amc.o -lm

test_amc.o: test_amc.c_
	$(CC) $(localCFLAGS) -c -o test_amc.o test_amc.c_
