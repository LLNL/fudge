# <<BEGIN-copyright>>
# Copyright (c) 2011, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the LLNL Computational Nuclear Physics group
#         (email: mattoon1@llnl.gov)
# LLNL-CODE-494171 All rights reserved.
# 
# This file is part of the FUDGE package (For Updating Data and 
#         Generating Evaluations)
# 
# 
#     Please also read this link - Our Notice and GNU General Public License.
# 
# This program is free software; you can redistribute it and/or modify it under 
# the terms of the GNU General Public License (as published by the Free Software
# Foundation) version 2, dated June 1991.
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of 
# the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with 
# this program; if not, write to 
# 
# the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA
# <<END-copyright>>

NUMERICALFUNCTIONS_PATH = ../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include

CFLAGS = -g -O0 -Wall
localCFLAGS = $(CFLAGS) $(NFUCFLAGS) -I$(NUMERICALFUNCTIONS_INCLUDE_PATH)

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
