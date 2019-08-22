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
localCFLAGS = $(CFLAGS) $(nf_LegendreCFLAGS) -I$(NUMERICALFUNCTIONS_INCLUDE_PATH)
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
