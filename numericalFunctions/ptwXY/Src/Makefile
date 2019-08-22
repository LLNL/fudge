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
localCFLAGS = $(CFLAGS) $(ptwCFLAGS) $(ptwXYCFLAGS) -I$(NUMERICALFUNCTIONS_INCLUDE_PATH)
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
	$(CC) $(localCFLAGS) -c -o ptwXY_integration.o ptwXY_integration.c

ptwXY_binaryOperators.o: ptwXY_binaryOperators.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_binaryOperators.o ptwXY_binaryOperators.c

ptwXY_unitaryOperators.o: ptwXY_unitaryOperators.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_unitaryOperators.o ptwXY_unitaryOperators.c

ptwXY_interpolation.o: ptwXY_interpolation.c ptwXY.h
	$(CC) $(localCFLAGS) -c -o ptwXY_interpolation.o ptwXY_interpolation.c

clean:
	rm -f *.o *.a

realclean: clean
