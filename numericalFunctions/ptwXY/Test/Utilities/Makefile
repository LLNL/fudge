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

.PHONY: check clean realclean

CC = gcc
CFLAGS = -g -O0 -Wall
ARFLAGS = rc

NUMERICALFUNCTIONS_PATH = ../../..
NUMERICALFUNCTIONS_INCLUDE_PATH = $(NUMERICALFUNCTIONS_PATH)/include

PTWXY_PATH = ../..
PTWXY_INCLUDE_PATH = $(PTWXY_PATH)/include

localCFLAGS = -I$(NUMERICALFUNCTIONS_INCLUDE_PATH) -I$(PTWXY_INCLUDE_PATH)

libptwXY_utilities.a: ptwXY_utilities.o
	$(AR) $(ARFLAGS) $@ $^
	ranlib $@

ptwXY_utilities.o: ptwXY_utilities.c ptwXY_utilities.h
	$(CC) $(CFLAGS) -c -o ptwXY_utilities.o ptwXY_utilities.c $(localCFLAGS)

check: libptwXY_utilities.a

clean:
	rm -f *.a *.o

realclean: clean
