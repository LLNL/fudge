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

CFLAGS = -g -O0 -Wall
localCFLAGS = $(CFLAGS) $(NFUCFLAGS)
objects = nf_utilities.o nf_stringToDoubles.o

.PHONY: default clean realclean

default: libnf_utilities.a

libnf_utilities.a: $(objects)
	$(AR) -cr libnf_utilities.a *.o

nf_utilities.o: nf_utilities.c nf_utilities.h
	$(CC) $(localCFLAGS) -c -o nf_utilities.o nf_utilities.c

nf_stringToDoubles.o: nf_stringToDoubles.c nf_utilities.h
	$(CC) $(localCFLAGS) -c -o nf_stringToDoubles.o nf_stringToDoubles.c

nf_stringToDoubles_main: nf_stringToDoubles_main.c nf_stringToDoubles.o nf_utilities.o nf_utilities.h
	$(CC) $(localCFLAGS) -o nf_stringToDoubles_main nf_stringToDoubles_main.c nf_stringToDoubles.o nf_utilities.o

clean:
	rm -f *.o *.a

realclean: clean
