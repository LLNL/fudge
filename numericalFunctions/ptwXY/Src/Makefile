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
# When citing FUDGE, please use the following reference:
#   C.M. Mattoon, B.R. Beck, N.R. Patel, N.C. Summers, G.W. Hedstrom, D.A. Brown, "Generalized Nuclear Data: A New Structure (with Supporting Infrastructure) for Handling Nuclear Data", Nuclear Data Sheets, Volume 113, Issue 12, December 2012, Pages 3145-3171, ISSN 0090-3752, http://dx.doi.org/10. 1016/j.nds.2012.11.008
# 
# 
#     Please also read this link - Our Notice and Modified BSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Lawrence Livermore National Security, LLC. nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
